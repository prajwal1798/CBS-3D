//=============================================================================
// CBS3D++_SI
//
// Production distributed-memory CBS loop with global convergence control and
// parallel VTU/PVTU/PVD output.
//
// Distributed finite-element sequence for every iteration:
//
//   forward old fields -> owned-element assembly -> reverse-add nodal residual
//   -> owner-only update/strong BC -> forward updated fields.
//
// The PETSc pressure topology, matrix, KSP and AMG hierarchy are built once and
// reused until the distributed calculation terminates.
//=============================================================================

#include "cbs/solver/Solver.hpp"

#include "cbs/assembly/EnergyAssembly.hpp"
#include "cbs/solver/Steps.hpp"
#include "cbs/turbulence/TurbulencePreprocess.hpp"
#include "cbs/assembly/MomentumAssembly.hpp"
#include "cbs/assembly/PressureAssembly.hpp"
#include "cbs/assembly/VelocityCorrectionAssembly.hpp"
#include "cbs/boundary/Boundary.hpp"
#include "cbs/io/DistributedPost.hpp"
#include "cbs/linalg/PetscPersistentDistributedPressureSystem.hpp"
#include "cbs/parallel/Coloring.hpp"
#include "cbs/parallel/HaloExchange.hpp"
#include "cbs/preprocess/Preprocess.hpp"
#include "cbs/timestep/TimeStep.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <streambuf>
#include <string>

#ifdef CBS3D_USE_MPI
#include <mpi.h>
#endif

namespace cbs
{
#if defined(CBS3D_USE_MPI) && defined(CBS3D_USE_PETSC)
    namespace
    {
        constexpr Real residual_epsilon = 1.0e-30;

        class NullStreamBuffer final : public std::streambuf
        {
        protected:
            int overflow(const int character) override
            {
                return traits_type::not_eof(character);
            }
        };


        class ScopedStdoutSilence final
        {
        public:
            explicit ScopedStdoutSilence(const bool suppress)
            {
                if (suppress)
                {
                    previous_ = std::cout.rdbuf(&null_buffer_);
                }
            }

            ~ScopedStdoutSilence()
            {
                if (previous_ != nullptr)
                {
                    std::cout.rdbuf(previous_);
                }
            }

            ScopedStdoutSilence(const ScopedStdoutSilence&) = delete;
            ScopedStdoutSilence& operator=(const ScopedStdoutSilence&) = delete;

        private:
            NullStreamBuffer null_buffer_;
            std::streambuf* previous_ = nullptr;
        };


        bool environment_flag_enabled(const char* name)
        {
            const char* value = std::getenv(name);

            return value != nullptr &&
                   value[0] != '\0' &&
                   std::string(value) != "0";
        }

        void check_mpi(const int error_code, const char* operation)
        {
            if (error_code != MPI_SUCCESS)
            {
                throw std::runtime_error(
                    std::string(
                        "Solver::runDistributedProductionLoop MPI failure in ")
                    + operation);
            }
        }

        bool touches_fluid(const CBSStateSI& s, const Int ip)
        {
            return
                (s.node_material_mask(ip) &
                 CBSStateSI::node_touches_fluid) != 0;
        }

        bool touches_solid(const CBSStateSI& s, const Int ip)
        {
            return
                (s.node_material_mask(ip) &
                 CBSStateSI::node_touches_solid) != 0;
        }

        bool velocity_active(const CBSStateSI& s, const Int ip)
        {
            return touches_fluid(s, ip) && !touches_solid(s, ip);
        }

        Real fixed_pressure_value(const CBSStateSI& s, const Int ip)
        {
            for (Int i = 1; i <= s.cfg.bc_fixed; ++i)
            {
                if (s.bc_list(i) == ip)
                {
                    const Real value = s.bc_values(i);

                    if (!std::isfinite(value))
                    {
                        throw std::runtime_error(
                            "Distributed fixed pressure is non-finite");
                    }

                    return value;
                }
            }

            throw std::runtime_error(
                "Distributed fixed pressure node is absent from bc_list");
        }

        bool temperature_dirichlet_value(
            const CBSStateSI& s,
            const Int bc,
            Real& value)
        {
            if (bc == s.cfg.bc_temperature_one_noslip)
            {
                value = 1.0;
                return true;
            }

            if (bc == s.cfg.bc_temperature_zero_noslip ||
                bc == s.cfg.bc_temperature_zero_prescribed_velocity ||
                bc == s.cfg.bc_parabolic_inlet)
            {
                value = 0.0;
                return true;
            }

            if (bc == s.cfg.bc_velocity_temperature_inlet ||
                bc == s.cfg.bc_massflow_temperature_inlet)
            {
                value = s.cfg.inlet_temperature;
                return true;
            }

            return false;
        }

        struct ThermalBoundaryState
        {
            explicit ThermalBoundaryState(const Int npoin)
                : fixed(npoin),
                  value(npoin),
                  heat_flux_load(npoin)
            {
                fixed.fill(0);
                value.fill(0.0);
                heat_flux_load.fill(0.0);
            }

            Array1D<Int> fixed;
            Array1D<Real> value;
            Array1D<Real> heat_flux_load;
        };

        ThermalBoundaryState build_thermal_boundary_state(CBSStateSI& s)
        {
            ThermalBoundaryState state(s.cfg.npoin);

            Array1D<Real> count(s.cfg.npoin);
            Array1D<Real> sum(s.cfg.npoin);
            Array1D<Real> sum_squared(s.cfg.npoin);

            count.fill(0.0);
            sum.fill(0.0);
            sum_squared.fill(0.0);

            for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
            {
                const Int bc = s.iside(s.cfg.bsid, ib);
                Real prescribed_temperature = 0.0;

                if (temperature_dirichlet_value(
                        s,
                        bc,
                        prescribed_temperature))
                {
                    if (!std::isfinite(prescribed_temperature))
                    {
                        throw std::runtime_error(
                            "Distributed prescribed temperature is non-finite");
                    }

                    for (Int in = 1; in <= s.cfg.nsidp; ++in)
                    {
                        const Int ip = s.iside(in, ib);

                        count(ip) += 1.0;
                        sum(ip) += prescribed_temperature;
                        sum_squared(ip) +=
                            prescribed_temperature * prescribed_temperature;
                    }
                }

                if (bc == s.cfg.bc_noslip_heatflux_wall &&
                    s.cfg.heat_flux_bc != 0.0)
                {
                    const Real area = s.face_norm(4, ib);

                    if (area <= 0.0 || !std::isfinite(area))
                    {
                        throw std::runtime_error(
                            "Distributed heat-flux boundary has invalid area");
                    }

                    const Real contribution =
                        s.cfg.heat_flux_bc * area /
                        static_cast<Real>(s.cfg.nsidp);

                    for (Int in = 1; in <= s.cfg.nsidp; ++in)
                    {
                        state.heat_flux_load(
                            s.iside(in, ib)) += contribution;
                    }
                }
            }

            HaloExchange::sumGhostContributionsToOwners(
                count,
                s.partition_metadata,
                MPI_COMM_WORLD);

            HaloExchange::sumGhostContributionsToOwners(
                sum,
                s.partition_metadata,
                MPI_COMM_WORLD);

            HaloExchange::sumGhostContributionsToOwners(
                sum_squared,
                s.partition_metadata,
                MPI_COMM_WORLD);

            HaloExchange::sumGhostContributionsToOwners(
                state.heat_flux_load,
                s.partition_metadata,
                MPI_COMM_WORLD);

            for (const Int ip : s.owned_nodes)
            {
                if (count(ip) <= 0.0)
                {
                    continue;
                }

                const Real mean = sum(ip) / count(ip);
                const Real variance =
                    std::max(
                        0.0,
                        sum_squared(ip) / count(ip) - mean * mean);

                const Real tolerance =
                    1.0e-12 * std::max(1.0, mean * mean);

                if (variance > tolerance)
                {
                    throw std::runtime_error(
                        "Conflicting distributed prescribed temperatures meet "
                        "at one global node");
                }

                state.fixed(ip) = 1;
                state.value(ip) = mean;
            }

            HaloExchange::broadcastOwnedToGhosts(
                state.fixed,
                s.partition_metadata,
                MPI_COMM_WORLD);

            HaloExchange::broadcastOwnedToGhosts(
                state.value,
                s.partition_metadata,
                MPI_COMM_WORLD);

            return state;
        }

        void apply_owned_temperature_constraints(
            CBSStateSI& s,
            const ThermalBoundaryState& boundary)
        {
            for (const Int ip : s.owned_nodes)
            {
                if (boundary.fixed(ip) != 0)
                {
                    s.temperature(ip) = boundary.value(ip);
                }
            }
        }

        void apply_owned_pressure_constraints(CBSStateSI& s)
        {
            for (const Int ip : s.owned_nodes)
            {
                if (s.node_pressure_fixed(ip) != 0)
                {
                    s.pres(ip) = fixed_pressure_value(s, ip);
                }
            }
        }

        void reject_unsupported_options(const CBSStateSI& s)
        {
            if (s.cfg.cbs_scheme != 1)
            {
                throw std::runtime_error(
                    "Distributed production loop requires semi-implicit CBS");
            }

            if (s.cfg.ntime < 1)
            {
                throw std::runtime_error(
                    "Distributed production loop requires ntime >= 1");
            }

            // Spalart-Allmaras is supported in the distributed loop.  The
            // wall-distance search gathers the global wall surface, the SA
            // classification is reduced over shared nodes, and the transport
            // step exchanges its element-assembled accumulators before the nodal
            // update.  Only the non-dimensional path is rejected: SA computes an
            // eddy viscosity from cfg.ani there, but MomentumAssembly ignores
            // mu_eff_e unless dimensional material properties are active, so the
            // turbulence would be computed and then silently discarded.
            if (s.cfg.turbulence_on > 0 &&
                !(s.cfg.dimensional_mode > 0 &&
                  s.cfg.material_properties_enabled > 0))
            {
                throw std::runtime_error(
                    "Spalart-Allmaras requires dimensional_mode and "
                    "material_properties_enabled; in non-dimensional mode the "
                    "momentum assembly ignores the turbulent viscosity");
            }

            if (s.cfg.step2_check > 0)
            {
                throw std::runtime_error(
                    "Distributed pressure-based timestep correction is not yet "
                    "implemented");
            }

            if (s.cfg.htype == 2)
            {
                throw std::runtime_error(
                    "Distributed htype=2 requires a shared-node minimum-length "
                    "halo operation");
            }

            for (const Int ip : s.owned_nodes)
            {
                if (s.node_symmetry(ip) != 0)
                {
                    throw std::runtime_error(
                        "Distributed symmetry requires a persistent nodal "
                        "symmetry-normal operator");
                }
            }
        }

        //=====================================================================
        // Establishes the timestep field for one distributed iteration.
        //
        // Two regimes are supported.
        //
        // Global uniform timestep.
        //     Every element and node is advanced with the smallest stable
        //     timestep anywhere in the mesh.  This is required whenever the
        //     pseudo-time path is physically meaningful, that is for a
        //     transient run, and it is what a fixed dtfix requests.
        //
        // Local timestep.
        //     Each element is advanced at its own stable rate.  For a
        //     steady-state calculation the path through pseudo-time carries no
        //     meaning and only the fixed point is sought, so advancing each
        //     element at its own rate is a legitimate relaxation and converges
        //     enormously faster on a stretched mesh: on a wall-resolved
        //     flat-plate grid the near-wall cells are hundreds of times smaller
        //     than the freestream cells, so a global timestep throttles the
        //     entire domain to the rate of its smallest element and the flow
        //     needs millions of iterations merely to cross the domain once.
        //
        // The distributed subtlety is that the nodal timestep deltp is shared
        // across partition interfaces.  TimeStep::computeTimeStep forms it from
        // locally visible elements only, so each rank's value at an interface
        // node is a minimum over a subset and is therefore too large.  Left
        // uncorrected the two ranks would advance the same node with different
        // dt, and since deltp enters the Step 2 pressure operator and the Step 3
        // correction, the assembled system would depend on the partitioning.
        // The nodal field is therefore reduced by minimum onto the owners and
        // broadcast back to the ghosts before use.
        //
        // Element timesteps need no exchange: an element belongs to exactly one
        // rank.
        //=====================================================================
        bool distributed_local_timestep_enabled(const CBSStateSI& s)
        {
            // Mirror the serial ilots semantics exactly rather than inventing a
            // second rule, so that a case behaves the same way in both paths.
            if (s.cfg.ilots != 1 && s.cfg.ilots != 2)
            {
                return false;
            }

            if (s.cfg.transient_on > 0 ||
                s.cfg.solver_opt > 1 ||
                s.cfg.dtfixed > 0 ||
                s.cfg.dtfixed == -1 ||
                s.cfg.ilots <= -1)
            {
                return false;
            }

            return true;
        }

        Real compute_communicator_timestep(
            CBSStateSI& s,
            const Int iteration)
        {
            TimeStep::computeTimeStep(s, iteration);

            if (distributed_local_timestep_enabled(s))
            {
                // computeTimeStep has just written the currently stable
                // timestep of every element and node.  Treat that as the
                // candidate field, reconcile the nodal part across partition
                // interfaces, and then decide whether the frozen field that the
                // pressure operator was built from is still safe.
                HaloExchange::minGhostContributionsToOwners(
                    s.deltp,
                    s.partition_metadata,
                    MPI_COMM_WORLD);

                HaloExchange::broadcastOwnedToGhosts(
                    s.deltp,
                    s.partition_metadata,
                    MPI_COMM_WORLD);

                HaloExchange::minGhostContributionsToOwners(
                    s.deltp2,
                    s.partition_metadata,
                    MPI_COMM_WORLD);

                HaloExchange::broadcastOwnedToGhosts(
                    s.deltp2,
                    s.partition_metadata,
                    MPI_COMM_WORLD);

                // Stability monitor.
                //
                // The frozen field is unsafe as soon as any element's frozen
                // timestep exceeds a margin fraction of the timestep that
                // element can currently sustain.  The margin is below one so
                // that the refresh happens before the frozen value actually
                // violates stability rather than after.  This is the mechanism
                // that keeps the method safe while nu_tilde grows and the SA
                // diffusion limit tightens, which is exactly the regime in which
                // a blindly frozen field would eventually diverge.
                Int refresh_needed = s.lts_frozen_valid ? 0 : 1;

                if (s.lts_frozen_valid)
                {
                    const Real margin = s.cfg.lts_refresh_margin > 0.0
                        ? s.cfg.lts_refresh_margin
                        : 0.8;

                    for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
                    {
                        if (s.lts_delte_frozen(ie) >
                            margin * s.delte(ie))
                        {
                            refresh_needed = 1;
                            break;
                        }
                    }
                }

                // Every rank must reach the same decision, or some would rebuild
                // the pressure operator while others did not and the next
                // collective solve would deadlock.
                Int global_refresh = 0;

                check_mpi(
                    MPI_Allreduce(
                        &refresh_needed,
                        &global_refresh,
                        1,
                        MPI_INT,
                        MPI_MAX,
                        MPI_COMM_WORLD),
                    "MPI_Allreduce local timestep refresh decision");

                if (global_refresh != 0)
                {
                    for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
                    {
                        s.lts_delte_frozen(ie) = s.delte(ie);
                    }

                    for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
                    {
                        s.lts_deltp_frozen(ip) = s.deltp(ip);
                        s.lts_deltp2_frozen(ip) = s.deltp2(ip);
                    }

                    s.lts_frozen_valid = true;
                    ++s.cfg.lts_refresh_count;
                }

                // Install the frozen field.  Everything downstream, including
                // the momentum time diagonals and the pressure operator, must
                // use the same timestep the operator was assembled with.
                for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
                {
                    s.delte(ie) = s.lts_delte_frozen(ie);
                }

                for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
                {
                    s.deltp(ip) = s.lts_deltp_frozen(ip);
                    s.deltp1(ip) = s.lts_deltp_frozen(ip);
                    s.deltp2(ip) = s.lts_deltp2_frozen(ip);
                }

                Real local_minimum =
                    std::numeric_limits<Real>::max();

                Real local_maximum = 0.0;

                for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
                {
                    const Real dt = s.delte(ie);

                    if (dt <= 0.0 || !std::isfinite(dt))
                    {
                        throw std::runtime_error(
                            "Distributed production loop obtained an invalid "
                            "local element timestep");
                    }

                    local_minimum = std::min(local_minimum, dt);
                    local_maximum = std::max(local_maximum, dt);
                }

                Real extrema[2] = { -local_minimum, local_maximum };
                Real global_extrema[2] = { 0.0, 0.0 };

                check_mpi(
                    MPI_Allreduce(
                        extrema,
                        global_extrema,
                        2,
                        MPI_DOUBLE,
                        MPI_MAX,
                        MPI_COMM_WORLD),
                    "MPI_Allreduce local timestep extrema");

                const Real global_minimum = -global_extrema[0];

                if (global_minimum <= 0.0 || !std::isfinite(global_minimum))
                {
                    throw std::runtime_error(
                        "Distributed production loop obtained invalid timestep");
                }

                s.cfg.dtreal = global_minimum;
                s.cfg.lts_dt_min = global_minimum;
                s.cfg.lts_dt_max = global_extrema[1];

                s.lts_operator_stale = global_refresh != 0;

                return global_minimum;
            }

            Real local_minimum =
                std::numeric_limits<Real>::max();

            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                local_minimum =
                    std::min(local_minimum, s.delte(ie));
            }

            Real global_minimum = 0.0;

            check_mpi(
                MPI_Allreduce(
                    &local_minimum,
                    &global_minimum,
                    1,
                    MPI_DOUBLE,
                    MPI_MIN,
                    MPI_COMM_WORLD),
                "MPI_Allreduce global timestep");

            if (global_minimum <= 0.0 ||
                !std::isfinite(global_minimum))
            {
                throw std::runtime_error(
                    "Distributed production loop obtained invalid timestep");
            }

            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                s.delte(ie) = global_minimum;
            }

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                s.deltp(ip) = global_minimum;
                s.deltp1(ip) = global_minimum;
                s.deltp2(ip) = global_minimum;
            }

            s.cfg.dtreal = global_minimum;
            return global_minimum;
        }

        void build_distributed_time_diagonals(CBSStateSI& s)
        {
            Array1D<Real> momentum_diagonal(s.cfg.npoin);
            Array1D<Real> thermal_diagonal(s.cfg.npoin);

            momentum_diagonal.fill(0.0);
            thermal_diagonal.fill(0.0);

            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                const Real dt = s.delte(ie);

                if (dt <= 0.0 || !std::isfinite(dt))
                {
                    throw std::runtime_error(
                        "Invalid element timestep in distributed diagonal assembly");
                }

                const bool fluid_element = s.mat_elem(ie) == 0;

                for (Int in = 1; in <= s.cfg.nep; ++in)
                {
                    const Int ip = s.intma(in, ie);
                    const Int coefficient_index =
                        (ie - 1) * s.cfg.nep + in;

                    const Real mass = s.elcoe_e(coefficient_index);

                    if (mass <= 0.0 || !std::isfinite(mass))
                    {
                        throw std::runtime_error(
                            "Invalid element-node mass in distributed diagonal "
                            "assembly");
                    }

                    if (fluid_element)
                    {
                        momentum_diagonal(ip) += mass / dt;
                    }

                    thermal_diagonal(ip) +=
                        s.rho_cp_e(ie) * mass / dt;
                }
            }

            HaloExchange::sumGhostContributionsToOwners(
                momentum_diagonal,
                s.partition_metadata,
                MPI_COMM_WORLD);

            HaloExchange::sumGhostContributionsToOwners(
                thermal_diagonal,
                s.partition_metadata,
                MPI_COMM_WORLD);

            s.elcoe2.fill(0.0);
            s.elcoe2p.fill(0.0);

            for (const Int ip : s.owned_nodes)
            {
                if (touches_fluid(s, ip))
                {
                    if (momentum_diagonal(ip) <= 0.0 ||
                        !std::isfinite(momentum_diagonal(ip)))
                    {
                        throw std::runtime_error(
                            "Non-positive owned momentum time diagonal");
                    }

                    s.elcoe2(ip) = 1.0 / momentum_diagonal(ip);
                }

                if (thermal_diagonal(ip) <= 0.0 ||
                    !std::isfinite(thermal_diagonal(ip)))
                {
                    throw std::runtime_error(
                        "Non-positive owned thermal time diagonal");
                }

                s.elcoe2p(ip) = 1.0 / thermal_diagonal(ip);
            }

            HaloExchange::broadcastOwnedToGhosts(
                s.elcoe2,
                s.partition_metadata,
                MPI_COMM_WORLD);

            HaloExchange::broadcastOwnedToGhosts(
                s.elcoe2p,
                s.partition_metadata,
                MPI_COMM_WORLD);
        }

        struct ContinuityMetrics
        {
            Real l2_squared = 0.0;
            Real max_abs = 0.0;
            long long count = 0;
        };

        ContinuityMetrics global_preprojection_continuity(
            const CBSStateSI& s)
        {
            ContinuityMetrics local;

            for (const Int ip : s.owned_nodes)
            {
                if (!touches_fluid(s, ip) ||
                    s.node_pressure_fixed(ip) != 0)
                {
                    continue;
                }

                const Real value = s.rhs1(ip);

                if (!std::isfinite(value))
                {
                    throw std::runtime_error(
                        "Non-finite distributed continuity residual");
                }

                local.l2_squared += value * value;
                local.max_abs =
                    std::max(local.max_abs, std::abs(value));
                ++local.count;
            }

            ContinuityMetrics global;

            check_mpi(
                MPI_Allreduce(
                    &local.l2_squared,
                    &global.l2_squared,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM,
                    MPI_COMM_WORLD),
                "MPI_Allreduce continuity L2");

            check_mpi(
                MPI_Allreduce(
                    &local.max_abs,
                    &global.max_abs,
                    1,
                    MPI_DOUBLE,
                    MPI_MAX,
                    MPI_COMM_WORLD),
                "MPI_Allreduce continuity maximum");

            check_mpi(
                MPI_Allreduce(
                    &local.count,
                    &global.count,
                    1,
                    MPI_LONG_LONG,
                    MPI_SUM,
                    MPI_COMM_WORLD),
                "MPI_Allreduce continuity count");

            return global;
        }

        struct IterationMetrics
        {
            Real relative_u = 0.0;
            Real relative_v = 0.0;
            Real relative_w = 0.0;
            Real relative_velocity = 0.0;
            Real relative_pressure = 0.0;
            Real relative_temperature = 0.0;
            Real maximum_velocity = 0.0;
            Real maximum_velocity_correction = 0.0;
            Real minimum_temperature = 0.0;
            Real maximum_temperature = 0.0;
        };

        Real rhs_scale_residual(
            const Real difference,
            const Real inverse_diagonal,
            const char* field_name)
        {
            if (inverse_diagonal <= 0.0 ||
                !std::isfinite(inverse_diagonal))
            {
                throw std::runtime_error(
                    std::string("Invalid inverse diagonal in distributed ")
                    + field_name + " convergence");
            }

            return difference / inverse_diagonal;
        }

        void finalise_residual_triplet(
            CBSStateSI& s,
            const std::array<Real, 15>& sums,
            const Int base)
        {
            s.hb[static_cast<std::size_t>(base)] =
                std::sqrt(
                    sums[static_cast<std::size_t>(base)] /
                    (sums[static_cast<std::size_t>(base + 1)] +
                     residual_epsilon));

            s.hb[static_cast<std::size_t>(base + 1)] =
                std::sqrt(sums[static_cast<std::size_t>(base + 1)]);

            s.hb[static_cast<std::size_t>(base + 2)] =
                std::sqrt(sums[static_cast<std::size_t>(base + 2)]);
        }

        IterationMetrics evaluate_global_convergence(
            CBSStateSI& s,
            const Real local_maximum_velocity_correction)
        {
            std::array<Real, 15> local_sums{};
            std::array<Real, 15> global_sums{};

            Real local_maximum_velocity = 0.0;
            Real local_minimum_temperature =
                std::numeric_limits<Real>::max();
            Real local_maximum_temperature =
                -std::numeric_limits<Real>::max();

            for (const Int ip : s.owned_nodes)
            {
                const Real u = s.unkno(1, ip);
                const Real v = s.unkno(2, ip);
                const Real w = s.unkno(3, ip);

                if (!std::isfinite(u) ||
                    !std::isfinite(v) ||
                    !std::isfinite(w))
                {
                    throw std::runtime_error(
                        "Distributed loop produced non-finite velocity");
                }

                if (velocity_active(s, ip))
                {
                    for (Int dim = 1; dim <= s.cfg.ndim; ++dim)
                    {
                        const Int base = (dim - 1) * 3;
                        const Real current = s.unkno(dim, ip);
                        const Real difference =
                            current - s.unkn1(dim, ip);

                        const Real rhs_residual =
                            rhs_scale_residual(
                                difference,
                                s.elcoe2(ip),
                                "velocity");

                        local_sums[static_cast<std::size_t>(base)] +=
                            difference * difference;

                        local_sums[static_cast<std::size_t>(base + 1)] +=
                            current * current;

                        local_sums[static_cast<std::size_t>(base + 2)] +=
                            rhs_residual * rhs_residual;
                    }

                    local_maximum_velocity =
                        std::max(
                            local_maximum_velocity,
                            std::sqrt(u * u + v * v + w * w));
                }

                if (touches_fluid(s, ip))
                {
                    const Real pressure = s.pres(ip);
                    const Real pressure_difference =
                        pressure - s.pres1(ip);

                    if (!std::isfinite(pressure))
                    {
                        throw std::runtime_error(
                            "Distributed loop produced non-finite pressure");
                    }

                    local_sums[9] +=
                        pressure_difference * pressure_difference;
                    local_sums[10] += pressure * pressure;
                    local_sums[11] +=
                        pressure_difference * pressure_difference;
                }

                if (s.cfg.temp_calc > 0)
                {
                    const Real temperature = s.temperature(ip);
                    const Real temperature_difference =
                        temperature - s.temperature1(ip);

                    if (!std::isfinite(temperature))
                    {
                        throw std::runtime_error(
                            "Distributed loop produced non-finite temperature");
                    }

                    const Real rhs_residual =
                        rhs_scale_residual(
                            temperature_difference,
                            s.elcoe2p(ip),
                            "temperature");

                    local_sums[12] +=
                        temperature_difference * temperature_difference;
                    local_sums[13] += temperature * temperature;
                    local_sums[14] += rhs_residual * rhs_residual;

                    local_minimum_temperature =
                        std::min(local_minimum_temperature, temperature);
                    local_maximum_temperature =
                        std::max(local_maximum_temperature, temperature);
                }
            }

            check_mpi(
                MPI_Allreduce(
                    local_sums.data(),
                    global_sums.data(),
                    static_cast<int>(global_sums.size()),
                    MPI_DOUBLE,
                    MPI_SUM,
                    MPI_COMM_WORLD),
                "MPI_Allreduce distributed convergence sums");

            finalise_residual_triplet(s, global_sums, 0);
            finalise_residual_triplet(s, global_sums, 3);
            finalise_residual_triplet(s, global_sums, 6);
            finalise_residual_triplet(s, global_sums, 9);

            if (s.cfg.temp_calc > 0)
            {
                finalise_residual_triplet(s, global_sums, 12);
            }
            else
            {
                s.hb[12] = 0.0;
                s.hb[13] = 0.0;
                s.hb[14] = 0.0;
            }

            IterationMetrics metrics;
            metrics.relative_u = s.hb[0];
            metrics.relative_v = s.hb[3];
            metrics.relative_w = s.hb[6];
            metrics.relative_velocity =
                std::max({
                    metrics.relative_u,
                    metrics.relative_v,
                    metrics.relative_w});
            metrics.relative_pressure = s.hb[9];
            metrics.relative_temperature = s.hb[12];

            check_mpi(
                MPI_Allreduce(
                    &local_maximum_velocity,
                    &metrics.maximum_velocity,
                    1,
                    MPI_DOUBLE,
                    MPI_MAX,
                    MPI_COMM_WORLD),
                "MPI_Allreduce maximum velocity");

            check_mpi(
                MPI_Allreduce(
                    &local_maximum_velocity_correction,
                    &metrics.maximum_velocity_correction,
                    1,
                    MPI_DOUBLE,
                    MPI_MAX,
                    MPI_COMM_WORLD),
                "MPI_Allreduce maximum velocity correction");

            if (s.cfg.temp_calc > 0)
            {
                check_mpi(
                    MPI_Allreduce(
                        &local_minimum_temperature,
                        &metrics.minimum_temperature,
                        1,
                        MPI_DOUBLE,
                        MPI_MIN,
                        MPI_COMM_WORLD),
                    "MPI_Allreduce minimum temperature");

                check_mpi(
                    MPI_Allreduce(
                        &local_maximum_temperature,
                        &metrics.maximum_temperature,
                        1,
                        MPI_DOUBLE,
                        MPI_MAX,
                        MPI_COMM_WORLD),
                    "MPI_Allreduce maximum temperature");
            }

            return metrics;
        }

        bool steady_state_reached(const CBSStateSI& s)
        {
            const Real velocity_residual =
                std::max({s.hb[0], s.hb[3], s.hb[6]});

            const bool velocity_ok =
                s.cfg.vel_check < 1 ||
                velocity_residual < s.cfg.l2norm_vel_tolerance;

            const bool temperature_ok =
                s.cfg.temp_check < 1 ||
                s.cfg.temp_calc < 1 ||
                s.hb[12] < s.cfg.l2norm_temp_tolerance;

            const bool pressure_ok =
                s.cfg.l2norm_pres_tolerance <= 0.0 ||
                s.hb[9] < s.cfg.l2norm_pres_tolerance;

            return velocity_ok && temperature_ok && pressure_ok;
        }

        bool transient_end_time_reached(const CBSStateSI& s)
        {
            if (s.cfg.transient_on < 1 || s.cfg.end_rtime <= 0.0)
            {
                return false;
            }

            const Real tolerance =
                1.0e-12 * std::max(1.0, std::abs(s.cfg.end_rtime));

            return s.cfg.rtime + tolerance >= s.cfg.end_rtime;
        }

        class ScopedHeatFluxSuppression
        {
        public:
            explicit ScopedHeatFluxSuppression(CBSStateSI& s)
                : state_(s),
                  saved_value_(s.cfg.heat_flux_bc)
            {
                state_.cfg.heat_flux_bc = 0.0;
            }

            ~ScopedHeatFluxSuppression()
            {
                state_.cfg.heat_flux_bc = saved_value_;
            }

            ScopedHeatFluxSuppression(
                const ScopedHeatFluxSuppression&) = delete;

            ScopedHeatFluxSuppression& operator=(
                const ScopedHeatFluxSuppression&) = delete;

        private:
            CBSStateSI& state_;
            Real saved_value_;
        };
    }
#endif

    //=========================================================================
    // Executes DD-4B/DD-4C production distributed loop.
    //=========================================================================
    void Solver::runDistributedProductionLoop()
    {
#if defined(CBS3D_USE_MPI) && defined(CBS3D_USE_PETSC)
        if (!s_.mpi_enabled)
        {
            throw std::runtime_error(
                "Solver::runDistributedProductionLoop requires more than one "
                "MPI rank");
        }

        const double setup_start = MPI_Wtime();

        ThermalBoundaryState thermal_boundary(s_.cfg.npoin);
        PetscPersistentDistributedPressureSystem pressure_system;

        {
            const ScopedStdoutSilence setup_stdout(
                s_.mpi_rank != 0 ||
                !environment_flag_enabled("CBS3D_VERBOSE"));

            runDistributedPreprocessing();

            Preprocess::classifyFaceEdges(s_);
            Preprocess::elementSize(s_);
            Coloring::build(s_);

            reject_unsupported_options(s_);

            PressureAssembly::buildElementPressureTerms(s_);

            if (s_.cfg.temp_calc > 0)
            {
                thermal_boundary = build_thermal_boundary_state(s_);
            }

            Boundary::applyOwnedVelocityConstraints(s_);
            apply_owned_pressure_constraints(s_);

            if (s_.cfg.temp_calc > 0)
            {
                apply_owned_temperature_constraints(s_, thermal_boundary);
            }

            HaloExchange::broadcastOwnedToGhosts(
                s_.unkno,
                s_.partition_metadata,
                MPI_COMM_WORLD);

            HaloExchange::broadcastOwnedToGhosts(
                s_.pres,
                s_.partition_metadata,
                MPI_COMM_WORLD);

            if (s_.cfg.temp_calc > 0)
            {
                HaloExchange::broadcastOwnedToGhosts(
                    s_.temperature,
                    s_.partition_metadata,
                    MPI_COMM_WORLD);
            }

            updateVelocityMagnitude();

            // Wall distance, SA node classification and the initial eddy
            // viscosity.  This must run before the first Step 1, because the
            // momentum assembly reads mu_eff_e.
            if (s_.mpi_rank == 0)
            {
                if (distributed_local_timestep_enabled(s_))
                {
                    std::cout
                        << "  Timestep regime : LOCAL (ilots="
                        << s_.cfg.ilots << ")\n"
                        << "    Each element advances at its own stable rate."
                           " Valid because this is a\n"
                           "    steady-state run, where the path through"
                           " pseudo-time carries no meaning\n"
                           "    and only the fixed point is sought. Reported dt"
                           " is the global minimum.\n";
                }
                else
                {
                    std::cout
                        << "  Timestep regime : GLOBAL UNIFORM\n"
                        << "    Every element advances at the smallest stable"
                           " timestep in the mesh.\n"
                           "    On a stretched wall-resolved mesh this throttles"
                           " the whole domain to the\n"
                           "    rate of its smallest cell. For a steady run set"
                           " ilots=1 and dtfixed=0\n"
                           "    to enable local timestepping.\n";
                }
            }

            TurbulencePreprocess::prepareSpalartAllmaras(s_);

            pressure_system.initialise(s_);
        }

        double local_setup_seconds = MPI_Wtime() - setup_start;
        double global_setup_seconds = 0.0;

        check_mpi(
            MPI_Allreduce(
                &local_setup_seconds,
                &global_setup_seconds,
                1,
                MPI_DOUBLE,
                MPI_MAX,
                MPI_COMM_WORLD),
            "MPI_Allreduce setup time");

        const std::string distributed_case_name =
            DistributedPost::distributedCaseName(case_name_);

        if (s_.mpi_rank == 0)
        {
            std::cout
                << "Run configuration\n"
                << "  mesh       : "
                << s_.partition_metadata.global_nelem << " tetrahedra, "
                << s_.partition_metadata.global_npoin << " nodes\n"
                << "  parallel   : " << s_.mpi_size << " MPI ranks\n"
                << "  iterations : " << s_.cfg.ntime << "\n"
                << "  energy     : "
                << (s_.cfg.temp_calc > 0 ? "ON" : "OFF") << "\n"
                << "  residuals  : "
                << (s_.cfg.residual_log_enabled > 0 ? "ON" : "OFF") << "\n"
                << "  VTU/PVD    : "
                << (s_.cfg.vtu_output_enabled > 0 ? "ON" : "OFF") << "\n"
                << "  setup time : " << global_setup_seconds << " s\n";
        }

        DistributedPost::initialise(s_, case_name_);

        const double loop_start = MPI_Wtime();

        Int last_iteration = 0;
        std::string stop_reason = "maximum iteration count reached";

        for (Int iteration = 1;
             iteration <= s_.cfg.ntime;
             ++iteration)
        {
            const double iteration_start = MPI_Wtime();
            last_iteration = iteration;

            s_.unkn1 = s_.unkno;
            s_.pres1 = s_.pres;
            s_.temperature1 = s_.temperature;
            s_.velocity_old = s_.velocity;

            HaloExchange::broadcastOwnedToGhosts(
                s_.unkn1,
                s_.partition_metadata,
                MPI_COMM_WORLD);

            if (s_.cfg.temp_calc > 0)
            {
                HaloExchange::broadcastOwnedToGhosts(
                    s_.temperature1,
                    s_.partition_metadata,
                    MPI_COMM_WORLD);
            }

            const Real global_dt =
                compute_communicator_timestep(s_, iteration);

            if (s_.cfg.transient_on > 0)
            {
                s_.cfg.rtime += global_dt;
            }

            // The pressure operator is assembled from the frozen element
            // timesteps, so it must be rewritten whenever that field changes.
            if (s_.lts_operator_stale)
            {
                pressure_system.rebuildOperator(s_);
                s_.lts_operator_stale = false;

                if (s_.mpi_rank == 0)
                {
                    std::cout
                        << "  [LTS] timestep field refreshed at iteration "
                        << iteration
                        << " (refresh " << s_.cfg.lts_refresh_count
                        << "); pressure operator and preconditioner rebuilt\n";
                }
            }

            build_distributed_time_diagonals(s_);

            // CBS Step 1: local owned-element momentum assembly, reverse-add,
            // owner-only predictor update and strong velocity constraints.
            MomentumAssembly::assembleStep1Rhs(s_);

            HaloExchange::sumGhostContributionsToOwners(
                s_.rhs,
                s_.partition_metadata,
                MPI_COMM_WORLD);

            for (const Int ip : s_.owned_nodes)
            {
                if (!touches_fluid(s_, ip))
                {
                    continue;
                }

                for (Int dim = 1; dim <= s_.cfg.ndim; ++dim)
                {
                    s_.unkno(dim, ip) =
                        s_.unkn1(dim, ip) +
                        s_.rhs(dim, ip) * s_.elcoe2(ip);
                }
            }

            Boundary::applyOwnedVelocityConstraints(s_);

            HaloExchange::broadcastOwnedToGhosts(
                s_.unkno,
                s_.partition_metadata,
                MPI_COMM_WORLD);

            // CBS Step 2: distributed continuity RHS and persistent PETSc solve.
            PressureAssembly::assembleStep2Rhs(s_);

            HaloExchange::sumGhostContributionsToOwners(
                s_.rhs1,
                s_.partition_metadata,
                MPI_COMM_WORLD);

            const ContinuityMetrics continuity =
                global_preprojection_continuity(s_);

            const ConjugateGradient::Result pressure_result =
                pressure_system.solve(s_);

            s_.last_cg_iterations = pressure_result.iterations;
            s_.last_cg_initial_l2 = pressure_result.initial_l2;
            s_.last_cg_final_l2 = pressure_result.final_l2;
            s_.last_cg_relative_l2 =
                pressure_result.final_relative_l2;
            s_.last_cg_max_abs = pressure_result.final_max_abs;

            if (!pressure_result.converged)
            {
                throw std::runtime_error(
                    "Persistent distributed pressure solve failed to converge");
            }

            HaloExchange::broadcastOwnedToGhosts(
                s_.pres,
                s_.partition_metadata,
                MPI_COMM_WORLD);

            // CBS Step 3: pressure-gradient correction, reverse-add, owner-only
            // velocity update and strong constraints.
            VelocityCorrectionAssembly::assembleStep3Rhs(s_);

            HaloExchange::sumGhostContributionsToOwners(
                s_.rhs,
                s_.partition_metadata,
                MPI_COMM_WORLD);

            Real local_maximum_velocity_correction = 0.0;

            for (const Int ip : s_.owned_nodes)
            {
                if (!touches_fluid(s_, ip))
                {
                    continue;
                }

                for (Int dim = 1; dim <= s_.cfg.ndim; ++dim)
                {
                    const Real correction =
                        s_.rhs(dim, ip) * s_.elcoe2(ip);

                    if (!std::isfinite(correction))
                    {
                        throw std::runtime_error(
                            "Distributed loop produced non-finite velocity "
                            "correction");
                    }

                    s_.unkno(dim, ip) += correction;

                    local_maximum_velocity_correction =
                        std::max(
                            local_maximum_velocity_correction,
                            std::abs(correction));
                }
            }

            Boundary::applyOwnedVelocityConstraints(s_);

            HaloExchange::broadcastOwnedToGhosts(
                s_.unkno,
                s_.partition_metadata,
                MPI_COMM_WORLD);

            // Step SA: Spalart-Allmaras transport.
            //
            // It runs after Step 3 because the SA advection and vorticity
            // production need the corrected velocity, and before Step 4 because
            // the energy assembly reads k_eff_e.  Steps::stepSpalartAllmaras
            // performs its own interface exchanges: the element-assembled
            // accumulators are summed onto their owners before the nodal update
            // and nu_tilde is broadcast back to the ghost layer afterwards.
            if (s_.cfg.turbulence_on > 0)
            {
                s_.nu_tilde1 = s_.nu_tilde;

                Steps::stepSpalartAllmaras(s_);
            }

            // CBS Step 4: thermal assembly over fluid and solid elements.
            if (s_.cfg.temp_calc > 0)
            {
                {
                    ScopedHeatFluxSuppression suppress_heat_flux(s_);
                    EnergyAssembly::assembleStep4Rhs(s_);
                }

                HaloExchange::sumGhostContributionsToOwners(
                    s_.rhs1,
                    s_.partition_metadata,
                    MPI_COMM_WORLD);

                for (const Int ip : s_.owned_nodes)
                {
                    s_.rhs1(ip) +=
                        thermal_boundary.heat_flux_load(ip);

                    s_.temperature(ip) =
                        s_.temperature1(ip) +
                        s_.rhs1(ip) * s_.elcoe2p(ip);
                }

                apply_owned_temperature_constraints(
                    s_,
                    thermal_boundary);

                HaloExchange::broadcastOwnedToGhosts(
                    s_.temperature,
                    s_.partition_metadata,
                    MPI_COMM_WORLD);
            }

            updateVelocityMagnitude();

            const IterationMetrics metrics =
                evaluate_global_convergence(
                    s_,
                    local_maximum_velocity_correction);

            // Collective, so it must be executed by every rank on every
            // iteration rather than only on the ranks that print.
            const Convergence::TurbulenceDiagnostics turbulence_diagnostics =
                Convergence::turbulenceDiagnostics(s_);

            const Real continuity_l2 =
                std::sqrt(continuity.l2_squared);

            const Real continuity_rms =
                continuity.count > 0
                    ? continuity_l2 /
                        std::sqrt(static_cast<Real>(continuity.count))
                    : 0.0;

            double local_iteration_seconds =
                MPI_Wtime() - iteration_start;
            double global_iteration_seconds = 0.0;

            check_mpi(
                MPI_Allreduce(
                    &local_iteration_seconds,
                    &global_iteration_seconds,
                    1,
                    MPI_DOUBLE,
                    MPI_MAX,
                    MPI_COMM_WORLD),
                "MPI_Allreduce iteration time");

            DistributedPost::writeResidualRow(
                s_,
                case_name_,
                iteration,
                continuity_rms,
                continuity.max_abs,
                metrics.maximum_velocity,
                metrics.maximum_velocity_correction,
                global_iteration_seconds,
                turbulence_diagnostics);

            if (DistributedPost::shouldWriteSolution(s_, iteration))
            {
                DistributedPost::writeSolution(
                    s_,
                    case_name_,
                    iteration);
            }

            bool stop_now = false;

            if (transient_end_time_reached(s_))
            {
                stop_reason = "transient end time reached";
                stop_now = true;
            }
            else if (s_.cfg.transient_on < 1 &&
                     iteration >= s_.cfg.steady_min_iterations &&
                     steady_state_reached(s_))
            {
                stop_reason = "steady-state convergence reached";
                stop_now = true;
            }

            const Int console_interval =
                std::max<Int>(1, s_.cfg.console_log_every);

            const bool print_iteration =
                iteration == 1 ||
                (iteration % console_interval) == 0 ||
                stop_now ||
                iteration == s_.cfg.ntime;

            if (s_.mpi_rank == 0 && print_iteration)
            {
                std::cout
                    << "Iteration " << iteration
                    << "/" << s_.cfg.ntime
                    << "  dt=" << global_dt
                    << "  CG=" << pressure_result.iterations
                    << "  CGrel=" << pressure_result.final_relative_l2
                    << "  RelU=" << metrics.relative_u
                    << "  RelV=" << metrics.relative_v
                    << "  RelW=" << metrics.relative_w
                    << "  RelP=" << metrics.relative_pressure;

                if (s_.cfg.temp_calc > 0)
                {
                    std::cout
                        << "  RelT=" << metrics.relative_temperature
                        << "  Tmin=" << metrics.minimum_temperature
                        << "  Tmax=" << metrics.maximum_temperature;
                }

                if (distributed_local_timestep_enabled(s_))
                {
                    const Real ratio = s_.cfg.lts_dt_min > 0.0
                        ? s_.cfg.lts_dt_max / s_.cfg.lts_dt_min
                        : 0.0;

                    std::cout
                        << "  LTSdt=[" << s_.cfg.lts_dt_min
                        << "," << s_.cfg.lts_dt_max << "]"
                        << "  LTSratio=" << ratio
                        << "  LTSrefresh=" << s_.cfg.lts_refresh_count;
                }

                if (s_.cfg.turbulence_on > 0)
                {
                    std::cout
                        << "  RelSA=" << turbulence_diagnostics.residual
                        << "  NuTilde=[" << turbulence_diagnostics.nu_tilde_min
                        << "," << turbulence_diagnostics.nu_tilde_max << "]"
                        << "  NuT=[" << turbulence_diagnostics.nu_t_min
                        << "," << turbulence_diagnostics.nu_t_max << "]"
                        << "  MuTmax=" << turbulence_diagnostics.mu_t_max
                        << "  MuEffmax=" << turbulence_diagnostics.mu_eff_max;
                }

                std::cout
                    << "  DivRMS=" << continuity_rms
                    << "  Umax=" << metrics.maximum_velocity
                    << "  dUmax="
                    << metrics.maximum_velocity_correction
                    << "  wall=" << global_iteration_seconds
                    << " s\n";
            }

            if (stop_now)
            {
                break;
            }
        }

        if (s_.cfg.vtu_output_enabled > 0 &&
            !DistributedPost::solutionAlreadyWritten(
                s_,
                last_iteration))
        {
            DistributedPost::writeSolution(
                s_,
                case_name_,
                last_iteration);
        }

        double local_loop_seconds = MPI_Wtime() - loop_start;
        double global_loop_seconds = 0.0;

        check_mpi(
            MPI_Allreduce(
                &local_loop_seconds,
                &global_loop_seconds,
                1,
                MPI_DOUBLE,
                MPI_MAX,
                MPI_COMM_WORLD),
            "MPI_Allreduce loop time");

        if (s_.mpi_rank == 0)
        {
            std::cout
                << "Run complete\n"
                << "  iterations    : " << last_iteration << "\n"
                << "  stop reason   : " << stop_reason << "\n"
                << "  physical time : " << s_.cfg.rtime << "\n"
                << "  loop time     : " << global_loop_seconds << " s\n"
                << "  output case   : " << distributed_case_name << "\n";
        }
#else
        throw std::runtime_error(
            "Solver::runDistributedProductionLoop requires MPI and PETSc support");
#endif
    }
}
