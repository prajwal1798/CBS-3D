//=============================================================================
// CBS3D++_SI
//
// Persistent distributed-memory CBS iteration loop.
//
// Preprocessing, pressure sparsity construction and AMG setup are performed
// once. Each iteration then executes owner/ghost-consistent CBS Steps 1 to 4.
//=============================================================================

#include "cbs/solver/Solver.hpp"

#include "cbs/assembly/EnergyAssembly.hpp"
#include "cbs/assembly/MomentumAssembly.hpp"
#include "cbs/assembly/PressureAssembly.hpp"
#include "cbs/assembly/VelocityCorrectionAssembly.hpp"
#include "cbs/boundary/Boundary.hpp"
#include "cbs/linalg/PetscPersistentDistributedPressureSystem.hpp"
#include "cbs/parallel/Coloring.hpp"
#include "cbs/parallel/HaloExchange.hpp"
#include "cbs/preprocess/Preprocess.hpp"
#include "cbs/timestep/TimeStep.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#ifdef CBS3D_USE_MPI
#include <mpi.h>
#endif

namespace cbs
{
#if defined(CBS3D_USE_MPI) && defined(CBS3D_USE_PETSC)
    namespace
    {
        void check_loop_mpi(
            const int error_code,
            const char* operation)
        {
            if (error_code != MPI_SUCCESS)
            {
                throw std::runtime_error(
                    std::string("Solver::runDistributedLoop MPI failure in ")
                    + operation);
            }
        }

        bool touches_fluid(
            const CBSStateSI& s,
            const Int ip)
        {
            return
                (s.node_material_mask(ip) &
                 CBSStateSI::node_touches_fluid) != 0;
        }

        Real fixed_pressure_value(
            const CBSStateSI& s,
            const Int ip)
        {
            for (Int i = 1; i <= s.cfg.bc_fixed; ++i)
            {
                if (s.bc_list(i) == ip)
                {
                    return s.bc_values(i);
                }
            }

            throw std::runtime_error(
                "Distributed loop fixed pressure node is absent from bc_list");
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

        ThermalBoundaryState build_thermal_boundary_state(
            CBSStateSI& s)
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
                            "Distributed thermal boundary contains a non-finite "
                            "prescribed temperature");
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
                            "Distributed thermal boundary contains an invalid "
                            "heat-flux face area");
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
                        "Conflicting prescribed temperatures meet at a "
                        "distributed boundary node");
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

        void reject_unsupported_distributed_options(
            const CBSStateSI& s)
        {
            if (s.cfg.cbs_scheme != 1)
            {
                throw std::runtime_error(
                    "Distributed loop currently requires semi-implicit CBS");
            }

            if (s.cfg.ntime < 1)
            {
                throw std::runtime_error(
                    "Distributed loop requires ntime >= 1");
            }

            if (s.cfg.turbulence_on > 0)
            {
                throw std::runtime_error(
                    "Distributed Spalart-Allmaras transport is not yet enabled");
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

        Real compute_communicator_timestep(
            CBSStateSI& s,
            const Int iteration)
        {
            TimeStep::computeTimeStep(s, iteration);

            Real local_minimum =
                std::numeric_limits<Real>::max();

            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                local_minimum =
                    std::min(local_minimum, s.delte(ie));
            }

            Real global_minimum = 0.0;

            check_loop_mpi(
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
                    "Distributed loop obtained an invalid global timestep");
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

        void build_distributed_time_diagonals(
            CBSStateSI& s)
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

                const bool fluid_element =
                    s.mat_elem(ie) == 0;

                for (Int in = 1; in <= s.cfg.nep; ++in)
                {
                    const Int ip = s.intma(in, ie);
                    const Int coefficient_index =
                        (ie - 1) * s.cfg.nep + in;

                    const Real mass =
                        s.elcoe_e(coefficient_index);

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

                    s.elcoe2(ip) =
                        1.0 / momentum_diagonal(ip);
                }

                if (thermal_diagonal(ip) <= 0.0 ||
                    !std::isfinite(thermal_diagonal(ip)))
                {
                    throw std::runtime_error(
                        "Non-positive owned thermal time diagonal");
                }

                s.elcoe2p(ip) =
                    1.0 / thermal_diagonal(ip);
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

            check_loop_mpi(
                MPI_Allreduce(
                    &local.l2_squared,
                    &global.l2_squared,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM,
                    MPI_COMM_WORLD),
                "MPI_Allreduce continuity L2");

            check_loop_mpi(
                MPI_Allreduce(
                    &local.max_abs,
                    &global.max_abs,
                    1,
                    MPI_DOUBLE,
                    MPI_MAX,
                    MPI_COMM_WORLD),
                "MPI_Allreduce continuity maximum");

            check_loop_mpi(
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
            Real relative_velocity = 0.0;
            Real relative_pressure = 0.0;
            Real relative_temperature = 0.0;
            Real maximum_velocity = 0.0;
            Real maximum_velocity_correction = 0.0;
            Real minimum_temperature = 0.0;
            Real maximum_temperature = 0.0;
        };

        IterationMetrics evaluate_iteration_metrics(
            const CBSStateSI& s,
            const Real local_maximum_velocity_correction)
        {
            Real local_sums[6] =
                {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

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

                if (touches_fluid(s, ip))
                {
                    for (Int dim = 1; dim <= s.cfg.ndim; ++dim)
                    {
                        const Real current = s.unkno(dim, ip);
                        const Real difference =
                            current - s.unkn1(dim, ip);

                        local_sums[0] += difference * difference;
                        local_sums[1] += current * current;
                    }
                }

                local_maximum_velocity =
                    std::max(
                        local_maximum_velocity,
                        std::sqrt(u * u + v * v + w * w));

                if (touches_fluid(s, ip))
                {
                    if (!std::isfinite(s.pres(ip)))
                    {
                        throw std::runtime_error(
                            "Distributed loop produced non-finite pressure");
                    }

                    const Real pressure_difference =
                        s.pres(ip) - s.pres1(ip);

                    local_sums[2] +=
                        pressure_difference * pressure_difference;

                    local_sums[3] +=
                        s.pres(ip) * s.pres(ip);
                }

                if (s.cfg.temp_calc > 0)
                {
                    const Real temperature = s.temperature(ip);
                    const Real difference =
                        temperature - s.temperature1(ip);

                    if (!std::isfinite(temperature))
                    {
                        throw std::runtime_error(
                            "Distributed loop produced non-finite temperature");
                    }

                    local_sums[4] += difference * difference;
                    local_sums[5] += temperature * temperature;

                    local_minimum_temperature =
                        std::min(local_minimum_temperature, temperature);

                    local_maximum_temperature =
                        std::max(local_maximum_temperature, temperature);
                }
            }

            Real global_sums[6] = {};

            check_loop_mpi(
                MPI_Allreduce(
                    local_sums,
                    global_sums,
                    6,
                    MPI_DOUBLE,
                    MPI_SUM,
                    MPI_COMM_WORLD),
                "MPI_Allreduce iteration residual sums");

            IterationMetrics metrics;

            const auto relative_change = [](
                const Real numerator,
                const Real denominator)
            {
                if (denominator > 1.0e-300)
                {
                    return std::sqrt(numerator / denominator);
                }

                return std::sqrt(numerator);
            };

            metrics.relative_velocity =
                relative_change(global_sums[0], global_sums[1]);

            metrics.relative_pressure =
                relative_change(global_sums[2], global_sums[3]);

            if (s.cfg.temp_calc > 0)
            {
                metrics.relative_temperature =
                    relative_change(global_sums[4], global_sums[5]);
            }

            check_loop_mpi(
                MPI_Allreduce(
                    &local_maximum_velocity,
                    &metrics.maximum_velocity,
                    1,
                    MPI_DOUBLE,
                    MPI_MAX,
                    MPI_COMM_WORLD),
                "MPI_Allreduce maximum velocity");

            check_loop_mpi(
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
                check_loop_mpi(
                    MPI_Allreduce(
                        &local_minimum_temperature,
                        &metrics.minimum_temperature,
                        1,
                        MPI_DOUBLE,
                        MPI_MIN,
                        MPI_COMM_WORLD),
                    "MPI_Allreduce minimum temperature");

                check_loop_mpi(
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
    // Executes the persistent distributed CBS iteration loop.
    //=========================================================================
    void Solver::runDistributedLoop()
    {
#if defined(CBS3D_USE_MPI) && defined(CBS3D_USE_PETSC)
        if (!s_.mpi_enabled)
        {
            throw std::runtime_error(
                "Solver::runDistributedLoop requires more than one MPI rank");
        }

        const double setup_start = MPI_Wtime();

        runDistributedPreprocessing();

        Preprocess::classifyFaceEdges(s_);
        Preprocess::elementSize(s_);
        Coloring::build(s_);

        reject_unsupported_distributed_options(s_);

        PressureAssembly::buildElementPressureTerms(s_);

        ThermalBoundaryState thermal_boundary(s_.cfg.npoin);

        if (s_.cfg.temp_calc > 0)
        {
            thermal_boundary =
                build_thermal_boundary_state(s_);
        }

        Boundary::applyOwnedVelocityConstraints(s_);
        apply_owned_pressure_constraints(s_);

        if (s_.cfg.temp_calc > 0)
        {
            apply_owned_temperature_constraints(
                s_,
                thermal_boundary);
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

        PetscPersistentDistributedPressureSystem pressure_system;
        pressure_system.initialise(s_);

        double local_setup_seconds =
            MPI_Wtime() - setup_start;

        double global_setup_seconds = 0.0;

        check_loop_mpi(
            MPI_Allreduce(
                &local_setup_seconds,
                &global_setup_seconds,
                1,
                MPI_DOUBLE,
                MPI_MAX,
                MPI_COMM_WORLD),
            "MPI_Allreduce setup time");

        if (s_.mpi_rank == 0)
        {
            std::cout
                << "============================================================\n"
                << "CBS3D PERSISTENT DISTRIBUTED LOOP READY\n"
                << "============================================================\n"
                << "MPI ranks                  : " << s_.mpi_size << "\n"
                << "iterations requested       : " << s_.cfg.ntime << "\n"
                << "energy equation            : "
                << (s_.cfg.temp_calc > 0 ? "ON" : "OFF") << "\n"
                << "pressure matrix builds     : 1\n"
                << "AMG hierarchy builds       : 1\n"
                << "maximum setup time [s]     : "
                << global_setup_seconds << "\n"
                << "============================================================\n";
        }

        const double loop_start = MPI_Wtime();

        for (Int iteration = 1;
             iteration <= s_.cfg.ntime;
             ++iteration)
        {
            const double iteration_start = MPI_Wtime();

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
                compute_communicator_timestep(
                    s_,
                    iteration);

            if (s_.cfg.transient_on > 0)
            {
                s_.cfg.rtime += global_dt;
            }

            build_distributed_time_diagonals(s_);

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

                for (Int dim = 1;
                     dim <= s_.cfg.ndim;
                     ++dim)
                {
                    s_.unkno(dim, ip) =
                        s_.unkn1(dim, ip)
                        + s_.rhs(dim, ip) * s_.elcoe2(ip);
                }
            }

            Boundary::applyOwnedVelocityConstraints(s_);

            HaloExchange::broadcastOwnedToGhosts(
                s_.unkno,
                s_.partition_metadata,
                MPI_COMM_WORLD);

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

                for (Int dim = 1;
                     dim <= s_.cfg.ndim;
                     ++dim)
                {
                    const Real correction =
                        s_.rhs(dim, ip) * s_.elcoe2(ip);

                    if (!std::isfinite(correction))
                    {
                        throw std::runtime_error(
                            "Distributed loop produced a non-finite velocity "
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
                        s_.temperature1(ip)
                        + s_.rhs1(ip) * s_.elcoe2p(ip);
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
                evaluate_iteration_metrics(
                    s_,
                    local_maximum_velocity_correction);

            const Real continuity_l2 =
                std::sqrt(continuity.l2_squared);

            const Real continuity_rms =
                continuity.count > 0
                    ? continuity_l2 /
                        std::sqrt(
                            static_cast<Real>(continuity.count))
                    : 0.0;

            double local_iteration_seconds =
                MPI_Wtime() - iteration_start;

            double global_iteration_seconds = 0.0;

            check_loop_mpi(
                MPI_Allreduce(
                    &local_iteration_seconds,
                    &global_iteration_seconds,
                    1,
                    MPI_DOUBLE,
                    MPI_MAX,
                    MPI_COMM_WORLD),
                "MPI_Allreduce iteration time");

            if (s_.mpi_rank == 0)
            {
                std::cout
                    << "DD iteration " << iteration
                    << "/" << s_.cfg.ntime
                    << "  dt=" << global_dt
                    << "  CG=" << pressure_result.iterations
                    << "  CGrel=" << pressure_result.final_relative_l2
                    << "  RelU=" << metrics.relative_velocity
                    << "  RelP=" << metrics.relative_pressure;

                if (s_.cfg.temp_calc > 0)
                {
                    std::cout
                        << "  RelT=" << metrics.relative_temperature
                        << "  Tmin=" << metrics.minimum_temperature
                        << "  Tmax=" << metrics.maximum_temperature;
                }

                std::cout
                    << "  DivRMS=" << continuity_rms
                    << "  Umax=" << metrics.maximum_velocity
                    << "  dUmax="
                    << metrics.maximum_velocity_correction
                    << "  wall=" << global_iteration_seconds
                    << " s\n";
            }
        }

        double local_loop_seconds =
            MPI_Wtime() - loop_start;

        double global_loop_seconds = 0.0;

        check_loop_mpi(
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
                << "============================================================\n"
                << "CBS3D PERSISTENT DISTRIBUTED LOOP COMPLETE\n"
                << "============================================================\n"
                << "MPI ranks                  : " << s_.mpi_size << "\n"
                << "completed iterations       : " << s_.cfg.ntime << "\n"
                << "final physical time        : " << s_.cfg.rtime << "\n"
                << "pressure matrix builds     : 1\n"
                << "AMG hierarchy builds       : 1\n"
                << "maximum loop time [s]      : "
                << global_loop_seconds << "\n"
                << "next development stage     : distributed output/convergence\n"
                << "============================================================\n";
        }
#else
        throw std::runtime_error(
            "Solver::runDistributedLoop requires MPI and PETSc support");
#endif
    }
}
