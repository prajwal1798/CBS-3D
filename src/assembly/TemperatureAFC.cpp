//=============================================================================
// CBS3D++_SI
//
// Distributed algebraic flux correction for the Step-4 temperature equation.
//=============================================================================

#include "cbs/assembly/TemperatureAFC.hpp"

#include "cbs/parallel/HaloExchange.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef CBS3D_USE_MPI
#include <mpi.h>
#endif

namespace cbs
{
#if defined(CBS3D_USE_MPI)
    namespace
    {
        struct EdgeFlux
        {
            Int node_i = 0;
            Int node_j = 0;

            // Antidiffusive residual flux added to node_i and removed from
            // node_j when its limiter coefficient is one.
            Real flux_to_i = 0.0;
        };


        void check_mpi(const int error_code, const char* operation)
        {
            if (error_code != MPI_SUCCESS)
            {
                throw std::runtime_error(
                    std::string("TemperatureAFC MPI failure in ") + operation);
            }
        }


        bool environment_flag_enabled(const char* name)
        {
            const char* value = std::getenv(name);

            return value != nullptr &&
                   value[0] != '\0' &&
                   std::string(value) != "0";
        }


        bool environment_real(
            const char* name,
            Real& result)
        {
            const char* value = std::getenv(name);

            if (value == nullptr || value[0] == '\0')
            {
                return false;
            }

            const std::string text(value);
            std::size_t parsed = 0;

            try
            {
                result = std::stod(text, &parsed);
            }
            catch (const std::exception&)
            {
                throw std::runtime_error(
                    std::string(name) + " is not a valid real number");
            }

            if (parsed != text.size() || !std::isfinite(result))
            {
                throw std::runtime_error(
                    std::string(name) + " must be finite");
            }

            return true;
        }


        Int environment_positive_int(
            const char* name,
            const Int fallback)
        {
            const char* value = std::getenv(name);

            if (value == nullptr || value[0] == '\0')
            {
                return fallback;
            }

            const std::string text(value);
            std::size_t parsed = 0;
            long long number = 0;

            try
            {
                number = std::stoll(text, &parsed, 10);
            }
            catch (const std::exception&)
            {
                throw std::runtime_error(
                    std::string(name) + " is not a valid integer");
            }

            if (parsed != text.size() ||
                number < 1 ||
                number > std::numeric_limits<Int>::max())
            {
                throw std::runtime_error(
                    std::string(name) + " must be a positive CBS3D Int");
            }

            return static_cast<Int>(number);
        }


        Int dNkdx_index(
            const CBSStateSI& s,
            const Int ie,
            const Int dim,
            const Int local_node)
        {
            return (ie - 1) * s.cfg.ndim * s.cfg.nep
                + (dim - 1) * s.cfg.nep
                + local_node;
        }


        Real grad(
            const CBSStateSI& s,
            const Int ie,
            const Int dim,
            const Int local_node)
        {
            return s.dNkdx(dNkdx_index(s, ie, dim, local_node));
        }


        bool is_fluid_element(
            const CBSStateSI& s,
            const Int ie)
        {
            return s.mat_elem(ie) == 0;
        }


        bool is_temperature_dirichlet_bc(
            const CBSStateSI& s,
            const Int bc)
        {
            return
                bc == s.cfg.bc_temperature_one_noslip ||
                bc == s.cfg.bc_temperature_zero_noslip ||
                bc == s.cfg.bc_temperature_zero_prescribed_velocity ||
                bc == s.cfg.bc_parabolic_inlet ||
                bc == s.cfg.bc_velocity_temperature_inlet ||
                bc == s.cfg.bc_massflow_temperature_inlet;
        }


        // Builds the exact elemental linear operator A_e used by the current
        // high-order Step-4 residual:
        //
        //     r_internal^(e) = A_e T^n.
        //
        // Volumetric and boundary sources are deliberately excluded because
        // they are present in both the low- and high-order schemes and therefore
        // generate no antidiffusive flux.
        void build_element_operator(
            const CBSStateSI& s,
            const Int ie,
            Real operator_matrix[5][5])
        {
            for (Int a = 1; a <= s.cfg.nep; ++a)
            {
                for (Int b = 1; b <= s.cfg.nep; ++b)
                {
                    operator_matrix[a][b] = 0.0;
                }
            }

            const Real volume = s.detJ(ie) / 6.0;

            if (volume <= 0.0 || !std::isfinite(volume))
            {
                throw std::runtime_error(
                    "TemperatureAFC encountered invalid element volume");
            }

            if (is_fluid_element(s, ie))
            {
                Real u_sum = 0.0;
                Real v_sum = 0.0;
                Real w_sum = 0.0;

                for (Int a = 1; a <= s.cfg.nep; ++a)
                {
                    const Int ip = s.intma(a, ie);
                    u_sum += s.unkno(1, ip);
                    v_sum += s.unkno(2, ip);
                    w_sum += s.unkno(3, ip);
                }

                const Real advection_factor =
                    s.rho_cp_e(ie) * s.detJ(ie) / 120.0;

                for (Int a = 1; a <= s.cfg.nep; ++a)
                {
                    const Int ip = s.intma(a, ie);

                    const Real u_weight = u_sum + s.unkno(1, ip);
                    const Real v_weight = v_sum + s.unkno(2, ip);
                    const Real w_weight = w_sum + s.unkno(3, ip);

                    for (Int b = 1; b <= s.cfg.nep; ++b)
                    {
                        operator_matrix[a][b] -=
                            advection_factor *
                            (u_weight * grad(s, ie, 1, b) +
                             v_weight * grad(s, ie, 2, b) +
                             w_weight * grad(s, ie, 3, b));
                    }
                }

                Real ubar = u_sum / static_cast<Real>(s.cfg.nep);
                Real vbar = v_sum / static_cast<Real>(s.cfg.nep);
                Real wbar = w_sum / static_cast<Real>(s.cfg.nep);

                const Real dt = s.delte(ie);

                if (dt <= 0.0 || !std::isfinite(dt))
                {
                    throw std::runtime_error(
                        "TemperatureAFC encountered invalid element timestep");
                }

                const Real stabilisation_factor =
                    0.5 * dt * s.rho_cp_e(ie) * volume;

                for (Int a = 1; a <= s.cfg.nep; ++a)
                {
                    const Real u_grad_a =
                        ubar * grad(s, ie, 1, a) +
                        vbar * grad(s, ie, 2, a) +
                        wbar * grad(s, ie, 3, a);

                    for (Int b = 1; b <= s.cfg.nep; ++b)
                    {
                        const Real u_grad_b =
                            ubar * grad(s, ie, 1, b) +
                            vbar * grad(s, ie, 2, b) +
                            wbar * grad(s, ie, 3, b);

                        operator_matrix[a][b] +=
                            stabilisation_factor * u_grad_a * u_grad_b;
                    }
                }
            }

            Real conductivity = s.k_e(ie);

            if (s.cfg.turbulence_on > 0 &&
                s.cfg.turbulent_thermal_diffusivity_on > 0 &&
                is_fluid_element(s, ie))
            {
                conductivity = s.k_eff_e(ie);
            }

            if (conductivity <= 0.0 || !std::isfinite(conductivity))
            {
                throw std::runtime_error(
                    "TemperatureAFC encountered invalid thermal conductivity");
            }

            for (Int a = 1; a <= s.cfg.nep; ++a)
            {
                for (Int b = 1; b <= s.cfg.nep; ++b)
                {
                    operator_matrix[a][b] -=
                        conductivity * volume *
                        (grad(s, ie, 1, a) * grad(s, ie, 1, b) +
                         grad(s, ie, 2, a) * grad(s, ie, 2, b) +
                         grad(s, ie, 3, a) * grad(s, ie, 3, b));
                }
            }
        }


        // Reverse min/max reductions for one scalar nodal array.  The existing
        // HaloExchange sum path cannot represent extrema, so AFC uses dedicated
        // tags and the same owner/ghost node maps.
        void reduce_extrema_to_owners(
            Array1D<Real>& values,
            const PartitionMetadata& metadata,
            const bool take_minimum,
            MPI_Comm communicator)
        {
            const int tag = take_minimum ? 930 : 931;

            std::vector<std::vector<Real>> send_buffers(
                metadata.neighbours.size());
            std::vector<std::vector<Real>> recv_buffers(
                metadata.neighbours.size());
            std::vector<MPI_Request> requests;
            requests.reserve(2U * metadata.neighbours.size());

            for (Size i = 0; i < metadata.neighbours.size(); ++i)
            {
                const PartitionNeighbour& neighbour = metadata.neighbours[i];
                recv_buffers[i].resize(neighbour.send_local_nodes.size());

                MPI_Request request{};
                check_mpi(
                    MPI_Irecv(
                        recv_buffers[i].data(),
                        static_cast<int>(recv_buffers[i].size()),
                        MPI_DOUBLE,
                        neighbour.rank,
                        tag,
                        communicator,
                        &request),
                    "MPI_Irecv extrema");
                requests.push_back(request);
            }

            for (Size i = 0; i < metadata.neighbours.size(); ++i)
            {
                const PartitionNeighbour& neighbour = metadata.neighbours[i];
                const std::vector<Int>& nodes = neighbour.recv_local_nodes;

                send_buffers[i].resize(nodes.size());

                for (Size j = 0; j < nodes.size(); ++j)
                {
                    send_buffers[i][j] = values(nodes[j]);
                }

                MPI_Request request{};
                check_mpi(
                    MPI_Isend(
                        send_buffers[i].data(),
                        static_cast<int>(send_buffers[i].size()),
                        MPI_DOUBLE,
                        neighbour.rank,
                        tag,
                        communicator,
                        &request),
                    "MPI_Isend extrema");
                requests.push_back(request);
            }

            if (!requests.empty())
            {
                check_mpi(
                    MPI_Waitall(
                        static_cast<int>(requests.size()),
                        requests.data(),
                        MPI_STATUSES_IGNORE),
                    "MPI_Waitall extrema");
            }

            for (Size i = 0; i < metadata.neighbours.size(); ++i)
            {
                const PartitionNeighbour& neighbour = metadata.neighbours[i];
                const std::vector<Int>& nodes = neighbour.send_local_nodes;

                for (Size j = 0; j < nodes.size(); ++j)
                {
                    if (take_minimum)
                    {
                        values(nodes[j]) =
                            std::min(values(nodes[j]), recv_buffers[i][j]);
                    }
                    else
                    {
                        values(nodes[j]) =
                            std::max(values(nodes[j]), recv_buffers[i][j]);
                    }
                }
            }
        }


        void build_temperature_stencil_bounds(
            const CBSStateSI& s,
            Array1D<Real>& lower,
            Array1D<Real>& upper)
        {
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                lower(ip) = s.temperature1(ip);
                upper(ip) = s.temperature1(ip);
            }

            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                Real element_minimum =
                    std::numeric_limits<Real>::max();
                Real element_maximum =
                    -std::numeric_limits<Real>::max();

                for (Int a = 1; a <= s.cfg.nep; ++a)
                {
                    const Real value =
                        s.temperature1(s.intma(a, ie));

                    element_minimum = std::min(element_minimum, value);
                    element_maximum = std::max(element_maximum, value);
                }

                for (Int a = 1; a <= s.cfg.nep; ++a)
                {
                    const Int ip = s.intma(a, ie);
                    lower(ip) = std::min(lower(ip), element_minimum);
                    upper(ip) = std::max(upper(ip), element_maximum);
                }
            }

            reduce_extrema_to_owners(
                lower,
                s.partition_metadata,
                true,
                MPI_COMM_WORLD);

            reduce_extrema_to_owners(
                upper,
                s.partition_metadata,
                false,
                MPI_COMM_WORLD);

            HaloExchange::broadcastOwnedToGhosts(
                lower,
                s.partition_metadata,
                MPI_COMM_WORLD);

            HaloExchange::broadcastOwnedToGhosts(
                upper,
                s.partition_metadata,
                MPI_COMM_WORLD);
        }


        Array1D<Int> build_fixed_temperature_mask(CBSStateSI& s)
        {
            Array1D<Int> fixed(s.cfg.npoin);
            fixed.fill(0);

            for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
            {
                const Int bc = s.iside(s.cfg.bsid, ib);

                if (!is_temperature_dirichlet_bc(s, bc))
                {
                    continue;
                }

                for (Int in = 1; in <= s.cfg.nsidp; ++in)
                {
                    fixed(s.iside(in, ib)) = 1;
                }
            }

            HaloExchange::orGhostMasksToOwners(
                fixed,
                s.partition_metadata,
                MPI_COMM_WORLD);

            HaloExchange::broadcastOwnedToGhosts(
                fixed,
                s.partition_metadata,
                MPI_COMM_WORLD);

            return fixed;
        }
    }
#endif


    bool TemperatureAFC::enabled()
    {
#if defined(CBS3D_USE_MPI)
        return environment_flag_enabled("CBS3D_TEMPERATURE_AFC");
#else
        return false;
#endif
    }


    TemperatureAFC::Diagnostics
    TemperatureAFC::limitDistributedResidual(CBSStateSI& s)
    {
#if !defined(CBS3D_USE_MPI)
        (void)s;
        throw std::runtime_error(
            "TemperatureAFC requires an MPI-enabled build");
#else
        if (!s.mpi_enabled)
        {
            throw std::runtime_error(
                "TemperatureAFC requires more than one MPI rank");
        }

        if (s.cfg.ndim != 3 || s.cfg.nep != 4)
        {
            throw std::runtime_error(
                "TemperatureAFC requires four-node 3D tetrahedra");
        }

        std::vector<EdgeFlux> edge_fluxes;
        edge_fluxes.reserve(
            static_cast<Size>(s.cfg.nelem) * 6U);

        // Add elemental graph viscosity to the high-order residual.  For each
        // pair a,b, d_ab=max(0,-A_ab,-A_ba) makes both low-order off-diagonal
        // coefficients non-negative.  The equal-and-opposite raw flux stored
        // below removes this added diffusion when alpha_ab=1.
        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            Real element_operator[5][5]{};
            build_element_operator(s, ie, element_operator);

            for (Int a = 1; a <= s.cfg.nep; ++a)
            {
                for (Int b = a + 1; b <= s.cfg.nep; ++b)
                {
                    const Real graph_viscosity =
                        std::max({
                            0.0,
                            -element_operator[a][b],
                            -element_operator[b][a]});

                    if (graph_viscosity <= 0.0)
                    {
                        continue;
                    }

                    const Int node_i = s.intma(a, ie);
                    const Int node_j = s.intma(b, ie);

                    const Real raw_flux =
                        graph_viscosity *
                        (s.temperature1(node_i) -
                         s.temperature1(node_j));

                    if (!std::isfinite(raw_flux))
                    {
                        throw std::runtime_error(
                            "TemperatureAFC generated a non-finite edge flux");
                    }

                    // Low-order residual = high-order residual + graph
                    // diffusion.  raw_flux is the antidiffusive flux that
                    // recovers the high-order residual.
                    s.rhs1(node_i) -= raw_flux;
                    s.rhs1(node_j) += raw_flux;

                    edge_fluxes.push_back(
                        EdgeFlux{node_i, node_j, raw_flux});
                }
            }
        }

        // Complete the distributed low-order residual on owners.
        HaloExchange::sumGhostContributionsToOwners(
            s.rhs1,
            s.partition_metadata,
            MPI_COMM_WORLD);

        Array1D<Real> stencil_minimum(s.cfg.npoin);
        Array1D<Real> stencil_maximum(s.cfg.npoin);
        build_temperature_stencil_bounds(
            s,
            stencil_minimum,
            stencil_maximum);

        Array1D<Int> fixed = build_fixed_temperature_mask(s);

        Array1D<Real> low_temperature(s.cfg.npoin);
        Array1D<Real> p_plus(s.cfg.npoin);
        Array1D<Real> p_minus(s.cfg.npoin);
        Array1D<Real> r_plus(s.cfg.npoin);
        Array1D<Real> r_minus(s.cfg.npoin);
        Array1D<Real> correction(s.cfg.npoin);

        low_temperature.fill(0.0);
        p_plus.fill(0.0);
        p_minus.fill(0.0);
        r_plus.fill(1.0);
        r_minus.fill(1.0);
        correction.fill(0.0);

        Real local_minimum_low =
            std::numeric_limits<Real>::max();

        for (const Int ip : s.owned_nodes)
        {
            if (s.elcoe2p(ip) <= 0.0 ||
                !std::isfinite(s.elcoe2p(ip)))
            {
                throw std::runtime_error(
                    "TemperatureAFC encountered invalid thermal diagonal");
            }

            low_temperature(ip) =
                s.temperature1(ip) +
                s.rhs1(ip) * s.elcoe2p(ip);

            if (!std::isfinite(low_temperature(ip)))
            {
                throw std::runtime_error(
                    "TemperatureAFC generated non-finite low-order temperature");
            }

            local_minimum_low =
                std::min(local_minimum_low, low_temperature(ip));
        }

        HaloExchange::broadcastOwnedToGhosts(
            low_temperature,
            s.partition_metadata,
            MPI_COMM_WORLD);

        for (const EdgeFlux& edge : edge_fluxes)
        {
            const Real flux = edge.flux_to_i;

            if (flux >= 0.0)
            {
                p_plus(edge.node_i) += flux;
                p_minus(edge.node_j) -= flux;
            }
            else
            {
                p_minus(edge.node_i) += flux;
                p_plus(edge.node_j) -= flux;
            }
        }

        HaloExchange::sumGhostContributionsToOwners(
            p_plus,
            s.partition_metadata,
            MPI_COMM_WORLD);

        HaloExchange::sumGhostContributionsToOwners(
            p_minus,
            s.partition_metadata,
            MPI_COMM_WORLD);

        Real physical_floor = 0.0;
        const bool use_physical_floor =
            environment_real(
                "CBS3D_AFC_MIN_TEMPERATURE",
                physical_floor);

        Real global_minimum_low = 0.0;
        check_mpi(
            MPI_Allreduce(
                &local_minimum_low,
                &global_minimum_low,
                1,
                MPI_DOUBLE,
                MPI_MIN,
                MPI_COMM_WORLD),
            "MPI_Allreduce AFC low-order Tmin");

        if (use_physical_floor)
        {
            const Real tolerance =
                1.0e-10 * std::max(1.0, std::abs(physical_floor));

            if (global_minimum_low < physical_floor - tolerance)
            {
                throw std::runtime_error(
                    "TemperatureAFC low-order stage violated "
                    "CBS3D_AFC_MIN_TEMPERATURE; reduce dt or inspect the "
                    "thermal operator");
            }
        }

        for (const Int ip : s.owned_nodes)
        {
            if (fixed(ip) != 0)
            {
                r_plus(ip) = 0.0;
                r_minus(ip) = 0.0;
                continue;
            }

            // Source-aware Zalesak bounds.  The low-order state is included so
            // a legitimate volumetric source is not clipped.  In the heating
            // blanket case the optional physical floor remains the strict lower
            // safety bound.
            Real lower_bound =
                std::min(stencil_minimum(ip), low_temperature(ip));
            const Real upper_bound =
                std::max(stencil_maximum(ip), low_temperature(ip));

            if (use_physical_floor)
            {
                lower_bound = std::max(lower_bound, physical_floor);
            }

            const Real q_plus =
                std::max(
                    0.0,
                    (upper_bound - low_temperature(ip)) /
                    s.elcoe2p(ip));

            const Real q_minus =
                std::min(
                    0.0,
                    (lower_bound - low_temperature(ip)) /
                    s.elcoe2p(ip));

            if (p_plus(ip) > 0.0)
            {
                r_plus(ip) =
                    std::min(1.0, q_plus / p_plus(ip));
            }

            if (p_minus(ip) < 0.0)
            {
                r_minus(ip) =
                    std::min(1.0, q_minus / p_minus(ip));
            }

            r_plus(ip) = std::clamp(r_plus(ip), 0.0, 1.0);
            r_minus(ip) = std::clamp(r_minus(ip), 0.0, 1.0);
        }

        HaloExchange::broadcastOwnedToGhosts(
            r_plus,
            s.partition_metadata,
            MPI_COMM_WORLD);

        HaloExchange::broadcastOwnedToGhosts(
            r_minus,
            s.partition_metadata,
            MPI_COMM_WORLD);

        long long local_limited_edges = 0;
        Real local_minimum_alpha = 1.0;

        for (const EdgeFlux& edge : edge_fluxes)
        {
            Real alpha = 1.0;

            if (edge.flux_to_i >= 0.0)
            {
                alpha =
                    std::min(
                        r_plus(edge.node_i),
                        r_minus(edge.node_j));
            }
            else
            {
                alpha =
                    std::min(
                        r_minus(edge.node_i),
                        r_plus(edge.node_j));
            }

            alpha = std::clamp(alpha, 0.0, 1.0);
            local_minimum_alpha =
                std::min(local_minimum_alpha, alpha);

            if (alpha < 1.0 - 1.0e-13)
            {
                ++local_limited_edges;
            }

            const Real limited_flux = alpha * edge.flux_to_i;
            correction(edge.node_i) += limited_flux;
            correction(edge.node_j) -= limited_flux;
        }

        HaloExchange::sumGhostContributionsToOwners(
            correction,
            s.partition_metadata,
            MPI_COMM_WORLD);

        Array1D<Real> limited_residual(s.cfg.npoin);
        limited_residual.fill(0.0);

        Real local_minimum_limited =
            std::numeric_limits<Real>::max();
        Real local_maximum_limited =
            -std::numeric_limits<Real>::max();
        Real local_correction_balance = 0.0;

        for (const Int ip : s.owned_nodes)
        {
            limited_residual(ip) = s.rhs1(ip) + correction(ip);

            const Real limited_temperature =
                s.temperature1(ip) +
                limited_residual(ip) * s.elcoe2p(ip);

            if (!std::isfinite(limited_temperature))
            {
                throw std::runtime_error(
                    "TemperatureAFC generated non-finite limited temperature");
            }

            local_minimum_limited =
                std::min(local_minimum_limited, limited_temperature);
            local_maximum_limited =
                std::max(local_maximum_limited, limited_temperature);
            local_correction_balance += correction(ip);
        }

        Real global_minimum_limited = 0.0;
        Real global_maximum_limited = 0.0;

        check_mpi(
            MPI_Allreduce(
                &local_minimum_limited,
                &global_minimum_limited,
                1,
                MPI_DOUBLE,
                MPI_MIN,
                MPI_COMM_WORLD),
            "MPI_Allreduce AFC limited Tmin");

        check_mpi(
            MPI_Allreduce(
                &local_maximum_limited,
                &global_maximum_limited,
                1,
                MPI_DOUBLE,
                MPI_MAX,
                MPI_COMM_WORLD),
            "MPI_Allreduce AFC limited Tmax");

        if (use_physical_floor)
        {
            const Real tolerance =
                1.0e-10 * std::max(1.0, std::abs(physical_floor));

            if (global_minimum_limited < physical_floor - tolerance)
            {
                throw std::runtime_error(
                    "TemperatureAFC limited state violated "
                    "CBS3D_AFC_MIN_TEMPERATURE");
            }
        }

        // The established production loop performs one more reverse-add after
        // this function returns.  Keep only globally assembled owner residuals
        // so that operation does not double-count AFC contributions.
        s.rhs1 = limited_residual;

        Diagnostics diagnostics;
        const long long local_edges =
            static_cast<long long>(edge_fluxes.size());

        check_mpi(
            MPI_Allreduce(
                &local_edges,
                &diagnostics.elemental_edges,
                1,
                MPI_LONG_LONG,
                MPI_SUM,
                MPI_COMM_WORLD),
            "MPI_Allreduce AFC edge count");

        check_mpi(
            MPI_Allreduce(
                &local_limited_edges,
                &diagnostics.limited_edges,
                1,
                MPI_LONG_LONG,
                MPI_SUM,
                MPI_COMM_WORLD),
            "MPI_Allreduce AFC limited-edge count");

        check_mpi(
            MPI_Allreduce(
                &local_minimum_alpha,
                &diagnostics.minimum_alpha,
                1,
                MPI_DOUBLE,
                MPI_MIN,
                MPI_COMM_WORLD),
            "MPI_Allreduce AFC minimum alpha");

        check_mpi(
            MPI_Allreduce(
                &local_correction_balance,
                &diagnostics.global_correction_balance,
                1,
                MPI_DOUBLE,
                MPI_SUM,
                MPI_COMM_WORLD),
            "MPI_Allreduce AFC correction balance");

        diagnostics.minimum_low_order_temperature =
            global_minimum_low;
        diagnostics.minimum_limited_temperature =
            global_minimum_limited;
        diagnostics.maximum_limited_temperature =
            global_maximum_limited;

        static long long call_count = 0;
        ++call_count;

        const Int log_every =
            environment_positive_int(
                "CBS3D_AFC_LOG_EVERY",
                1000);

        if (s.mpi_rank == 0 &&
            (call_count == 1 ||
             (call_count % log_every) == 0))
        {
            std::cout
                << "Temperature AFC"
                << "  call=" << call_count
                << "  edges=" << diagnostics.elemental_edges
                << "  limited=" << diagnostics.limited_edges
                << "  alpha_min=" << diagnostics.minimum_alpha
                << "  Tlow_min="
                << diagnostics.minimum_low_order_temperature
                << "  Tafc_min="
                << diagnostics.minimum_limited_temperature
                << "  Tafc_max="
                << diagnostics.maximum_limited_temperature
                << "  correction_balance="
                << diagnostics.global_correction_balance
                << '\n';
        }

        return diagnostics;
#endif
    }
}
