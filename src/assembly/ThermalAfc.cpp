//=============================================================================
// CBS3D++_SI
//
// Local-extremum-preserving algebraic flux correction for CBS Step 4.
//
// The existing EnergyAssembly residual is retained as the high-order target.
// For each tetrahedral node pair, symmetric algebraic diffusion is added using
//
//     d_ab = max(0, A_ab, A_ba).
//
// This makes every low-order off-diagonal coefficient non-positive. The raw
// antidiffusive pair flux that recovers the high-order residual is then limited
// with nodal Zalesak factors. Element-edge contributions are intentionally not
// deduplicated: their sum is conservative, monotone and naturally compatible
// with partitioned finite-element assembly.
//=============================================================================

#include "cbs/assembly/ThermalAfc.hpp"

#include "cbs/assembly/EnergyAssembly.hpp"
#include "cbs/parallel/HaloExchange.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace cbs
{
    namespace
    {
        constexpr Real numerical_tolerance = 1.0e-12;

        struct AfcWork
        {
            explicit AfcWork(const Int npoin)
                : low_rhs(npoin),
                  low_diagonal(npoin),
                  lower_bound(npoin),
                  upper_bound(npoin),
                  positive_flux(npoin),
                  negative_flux(npoin),
                  low_temperature(npoin),
                  positive_ratio(npoin),
                  negative_ratio(npoin),
                  correction(npoin)
            {
                low_rhs.fill(0.0);
                low_diagonal.fill(0.0);
                lower_bound.fill(0.0);
                upper_bound.fill(0.0);
                positive_flux.fill(0.0);
                negative_flux.fill(0.0);
                low_temperature.fill(0.0);
                positive_ratio.fill(1.0);
                negative_ratio.fill(1.0);
                correction.fill(0.0);
            }

            Array1D<Real> low_rhs;
            Array1D<Real> low_diagonal;
            Array1D<Real> lower_bound;
            Array1D<Real> upper_bound;
            Array1D<Real> positive_flux;
            Array1D<Real> negative_flux;
            Array1D<Real> low_temperature;
            Array1D<Real> positive_ratio;
            Array1D<Real> negative_ratio;
            Array1D<Real> correction;
        };

        constexpr std::array<std::array<Int, 2>, 6> local_edges =
        {{
            {{1, 2}},
            {{1, 3}},
            {{1, 4}},
            {{2, 3}},
            {{2, 4}},
            {{3, 4}}
        }};

        Real pair_diffusion(
            const Real matrix[5][5],
            const Int a,
            const Int b)
        {
            return std::max({0.0, matrix[a][b], matrix[b][a]});
        }

        void add_signed_flux(
            Array1D<Real>& positive,
            Array1D<Real>& negative,
            const Int node,
            const Real flux)
        {
            if (flux >= 0.0)
            {
                positive(node) += flux;
            }
            else
            {
                negative(node) += flux;
            }
        }

        void initialise_work(
            const CBSStateSI& s,
            AfcWork& work)
        {
            for (Int node = 1; node <= s.cfg.npoin; ++node)
            {
                const Real temperature = s.temperature1(node);

                if (!std::isfinite(temperature))
                {
                    throw std::runtime_error(
                        "ThermalAfc - non-finite old temperature at node "
                        + std::to_string(node));
                }

                work.low_rhs(node) = s.rhs1(node);
                work.lower_bound(node) = temperature;
                work.upper_bound(node) = temperature;
            }
        }

        void assemble_low_order_data(
            const CBSStateSI& s,
            AfcWork& work)
        {
            Real matrix[5][5]{};

            for (Int element = 1; element <= s.cfg.nelem; ++element)
            {
                EnergyAssembly::buildElementTransportMatrix(
                    s,
                    element,
                    matrix);

                for (Int a = 1; a <= s.cfg.nep; ++a)
                {
                    const Int node = s.intma(a, element);
                    work.low_diagonal(node) += matrix[a][a];
                }

                for (const auto& edge : local_edges)
                {
                    const Int a = edge[0];
                    const Int b = edge[1];
                    const Int node_a = s.intma(a, element);
                    const Int node_b = s.intma(b, element);
                    const Real diffusion = pair_diffusion(matrix, a, b);

                    const Real temperature_a = s.temperature1(node_a);
                    const Real temperature_b = s.temperature1(node_b);

                    work.lower_bound(node_a) =
                        std::min(work.lower_bound(node_a), temperature_b);
                    work.upper_bound(node_a) =
                        std::max(work.upper_bound(node_a), temperature_b);
                    work.lower_bound(node_b) =
                        std::min(work.lower_bound(node_b), temperature_a);
                    work.upper_bound(node_b) =
                        std::max(work.upper_bound(node_b), temperature_a);

                    if (diffusion <= 0.0)
                    {
                        continue;
                    }

                    work.low_diagonal(node_a) += diffusion;
                    work.low_diagonal(node_b) += diffusion;

                    // D*T pair contribution. The low-order residual is
                    //
                    //     r_L = r_H - D*T.
                    const Real antidiffusive_flux =
                        diffusion * (temperature_a - temperature_b);

                    work.low_rhs(node_a) -= antidiffusive_flux;
                    work.low_rhs(node_b) += antidiffusive_flux;

                    // The opposite pair is the raw correction that recovers the
                    // high-order residual from the low-order residual.
                    add_signed_flux(
                        work.positive_flux,
                        work.negative_flux,
                        node_a,
                        antidiffusive_flux);
                    add_signed_flux(
                        work.positive_flux,
                        work.negative_flux,
                        node_b,
                        -antidiffusive_flux);
                }
            }
        }

        Real bounded_ratio(const Real value)
        {
            if (!std::isfinite(value))
            {
                throw std::runtime_error(
                    "ThermalAfc - non-finite limiter ratio");
            }

            return std::clamp(value, 0.0, 1.0);
        }

        void calculate_limiter_ratios(
            const CBSStateSI& s,
            const Array1D<Int>& fixed,
            AfcWork& work,
            const std::vector<Int>* active_nodes)
        {
            const auto process_node = [&](const Int node)
            {
                if (fixed(node) != 0)
                {
                    work.positive_ratio(node) = 0.0;
                    work.negative_ratio(node) = 0.0;
                    return;
                }

                const Real inverse_diagonal = s.elcoe2p(node);
                if (inverse_diagonal <= 0.0
                    || !std::isfinite(inverse_diagonal))
                {
                    throw std::runtime_error(
                        "ThermalAfc - invalid inverse thermal diagonal at node "
                        + std::to_string(node));
                }

                // External heating/cooling is already present in the monotone
                // low-order state. Include that state in the admissible interval
                // so that the limiter constrains only the antidiffusive recovery.
                const Real lower = std::min(
                    work.lower_bound(node),
                    work.low_temperature(node));
                const Real upper = std::max(
                    work.upper_bound(node),
                    work.low_temperature(node));

                const Real positive_budget =
                    (upper - work.low_temperature(node)) / inverse_diagonal;
                const Real negative_budget =
                    (lower - work.low_temperature(node)) / inverse_diagonal;

                work.positive_ratio(node) =
                    work.positive_flux(node) > 0.0
                    ? bounded_ratio(
                        positive_budget / work.positive_flux(node))
                    : 1.0;

                work.negative_ratio(node) =
                    work.negative_flux(node) < 0.0
                    ? bounded_ratio(
                        negative_budget / work.negative_flux(node))
                    : 1.0;
            };

            if (active_nodes != nullptr)
            {
                for (const Int node : *active_nodes)
                {
                    process_node(node);
                }
                return;
            }

            for (Int node = 1; node <= s.cfg.npoin; ++node)
            {
                process_node(node);
            }
        }

        void assemble_limited_correction(
            const CBSStateSI& s,
            AfcWork& work)
        {
            work.correction.fill(0.0);
            Real matrix[5][5]{};

            for (Int element = 1; element <= s.cfg.nelem; ++element)
            {
                EnergyAssembly::buildElementTransportMatrix(
                    s,
                    element,
                    matrix);

                for (const auto& edge : local_edges)
                {
                    const Int a = edge[0];
                    const Int b = edge[1];
                    const Real diffusion = pair_diffusion(matrix, a, b);

                    if (diffusion <= 0.0)
                    {
                        continue;
                    }

                    const Int node_a = s.intma(a, element);
                    const Int node_b = s.intma(b, element);
                    const Real raw_flux = diffusion *
                        (s.temperature1(node_a) - s.temperature1(node_b));

                    const Real alpha = raw_flux >= 0.0
                        ? std::min(
                            work.positive_ratio(node_a),
                            work.negative_ratio(node_b))
                        : std::min(
                            work.negative_ratio(node_a),
                            work.positive_ratio(node_b));

                    const Real limited_flux = alpha * raw_flux;
                    work.correction(node_a) += limited_flux;
                    work.correction(node_b) -= limited_flux;
                }
            }
        }

        Real temperature_dirichlet_value(
            const CBSStateSI& s,
            const Int bc,
            bool& prescribed)
        {
            prescribed = true;

            if (bc == s.cfg.bc_temperature_one_noslip)
            {
                return 1.0;
            }

            if (bc == s.cfg.bc_temperature_zero_noslip
                || bc == s.cfg.bc_temperature_zero_prescribed_velocity
                || bc == s.cfg.bc_parabolic_inlet)
            {
                return 0.0;
            }

            if (bc == s.cfg.bc_velocity_temperature_inlet
                || bc == s.cfg.bc_massflow_temperature_inlet)
            {
                return s.cfg.inlet_temperature;
            }

            prescribed = false;
            return 0.0;
        }

        void build_serial_fixed_temperature(
            const CBSStateSI& s,
            Array1D<Int>& fixed,
            Array1D<Real>& value)
        {
            Array1D<Real> count(s.cfg.npoin);
            Array1D<Real> sum(s.cfg.npoin);
            Array1D<Real> sum_squared(s.cfg.npoin);

            fixed.fill(0);
            value.fill(0.0);
            count.fill(0.0);
            sum.fill(0.0);
            sum_squared.fill(0.0);

            for (Int face = 1; face <= s.cfg.nboun; ++face)
            {
                bool prescribed = false;
                const Real boundary_value = temperature_dirichlet_value(
                    s,
                    s.iside(s.cfg.bsid, face),
                    prescribed);

                if (!prescribed)
                {
                    continue;
                }

                for (Int local_node = 1;
                     local_node <= s.cfg.nsidp;
                     ++local_node)
                {
                    const Int node = s.iside(local_node, face);
                    count(node) += 1.0;
                    sum(node) += boundary_value;
                    sum_squared(node) += boundary_value * boundary_value;
                }
            }

            for (Int node = 1; node <= s.cfg.npoin; ++node)
            {
                if (count(node) <= 0.0)
                {
                    continue;
                }

                const Real mean = sum(node) / count(node);
                const Real variance = std::max(
                    0.0,
                    sum_squared(node) / count(node) - mean * mean);

                if (variance > numerical_tolerance
                    * std::max(1.0, mean * mean))
                {
                    throw std::runtime_error(
                        "ThermalAfc - conflicting prescribed temperatures at "
                        "node " + std::to_string(node));
                }

                fixed(node) = 1;
                value(node) = mean;
            }
        }

        void check_serial_positivity(
            const CBSStateSI& s,
            const Array1D<Int>& fixed,
            const AfcWork& work)
        {
            Real maximum_ratio = 0.0;

            for (Int node = 1; node <= s.cfg.npoin; ++node)
            {
                if (fixed(node) != 0)
                {
                    continue;
                }

                maximum_ratio = std::max(
                    maximum_ratio,
                    s.elcoe2p(node) * work.low_diagonal(node));
            }

            if (maximum_ratio > 1.0 + numerical_tolerance)
            {
                throw std::runtime_error(
                    "ThermalAfc - Step-4 timestep violates the explicit "
                    "low-order positivity limit; maximum dt*Aii/M = "
                    + std::to_string(maximum_ratio));
            }
        }

#ifdef CBS3D_USE_MPI
        void check_mpi(const int error_code, const char* operation)
        {
            if (error_code != MPI_SUCCESS)
            {
                throw std::runtime_error(
                    std::string("ThermalAfc MPI failure in ") + operation);
            }
        }
#endif
    }

    void ThermalAfc::applySerial(CBSStateSI& s)
    {
        AfcWork work(s.cfg.npoin);
        Array1D<Int> fixed(s.cfg.npoin);
        Array1D<Real> fixed_value(s.cfg.npoin);

        build_serial_fixed_temperature(s, fixed, fixed_value);
        initialise_work(s, work);
        assemble_low_order_data(s, work);
        check_serial_positivity(s, fixed, work);

        for (Int node = 1; node <= s.cfg.npoin; ++node)
        {
            work.low_temperature(node) =
                s.temperature1(node)
                + s.elcoe2p(node) * work.low_rhs(node);

            if (fixed(node) != 0)
            {
                work.low_temperature(node) = fixed_value(node);
            }
        }

        calculate_limiter_ratios(s, fixed, work, nullptr);
        assemble_limited_correction(s, work);

        for (Int node = 1; node <= s.cfg.npoin; ++node)
        {
            s.temperature(node) =
                work.low_temperature(node)
                + s.elcoe2p(node) * work.correction(node);

            if (fixed(node) != 0)
            {
                s.temperature(node) = fixed_value(node);
            }

            if (!std::isfinite(s.temperature(node)))
            {
                throw std::runtime_error(
                    "ThermalAfc - non-finite corrected temperature at node "
                    + std::to_string(node));
            }

            s.rhs1(node) =
                (s.temperature(node) - s.temperature1(node))
                / s.elcoe2p(node);
        }
    }

#ifdef CBS3D_USE_MPI
    void ThermalAfc::applyDistributed(
        CBSStateSI& s,
        const Array1D<Int>& fixed_temperature,
        const Array1D<Real>& fixed_value,
        const Array1D<Real>& external_heat_flux_load,
        MPI_Comm communicator)
    {
        AfcWork work(s.cfg.npoin);

        initialise_work(s, work);
        assemble_low_order_data(s, work);

        HaloExchange::sumGhostContributionsToOwners(
            work.low_rhs,
            s.partition_metadata,
            communicator);
        HaloExchange::sumGhostContributionsToOwners(
            work.low_diagonal,
            s.partition_metadata,
            communicator);
        HaloExchange::sumGhostContributionsToOwners(
            work.positive_flux,
            s.partition_metadata,
            communicator);
        HaloExchange::sumGhostContributionsToOwners(
            work.negative_flux,
            s.partition_metadata,
            communicator);
        HaloExchange::minGhostContributionsToOwners(
            work.lower_bound,
            s.partition_metadata,
            communicator);
        HaloExchange::maxGhostContributionsToOwners(
            work.upper_bound,
            s.partition_metadata,
            communicator);

        Real local_maximum_ratio = 0.0;

        for (const Int node : s.owned_nodes)
        {
            work.low_rhs(node) += external_heat_flux_load(node);

            if (fixed_temperature(node) == 0)
            {
                local_maximum_ratio = std::max(
                    local_maximum_ratio,
                    s.elcoe2p(node) * work.low_diagonal(node));
            }

            work.low_temperature(node) =
                s.temperature1(node)
                + s.elcoe2p(node) * work.low_rhs(node);

            if (fixed_temperature(node) != 0)
            {
                work.low_temperature(node) = fixed_value(node);
            }
        }

        Real global_maximum_ratio = 0.0;
        check_mpi(
            MPI_Allreduce(
                &local_maximum_ratio,
                &global_maximum_ratio,
                1,
                MPI_DOUBLE,
                MPI_MAX,
                communicator),
            "MPI_Allreduce positivity ratio");

        if (global_maximum_ratio > 1.0 + numerical_tolerance)
        {
            throw std::runtime_error(
                "ThermalAfc - distributed Step-4 timestep violates the "
                "explicit low-order positivity limit; maximum dt*Aii/M = "
                + std::to_string(global_maximum_ratio));
        }

        HaloExchange::broadcastOwnedToGhosts(
            work.low_temperature,
            s.partition_metadata,
            communicator);
        HaloExchange::broadcastOwnedToGhosts(
            work.lower_bound,
            s.partition_metadata,
            communicator);
        HaloExchange::broadcastOwnedToGhosts(
            work.upper_bound,
            s.partition_metadata,
            communicator);

        calculate_limiter_ratios(
            s,
            fixed_temperature,
            work,
            &s.owned_nodes);

        HaloExchange::broadcastOwnedToGhosts(
            work.positive_ratio,
            s.partition_metadata,
            communicator);
        HaloExchange::broadcastOwnedToGhosts(
            work.negative_ratio,
            s.partition_metadata,
            communicator);

        assemble_limited_correction(s, work);

        HaloExchange::sumGhostContributionsToOwners(
            work.correction,
            s.partition_metadata,
            communicator);

        for (const Int node : s.owned_nodes)
        {
            s.temperature(node) =
                work.low_temperature(node)
                + s.elcoe2p(node) * work.correction(node);

            if (fixed_temperature(node) != 0)
            {
                s.temperature(node) = fixed_value(node);
            }

            if (!std::isfinite(s.temperature(node)))
            {
                throw std::runtime_error(
                    "ThermalAfc - distributed update produced non-finite "
                    "temperature at node " + std::to_string(node));
            }

            s.rhs1(node) =
                (s.temperature(node) - s.temperature1(node))
                / s.elcoe2p(node);
        }

        HaloExchange::broadcastOwnedToGhosts(
            s.temperature,
            s.partition_metadata,
            communicator);
    }
#endif
}
