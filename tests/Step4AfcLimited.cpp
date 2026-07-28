//=============================================================================
// CBS3D++_SI
//
// Limited algebraic-flux-correction validation for Step 4.
//
// The production high-order thermal matrix A_H is recovered by repeatedly
// calling EnergyAssembly::assembleStep4Rhs(). A monotone low-order matrix is
// formed as
//
//     A_L = A_H + D,
//
// where D is symmetric, conservative algebraic diffusion. The raw pairwise
// antidiffusive flux that recovers the high-order explicit update is
//
//     f_ij = d_ij (T_i^n - T_j^n),   f_ji = -f_ij.
//
// A Zalesak-style nodal limiter constrains these pairwise fluxes using local
// old-time neighbour bounds. The corrected update is therefore conservative,
// recovers part of the high-order transport, and introduces no new extrema.
//
// This executable is a validation candidate only. It does not modify the
// production EnergyAssembly.cpp implementation.
//=============================================================================

#define main cbs3d_step4_bounded_transport_baseline_main
#include "Step4BoundedTransport.cpp"
#undef main

#include <utility>

namespace
{
    struct FctSparseEntry
    {
        Int column = 0;
        Real value = 0.0;
    };

    struct FctEdge
    {
        Int i = 0;
        Int j = 0;
        Real diffusion = 0.0;
    };

    using FctSparseRows = std::vector<std::vector<FctSparseEntry>>;

    std::size_t fct_matrix_index(
        const Int nodes,
        const Int row,
        const Int column)
    {
        const std::size_t stride = static_cast<std::size_t>(nodes + 1);
        return static_cast<std::size_t>(row) * stride
             + static_cast<std::size_t>(column);
    }

    std::vector<Real> recover_high_order_matrix(
        const Mesh& mesh,
        const Real dt)
    {
        CBSStateSI state = make_state(mesh, dt);
        const Int nodes = state.cfg.npoin;

        std::vector<Real> matrix(
            static_cast<std::size_t>(nodes + 1)
          * static_cast<std::size_t>(nodes + 1),
            0.0);

        // For this prescribed-velocity, constant-property test:
        //
        //     rhs(T) = -A_H T.
        //
        // Applying the unmodified production residual to every nodal basis
        // vector recovers A_H without duplicating its element formula.
        for (Int column = 1; column <= nodes; ++column)
        {
            state.temperature.fill(0.0);
            state.temperature1.fill(0.0);
            state.temperature(column) = 1.0;
            state.temperature1(column) = 1.0;

            cbs::EnergyAssembly::assembleStep4Rhs(state);

            for (Int row = 1; row <= nodes; ++row)
            {
                matrix[fct_matrix_index(nodes, row, column)] =
                    -state.rhs1(row);
            }
        }

        return matrix;
    }

    std::vector<Real> build_low_order_matrix(
        const std::vector<Real>& high_order,
        const Int nodes)
    {
        std::vector<Real> low_order = high_order;

        for (Int i = 1; i <= nodes; ++i)
        {
            for (Int j = i + 1; j <= nodes; ++j)
            {
                const Real aij =
                    high_order[fct_matrix_index(nodes, i, j)];
                const Real aji =
                    high_order[fct_matrix_index(nodes, j, i)];
                const Real diffusion = std::max({0.0, aij, aji});

                if (diffusion <= 0.0)
                {
                    continue;
                }

                low_order[fct_matrix_index(nodes, i, j)] -= diffusion;
                low_order[fct_matrix_index(nodes, j, i)] -= diffusion;
                low_order[fct_matrix_index(nodes, i, i)] += diffusion;
                low_order[fct_matrix_index(nodes, j, j)] += diffusion;
            }
        }

        return low_order;
    }

    FctSparseRows make_sparse_rows(
        const std::vector<Real>& matrix,
        const Int nodes)
    {
        FctSparseRows rows(static_cast<std::size_t>(nodes + 1));
        constexpr Real storage_tolerance = 1.0e-24;

        for (Int i = 1; i <= nodes; ++i)
        {
            auto& row = rows[static_cast<std::size_t>(i)];
            for (Int j = 1; j <= nodes; ++j)
            {
                const Real value = matrix[fct_matrix_index(nodes, i, j)];
                if (std::abs(value) > storage_tolerance)
                {
                    row.push_back({j, value});
                }
            }
        }

        return rows;
    }

    std::vector<FctEdge> make_antidiffusive_edges(
        const std::vector<Real>& high_order,
        const std::vector<Real>& low_order,
        const Int nodes)
    {
        std::vector<FctEdge> edges;
        constexpr Real storage_tolerance = 1.0e-24;

        for (Int i = 1; i <= nodes; ++i)
        {
            for (Int j = i + 1; j <= nodes; ++j)
            {
                const Real diffusion_ij =
                    high_order[fct_matrix_index(nodes, i, j)]
                  - low_order[fct_matrix_index(nodes, i, j)];
                const Real diffusion_ji =
                    high_order[fct_matrix_index(nodes, j, i)]
                  - low_order[fct_matrix_index(nodes, j, i)];
                const Real diffusion =
                    0.5 * (diffusion_ij + diffusion_ji);

                if (diffusion > storage_tolerance)
                {
                    edges.push_back({i, j, diffusion});
                }
            }
        }

        return edges;
    }

    Real maximum_absolute_row_sum_fct(
        const std::vector<Real>& matrix,
        const Int nodes)
    {
        Real maximum = 0.0;
        for (Int i = 1; i <= nodes; ++i)
        {
            Real row_sum = 0.0;
            for (Int j = 1; j <= nodes; ++j)
            {
                row_sum += matrix[fct_matrix_index(nodes, i, j)];
            }
            maximum = std::max(maximum, std::abs(row_sum));
        }
        return maximum;
    }

    Real maximum_positive_off_diagonal_fct(
        const std::vector<Real>& matrix,
        const Int nodes)
    {
        Real maximum = 0.0;
        for (Int i = 1; i <= nodes; ++i)
        {
            for (Int j = 1; j <= nodes; ++j)
            {
                if (i != j)
                {
                    maximum = std::max(
                        maximum,
                        matrix[fct_matrix_index(nodes, i, j)]);
                }
            }
        }
        return maximum;
    }

    Real positivity_limit_fct(
        const std::vector<Real>& low_order,
        const std::vector<Real>& capacity,
        const Int nodes)
    {
        Real limit = std::numeric_limits<Real>::max();

        for (Int i = 1; i <= nodes; ++i)
        {
            const Real diagonal =
                low_order[fct_matrix_index(nodes, i, i)];
            if (diagonal > 0.0)
            {
                limit = std::min(
                    limit,
                    capacity[static_cast<std::size_t>(i)] / diagonal);
            }
        }

        if (!std::isfinite(limit)
            || limit == std::numeric_limits<Real>::max())
        {
            throw std::runtime_error(
                "Step4AfcLimited - unable to determine positivity timestep");
        }

        return limit;
    }

    void apply_low_order_residual(
        const FctSparseRows& rows,
        const std::vector<Real>& temperature,
        std::vector<Real>& residual,
        const Int nodes)
    {
        for (Int i = 1; i <= nodes; ++i)
        {
            Real matrix_temperature = 0.0;
            for (const FctSparseEntry& entry :
                 rows[static_cast<std::size_t>(i)])
            {
                matrix_temperature +=
                    entry.value
                    * temperature[static_cast<std::size_t>(entry.column)];
            }
            residual[static_cast<std::size_t>(i)] = -matrix_temperature;
        }
    }

    void enforce_inlet_vector(
        const CBSStateSI& state,
        std::vector<Real>& temperature)
    {
        for (Int node = 1; node <= state.cfg.npoin; ++node)
        {
            if (std::abs(state.coord(1, node)) <= 1.0e-14)
            {
                temperature[static_cast<std::size_t>(node)] =
                    reference_temperature;
            }
        }
    }

    struct LimiterTotals
    {
        Real raw_absolute_flux = 0.0;
        Real accepted_absolute_flux = 0.0;
        Real minimum_alpha = 1.0;
        std::uint64_t edge_visits = 0;
        std::uint64_t limited_edge_visits = 0;
    };

    bool run_limited_afc_case(const Mesh& mesh, const Real dt)
    {
        CBSStateSI state = make_state(mesh, dt);
        const std::vector<Real> capacity =
            assemble_lumped_capacity(state, dt);
        const Int nodes = state.cfg.npoin;

        const std::vector<Real> high_order =
            recover_high_order_matrix(mesh, dt);
        const std::vector<Real> low_order =
            build_low_order_matrix(high_order, nodes);
        const FctSparseRows sparse_low =
            make_sparse_rows(low_order, nodes);
        const std::vector<FctEdge> edges =
            make_antidiffusive_edges(high_order, low_order, nodes);

        const Real high_row_sum =
            maximum_absolute_row_sum_fct(high_order, nodes);
        const Real low_row_sum =
            maximum_absolute_row_sum_fct(low_order, nodes);
        const Real maximum_positive_offdiag =
            maximum_positive_off_diagonal_fct(low_order, nodes);
        const Real positivity_dt =
            positivity_limit_fct(low_order, capacity, nodes);

        const Int steps = static_cast<Int>(std::llround(final_time / dt));
        const Real initial_energy = excess_energy(state, capacity);

        std::vector<Real> old_temperature(
            static_cast<std::size_t>(nodes + 1), reference_temperature);
        std::vector<Real> low_temperature(
            static_cast<std::size_t>(nodes + 1), reference_temperature);
        std::vector<Real> residual(
            static_cast<std::size_t>(nodes + 1), 0.0);
        std::vector<Real> raw_flux(edges.size(), 0.0);
        std::vector<Real> positive_sum(
            static_cast<std::size_t>(nodes + 1), 0.0);
        std::vector<Real> negative_sum(
            static_cast<std::size_t>(nodes + 1), 0.0);
        std::vector<Real> upper_bound(
            static_cast<std::size_t>(nodes + 1), reference_temperature);
        std::vector<Real> lower_bound(
            static_cast<std::size_t>(nodes + 1), reference_temperature);
        std::vector<Real> positive_ratio(
            static_cast<std::size_t>(nodes + 1), 1.0);
        std::vector<Real> negative_ratio(
            static_cast<std::size_t>(nodes + 1), 1.0);
        std::vector<Real> correction(
            static_cast<std::size_t>(nodes + 1), 0.0);

        LimiterTotals limiter_totals;
        Int first_violation_step = -1;

        for (Int step = 1; step <= steps; ++step)
        {
            for (Int node = 1; node <= nodes; ++node)
            {
                const Real temperature = state.temperature(node);
                old_temperature[static_cast<std::size_t>(node)] = temperature;
                state.temperature1(node) = temperature;
                lower_bound[static_cast<std::size_t>(node)] = temperature;
                upper_bound[static_cast<std::size_t>(node)] = temperature;
                positive_sum[static_cast<std::size_t>(node)] = 0.0;
                negative_sum[static_cast<std::size_t>(node)] = 0.0;
                correction[static_cast<std::size_t>(node)] = 0.0;
            }

            // Local old-time neighbour bounds used by the nodal limiter.
            for (const FctEdge& edge : edges)
            {
                const Real ti =
                    old_temperature[static_cast<std::size_t>(edge.i)];
                const Real tj =
                    old_temperature[static_cast<std::size_t>(edge.j)];

                lower_bound[static_cast<std::size_t>(edge.i)] = std::min(
                    lower_bound[static_cast<std::size_t>(edge.i)], tj);
                upper_bound[static_cast<std::size_t>(edge.i)] = std::max(
                    upper_bound[static_cast<std::size_t>(edge.i)], tj);
                lower_bound[static_cast<std::size_t>(edge.j)] = std::min(
                    lower_bound[static_cast<std::size_t>(edge.j)], ti);
                upper_bound[static_cast<std::size_t>(edge.j)] = std::max(
                    upper_bound[static_cast<std::size_t>(edge.j)], ti);
            }

            apply_low_order_residual(
                sparse_low, old_temperature, residual, nodes);

            for (Int node = 1; node <= nodes; ++node)
            {
                low_temperature[static_cast<std::size_t>(node)] =
                    old_temperature[static_cast<std::size_t>(node)]
                    + dt * residual[static_cast<std::size_t>(node)]
                    / capacity[static_cast<std::size_t>(node)];
            }
            enforce_inlet_vector(state, low_temperature);

            // Raw conservative antidiffusive fluxes that recover A_H from A_L.
            for (std::size_t edge_index = 0;
                 edge_index < edges.size();
                 ++edge_index)
            {
                const FctEdge& edge = edges[edge_index];
                const Real flux = edge.diffusion
                    * (old_temperature[static_cast<std::size_t>(edge.i)]
                     - old_temperature[static_cast<std::size_t>(edge.j)]);
                raw_flux[edge_index] = flux;

                if (flux >= 0.0)
                {
                    positive_sum[static_cast<std::size_t>(edge.i)] += flux;
                    negative_sum[static_cast<std::size_t>(edge.j)] -= flux;
                }
                else
                {
                    negative_sum[static_cast<std::size_t>(edge.i)] += flux;
                    positive_sum[static_cast<std::size_t>(edge.j)] -= flux;
                }
            }

            for (Int node = 1; node <= nodes; ++node)
            {
                const std::size_t index = static_cast<std::size_t>(node);
                const Real allowable_positive =
                    capacity[index]
                    * (upper_bound[index] - low_temperature[index]) / dt;
                const Real allowable_negative =
                    capacity[index]
                    * (lower_bound[index] - low_temperature[index]) / dt;

                positive_ratio[index] =
                    (positive_sum[index] > 0.0)
                    ? std::min(
                        1.0,
                        std::max(0.0, allowable_positive)
                        / positive_sum[index])
                    : 1.0;

                negative_ratio[index] =
                    (negative_sum[index] < 0.0)
                    ? std::min(
                        1.0,
                        std::max(0.0, -allowable_negative)
                        / (-negative_sum[index]))
                    : 1.0;
            }

            for (std::size_t edge_index = 0;
                 edge_index < edges.size();
                 ++edge_index)
            {
                const FctEdge& edge = edges[edge_index];
                const Real flux = raw_flux[edge_index];
                Real alpha = 1.0;

                if (flux > 0.0)
                {
                    alpha = std::min(
                        positive_ratio[static_cast<std::size_t>(edge.i)],
                        negative_ratio[static_cast<std::size_t>(edge.j)]);
                }
                else if (flux < 0.0)
                {
                    alpha = std::min(
                        negative_ratio[static_cast<std::size_t>(edge.i)],
                        positive_ratio[static_cast<std::size_t>(edge.j)]);
                }

                const Real limited_flux = alpha * flux;
                correction[static_cast<std::size_t>(edge.i)] += limited_flux;
                correction[static_cast<std::size_t>(edge.j)] -= limited_flux;

                limiter_totals.raw_absolute_flux += std::abs(flux);
                limiter_totals.accepted_absolute_flux +=
                    std::abs(limited_flux);
                limiter_totals.minimum_alpha = std::min(
                    limiter_totals.minimum_alpha, alpha);
                ++limiter_totals.edge_visits;
                if (alpha < 1.0 - 1.0e-12 && std::abs(flux) > 1.0e-24)
                {
                    ++limiter_totals.limited_edge_visits;
                }
            }

            for (Int node = 1; node <= nodes; ++node)
            {
                const std::size_t index = static_cast<std::size_t>(node);
                state.temperature(node) =
                    low_temperature[index]
                    + dt * correction[index] / capacity[index];
            }
            enforce_inlet_temperature(state);

            if (first_violation_step < 0)
            {
                const Metrics current = measure(state);
                if (current.below_count > 0 || current.above_count > 0)
                {
                    first_violation_step = step;
                }
            }
        }

        const Metrics final_metrics = measure(state);
        const Real final_energy = excess_energy(state, capacity);
        const Real relative_energy_drift =
            (initial_energy != 0.0)
            ? (final_energy - initial_energy) / initial_energy
            : 0.0;
        const Real accepted_fraction =
            (limiter_totals.raw_absolute_flux > 0.0)
            ? limiter_totals.accepted_absolute_flux
              / limiter_totals.raw_absolute_flux
            : 0.0;
        const Real limited_edge_fraction =
            (limiter_totals.edge_visits > 0)
            ? static_cast<Real>(limiter_totals.limited_edge_visits)
              / static_cast<Real>(limiter_totals.edge_visits)
            : 0.0;

        const bool low_matrix_is_monotone =
            maximum_positive_offdiag <= 1.0e-12;
        const bool timestep_is_admissible =
            dt <= positivity_dt * (1.0 + 1.0e-12);
        const bool bounds_preserved =
            final_metrics.below_count == 0
            && final_metrics.above_count == 0
            && first_violation_step < 0;
        const bool antidiffusion_recovered =
            limiter_totals.accepted_absolute_flux > 0.0
            && accepted_fraction > 1.0e-6;
        const bool energy_conserved =
            std::abs(relative_energy_drift) <= 1.0e-4;

        std::cout << "\nCASE dt = " << std::scientific << dt << " s\n";
        std::cout << "  steps                         : " << steps << "\n";
        std::cout << "  production-matrix row sum    : "
                  << high_row_sum << "\n";
        std::cout << "  low-order-matrix row sum     : "
                  << low_row_sum << "\n";
        std::cout << "  maximum positive offdiag     : "
                  << maximum_positive_offdiag << "\n";
        std::cout << "  antidiffusive graph edges    : "
                  << edges.size() << "\n";
        std::cout << "  explicit positivity dt limit : "
                  << positivity_dt << " s\n";
        std::cout << "  dt / positivity limit        : "
                  << dt / positivity_dt << "\n";
        std::cout << "  accepted antidiffusion       : "
                  << accepted_fraction << "\n";
        std::cout << "  limited edge-visit fraction  : "
                  << limited_edge_fraction << "\n";
        std::cout << "  minimum limiter alpha        : "
                  << limiter_totals.minimum_alpha << "\n";
        std::cout << "  Tmin                          : "
                  << std::setprecision(12) << final_metrics.minimum << " K\n";
        std::cout << "  Tmax                          : "
                  << final_metrics.maximum << " K\n";
        std::cout << "  nodes below 300 K             : "
                  << final_metrics.below_count << "\n";
        std::cout << "  nodes above 400 K             : "
                  << final_metrics.above_count << "\n";
        std::cout << "  first violation step          : "
                  << first_violation_step << "\n";
        std::cout << "  relative excess-energy drift  : "
                  << relative_energy_drift << "\n";
        std::cout << "  LOW-ORDER M-MATRIX CHECK      : "
                  << (low_matrix_is_monotone ? "PASS" : "FAIL") << "\n";
        std::cout << "  TIMESTEP CHECK                : "
                  << (timestep_is_admissible ? "PASS" : "FAIL") << "\n";
        std::cout << "  ANTIDIFFUSION RECOVERY CHECK  : "
                  << (antidiffusion_recovered ? "PASS" : "FAIL") << "\n";
        std::cout << "  ENERGY CHECK                  : "
                  << (energy_conserved ? "PASS" : "FAIL") << "\n";
        std::cout << "  BOUND CHECK                   : "
                  << (bounds_preserved ? "PASS" : "FAIL") << "\n";

        return low_matrix_is_monotone
            && timestep_is_admissible
            && antidiffusion_recovered
            && energy_conserved
            && bounds_preserved;
    }
}

int main()
{
    const Mesh mesh = make_mesh();

    std::cout
        << "CBS3D++_SI Step-4 limited AFC validation\n"
        << "Production operator recovered through EnergyAssembly::assembleStep4Rhs\n"
        << "Limiter          : local-extremum-preserving pairwise Zalesak AFC\n"
        << "Mesh nodes      : " << mesh.coordinates.size() << "\n"
        << "Mesh tetrahedra : " << mesh.tetrahedra.size() << "\n"
        << "Initial bounds  : [300,400] K\n"
        << "Prescribed U    : (0.05,0,0) m/s\n";

    int failures = 0;
    failures += run_limited_afc_case(mesh, 1.0e-4) ? 0 : 1;
    failures += run_limited_afc_case(mesh, 1.0e-5) ? 0 : 1;

    std::cout << "\n============================================================\n";
    if (failures == 0)
    {
        std::cout
            << "LIMITED AFC VALIDATION: ALL CASES PASSED\n"
            << "The bounded low-order backbone now recovers limited high-order\n"
            << "antidiffusive flux without creating new temperature extrema.\n";
        return 0;
    }

    std::cout << "LIMITED AFC VALIDATION: "
              << failures << " CASE(S) FAILED\n";
    return 1;
}
