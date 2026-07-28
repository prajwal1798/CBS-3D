//=============================================================================
// CBS3D++_SI
//
// Algebraic low-order boundedness validation for Step 4.
//
// The current production thermal operator is first recovered as a global
// matrix by repeatedly calling EnergyAssembly::assembleStep4Rhs().  No
// duplicate finite-element formula is used.  A symmetric algebraic diffusion
// operator is then added so that every off-diagonal entry of the resulting
// transport matrix is non-positive:
//
//     d_ij = -max(0, a_ij, a_ji),  i != j,
//     d_ii = -sum_{j!=i} d_ij,
//     A_L  = A_H + D.
//
// With lumped thermal capacity M_L and an explicit step satisfying
//
//     dt <= min_i M_i / (A_L)_ii,
//
// the update is a convex combination of neighbouring nodal temperatures.
// This is the monotone low-order backbone required by AFC/FCT.  It is
// intentionally diffusive and is not proposed as the final production scheme.
//=============================================================================

#define main cbs3d_step4_bounded_transport_baseline_main
#include "Step4BoundedTransport.cpp"
#undef main

#include <utility>

namespace
{
    struct SparseEntry
    {
        Int column = 0;
        Real value = 0.0;
    };

    using SparseRows = std::vector<std::vector<SparseEntry>>;

    std::size_t matrix_index(
        const Int nodes,
        const Int row,
        const Int column)
    {
        const std::size_t stride = static_cast<std::size_t>(nodes + 1);
        return static_cast<std::size_t>(row) * stride
             + static_cast<std::size_t>(column);
    }

    std::vector<Real> recover_production_matrix(
        const Mesh& mesh,
        const Real dt)
    {
        CBSStateSI state = make_state(mesh, dt);
        const Int nodes = state.cfg.npoin;

        std::vector<Real> matrix(
            static_cast<std::size_t>(nodes + 1)
          * static_cast<std::size_t>(nodes + 1),
            0.0);

        // The production residual is linear for this constant-property,
        // prescribed-velocity test:
        //
        //     rhs(T) = -A_H T.
        //
        // Applying the residual to each nodal basis vector recovers one matrix
        // column without reproducing the element formula in the test code.
        for (Int column = 1; column <= nodes; ++column)
        {
            state.temperature.fill(0.0);
            state.temperature1.fill(0.0);
            state.temperature(column) = 1.0;
            state.temperature1(column) = 1.0;

            cbs::EnergyAssembly::assembleStep4Rhs(state);

            for (Int row = 1; row <= nodes; ++row)
            {
                matrix[matrix_index(nodes, row, column)] =
                    -state.rhs1(row);
            }
        }

        return matrix;
    }

    std::vector<Real> make_algebraic_low_order_matrix(
        const std::vector<Real>& high_order,
        const Int nodes)
    {
        std::vector<Real> low_order = high_order;

        for (Int i = 1; i <= nodes; ++i)
        {
            for (Int j = i + 1; j <= nodes; ++j)
            {
                const Real aij =
                    high_order[matrix_index(nodes, i, j)];
                const Real aji =
                    high_order[matrix_index(nodes, j, i)];

                const Real diffusion = std::max({0.0, aij, aji});
                if (diffusion == 0.0)
                {
                    continue;
                }

                // Add the symmetric zero-row-sum diffusion matrix D.
                low_order[matrix_index(nodes, i, j)] -= diffusion;
                low_order[matrix_index(nodes, j, i)] -= diffusion;
                low_order[matrix_index(nodes, i, i)] += diffusion;
                low_order[matrix_index(nodes, j, j)] += diffusion;
            }
        }

        return low_order;
    }

    Real maximum_absolute_row_sum(
        const std::vector<Real>& matrix,
        const Int nodes)
    {
        Real maximum = 0.0;
        for (Int i = 1; i <= nodes; ++i)
        {
            Real row_sum = 0.0;
            for (Int j = 1; j <= nodes; ++j)
            {
                row_sum += matrix[matrix_index(nodes, i, j)];
            }
            maximum = std::max(maximum, std::abs(row_sum));
        }
        return maximum;
    }

    Real maximum_positive_off_diagonal(
        const std::vector<Real>& matrix,
        const Int nodes)
    {
        Real maximum = 0.0;
        for (Int i = 1; i <= nodes; ++i)
        {
            for (Int j = 1; j <= nodes; ++j)
            {
                if (i == j)
                {
                    continue;
                }
                maximum = std::max(
                    maximum,
                    matrix[matrix_index(nodes, i, j)]);
            }
        }
        return maximum;
    }

    SparseRows to_sparse_rows(
        const std::vector<Real>& matrix,
        const Int nodes)
    {
        SparseRows rows(static_cast<std::size_t>(nodes + 1));
        constexpr Real storage_tolerance = 1.0e-24;

        for (Int i = 1; i <= nodes; ++i)
        {
            auto& row = rows[static_cast<std::size_t>(i)];
            for (Int j = 1; j <= nodes; ++j)
            {
                const Real value = matrix[matrix_index(nodes, i, j)];
                if (std::abs(value) > storage_tolerance)
                {
                    row.push_back({j, value});
                }
            }
        }

        return rows;
    }

    Real explicit_positivity_limit(
        const std::vector<Real>& low_order,
        const std::vector<Real>& capacity,
        const Int nodes)
    {
        Real limit = std::numeric_limits<Real>::max();

        for (Int i = 1; i <= nodes; ++i)
        {
            const Real diagonal =
                low_order[matrix_index(nodes, i, i)];

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
                "Step4AfcLowOrder - unable to determine positivity timestep");
        }

        return limit;
    }

    Int nonzero_count(const SparseRows& rows)
    {
        Int count = 0;
        for (std::size_t i = 1; i < rows.size(); ++i)
        {
            count += static_cast<Int>(rows[i].size());
        }
        return count;
    }

    bool run_low_order_case(const Mesh& mesh, const Real dt)
    {
        CBSStateSI state = make_state(mesh, dt);
        const std::vector<Real> capacity =
            assemble_lumped_capacity(state, dt);

        const Int nodes = state.cfg.npoin;
        const std::vector<Real> high_order =
            recover_production_matrix(mesh, dt);
        const std::vector<Real> low_order =
            make_algebraic_low_order_matrix(high_order, nodes);
        const SparseRows sparse_low_order =
            to_sparse_rows(low_order, nodes);

        const Real high_row_sum =
            maximum_absolute_row_sum(high_order, nodes);
        const Real low_row_sum =
            maximum_absolute_row_sum(low_order, nodes);
        const Real maximum_positive_offdiag =
            maximum_positive_off_diagonal(low_order, nodes);
        const Real positivity_dt =
            explicit_positivity_limit(low_order, capacity, nodes);

        const Int steps = static_cast<Int>(std::llround(final_time / dt));
        const Real initial_energy = excess_energy(state, capacity);

        std::vector<Real> residual(
            static_cast<std::size_t>(nodes + 1), 0.0);
        Int first_violation_step = -1;

        for (Int step = 1; step <= steps; ++step)
        {
            for (Int node = 1; node <= nodes; ++node)
            {
                state.temperature1(node) = state.temperature(node);
            }

            for (Int i = 1; i <= nodes; ++i)
            {
                Real matrix_temperature = 0.0;
                for (const SparseEntry& entry :
                     sparse_low_order[static_cast<std::size_t>(i)])
                {
                    matrix_temperature +=
                        entry.value * state.temperature1(entry.column);
                }
                residual[static_cast<std::size_t>(i)] =
                    -matrix_temperature;
            }

            for (Int node = 1; node <= nodes; ++node)
            {
                state.temperature(node) =
                    state.temperature1(node)
                    + dt
                    * residual[static_cast<std::size_t>(node)]
                    / capacity[static_cast<std::size_t>(node)];
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

        const bool matrix_is_monotone =
            maximum_positive_offdiag <= 1.0e-12;
        const bool timestep_is_admissible =
            dt <= positivity_dt * (1.0 + 1.0e-12);
        const bool bounds_preserved =
            final_metrics.below_count == 0
            && final_metrics.above_count == 0
            && first_violation_step < 0;

        std::cout << "\nCASE dt = " << std::scientific << dt << " s\n";
        std::cout << "  steps                         : " << steps << "\n";
        std::cout << "  production-matrix row sum    : "
                  << high_row_sum << "\n";
        std::cout << "  low-order-matrix row sum     : "
                  << low_row_sum << "\n";
        std::cout << "  maximum positive offdiag     : "
                  << maximum_positive_offdiag << "\n";
        std::cout << "  low-order sparse nnz         : "
                  << nonzero_count(sparse_low_order) << "\n";
        std::cout << "  explicit positivity dt limit : "
                  << positivity_dt << " s\n";
        std::cout << "  dt / positivity limit        : "
                  << dt / positivity_dt << "\n";
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
        std::cout << "  M-MATRIX CHECK                : "
                  << (matrix_is_monotone ? "PASS" : "FAIL") << "\n";
        std::cout << "  TIMESTEP CHECK                : "
                  << (timestep_is_admissible ? "PASS" : "FAIL") << "\n";
        std::cout << "  BOUND CHECK                   : "
                  << (bounds_preserved ? "PASS" : "FAIL") << "\n";

        return matrix_is_monotone
            && timestep_is_admissible
            && bounds_preserved;
    }
}

int main()
{
    const Mesh mesh = make_mesh();

    std::cout
        << "CBS3D++_SI Step-4 algebraic low-order validation\n"
        << "Production operator recovered through EnergyAssembly::assembleStep4Rhs\n"
        << "Mesh nodes      : " << mesh.coordinates.size() << "\n"
        << "Mesh tetrahedra : " << mesh.tetrahedra.size() << "\n"
        << "Initial bounds  : [300,400] K\n"
        << "Prescribed U    : (0.05,0,0) m/s\n";

    int failures = 0;
    failures += run_low_order_case(mesh, 1.0e-4) ? 0 : 1;
    failures += run_low_order_case(mesh, 1.0e-5) ? 0 : 1;

    std::cout << "\n============================================================\n";
    if (failures == 0)
    {
        std::cout
            << "ALGEBRAIC LOW-ORDER VALIDATION: ALL CASES PASSED\n"
            << "This establishes a bounded AFC/FCT low-order backbone.\n"
            << "The next stage is limited antidiffusive-flux recovery.\n";
        return 0;
    }

    std::cout << "ALGEBRAIC LOW-ORDER VALIDATION: "
              << failures << " CASE(S) FAILED\n";
    return 1;
}
