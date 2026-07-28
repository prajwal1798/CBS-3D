//=============================================================================
// CBS3D++_SI
//
// End-to-end boundedness validation of the production ThermalAfc class.
// The test uses the same structured tetrahedral pulse problem as the earlier
// diagnostics but now executes:
//
//     EnergyAssembly::assembleStep4Rhs()
//     ThermalAfc::applySerial()
//
// at every step. This is the code path integrated into the solver.
//=============================================================================

#define main cbs3d_step4_bounded_transport_baseline_main
#include "Step4BoundedTransport.cpp"
#undef main

#include "cbs/assembly/ThermalAfc.hpp"

namespace
{
    bool run_integrated_case(const Mesh& mesh, const Real dt)
    {
        CBSStateSI state = make_state(mesh, dt);
        const std::vector<Real> capacity = assemble_lumped_capacity(state, dt);
        const Int steps = static_cast<Int>(std::llround(final_time / dt));
        const Real initial_energy = excess_energy(state, capacity);

        Int first_violation_step = -1;

        for (Int step = 1; step <= steps; ++step)
        {
            for (Int node = 1; node <= state.cfg.npoin; ++node)
            {
                state.temperature1(node) = state.temperature(node);
            }

            cbs::EnergyAssembly::assembleStep4Rhs(state);
            cbs::ThermalAfc::applySerial(state);
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
            (final_energy - initial_energy) / initial_energy;

        const bool bounds_preserved =
            first_violation_step < 0
            && final_metrics.below_count == 0
            && final_metrics.above_count == 0;

        const bool energy_preserved =
            std::abs(relative_energy_drift) <= 1.0e-5;

        std::cout << "\nCASE dt = " << std::scientific << dt << " s\n";
        std::cout << "  steps                         : " << steps << "\n";
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
        std::cout << "  ENERGY CHECK                  : "
                  << (energy_preserved ? "PASS" : "FAIL") << "\n";
        std::cout << "  BOUND CHECK                   : "
                  << (bounds_preserved ? "PASS" : "FAIL") << "\n";

        return bounds_preserved && energy_preserved;
    }
}

int main()
{
    const Mesh mesh = make_mesh();

    std::cout
        << "CBS3D++_SI integrated production AFC validation\n"
        << "Path: EnergyAssembly -> ThermalAfc::applySerial\n"
        << "Mesh nodes      : " << mesh.coordinates.size() << "\n"
        << "Mesh tetrahedra : " << mesh.tetrahedra.size() << "\n";

    int failures = 0;
    failures += run_integrated_case(mesh, 1.0e-4) ? 0 : 1;
    failures += run_integrated_case(mesh, 1.0e-5) ? 0 : 1;

    std::cout << "\n============================================================\n";
    if (failures == 0)
    {
        std::cout << "INTEGRATED PRODUCTION AFC VALIDATION: ALL CASES PASSED\n";
        return 0;
    }

    std::cout << "INTEGRATED PRODUCTION AFC VALIDATION: "
              << failures << " CASE(S) FAILED\n";
    return 1;
}
