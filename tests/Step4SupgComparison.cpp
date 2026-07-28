//=============================================================================
// CBS3D++_SI
//
// Diagnostic comparison of three Step-4 thermal stabilisation choices:
//
//   1. Current production characteristic correction, tau = dt/2.
//   2. No characteristic correction (operator-isolation reference only).
//   3. A Péclet-based streamline-upwind candidate,
//
//          tau = h_u/(2|u|) [coth(Pe) - 1/Pe],
//          Pe  = |u| h_u/(2 alpha),
//          h_u = 2|u| / sum_a |u . grad(N_a)|.
//
// The production EnergyAssembly::assembleStep4Rhs() routine is called first at
// every step. For cases 2 and 3 this diagnostic then removes the production
// tau=dt/2 contribution and replaces it with the requested candidate. No
// production source file is changed by this experiment.
//=============================================================================

#define main cbs3d_step4_bounded_transport_baseline_main
#include "Step4BoundedTransport.cpp"
#undef main

#include <string>

namespace
{
    enum class StabilisationMode
    {
        production_characteristic,
        no_characteristic,
        peclet_supg
    };

    const char* mode_name(const StabilisationMode mode)
    {
        switch (mode)
        {
        case StabilisationMode::production_characteristic:
            return "production characteristic tau=dt/2";
        case StabilisationMode::no_characteristic:
            return "no characteristic correction (diagnostic)";
        case StabilisationMode::peclet_supg:
            return "Peclet-based streamline SUPG candidate";
        }
        return "unknown";
    }

    Real element_gradient(
        const CBSStateSI& state,
        const Int element,
        const Int dimension,
        const Int local_node)
    {
        return state.dNkdx(
            gradient_index(state, element, dimension, local_node));
    }

    Real peclet_supg_tau(
        const CBSStateSI& state,
        const Int element,
        const Real ubar,
        const Real vbar,
        const Real wbar)
    {
        const Real speed = std::sqrt(
            ubar * ubar + vbar * vbar + wbar * wbar);

        if (speed <= 1.0e-14)
        {
            return 0.0;
        }

        Real directional_gradient_sum = 0.0;
        for (Int local_node = 1; local_node <= state.cfg.nep; ++local_node)
        {
            const Real u_grad_shape =
                ubar * element_gradient(state, element, 1, local_node)
              + vbar * element_gradient(state, element, 2, local_node)
              + wbar * element_gradient(state, element, 3, local_node);

            directional_gradient_sum += std::abs(u_grad_shape);
        }

        if (directional_gradient_sum <= 1.0e-30)
        {
            return 0.0;
        }

        const Real streamline_length =
            2.0 * speed / directional_gradient_sum;

        const Real alpha =
            state.k_e(element) / state.rho_cp_e(element);

        if (alpha <= 0.0 || !std::isfinite(alpha))
        {
            throw std::runtime_error(
                "Step4SupgComparison - invalid thermal diffusivity");
        }

        const Real peclet =
            speed * streamline_length / (2.0 * alpha);

        Real xi = 0.0;
        if (peclet < 1.0e-3)
        {
            const Real pe2 = peclet * peclet;
            xi = peclet / 3.0 - peclet * pe2 / 45.0;
        }
        else if (peclet > 20.0)
        {
            xi = 1.0 - 1.0 / peclet;
        }
        else
        {
            xi = 1.0 / std::tanh(peclet) - 1.0 / peclet;
        }

        const Real tau =
            streamline_length * xi / (2.0 * speed);

        if (tau < 0.0 || !std::isfinite(tau))
        {
            throw std::runtime_error(
                "Step4SupgComparison - invalid Péclet-based SUPG tau");
        }

        return tau;
    }

    void replace_production_stabilisation(
        CBSStateSI& state,
        const StabilisationMode mode)
    {
        if (mode == StabilisationMode::production_characteristic)
        {
            return;
        }

        for (Int element = 1; element <= state.cfg.nelem; ++element)
        {
            if (state.mat_elem(element) != 0)
            {
                continue;
            }

            Real dTdx = 0.0;
            Real dTdy = 0.0;
            Real dTdz = 0.0;
            Real ubar = 0.0;
            Real vbar = 0.0;
            Real wbar = 0.0;

            for (Int local_node = 1;
                 local_node <= state.cfg.nep;
                 ++local_node)
            {
                const Int node = state.intma(local_node, element);
                const Real temperature = state.temperature1(node);

                dTdx += element_gradient(
                    state, element, 1, local_node) * temperature;
                dTdy += element_gradient(
                    state, element, 2, local_node) * temperature;
                dTdz += element_gradient(
                    state, element, 3, local_node) * temperature;

                ubar += state.unkno(1, node);
                vbar += state.unkno(2, node);
                wbar += state.unkno(3, node);
            }

            const Real inverse_nodes =
                1.0 / static_cast<Real>(state.cfg.nep);
            ubar *= inverse_nodes;
            vbar *= inverse_nodes;
            wbar *= inverse_nodes;

            const Real advective_temperature_gradient =
                ubar * dTdx + vbar * dTdy + wbar * dTdz;

            const Real old_tau = 0.5 * state.delte(element);
            const Real new_tau =
                (mode == StabilisationMode::peclet_supg)
                ? peclet_supg_tau(state, element, ubar, vbar, wbar)
                : 0.0;

            const Real volume = state.detJ(element) / 6.0;
            const Real delta_factor =
                (new_tau - old_tau)
                * state.rho_cp_e(element)
                * volume;

            for (Int local_node = 1;
                 local_node <= state.cfg.nep;
                 ++local_node)
            {
                const Int node = state.intma(local_node, element);
                const Real u_grad_shape =
                    ubar * element_gradient(state, element, 1, local_node)
                  + vbar * element_gradient(state, element, 2, local_node)
                  + wbar * element_gradient(state, element, 3, local_node);

                state.rhs1(node) +=
                    delta_factor
                    * u_grad_shape
                    * advective_temperature_gradient;
            }
        }
    }

    void report_tau_range(
        const CBSStateSI& state,
        const StabilisationMode mode)
    {
        Real minimum_tau = std::numeric_limits<Real>::max();
        Real maximum_tau = 0.0;

        for (Int element = 1; element <= state.cfg.nelem; ++element)
        {
            Real ubar = 0.0;
            Real vbar = 0.0;
            Real wbar = 0.0;

            for (Int local_node = 1;
                 local_node <= state.cfg.nep;
                 ++local_node)
            {
                const Int node = state.intma(local_node, element);
                ubar += state.unkno(1, node);
                vbar += state.unkno(2, node);
                wbar += state.unkno(3, node);
            }

            const Real inverse_nodes =
                1.0 / static_cast<Real>(state.cfg.nep);
            ubar *= inverse_nodes;
            vbar *= inverse_nodes;
            wbar *= inverse_nodes;

            Real tau = 0.0;
            if (mode == StabilisationMode::production_characteristic)
            {
                tau = 0.5 * state.delte(element);
            }
            else if (mode == StabilisationMode::peclet_supg)
            {
                tau = peclet_supg_tau(
                    state, element, ubar, vbar, wbar);
            }

            minimum_tau = std::min(minimum_tau, tau);
            maximum_tau = std::max(maximum_tau, tau);
        }

        std::cout << "  tau range                     : ["
                  << minimum_tau << ", " << maximum_tau << "] s\n";
    }

    void run_comparison_case(
        const Mesh& mesh,
        const Real dt,
        const StabilisationMode mode)
    {
        const Int steps =
            static_cast<Int>(std::llround(final_time / dt));

        CBSStateSI state = make_state(mesh, dt);
        const std::vector<Real> capacity =
            assemble_lumped_capacity(state, dt);

        const Real initial_energy = excess_energy(state, capacity);
        Int first_violation_step = -1;

        for (Int step = 1; step <= steps; ++step)
        {
            for (Int node = 1; node <= state.cfg.npoin; ++node)
            {
                state.temperature1(node) = state.temperature(node);
            }

            cbs::EnergyAssembly::assembleStep4Rhs(state);
            replace_production_stabilisation(state, mode);

            for (Int node = 1; node <= state.cfg.npoin; ++node)
            {
                state.temperature(node) =
                    state.temperature1(node)
                    + state.rhs1(node) * state.elcoe2p(node);
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

        const Real peclet =
            velocity_x * spacing / (2.0 * thermal_diffusivity);
        const Real courant = velocity_x * dt / spacing;

        std::cout << "\nMODE: " << mode_name(mode) << "\n";
        std::cout << "  dt                            : "
                  << std::scientific << dt << " s\n";
        std::cout << "  steps                         : " << steps << "\n";
        std::cout << "  element Peclet U*h/(2*alpha) : " << peclet << "\n";
        std::cout << "  Courant U*dt/h                : " << courant << "\n";
        report_tau_range(state, mode);
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
        if (first_violation_step >= 0)
        {
            std::cout << "  first violation time          : "
                      << static_cast<Real>(first_violation_step) * dt
                      << " s\n";
        }
        std::cout << "  relative excess-energy drift  : "
                  << relative_energy_drift << "\n";
        std::cout << "  BOUND CHECK                   : "
                  << ((final_metrics.below_count == 0
                       && final_metrics.above_count == 0)
                      ? "PASS" : "FAIL")
                  << "\n";
    }
}

int main()
{
    const Mesh mesh = make_mesh();

    std::cout
        << "CBS3D++_SI Step-4 stabilisation comparison\n"
        << "Production residual assembled before each diagnostic replacement\n"
        << "Mesh nodes      : " << mesh.coordinates.size() << "\n"
        << "Mesh tetrahedra : " << mesh.tetrahedra.size() << "\n"
        << "Initial bounds  : [300,400] K\n"
        << "Prescribed U    : (0.05,0,0) m/s\n";

    for (const Real dt : {1.0e-4, 1.0e-5})
    {
        std::cout << "\n============================================================\n";
        std::cout << "COMPARISON AT dt = "
                  << std::scientific << dt << " s\n";

        run_comparison_case(
            mesh,
            dt,
            StabilisationMode::production_characteristic);
        run_comparison_case(
            mesh,
            dt,
            StabilisationMode::no_characteristic);
        run_comparison_case(
            mesh,
            dt,
            StabilisationMode::peclet_supg);
    }

    std::cout << "\n============================================================\n";
    std::cout
        << "The Péclet-SUPG case is a validation candidate, not a production\n"
        << "change. A passing result here must still be followed by smooth\n"
        << "accuracy, CHT-interface, timestep, and MPI-equivalence tests.\n";

    return 0;
}
