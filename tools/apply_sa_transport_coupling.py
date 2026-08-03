from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return path.read_text(encoding="utf-8")


def write(path, text):
    path.write_text(text, encoding="utf-8")


def replace_once(text, old, new, label):
    if new in text:
        print(f"[skip] {label}")
        return text
    if old not in text:
        raise SystemExit(f"Patch anchor not found: {label}")
    print(f"[patch] {label}")
    return text.replace(old, new, 1)


def patch_solver():
    path = ROOT / "src" / "solver" / "Solver.cpp"
    text = read(path)

    old_copy = """        s_.unkn1 = s_.unkno;
        s_.pres1 = s_.pres;
        s_.temperature1 = s_.temperature;
"""
    new_copy = """        s_.unkn1 = s_.unkno;
        s_.pres1 = s_.pres;
        s_.temperature1 = s_.temperature;

        // Store the SA working variable at the beginning of the current
        // iteration.  The SA residual assembly is explicit in nu_tilde and
        // therefore uses nu_tilde1 as q^n.
        if (s_.cfg.turbulence_on > 0)
        {
            s_.nu_tilde1 = s_.nu_tilde;
        }
"""
    text = replace_once(text, old_copy, new_copy, "Solver.cpp: copy nu_tilde1")

    old_step4 = """        // CBS Step 4: advance the energy equation when temperature is enabled.
        {
            auto timer = profiler_.time(SolverProfiler::Section::Step4Energy);
            Steps::step4(s_);
        }
"""
    new_step4 = """        // Optional Spalart-Allmaras transport step.  It is placed after
        // Step 3 because the SA advection and production terms use the corrected
        // velocity field.  It is placed before Step 4 so that the energy equation
        // can use the updated turbulent thermal conductivity.
        if (s_.cfg.turbulence_on > 0)
        {
            auto timer = profiler_.time(SolverProfiler::Section::StepSpalartAllmaras);
            Steps::stepSpalartAllmaras(s_);
        }

        // CBS Step 4: advance the energy equation when temperature is enabled.
        {
            auto timer = profiler_.time(SolverProfiler::Section::Step4Energy);
            Steps::step4(s_);
        }
"""
    text = replace_once(text, old_step4, new_step4, "Solver.cpp: insert Step SA")

    write(path, text)


def patch_momentum():
    path = ROOT / "src" / "assembly" / "MomentumAssembly.cpp"
    text = read(path)

    old = """        // Returns the momentum diffusivity used in the viscous term.
        //
        // Dimensional material mode:
        //
        //     nu_e = mu_e / rho_e
        //
        // Non-dimensional mode:
        //
        //     nu_e = ani
        Real momentum_diffusivity(
            const CBSStateSI& s,
            Int ie)
        {
            if (s.cfg.dimensional_mode > 0 && s.cfg.material_properties_enabled > 0)
            {
                if (s.rho_e(ie) <= 0.0 || s.mu_e(ie) < 0.0)
                {
                    throw std::runtime_error(
                        "MomentumAssembly::assembleStep1Rhs - invalid fluid rho/mu at element "
                        + std::to_string(ie));
                }

                return s.mu_e(ie) / s.rho_e(ie);
            }

            return s.cfg.ani;
        }
"""
    new = """        // Returns the kinematic viscosity used in the Step-1 viscous term.
        //
        // Laminar dimensional mode:
        //
        //     nu_e = mu_e / rho_e
        //
        // Spalart-Allmaras dimensional mode:
        //
        //     nu_eff,e = mu_eff_e / rho_e
        //
        // where:
        //
        //     mu_eff_e = mu_e + mu_t_e
        //
        // Non-dimensional laminar mode continues to use ani.  The first SA
        // implementation is intended for dimensional material cases, where rho
        // and mu are available element by element.
        Real momentum_diffusivity(
            const CBSStateSI& s,
            Int ie)
        {
            if (s.cfg.dimensional_mode > 0 && s.cfg.material_properties_enabled > 0)
            {
                if (s.rho_e(ie) <= 0.0 || !std::isfinite(s.rho_e(ie)))
                {
                    throw std::runtime_error(
                        "MomentumAssembly::assembleStep1Rhs - invalid fluid density at element "
                        + std::to_string(ie));
                }

                Real mu = s.mu_e(ie);

                if (s.cfg.turbulence_on > 0)
                {
                    mu = s.mu_eff_e(ie);
                }

                if (mu < 0.0 || !std::isfinite(mu))
                {
                    throw std::runtime_error(
                        "MomentumAssembly::assembleStep1Rhs - invalid effective viscosity at element "
                        + std::to_string(ie));
                }

                return mu / s.rho_e(ie);
            }

            return s.cfg.ani;
        }
"""
    text = replace_once(text, old, new, "MomentumAssembly.cpp: effective viscosity")
    write(path, text)


def patch_energy():
    path = ROOT / "src" / "assembly" / "EnergyAssembly.cpp"
    text = read(path)

    old = """            const Real volume = s.detJ(ie) / 6.0;
            const Real k = s.k_e(ie);

            for (Int a = 1; a <= s.cfg.nep; ++a)
"""
    new = """            const Real volume = s.detJ(ie) / 6.0;

            // Laminar and solid calculation:
            //
            //     k_used = k_e
            //
            // Turbulent fluid heat-transfer calculation:
            //
            //     k_used = k_eff_e
            //
            // where:
            //
            //     k_eff_e = k_e + rho cp nu_t / Pr_t
            //
            // The turbulent addition is never applied in solid elements.
            Real k = s.k_e(ie);

            if (s.cfg.turbulence_on > 0 &&
                s.cfg.turbulent_thermal_diffusivity_on > 0 &&
                is_fluid_element(s, ie))
            {
                k = s.k_eff_e(ie);
            }

            if (k <= 0.0 || !std::isfinite(k))
            {
                throw std::runtime_error(
                    "EnergyAssembly::assembleStep4Rhs - invalid effective thermal conductivity at element "
                    + std::to_string(ie));
            }

            for (Int a = 1; a <= s.cfg.nep; ++a)
"""
    text = replace_once(text, old, new, "EnergyAssembly.cpp: effective thermal conductivity")
    write(path, text)


def patch_timestep():
    path = ROOT / "src" / "timestep" / "TimeStep.cpp"
    text = read(path)

    old = """                    const Real nu = s.mu_e(ie) / s.rho_e(ie);
                    diff = std::max(diff, nu);
"""
    new = """                    Real nu = s.mu_e(ie) / s.rho_e(ie);

                    if (s.cfg.turbulence_on > 0)
                    {
                        if (s.mu_eff_e(ie) < 0.0 || !std::isfinite(s.mu_eff_e(ie)))
                        {
                            throw std::runtime_error(
                                "TimeStep - invalid effective viscosity at element "
                                + std::to_string(ie));
                        }

                        nu = s.mu_eff_e(ie) / s.rho_e(ie);

                        if (s.cfg.turbulent_thermal_diffusivity_on > 0)
                        {
                            if (s.k_eff_e(ie) <= 0.0 || !std::isfinite(s.k_eff_e(ie)))
                            {
                                throw std::runtime_error(
                                    "TimeStep - invalid effective thermal conductivity at element "
                                    + std::to_string(ie));
                            }

                            const Real alpha_eff = s.k_eff_e(ie) / s.rho_cp_e(ie);
                            diff = std::max(diff, alpha_eff);
                        }
                    }

                    diff = std::max(diff, nu);
"""
    text = replace_once(text, old, new, "TimeStep.cpp: effective turbulent diffusivity")
    write(path, text)


def patch_sa_assembly():
    path = ROOT / "src" / "assembly" / "SpalartAllmarasAssembly.cpp"
    text = read(path)

    old_decl = """                Real local_rhs[5] = {};
                Real local_production[5] = {};
                Real local_destruction[5] = {};
                Real local_diffusion[5] = {};
                Real local_source[5] = {};
"""
    new_decl = """                Real local_rhs[5] = {};
                Real local_production[5] = {};
                Real local_destruction[5] = {};
                Real local_diffusion[5] = {};
                Real local_source[5] = {};

                // local_destruction_lhs stores the nodal coefficient added to
                // the left-hand side when positive destruction is treated
                // semi-implicitly:
                //
                //     D = C_D q^2  ->  C_D q_old q_new
                //
                // The global temporary storage is sa_residual during assembly.
                // It is overwritten by the actual nodal update after updateNuTilde().
                Real local_destruction_lhs[5] = {};
"""
    text = replace_once(text, old_decl, new_decl, "SpalartAllmarasAssembly.cpp: local implicit destruction array")

    old_terms = """                    const Real production_rhs = nodal_volume * production;
                    const Real destruction_rhs = nodal_volume * destruction;
                    const Real nonlinear_rhs = nodal_volume * nonlinear_gradient_source;

                    local_rhs[a] += advection;
                    local_rhs[a] += diffusion;
                    local_rhs[a] += production_rhs;
                    local_rhs[a] -= destruction_rhs;
                    local_rhs[a] += nonlinear_rhs;

                    local_production[a] += production_rhs;
                    local_destruction[a] += destruction_rhs;
                    local_diffusion[a] += diffusion + nonlinear_rhs;
                    local_source[a] += production_rhs - destruction_rhs;
"""
    new_terms = """                    const Real production_rhs = nodal_volume * production;
                    const Real destruction_rhs = nodal_volume * destruction;
                    const Real nonlinear_rhs = nodal_volume * nonlinear_gradient_source;

                    local_rhs[a] += advection;
                    local_rhs[a] += diffusion;
                    local_rhs[a] += production_rhs;
                    local_rhs[a] += nonlinear_rhs;

                    // Destruction is stiff near the wall because it contains
                    // (nu_tilde/d)^2.  If sa_implicit_destruction is enabled and
                    // the destruction coefficient is positive, use the standard
                    // first-order linearisation:
                    //
                    //     C_D q^2  ->  C_D q_old q_new
                    //
                    // The term C_D q_old is therefore added to the nodal
                    // mass/time diagonal.  Negative C_D is not a destruction
                    // sink; it is kept explicit in the RHS with the correct sign.
                    if (s.cfg.sa_implicit_destruction > 0 &&
                        destruction_coefficient > 0.0 &&
                        q_average > 0.0)
                    {
                        local_destruction_lhs[a] +=
                            nodal_volume * destruction_coefficient * q_average;
                    }
                    else
                    {
                        local_rhs[a] -= destruction_rhs;
                    }

                    local_production[a] += production_rhs;
                    local_destruction[a] += destruction_rhs;
                    local_diffusion[a] += diffusion + nonlinear_rhs;
                    local_source[a] += production_rhs - destruction_rhs;
"""
    text = replace_once(text, old_terms, new_terms, "SpalartAllmarasAssembly.cpp: semi-implicit destruction")

    old_scatter = """                    s.sa_destruction(ip) += local_destruction[a];
                    s.sa_diffusion(ip) += local_diffusion[a];
                    s.sa_source(ip) += local_source[a];
"""
    new_scatter = """                    s.sa_destruction(ip) += local_destruction[a];
                    s.sa_diffusion(ip) += local_diffusion[a];
                    s.sa_source(ip) += local_source[a];
                    s.sa_residual(ip) += local_destruction_lhs[a];
"""
    text = replace_once(text, old_scatter, new_scatter, "SpalartAllmarasAssembly.cpp: scatter implicit destruction coefficient")

    old_update = """            const Real old_value = s.nu_tilde1(ip);
            Real new_value = old_value + s.elcoe2(ip) * s.sa_rhs(ip);

            if (!std::isfinite(new_value))
            {
                new_value = old_value;
            }

            new_value = std::max(s.cfg.sa_nu_tilde_floor, new_value);

            s.nu_tilde(ip) = new_value;
            s.sa_residual(ip) = new_value - old_value;
"""
    new_update = """            const Real old_value = s.nu_tilde1(ip);

            Real denominator = 1.0;
            if (s.cfg.sa_implicit_destruction > 0)
            {
                denominator += s.elcoe2(ip) * s.sa_residual(ip);
            }

            if (denominator <= 0.0 || !std::isfinite(denominator))
            {
                denominator = 1.0;
            }

            Real new_value =
                (old_value + s.elcoe2(ip) * s.sa_rhs(ip)) / denominator;

            if (!std::isfinite(new_value))
            {
                new_value = old_value;
            }

            new_value = std::max(s.cfg.sa_nu_tilde_floor, new_value);

            s.nu_tilde(ip) = new_value;
            s.sa_residual(ip) = new_value - old_value;
"""
    text = replace_once(text, old_update, new_update, "SpalartAllmarasAssembly.cpp: implicit update denominator")

    write(path, text)


def main():
    patch_solver()
    patch_momentum()
    patch_energy()
    patch_timestep()
    patch_sa_assembly()
    print("SA transport coupling patch applied.")


if __name__ == "__main__":
    main()
