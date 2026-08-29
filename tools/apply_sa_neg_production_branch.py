#!/usr/bin/env python3
"""Apply the NASA/TMR SA-neg branch to the Liu/Nithiarasu SA production path.

The patch is intentionally narrow:

* turbulence_model=0 remains the existing Liu/SA-noft2 implementation;
* turbulence_model=1 becomes the standard NASA/TMR SA-neg model;
* positive SA-neg states use standard SA including ft2 and the AJS S-tilde guard;
* negative SA-neg states use the Allmaras-Johnson-Spalart negative production,
  recovery/destruction and diffusion coefficient;
* negative nu_tilde is preserved by the update instead of being clipped to zero;
* momentum still receives zero eddy viscosity from negative SA states;
* the SA explicit timestep uses (nu + nu_tilde*f_n)/sigma on the negative branch.

The CBS conservative scalar transport/characteristic discretisation is retained.
For incompressible flow it is the conservation-form counterpart of the standard
SA transport equation, and the existing frozen-field audit showed its discrete
divergence correction is far too small to explain the flat-plate error.
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ASSEMBLY = ROOT / "src" / "assembly" / "SpalartAllmarasAssemblyLiu.cpp"
TIMESTEP = ROOT / "src" / "timestep" / "TimeStep.cpp"
CMAKE = ROOT / "CMakeLists.txt"


def replace_once(text: str, old: str, new: str, description: str) -> str:
    if new in text:
        print("[skip] {} already applied".format(description))
        return text
    count = text.count(old)
    if count != 1:
        raise SystemExit(
            "FATAL: expected exactly one anchor for {}; found {}".format(
                description, count
            )
        )
    print("[patch] {}".format(description))
    return text.replace(old, new, 1)


def write_text(path, text):
    with path.open("w", encoding="utf-8", newline="\n") as output:
        output.write(text)


def patch_assembly() -> None:
    text = ASSEMBLY.read_text(encoding="utf-8")

    old = '''                    q_average += s.sa_wall_node(ip) != 0
                        ? 0.0
                        : std::max(s.cfg.sa_nu_tilde_floor, s.nu_tilde(ip));
'''
    new = '''                    const Real q_node = s.sa_wall_node(ip) != 0
                        ? 0.0
                        : s.nu_tilde(ip);
                    q_average += s.cfg.turbulence_model == 1
                        ? q_node
                        : std::max(s.cfg.sa_nu_tilde_floor, q_node);
'''
    text = replace_once(
        text, old, new,
        "preserve signed nu_tilde in SA-neg eddy-viscosity preprocessing")

    old = '''                    q[a] = s.sa_wall_node(ip) != 0
                        ? 0.0
                        : std::max(s.cfg.sa_nu_tilde_floor, s.nu_tilde1(ip));
'''
    new = '''                    const Real q_node = s.sa_wall_node(ip) != 0
                        ? 0.0
                        : s.nu_tilde1(ip);
                    q[a] = s.cfg.turbulence_model == 1
                        ? q_node
                        : std::max(s.cfg.sa_nu_tilde_floor, q_node);
'''
    text = replace_once(
        text, old, new,
        "preserve signed transported field in SA-neg assembly")

    old = '''                const Real omega = vorticity_magnitude(s, ie);
                const Real s_bar = q_average > 0.0
                    ? turbulence::sBar(
                        q_average,
                        molecular_nu,
                        wall_distance,
                        constants)
                    : 0.0;

                // Liu Eq. 3.26.  No ft2 term and no AJS S-tilde branch appears
                // in the thesis formulation.  The tiny floor is only a division
                // guard for r and does not act in ordinary attached flow.
                const Real s_tilde = std::max(
                    omega + s_bar,
                    s.cfg.sa_min_stilde);

                const Real r_value = turbulence::rFunction(
                    q_average,
                    s_tilde,
                    wall_distance,
                    constants);
                const Real fw_value = turbulence::fw(r_value, constants);

                // Passing ft2=0 selects the fully turbulent Liu/SA-noft2 model
                // while retaining the well-tested scalar helper functions.
                const Real production = turbulence::productionTerm(
                    q_average,
                    s_tilde,
                    0.0,
                    constants);
                const Real destruction_coefficient =
                    turbulence::destructionCoefficient(
                        wall_distance,
                        fw_value,
                        0.0,
                        constants);

                const Real diffusion_coefficient =
                    (molecular_nu + q_average) / constants.sigma;
                const Real nonlinear_source =
                    constants.cb2 * grad_q_squared / constants.sigma;
'''
    new = '''                const Real omega = vorticity_magnitude(s, ie);
                const bool sa_neg_model = s.cfg.turbulence_model == 1;
                const bool negative_branch = sa_neg_model && q_average < 0.0;

                Real production_coefficient = 0.0;
                Real destruction_coefficient = 0.0;
                Real diffusion_coefficient = 0.0;

                if (negative_branch)
                {
                    // NASA/TMR SA-neg, Allmaras-Johnson-Spalart (2012):
                    //
                    //   P- = cb1 (1-ct3) Omega nu_tilde
                    //   D- = +cw1 (nu_tilde/d)^2
                    //   K- = (nu + nu_tilde f_n) / sigma
                    //
                    // D- is a recovery source on the negative branch; the plus
                    // sign is intentional and is handled below when assembling
                    // the source vector.
                    const Real chi_value =
                        turbulence::chi(q_average, molecular_nu);
                    const Real fn_value =
                        turbulence::negativeBranchFn(chi_value, constants);

                    production_coefficient =
                        constants.cb1 * (1.0 - constants.ct3) * omega;
                    destruction_coefficient =
                        constants.cw1() /
                        (wall_distance * wall_distance);
                    diffusion_coefficient =
                        (molecular_nu + q_average * fn_value) /
                        constants.sigma;
                }
                else
                {
                    const Real s_bar = q_average > 0.0
                        ? turbulence::sBar(
                            q_average,
                            molecular_nu,
                            wall_distance,
                            constants)
                        : 0.0;

                    Real s_tilde = 0.0;
                    Real ft2_value = 0.0;

                    if (sa_neg_model)
                    {
                        // SA-neg is identical to standard SA for nu_tilde >= 0.
                        // The AJS protected S-tilde form is considered standard
                        // for SA-neg and avoids clipping S-tilde itself.
                        s_tilde = turbulence::limitedSTilde(
                            omega,
                            s_bar,
                            constants);
                        ft2_value = turbulence::ft2(
                            turbulence::chi(q_average, molecular_nu),
                            constants);
                    }
                    else
                    {
                        // Preserve the established Liu/SA-noft2 branch exactly.
                        s_tilde = std::max(
                            omega + s_bar,
                            s.cfg.sa_min_stilde);
                        ft2_value = 0.0;
                    }

                    const Real r_value = turbulence::rFunction(
                        q_average,
                        s_tilde,
                        wall_distance,
                        constants);
                    const Real fw_value = turbulence::fw(r_value, constants);

                    production_coefficient =
                        constants.cb1 * (1.0 - ft2_value) * s_tilde;
                    destruction_coefficient =
                        turbulence::destructionCoefficient(
                            wall_distance,
                            fw_value,
                            ft2_value,
                            constants);
                    diffusion_coefficient =
                        (molecular_nu + q_average) / constants.sigma;
                }

                if (!(diffusion_coefficient > 0.0) ||
                    !std::isfinite(diffusion_coefficient))
                {
#ifdef CBS3D_USE_OPENMP
#pragma omp atomic write
#endif
                    bad_material = 1;
                    continue;
                }

                const Real nonlinear_source =
                    constants.cb2 * grad_q_squared / constants.sigma;
'''
    text = replace_once(
        text, old, new,
        "add positive and negative SA-neg source/diffusion closures")

    old = '''                    const Real production_rhs =
                        (q_average > 0.0 ? production / q_average : 0.0)
                        * q_load;
                    const Real destruction_rhs =
                        destruction_coefficient * q_average * q_load;

                    local_rhs[a] =
                        advection
                        + characteristic
                        + diffusion
                        + nonlinear
                        + production_rhs
                        - destruction_rhs;

                    local_production[a] = production_rhs;
                    local_destruction[a] = destruction_rhs;
                    local_diffusion[a] = diffusion + nonlinear;
                    local_source[a] = production_rhs - destruction_rhs;
'''
    new = '''                    const Real production_rhs =
                        production_coefficient * q_load;
                    const Real destruction_rhs =
                        destruction_coefficient * q_average * q_load;
                    const Real source_rhs = negative_branch
                        ? production_rhs + destruction_rhs
                        : production_rhs - destruction_rhs;

                    local_rhs[a] =
                        advection
                        + characteristic
                        + diffusion
                        + nonlinear
                        + source_rhs;

                    local_production[a] = production_rhs;
                    // On the SA-neg negative branch this array stores the
                    // positive recovery term cw1*(nu_tilde/d)^2 even though the
                    // historical diagnostic name is "destruction".
                    local_destruction[a] = destruction_rhs;
                    local_diffusion[a] = diffusion + nonlinear;
                    local_source[a] = source_rhs;
'''
    text = replace_once(
        text, old, new,
        "assemble SA-neg recovery with the correct positive sign")

    old = '''            // Liu's Step 4 is explicit.  sa_implicit_destruction is deliberately
            // not used by this production assembly.  The non-negative projection
            // remains a guard for the present non-SA-neg branch and must be
            // monitored during verification rather than mistaken for physics.
            new_value = std::max(s.cfg.sa_nu_tilde_floor, new_value);
            s.nu_tilde(ip) = new_value;
            s.sa_residual(ip) = new_value - old_value;
'''
    new = '''            // Liu's Step 4 is explicit.  For the legacy/model-0
            // SA-noft2 path retain the historical non-negative safeguard.  A
            // genuine SA-neg calculation must preserve finite negative values so
            // that the negative recovery equation can act on the next iteration.
            if (s.cfg.turbulence_model != 1)
            {
                new_value = std::max(s.cfg.sa_nu_tilde_floor, new_value);
            }

            s.nu_tilde(ip) = new_value;
            s.sa_residual(ip) = new_value - old_value;
'''
    text = replace_once(
        text, old, new,
        "remove hard zero projection from turbulence_model=1")

    write_text(ASSEMBLY, text)


def patch_timestep() -> None:
    text = TIMESTEP.read_text(encoding="utf-8")

    old = '''                        const Real nu_tilde_e = s.nu_tilde_e(ie);

                        if (nu_tilde_e < 0.0 || !std::isfinite(nu_tilde_e))
                        {
                            throw std::runtime_error(
                                "TimeStep - invalid nu_tilde_e at element "
                                + std::to_string(ie));
                        }

                        const Real molecular_nu = s.mu_e(ie) / s.rho_e(ie);

                        const Real sa_diffusivity =
                            (molecular_nu + nu_tilde_e) / sa_constants.sigma;

                        diff = std::max(diff, sa_diffusivity);
'''
    new = '''                        const Real nu_tilde_e = s.nu_tilde_e(ie);

                        if (!std::isfinite(nu_tilde_e))
                        {
                            throw std::runtime_error(
                                "TimeStep - non-finite nu_tilde_e at element "
                                + std::to_string(ie));
                        }

                        const Real molecular_nu = s.mu_e(ie) / s.rho_e(ie);
                        Real sa_transport_viscosity = nu_tilde_e;

                        if (s.cfg.turbulence_model == 1 && nu_tilde_e < 0.0)
                        {
                            const Real chi_value = turbulence::chi(
                                nu_tilde_e,
                                molecular_nu);
                            const Real fn_value = turbulence::negativeBranchFn(
                                chi_value,
                                sa_constants);
                            sa_transport_viscosity = nu_tilde_e * fn_value;
                        }
                        else if (nu_tilde_e < 0.0)
                        {
                            throw std::runtime_error(
                                "TimeStep - negative nu_tilde_e in non-SA-neg model at element "
                                + std::to_string(ie));
                        }

                        const Real sa_diffusivity =
                            (molecular_nu + sa_transport_viscosity) /
                            sa_constants.sigma;

                        if (!(sa_diffusivity > 0.0) ||
                            !std::isfinite(sa_diffusivity))
                        {
                            throw std::runtime_error(
                                "TimeStep - invalid SA transport diffusivity at element "
                                + std::to_string(ie));
                        }

                        diff = std::max(diff, sa_diffusivity);
'''
    text = replace_once(
        text, old, new,
        "use SA-neg negative diffusion coefficient in timestep stability bound")

    write_text(TIMESTEP, text)


def patch_cmake() -> None:
    text = CMAKE.read_text(encoding="utf-8")
    marker = "# SA-neg deterministic negative-branch verification"
    if marker in text:
        print("[skip] SA-neg CMake verification target already present")
        return

    addition = r'''

#==============================================================================
# SA-neg deterministic negative-branch verification
#==============================================================================
if (CBS3D_BUILD_TESTS)
    add_executable(sa_neg_element_verification
        tests/SANegElementVerification.cpp
        src/assembly/SpalartAllmarasAssemblyLiu.cpp
        src/turbulence/SpalartAllmaras.cpp
    )

    target_include_directories(sa_neg_element_verification PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
    )

    target_compile_features(sa_neg_element_verification PRIVATE cxx_std_20)

    if (MSVC)
        target_compile_options(sa_neg_element_verification PRIVATE
            /W4
            /permissive-
        )
    else()
        target_compile_options(sa_neg_element_verification PRIVATE
            -O2
            -Wall
            -Wextra
            -pedantic
        )
    endif()

    add_test(
        NAME sa_neg_negative_branch
        COMMAND sa_neg_element_verification
    )
endif()
'''
    print("[patch] add SA-neg deterministic CMake verification target")
    write_text(CMAKE, text.rstrip() + addition + "\n")


def main() -> int:
    patch_assembly()
    patch_timestep()
    patch_cmake()
    print("SA-neg production patch applied successfully.")
    print("Review with: git diff --check && git diff")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
