//=============================================================================
// CBS3D++_SI
// INDEPENDENT P1 TETRAHEDRAL STEP-4 OPERATOR VERIFICATION.
//
// THE EXPECTED FLUID RHS IS EVALUATED WITH A DEGREE-2 TETRAHEDRAL
// QUADRATURE RULE, INDEPENDENT OF ENERGYASSEMBLY'S CLOSED-FORM
// COEFFICIENTS. NO PRODUCTION SOLVER CODE IS CHANGED.
//=============================================================================

#include "cbs/assembly/EnergyAssembly.hpp"
#include "cbs/core/CBSStateSI.hpp"

#include <array>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

namespace
{
    using cbs::CBSStateSI;
    using cbs::Int;
    using cbs::Real;

    using Vector = std::array<Real, 3>;
    using NodalVector = std::array<Vector, 4>;
    using NodalScalar = std::array<Real, 4>;

    const NodalVector gradients = {{
        {{-1.0, -1.0, -1.0}},
        {{ 1.0,  0.0,  0.0}},
        {{ 0.0,  1.0,  0.0}},
        {{ 0.0,  0.0,  1.0}}
    }};

    const NodalVector velocities = {{
        {{ 0.70, -0.10,  0.20}},
        {{ 0.60,  0.10,  0.00}},
        {{ 0.20,  0.10, -0.20}},
        {{ 0.30, -0.40,  0.50}}
    }};

    const NodalScalar temperatures = {{301.0, 312.0, 323.0, 333.0}};

    Real dot(const Vector& a, const Vector& b)
    {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    bool close_to(const Real actual, const Real expected,
                  const Real relative_tolerance = 2.0e-10)
    {
        return std::abs(actual - expected) <=
            relative_tolerance * std::max(1.0, std::abs(expected));
    }

    void print_result(const std::string& title, const bool passed)
    {
        std::cout << title << ": " << (passed ? "PASS" : "FAIL") << "\n";
    }

    CBSStateSI reference_tetrahedron(const Int material_id,
                                    const bool with_boundary = false)
    {
        CBSStateSI s;
        s.set_problem_sizes(1, 4, with_boundary ? 1 : 0, 0);
        s.cfg.temp_calc = 1;
        s.cfg.source_solid = 0.0;
        s.cfg.heat_flux_bc = 0.0;
        s.cfg.turbulence_on = 0;
        s.cfg.turbulent_thermal_diffusivity_on = 0;
        s.cfg.inlet_temperature = 298.15;

        for (Int a = 1; a <= 4; ++a)
        {
            s.intma(a, 1) = a;
            for (Int d = 1; d <= 3; ++d)
            {
                s.dNkdx((d - 1) * 4 + a) =
                    gradients[static_cast<std::size_t>(a - 1)]
                             [static_cast<std::size_t>(d - 1)];
            }
        }

        s.detJ(1) = 1.0;  // UNIT REFERENCE TETRAHEDRON HAS V=1/6.
        s.delte(1) = 0.017;
        s.mat_elem(1) = material_id;
        s.rho_e(1) = 1.0;
        s.cp_e(1) = 1800.0;
        s.rho_cp_e(1) = 1800.0;
        s.k_e(1) = 5.7;
        s.k_eff_e(1) = 5.7;
        s.Qvol_e(1) = 0.0;
        s.temperature.fill(298.15);
        s.temperature1.fill(298.15);
        s.unkno.fill(0.0);
        s.rhs1.fill(0.0);
        return s;
    }

    void set_field(CBSStateSI& s)
    {
        for (Int a = 1; a <= 4; ++a)
        {
            const std::size_t i = static_cast<std::size_t>(a - 1);
            s.temperature1(a) = temperatures[i];
            for (Int d = 1; d <= 3; ++d)
            {
                s.unkno(d, a) = velocities[i][static_cast<std::size_t>(d - 1)];
            }
        }
    }

    // EXACT DEGREE-2 QUADRATURE FOR THE P1-INTERPOLATED CONVECTIVE
    // INTEGRAND. THE CHARACTERISTIC TERM USES THE CENTROIDAL VELOCITY
    // BY DEFINITION OF THE CODE'S CHOSEN APPROXIMATION.
    NodalScalar independent_fluid_rhs(const CBSStateSI& s)
    {
        NodalScalar result = {{0.0, 0.0, 0.0, 0.0}};
        Vector grad_t = {{0.0, 0.0, 0.0}};
        Vector u_bar = {{0.0, 0.0, 0.0}};
        Real divergence = 0.0;

        for (std::size_t a = 0; a < 4; ++a)
        {
            const Real t = s.temperature1(static_cast<Int>(a + 1));
            Vector velocity = {{
                s.unkno(1, static_cast<Int>(a + 1)),
                s.unkno(2, static_cast<Int>(a + 1)),
                s.unkno(3, static_cast<Int>(a + 1))
            }};

            for (std::size_t d = 0; d < 3; ++d)
            {
                grad_t[d] += t * gradients[a][d];
                u_bar[d] += velocity[d] / 4.0;
            }
            divergence += dot(velocity, gradients[a]);
        }

        const Real alpha = (5.0 + 3.0 * std::sqrt(5.0)) / 20.0;
        const Real beta = (5.0 - std::sqrt(5.0)) / 20.0;
        const Real volume = 1.0 / 6.0;
        const Real capacity = s.rho_cp_e(1);
        const Real conductivity = s.k_e(1);
        const Real tref = s.cfg.inlet_temperature;
        const Real tau = s.delte(1) / 2.0;

        for (std::size_t p = 0; p < 4; ++p)
        {
            NodalScalar shape = {{beta, beta, beta, beta}};
            shape[p] = alpha;
            Vector u_q = {{0.0, 0.0, 0.0}};
            Real t_q = 0.0;

            for (std::size_t b = 0; b < 4; ++b)
            {
                const Int node = static_cast<Int>(b + 1);
                t_q += shape[b] * s.temperature1(node);
                for (std::size_t d = 0; d < 3; ++d)
                {
                    u_q[d] += shape[b] * s.unkno(static_cast<Int>(d + 1), node);
                }
            }

            for (std::size_t a = 0; a < 4; ++a)
            {
                result[a] -= capacity * (volume / 4.0) * shape[a] *
                    (dot(u_q, grad_t) + (t_q - tref) * divergence);
            }
        }

        for (std::size_t a = 0; a < 4; ++a)
        {
            result[a] -= conductivity * volume * dot(gradients[a], grad_t);
            result[a] -= tau * capacity * volume *
                dot(u_bar, gradients[a]) * dot(u_bar, grad_t);
        }
        return result;
    }

    bool test_uniform_temperature_invariance()
    {
        CBSStateSI s = reference_tetrahedron(0);
        set_field(s);
        s.temperature1.fill(s.cfg.inlet_temperature);
        cbs::EnergyAssembly::assembleStep4Rhs(s);

        bool passed = true;
        for (Int a = 1; a <= 4; ++a)
        {
            passed = passed && close_to(s.rhs1(a), 0.0);
        }
        print_result("FLUID UNIFORM REFERENCE TEMPERATURE", passed);
        return passed;
    }

    bool test_fluid_variational_integrals()
    {
        CBSStateSI s = reference_tetrahedron(0);
        set_field(s);
        const NodalScalar expected = independent_fluid_rhs(s);
        cbs::EnergyAssembly::assembleStep4Rhs(s);

        bool passed = true;
        for (Int a = 1; a <= 4; ++a)
        {
            const Real actual = s.rhs1(a);
            const Real target = expected[static_cast<std::size_t>(a - 1)];
            if (!close_to(actual, target))
            {
                std::cout << "  NODE " << a << " actual=" << actual
                          << " expected=" << target << "\n";
                passed = false;
            }
        }
        print_result("FLUID P1 QUADRATURE VS PRODUCTION RESIDUAL", passed);
        return passed;
    }

    bool test_solid_conduction()
    {
        CBSStateSI s = reference_tetrahedron(1);
        set_field(s);  // NONZERO VELOCITY MUST NOT ENTER SOLID RESIDUAL.
        cbs::EnergyAssembly::assembleStep4Rhs(s);

        Vector grad_t = {{0.0, 0.0, 0.0}};
        for (std::size_t a = 0; a < 4; ++a)
        {
            for (std::size_t d = 0; d < 3; ++d)
            {
                grad_t[d] += temperatures[a] * gradients[a][d];
            }
        }

        bool passed = true;
        Real residual_sum = 0.0;
        for (Int a = 1; a <= 4; ++a)
        {
            const Real expected = -s.k_e(1) / 6.0 *
                dot(gradients[static_cast<std::size_t>(a - 1)], grad_t);
            passed = passed && close_to(s.rhs1(a), expected);
            residual_sum += s.rhs1(a);
        }
        passed = passed && close_to(residual_sum, 0.0);
        print_result("SOLID CONDUCTION AND CONSERVATIVE ROW SUM", passed);
        return passed;
    }

    bool test_heat_source_and_neumann_flux()
    {
        CBSStateSI s = reference_tetrahedron(1, true);
        s.Qvol_e(1) = 24.0;
        s.cfg.heat_flux_bc = 12.0;

        s.iside(1, 1) = 1;
        s.iside(2, 1) = 2;
        s.iside(3, 1) = 3;
        s.iside(s.cfg.bsid, 1) = s.cfg.bc_noslip_heatflux_wall;
        s.face_norm(4, 1) = 0.5;

        cbs::EnergyAssembly::assembleStep4Rhs(s);

        // SOURCE: 24*(1/6)/4 = 1 W ON EVERY NODE.
        // FLUX: 12*(1/2)/3 = 2 W ON FACE NODES 1,2,3.
        bool passed = true;
        for (Int a = 1; a <= 4; ++a)
        {
            const Real expected = (a <= 3) ? 3.0 : 1.0;
            passed = passed && close_to(s.rhs1(a), expected);
        }
        print_result("VOLUMETRIC SOURCE + PRESCRIBED NEUMANN FLUX", passed);
        return passed;
    }

    bool test_element_matrix_matches_residual()
    {
        CBSStateSI s = reference_tetrahedron(0);
        set_field(s);
        Real matrix[5][5] = {{0.0}};
        cbs::EnergyAssembly::elementOperator(s, 1, matrix);
        cbs::EnergyAssembly::assembleStep4Rhs(s);

        NodalScalar original_rhs = {{
            s.rhs1(1), s.rhs1(2), s.rhs1(3), s.rhs1(4)
        }};

        bool passed = true;
        const Real epsilon = 1.0e-4;
        for (Int b = 1; b <= 4; ++b)
        {
            s.temperature1(b) += epsilon;
            cbs::EnergyAssembly::assembleStep4Rhs(s);
            for (Int a = 1; a <= 4; ++a)
            {
                const Real derivative =
                    (s.rhs1(a) - original_rhs[static_cast<std::size_t>(a - 1)]) /
                    epsilon;
                if (!close_to(derivative, matrix[a][b], 2.0e-6))
                {
                    std::cout << "  MATRIX " << a << "," << b
                              << " finite_difference=" << derivative
                              << " element_operator=" << matrix[a][b] << "\n";
                    passed = false;
                }
            }
            s.temperature1(b) -= epsilon;
        }
        print_result("ELEMENT MATRIX VS RESIDUAL DIRECTIONAL DERIVATIVE", passed);
        return passed;
    }
}

int main()
{
    int failures = 0;
    failures += test_uniform_temperature_invariance() ? 0 : 1;
    failures += test_fluid_variational_integrals() ? 0 : 1;
    failures += test_solid_conduction() ? 0 : 1;
    failures += test_heat_source_and_neumann_flux() ? 0 : 1;
    failures += test_element_matrix_matches_residual() ? 0 : 1;

    if (failures != 0)
    {
        std::cout << "STEP-4 OPERATOR VERIFICATION FAILED: "
                  << failures << " TEST(S)\n";
        return 1;
    }
    std::cout << "STEP-4 OPERATOR VERIFICATION PASS\n";
    return 0;
}
