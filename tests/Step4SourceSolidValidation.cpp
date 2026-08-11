//=============================================================================
// CBS3D++_SI
//
// Production validation for the solid volumetric-source input contract.
//=============================================================================

#include "cbs/assembly/EnergyAssembly.hpp"
#include "cbs/core/CBSStateSI.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{
    using cbs::CBSStateSI;
    using cbs::Int;
    using cbs::Real;

    constexpr Real tolerance = 1.0e-12;

    Int gradient_index(const Int dim, const Int local_node)
    {
        return (dim - 1) * 4 + local_node;
    }

    CBSStateSI make_reference_tetrahedron(const Int material_id)
    {
        CBSStateSI s;
        s.set_problem_sizes(1, 4, 0, 0);

        s.cfg.temp_calc = 1;
        s.cfg.source_solid = 0.0;
        s.cfg.heat_flux_bc = 0.0;
        s.cfg.turbulence_on = 0;
        s.cfg.turbulent_thermal_diffusivity_on = 0;

        for (Int a = 1; a <= 4; ++a)
        {
            s.intma(a, 1) = a;
        }

        const std::array<Real, 4> dNdx = {-1.0, 1.0, 0.0, 0.0};
        const std::array<Real, 4> dNdy = {-1.0, 0.0, 1.0, 0.0};
        const std::array<Real, 4> dNdz = {-1.0, 0.0, 0.0, 1.0};

        for (Int a = 1; a <= 4; ++a)
        {
            s.dNkdx(gradient_index(1, a)) = dNdx[static_cast<std::size_t>(a - 1)];
            s.dNkdx(gradient_index(2, a)) = dNdy[static_cast<std::size_t>(a - 1)];
            s.dNkdx(gradient_index(3, a)) = dNdz[static_cast<std::size_t>(a - 1)];
        }

        s.detJ(1) = 1.0;
        s.delte(1) = 1.0e-2;
        s.mat_elem(1) = material_id;

        s.rho_e(1) = 1.0;
        s.cp_e(1) = 1.0;
        s.rho_cp_e(1) = 1.0;
        s.k_e(1) = 1.0;
        s.k_eff_e(1) = 1.0;
        s.Qvol_e(1) = 0.0;

        s.temperature.fill(300.0);
        s.temperature1.fill(300.0);
        s.unkno.fill(0.0);
        s.rhs1.fill(0.0);

        return s;
    }

    bool check_rhs(
        const std::string& name,
        CBSStateSI& s,
        const Real expected)
    {
        cbs::EnergyAssembly::assembleStep4Rhs(s);

        bool passed = true;

        for (Int a = 1; a <= 4; ++a)
        {
            passed = passed &&
                std::abs(s.rhs1(a) - expected) <= tolerance;
        }

        std::cout << name << ": " << (passed ? "PASS" : "FAIL") << "\n";
        return passed;
    }

    bool test_par_source_on_solid()
    {
        CBSStateSI s = make_reference_tetrahedron(1);
        s.cfg.source_solid = 24.0;

        // detJ=1 => V=1/6 and Q*V/4 = Q*detJ/24 = 1.
        return check_rhs(".par source_solid on solid", s, 1.0);
    }

    bool test_matprop_fallback_on_solid()
    {
        CBSStateSI s = make_reference_tetrahedron(1);
        s.Qvol_e(1) = 24.0;

        return check_rhs(".matprop Qvol fallback on solid", s, 1.0);
    }

    bool test_par_source_not_applied_to_fluid()
    {
        CBSStateSI s = make_reference_tetrahedron(0);
        s.cfg.source_solid = 24.0;

        return check_rhs(".par source_solid ignored by fluid", s, 0.0);
    }

    bool test_double_specification_rejected()
    {
        CBSStateSI s = make_reference_tetrahedron(1);
        s.cfg.source_solid = 24.0;
        s.Qvol_e(1) = 12.0;

        try
        {
            cbs::EnergyAssembly::assembleStep4Rhs(s);
        }
        catch (const std::runtime_error& error)
        {
            const std::string message = error.what();
            const bool passed =
                message.find("specified in both .par source_solid and .matprop Qvol")
                != std::string::npos;

            std::cout
                << "double source specification rejected: "
                << (passed ? "PASS" : "FAIL") << "\n";

            return passed;
        }

        std::cout << "double source specification rejected: FAIL\n";
        return false;
    }
}

int main()
{
    int failures = 0;

    failures += test_par_source_on_solid() ? 0 : 1;
    failures += test_matprop_fallback_on_solid() ? 0 : 1;
    failures += test_par_source_not_applied_to_fluid() ? 0 : 1;
    failures += test_double_specification_rejected() ? 0 : 1;

    if (failures == 0)
    {
        std::cout << "SOLID VOLUMETRIC SOURCE VALIDATION: ALL TESTS PASSED\n";
        return 0;
    }

    std::cout
        << "SOLID VOLUMETRIC SOURCE VALIDATION: "
        << failures << " TEST(S) FAILED\n";

    return 1;
}
