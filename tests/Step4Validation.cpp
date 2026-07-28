//=============================================================================
// CBS3D++_SI
//
// Standalone verification of the production Step-4 thermal residual assembly
// on one reference P1 tetrahedron.
//
// The test calls EnergyAssembly::assembleStep4Rhs() directly.  It therefore
// checks the production implementation rather than a duplicated test formula.
//=============================================================================

#include "cbs/assembly/EnergyAssembly.hpp"
#include "cbs/core/CBSStateSI.hpp"

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

namespace
{
    using cbs::CBSStateSI;
    using cbs::Int;
    using cbs::Real;

    constexpr Real absolute_tolerance = 1.0e-12;
    constexpr Real relative_tolerance = 1.0e-12;

    Int gradient_index(const Int dim, const Int local_node)
    {
        // One element only.  EnergyAssembly stores the gradients in dim-major
        // order: dN1/dx..dN4/dx, dN1/dy..dN4/dy, dN1/dz..dN4/dz.
        return (dim - 1) * 4 + local_node;
    }

    bool nearly_equal(const Real actual, const Real expected)
    {
        const Real scale = 1.0 + std::abs(expected);
        return std::abs(actual - expected)
            <= absolute_tolerance + relative_tolerance * scale;
    }

    CBSStateSI make_reference_tetrahedron(const Int nboun = 0)
    {
        CBSStateSI s;
        s.set_problem_sizes(1, 4, nboun, 0);

        s.cfg.temp_calc = 1;
        s.cfg.heat_flux_bc = 0.0;
        s.cfg.turbulence_on = 0;
        s.cfg.turbulent_thermal_diffusivity_on = 0;

        // Reference tetrahedron:
        //
        //   node 1 = (0,0,0)
        //   node 2 = (1,0,0)
        //   node 3 = (0,1,0)
        //   node 4 = (0,0,1)
        //
        // det(J)=1 and V=1/6.
        s.coord(1, 1) = 0.0; s.coord(2, 1) = 0.0; s.coord(3, 1) = 0.0;
        s.coord(1, 2) = 1.0; s.coord(2, 2) = 0.0; s.coord(3, 2) = 0.0;
        s.coord(1, 3) = 0.0; s.coord(2, 3) = 1.0; s.coord(3, 3) = 0.0;
        s.coord(1, 4) = 0.0; s.coord(2, 4) = 0.0; s.coord(3, 4) = 1.0;

        for (Int a = 1; a <= 4; ++a)
        {
            s.intma(a, 1) = a;
        }

        // N1=1-x-y-z, N2=x, N3=y, N4=z.
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
        s.mat_elem(1) = 0;

        s.rho_e(1) = 1.0;
        s.cp_e(1) = 3.0;
        s.rho_cp_e(1) = 3.0;
        s.k_e(1) = 5.0;
        s.Qvol_e(1) = 0.0;

        s.temperature.fill(300.0);
        s.temperature1.fill(300.0);
        s.unkno.fill(0.0);
        s.rhs1.fill(0.0);

        return s;
    }

    void set_linear_temperature_x(CBSStateSI& s)
    {
        // T=300+x, hence grad(T)=(1,0,0) exactly on this tetrahedron.
        for (Int ip = 1; ip <= 4; ++ip)
        {
            const Real temperature = 300.0 + s.coord(1, ip);
            s.temperature(ip) = temperature;
            s.temperature1(ip) = temperature;
        }
    }

    bool check_vector(
        const std::string& test_name,
        const std::array<Real, 4>& actual,
        const std::array<Real, 4>& expected)
    {
        bool passed = true;

        std::cout << "\n" << test_name << "\n";
        std::cout << "node              actual              expected             error\n";

        for (Int a = 1; a <= 4; ++a)
        {
            const Real value = actual[static_cast<std::size_t>(a - 1)];
            const Real target = expected[static_cast<std::size_t>(a - 1)];
            const Real error = value - target;
            const bool node_passed = nearly_equal(value, target);
            passed = passed && node_passed;

            std::cout
                << std::setw(4) << a
                << std::setw(20) << std::setprecision(12) << std::scientific << value
                << std::setw(22) << target
                << std::setw(20) << error
                << (node_passed ? "  PASS" : "  FAIL")
                << "\n";
        }

        std::cout << (passed ? "RESULT: PASS\n" : "RESULT: FAIL\n");
        return passed;
    }

    std::array<Real, 4> assembled_rhs(CBSStateSI& s)
    {
        cbs::EnergyAssembly::assembleStep4Rhs(s);

        return {
            s.rhs1(1),
            s.rhs1(2),
            s.rhs1(3),
            s.rhs1(4)
        };
    }

    bool test_uniform_temperature_zero_residual()
    {
        CBSStateSI s = make_reference_tetrahedron();

        for (Int ip = 1; ip <= 4; ++ip)
        {
            s.unkno(1, ip) = 2.0;
            s.unkno(2, ip) = -1.0;
            s.unkno(3, ip) = 0.5;
        }

        return check_vector(
            "TEST 1: uniform temperature with nonzero velocity",
            assembled_rhs(s),
            {0.0, 0.0, 0.0, 0.0});
    }

    bool test_linear_diffusion_element_vector()
    {
        CBSStateSI s = make_reference_tetrahedron();
        set_linear_temperature_x(s);

        // r_diff,a = -k V grad(N_a).grad(T), k=5, V=1/6.
        return check_vector(
            "TEST 2: exact linear-temperature diffusion vector",
            assembled_rhs(s),
            {5.0 / 6.0, -5.0 / 6.0, 0.0, 0.0});
    }

    bool test_convection_and_characteristic_increment()
    {
        CBSStateSI stationary = make_reference_tetrahedron();
        CBSStateSI moving = make_reference_tetrahedron();

        set_linear_temperature_x(stationary);
        set_linear_temperature_x(moving);

        for (Int ip = 1; ip <= 4; ++ip)
        {
            moving.unkno(1, ip) = 2.0;
        }

        const std::array<Real, 4> stationary_rhs = assembled_rhs(stationary);
        const std::array<Real, 4> moving_rhs = assembled_rhs(moving);

        std::array<Real, 4> velocity_dependent_increment{};
        for (std::size_t i = 0; i < velocity_dependent_increment.size(); ++i)
        {
            // Subtraction cancels the identical diffusion contribution exactly,
            // isolating convection plus the characteristic correction while
            // still calling the unmodified production assembly routine.
            velocity_dependent_increment[i] = moving_rhs[i] - stationary_rhs[i];
        }

        // Exact values for rhoCp=3, detJ=1, dt=0.01, u=(2,0,0), gradT=(1,0,0):
        //
        // Galerkin convection = -0.25 at every node.
        // Characteristic term = {-0.01,+0.01,0,0}.
        return check_vector(
            "TEST 3: exact convection plus characteristic increment",
            velocity_dependent_increment,
            {-0.26, -0.24, -0.25, -0.25});
    }

    bool test_uniform_volumetric_source()
    {
        CBSStateSI s = make_reference_tetrahedron();
        s.Qvol_e(1) = 24.0;

        // Q V/4 = Q detJ/24 = 1 at every element node.
        return check_vector(
            "TEST 4: exact uniform volumetric-source vector",
            assembled_rhs(s),
            {1.0, 1.0, 1.0, 1.0});
    }

    bool test_prescribed_heat_flux_face()
    {
        CBSStateSI s = make_reference_tetrahedron(1);

        s.cfg.heat_flux_bc = 12.0;
        s.iside(1, 1) = 1;
        s.iside(2, 1) = 2;
        s.iside(3, 1) = 3;
        s.iside(s.cfg.bsid, 1) = s.cfg.bc_noslip_heatflux_wall;
        s.face_norm(4, 1) = 0.5;

        // q'' A/3 = 12*0.5/3 = 2 at each node of the triangular face.
        return check_vector(
            "TEST 5: exact prescribed heat-flux face vector",
            assembled_rhs(s),
            {2.0, 2.0, 2.0, 0.0});
    }
}

int main()
{
    std::cout << "CBS3D++_SI Step-4 element validation\n";
    std::cout << "Production routine: EnergyAssembly::assembleStep4Rhs\n";

    int failures = 0;

    failures += test_uniform_temperature_zero_residual() ? 0 : 1;
    failures += test_linear_diffusion_element_vector() ? 0 : 1;
    failures += test_convection_and_characteristic_increment() ? 0 : 1;
    failures += test_uniform_volumetric_source() ? 0 : 1;
    failures += test_prescribed_heat_flux_face() ? 0 : 1;

    std::cout << "\n============================================================\n";
    if (failures == 0)
    {
        std::cout << "STEP-4 ELEMENT VALIDATION: ALL TESTS PASSED\n";
        return 0;
    }

    std::cout << "STEP-4 ELEMENT VALIDATION: " << failures
              << " TEST(S) FAILED\n";
    return 1;
}
