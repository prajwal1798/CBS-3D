//=============================================================================
// CBS3D++_SI
// Thermal lumped-mass/time-diagonal verification, including a conformal
// solid-fluid interface. Calls the production TimeStep code unmodified.
//
// Independent P1 tetrahedral reference:
//   M_L,a^(e) = rho_cp_e * V_e/4
//   H_i = sum_(e incident to i) M_L,i^(e) / dt_e
//   elcoe2p_i = 1 / H_i
// Momentum: fluid-only M_L/step; no solid contribution.
//=============================================================================

#include "cbs/timestep/TimeStep.hpp"
#include "cbs/core/CBSStateSI.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace
{
    using cbs::CBSStateSI;
    using cbs::Int;
    using cbs::Real;

    bool approximately(const Real value, const Real reference)
    {
        return std::abs(value - reference) <=
            1.0e-12 * std::max(1.0, std::abs(reference));
    }

    CBSStateSI two_material_reference()
    {
        CBSStateSI s;
        // Two TET4 elements sharing node 1. Other nodes are disjoint:
        //  element 1: fluid nodes 1,2,3,4
        //  element 2: solid nodes 1,5,6,7
        s.set_problem_sizes(2, 7, 0, 0);
        s.cfg.bc_fixed = 0;
        for (Int a = 1; a <= 4; ++a)
        {
            s.intma(a, 1) = a;
            s.intma(a, 2) = a == 1 ? 1 : a + 3;
        }

        s.mat_elem(1) = 0;
        s.mat_elem(2) = 1;

        // Tetrahedral volumes 1 and 2; lumped element nodal masses
        // are V/4 = 0.25 and 0.5 respectively.
        s.detJ(1) = 6.0;
        s.detJ(2) = 12.0;
        s.delte(1) = 0.1;
        s.delte(2) = 0.2;
        s.rho_cp_e(1) = 10.0;
        s.rho_cp_e(2) = 100.0;

        for (Int a = 1; a <= 4; ++a)
        {
            s.elcoe_e(a) = 0.25;
            s.elcoe_e(4 + a) = 0.5;
            // Nonzero pressure diagonal for all fluid-connected nodes.
            s.pdiagE(a) = 2.0;
            s.pdiagE(4 + a) = 0.0;
        }

        s.gstifE.fill(0.0);
        s.pdiag.fill(0.0);
        s.rhs1.fill(0.0);
        s.pres.fill(0.0);
        return s;
    }
}

int main()
{
    CBSStateSI s = two_material_reference();
    cbs::TimeStep::updateLhsDiagonal(s);

    // At node 1, incident heat capacities add across the two materials:
    //  10*(1/4)/0.1 + 100*(2/4)/0.2 = 25 + 250 = 275.
    bool ok = approximately(s.elcoe2p(1), 1.0 / 275.0);

    // Fluid-only nodes 2..4: 10*(1/4)/0.1 = 25.
    for (Int n = 2; n <= 4; ++n)
    {
        ok = ok && approximately(s.elcoe2p(n), 1.0 / 25.0);
    }

    // Solid-only nodes 5..7: 100*(2/4)/0.2 = 250.
    for (Int n = 5; n <= 7; ++n)
    {
        ok = ok && approximately(s.elcoe2p(n), 1.0 / 250.0);
    }

    std::cout << "MIXED-MATERIAL THERMAL DIAGONAL: "
              << (ok ? "PASS" : "FAIL") << "\n";

    // Fluid momentum capacity/mass is assembled from fluid elements only.
    // Fluid lumped diagonal 0.25/0.1 = 2.5, inverse 0.4.
    bool momentum_ok = true;
    for (Int n = 1; n <= 4; ++n)
    {
        momentum_ok = momentum_ok && approximately(s.elcoe2(n), 0.4);
    }
    for (Int n = 5; n <= 7; ++n)
    {
        momentum_ok = momentum_ok && approximately(s.elcoe2(n), 0.0);
    }

    std::cout << "FLUID-ONLY MOMENTUM DIAGONAL: "
              << (momentum_ok ? "PASS" : "FAIL") << "\n";

    // Pressure matrix is scaled by fluid element dt only; solid pressure
    // diagonals must remain inactive.
    bool pressure_ok = true;
    for (Int n = 1; n <= 4; ++n)
    {
        pressure_ok = pressure_ok && approximately(s.pdiag(n), 0.2);
    }
    for (Int n = 5; n <= 7; ++n)
    {
        pressure_ok = pressure_ok && approximately(s.pdiag(n), 0.0);
    }

    std::cout << "PRESSURE FLUID/SOLID MASK: "
              << (pressure_ok ? "PASS" : "FAIL") << "\n";

    if (!(ok && momentum_ok && pressure_ok))
    {
        std::cerr << "LUMPED MASS VERIFICATION FAILED\n";
        return 1;
    }
    std::cout << "LUMPED MASS VERIFICATION PASS\n";
    return 0;
}
