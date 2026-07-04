#pragma once

//=============================================================================
// CBS3D++_SI
//
// Pressure assembly for the semi-implicit three-dimensional CBS solver.
//
// The pressure-correction equation is assembled only over fluid elements:
//
//     mat_elem(e) = 0
//
// For a four-node linear tetrahedral element, the pressure stiffness matrix is
//
//     H_ij^(e) = integral(V_e) grad(N_i) . grad(N_j) dV
//
// Because grad(N_i) is constant inside a P1 tetrahedron,
//
//     H_ij^(e) = V_e grad(N_i) . grad(N_j)
//
// The element matrix is stored in compact form:
//
//     pdiagE : four diagonal entries
//     gstifE : six upper-triangular off-diagonal entries
//
// Off-diagonal ordering:
//
//     1 -> (1,2)
//     2 -> (1,3)
//     3 -> (1,4)
//     4 -> (2,3)
//     5 -> (2,4)
//     6 -> (3,4)
//
// The time-step-scaled global pressure operator is completed later by
// TimeStep::updateLhsDiagonal().
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"

namespace cbs
{
    class PressureAssembly
    {
    public:
        // Calculates the compact element pressure-stiffness coefficients.
        static void buildElementPressureTerms(CBSStateSI& s);

        // Initialises the global pressure arrays before time-step scaling.
        static void buildGlobalPressureTerms(CBSStateSI& s);

        // Assembles the weak-divergence right-hand side for CBS Step 2.
        static void assembleStep2Rhs(CBSStateSI& s);

        // Compatibility wrapper for the Step 2 pressure RHS assembly.
        static void assembleRhs(CBSStateSI& s);
    };
}
