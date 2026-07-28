#pragma once

//=============================================================================
// CBS3D++_SI
//
// Energy-equation residual assembly for three-dimensional conjugate heat
// transfer.
//
// The dimensional energy equation is:
//
//     rho cp [ dT/dt + u . grad(T) ]
//         = div(k grad(T)) + Q
//
// The semi-implicit nodal update is:
//
//     T^(n+1) = T^n + D_T^(-1) r_T
//
// where:
//
//     D_T^(-1) = elcoe2p
//
// is the inverse lumped thermal-capacitance/time diagonal and r_T is assembled
// into rhs1.
//
// Element treatment:
//
//     Fluid element:
//         convection
//         CBS/SUPG-style convection stabilisation
//         thermal diffusion
//         volumetric heat source
//
//     Solid element:
//         thermal diffusion
//         volumetric heat source
//
// Boundary treatment:
//
//     BC 532:
//         prescribed external heat flux
//
//     BC 901:
//         conformal fluid-solid interface
//
// No separate interface RHS is required for BC 901 because fluid and solid
// tetrahedra share the same interface nodes. Temperature continuity is
// therefore imposed naturally by the conformal finite-element mesh.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"

namespace cbs
{
    class EnergyAssembly
    {
    public:
        // Assembles the complete CBS Step 4 thermal residual.
        static void assembleStep4Rhs(CBSStateSI& s);

        // Builds the four-node element transport matrix A_e defined by:
        //
        //     r_e(T) = -A_e T_e + b_e,
        //
        // where b_e contains volumetric and boundary loads and is deliberately
        // excluded. The matrix contains exactly the same convection,
        // characteristic and diffusion coefficients as assembleStep4Rhs().
        // Positions 1..4 are active; row/column zero are unused.
        static void buildElementTransportMatrix(
            const CBSStateSI& s,
            Int ie,
            Real matrix[5][5]);

        // Reserved hook for future real-time or BDF thermal-history terms.
        static void applyRealTimeEnergyTerm(CBSStateSI& s);
    };
}
