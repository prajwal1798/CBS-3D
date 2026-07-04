#pragma once

//=============================================================================
// CBS3D++_SI
//
// Momentum residual assembly for CBS Step 1.
//
// The semi-implicit momentum predictor is written in the nodal form:
//
//     M_L / dt * (u* - u^n) = r_m
//
// or:
//
//     u* = u^n + D_u^(-1) r_m
//
// where:
//
//     M_L        lumped momentum mass matrix
//     D_u^(-1)   inverse momentum mass/time diagonal stored in elcoe2
//     r_m        Step 1 momentum residual stored in rhs
//
// The assembled residual contains:
//
//     1. Galerkin convection
//     2. CBS characteristic correction
//     3. viscous diffusion
//
// Pressure correction is not performed here. It is applied later in CBS
// Step 3.
//
// Momentum is assembled only over fluid elements:
//
//     mat_elem(e) = 0
//
// The nodal velocity update itself is performed in Steps::step1SemiImplicit().
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"

namespace cbs
{
    class MomentumAssembly
    {
    public:
        // Assembles the complete CBS Step 1 momentum residual.
        static void assembleStep1Rhs(CBSStateSI& s);

        // Reserved hook for additional real-time momentum terms.
        static void applyRealTimeMomentumTerm(CBSStateSI& s);
    };
}
