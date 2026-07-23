#pragma once

//=============================================================================
// CBS3D++_SI
//
// Finite-element assembly of the CBS Step 3 pressure-gradient correction.
//
// For each fluid P1 tetrahedron, the element contribution is
//
//     r_{p,a}^{(e)} = - (V_e / 4) grad(p_h)|_e
//
// with
//
//     grad(p_h)|_e = sum_b p_b grad(N_b).
//
// The routine performs rank-local element assembly only. In distributed-memory
// execution, shared-node contributions must subsequently be reverse-summed to
// the unique node owner before the nodal velocity update is applied.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"

namespace cbs
{
    class VelocityCorrectionAssembly
    {
    public:
        // Assembles the rank-local Step 3 pressure-gradient vector into s.rhs.
        static void assembleStep3Rhs(CBSStateSI& s);
    };
}
