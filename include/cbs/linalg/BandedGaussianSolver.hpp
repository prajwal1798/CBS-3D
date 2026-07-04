#pragma once

//=============================================================================
// CBS3D++_SI
//
// Direct linear solver for a real square banded system:
//
//     A x = b
//
// The algorithm performs:
//
//     1. Gaussian forward elimination inside the stored lower and upper bands.
//     2. Back substitution inside the resulting upper band.
//
// For elimination column k:
//
//     m_ik = A(i,k) / A(k,k)
//
//     A(i,j) <- A(i,j) - m_ik A(k,j)
//
//     b_i    <- b_i    - m_ik b_k
//
// Back substitution uses:
//
//     x_i =
//         [ b_i - sum_(j=i+1)^(min(n,i+q)) A(i,j)x_j ]
//         / A(i,i)
//
// where q is the upper bandwidth.
//
// Important properties:
//
//     - No row pivoting is performed.
//     - A is overwritten by the elimination process.
//     - b is preserved because the routine works on an internal RHS copy.
//     - x is resized and overwritten with the solution.
//     - Pivots with magnitude <= 1.0e-30 are treated as singular.
//
// This module is retained for the inactive banded pressure-solver path. The
// active CBS3D++_SI pressure solve currently uses native Conjugate Gradient or
// PETSc CG with AMG.
//=============================================================================

#include "cbs/core/Array.hpp"
#include "cbs/core/Types.hpp"
#include "cbs/linalg/BandedMatrix.hpp"

namespace cbs
{
    class BandedGaussianSolver
    {
    public:
        // Solves A x = b by in-place banded Gaussian elimination followed by
        // back substitution.
        static void solve(
            BandedMatrix& A,
            Array1D<Real>& b,
            Array1D<Real>& x);
    };
}
