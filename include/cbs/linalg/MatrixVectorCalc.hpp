#pragma once

//=============================================================================
// CBS3D++_SI
//
// Matrix-vector products for the compact tetrahedral pressure operator.
//
// The pressure matrix is not stored as a global dense or CSR matrix. Instead,
// its coefficients are kept in the original finite-element form:
//
//     pdiag(i)   assembled diagonal coefficient A_ii
//
//     gstif(6(e-1)+1) = A_12^(e)
//     gstif(6(e-1)+2) = A_13^(e)
//     gstif(6(e-1)+3) = A_14^(e)
//     gstif(6(e-1)+4) = A_23^(e)
//     gstif(6(e-1)+5) = A_24^(e)
//     gstif(6(e-1)+6) = A_34^(e)
//
// For one four-node tetrahedron:
//
//           [ A11  A12  A13  A14 ]
//     A_e = [ A12  A22  A23  A24 ]
//           [ A13  A23  A33  A34 ]
//           [ A14  A24  A34  A44 ]
//
// The off-diagonal terms are applied element by element. The assembled nodal
// diagonal is applied afterwards through pdiag.
//
// Pressure exists only in the fluid domain:
//
//     mat_elem(e) = 0
//
// Prescribed-pressure rows are treated as identity equations:
//
//     y_i = x_i
//
// This matrix-free form is used by the native Pressure CG solver.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"

namespace cbs
{
    class MatrixVectorCalc
    {
    public:
        // Calculates:
        //
        //     y = A x
        //
        // for a scalar pressure vector.
        static void pressureMultiply(
            const CBSStateSI& s,
            const Array1D<Real>& x,
            Array1D<Real>& y);

        // Applies the same pressure operator independently to every component
        // of a nodal vector field.
        static void pressureMultiplyVector(
            const CBSStateSI& s,
            const Array2D<Real>& x,
            Array2D<Real>& y);
    };
}
