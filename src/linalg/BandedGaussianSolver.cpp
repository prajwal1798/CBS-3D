//=============================================================================
// CBS3D++_SI
//
// Direct solution of a banded linear system:
//
//     A x = b
//
// Forward elimination:
//
//     m_ik = A(i,k) / A(k,k)
//
//     A(i,j) <- A(i,j) - m_ik A(k,j)
//
//     b_i    <- b_i    - m_ik b_k
//
// Back substitution:
//
//     x_i =
//         [ b_i - sum_j A(i,j)x_j ]
//         / A(i,i)
//
// Only coefficients inside the stored matrix band are visited. No pivoting is
// performed, so the method requires acceptable diagonal pivots in the supplied
// matrix ordering.
//=============================================================================

#include "cbs/linalg/BandedGaussianSolver.hpp"

#include <cmath>
#include <stdexcept>

namespace cbs
{
    //=========================================================================
    // Solves the banded system A x = b.
    //
    // Inputs:
    //     A   square banded coefficient matrix
    //     b   right-hand side
    //
    // Output:
    //     x   solution vector
    //
    // Modified data:
    //     A   overwritten by its upper-triangular elimination form
    //
    // Preserved data:
    //     b   copied to the local vector rhs before elimination
    //
    // The algorithm does not perform partial or complete pivoting.
    //=========================================================================
    void BandedGaussianSolver::solve(BandedMatrix& A, Array1D<Real>& b, Array1D<Real>& x)
    {
        const Int n = A.size();

        if (b.size() != n)
        {
            throw std::runtime_error("BandedGaussianSolver::solve - RHS size mismatch");
        }

        const Int ml = A.lowerBandwidth();
        const Int mu = A.upperBandwidth();

        // Preserve the caller's RHS while applying elimination to a local copy.
        Array1D<Real> rhs(n);
        for (Int i = 1; i <= n; ++i)
        {
            rhs(i) = b(i);
        }

        // --------------------------------------------------------------------
        // Forward elimination
        //
        // At pivot column k, only rows:
        //
        //     k + 1 <= i <= min(n, k + ml)
        //
        // can contain a stored entry A(i,k).
        // --------------------------------------------------------------------
        for (Int k = 1; k <= n - 1; ++k)
        {
            const Real akk = A(k, k);

            if (std::abs(akk) <= 1.0e-30)
            {
                throw std::runtime_error("BandedGaussianSolver::solve - near-zero pivot detected");
            }

            const Int i_max = (k + ml <= n) ? (k + ml) : n;

            for (Int i = k + 1; i <= i_max; ++i)
            {
                if (!A.inBand(i, k))
                {
                    continue;
                }

                // Elimination multiplier:
                //
                //     m_ik = A(i,k) / A(k,k)
                const Real factor = A(i, k) / akk;
                A(i, k) = 0.0;

                // Only columns in the pivot row's upper band can contribute:
                //
                //     k + 1 <= j <= min(n, k + mu)
                const Int j_max = (k + mu <= n) ? (k + mu) : n;
                for (Int j = k + 1; j <= j_max; ++j)
                {
                    if (A.inBand(i, j) && A.inBand(k, j))
                    {
                        A(i, j) -= factor * A(k, j);
                    }
                }

                rhs(i) -= factor * rhs(k);
            }
        }

        // --------------------------------------------------------------------
        // Back substitution
        // --------------------------------------------------------------------
        x.resize(n);
        x.fill(0.0);

        if (std::abs(A(n, n)) <= 1.0e-30)
        {
            throw std::runtime_error("BandedGaussianSolver::solve - near-zero pivot detected at back substitution");
        }

        x(n) = rhs(n) / A(n, n);

        for (Int i = n - 1; i >= 1; --i)
        {
            Real sum = rhs(i);

            const Int j_max = (i + mu <= n) ? (i + mu) : n;
            for (Int j = i + 1; j <= j_max; ++j)
            {
                if (A.inBand(i, j))
                {
                    sum -= A(i, j) * x(j);
                }
            }

            if (std::abs(A(i, i)) <= 1.0e-30)
            {
                throw std::runtime_error("BandedGaussianSolver::solve - near-zero diagonal in back substitution");
            }

            x(i) = sum / A(i, i);
        }
    }
}
