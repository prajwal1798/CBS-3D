#pragma once

//=============================================================================
// CBS3D++_SI
//
// One-based storage container for a real square banded matrix.
//
// Let A be an n x n matrix with:
//
//     lower bandwidth = p
//     upper bandwidth = q
//
// A coefficient A(i,j) is inside the stored band when:
//
//     -p <= j - i <= q
//
// The matrix is stored in an Array2D with:
//
//     p + q + 1 rows
//     n columns
//
// The storage row corresponding to A(i,j) is:
//
//     band_row(i,j) = q + 1 + i - j
//
// Therefore:
//
//     main diagonal         band row q + 1
//     first superdiagonal  band row q
//     first subdiagonal    band row q + 2
//
// Example for p = 2 and q = 1:
//
//     stored row 1   A(1,2), A(2,3), A(3,4), ...
//     stored row 2   A(1,1), A(2,2), A(3,3), ...
//     stored row 3   A(2,1), A(3,2), A(4,3), ...
//     stored row 4   A(3,1), A(4,2), A(5,3), ...
//
// Matrix indices use the solver's one-based convention.
//
// This container is retained for the banded pressure-solver path. The current
// active pressure path uses either native Conjugate Gradient or PETSc.
//=============================================================================

#include "cbs/core/Array.hpp"
#include "cbs/core/Types.hpp"

namespace cbs
{
    class BandedMatrix
    {
    public:
        BandedMatrix() = default;

        // Removes all matrix storage and resets every dimension to zero.
        void clear();

        // Allocates an n x n banded matrix with the requested lower and upper
        // bandwidths. All stored coefficients are initialised to zero.
        void resize(
            Int n,
            Int lower_bw,
            Int upper_bw);

        // Assigns the same value to every stored band coefficient.
        void fill(Real value);

        Int size() const;
        Int lowerBandwidth() const;
        Int upperBandwidth() const;

        // Returns true when A(i,j) lies inside the matrix dimensions and inside
        // the stored band.
        bool inBand(
            Int i,
            Int j) const;

        // Provides checked read/write access to one stored coefficient.
        //
        // Access outside the matrix dimensions or outside the stored band
        // throws std::out_of_range.
        Real& operator()(
            Int i,
            Int j);

        const Real& operator()(
            Int i,
            Int j) const;

        // Adds value to A(i,j). Contributions outside the stored band are
        // ignored rather than treated as an error.
        void add(
            Int i,
            Int j,
            Real value);

        // Exposes the compact band storage when a linear solver needs direct
        // access to the underlying Array2D.
        const Array2D<Real>& data() const;
        Array2D<Real>& data();

    private:
        // Converts full-matrix indices to the compact band-storage row:
        //
        //     band_row = upper_bw + 1 + i - j
        Int bandRow(
            Int i,
            Int j) const;

        // Verifies full-matrix bounds and band membership.
        void checkIndex(
            Int i,
            Int j) const;

    private:
        Int n_ = 0;
        Int lower_bw_ = 0;
        Int upper_bw_ = 0;
        Int band_rows_ = 0;

        // Compact storage:
        //
        //     a_(band_row, column)
        Array2D<Real> a_;
    };
}
