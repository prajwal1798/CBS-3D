#pragma once

//=============================================================================
// CBS3D++_SI
//
// Compressed Sparse Row storage for a real rectangular matrix.
//
// A sparse matrix A with nrows rows and ncols columns is represented by:
//
//     row_ptr[i-1] <= k < row_ptr[i]
//
// for all stored coefficients in matrix row i. The corresponding column and
// value are:
//
//     j      = col_ind[k]
//     A(i,j) = values[k]
//
// Indexing conventions:
//
//     matrix rows and columns     one-based
//     Array1D vectors             one-based
//     std::vector storage index   zero-based
//
// Therefore row_ptr has nrows + 1 entries:
//
//     row_ptr[0]       = 0
//     row_ptr[nrows]   = number of stored coefficients
//
// Triplet construction accepts entries:
//
//     (row[k], column[k], value[k])
//
// The entries are sorted by row and then column. Repeated coordinates are
// combined by addition:
//
//     A(i,j) = sum of all supplied values at coordinate (i,j)
//
// Matrix-vector multiplication evaluates:
//
//     y_i = sum_j A(i,j) x_j
//
// This class is currently retained as a general sparse-matrix utility. The
// active CBS pressure operator uses its specialised compact tetrahedral storage
// or PETSc matrix assembly.
//=============================================================================

#include "cbs/core/Array.hpp"
#include "cbs/core/Types.hpp"

#include <vector>

namespace cbs
{
    class CSRMatrix
    {
    public:
        CSRMatrix() = default;

        // Removes the dimensions and all sparse storage.
        void clear();

        // Sets the matrix dimensions and creates an empty CSR row-pointer
        // structure. No coefficients are stored.
        void setSize(
            Int nrows,
            Int ncols);

        // Builds the complete CSR matrix from one-based coordinate triplets.
        // Duplicate coordinates are summed.
        void setFromTriplets(
            Int nrows,
            Int ncols,
            const std::vector<Int>& rows,
            const std::vector<Int>& cols,
            const std::vector<Real>& vals);

        Int nrows() const;
        Int ncols() const;
        Int nnz() const;

        // Read-only access to the CSR arrays.
        const std::vector<Int>& rowPtr() const;
        const std::vector<Int>& colInd() const;
        const std::vector<Real>& values() const;

        // Writable access is retained for external sparse-matrix assembly.
        // The caller is responsible for preserving all CSR invariants.
        std::vector<Int>& rowPtr();
        std::vector<Int>& colInd();
        std::vector<Real>& values();

        // Calculates y = A x.
        void matVec(
            const Array1D<Real>& x,
            Array1D<Real>& y) const;

    private:
        Int nrows_ = 0;
        Int ncols_ = 0;

        // CSR row offsets. Size is nrows + 1.
        std::vector<Int> row_ptr_;

        // One-based matrix column number for each stored coefficient.
        std::vector<Int> col_ind_;

        // Numerical value of each stored coefficient.
        std::vector<Real> values_;
    };
}
