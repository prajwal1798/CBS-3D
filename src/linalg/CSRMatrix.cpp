//=============================================================================
// CBS3D++_SI
//
// Compressed Sparse Row matrix construction and multiplication.
//
// For row i, stored coefficients occupy:
//
//     row_ptr[i-1] <= k < row_ptr[i]
//
// and the matrix-vector product is:
//
//     y_i = sum_k values[k] x(col_ind[k])
//
// Matrix row and column numbers remain one-based, while the underlying
// std::vector positions are zero-based.
//=============================================================================

#include "cbs/linalg/CSRMatrix.hpp"

#include <algorithm>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace cbs
{
    //=========================================================================
    // Removes all dimensions and sparse coefficients.
    //=========================================================================
    void CSRMatrix::clear()
    {
        nrows_ = 0;
        ncols_ = 0;
        row_ptr_.clear();
        col_ind_.clear();
        values_.clear();
    }


    //=========================================================================
    // Creates an empty nrows x ncols CSR matrix.
    //
    // The row-pointer array is initialised as:
    //
    //     row_ptr[0:nrows] = 0
    //
    // and no column indices or numerical values are stored.
    //=========================================================================
    void CSRMatrix::setSize(const Int nrows, const Int ncols)
    {
        if (nrows < 0 || ncols < 0)
        {
            throw std::runtime_error("CSRMatrix::setSize - negative dimension");
        }

        nrows_ = nrows;
        ncols_ = ncols;

        row_ptr_.assign(static_cast<Size>(nrows_) + 1U, 0);
        col_ind_.clear();
        values_.clear();
    }


    //=========================================================================
    // Builds CSR storage from one-based coordinate triplets.
    //
    // Input arrays:
    //
    //     rows[k]   matrix row
    //     cols[k]   matrix column
    //     vals[k]   coefficient value
    //
    // Construction sequence:
    //
    //     1. Validate dimensions and triplet indices.
    //     2. Sort entries by row and then column.
    //     3. Sum repeated coordinates.
    //     4. Count entries in every row.
    //     5. Prefix-sum the row counts to obtain row_ptr.
    //     6. Copy the sorted columns and values into CSR storage.
    //
    // Repeated coordinates satisfy:
    //
    //     A(i,j) = sum_m value_m(i,j)
    //=========================================================================
    void CSRMatrix::setFromTriplets(const Int nrows,
        const Int ncols,
        const std::vector<Int>& rows,
        const std::vector<Int>& cols,
        const std::vector<Real>& vals)
    {
        if (rows.size() != cols.size() || rows.size() != vals.size())
        {
            throw std::runtime_error("CSRMatrix::setFromTriplets - inconsistent triplet sizes");
        }

        if (nrows < 0 || ncols < 0)
        {
            throw std::runtime_error("CSRMatrix::setFromTriplets - negative dimension");
        }

        nrows_ = nrows;
        ncols_ = ncols;

        std::vector<std::tuple<Int, Int, Real>> triplets;
        triplets.reserve(rows.size());

        // Validate and collect one-based coordinate entries.
        for (Size k = 0; k < rows.size(); ++k)
        {
            const Int r = rows[k];
            const Int c = cols[k];
            const Real v = vals[k];

            if (r < 1 || r > nrows_ || c < 1 || c > ncols_)
            {
                throw std::runtime_error("CSRMatrix::setFromTriplets - triplet index out of bounds");
            }

            triplets.emplace_back(r, c, v);
        }

        // CSR requires entries to be grouped by row. Sorting by column inside
        // each row also makes duplicate coordinates adjacent.
        std::sort(triplets.begin(), triplets.end(),
            [](const auto& a, const auto& b)
            {
                if (std::get<0>(a) != std::get<0>(b))
                {
                    return std::get<0>(a) < std::get<0>(b);
                }
                return std::get<1>(a) < std::get<1>(b);
            });

        std::vector<Int> rows_u;
        std::vector<Int> cols_u;
        std::vector<Real> vals_u;

        rows_u.reserve(triplets.size());
        cols_u.reserve(triplets.size());
        vals_u.reserve(triplets.size());

        // Merge repeated matrix coordinates by summing their values.
        for (const auto& t : triplets)
        {
            const Int r = std::get<0>(t);
            const Int c = std::get<1>(t);
            const Real v = std::get<2>(t);

            if (!rows_u.empty() && rows_u.back() == r && cols_u.back() == c)
            {
                vals_u.back() += v;
            }
            else
            {
                rows_u.push_back(r);
                cols_u.push_back(c);
                vals_u.push_back(v);
            }
        }

        row_ptr_.assign(static_cast<Size>(nrows_) + 1U, 0);

        // Count the number of stored coefficients in every matrix row.
        //
        // Count for one-based row r is temporarily stored at row_ptr[r].
        for (const Int r : rows_u)
        {
            ++row_ptr_[static_cast<Size>(r)];
        }

        // Prefix sum:
        //
        //     row_ptr[i] =
        //         total number of stored coefficients in rows 1 ... i
        for (Int i = 1; i <= nrows_; ++i)
        {
            row_ptr_[static_cast<Size>(i)] += row_ptr_[static_cast<Size>(i - 1)];
        }

        col_ind_.assign(cols_u.begin(), cols_u.end());
        values_.assign(vals_u.begin(), vals_u.end());
    }


    //=========================================================================
    // Returns the number of matrix rows.
    //=========================================================================
    Int CSRMatrix::nrows() const
    {
        return nrows_;
    }


    //=========================================================================
    // Returns the number of matrix columns.
    //=========================================================================
    Int CSRMatrix::ncols() const
    {
        return ncols_;
    }


    //=========================================================================
    // Returns the number of coefficients stored in CSR form.
    //=========================================================================
    Int CSRMatrix::nnz() const
    {
        return static_cast<Int>(values_.size());
    }


    //=========================================================================
    // Returns read-only access to the CSR row-pointer array.
    //=========================================================================
    const std::vector<Int>& CSRMatrix::rowPtr() const
    {
        return row_ptr_;
    }


    //=========================================================================
    // Returns read-only access to the stored one-based column indices.
    //=========================================================================
    const std::vector<Int>& CSRMatrix::colInd() const
    {
        return col_ind_;
    }


    //=========================================================================
    // Returns read-only access to the stored coefficient values.
    //=========================================================================
    const std::vector<Real>& CSRMatrix::values() const
    {
        return values_;
    }


    //=========================================================================
    // Returns writable access to the CSR row-pointer array.
    //
    // The caller must preserve monotonic row offsets and the final nnz value.
    //=========================================================================
    std::vector<Int>& CSRMatrix::rowPtr()
    {
        return row_ptr_;
    }


    //=========================================================================
    // Returns writable access to the stored column indices.
    //
    // Every column number must remain in the one-based range 1 ... ncols.
    //=========================================================================
    std::vector<Int>& CSRMatrix::colInd()
    {
        return col_ind_;
    }


    //=========================================================================
    // Returns writable access to the stored numerical coefficients.
    //=========================================================================
    std::vector<Real>& CSRMatrix::values()
    {
        return values_;
    }


    //=========================================================================
    // Calculates the sparse matrix-vector product:
    //
    //     y = A x
    //
    // For matrix row i:
    //
    //     y_i =
    //         sum_(k=row_ptr[i-1])^(row_ptr[i]-1)
    //         values[k] x(col_ind[k])
    //
    // Input:
    //     x   one-based vector of length ncols
    //
    // Output:
    //     y   one-based vector resized to nrows
    //=========================================================================
    void CSRMatrix::matVec(const Array1D<Real>& x, Array1D<Real>& y) const
    {
        if (x.size() != ncols_)
        {
            throw std::runtime_error("CSRMatrix::matVec - input vector size mismatch");
        }

        y.resize(nrows_);
        y.fill(0.0);

        for (Int i = 1; i <= nrows_; ++i)
        {
            const Int k0 = row_ptr_[static_cast<Size>(i - 1)];
            const Int k1 = row_ptr_[static_cast<Size>(i)];

            Real sum = 0.0;

            for (Int k = k0; k < k1; ++k)
            {
                const Int j = col_ind_[static_cast<Size>(k)];
                sum += values_[static_cast<Size>(k)] * x(j);
            }

            y(i) = sum;
        }
    }
}
