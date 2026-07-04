//=============================================================================
// CBS3D++_SI
//
// One-based square banded-matrix storage.
//
// For a matrix with upper bandwidth q, coefficient A(i,j) is stored as:
//
//     a(q + 1 + i - j, j)
//
// Only coefficients satisfying:
//
//     -lower_bandwidth <= j - i <= upper_bandwidth
//
// are represented.
//=============================================================================

#include "cbs/linalg/BandedMatrix.hpp"

#include <stdexcept>

namespace cbs
{
    //=========================================================================
    // Removes all matrix storage and resets the dimensions and bandwidths.
    //=========================================================================
    void BandedMatrix::clear()
    {
        n_ = 0;
        lower_bw_ = 0;
        upper_bw_ = 0;
        band_rows_ = 0;
        a_ = Array2D<Real>();
    }


    //=========================================================================
    // Allocates the compact storage for an n x n banded matrix.
    //
    // Number of stored rows:
    //
    //     n_band_rows = lower_bw + upper_bw + 1
    //
    // Storage dimensions:
    //
    //     a_(1:n_band_rows, 1:n)
    //
    // Every stored coefficient is initialised to zero.
    //=========================================================================
    void BandedMatrix::resize(
        const Int n,
        const Int lower_bw,
        const Int upper_bw)
    {
        if (n < 0 || lower_bw < 0 || upper_bw < 0)
        {
            throw std::runtime_error("BandedMatrix::resize - negative dimension or bandwidth");
        }

        n_ = n;
        lower_bw_ = lower_bw;
        upper_bw_ = upper_bw;
        band_rows_ = lower_bw_ + upper_bw_ + 1;

        a_.resize(band_rows_, n_);
        a_.fill(0.0);
    }


    //=========================================================================
    // Assigns value to every coefficient represented by the compact band
    // storage.
    //=========================================================================
    void BandedMatrix::fill(const Real value)
    {
        a_.fill(value);
    }


    //=========================================================================
    // Returns the full matrix dimension n.
    //=========================================================================
    Int BandedMatrix::size() const
    {
        return n_;
    }


    //=========================================================================
    // Returns the number of stored subdiagonals.
    //=========================================================================
    Int BandedMatrix::lowerBandwidth() const
    {
        return lower_bw_;
    }


    //=========================================================================
    // Returns the number of stored superdiagonals.
    //=========================================================================
    Int BandedMatrix::upperBandwidth() const
    {
        return upper_bw_;
    }


    //=========================================================================
    // Tests whether coefficient A(i,j) is represented.
    //
    // A valid stored coefficient must satisfy:
    //
    //     1 <= i,j <= n
    //
    // and:
    //
    //     j - i <= upper_bw
    //
    //     i - j <= lower_bw
    //=========================================================================
    bool BandedMatrix::inBand(
        const Int i,
        const Int j) const
    {
        if (i < 1 || i > n_ || j < 1 || j > n_)
        {
            return false;
        }

        const Int d = j - i;
        return (d <= upper_bw_ && -d <= lower_bw_);
    }


    //=========================================================================
    // Returns checked writable access to A(i,j).
    //=========================================================================
    Real& BandedMatrix::operator()(
        const Int i,
        const Int j)
    {
        checkIndex(i, j);
        return a_(bandRow(i, j), j);
    }


    //=========================================================================
    // Returns checked read-only access to A(i,j).
    //=========================================================================
    const Real& BandedMatrix::operator()(
        const Int i,
        const Int j) const
    {
        checkIndex(i, j);
        return a_(bandRow(i, j), j);
    }


    //=========================================================================
    // Adds one contribution to A(i,j).
    //
    // Contributions outside the represented band are ignored. This supports
    // finite-element assembly when the caller has already chosen a fixed
    // global bandwidth.
    //=========================================================================
    void BandedMatrix::add(
        const Int i,
        const Int j,
        const Real value)
    {
        if (!inBand(i, j))
        {
            return;
        }

        a_(bandRow(i, j), j) += value;
    }


    //=========================================================================
    // Returns read-only access to the compact Array2D storage.
    //=========================================================================
    const Array2D<Real>& BandedMatrix::data() const
    {
        return a_;
    }


    //=========================================================================
    // Returns writable access to the compact Array2D storage.
    //=========================================================================
    Array2D<Real>& BandedMatrix::data()
    {
        return a_;
    }


    //=========================================================================
    // Converts full-matrix coordinates A(i,j) to the compact storage row:
    //
    //     band_row(i,j) = upper_bw + 1 + i - j
    //
    // The compact storage column is the original matrix column j.
    //=========================================================================
    Int BandedMatrix::bandRow(
        const Int i,
        const Int j) const
    {
        return upper_bw_ + 1 + i - j;
    }


    //=========================================================================
    // Verifies matrix bounds and band membership before checked coefficient
    // access.
    //=========================================================================
    void BandedMatrix::checkIndex(
        const Int i,
        const Int j) const
    {
        if (i < 1 || i > n_ || j < 1 || j > n_)
        {
            throw std::out_of_range("BandedMatrix index out of bounds");
        }

        if (!inBand(i, j))
        {
            throw std::out_of_range("BandedMatrix access outside band");
        }
    }
}
