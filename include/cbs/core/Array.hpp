#pragma once

//=============================================================================
// CBS3D++_SI
//
// Small array containers with one-based scientific indexing.
//
// The solver equations and legacy CBS data structures use indices beginning at
// one. These containers preserve that convention while storing their entries
// in contiguous std::vector memory.
//
// Storage order:
//
//     Array1D:
//         A(i)
//
//     Array2D, column-major:
//         offset(i,j) = (j-1)n1 + i
//
//     Array3D, column-major:
//         offset(i,j,k) = [(k-1)n2 + (j-1)]n1 + i
//
// Position zero in the underlying vector is intentionally unused. This allows
// element, node and boundary loops to follow the mathematical indexing used by
// the finite-element formulation.
//
// Bounds checking is active in debug builds and omitted when NDEBUG is defined.
// Allocation offsets are calculated with Size to avoid integer overflow for
// large three-dimensional meshes.
//=============================================================================

#include "cbs/core/Types.hpp"

#include <stdexcept>
#include <vector>

namespace cbs
{
    template <typename T>
    class Array1D
    {
    public:
        Array1D() = default;

        explicit Array1D(const Int n) { resize(n); }

        // Allocates n usable entries plus the unused position zero.
        void resize(const Int n)
        {
            if (n < 0) throw std::runtime_error("Array1D::resize - negative size");
            n_ = n;
            data_.assign(static_cast<Size>(n_) + 1U, T{});
        }

        Int size() const { return n_; }

        // Assigns one value to every active array entry.
        // Assigns one value to every active array entry.
        void fill(const T& value)
        {
            for (Int i = 1; i <= n_; ++i) data_[static_cast<Size>(i)] = value;
        }

        T& operator()(const Int i)
        {
            check_index(i);
            return data_[static_cast<Size>(i)];
        }
        const T& operator()(const Int i) const
        {
            check_index(i);
            return data_[static_cast<Size>(i)];
        }

        std::vector<T>& raw() { return data_; }
        const std::vector<T>& raw() const { return data_; }

    private:
        void check_index([[maybe_unused]] const Int i) const
        {
#ifndef NDEBUG
            if (i < 1 || i > n_) throw std::out_of_range("Array1D index out of bounds");
#endif
        }

        Int n_ = 0;
        std::vector<T> data_;
    };

    template <typename T>
    class Array2D
    {
    public:
        Array2D() = default;

        Array2D(const Int n1, const Int n2) { resize(n1, n2); }

        // Allocates a column-major n1 x n2 array plus position zero.
        void resize(const Int n1, const Int n2)
        {
            if (n1 < 0 || n2 < 0) throw std::runtime_error("Array2D::resize - negative size");
            n1_ = n1;
            n2_ = n2;
            data_.assign(static_cast<Size>(n1_) * static_cast<Size>(n2_) + 1U, T{});
        }

        Int dim1() const { return n1_; }
        Int dim2() const { return n2_; }

        void fill(const T& value)
        {
            for (Int j = 1; j <= n2_; ++j)
                for (Int i = 1; i <= n1_; ++i) (*this)(i, j) = value;
        }

        T& operator()(const Int i, const Int j)
        {
            check_index(i, j);
            return data_[offset(i, j)];
        }
        const T& operator()(const Int i, const Int j) const
        {
            check_index(i, j);
            return data_[offset(i, j)];
        }

        std::vector<T>& raw() { return data_; }
        const std::vector<T>& raw() const { return data_; }

    private:
        Size offset(const Int i, const Int j) const
        {
            return (static_cast<Size>(j) - 1U) * static_cast<Size>(n1_) + static_cast<Size>(i);
        }
        void check_index([[maybe_unused]] const Int i, [[maybe_unused]] const Int j) const
        {
#ifndef NDEBUG
            if (i < 1 || i > n1_ || j < 1 || j > n2_)
                throw std::out_of_range("Array2D index out of bounds");
#endif
        }

        Int n1_ = 0;
        Int n2_ = 0;
        std::vector<T> data_;
    };

    template <typename T>
    class Array3D
    {
        // 1-based, column-major: offset = ((k-1)*n2 + (j-1))*n1 + i
        // Mirrors legacy 3-index arrays such as annxf(ndim1,nsid,nelem).
    public:
        Array3D() = default;

        Array3D(const Int n1, const Int n2, const Int n3) { resize(n1, n2, n3); }

        // Allocates a column-major n1 x n2 x n3 array plus position zero.
        void resize(const Int n1, const Int n2, const Int n3)
        {
            if (n1 < 0 || n2 < 0 || n3 < 0) throw std::runtime_error("Array3D::resize - negative size");
            n1_ = n1;
            n2_ = n2;
            n3_ = n3;
            data_.assign(static_cast<Size>(n1_) * static_cast<Size>(n2_) * static_cast<Size>(n3_) + 1U, T{});
        }

        Int dim1() const { return n1_; }
        Int dim2() const { return n2_; }
        Int dim3() const { return n3_; }

        void fill(const T& value)
        {
            for (Size m = 1; m < data_.size(); ++m) data_[m] = value;
        }

        T& operator()(const Int i, const Int j, const Int k)
        {
            check_index(i, j, k);
            return data_[offset(i, j, k)];
        }
        const T& operator()(const Int i, const Int j, const Int k) const
        {
            check_index(i, j, k);
            return data_[offset(i, j, k)];
        }

        std::vector<T>& raw() { return data_; }
        const std::vector<T>& raw() const { return data_; }

    private:
        Size offset(const Int i, const Int j, const Int k) const
        {
            return ((static_cast<Size>(k) - 1U) * static_cast<Size>(n2_)
                + (static_cast<Size>(j) - 1U)) * static_cast<Size>(n1_)
                + static_cast<Size>(i);
        }
        void check_index([[maybe_unused]] const Int i,
            [[maybe_unused]] const Int j,
            [[maybe_unused]] const Int k) const
        {
#ifndef NDEBUG
            if (i < 1 || i > n1_ || j < 1 || j > n2_ || k < 1 || k > n3_)
                throw std::out_of_range("Array3D index out of bounds");
#endif
        }

        Int n1_ = 0;
        Int n2_ = 0;
        Int n3_ = 0;
        std::vector<T> data_;
    };
}