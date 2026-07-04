//=============================================================================
// CBS3D++_SI
//
// Matrix-vector products for the compact tetrahedral pressure operator.
//
// The scalar operation is:
//
//     y = A x
//
// where A is assembled from:
//
//     1. six off-diagonal coefficients stored for each tetrahedron;
//     2. one assembled diagonal coefficient stored at each mesh node.
//
// For one tetrahedral element:
//
//     y_1^(e) = A_12 x_2 + A_13 x_3 + A_14 x_4
//
//     y_2^(e) = A_12 x_1 + A_23 x_3 + A_24 x_4
//
//     y_3^(e) = A_13 x_1 + A_23 x_2 + A_34 x_4
//
//     y_4^(e) = A_14 x_1 + A_24 x_2 + A_34 x_3
//
// After all element contributions are assembled:
//
//     y_i = y_i + A_ii x_i
//
// Pressure is active only in fluid elements. Prescribed-pressure rows are
// replaced by identity rows so that:
//
//     y_i = x_i
//
// on every fixed-pressure node.
//=============================================================================

#include "cbs/linalg/MatrixVectorCalc.hpp"

#include <stdexcept>
#include <string>

namespace cbs
{
    namespace
    {
        // Checks the element dimensions required by the compact four-node
        // tetrahedral pressure operator.
        void validate_pressure_operator_dimensions(const CBSStateSI& s)
        {
            if (s.cfg.nep != 4 || s.cfg.gsdim != 6)
            {
                throw std::runtime_error(
                    "MatrixVectorCalc - pressureMultiply requires tetrahedral nep=4 and gsdim=6");
            }
        }


        // Replaces every prescribed-pressure row by an identity row:
        //
        //     y_i = x_i
        //
        // The fixed-node list is used directly, avoiding a search through all
        // mesh nodes during every Pressure CG matrix-vector product.
        void enforce_fixed_pressure_rows(
            const CBSStateSI& s,
            const Array1D<Real>& x,
            Array1D<Real>& y)
        {
            // Fixed pressure/reference rows are identity rows.
            //
            // The old implementation scanned every node and, for each node,
            // searched the whole bc_list.  On large CHT meshes this made every
            // pressure matvec perform O(npoin * bc_fixed) unnecessary work.
            // The fixed pressure list is already available, so apply only those
            // rows directly.
            for (Int i = 1; i <= s.cfg.bc_fixed; ++i)
            {
                const Int ip = s.bc_list(i);

                if (ip < 1 || ip > s.cfg.npoin)
                {
                    throw std::runtime_error(
                        "MatrixVectorCalc::pressureMultiply - fixed pressure node out of range");
                }

                y(ip) = x(ip);
            }
        }
    }


    //=========================================================================
    // Calculates the scalar pressure matrix-vector product:
    //
    //     y = A x
    //
    // The operation is performed in three stages:
    //
    //     1. Add element off-diagonal contributions.
    //     2. Add the assembled nodal diagonal contribution.
    //     3. Impose prescribed-pressure identity rows.
    //
    // The compact off-diagonal ordering is:
    //
    //     1 -> (1,2)
    //     2 -> (1,3)
    //     3 -> (1,4)
    //     4 -> (2,3)
    //     5 -> (2,4)
    //     6 -> (3,4)
    //
    // Elements are processed colour by colour. Within one colour, no two
    // tetrahedra share a node, so the nodal scatter is race-free under OpenMP.
    //
    // Inputs:
    //     x          input pressure vector
    //     gstif      element off-diagonal pressure coefficients
    //     pdiag      assembled pressure diagonal
    //     intma      tetrahedral connectivity
    //     mat_elem   material number
    //
    // Output:
    //     y          pressure matrix-vector product
    //=========================================================================
    void MatrixVectorCalc::pressureMultiply(
        const CBSStateSI& s,
        const Array1D<Real>& x,
        Array1D<Real>& y)
    {
        validate_pressure_operator_dimensions(s);

        y.fill(0.0);

        // ---------------------------------------------------------------------
        // Element off-diagonal action for the pressure Poisson operator.
        //
        // IMPORTANT CHT rule:
        //     Pressure/momentum exist only in the fluid domain.
        //
        // The CHT mesh contains many solid tetrahedra.  Even if their pressure
        // stiffness entries are zeroed elsewhere, looping over solid elements in
        // every CG matvec is a major waste.  The pressure operator must therefore
        // skip mat_elem != 0 directly inside the hot matvec loop.
        //
        // Compact pair ordering inherited from legacy multi4nPdiag for ndim=3:
        //   1 : (1,2)
        //   2 : (1,3)
        //   3 : (1,4)
        //   4 : (2,3)
        //   5 : (2,4)
        //   6 : (3,4)
        // ---------------------------------------------------------------------
        for (Int c = 0; c < s.ncolor; ++c)
        {
            const Int beg = s.color_ptr[static_cast<std::size_t>(c)];
            const Int end = s.color_ptr[static_cast<std::size_t>(c) + 1];

#pragma omp parallel for schedule(static)
            for (Int k = beg; k < end; ++k)
            {
                const Int ie = s.color_elem[static_cast<std::size_t>(k)];

                if (s.mat_elem(ie) != 0)
                {
                    continue;
                }

                const Int i1 = s.intma(1, ie);
                const Int i2 = s.intma(2, ie);
                const Int i3 = s.intma(3, ie);
                const Int i4 = s.intma(4, ie);

                const Real x1 = x(i1);
                const Real x2 = x(i2);
                const Real x3 = x(i3);
                const Real x4 = x(i4);

                const Int base = (ie - 1) * s.cfg.gsdim;

                const Real a12 = s.gstif(base + 1);
                const Real a13 = s.gstif(base + 2);
                const Real a14 = s.gstif(base + 3);
                const Real a23 = s.gstif(base + 4);
                const Real a24 = s.gstif(base + 5);
                const Real a34 = s.gstif(base + 6);

                // Race-free: within a color no two elements share a node.
                y(i1) += a12 * x2 + a13 * x3 + a14 * x4;
                y(i2) += a12 * x1 + a23 * x3 + a24 * x4;
                y(i3) += a13 * x1 + a23 * x2 + a34 * x4;
                y(i4) += a14 * x1 + a24 * x2 + a34 * x3;
            }
        }

        // Diagonal action is assembled nodally in pdiag.  Solid-only nodes have
        // pdiag = 0 in the corrected CHT preprocessing/time-step path, so this
        // loop is safe over all nodes.
#pragma omp parallel for schedule(static)
        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            y(ip) += s.pdiag(ip) * x(ip);
        }

        enforce_fixed_pressure_rows(s, x, y);
    }


    //=========================================================================
    // Applies the same compact pressure operator to every component of a nodal
    // vector field:
    //
    //     y_d = A x_d,    d = 1, 2, 3
    //
    // The element off-diagonal and nodal diagonal actions are identical to the
    // scalar pressure operation and are repeated independently for each
    // Cartesian component.
    //
    // This routine is not the principal Pressure CG matrix-vector kernel.
    //=========================================================================
    void MatrixVectorCalc::pressureMultiplyVector(
        const CBSStateSI& s,
        const Array2D<Real>& x,
        Array2D<Real>& y)
    {
        validate_pressure_operator_dimensions(s);

        y.fill(0.0);

        // Vector form of the same compact operator.  This is not the main
        // pressure-CG bottleneck, but it must obey the same CHT rule: pressure-
        // style matrix action is fluid-element-only.
        for (Int c = 0; c < s.ncolor; ++c)
        {
            const Int beg = s.color_ptr[static_cast<std::size_t>(c)];
            const Int end = s.color_ptr[static_cast<std::size_t>(c) + 1];

#pragma omp parallel for schedule(static)
            for (Int k = beg; k < end; ++k)
            {
                const Int ie = s.color_elem[static_cast<std::size_t>(k)];

                if (s.mat_elem(ie) != 0)
                {
                    continue;
                }

                const Int i1 = s.intma(1, ie);
                const Int i2 = s.intma(2, ie);
                const Int i3 = s.intma(3, ie);
                const Int i4 = s.intma(4, ie);

                const Int base = (ie - 1) * s.cfg.gsdim;

                const Real a12 = s.gstif(base + 1);
                const Real a13 = s.gstif(base + 2);
                const Real a14 = s.gstif(base + 3);
                const Real a23 = s.gstif(base + 4);
                const Real a24 = s.gstif(base + 5);
                const Real a34 = s.gstif(base + 6);

                for (Int idim = 1; idim <= s.cfg.ndim; ++idim)
                {
                    const Real x1 = x(idim, i1);
                    const Real x2 = x(idim, i2);
                    const Real x3 = x(idim, i3);
                    const Real x4 = x(idim, i4);

                    y(idim, i1) += a12 * x2 + a13 * x3 + a14 * x4;
                    y(idim, i2) += a12 * x1 + a23 * x3 + a24 * x4;
                    y(idim, i3) += a13 * x1 + a23 * x2 + a34 * x4;
                    y(idim, i4) += a14 * x1 + a24 * x2 + a34 * x3;
                }
            }
        }

#pragma omp parallel for schedule(static)
        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            for (Int idim = 1; idim <= s.cfg.ndim; ++idim)
            {
                y(idim, ip) += s.pdiag(ip) * x(idim, ip);
            }
        }
    }
}
