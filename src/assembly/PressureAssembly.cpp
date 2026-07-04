//=============================================================================
// CBS3D++_SI
//
// Pressure assembly for the semi-implicit CBS pressure-correction step.
//
// This file performs two operations:
//
//     1. Calculate the element pressure stiffness matrix.
//
//            H_ij^(e) = V_e grad(N_i) . grad(N_j)
//
//     2. Assemble the Step 2 pressure right-hand side from the weak divergence
//        of the predicted velocity field.
//
// Pressure is active only in fluid elements, identified by:
//
//     mat_elem(e) = 0
//=============================================================================

#include "cbs/assembly/PressureAssembly.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace cbs
{
    namespace
    {
        // Returns the one-dimensional storage position of
        //
        //     dN_local_node / dx_dim
        //
        // for element ie.
        Int dNkdx_index(
            const CBSStateSI& s,
            Int ie,
            Int dim,
            Int local_node)
        {
            return (ie - 1) * s.cfg.ndim * s.cfg.nep
                + (dim - 1) * s.cfg.nep
                + local_node;
        }


        // Returns the compact storage position of a local tetrahedral node.
        Int element_node_index(
            const CBSStateSI& s,
            Int ie,
            Int local_node)
        {
            return (ie - 1) * s.cfg.nep + local_node;
        }


        // Returns the compact storage position of one off-diagonal
        // tetrahedral node pair.
        Int offdiag_index(
            const CBSStateSI& s,
            Int ie,
            Int pair_id)
        {
            return (ie - 1) * s.cfg.gsdim + pair_id;
        }


        // Pressure is assembled only on fluid elements.
        bool is_fluid_element(
            const CBSStateSI& s,
            Int ie)
        {
            return s.mat_elem(ie) == 0;
        }


        // Checks the fixed dimensions required by the CBS3D P1 tetrahedral
        // pressure formulation.
        void validate_pressure_dimensions(const CBSStateSI& s)
        {
            if (s.cfg.ndim != 3 ||
                s.cfg.nep != 4 ||
                s.cfg.gsdim != 6 ||
                s.cfg.gdim != 13)
            {
                throw std::runtime_error(
                    "PressureAssembly - CBS3D pressure assembly requires ndim=3, nep=4, gsdim=6, gdim=13");
            }
        }


        // Calculates
        //
        //     grad(N_a) . grad(N_b)
        //
        // for one linear tetrahedral element.
        Real grad_dot(
            const CBSStateSI& s,
            Int ie,
            Int a,
            Int b)
        {
            return
                s.dNkdx(dNkdx_index(s, ie, 1, a)) * s.dNkdx(dNkdx_index(s, ie, 1, b)) +
                s.dNkdx(dNkdx_index(s, ie, 2, a)) * s.dNkdx(dNkdx_index(s, ie, 2, b)) +
                s.dNkdx(dNkdx_index(s, ie, 3, a)) * s.dNkdx(dNkdx_index(s, ie, 3, b));
        }
    }


    //=========================================================================
    // Calculates the compact pressure-stiffness matrix for every fluid
    // tetrahedron.
    //
    // Governing element expression:
    //
    //     H_ab^(e) = integral(V_e) grad(N_a) . grad(N_b) dV
    //
    // For a P1 tetrahedron, grad(N_a) is constant. Therefore:
    //
    //     H_ab^(e) = V_e grad(N_a) . grad(N_b)
    //
    // The tetrahedral volume is:
    //
    //     V_e = det(J_e) / 6
    //
    // Inputs:
    //     detJ(e)       tetrahedral Jacobian determinant
    //     dNkdx         Cartesian shape-function derivatives
    //     mat_elem(e)   material number
    //
    // Outputs:
    //     pdiagE        four element diagonal entries
    //     gstifE        six element off-diagonal entries
    //=========================================================================
    void PressureAssembly::buildElementPressureTerms(CBSStateSI& s)
    {
        validate_pressure_dimensions(s);

        s.pdiagE.fill(0.0);
        s.gstifE.fill(0.0);

        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            if (!is_fluid_element(s, ie))
            {
                continue;
            }

            if (s.detJ(ie) <= 0.0 || !std::isfinite(s.detJ(ie)))
            {
                throw std::runtime_error(
                    "PressureAssembly::buildElementPressureTerms - invalid detJ at element "
                    + std::to_string(ie));
            }

            const Real volume = s.detJ(ie) / 6.0;

            for (Int a = 1; a <= s.cfg.nep; ++a)
            {
                s.pdiagE(element_node_index(s, ie, a)) =
                    volume * grad_dot(s, ie, a, a);
            }

            // Compact upper-triangular ordering:
            //
            //     (1,2), (1,3), (1,4), (2,3), (2,4), (3,4)
            //
            // This ordering must remain consistent with the pressure
            // matrix-vector multiplication.
            s.gstifE(offdiag_index(s, ie, 1)) = volume * grad_dot(s, ie, 1, 2);
            s.gstifE(offdiag_index(s, ie, 2)) = volume * grad_dot(s, ie, 1, 3);
            s.gstifE(offdiag_index(s, ie, 3)) = volume * grad_dot(s, ie, 1, 4);
            s.gstifE(offdiag_index(s, ie, 4)) = volume * grad_dot(s, ie, 2, 3);
            s.gstifE(offdiag_index(s, ie, 5)) = volume * grad_dot(s, ie, 2, 4);
            s.gstifE(offdiag_index(s, ie, 6)) = volume * grad_dot(s, ie, 3, 4);
        }
    }


    //=========================================================================
    // Initialises the global pressure arrays.
    //
    // The final active pressure operator is not assembled here. The element
    // coefficients pdiagE and gstifE are multiplied by the local CBS time step
    // inside TimeStep::updateLhsDiagonal().
    //
    // No CSR or banded matrix is constructed in the matrix-free pressure path.
    //=========================================================================
    void PressureAssembly::buildGlobalPressureTerms(CBSStateSI& s)
    {
        s.pdiag.fill(0.0);
        s.gstif.fill(0.0);
    }


    //=========================================================================
    // Assembles the pressure right-hand side for CBS Step 2.
    //
    // The volume contribution represents the weak divergence of the predicted
    // velocity:
    //
    //     b_a^(e) = integral(V_e) grad(N_a) . u* dV
    //
    // Using the P1 tetrahedral nodal quadrature employed by the legacy method:
    //
    //     integral(V_e) u* dV
    //         approximately (V_e / 4) sum_b u_b*
    //
    // giving:
    //
    //     b_a^(e) =
    //         grad(N_a) . [(V_e / 4) sum_b u_b*]
    //
    // Boundary-face terms use the previous velocity and the area-weighted
    // outward normal. For a triangular P1 face, fcon[2] = 1/12 is the
    // consistent integration factor.
    //
    // Important time-step convention:
    //
    //     RHS      = raw weak divergence
    //     operator = dt_e K
    //     Step 3   = u* - dt_e M_L^(-1) G delta_p
    //
    // Therefore the RHS must not be divided by dtreal in this routine.
    //
    // Inputs:
    //     unkno    predicted velocity
    //     unkn1    previous velocity used by boundary-face terms
    //     dNkdx    shape-function gradients
    //     annxf    area-weighted boundary normals
    //
    // Output:
    //     rhs1     assembled nodal pressure RHS
    //=========================================================================
    void PressureAssembly::assembleStep2Rhs(CBSStateSI& s)
    {
        validate_pressure_dimensions(s);

        if (s.cfg.dtreal <= 0.0 || !std::isfinite(s.cfg.dtreal))
        {
            throw std::runtime_error(
                "PressureAssembly::assembleStep2Rhs - dtreal must be positive before Step 2 RHS assembly");
        }

        s.rhs1.fill(0.0);

        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            if (!is_fluid_element(s, ie))
            {
                continue;
            }

            if (s.detJ(ie) <= 0.0 || !std::isfinite(s.detJ(ie)))
            {
                throw std::runtime_error(
                    "PressureAssembly::assembleStep2Rhs - invalid detJ at element "
                    + std::to_string(ie));
            }

            const Real vol4 = s.detJ(ie) * s.cfg.fcon[1];

            Real u_sum[4] = { 0.0, 0.0, 0.0, 0.0 };

            for (Int in = 1; in <= s.cfg.nep; ++in)
            {
                const Int ip = s.intma(in, ie);

                u_sum[1] += s.unkno(1, ip);
                u_sum[2] += s.unkno(2, ip);
                u_sum[3] += s.unkno(3, ip);
            }

            u_sum[1] *= vol4;
            u_sum[2] *= vol4;
            u_sum[3] *= vol4;

            Real lrhs[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };

            for (Int in = 1; in <= s.cfg.nep; ++in)
            {
                lrhs[in] =
                    s.dNkdx(dNkdx_index(s, ie, 1, in)) * u_sum[1] +
                    s.dNkdx(dNkdx_index(s, ie, 2, in)) * u_sum[2] +
                    s.dNkdx(dNkdx_index(s, ie, 3, in)) * u_sum[3];
            }

            // Add the boundary-face contribution on marked element faces.
            for (Int is = 1; is <= s.cfg.nsid; ++is)
            {
                if (s.fedge(is, ie) != 1)
                {
                    continue;
                }

                Real face_u_sum[4] = { 0.0, 0.0, 0.0, 0.0 };

                for (Int i = 1; i <= s.cfg.nsidp; ++i)
                {
                    const Int local_node = s.ippn1(is, i);
                    const Int ip = s.intma(local_node, ie);

                    face_u_sum[1] += s.unkn1(1, ip);
                    face_u_sum[2] += s.unkn1(2, ip);
                    face_u_sum[3] += s.unkn1(3, ip);
                }

                for (Int idim = 1; idim <= s.cfg.ndim; ++idim)
                {
                    const Real normal_factor =
                        s.annxf(idim, is, ie) * s.cfg.fcon[2];

                    for (Int i = 1; i <= s.cfg.nsidp; ++i)
                    {
                        const Int local_node = s.ippn1(is, i);
                        const Int ip = s.intma(local_node, ie);

                        const Real node_u_old = s.unkn1(idim, ip);

                        lrhs[local_node] -=
                            (face_u_sum[idim] + node_u_old) * normal_factor;
                    }
                }
            }

            // Assemble the local tetrahedral contribution into the global
            // nodal pressure RHS.
            for (Int in = 1; in <= s.cfg.nep; ++in)
            {
                const Int ip = s.intma(in, ie);
                s.rhs1(ip) += lrhs[in];
            }
        }
    }


    //=========================================================================
    // Compatibility wrapper retained for callers using the older routine name.
    //=========================================================================
    void PressureAssembly::assembleRhs(CBSStateSI& s)
    {
        assembleStep2Rhs(s);
    }
}
