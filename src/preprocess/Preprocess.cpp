//=============================================================================
// CBS3D++_SI
//
// Finite-element preprocessing for the three-dimensional semi-implicit
// Characteristic-Based Split solver.
//
// The routines in this file prepare the mesh geometry, mass coefficients,
// boundary-face data and nodal boundary lists required before the CBS
// iterations begin. The spatial discretisation uses four-node linear
// tetrahedral elements with one-based element and node numbering.
//=============================================================================

#include "cbs/preprocess/Preprocess.hpp"

#ifdef CBS3D_USE_MPI
#include "cbs/parallel/HaloExchange.hpp"
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace cbs
{
    namespace
    {
        // Boundary identifiers recognised by the current solver. The identifiers
        // below 511 are retained for legacy validation cases. The CHT blanket
        // cases primarily use 511, 520, 530, 532 and 901.
        //
        // BC 902 is a legacy heat-flux marker. It is not a fluid-solid interface.
        constexpr Int BC_MOVING_WALL = 500;
        constexpr Int BC_T_HOT_WALL  = 501;
        constexpr Int BC_T_COLD_WALL = 502;
        constexpr Int BC_LEGACY_VEL  = 503;
        constexpr Int BC_PRESSURE    = 504;
        constexpr Int BC_SYMMETRY    = 506;
        constexpr Int BC_BFS_INLET   = 507;
        constexpr Int BC_PARAB_INLET = 508;
        constexpr Int BC_VEL_TEMP    = 510;
        constexpr Int BC_INLET       = 511;
        constexpr Int BC_OUTLET      = 520;
        constexpr Int BC_ADIABATIC   = 530;
        constexpr Int BC_HEATFLUX    = 532;
        constexpr Int BC_INTERFACE   = 901;
        constexpr Int BC_HFLUX_MARK  = 902;

        //-------------------------------------------------------------------------
        // Returns true when the boundary identifier is recognised by the solver.
        //
        // This function only classifies identifiers. The physical action
        // associated with each identifier is applied later by Boundary and the
        // CBS assembly routines.
        //-------------------------------------------------------------------------
        bool is_supported_cbs3d_bc(Int bc)
        {
            switch (bc)
            {
                case BC_MOVING_WALL:
                case BC_T_HOT_WALL:
                case BC_T_COLD_WALL:
                case BC_LEGACY_VEL:
                case BC_PRESSURE:
                case BC_SYMMETRY:
                case BC_BFS_INLET:
                case BC_PARAB_INLET:
                case BC_VEL_TEMP:
                case BC_INLET:
                case BC_OUTLET:
                case BC_ADIABATIC:
                case BC_HEATFLUX:
                case BC_INTERFACE:
                case BC_HFLUX_MARK:
                    return true;

                default:
                    return false;
            }
        }

        // Returns true for the pressure-outlet boundary used to prescribe
        // the pressure correction/reference value.
        bool is_pressure_outlet_bc(Int bc)
        {
            return bc == BC_OUTLET;
        }

        // Returns true for the inlet where the velocity magnitude is obtained
        // from the prescribed mass-flow rate.
        bool is_mass_flow_inlet_bc(Int bc)
        {
            return bc == BC_INLET;
        }

        //-------------------------------------------------------------------------
        // Identifies boundary faces whose nodes belong to the strict no-slip list.
        //
        // Moving-wall BC 500 is deliberately excluded because its non-zero
        // velocity is imposed separately by Boundary::applyVelocity().
        //-------------------------------------------------------------------------
        bool is_no_slip_wall_bc(Int bc)
        {
            return bc == BC_T_HOT_WALL ||
                   bc == BC_T_COLD_WALL ||
                   bc == BC_ADIABATIC ||
                   bc == BC_HEATFLUX ||
                   bc == BC_INTERFACE;
        }

        //-------------------------------------------------------------------------
        // Identifies faces on which velocity is imposed strongly.
        //
        // fedge = 1  prescribed-velocity face
        // fedge = 2  other exterior face
        // fedge = 0  interior face
        //-------------------------------------------------------------------------
        bool is_prescribed_velocity_face_bc(Int bc)
        {
            return bc == BC_MOVING_WALL ||
                   bc == BC_T_HOT_WALL ||
                   bc == BC_T_COLD_WALL ||
                   bc == BC_LEGACY_VEL ||
                   bc == BC_BFS_INLET ||
                   bc == BC_PARAB_INLET ||
                   bc == BC_VEL_TEMP ||
                   bc == BC_INLET ||
                   bc == BC_ADIABATIC ||
                   bc == BC_HEATFLUX ||
                   bc == BC_INTERFACE;
        }

        //-------------------------------------------------------------------------
        // Converts (element, Cartesian direction, local node) into the one-based
        // storage position used by dNkdx.
        //
        // For each element, the array stores:
        //
        //     dN_1/dx ... dN_4/dx
        //     dN_1/dy ... dN_4/dy
        //     dN_1/dz ... dN_4/dz
        //-------------------------------------------------------------------------
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

        // Converts (element, local node) into the one-based position used by
        // element-node arrays such as elcoe_e.
        Int element_node_index(
            const CBSStateSI& s,
            Int ie,
            Int local_node)
        {
            return (ie - 1) * s.cfg.nep + local_node;
        }

        // Returns the three global node numbers of a triangular face in
        // ascending order so that face matching is independent of orientation.
        std::array<Int, 3> sorted_face_nodes(
            Int a,
            Int b,
            Int c)
        {
            std::array<Int, 3> nodes = { a, b, c };
            std::sort(nodes.begin(), nodes.end());
            return nodes;
        }

        // Returns true when two sorted triangular faces contain the same nodes.
        bool same_face_nodes(
            const std::array<Int, 3>& a,
            const std::array<Int, 3>& b)
        {
            return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
        }

        //-------------------------------------------------------------------------
        // Calculates the determinant of a 3 x 3 matrix:
        //
        //     |a11 a12 a13|
        //     |a21 a22 a23|
        //     |a31 a32 a33|
        //
        // The determinant is used to calculate det(J) for a tetrahedral element.
        //-------------------------------------------------------------------------
        Real determinant3(
            Real a11, Real a12, Real a13,
            Real a21, Real a22, Real a23,
            Real a31, Real a32, Real a33)
        {
            return a11 * (a22 * a33 - a23 * a32)
                - a12 * (a21 * a33 - a23 * a31)
                + a13 * (a21 * a32 - a22 * a31);
        }

        //-------------------------------------------------------------------------
        // Verifies the fixed topology used by the present discretisation:
        //
        //     spatial dimensions             ndim   = 3
        //     nodes per tetrahedron           nep    = 4
        //     faces per tetrahedron           nsid   = 4
        //     nodes per triangular face       nsidp  = 3
        //     face-normal entries             ndim1  = 4
        //     node pairs per tetrahedron      gsdim  = 6
        //     boundary-face record entries    bsid   = 6
        //-------------------------------------------------------------------------
        void validate_core_dimensions(const CBSStateSI& s)
        {
            if (s.cfg.ndim != 3 ||
                s.cfg.nep != 4 ||
                s.cfg.nsid != 4 ||
                s.cfg.nsidp != 3 ||
                s.cfg.ndim1 != 4 ||
                s.cfg.gsdim != 6 ||
                s.cfg.bsid != 6)
            {
                throw std::runtime_error(
                    "Preprocess - CBS3D CHT preprocess requires ndim=3, nep=4, nsid=4, nsidp=3, ndim1=4, gsdim=6, bsid=6");
            }
        }

        // Adds a node to a one-based list only if it has not already been added.
        void add_unique_node(
            Array1D<Int>& node_list,
            Int& count,
            std::vector<Int>& marker,
            Int ip)
        {
            if (ip < 1 || ip >= static_cast<Int>(marker.size()))
            {
                throw std::runtime_error(
                    "Preprocess - node index out of range while building unique node list");
            }

            if (marker[static_cast<Size>(ip)] != 0)
            {
                return;
            }

            ++count;
            node_list(count) = ip;
            marker[static_cast<Size>(ip)] = 1;
        }

        //-------------------------------------------------------------------------
        // Returns the validated persistent material-connectivity mask for one
        // node.
        //
        // The mask must already have been built by
        // Preprocess::buildMaterialNodeMasks(). In an MPI calculation that
        // routine reconciles owner and ghost copies before these predicates are
        // used.
        //-------------------------------------------------------------------------
        Int material_node_mask(
            const CBSStateSI& s,
            Int ip,
            const char* context)
        {
            if (ip < 1 || ip > s.cfg.npoin)
            {
                throw std::runtime_error(
                    std::string(context) + " - node index out of range");
            }

            const Int mask = s.node_material_mask(ip);
            const Int valid_mask =
                CBSStateSI::node_touches_fluid |
                CBSStateSI::node_touches_solid;

            if (mask < CBSStateSI::node_touches_fluid ||
                mask > valid_mask)
            {
                throw std::runtime_error(
                    std::string(context)
                    + " - material node mask has not been built or is invalid");
            }

            return mask;
        }

        // Returns true when the node belongs to at least one fluid element.
        bool touches_fluid_domain(
            const CBSStateSI& s,
            Int ip)
        {
            return
                (material_node_mask(
                    s,
                    ip,
                    "Preprocess::touches_fluid_domain")
                 & CBSStateSI::node_touches_fluid) != 0;
        }

        // Returns true when the node is shared by fluid and solid elements.
        bool is_conformal_fluid_solid_interface_node(
            const CBSStateSI& s,
            Int ip)
        {
            const Int interface_mask =
                CBSStateSI::node_touches_fluid |
                CBSStateSI::node_touches_solid;

            return
                material_node_mask(
                    s,
                    ip,
                    "Preprocess::is_conformal_fluid_solid_interface_node")
                == interface_mask;
        }
    }

    //=========================================================================
    // Validates the fixed tetrahedral dimensions and checks every boundary ID.
    //
    // Unknown identifiers are reported once and are left without a strong
    // boundary constraint. They therefore receive the natural finite-element
    // boundary treatment in the downstream assembly.
    //
    // Input:
    //     s.iside(bsid, ib)   solver boundary identifier of boundary face ib
    //
    // Output:
    //     No solver array is modified.
    //=========================================================================
    void Preprocess::validateBoundaryFlags(CBSStateSI& s)
    {
        validate_core_dimensions(s);

        for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
        {
            const Int bc = s.iside(s.cfg.bsid, ib);

            if (!is_supported_cbs3d_bc(bc))
            {
                // Report each unknown identifier once. The downstream solver
                // applies no strong boundary condition to an unknown identifier,
                // so it receives the natural finite-element boundary treatment.
                static std::set<Int> warned;
                if (warned.insert(bc).second)
                {
                    std::cerr << "WARNING: Preprocess::validateBoundaryFlags - unrecognized BC_ID "
                              << bc << " (first seen at boundary face " << ib
                              << "). It will receive the default natural boundary treatment. "
                                 "Recognized CBS3D BC_IDs: 500, 501, 502, 503, 504, 506, 507, 508, "
                                 "510, 511, 520, 530, 532, 901, 902.\n";
                }
            }
        }
    }

    //=========================================================================
    // Calculates the shape-function gradients and Jacobian determinant for
    // every four-node linear tetrahedral element.
    //
    // The physical coordinates are obtained from the reference tetrahedron by
    //
    //     x = x_1 + J [xi, eta, zeta]^T
    //
    // with
    //
    //         |x_2-x_1  x_3-x_1  x_4-x_1|
    //     J = |y_2-y_1  y_3-y_1  y_4-y_1|
    //         |z_2-z_1  z_3-z_1  z_4-z_1|
    //
    // The physical gradients of the P1 shape functions are constant within
    // each tetrahedron:
    //
    //     grad(N_2) = row 1 of J^(-1)
    //     grad(N_3) = row 2 of J^(-1)
    //     grad(N_4) = row 3 of J^(-1)
    //     grad(N_1) = -grad(N_2) - grad(N_3) - grad(N_4)
    //
    // The element volume is
    //
    //     V_e = det(J_e) / 6
    //
    // A positive det(J) is required, so the tetrahedral connectivity must use
    // a consistent positive orientation.
    //
    // Output:
    //     s.dNkdx   Cartesian shape-function gradients
    //     s.detJ    Jacobian determinant, equal to 6 V_e
    //=========================================================================
    void Preprocess::shapeFunctionDerivatives(CBSStateSI& s)
    {
        validate_core_dimensions(s);

        Real max_volume = 0.0;
        Real min_volume = std::numeric_limits<Real>::max();

        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            const Int i1 = s.intma(1, ie);
            const Int i2 = s.intma(2, ie);
            const Int i3 = s.intma(3, ie);
            const Int i4 = s.intma(4, ie);

            const Real x1 = s.coord(1, i1);
            const Real y1 = s.coord(2, i1);
            const Real z1 = s.coord(3, i1);

            const Real x2 = s.coord(1, i2);
            const Real y2 = s.coord(2, i2);
            const Real z2 = s.coord(3, i2);

            const Real x3 = s.coord(1, i3);
            const Real y3 = s.coord(2, i3);
            const Real z3 = s.coord(3, i3);

            const Real x4 = s.coord(1, i4);
            const Real y4 = s.coord(2, i4);
            const Real z4 = s.coord(3, i4);

            const Real j11 = x2 - x1;
            const Real j21 = y2 - y1;
            const Real j31 = z2 - z1;

            const Real j12 = x3 - x1;
            const Real j22 = y3 - y1;
            const Real j32 = z3 - z1;

            const Real j13 = x4 - x1;
            const Real j23 = y4 - y1;
            const Real j33 = z4 - z1;

            const Real detJ = determinant3(
                j11, j12, j13,
                j21, j22, j23,
                j31, j32, j33);

            if (detJ <= 0.0 || !std::isfinite(detJ))
            {
                throw std::runtime_error(
                    "Preprocess::shapeFunctionDerivatives - non-positive tetrahedron volume at element "
                    + std::to_string(ie)
                    + ". Tetrahedral connectivity must be consistently positively oriented.");
            }

            const Real inv_det = 1.0 / detJ;

            const Real inv11 = (j22 * j33 - j23 * j32) * inv_det;
            const Real inv12 = (j13 * j32 - j12 * j33) * inv_det;
            const Real inv13 = (j12 * j23 - j13 * j22) * inv_det;

            const Real inv21 = (j23 * j31 - j21 * j33) * inv_det;
            const Real inv22 = (j11 * j33 - j13 * j31) * inv_det;
            const Real inv23 = (j13 * j21 - j11 * j23) * inv_det;

            const Real inv31 = (j21 * j32 - j22 * j31) * inv_det;
            const Real inv32 = (j12 * j31 - j11 * j32) * inv_det;
            const Real inv33 = (j11 * j22 - j12 * j21) * inv_det;

            const Real gx2 = inv11;
            const Real gy2 = inv12;
            const Real gz2 = inv13;

            const Real gx3 = inv21;
            const Real gy3 = inv22;
            const Real gz3 = inv23;

            const Real gx4 = inv31;
            const Real gy4 = inv32;
            const Real gz4 = inv33;

            const Real gx1 = -(gx2 + gx3 + gx4);
            const Real gy1 = -(gy2 + gy3 + gy4);
            const Real gz1 = -(gz2 + gz3 + gz4);

            s.dNkdx(dNkdx_index(s, ie, 1, 1)) = gx1;
            s.dNkdx(dNkdx_index(s, ie, 2, 1)) = gy1;
            s.dNkdx(dNkdx_index(s, ie, 3, 1)) = gz1;

            s.dNkdx(dNkdx_index(s, ie, 1, 2)) = gx2;
            s.dNkdx(dNkdx_index(s, ie, 2, 2)) = gy2;
            s.dNkdx(dNkdx_index(s, ie, 3, 2)) = gz2;

            s.dNkdx(dNkdx_index(s, ie, 1, 3)) = gx3;
            s.dNkdx(dNkdx_index(s, ie, 2, 3)) = gy3;
            s.dNkdx(dNkdx_index(s, ie, 3, 3)) = gz3;

            s.dNkdx(dNkdx_index(s, ie, 1, 4)) = gx4;
            s.dNkdx(dNkdx_index(s, ie, 2, 4)) = gy4;
            s.dNkdx(dNkdx_index(s, ie, 3, 4)) = gz4;

            s.detJ(ie) = detJ;

            const Real volume = detJ / 6.0;
            max_volume = std::max(max_volume, volume);
            min_volume = std::min(min_volume, volume);
        }

        std::cout << "Maximum Volume      " << max_volume << "\n";
        std::cout << "Minimum Volume      " << min_volume << "\n";

        if (min_volume > 0.0)
        {
            std::cout << "Ratio of Max to Min " << (max_volume / min_volume) << "\n";
        }
    }

    //=========================================================================
    // Matches each triangular boundary face to one local face of its parent
    // tetrahedral element.
    //
    // The three boundary node numbers are sorted and compared with the four
    // local tetrahedral faces defined by ippn1. Sorting removes any dependence
    // on the node order used in the input boundary-face record.
    //
    // Input:
    //     s.iside(1:3, ib)     boundary-face node numbers
    //     s.iside(nsidpe, ib)  parent tetrahedral element
    //     s.intma              tetrahedral connectivity
    //
    // Output:
    //     s.iside(nsidpl, ib)  matched local face number, 1 to 4
    //=========================================================================
    void Preprocess::assignBoundaryFaceNumbers(CBSStateSI& s)
    {
        validate_core_dimensions(s);
        validateBoundaryFlags(s);

        for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
        {
            const Int parent = s.iside(s.cfg.nsidpe, ib);

            if (parent < 1 || parent > s.cfg.nelem)
            {
                throw std::runtime_error(
                    "Preprocess::assignBoundaryFaceNumbers - invalid parent tetrahedron at boundary face "
                    + std::to_string(ib));
            }

            const std::array<Int, 3> boundary_face =
                sorted_face_nodes(
                    s.iside(1, ib),
                    s.iside(2, ib),
                    s.iside(3, ib));

            Int matched_face = 0;

            for (Int is = 1; is <= s.cfg.nsid; ++is)
            {
                const Int p1 = s.intma(s.ippn1(is, 1), parent);
                const Int p2 = s.intma(s.ippn1(is, 2), parent);
                const Int p3 = s.intma(s.ippn1(is, 3), parent);

                const std::array<Int, 3> element_face =
                    sorted_face_nodes(p1, p2, p3);

                if (same_face_nodes(boundary_face, element_face))
                {
                    matched_face = is;
                    break;
                }
            }

            if (matched_face == 0)
            {
                throw std::runtime_error(
                    "Preprocess::assignBoundaryFaceNumbers - boundary face "
                    + std::to_string(ib)
                    + " does not match any face of parent tetrahedron "
                    + std::to_string(parent));
            }

            s.iside(s.cfg.nsidpl, ib) = matched_face;
        }
    }

    //=========================================================================
    // Calculates outward area-weighted normals and areas for tetrahedral faces.
    //
    // For face nodes x_1, x_2 and x_3:
    //
    //     a = x_2 - x_1
    //     b = x_3 - x_1
    //     c = a x b
    //
    // The triangle area is
    //
    //     A_f = |c| / 2
    //
    // and the area-weighted normal is
    //
    //     n_A = c / 2 = A_f n
    //
    // where n is the unit normal. The sign is reversed when c points towards
    // the tetrahedral node opposite the face, ensuring an outward direction.
    //
    // Output:
    //     s.annxf(1:3, face, element)  area-weighted outward normal
    //     s.annxf(4,   face, element)  face area
    //     s.face_norm(:, boundary)     corresponding boundary-face values
    //=========================================================================
    void Preprocess::getNormals(CBSStateSI& s)
    {
        validate_core_dimensions(s);

        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            for (Int is = 1; is <= s.cfg.nsid; ++is)
            {
                const Int p1 = s.intma(s.ippn1(is, 1), ie);
                const Int p2 = s.intma(s.ippn1(is, 2), ie);
                const Int p3 = s.intma(s.ippn1(is, 3), ie);
                const Int popp = s.intma(is, ie);

                const Real e1x = s.coord(1, p2) - s.coord(1, p1);
                const Real e1y = s.coord(2, p2) - s.coord(2, p1);
                const Real e1z = s.coord(3, p2) - s.coord(3, p1);

                const Real e2x = s.coord(1, p3) - s.coord(1, p1);
                const Real e2y = s.coord(2, p3) - s.coord(2, p1);
                const Real e2z = s.coord(3, p3) - s.coord(3, p1);

                Real ax = e1y * e2z - e1z * e2y;
                Real ay = e1z * e2x - e1x * e2z;
                Real az = e1x * e2y - e1y * e2x;

                const Real amag = std::sqrt(ax * ax + ay * ay + az * az);

                if (amag <= 0.0 || !std::isfinite(amag))
                {
                    throw std::runtime_error(
                        "Preprocess::getNormals - zero-area face at element "
                        + std::to_string(ie)
                        + ", face "
                        + std::to_string(is));
                }

                const Real fx =
                    (s.coord(1, p1) + s.coord(1, p2) + s.coord(1, p3)) / 3.0;
                const Real fy =
                    (s.coord(2, p1) + s.coord(2, p2) + s.coord(2, p3)) / 3.0;
                const Real fz =
                    (s.coord(3, p1) + s.coord(3, p2) + s.coord(3, p3)) / 3.0;

                const Real vx = s.coord(1, popp) - fx;
                const Real vy = s.coord(2, popp) - fy;
                const Real vz = s.coord(3, popp) - fz;

                if (ax * vx + ay * vy + az * vz > 0.0)
                {
                    ax = -ax;
                    ay = -ay;
                    az = -az;
                }

                s.annxf(1, is, ie) = 0.5 * ax;
                s.annxf(2, is, ie) = 0.5 * ay;
                s.annxf(3, is, ie) = 0.5 * az;
                s.annxf(4, is, ie) = 0.5 * amag;
            }
        }

        for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
        {
            const Int parent = s.iside(s.cfg.nsidpe, ib);
            const Int local_face = s.iside(s.cfg.nsidpl, ib);

            if (parent < 1 || parent > s.cfg.nelem ||
                local_face < 1 || local_face > s.cfg.nsid)
            {
                throw std::runtime_error(
                    "Preprocess::getNormals - boundary local face numbers must be assigned before boundary normal extraction");
            }

            for (Int idim = 1; idim <= s.cfg.ndim1; ++idim)
            {
                s.face_norm(idim, ib) = s.annxf(idim, local_face, parent);
            }
        }
    }

    //=========================================================================
    // Assembles nodal lumped mass, thermal capacitance and the correction
    // between lumped and consistent P1 tetrahedral mass matrices.
    //
    // The consistent scalar element mass matrix is
    //
    //     M_ij^(e) = integral(V_e) N_i N_j dV
    //
    // For a four-node linear tetrahedron:
    //
    //     M_ii^(e) = V_e / 10 = det(J_e) / 60
    //     M_ij^(e) = V_e / 20 = det(J_e) / 120,  i != j
    //
    // The code stores the lumped nodal contribution as
    //
    //     m_i^(e) = det(J_e) * mass_factor
    //
    // With the standard value mass_factor = 1/24, this gives V_e/4.
    //
    // Thermal capacitance is
    //
    //     c_i^(e) = (rho cp)_e m_i^(e)
    //
    // The arrays M_diag and Mconsist store the matrix correction
    //
    //     M_L - M_C
    //
    // used by the characteristic formulation:
    //
    //     diagonal     m_i^(e) - V_e/10
    //     off-diagonal            -V_e/20
    //
    // Output:
    //     s.Mdiag_real  assembled lumped nodal mass
    //     s.elcoe_e     element-node lumped mass contribution
    //     s.elcoe2      inverse lumped nodal mass
    //     s.elcoe2p     inverse lumped nodal thermal capacitance
    //     s.M_diag      diagonal of M_L - M_C
    //     s.Mconsist    off-diagonals of M_L - M_C
    //=========================================================================
    void Preprocess::massMatrix(CBSStateSI& s)
    {
        validate_core_dimensions(s);

        s.Mdiag_real.fill(0.0);
        s.elcoe_e.fill(0.0);
        s.elcoe2.fill(0.0);
        s.elcoe2p.fill(0.0);
        s.M_diag.fill(0.0);
        s.Mconsist.fill(0.0);

        Array1D<Real> thermal_lumped;
        thermal_lumped.resize(s.cfg.npoin);
        thermal_lumped.fill(0.0);

        // Every rank stores only owned tetrahedra. Therefore each element
        // contribution is assembled exactly once globally.
        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            if (s.detJ(ie) <= 0.0 || !std::isfinite(s.detJ(ie)))
            {
                throw std::runtime_error(
                    "Preprocess::massMatrix - invalid detJ at element "
                    + std::to_string(ie));
            }

            if (s.rho_cp_e(ie) <= 0.0 || !std::isfinite(s.rho_cp_e(ie)))
            {
                throw std::runtime_error(
                    "Preprocess::massMatrix - invalid rho*cp at element "
                    + std::to_string(ie));
            }

            const Real nodal_mass =
                s.detJ(ie) * s.cfg.mass_factor;

            const Real nodal_capacity =
                s.rho_cp_e(ie) * nodal_mass;

            const Real consistent_diag =
                s.detJ(ie) / 60.0;

            const Real correction_off =
                -s.detJ(ie) / 120.0;

            for (Int in = 1; in <= s.cfg.nep; ++in)
            {
                const Int ip = s.intma(in, ie);

                s.Mdiag_real(ip) += nodal_mass;
                thermal_lumped(ip) += nodal_capacity;
                s.M_diag(ip) += consistent_diag;

                s.elcoe_e(
                    element_node_index(s, ie, in)) = nodal_mass;
            }

            const Int first_pair =
                (ie - 1) * s.cfg.gsdim + 1;

            for (Int ig = 0; ig < s.cfg.gsdim; ++ig)
            {
                s.Mconsist(first_pair + ig) = correction_off;
            }
        }

#ifdef CBS3D_USE_MPI
        if (s.mpi_enabled)
        {
            // A rank-owned tetrahedron can contribute to a node owned by a
            // neighbouring rank. Such contributions are stored temporarily
            // in the local ghost-node entries and must be added to the owner.
            HaloExchange::sumGhostContributionsToOwners(
                s.Mdiag_real,
                s.partition_metadata);

            HaloExchange::sumGhostContributionsToOwners(
                thermal_lumped,
                s.partition_metadata);

            HaloExchange::sumGhostContributionsToOwners(
                s.M_diag,
                s.partition_metadata);

            // The owner now contains the complete shared-node coefficient.
            // Broadcast it back so all ghost copies contain identical values.
            HaloExchange::broadcastOwnedToGhosts(
                s.Mdiag_real,
                s.partition_metadata);

            HaloExchange::broadcastOwnedToGhosts(
                thermal_lumped,
                s.partition_metadata);

            HaloExchange::broadcastOwnedToGhosts(
                s.M_diag,
                s.partition_metadata);
        }
#else
        if (s.mpi_enabled)
        {
            throw std::runtime_error(
                "Preprocess::massMatrix - MPI state requires an MPI-enabled build");
        }
#endif

        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            if (s.Mdiag_real(ip) <= 0.0 ||
                !std::isfinite(s.Mdiag_real(ip)))
            {
                throw std::runtime_error(
                    "Preprocess::massMatrix - non-positive lumped mass at node "
                    + std::to_string(ip));
            }

            if (thermal_lumped(ip) <= 0.0 ||
                !std::isfinite(thermal_lumped(ip)))
            {
                throw std::runtime_error(
                    "Preprocess::massMatrix - non-positive thermal capacitance at node "
                    + std::to_string(ip));
            }

            s.elcoe2(ip) =
                1.0 / s.Mdiag_real(ip);

            s.elcoe2p(ip) =
                1.0 / thermal_lumped(ip);

            s.M_diag(ip) =
                s.Mdiag_real(ip) - s.M_diag(ip);
        }
    }

    //=========================================================================
    // Builds the persistent nodal material-connectivity mask.
    //
    // Each owned tetrahedron contributes one bit to all four of its nodes:
    //
    //     node_touches_fluid = 1   mat_elem(e) == 0
    //     node_touches_solid = 2   mat_elem(e) > 0
    //
    // A conformal fluid-solid interface node therefore has mask 3.
    //
    // In distributed memory, a shared node may see its fluid tetrahedra on
    // one rank and its solid tetrahedra on another rank. Ghost contributions
    // are therefore combined on the owner with bitwise OR and the completed
    // owner value is subsequently broadcast to all ghost copies.
    //=========================================================================
    void Preprocess::buildMaterialNodeMasks(CBSStateSI& s)
    {
        validate_core_dimensions(s);

        s.node_material_mask.fill(0);

        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            const Int material_bit =
                s.mat_elem(ie) == 0
                    ? CBSStateSI::node_touches_fluid
                    : CBSStateSI::node_touches_solid;

            for (Int in = 1; in <= s.cfg.nep; ++in)
            {
                const Int ip = s.intma(in, ie);

                if (ip < 1 || ip > s.cfg.npoin)
                {
                    throw std::runtime_error(
                        "Preprocess::buildMaterialNodeMasks - "
                        "element node is outside the local node range");
                }

                s.node_material_mask(ip) |= material_bit;
            }
        }

#ifdef CBS3D_USE_MPI
        if (s.mpi_enabled)
        {
            HaloExchange::orGhostMasksToOwners(
                s.node_material_mask,
                s.partition_metadata);

            HaloExchange::broadcastOwnedToGhosts(
                s.node_material_mask,
                s.partition_metadata);
        }
#endif

        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            const Int mask = s.node_material_mask(ip);

            if (mask < CBSStateSI::node_touches_fluid ||
                mask >
                    (CBSStateSI::node_touches_fluid |
                     CBSStateSI::node_touches_solid))
            {
                throw std::runtime_error(
                    "Preprocess::buildMaterialNodeMasks - "
                    "invalid reconciled material mask at node "
                    + std::to_string(ip));
            }
        }
    }


    //=========================================================================
    // Classifies every tetrahedral face for the CBS boundary correction terms.
    //
    //     fedge = 0   interior face
    //     fedge = 1   exterior face with strongly prescribed velocity
    //     fedge = 2   other exterior face
    //
    // Only faces listed in the boundary-face array are changed from zero.
    //=========================================================================
    void Preprocess::classifyFaceEdges(CBSStateSI& s)
    {
        validate_core_dimensions(s);
        validateBoundaryFlags(s);

        s.fedge.fill(0);

        for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
        {
            const Int parent = s.iside(s.cfg.nsidpe, ib);
            const Int local_face = s.iside(s.cfg.nsidpl, ib);
            const Int bc = s.iside(s.cfg.bsid, ib);

            if (parent < 1 || parent > s.cfg.nelem ||
                local_face < 1 || local_face > s.cfg.nsid)
            {
                throw std::runtime_error(
                    "Preprocess::classifyFaceEdges - boundary face number not assigned");
            }

            s.fedge(local_face, parent) =
                is_prescribed_velocity_face_bc(bc) ? 1 : 2;
        }
    }

    //=========================================================================
    // Calculates the characteristic length of every tetrahedral element.
    //
    // For each face f, the altitude from the opposite vertex is
    //
    //     h_f = 3 V_e / A_f
    //
    // Since det(J_e) = 6 V_e, the code evaluates
    //
    //     h_f = [det(J_e) / 2] / A_f
    //
    // The element length is the minimum of the four altitudes:
    //
    //     h_e = min_f(h_f)
    //
    // Output:
    //     s.alen_e(e)  minimum tetrahedral altitude
    //=========================================================================
    void Preprocess::elementSize(CBSStateSI& s)
    {
        validate_core_dimensions(s);

        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            if (s.detJ(ie) <= 0.0 || !std::isfinite(s.detJ(ie)))
            {
                throw std::runtime_error(
                    "Preprocess::elementSize - invalid element volume at element "
                    + std::to_string(ie));
            }

            const Real three_volume = 0.5 * s.detJ(ie);
            Real hmin = std::numeric_limits<Real>::max();

            for (Int is = 1; is <= s.cfg.nsid; ++is)
            {
                const Real face_area = s.annxf(4, is, ie);

                if (face_area <= 0.0 || !std::isfinite(face_area))
                {
                    throw std::runtime_error(
                        "Preprocess::elementSize - invalid face area at element "
                        + std::to_string(ie)
                        + ", face "
                        + std::to_string(is));
                }

                const Real h = three_volume / face_area;
                hmin = std::min(hmin, h);
            }

            if (hmin <= 0.0 || !std::isfinite(hmin))
            {
                throw std::runtime_error(
                    "Preprocess::elementSize - invalid characteristic length at element "
                    + std::to_string(ie));
            }

            s.alen_e(ie) = hmin;
        }
    }

    //=========================================================================
    // Builds the nodal no-slip list from physical wall faces and the conformal
    // fluid-solid material interface.
    //
    // Boundary-wall nodes are collected from no-slip boundary faces. Their
    // area-weighted face normals are summed and then normalised:
    //
    //     n_i = [sum_f A_f n_f] / |sum_f A_f n_f|
    //
    // A conformal CHT interface is an internal mesh surface and therefore may
    // not appear in the external boundary-face list. Interface nodes are also
    // detected from material adjacency:
    //
    //     interface node = touches fluid AND touches solid
    //
    // These nodes are added to wall_node_list so that the fluid velocity is
    // constrained to zero at the solid surface. Internal interface nodes do
    // not have one unique external wall normal, so only their no-slip status is
    // required here.
    //
    // Output:
    //     s.wall_node_list
    //     s.wall_node_norm
    //     s.cfg.npoin_wall
    //=========================================================================
    void Preprocess::wallDetermination(CBSStateSI& s)
    {
        validate_core_dimensions(s);
        validateBoundaryFlags(s);

        const Int physical_wall_bit =
            CBSStateSI::node_on_physical_wall;

        const Int material_interface_bit =
            CBSStateSI::node_on_material_interface;

        const Int valid_wall_mask =
            physical_wall_bit |
            material_interface_bit;

        s.cfg.npoin_wall = 0;
        s.wall_node_list.fill(0);
        s.wall_node_norm.fill(0.0);
        s.node_wall_mask.fill(0);
        s.node_wall_normal_sum.fill(0.0);

        // -------------------------------------------------------------
        // Rank-local physical-boundary contribution.
        //
        // Physical boundary faces occur exactly once in the distributed
        // mesh. Their nodal flags and area-weighted normal contributions
        // are initially accumulated on the rank that owns the face.
        // -------------------------------------------------------------
        for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
        {
            const Int bc = s.iside(s.cfg.bsid, ib);

            if (!is_no_slip_wall_bc(bc))
            {
                continue;
            }

            const Real area = s.face_norm(4, ib);

            if (area <= 0.0 || !std::isfinite(area))
            {
                throw std::runtime_error(
                    "Preprocess::wallDetermination - "
                    "invalid physical-wall face area");
            }

            for (Int in = 1; in <= s.cfg.nsidp; ++in)
            {
                const Int ip = s.iside(in, ib);

                if (ip < 1 || ip > s.cfg.npoin)
                {
                    throw std::runtime_error(
                        "Preprocess::wallDetermination - "
                        "physical-wall node out of range");
                }

                s.node_wall_mask(ip) |= physical_wall_bit;

                s.node_wall_normal_sum(1, ip) +=
                    s.face_norm(1, ib);

                s.node_wall_normal_sum(2, ip) +=
                    s.face_norm(2, ib);

                s.node_wall_normal_sum(3, ip) +=
                    s.face_norm(3, ib);
            }
        }

        // -------------------------------------------------------------
        // Conformal CHT interfaces are detected from the persistent
        // material-connectivity mask. This mask was reconciled before
        // wallDetermination() was called.
        // -------------------------------------------------------------
        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            if (is_conformal_fluid_solid_interface_node(s, ip))
            {
                s.node_wall_mask(ip) |= material_interface_bit;
            }
        }

#ifdef CBS3D_USE_MPI
        if (s.mpi_enabled)
        {
            // A shared node may be marked by a physical wall face or material
            // interface detected on a neighbouring rank. Combine all wall bits
            // on the owner and then broadcast the final classification.
            HaloExchange::orGhostMasksToOwners(
                s.node_wall_mask,
                s.partition_metadata);

            HaloExchange::broadcastOwnedToGhosts(
                s.node_wall_mask,
                s.partition_metadata);

            // Sum all area-weighted physical-wall normal contributions on the
            // owner and copy the complete vector back to every ghost.
            HaloExchange::sumGhostContributionsToOwners(
                s.node_wall_normal_sum,
                s.partition_metadata);

            HaloExchange::broadcastOwnedToGhosts(
                s.node_wall_normal_sum,
                s.partition_metadata);
        }
#endif

        Int material_interface_nodes_added = 0;

        // -------------------------------------------------------------
        // Rebuild the rank-local wall list after owner/ghost
        // reconciliation. Every local copy of a shared node therefore receives
        // the same wall classification and physical-wall normal.
        // -------------------------------------------------------------
        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            const Int wall_mask = s.node_wall_mask(ip);

            if (wall_mask < 0 || wall_mask > valid_wall_mask)
            {
                throw std::runtime_error(
                    "Preprocess::wallDetermination - "
                    "invalid reconciled wall mask at node "
                    + std::to_string(ip));
            }

            if (wall_mask == 0)
            {
                continue;
            }

            ++s.cfg.npoin_wall;

            const Int iw = s.cfg.npoin_wall;
            s.wall_node_list(iw) = ip;

            const bool physical_wall =
                (wall_mask & physical_wall_bit) != 0;

            const bool material_interface =
                (wall_mask & material_interface_bit) != 0;

            if (material_interface && !physical_wall)
            {
                ++material_interface_nodes_added;
            }

            const Real nx = s.node_wall_normal_sum(1, ip);
            const Real ny = s.node_wall_normal_sum(2, ip);
            const Real nz = s.node_wall_normal_sum(3, ip);

            if (!std::isfinite(nx) ||
                !std::isfinite(ny) ||
                !std::isfinite(nz))
            {
                throw std::runtime_error(
                    "Preprocess::wallDetermination - "
                    "non-finite reconciled wall normal");
            }

            const Real length =
                std::sqrt(nx * nx + ny * ny + nz * nz);

            if (physical_wall && length > 0.0)
            {
                s.wall_node_norm(1, iw) = nx / length;
                s.wall_node_norm(2, iw) = ny / length;
                s.wall_node_norm(3, iw) = nz / length;
            }
            else
            {
                // Interface-only nodes have no external boundary normal.
                // At geometric corners, physical-face contributions can also
                // cancel. Preserve the established zero-vector behaviour.
                s.wall_node_norm(1, iw) = 0.0;
                s.wall_node_norm(2, iw) = 0.0;
                s.wall_node_norm(3, iw) = 0.0;
            }
        }

        // These are meaningful global values only in a serial calculation.
        // Distributed global totals are printed
        if (!s.mpi_enabled)
        {
            std::cout
                << "Wall/interface no-slip nodes detected: "
                << s.cfg.npoin_wall << "\n"
                << "Conformal material-interface nodes added to "
                   "no-slip list: "
                << material_interface_nodes_added << "\n";
        }
    }


    //=========================================================================
    // Converts a prescribed inlet mass-flow rate into velocity magnitude.
    //
    // The total inlet area is the sum of all BC 511 face areas:
    //
    //     A_in = sum_f A_f
    //
    // Conservation of mass gives
    //
    //     m_dot = rho_in A_in U_in
    //
    // and therefore
    //
    //     U_in = m_dot / (rho_in A_in)
    //
    // This routine calculates only the magnitude. The inlet direction and
    // nodal velocity components are imposed later by Boundary::applyVelocity().
    //
    // Output:
    //     s.cfg.inlet_u_from_massflow
    //=========================================================================
    void Preprocess::computeMassFlowInletVelocity(CBSStateSI& s)
    {
        validateBoundaryFlags(s);

        if (s.cfg.mass_flow_inlet_enabled < 1)
        {
            return;
        }

        if (s.cfg.inlet_density <= 0.0 ||
            !std::isfinite(s.cfg.inlet_density))
        {
            throw std::runtime_error(
                "Preprocess::computeMassFlowInletVelocity - "
                "inlet density must be positive");
        }

        Real local_inlet_area = 0.0;
        Int local_inlet_faces = 0;

        // Each physical inlet face occurs on exactly one rank-local
        // partition. Artificial partition boundaries are not present in
        // the physical boundary-face file.
        for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
        {
            const Int bc = s.iside(s.cfg.bsid, ib);

            if (!is_mass_flow_inlet_bc(bc))
            {
                continue;
            }

            const Real area = s.face_norm(4, ib);

            if (area <= 0.0 || !std::isfinite(area))
            {
                throw std::runtime_error(
                    "Preprocess::computeMassFlowInletVelocity - "
                    "invalid inlet face area");
            }

            local_inlet_area += area;
            ++local_inlet_faces;
        }

        Real inlet_area = local_inlet_area;
        Int inlet_faces = local_inlet_faces;

#ifdef CBS3D_USE_MPI
        if (s.mpi_enabled)
        {
            const int area_error =
                MPI_Allreduce(
                    &local_inlet_area,
                    &inlet_area,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM,
                    MPI_COMM_WORLD);

            if (area_error != MPI_SUCCESS)
            {
                throw std::runtime_error(
                    "Preprocess::computeMassFlowInletVelocity - "
                    "MPI_Allreduce failed for inlet area");
            }

            const int face_error =
                MPI_Allreduce(
                    &local_inlet_faces,
                    &inlet_faces,
                    1,
                    MPI_INT,
                    MPI_SUM,
                    MPI_COMM_WORLD);

            if (face_error != MPI_SUCCESS)
            {
                throw std::runtime_error(
                    "Preprocess::computeMassFlowInletVelocity - "
                    "MPI_Allreduce failed for inlet face count");
            }
        }
#endif

        if (inlet_faces < 1)
        {
            throw std::runtime_error(
                "Preprocess::computeMassFlowInletVelocity - "
                "mass-flow inlet enabled but no BC_ID 511 faces were found");
        }

        if (inlet_area <= 0.0 ||
            !std::isfinite(inlet_area))
        {
            throw std::runtime_error(
                "Preprocess::computeMassFlowInletVelocity - "
                "global inlet area is invalid");
        }

        s.cfg.inlet_u_from_massflow =
            s.cfg.inlet_mass_flow_rate /
            (s.cfg.inlet_density * inlet_area);

        if (!std::isfinite(s.cfg.inlet_u_from_massflow))
        {
            throw std::runtime_error(
                "Preprocess::computeMassFlowInletVelocity - "
                "computed inlet velocity is not finite");
        }

        // Avoid duplicate MPI output. Every rank nevertheless stores the same
        // globally calculated inlet velocity.
        if (!s.mpi_enabled || s.mpi_rank == 0)
        {
            std::cout
                << "Mass-flow inlet faces: "
                << inlet_faces << "\n"
                << "Mass-flow inlet area : "
                << inlet_area << "\n"
                << "Mass-flow velocity magnitude: "
                << s.cfg.inlet_u_from_massflow << "\n";
        }
    }


    //=========================================================================
    // Builds the persistent distributed nodal velocity-boundary state.
    //
    // Rank-local physical boundary faces first contribute integer
    // classification bits and, for BC 511, area-weighted outward normals:
    //
    //     N_i = sum_{f incident on i} A_f n_f
    //
    // Shared-node integer classifications are reconciled with bitwise OR.
    // Shared-node normal contributions are reconciled with summation.
    //
    // The final strong velocity priority is:
    //
    //     ordinary prescribed velocity
    //         -> physical no-slip wall
    //         -> moving-wall BC 500
    //         -> material-solid/interface no-slip
    //
    // Pressure-outlet and symmetry membership are retained separately because
    // neither represents a fixed three-component velocity value.
    //=========================================================================
    void Preprocess::buildVelocityBoundaryState(CBSStateSI& s)
    {
        validate_core_dimensions(s);
        validateBoundaryFlags(s);

        // Rank-independent boundary candidate bits. These are combined on
        // shared-node owners with bitwise OR.
        constexpr Int candidate_moving_wall = 1 << 0;
        constexpr Int candidate_bc_503 = 1 << 1;
        constexpr Int candidate_bc_507 = 1 << 2;
        constexpr Int candidate_bc_508 = 1 << 3;
        constexpr Int candidate_bc_510 = 1 << 4;
        constexpr Int candidate_bc_511 = 1 << 5;
        constexpr Int candidate_pressure_outlet = 1 << 6;
        constexpr Int candidate_symmetry = 1 << 7;

        const Int ordinary_prescribed_mask =
            candidate_bc_503 |
            candidate_bc_507 |
            candidate_bc_508 |
            candidate_bc_510 |
            candidate_bc_511;

        Array1D<Int> candidate_mask;
        candidate_mask.resize(s.cfg.npoin);
        candidate_mask.fill(0);

        s.node_velocity_bc_type.fill(
            CBSStateSI::velocity_bc_free);

        s.node_velocity_bc_priority.fill(
            CBSStateSI::velocity_priority_free);

        s.node_velocity_bc_value.fill(0.0);
        s.node_inlet_normal_sum.fill(0.0);
        s.node_inlet_normal.fill(0.0);

        s.node_massflow_inlet.fill(0);
        s.node_pressure_outlet.fill(0);
        s.node_symmetry.fill(0);

        // -------------------------------------------------------------
        // Rank-local physical-boundary contributions.
        //
        // Every physical boundary triangle occurs on exactly one rank.
        // Its nodes may nevertheless be owned by another rank.
        // -------------------------------------------------------------
        for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
        {
            const Int bc = s.iside(s.cfg.bsid, ib);

            Int boundary_bit = 0;

            if (bc == s.cfg.bc_adiabatic_prescribed_velocity)
            {
                boundary_bit = candidate_moving_wall;
            }
            else if (
                bc == s.cfg.bc_temperature_zero_prescribed_velocity)
            {
                boundary_bit = candidate_bc_503;
            }
            else if (bc == s.cfg.bc_bfs_parabolic_inlet)
            {
                boundary_bit = candidate_bc_507;
            }
            else if (bc == s.cfg.bc_parabolic_inlet)
            {
                boundary_bit = candidate_bc_508;
            }
            else if (bc == s.cfg.bc_velocity_temperature_inlet)
            {
                boundary_bit = candidate_bc_510;
            }
            else if (bc == s.cfg.bc_massflow_temperature_inlet)
            {
                boundary_bit = candidate_bc_511;
            }
            else if (bc == s.cfg.bc_pressure_outlet)
            {
                boundary_bit = candidate_pressure_outlet;
            }
            else if (bc == s.cfg.bc_symmetry_no_flux)
            {
                boundary_bit = candidate_symmetry;
            }
            else
            {
                // No-slip state is already represented by node_wall_mask.
                // Material-solid state is represented by node_material_mask.
                // BC 504 and BC 902 do not impose strong velocity values.
                continue;
            }

            const bool massflow_face =
                bc == s.cfg.bc_massflow_temperature_inlet;

            if (massflow_face)
            {
                const Real area = s.face_norm(4, ib);

                if (area <= 0.0 ||
                    !std::isfinite(area) ||
                    !std::isfinite(s.face_norm(1, ib)) ||
                    !std::isfinite(s.face_norm(2, ib)) ||
                    !std::isfinite(s.face_norm(3, ib)))
                {
                    throw std::runtime_error(
                        "Preprocess::buildVelocityBoundaryState - "
                        "invalid BC 511 face normal");
                }
            }

            for (Int in = 1; in <= s.cfg.nsidp; ++in)
            {
                const Int ip = s.iside(in, ib);

                if (ip < 1 || ip > s.cfg.npoin)
                {
                    throw std::runtime_error(
                        "Preprocess::buildVelocityBoundaryState - "
                        "boundary node is outside the local node range");
                }

                candidate_mask(ip) |= boundary_bit;

                if (massflow_face)
                {
                    // face_norm(1:3,ib) already stores A_f n_f.
                    s.node_inlet_normal_sum(1, ip) +=
                        s.face_norm(1, ib);

                    s.node_inlet_normal_sum(2, ip) +=
                        s.face_norm(2, ib);

                    s.node_inlet_normal_sum(3, ip) +=
                        s.face_norm(3, ib);
                }
            }
        }

#ifdef CBS3D_USE_MPI
        if (s.mpi_enabled)
        {
            // Reconcile boundary membership on shared-node owners.
            HaloExchange::orGhostMasksToOwners(
                candidate_mask,
                s.partition_metadata);

            HaloExchange::broadcastOwnedToGhosts(
                candidate_mask,
                s.partition_metadata);

            // Reconcile all BC 511 area-weighted normal contributions.
            HaloExchange::sumGhostContributionsToOwners(
                s.node_inlet_normal_sum,
                s.partition_metadata);

            HaloExchange::broadcastOwnedToGhosts(
                s.node_inlet_normal_sum,
                s.partition_metadata);
        }
#else
        if (s.mpi_enabled)
        {
            throw std::runtime_error(
                "Preprocess::buildVelocityBoundaryState - "
                "MPI state requires an MPI-enabled build");
        }
#endif

        const auto count_bits = [](Int value) -> Int
        {
            Int count = 0;

            while (value != 0)
            {
                count += value & 1;
                value >>= 1;
            }

            return count;
        };

        const auto bfs_profile = [](Real y) -> Real
        {
            return
                0.6624 * std::pow(y, 6)
                - 7.5547 * std::pow(y, 5)
                + 33.9 * std::pow(y, 4)
                - 75.283 * std::pow(y, 3)
                + 83.368 * std::pow(y, 2)
                - 37.793 * y
                + 2.6959;
        };

        const auto rectangular_profile = [](Real y) -> Real
        {
            return 6.0 * y * (1.0 - y);
        };

        // -------------------------------------------------------------
        // Derive the complete nodal state from reconciled input masks.
        // -------------------------------------------------------------
        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            const Int candidates = candidate_mask(ip);

            const bool has_moving_wall =
                (candidates & candidate_moving_wall) != 0;

            const bool has_bc_503 =
                (candidates & candidate_bc_503) != 0;

            const bool has_bc_507 =
                (candidates & candidate_bc_507) != 0;

            const bool has_bc_508 =
                (candidates & candidate_bc_508) != 0;

            const bool has_bc_510 =
                (candidates & candidate_bc_510) != 0;

            const bool has_bc_511 =
                (candidates & candidate_bc_511) != 0;

            const bool has_pressure_outlet =
                (candidates & candidate_pressure_outlet) != 0;

            const bool has_symmetry =
                (candidates & candidate_symmetry) != 0;

            s.node_massflow_inlet(ip) =
                has_bc_511 ? 1 : 0;

            s.node_pressure_outlet(ip) =
                has_pressure_outlet ? 1 : 0;

            s.node_symmetry(ip) =
                has_symmetry ? 1 : 0;

            // ---------------------------------------------------------
            // Normalise the globally assembled BC 511 nodal normal.
            // ---------------------------------------------------------
            if (has_bc_511)
            {
                const Real nx =
                    s.node_inlet_normal_sum(1, ip);

                const Real ny =
                    s.node_inlet_normal_sum(2, ip);

                const Real nz =
                    s.node_inlet_normal_sum(3, ip);

                const Real magnitude =
                    std::sqrt(
                        nx * nx +
                        ny * ny +
                        nz * nz);

                if (!std::isfinite(magnitude) ||
                    magnitude <= 1.0e-14)
                {
                    const Size local_index =
                        static_cast<Size>(ip);

                    long long global_node = -1;
                    Int owner_rank = -1;

                    if (local_index <
                        s.partition_metadata.local_to_global_node.size())
                    {
                        global_node = static_cast<long long>(
                            s.partition_metadata
                                .local_to_global_node[local_index]);
                    }

                    if (local_index <
                        s.partition_metadata.node_owner_rank.size())
                    {
                        owner_rank =
                            s.partition_metadata
                                .node_owner_rank[local_index];
                    }

                    std::cerr
                        << std::scientific
                        << std::setprecision(17)
                        << "DD-2D1A-2 BC 511 normal diagnostic:"
                        << " rank=" << s.mpi_rank
                        << " local_node=" << ip
                        << " global_node=" << global_node
                        << " owner_rank=" << owner_rank
                        << " candidate_mask=" << candidates
                        << " normal_sum=("
                        << nx << ","
                        << ny << ","
                        << nz << ")"
                        << " magnitude=" << magnitude
                        << " threshold=1.00000000000000000e-14"
                        << "\n";

                    throw std::runtime_error(
                        "Preprocess::buildVelocityBoundaryState - "
                        "BC 511 node has a zero or invalid resultant normal");
                }

                s.node_inlet_normal(1, ip) = nx / magnitude;
                s.node_inlet_normal(2, ip) = ny / magnitude;
                s.node_inlet_normal(3, ip) = nz / magnitude;
            }

            const Int material_mask =
                s.node_material_mask(ip);

            const Int valid_material_mask =
                CBSStateSI::node_touches_fluid |
                CBSStateSI::node_touches_solid;

            if (material_mask < CBSStateSI::node_touches_fluid ||
                material_mask > valid_material_mask)
            {
                throw std::runtime_error(
                    "Preprocess::buildVelocityBoundaryState - "
                    "invalid material node mask");
            }

            const Int wall_mask =
                s.node_wall_mask(ip);

            const Int valid_wall_mask =
                CBSStateSI::node_on_physical_wall |
                CBSStateSI::node_on_material_interface;

            if (wall_mask < 0 || wall_mask > valid_wall_mask)
            {
                throw std::runtime_error(
                    "Preprocess::buildVelocityBoundaryState - "
                    "invalid wall node mask");
            }

            const bool touches_solid =
                (material_mask &
                 CBSStateSI::node_touches_solid) != 0;

            const bool physical_wall =
                (wall_mask &
                 CBSStateSI::node_on_physical_wall) != 0;

            const Int ordinary_count =
                count_bits(
                    candidates &
                    ordinary_prescribed_mask);

            // ---------------------------------------------------------
            // Priority 4: solid-only and fluid-solid interface nodes.
            // ---------------------------------------------------------
            if (touches_solid)
            {
                s.node_velocity_bc_type(ip) =
                    CBSStateSI::velocity_bc_noslip;

                s.node_velocity_bc_priority(ip) =
                    CBSStateSI::velocity_priority_material_solid;
            }

            // ---------------------------------------------------------
            // Priority 3: moving wall BC 500.

            else if (has_moving_wall)
            {
                s.node_velocity_bc_type(ip) =
                    CBSStateSI::velocity_bc_moving_wall;

                s.node_velocity_bc_priority(ip) =
                    CBSStateSI::velocity_priority_moving_wall;

                s.node_velocity_bc_value(1, ip) =
                    s.cfg.inlet_u;

                s.node_velocity_bc_value(2, ip) =
                    s.cfg.inlet_v;

                s.node_velocity_bc_value(3, ip) =
                    s.cfg.inlet_w;
            }

            // ---------------------------------------------------------
            // Priority 2: physical no-slip wall.
            // ---------------------------------------------------------
            else if (physical_wall)
            {
                s.node_velocity_bc_type(ip) =
                    CBSStateSI::velocity_bc_noslip;

                s.node_velocity_bc_priority(ip) =
                    CBSStateSI::velocity_priority_physical_wall;
            }

            // ---------------------------------------------------------
            // Priority 1: ordinary prescribed velocity.
            // ---------------------------------------------------------
            else if (ordinary_count > 0)
            {
                if (ordinary_count != 1)
                {
                    throw std::runtime_error(
                        "Preprocess::buildVelocityBoundaryState - "
                        "multiple ordinary prescribed-velocity families "
                        "remain active at one unconstrained node");
                }

                s.node_velocity_bc_type(ip) =
                    CBSStateSI::velocity_bc_prescribed;

                s.node_velocity_bc_priority(ip) =
                    CBSStateSI::velocity_priority_prescribed;

                if (has_bc_503)
                {
                    s.node_velocity_bc_value(1, ip) = 1.0;
                }
                else if (has_bc_507)
                {
                    s.node_velocity_bc_value(1, ip) =
                        bfs_profile(s.coord(2, ip));
                }
                else if (has_bc_508)
                {
                    s.node_velocity_bc_value(1, ip) =
                        rectangular_profile(s.coord(2, ip));
                }
                else if (has_bc_510)
                {
                    s.node_velocity_bc_value(1, ip) =
                        s.cfg.inlet_u;

                    s.node_velocity_bc_value(2, ip) =
                        s.cfg.inlet_v;

                    s.node_velocity_bc_value(3, ip) =
                        s.cfg.inlet_w;
                }
                else if (has_bc_511)
                {
                    if (s.cfg.mass_flow_inlet_enabled > 0)
                    {
                        s.node_velocity_bc_value(1, ip) =
                            -s.cfg.inlet_u_from_massflow *
                            s.node_inlet_normal(1, ip);

                        s.node_velocity_bc_value(2, ip) =
                            -s.cfg.inlet_u_from_massflow *
                            s.node_inlet_normal(2, ip);

                        s.node_velocity_bc_value(3, ip) =
                            -s.cfg.inlet_u_from_massflow *
                            s.node_inlet_normal(3, ip);
                    }
                    else
                    {
                        s.node_velocity_bc_value(1, ip) =
                            s.cfg.inlet_u;

                        s.node_velocity_bc_value(2, ip) =
                            s.cfg.inlet_v;

                        s.node_velocity_bc_value(3, ip) =
                            s.cfg.inlet_w;
                    }
                }
            }

            for (Int idim = 1;
                 idim <= s.cfg.ndim;
                 ++idim)
            {
                if (!std::isfinite(
                        s.node_velocity_bc_value(idim, ip)) ||
                    !std::isfinite(
                        s.node_inlet_normal(idim, ip)) ||
                    !std::isfinite(
                        s.node_inlet_normal_sum(idim, ip)))
                {
                    throw std::runtime_error(
                        "Preprocess::buildVelocityBoundaryState - "
                        "non-finite persistent nodal boundary state");
                }
            }
        }

#ifdef CBS3D_USE_MPI
        if (s.mpi_enabled)
        {
            // The owner is authoritative for the completed persistent state.
            HaloExchange::broadcastOwnedToGhosts(
                s.node_velocity_bc_type,
                s.partition_metadata);

            HaloExchange::broadcastOwnedToGhosts(
                s.node_velocity_bc_priority,
                s.partition_metadata);

            HaloExchange::broadcastOwnedToGhosts(
                s.node_velocity_bc_value,
                s.partition_metadata);

            HaloExchange::broadcastOwnedToGhosts(
                s.node_inlet_normal,
                s.partition_metadata);

            HaloExchange::broadcastOwnedToGhosts(
                s.node_massflow_inlet,
                s.partition_metadata);

            HaloExchange::broadcastOwnedToGhosts(
                s.node_pressure_outlet,
                s.partition_metadata);

            HaloExchange::broadcastOwnedToGhosts(
                s.node_symmetry,
                s.partition_metadata);
        }
#endif
    }


    //=========================================================================
    // Initialises the scalar nodal velocity magnitude:
    //
    //     |u_i| = sqrt(u_i^2 + v_i^2 + w_i^2 + epsilon)
    //
    // where epsilon = 10^(-16) prevents an exactly zero square-root argument.
    // The same value is copied to velocity_old for the first residual update.
    //=========================================================================
    void Preprocess::initialiseVelocityMagnitude(CBSStateSI& s)
    {
        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            s.velocity(ip) = std::sqrt(
                s.unkno(1, ip) * s.unkno(1, ip) +
                s.unkno(2, ip) * s.unkno(2, ip) +
                s.unkno(3, ip) * s.unkno(3, ip) +
                1.0e-16);

            s.velocity_old(ip) = s.velocity(ip);
        }
    }

    //=========================================================================
    // Builds the list of nodes with prescribed pressure.
    //
    // All nodes on BC 520 pressure-outlet faces are added once, provided they
    // are connected to the fluid domain. Their prescribed value is
    //
    //     p_i = outlet_pressure_gauge
    //
    // If the mesh contains no pressure-outlet face, one fluid-connected node
    // is selected to remove the constant-pressure null space. The requested
    // pnode is used when valid; otherwise the first fluid-connected node is
    // selected.
    //
    // Output:
    //     s.bc_list
    //     s.bc_values
    //     s.cfg.bc_fixed
    //=========================================================================
    void Preprocess::detectPressureBoundaryNodes(CBSStateSI& s)
    {
        validateBoundaryFlags(s);

        s.node_pressure_fixed.fill(0);
        s.cfg.bc_fixed = 0;
        s.bc_list.fill(0);
        s.bc_values.fill(0.0);

        // -------------------------------------------------------------
        // Mark pressure-outlet nodes from this rank's physical BC 520
        // boundary faces.
        // -------------------------------------------------------------
        for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
        {
            const Int bc = s.iside(s.cfg.bsid, ib);

            if (!is_pressure_outlet_bc(bc))
            {
                continue;
            }

            for (Int in = 1; in <= s.cfg.nsidp; ++in)
            {
                const Int ip = s.iside(in, ib);

                if (ip < 1 || ip > s.cfg.npoin)
                {
                    throw std::runtime_error(
                        "Preprocess::detectPressureBoundaryNodes - "
                        "pressure node out of range");
                }

                if (!touches_fluid_domain(s, ip))
                {
                    throw std::runtime_error(
                        "Preprocess::detectPressureBoundaryNodes - "
                        "pressure outlet node "
                        + std::to_string(ip)
                        + " is not connected to a fluid element");
                }

                s.node_pressure_fixed(ip) = 1;
            }
        }

#ifdef CBS3D_USE_MPI
        if (s.mpi_enabled)
        {
            // A shared pressure-outlet node may be visible on only the rank
            // carrying the physical boundary face. Combine flags on the owner
            // and broadcast the final value to every ghost copy.
            HaloExchange::orGhostMasksToOwners(
                s.node_pressure_fixed,
                s.partition_metadata);

            HaloExchange::broadcastOwnedToGhosts(
                s.node_pressure_fixed,
                s.partition_metadata);
        }
#endif

        Int global_owned_fixed_nodes = 0;
        bool fallback_used = false;
        Int fallback_global_node = 0;

#ifdef CBS3D_USE_MPI
        if (s.mpi_enabled)
        {
            Int local_owned_fixed_nodes = 0;

            for (const Int ip : s.owned_nodes)
            {
                if (s.node_pressure_fixed(ip) != 0)
                {
                    ++local_owned_fixed_nodes;
                }
            }

            const int count_error =
                MPI_Allreduce(
                    &local_owned_fixed_nodes,
                    &global_owned_fixed_nodes,
                    1,
                    MPI_INT,
                    MPI_SUM,
                    MPI_COMM_WORLD);

            if (count_error != MPI_SUCCESS)
            {
                throw std::runtime_error(
                    "Preprocess::detectPressureBoundaryNodes - "
                    "MPI_Allreduce failed for pressure-fixed count");
            }

            // ---------------------------------------------------------
            // No explicit outlet exists. Select one deterministic global
            // fluid-connected reference node.
            //
            // In distributed input, cfg.pnode is interpreted as a global
            // mesh-node ID. If it is invalid, use the minimum global ID
            // among all owned fluid-connected nodes.
            // ---------------------------------------------------------
            if (global_owned_fixed_nodes < 1)
            {
                fallback_used = true;

                Int requested_min = 0;
                Int requested_max = 0;
                const Int requested_global_node = s.cfg.pnode;

                const int requested_min_error =
                    MPI_Allreduce(
                        &requested_global_node,
                        &requested_min,
                        1,
                        MPI_INT,
                        MPI_MIN,
                        MPI_COMM_WORLD);

                const int requested_max_error =
                    MPI_Allreduce(
                        &requested_global_node,
                        &requested_max,
                        1,
                        MPI_INT,
                        MPI_MAX,
                        MPI_COMM_WORLD);

                if (requested_min_error != MPI_SUCCESS ||
                    requested_max_error != MPI_SUCCESS)
                {
                    throw std::runtime_error(
                        "Preprocess::detectPressureBoundaryNodes - "
                        "MPI_Allreduce failed for requested pressure node");
                }

                if (requested_min != requested_max)
                {
                    throw std::runtime_error(
                        "Preprocess::detectPressureBoundaryNodes - "
                        "inconsistent pnode values across MPI ranks");
                }

                const Int global_npoin =
                    static_cast<Int>(
                        s.partition_metadata.global_npoin);

                Int local_requested_valid = 0;

                if (requested_global_node >= 1 &&
                    requested_global_node <= global_npoin)
                {
                    for (const Int ip : s.owned_nodes)
                    {
                        const Size local_index =
                            static_cast<Size>(ip);

                        if (local_index >=
                            s.local_to_global_node.size())
                        {
                            throw std::runtime_error(
                                "Preprocess::detectPressureBoundaryNodes - "
                                "incomplete local-to-global node map");
                        }

                        if (s.local_to_global_node[local_index] ==
                                requested_global_node &&
                            touches_fluid_domain(s, ip))
                        {
                            local_requested_valid = 1;
                            break;
                        }
                    }
                }

                Int global_requested_valid = 0;

                const int requested_valid_error =
                    MPI_Allreduce(
                        &local_requested_valid,
                        &global_requested_valid,
                        1,
                        MPI_INT,
                        MPI_MAX,
                        MPI_COMM_WORLD);

                if (requested_valid_error != MPI_SUCCESS)
                {
                    throw std::runtime_error(
                        "Preprocess::detectPressureBoundaryNodes - "
                        "MPI_Allreduce failed for pnode validation");
                }

                if (global_requested_valid != 0)
                {
                    fallback_global_node =
                        requested_global_node;
                }
                else
                {
                    Int local_minimum_global =
                        std::numeric_limits<Int>::max();

                    for (const Int ip : s.owned_nodes)
                    {
                        if (!touches_fluid_domain(s, ip))
                        {
                            continue;
                        }

                        const Size local_index =
                            static_cast<Size>(ip);

                        if (local_index >=
                            s.local_to_global_node.size())
                        {
                            throw std::runtime_error(
                                "Preprocess::detectPressureBoundaryNodes - "
                                "incomplete local-to-global node map");
                        }

                        local_minimum_global =
                            std::min(
                                local_minimum_global,
                                s.local_to_global_node[local_index]);
                    }

                    const int fallback_error =
                        MPI_Allreduce(
                            &local_minimum_global,
                            &fallback_global_node,
                            1,
                            MPI_INT,
                            MPI_MIN,
                            MPI_COMM_WORLD);

                    if (fallback_error != MPI_SUCCESS)
                    {
                        throw std::runtime_error(
                            "Preprocess::detectPressureBoundaryNodes - "
                            "MPI_Allreduce failed for fallback node");
                    }
                }

                if (fallback_global_node < 1 ||
                    fallback_global_node > global_npoin)
                {
                    throw std::runtime_error(
                        "Preprocess::detectPressureBoundaryNodes - "
                        "no global fluid-connected pressure reference exists");
                }

                if (static_cast<Size>(fallback_global_node) <
                    s.global_to_local_node.size())
                {
                    const Int local_fallback_node =
                        s.global_to_local_node[
                            static_cast<Size>(fallback_global_node)];

                    if (local_fallback_node >= 1 &&
                        local_fallback_node <= s.cfg.npoin)
                    {
                        s.node_pressure_fixed(
                            local_fallback_node) = 1;
                    }
                }

                HaloExchange::orGhostMasksToOwners(
                    s.node_pressure_fixed,
                    s.partition_metadata);

                HaloExchange::broadcastOwnedToGhosts(
                    s.node_pressure_fixed,
                    s.partition_metadata);

                local_owned_fixed_nodes = 0;

                for (const Int ip : s.owned_nodes)
                {
                    if (s.node_pressure_fixed(ip) != 0)
                    {
                        ++local_owned_fixed_nodes;
                    }
                }

                const int fallback_count_error =
                    MPI_Allreduce(
                        &local_owned_fixed_nodes,
                        &global_owned_fixed_nodes,
                        1,
                        MPI_INT,
                        MPI_SUM,
                        MPI_COMM_WORLD);

                if (fallback_count_error != MPI_SUCCESS)
                {
                    throw std::runtime_error(
                        "Preprocess::detectPressureBoundaryNodes - "
                        "MPI_Allreduce failed after fallback selection");
                }

                if (global_owned_fixed_nodes != 1)
                {
                    throw std::runtime_error(
                        "Preprocess::detectPressureBoundaryNodes - "
                        "distributed fallback did not select exactly one "
                        "owned pressure node");
                }
            }
        }
        else
#endif
        {
            // Serial behaviour: use the requested local node, otherwise the
            // first fluid-connected local node.
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                if (s.node_pressure_fixed(ip) != 0)
                {
                    ++global_owned_fixed_nodes;
                }
            }

            if (global_owned_fixed_nodes < 1)
            {
                fallback_used = true;

                Int ip = s.cfg.pnode;

                if (ip < 1 ||
                    ip > s.cfg.npoin ||
                    !touches_fluid_domain(s, ip))
                {
                    ip = 0;

                    for (Int candidate = 1;
                         candidate <= s.cfg.npoin;
                         ++candidate)
                    {
                        if (touches_fluid_domain(s, candidate))
                        {
                            ip = candidate;
                            break;
                        }
                    }
                }

                if (ip < 1 || ip > s.cfg.npoin)
                {
                    throw std::runtime_error(
                        "Preprocess::detectPressureBoundaryNodes - "
                        "no fluid-connected pressure reference exists");
                }

                s.node_pressure_fixed(ip) = 1;
                global_owned_fixed_nodes = 1;
                fallback_global_node = ip;
            }
        }

        // -------------------------------------------------------------
        // Rebuild the rank-local prescribed-pressure list after
        // reconciliation. Shared nodes appear in the list on every rank
        // carrying a local copy.
        // -------------------------------------------------------------
        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            const Int fixed_flag =
                s.node_pressure_fixed(ip);

            if (fixed_flag != 0 && fixed_flag != 1)
            {
                throw std::runtime_error(
                    "Preprocess::detectPressureBoundaryNodes - "
                    "invalid reconciled pressure-fixed flag");
            }

            if (fixed_flag == 0)
            {
                continue;
            }

            ++s.cfg.bc_fixed;
            s.bc_list(s.cfg.bc_fixed) = ip;
            s.bc_values(s.cfg.bc_fixed) =
                s.cfg.outlet_pressure_gauge;
        }

        if (!s.mpi_enabled || s.mpi_rank == 0)
        {
            std::cout
                << "Global pressure outlet/fixed nodes: "
                << global_owned_fixed_nodes << "\n";

            if (fallback_used)
            {
                std::cout
                    << "Pressure reference fallback global node: "
                    << fallback_global_node << "\n";
            }
            else
            {
                std::cout
                    << "Pressure reference fallback: NOT USED\n";
            }
        }
    }


}
