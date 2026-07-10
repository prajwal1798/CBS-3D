#pragma once

//=============================================================================
// CBS3D++_SI
//
// Complete numerical state of the three-dimensional semi-implicit CBS solver.
//
// The structure stores:
//
//     mesh connectivity and coordinates
//     tetrahedral and boundary-face topology
//     velocity, pressure and temperature fields
//     element geometry and mass coefficients
//     momentum, pressure and energy residuals
//     material properties for conjugate heat transfer
//     Pressure CG and banded-solver work arrays
//     OpenMP element-colouring data
//     convergence and output bookkeeping
//
// All Array1D, Array2D and Array3D objects use one-based scientific indexing.
//
// Important time-level convention used by Solver::advanceOneStep():
//
//     unkno, pres, temperature
//         current solution being advanced
//
//     unkn1, pres1, temperature1
//         copy of the solution at the beginning of the current CBS iteration
//
// Therefore:
//
//     delta_u = unkno - unkn1
//     delta_p = pres - pres1
//     delta_T = temperature - temperature1
//
// Element material convention:
//
//     mat_elem(e) = 0      fluid tetrahedron
//     mat_elem(e) > 0      solid/material tetrahedron
//=============================================================================

#include "cbs/core/Array.hpp"
#include "cbs/core/RunConfig.hpp"
#include "cbs/core/Types.hpp"
#include "cbs/parallel/PartitionMetadata.hpp"

#include <array>
#include <vector>

namespace cbs
{
    struct CBSStateSI
    {
        // --------------------------------------------------------------------
        // Problem configuration
        // --------------------------------------------------------------------
        RunConfig cfg;

        // --------------------------------------------------------------------
        // Mesh and local topology
        // --------------------------------------------------------------------
        Array2D<Int> intma;       // Element connectivity: (nep, nelem)
        Array2D<Int> iside;       // Boundary-face data: (bsid, nboun)
        Array2D<Int> ippn1;       // Local face-to-node map: (nsid, nsidp)
        Array2D<Int> GiD;         // Off-diagonal node pairs: (2, gsdim)
        Array2D<Int> flag_list;   // Boundary flag to BC-number map: (2, nflag)

        Array1D<Int> wall_node_list;
        Array1D<Int> bc_list;

        Array2D<Real> coord;      // Nodal coordinates: (ndim, npoin)
        Array1D<Int> mat_elem;    // Element material: 0 fluid, greater than 0 solid

        // --------------------------------------------------------------------
        // Primary solution variables
        // --------------------------------------------------------------------
        Array2D<Real> unkno;      // Current velocity: (ndim, npoin)
        Array2D<Real> unkn1;      // Velocity copied at iteration start

        Array1D<Real> pres;       // Current pressure
        Array1D<Real> pres1;      // Pressure copied at iteration start

        Array1D<Real> temperature;   // Current nodal temperature
        Array1D<Real> temperature1;  // Temperature copied at iteration start

        // --------------------------------------------------------------------
        // Previous physical-time solutions
        // --------------------------------------------------------------------
        Array2D<Real> unknn1;
        Array2D<Real> unknn2;

        Array1D<Real> tempert1;
        Array1D<Real> tempert2;

        // --------------------------------------------------------------------
        // Temporary interpolated solutions
        // --------------------------------------------------------------------
        Array2D<Real> unkno_temp;
        Array1D<Real> pres_temp;
        Array1D<Real> temperature_temp;

        // --------------------------------------------------------------------
        // Element geometry and mass data
        // --------------------------------------------------------------------
        // Flattened shape-function gradients:
        //
        //     dN_a/dx_j,  a = 1..4,  j = 1..3
        Array1D<Real> dNkdx;
        Array1D<Real> detJ;       // Element Jacobian determinant = 6 * volume
        Array1D<Real> delte;      // Element pseudo-time step
        // Lumped element nodal mass:
        //
        //     m_a^(e) = V_e/4 = det(J_e)/24
        Array1D<Real> elcoe_e;
        // Inverse assembled time diagonals:
        //
        //     elcoe2  = [M_u/dt]^-1
        //     elcoe2p = [C_T/dt]^-1
        Array1D<Real> elcoe2;
        Array1D<Real> elcoe2p;
        Array1D<Real> alen_e;     // Element characteristic length

        // --------------------------------------------------------------------
        // Nodal time-step and stabilisation data
        // --------------------------------------------------------------------
        Array1D<Real> deltp;
        Array1D<Real> deltp1;
        Array1D<Real> deltp2;
        Array1D<Real> beta1;
        Array1D<Real> alen;
        Array1D<Real> vvis;
        Array1D<Real> velcp;
        Array1D<Int> icflag;
        Array2D<Real> unkno_unit_vec;

        // --------------------------------------------------------------------
        // Boundary-face and wall geometry
        // --------------------------------------------------------------------
        Array2D<Real> face_norm;       // Boundary normal components and face area
        Array2D<Real> wall_node_norm;  // Averaged wall normal at each node

        // --------------------------------------------------------------------
        // Nodal diagnostics and coefficients
        // --------------------------------------------------------------------
        Array1D<Real> velocity;
        Array1D<Real> velocity_old;
        Array1D<Real> Mdiag_real;
        Array1D<Real> beta;

        // --------------------------------------------------------------------
        // Three-dimensional element-face data
        // --------------------------------------------------------------------
        // Area-weighted outward element-face normals:
        //
        //     annxf(1:3, face, element) = A_f n_f
        //     annxf(4,   face, element) = A_f
        Array3D<Real> annxf;
        Array2D<Int> fedge;       // Face type: 0 interior, 1 velocity BC, 2 exterior
        Array1D<Real> Mconsist;   // Consistent-mass off-diagonal terms
        Array1D<Real> M_diag;     // Consistent-mass diagonal terms

        // --------------------------------------------------------------------
        // Element material properties for conjugate heat transfer
        // --------------------------------------------------------------------
        Array1D<Real> rho_e;
        Array1D<Real> cp_e;
        Array1D<Real> k_e;
        Array1D<Real> mu_e;
        // Derived thermal properties:
        //
        //     rho_cp_e = rho_e cp_e
        //     alpha_e  = k_e / (rho_e cp_e)
        Array1D<Real> rho_cp_e;
        Array1D<Real> alpha_e;
        Array1D<Real> Qvol_e;


        // --------------------------------------------------------------------
        // Spalart-Allmaras turbulence state
        // --------------------------------------------------------------------
        // nu_tilde is the transported SA working variable.  nu_t and mu_t are
        // derived quantities and must never be prescribed independently.
        Array1D<Real> nu_tilde;
        Array1D<Real> nu_tilde1;
        Array1D<Real> nu_t;
        Array1D<Real> mu_t;

        // True minimum distance from a fluid node to the nearest physical
        // no-slip wall triangle.  This is geometric data, not a mesh-line
        // index or nearest-wall-node approximation.
        Array1D<Real> wall_distance;

        // SA nodal assembly and diagnostics.
        Array1D<Real> sa_rhs;
        Array1D<Real> sa_source;
        Array1D<Real> sa_production;
        Array1D<Real> sa_destruction;
        Array1D<Real> sa_diffusion;
        Array1D<Real> sa_residual;

        // Nodal turbulence classification flags.
        Array1D<Int> sa_active_node;
        Array1D<Int> sa_wall_node;
        Array1D<Int> sa_inlet_node;

        // Element-averaged turbulence quantities and effective properties used
        // by momentum and energy assemblies.  Molecular values remain stored in
        // mu_e and k_e and are not overwritten.
        Array1D<Real> nu_tilde_e;
        Array1D<Real> nu_t_e;
        Array1D<Real> mu_t_e;
        Array1D<Real> mu_eff_e;
        Array1D<Real> k_eff_e;

        // --------------------------------------------------------------------
        // Right-hand sides and correction arrays
        // --------------------------------------------------------------------
        Array2D<Real> rhs;
        Array1D<Real> rhs1;
        Array1D<Real> bc_values;

        Array2D<Real> Dunknn1;
        Array1D<Real> Dtempert1;

        Array2D<Real> Dudtau;
        Array1D<Real> DTdtau;

        // --------------------------------------------------------------------
        // Pressure-system coefficients
        // --------------------------------------------------------------------
        // Element pressure stiffness coefficients:
        //
        //     H_ab^(e) = V_e grad(N_a) . grad(N_b)
        Array1D<Real> pdiagE;
        Array1D<Real> gstifE;

        // Time-step-scaled pressure operator used by the linear solver.
        Array1D<Real> pdiag;
        Array1D<Real> gstif;

        // --------------------------------------------------------------------
        // Conjugate Gradient work arrays
        // --------------------------------------------------------------------
        Array1D<Real> apold;
        Array1D<Real> press;
        Array1D<Real> residual;
        Array1D<Real> Mmatrix;
        Array1D<Real> Zip;
        Array1D<Real> pvect;

        // --------------------------------------------------------------------
        // Banded-solver work arrays
        // --------------------------------------------------------------------
        Array2D<Real> gsm;
        Array2D<Real> gasfac;
        Array2D<Real> gsm1;
        Array2D<Real> gsminv;

        Array1D<Real> gaspiv;
        Array1D<Real> glm;

        // --------------------------------------------------------------------
        // Heat-transfer diagnostics
        // --------------------------------------------------------------------
        Array1D<Real> nusselt_local;
        Array1D<Real> nusselt_point_list;

        // Convergence measures:
        //
        //     hb[0:2]    u residual, field norm, RHS-scale norm
        //     hb[3:5]    v residual, field norm, RHS-scale norm
        //     hb[6:8]    w residual, field norm, RHS-scale norm
        //     hb[9:11]   pressure residual, field norm, absolute update
        //     hb[12:14]  temperature residual, field norm, RHS-scale norm
        std::array<Real, 15> hb{};

        // --------------------------------------------------------------------
        // Solver bookkeeping
        // --------------------------------------------------------------------
        Int elemConSize = 0;
        Int iiter = 0;
        Int iiter_total = 0;

        // --------------------------------------------------------------------
        // MPI and domain-decomposition metadata
        // --------------------------------------------------------------------
        // These fields are intentionally non-invasive while the validated
        // serial/OpenMP solver still owns the complete global arrays.
        //
        // mpi_rank
        //     Rank of the current process in the solver communicator.
        //
        // mpi_size
        //     Number of processes in the solver communicator.
        //
        // mpi_enabled
        //     True when the calculation uses more than one MPI rank.
        //
        // The ownership and global/local maps will be populated by the
        // domain-decomposition stage. Global mesh numbering remains one-based.
        Int mpi_rank = 0;
        Int mpi_size = 1;
        bool mpi_enabled = false;

        // Rank-local mesh ownership, global numbering and neighbour maps
        // read from the .mpi partition file.
        PartitionMetadata partition_metadata;

        std::vector<Int> owned_elements;
        std::vector<Int> ghost_elements;

        std::vector<Int> owned_nodes;
        std::vector<Int> ghost_nodes;

        // One-based global element/node id -> owning MPI rank.
        std::vector<Int> element_owner_rank;
        std::vector<Int> node_owner_rank;

        // One-based global node id -> compact local id.
        //
        // A value of -1 means that the global node is absent on this rank.
        std::vector<Int> global_to_local_node;

        // Compact local id -> one-based global node/element id.
        std::vector<Int> local_to_global_node;
        std::vector<Int> local_to_global_element;

        // Pressure-active ownership and halo maps.
        std::vector<Int> owned_pressure_nodes;
        std::vector<Int> ghost_pressure_nodes;
        std::vector<Int> pressure_global_to_local;

        // Element colouring prevents simultaneous OpenMP assembly to the same
        // node. Elements belonging to one colour do not share any node.
        Int ncolor = 0;
        std::vector<Int> color_ptr;
        std::vector<Int> color_elem;

        // Pressure colouring contains only pressure-active fluid elements.
        Int pressure_ncolor = 0;
        Int pressure_nelem = 0;
        std::vector<Int> pressure_color_ptr;
        std::vector<Int> pressure_color_elem;

        Real ani = 0.0;
        Real deltr6 = 0.0;
        Real deltr116 = 0.0;
        Real bmax = 0.0;
        Real bmax0 = 0.0;
        Real bmax1 = 0.0;

        // Last pressure-solver convergence information.
        Int last_cg_iterations = 0;
        Real last_cg_initial_l2 = 0.0;
        Real last_cg_final_l2 = 0.0;
        Real last_cg_relative_l2 = 0.0;
        Real last_cg_max_abs = 0.0;

        // Requested solution-output iterations and physical times.
        std::vector<Int> output_iterations;
        std::vector<Real> output_times;
        Real next_vtu_output_time = 0.0;

        //=====================================================================
        // Defines the local topology of the four-node tetrahedral element.
        //
        // ippn1 stores the three local nodes on each face. Face f is opposite
        // local node f:
        //
        //     face 1: local nodes 3, 2, 4
        //     face 2: local nodes 3, 4, 1
        //     face 3: local nodes 1, 4, 2
        //     face 4: local nodes 1, 2, 3
        //
        // GiD stores the six upper-triangular off-diagonal pairs of a symmetric
        // 4 x 4 element matrix:
        //
        //     (1,2), (1,3), (1,4), (2,3), (2,4), (3,4)
        //=====================================================================
        void initialise_local_topology()
        {
            cfg.bsid = cfg.nsidp + 3;
            cfg.nsidpl = cfg.nsidp + 1;
            cfg.nsidpe = cfg.nsidp + 2;

            // Three local nodes belonging to each tetrahedral face.
            ippn1.resize(cfg.nsid, cfg.nsidp);

            ippn1(1, 1) = 3;
            ippn1(1, 2) = 2;
            ippn1(1, 3) = 4;

            ippn1(2, 1) = 3;
            ippn1(2, 2) = 4;
            ippn1(2, 3) = 1;

            ippn1(3, 1) = 1;
            ippn1(3, 2) = 4;
            ippn1(3, 3) = 2;

            ippn1(4, 1) = 1;
            ippn1(4, 2) = 2;
            ippn1(4, 3) = 3;

            // Off-diagonal node pairs of a symmetric 4 x 4 element matrix.
            GiD.resize(2, cfg.gsdim);

            GiD(1, 1) = 1;
            GiD(2, 1) = 2;

            GiD(1, 2) = 1;
            GiD(2, 2) = 3;

            GiD(1, 3) = 1;
            GiD(2, 3) = 4;

            GiD(1, 4) = 2;
            GiD(2, 4) = 3;

            GiD(1, 5) = 2;
            GiD(2, 5) = 4;

            GiD(1, 6) = 3;
            GiD(2, 6) = 4;
        }

        //=====================================================================
        // Stores the problem sizes and allocates every mesh-dependent array.
        //
        // Inputs:
        //
        //     nelem_in   number of tetrahedral elements
        //     npoin_in   number of mesh nodes
        //     nboun_in   number of triangular boundary faces
        //     nflag_in   number of boundary-flag mappings
        //
        // Allocation convention:
        //
        //     element arrays       size nelem
        //     nodal arrays         size npoin
        //     boundary arrays      size nboun
        //     connectivity         (nep, nelem)
        //     coordinates          (ndim, npoin)
        //
        // Default material values describe a single fluid with:
        //
        //     rho = 1, cp = 1, k = 1, mu = 0
        //
        // until case-specific files overwrite them.
        //=====================================================================
        void set_problem_sizes(
            Int nelem_in,
            Int npoin_in,
            Int nboun_in,
            Int nflag_in)
        {
            cfg.nelem = nelem_in;
            cfg.npoin = npoin_in;
            cfg.nboun = nboun_in;
            cfg.nflag = nflag_in;

            elemConSize = cfg.nep * cfg.nelem;

            // Mesh and boundary data.
            intma.resize(cfg.nep, cfg.nelem);
            iside.resize(cfg.bsid, cfg.nboun);
            flag_list.resize(2, cfg.nflag);
            coord.resize(cfg.ndim, cfg.npoin);

            // Primary and previous solution fields.
            unkno.resize(cfg.ndim, cfg.npoin);
            unkn1.resize(cfg.ndim, cfg.npoin);
            pres.resize(cfg.npoin);
            pres1.resize(cfg.npoin);
            temperature.resize(cfg.npoin);
            temperature1.resize(cfg.npoin);

            unknn1.resize(cfg.ndim, cfg.npoin);
            unknn2.resize(cfg.ndim, cfg.npoin);
            tempert1.resize(cfg.npoin);
            tempert2.resize(cfg.npoin);

            unkno_temp.resize(cfg.ndim, cfg.npoin);
            pres_temp.resize(cfg.npoin);
            temperature_temp.resize(cfg.npoin);

            // Right-hand sides and iterative corrections.
            rhs.resize(cfg.ndim, cfg.npoin);
            rhs1.resize(cfg.npoin);
            Dunknn1.resize(cfg.ndim, cfg.npoin);
            Dtempert1.resize(cfg.npoin);
            Dudtau.resize(cfg.ndim, cfg.npoin);
            DTdtau.resize(cfg.npoin);

            // Nodal coefficients and diagnostics.
            velocity.resize(cfg.npoin);
            velocity_old.resize(cfg.npoin);
            Mdiag_real.resize(cfg.npoin);
            elcoe2.resize(cfg.npoin);
            elcoe2p.resize(cfg.npoin);
            wall_node_norm.resize(cfg.ndim, cfg.npoin);
            wall_node_list.resize(cfg.npoin);
            bc_values.resize(cfg.npoin);
            bc_list.resize(cfg.npoin);

            // Element geometry and coefficients.
            dNkdx.resize(cfg.ndim * cfg.nep * cfg.nelem);
            detJ.resize(cfg.nelem);
            delte.resize(cfg.nelem);
            elcoe_e.resize(cfg.nep * cfg.nelem);
            alen_e.resize(cfg.nelem);
            beta.resize(cfg.nelem);

            // Nodal time-step and stabilisation data.
            deltp.resize(cfg.npoin);
            deltp1.resize(cfg.npoin);
            deltp2.resize(cfg.npoin);
            beta1.resize(cfg.npoin);
            alen.resize(cfg.npoin);
            vvis.resize(cfg.npoin);
            velcp.resize(cfg.npoin);
            icflag.resize(cfg.npoin);
            unkno_unit_vec.resize(cfg.ndim, cfg.npoin);

            // Face geometry and element-face classification.
            face_norm.resize(cfg.ndim1, cfg.nboun);
            annxf.resize(cfg.ndim1, cfg.nsid, cfg.nelem);
            fedge.resize(cfg.nsid, cfg.nelem);
            Mconsist.resize(cfg.gsdim * cfg.nelem);
            M_diag.resize(cfg.npoin);

            // Pressure-system coefficients and solver work arrays.
            pdiagE.resize(cfg.nep * cfg.nelem);
            gstifE.resize(cfg.gsdim * cfg.nelem);
            pdiag.resize(cfg.npoin);
            gstif.resize(cfg.gsdim * cfg.nelem);

            apold.resize(cfg.npoin);
            press.resize(cfg.npoin);
            residual.resize(cfg.npoin);
            Mmatrix.resize(cfg.npoin);
            Zip.resize(cfg.npoin);
            pvect.resize(cfg.npoin);

            // Banded-solver work arrays.
            gsm.resize(cfg.npoin, 1);
            gasfac.resize(cfg.npoin, 1);
            gsm1.resize(cfg.npoin, 1);
            gsminv.resize(cfg.npoin, 1);
            gaspiv.resize(cfg.npoin);
            glm.resize(cfg.npoin);

            // Heat-transfer diagnostics and material data.
            nusselt_local.resize(cfg.npoin);
            nusselt_point_list.resize(cfg.npoin);

            mat_elem.resize(cfg.nelem);
            rho_e.resize(cfg.nelem);
            cp_e.resize(cfg.nelem);
            k_e.resize(cfg.nelem);
            mu_e.resize(cfg.nelem);
            rho_cp_e.resize(cfg.nelem);
            alpha_e.resize(cfg.nelem);
            Qvol_e.resize(cfg.nelem);


            // Spalart-Allmaras nodal fields.
            nu_tilde.resize(cfg.npoin);
            nu_tilde1.resize(cfg.npoin);
            nu_t.resize(cfg.npoin);
            mu_t.resize(cfg.npoin);
            wall_distance.resize(cfg.npoin);

            sa_rhs.resize(cfg.npoin);
            sa_source.resize(cfg.npoin);
            sa_production.resize(cfg.npoin);
            sa_destruction.resize(cfg.npoin);
            sa_diffusion.resize(cfg.npoin);
            sa_residual.resize(cfg.npoin);

            sa_active_node.resize(cfg.npoin);
            sa_wall_node.resize(cfg.npoin);
            sa_inlet_node.resize(cfg.npoin);

            // Spalart-Allmaras element fields.
            nu_tilde_e.resize(cfg.nelem);
            nu_t_e.resize(cfg.nelem);
            mu_t_e.resize(cfg.nelem);
            mu_eff_e.resize(cfg.nelem);
            k_eff_e.resize(cfg.nelem);

            // Default values before the input files are read.
            fedge.fill(0);
            mat_elem.fill(0);

            rho_e.fill(1.0);
            cp_e.fill(1.0);
            k_e.fill(1.0);
            mu_e.fill(0.0);
            rho_cp_e.fill(1.0);
            alpha_e.fill(1.0);
            Qvol_e.fill(0.0);


            nu_tilde.fill(0.0);
            nu_tilde1.fill(0.0);
            nu_t.fill(0.0);
            mu_t.fill(0.0);
            wall_distance.fill(1.0e300);

            sa_rhs.fill(0.0);
            sa_source.fill(0.0);
            sa_production.fill(0.0);
            sa_destruction.fill(0.0);
            sa_diffusion.fill(0.0);
            sa_residual.fill(0.0);

            sa_active_node.fill(0);
            sa_wall_node.fill(0);
            sa_inlet_node.fill(0);

            nu_tilde_e.fill(0.0);
            nu_t_e.fill(0.0);
            mu_t_e.fill(0.0);
            mu_eff_e.fill(0.0);
            k_eff_e.fill(1.0);

            deltp.fill(0.0);
            deltp1.fill(0.0);
            deltp2.fill(0.0);
            beta1.fill(0.0);
            alen.fill(0.0);
            vvis.fill(0.0);
            velcp.fill(0.0);
            icflag.fill(0);
            unkno_unit_vec.fill(0.0);

            // Reset only mesh-dependent partition data. The MPI communicator
            // context itself must survive mesh allocation because
            // Solver::setMpiContext() may be called before MeshIO::readAll().
            owned_elements.clear();
            ghost_elements.clear();
            owned_nodes.clear();
            ghost_nodes.clear();
            local_to_global_node.clear();
            local_to_global_element.clear();
            owned_pressure_nodes.clear();
            ghost_pressure_nodes.clear();

            element_owner_rank.assign(
                static_cast<Size>(cfg.nelem) + 1U,
                0);

            node_owner_rank.assign(
                static_cast<Size>(cfg.npoin) + 1U,
                0);

            global_to_local_node.assign(
                static_cast<Size>(cfg.npoin) + 1U,
                -1);

            pressure_global_to_local.assign(
                static_cast<Size>(cfg.npoin) + 1U,
                -1);
        }
    };
}
