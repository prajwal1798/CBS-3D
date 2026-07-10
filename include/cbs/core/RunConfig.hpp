#pragma once

//=============================================================================
// CBS3D++_SI
//
// Numerical, physical and output controls for one CBS3D calculation.
//
// The structure stores values read from the case parameter file together with
// fixed topology constants for a four-node linear tetrahedral element.
//
// Important conventions:
//
//     ndim   = 3     spatial dimensions
//     nep    = 4     nodes per tetrahedron
//     nsid   = 4     triangular faces per tetrahedron
//     nsidp  = 3     nodes per triangular face
//     gsdim  = 6     off-diagonal entries of a symmetric 4 x 4 matrix
//
// Free-stream ordering:
//
//     cinf[0] = Ux
//     cinf[1] = Uy
//     cinf[2] = Uz
//     cinf[3] = pressure
//     cinf[4] = temperature
//
// Index zero is retained in fcon and fdif so that their active entries match
// the one-based notation used in the numerical derivation.
//=============================================================================

#include "cbs/core/Types.hpp"

#include <array>
#include <string>

namespace cbs
{
    struct RunConfig
    {
        // Base name of the current problem.
        std::string case_name;

        // --------------------------------------------------------------------
        // Linear tetrahedral topology
        // --------------------------------------------------------------------
        Int ndim = 3;
        Int nep = 4;
        Int nsid = 4;
        Int nsidp = 3;

        // annxf stores three normal components and one face area.
        Int ndim1 = 4;

        // Per-element gradient storage:
        //
        //     ndim * nep + 1 = 3 * 4 + 1 = 13
        Int gdim = 13;

        // Six off-diagonal pairs of a symmetric four-node element matrix.
        Int gsdim = 6;

        // Boundary-face record:
        //
        //     1..3  face nodes
        //     4     local face number
        //     5     parent element
        //     6     boundary identifier
        Int bsid = 6;
        Int nsidpl = 4;
        Int nsidpe = 5;

        // --------------------------------------------------------------------
        // Mesh sizes
        // --------------------------------------------------------------------
        Int npoin = 0;
        Int nelem = 0;
        Int nboun = 0;
        Int nbw = 0;
        Int nflag = 0;

        // --------------------------------------------------------------------
        // Main solver controls
        // --------------------------------------------------------------------
        Int solver_opt = 1;
        Int restart_opt = 0;
        Int temp_calc = 0;
        Int ntime = 0;
        Int transient_on = 0;
        Int dtfixed = 0;
        Int iwrite = 1;
        Int convection_type = 0;
        Int pnode = 1;
        Int runtime_mod = 0;

        // Legacy output and diagnostic controls retained by the input format.
        Int write_output = 0;
        Int paraview_output = 0;
        Int tecplot_output = 0;
        Int nusselt_calc = 0;
        Int nusselt_flag = 0;

        // Enabled steady-state checks.
        Int vel_check = 0;
        Int temp_check = 0;

        // Restart, wall and iteration bookkeeping.
        Int istart = 1;
        Int iitime_start = 1;
        Int npoin_wall = 0;
        Int itrail = 0;
        Int bc_fixed = 0;
        Int first_access = 1;
        Int iiter = 0;
        Int iiter_total = 0;

        // Free-stream values: Ux, Uy, Uz, pressure and temperature.
        std::array<Real, 5> cinf =
        {
            0.0, 0.0, 0.0, 0.0, 0.0
        };

        // CBS time-weighting parameters.
        std::array<Real, 2> theta =
        {
            0.0, 1.0
        };

        Real dtfix = 0.0;
        Real dtreal = 0.0;
        Real csafm = 0.0;

        // --------------------------------------------------------------------
        // CBS formulation and time-step controls
        // --------------------------------------------------------------------
        // cbs_scheme:
        //
        //     1  semi-implicit pressure-correction CBS
        //     0  explicit CBS path
        //
        // ilots:
        //
        //     <= -1  fixed global dtfix
        //        1   local nodal and element pseudo-time step
        //        2   local pseudo-time with additional cap logic
        //     other  calculate local values, then use the global minimum
        //
        // htype:
        //
        //     1  standard geometric element length
        //     2  nodal/SUPG-style characteristic length
        //     3  geometric length with centroidal velocity
        Int cbs_scheme = 1;
        Int ilots = 0;
        Int htype = 1;
        Int step2_check = 0;
        Int rem_deltp = 0;
        Int beta_opt = 1;
        Int dtfix_end = 0;

        // --------------------------------------------------------------------
        // Residual, console and VTU output
        // --------------------------------------------------------------------
        Int residual_log_enabled = 1;
        Int residual_log_every = 1;
        Int console_log_every = 10;
        Int live_residual_plot = 0;

        Int vtu_output_enabled = 1;
        Int vtu_output_every_iterations = 100;
        Real vtu_output_every_sim_time = 0.0;

        Int write_boundary_debug_arrays = 1;
        Int steady_min_iterations = 5;

        // Artificial-diffusion switch retained for the numerical controls.
        Int art_diff = 0;

        // --------------------------------------------------------------------
        // Pressure Conjugate Gradient controls
        // --------------------------------------------------------------------
        // cg_preconditioner:
        //
        //     0  no preconditioner
        //     1  Jacobi diagonal preconditioner
        //
        // cg_conv_test:
        //
        //     1  maximum absolute residual
        //     2  relative L2 residual
        //     3  either absolute or relative criterion
        Int cg_preconditioner = 1;
        Int cg_conv_test = 3;
        Int cg_max_iter = 100000;

        Real csafm2 = 2.0;
        Real epsilon1 = 1.0e-6;
        Real deltr = 1.0e20;

        // --------------------------------------------------------------------
        // Dimensionless physical parameters
        // --------------------------------------------------------------------
        Real re = 0.0;
        Real pr = 0.0;
        Real ra = 0.0;
        Real ri = 0.0;

        // --------------------------------------------------------------------
        // Physical and pseudo-time state
        // --------------------------------------------------------------------
        Real rtime = 0.0;
        Real rtime_output = 0.0;
        Real frac = 0.0;
        Real end_rtime = 0.0;

        // Non-dimensional momentum diffusivity.
        Real ani = 0.0;

        // --------------------------------------------------------------------
        // Thermal and source controls
        // --------------------------------------------------------------------
        Real alpha_sf = 1.0;
        Real k_ratio = 1.0;
        Real source_solid = 1.0;
        Real heat_flux_bc = 0.0;


        // --------------------------------------------------------------------
        // Spalart-Allmaras turbulence controls
        // --------------------------------------------------------------------
        // turbulence_on:
        //     0  laminar solver path; all turbulence arrays remain passive
        //     1  solve and apply the selected one-equation turbulence model
        //
        // turbulence_model:
        //     0  standard Spalart-Allmaras model
        //     1  SA-neg compatible branch, reserved for the robust variant
        //
        // turbulent_thermal_diffusivity_on:
        //     0  turbulent viscosity affects momentum only
        //     1  fluid thermal conductivity receives rho*cp*nu_t/Pr_t
        Int turbulence_on = 0;
        Int turbulence_model = 0;
        Int turbulent_thermal_diffusivity_on = 0;

        // Fully turbulent inlet/farfield working-variable ratio:
        //     nu_tilde_inlet = sa_inlet_ratio * nu_molecular
        Real sa_inlet_ratio = 3.0;

        // Turbulent Prandtl number used only for fluid heat-transfer coupling.
        Real sa_prandtl_t = 0.90;

        // Numerical protections for the first non-negative SA implementation.
        Real sa_min_wall_distance = 1.0e-14;
        Real sa_min_stilde = 1.0e-14;
        Real sa_nu_tilde_floor = 0.0;

        Int sa_use_stilde_limiter = 1;
        Int sa_implicit_destruction = 1;

        Real nusselt_Tinf = 25.0;
        Real nusselt_Tref = 100.0;
        Real nusselt_diameter = 1.0;

        // --------------------------------------------------------------------
        // Pressure references and dimensional material mode
        // --------------------------------------------------------------------
        Real pboundary = 0.0;
        Real outlet_pressure_gauge = 0.0;

        Int dimensional_mode = 0;
        Int material_properties_enabled = 0;
        Int mass_flow_inlet_enabled = 0;

        Real p_ref = 15.5e6;
        Real model_depth = 1.0;

        // --------------------------------------------------------------------
        // Inlet conditions
        // --------------------------------------------------------------------
        // Mass-flow velocity magnitude:
        //
        //     U_in = mass_flow_rate / (rho_in A_in)
        Real inlet_mass_flow_rate = 0.0;
        Real inlet_density = 1.0;
        Real inlet_u = 0.0;
        Real inlet_v = 0.0;
        Real inlet_w = 0.0;
        Real inlet_u_from_massflow = 0.0;
        Real inlet_temperature = 568.15;

        // --------------------------------------------------------------------
        // Solver boundary-condition identifiers
        // --------------------------------------------------------------------
        Int bc_adiabatic_prescribed_velocity = 500;
        Int bc_temperature_one_noslip = 501;
        Int bc_temperature_zero_noslip = 502;
        Int bc_temperature_zero_prescribed_velocity = 503;
        Int bc_pressure = 504;
        Int bc_symmetry_no_flux = 506;
        Int bc_bfs_parabolic_inlet = 507;
        Int bc_parabolic_inlet = 508;

        Int bc_velocity_temperature_inlet = 510;
        Int bc_massflow_temperature_inlet = 511;
        Int bc_pressure_outlet = 520;
        Int bc_noslip_adiabatic_wall = 530;
        Int bc_noslip_heatflux_wall = 532;

        // Conformal CHT interface.
        Int bc_cht_interface = 901;

        // Legacy heat-flux marker. This is not the CHT interface identifier.
        Int bc_heatflux_marker = 902;

        // --------------------------------------------------------------------
        // Output-time and convergence tolerances
        // --------------------------------------------------------------------
        Real write_time_output = 0.0;
        Real time_output_interval = 0.0;

        Real relToler = 0.0;
        Real absToler = 0.0;

        Real l2norm_vel_tolerance = 0.0;
        Real l2norm_pres_tolerance = 0.0;
        Real l2norm_temp_tolerance = 0.0;

        // --------------------------------------------------------------------
        // P1 tetrahedral integration constants
        // --------------------------------------------------------------------
        // For det(J) = 6V:
        //
        //     fcon[1] = 1/24  gives V/4
        //     fcon[2] = 1/12  triangular consistent face factor
        //     fcon[3] = 1/6   gives V
        //     fcon[4] = 1/3   triangular face factor
        std::array<Real, 5> fcon =
        {
            0.0,
            1.0 / 24.0,
            1.0 / 12.0,
            1.0 / 6.0,
            1.0 / 3.0
        };

        std::array<Real, 3> fdif =
        {
            0.0,
            1.0 / 6.0,
            1.0 / 3.0
        };

        // Lumped nodal mass for one P1 tetrahedron:
        //
        //     m_a^(e) = V_e/4 = det(J_e)/24
        Real mass_factor = 1.0 / 24.0;
    };
}
