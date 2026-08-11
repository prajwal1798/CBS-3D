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
        std::string case_name;

        Int ndim = 3;
        Int nep = 4;
        Int nsid = 4;
        Int nsidp = 3;
        Int ndim1 = 4;
        Int gdim = 13;
        Int gsdim = 6;
        Int bsid = 6;
        Int nsidpl = 4;
        Int nsidpe = 5;

        Int npoin = 0;
        Int nelem = 0;
        Int nboun = 0;
        Int nbw = 0;
        Int nflag = 0;

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

        Int write_output = 0;
        Int paraview_output = 0;
        Int tecplot_output = 0;
        Int nusselt_calc = 0;
        Int nusselt_flag = 0;

        Int vel_check = 0;
        Int temp_check = 0;

        Int istart = 1;
        Int iitime_start = 1;
        Int npoin_wall = 0;
        Int itrail = 0;
        Int bc_fixed = 0;
        Int first_access = 1;
        Int iiter = 0;
        Int iiter_total = 0;

        std::array<Real, 5> cinf =
        {
            0.0, 0.0, 0.0, 0.0, 0.0
        };

        std::array<Real, 2> theta =
        {
            0.0, 1.0
        };

        Real dtfix = 0.0;
        Real dtreal = 0.0;
        Real csafm = 0.0;

        Int cbs_scheme = 1;
        Int ilots = 0;
        Int htype = 1;
        Int step2_check = 0;
        Int rem_deltp = 0;
        Int beta_opt = 1;
        Int dtfix_end = 0;

        Int residual_log_enabled = 1;
        Int residual_log_every = 1;
        Int console_log_every = 10;
        Int live_residual_plot = 0;
        Int vtu_output_enabled = 1;
        Int vtu_output_every_iterations = 100;
        Real vtu_output_every_sim_time = 0.0;
        Int write_boundary_debug_arrays = 1;
        Int steady_min_iterations = 5;

        Int art_diff = 0;

        Int cg_preconditioner = 1;
        Int cg_conv_test = 3;
        Int cg_max_iter = 100000;
        Real csafm2 = 2.0;
        Real epsilon1 = 1.0e-6;
        Real deltr = 1.0e20;

        Real re = 0.0;
        Real pr = 0.0;
        Real ra = 0.0;
        Real ri = 0.0;

        Real rtime = 0.0;
        Real rtime_output = 0.0;
        Real frac = 0.0;
        Real end_rtime = 0.0;
        Real ani = 0.0;

        // --------------------------------------------------------------------
        // Thermal and source controls
        // --------------------------------------------------------------------
        Real alpha_sf = 1.0;
        Real k_ratio = 1.0;

        // Uniform volumetric heat generation applied to every solid element.
        // Dimensional-mode unit: W/m^3. Zero disables the .par-level source
        // and leaves per-material .matprop Qvol available.
        Real source_solid = 0.0;

        // Positive value denotes heat entering BC 532.
        Real heat_flux_bc = 0.0;

        Int turbulence_on = 0;
        Int turbulence_model = 0;
        Int turbulent_thermal_diffusivity_on = 0;
        Real sa_inlet_ratio = 3.0;
        Real sa_prandtl_t = 0.90;
        Real sa_min_wall_distance = 1.0e-14;
        Real sa_min_stilde = 1.0e-14;
        Real sa_nu_tilde_floor = 0.0;
        Int sa_use_stilde_limiter = 1;
        Int sa_implicit_destruction = 1;

        Real nusselt_Tinf = 25.0;
        Real nusselt_Tref = 100.0;
        Real nusselt_diameter = 1.0;

        Real pboundary = 0.0;
        Real outlet_pressure_gauge = 0.0;

        Int dimensional_mode = 0;
        Int material_properties_enabled = 0;
        Int mass_flow_inlet_enabled = 0;

        Real p_ref = 15.5e6;
        Real model_depth = 1.0;

        Real inlet_mass_flow_rate = 0.0;
        Real inlet_density = 1.0;
        Real inlet_u = 0.0;
        Real inlet_v = 0.0;
        Real inlet_w = 0.0;
        Real inlet_u_from_massflow = 0.0;
        Real inlet_temperature = 568.15;

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
        Int bc_cht_interface = 901;
        Int bc_heatflux_marker = 902;

        Real write_time_output = 0.0;
        Real time_output_interval = 0.0;
        Real relToler = 0.0;
        Real absToler = 0.0;
        Real l2norm_vel_tolerance = 0.0;
        Real l2norm_pres_tolerance = 0.0;
        Real l2norm_temp_tolerance = 0.0;

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

        Real mass_factor = 1.0 / 24.0;
    };
}
