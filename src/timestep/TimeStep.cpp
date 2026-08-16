//=============================================================================
// CBS3D++_SI
//
// Time-step calculation for the three-dimensional semi-implicit CBS solver.
//
// The local stability time step is based on the advective and diffusive
// restrictions:
//
//     dt_adv  = h / |u|
//
//     dt_diff = h^2 / (2 D)
//
//     dt_e = C_safety min(dt_adv, dt_diff)
//
// The routine also constructs the inverse momentum and thermal time diagonals
// and scales the pressure stiffness coefficients by the local element time
// step.
//=============================================================================

#include "cbs/timestep/TimeStep.hpp"

#include "cbs/turbulence/SpalartAllmaras.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace cbs
{
    namespace
    {
        // Returns the compact one-dimensional storage position of local
        // node local_node in tetrahedral element ie.
        Int element_node_index(
            const CBSStateSI& s,
            const Int ie,
            const Int local_node)
        {
            return (ie - 1) * s.cfg.nep + local_node;
        }

        // Returns the compact storage position of one off-diagonal
        // pressure coefficient in tetrahedral element ie.
        Int pressure_offdiag_index(
            const CBSStateSI& s,
            const Int ie,
            const Int ig)
        {
            return (ie - 1) * s.cfg.gsdim + ig;
        }

        // Calculates the magnitude of the nodal velocity:
        //
        //     |u_i| = sqrt(u_i^2 + v_i^2 + w_i^2)
        Real node_speed(const CBSStateSI& s, const Int ip)
        {
            return std::sqrt(
                s.unkno(1, ip) * s.unkno(1, ip) +
                s.unkno(2, ip) * s.unkno(2, ip) +
                s.unkno(3, ip) * s.unkno(3, ip));
        }

        // Returns the maximum nodal velocity magnitude in one element:
        //
        //     U_e = max_a |u_a|
        Real element_max_speed(const CBSStateSI& s, const Int ie)
        {
            Real vmax = 0.0;

            for (Int in = 1; in <= s.cfg.nep; ++in)
            {
                const Int ip = s.intma(in, ie);
                vmax = std::max(vmax, node_speed(s, ip));
            }

            return vmax;
        }

        // Calculates the magnitude of the element-average velocity:
        //
        //     u_bar_e = (1 / n_e) sum_a u_a
        //
        //     U_e = |u_bar_e|
        Real element_centroid_speed(const CBSStateSI& s, const Int ie)
        {
            Real u = 0.0;
            Real v = 0.0;
            Real w = 0.0;

            for (Int in = 1; in <= s.cfg.nep; ++in)
            {
                const Int ip = s.intma(in, ie);
                u += s.unkno(1, ip);
                v += s.unkno(2, ip);
                w += s.unkno(3, ip);
            }

            u /= static_cast<Real>(s.cfg.nep);
            v /= static_cast<Real>(s.cfg.nep);
            w /= static_cast<Real>(s.cfg.nep);

            return std::sqrt(u * u + v * v + w * w);
        }

        // Prevents division by zero in the advective stability limit.
        Real safe_speed(const Real value)
        {
            return std::max(value, 1.0e-14);
        }

        // Returns the standard geometric element length calculated
        // during mesh preprocessing.
        Real standard_element_length(const CBSStateSI& s, const Int ie)
        {
            const Real h = s.alen_e(ie);

            if (h <= 0.0 || !std::isfinite(h))
            {
                throw std::runtime_error(
                    "TimeStep - invalid alen_e at element "
                    + std::to_string(ie)
                    + ". Run Preprocess::elementSize before timestep computation.");
            }

            return h;
        }

        // Calculates the arithmetic mean of the characteristic lengths
        // stored at the four nodes of one tetrahedron:
        //
        //     h_bar_e = (1 / 4) sum_a h_a
        Real nodal_length_average(const CBSStateSI& s, const Int ie)
        {
            Real h = 0.0;

            for (Int in = 1; in <= s.cfg.nep; ++in)
            {
                const Int ip = s.intma(in, ie);

                if (s.alen(ip) <= 0.0 || !std::isfinite(s.alen(ip)))
                {
                    return standard_element_length(s, ie);
                }

                h += s.alen(ip);
            }

            return h / static_cast<Real>(s.cfg.nep);
        }

        // Selects the characteristic element length according to htype.
        //
        //     htype = 1 : standard geometric length
        //     htype = 2 : minimum of the geometric and mean nodal lengths
        //     htype = 3 : standard length with centroidal velocity
        Real timestep_length(const CBSStateSI& s, const Int ie)
        {
            if (s.cfg.htype == 1)
            {
                return standard_element_length(s, ie);
            }

            if (s.cfg.htype == 2)
            {
                // Legacy htype=2: SUPG by checking each node of the element.
                // A full directional SUPG length needs the velocity-direction
                // projection from findh_calc.  Until that full kernel is ported,
                // use the nodal characteristic length if available; otherwise
                // fall back to the conservative element height.
                return std::min(standard_element_length(s, ie), nodal_length_average(s, ie));
            }

            if (s.cfg.htype == 3)
            {
                // Legacy htype=3: SUPG using centroidal/element-average velocity.
                // The geometric h is retained here; the speed is taken from the
                // centroidal average in timestep_speed().
                return standard_element_length(s, ie);
            }

            throw std::runtime_error(
                "TimeStep - htype must be 1, 2, or 3");
        }

        // Selects the characteristic velocity used in dt_adv = h / U.
        Real timestep_speed(const CBSStateSI& s, const Int ie)
        {
            if (s.cfg.htype == 3)
            {
                return element_centroid_speed(s, ie);
            }

            return element_max_speed(s, ie);
        }

        // Selects the diffusivity controlling the element stability limit.
        //
        // In non-dimensional mode:
        //
        //     D_e = ani
        //
        // In dimensional CHT mode:
        //
        //     D_e = alpha_e                         solid
        //
        //     D_e = max(alpha_e, nu_e)             fluid
        //
        // where:
        //
        //     alpha_e = k_e / (rho_e cp_e)
        //     nu_e    = mu_e / rho_e
        Real element_diffusivity_for_timestep(
            const CBSStateSI& s,
            const Int ie)
        {
            Real diff = s.cfg.ani;

            if (s.cfg.dimensional_mode > 0 && s.cfg.material_properties_enabled > 0)
            {
                if (s.alpha_e(ie) <= 0.0 || !std::isfinite(s.alpha_e(ie)))
                {
                    throw std::runtime_error(
                        "TimeStep - invalid thermal diffusivity alpha_e at element "
                        + std::to_string(ie));
                }

                diff = s.alpha_e(ie);

                if (s.mat_elem(ie) == 0)
                {
                    if (s.rho_e(ie) <= 0.0 || s.mu_e(ie) < 0.0)
                    {
                        throw std::runtime_error(
                            "TimeStep - invalid fluid rho/mu at element "
                            + std::to_string(ie));
                    }

                    Real nu = s.mu_e(ie) / s.rho_e(ie);

                    if (s.cfg.turbulence_on > 0)
                    {
                        if (s.mu_eff_e(ie) < 0.0 || !std::isfinite(s.mu_eff_e(ie)))
                        {
                            throw std::runtime_error(
                                "TimeStep - invalid effective viscosity at element "
                                + std::to_string(ie));
                        }

                        nu = s.mu_eff_e(ie) / s.rho_e(ie);

                        if (s.cfg.temp_calc > 0 &&
                            s.cfg.turbulent_thermal_diffusivity_on > 0)
                        {
                            if (s.k_eff_e(ie) <= 0.0 || !std::isfinite(s.k_eff_e(ie)))
                            {
                                throw std::runtime_error(
                                    "TimeStep - invalid effective thermal conductivity at element "
                                    + std::to_string(ie));
                            }

                            const Real alpha_eff = s.k_eff_e(ie) / s.rho_cp_e(ie);
                            diff = std::max(diff, alpha_eff);
                        }
                    }

                    // Thermal diffusion constrains the pseudo-time only
                    // when the energy equation is actually being advanced.
                    // For an isothermal fluid calculation, the active
                    // restrictions are momentum and, when enabled, SA.
                    if (s.cfg.temp_calc > 0)
                    {
                        diff = std::max(diff, nu);
                    }
                    else
                    {
                        diff = nu;
                    }

                    // Spalart-Allmaras diffusion limit.
                    //
                    // The SA transport equation is advanced explicitly and its
                    // diffusion coefficient is
                    //
                    //     D_SA = (nu + nu_tilde) / sigma
                    //
                    // which is not nu_eff.  Since nu_t = nu_tilde * fv1 with
                    // fv1 <= 1, the transported variable exceeds the eddy
                    // viscosity it produces, and sigma = 2/3 divides rather than
                    // multiplies, so D_SA can exceed nu_eff by a factor of order
                    // 1/sigma even where the model is behaving.  Sizing the
                    // timestep on nu_eff alone therefore leaves the SA equation
                    // unprotected, which matters because SA has already been
                    // observed to destabilise on this solver while the momentum
                    // field still looked healthy.
                    if (s.cfg.turbulence_on > 0)
                    {
                        const turbulence::SpalartAllmarasConstants sa_constants;

                        const Real nu_tilde_e = s.nu_tilde_e(ie);

                        if (nu_tilde_e < 0.0 || !std::isfinite(nu_tilde_e))
                        {
                            throw std::runtime_error(
                                "TimeStep - invalid nu_tilde_e at element "
                                + std::to_string(ie));
                        }

                        const Real molecular_nu = s.mu_e(ie) / s.rho_e(ie);

                        const Real sa_diffusivity =
                            (molecular_nu + nu_tilde_e) / sa_constants.sigma;

                        diff = std::max(diff, sa_diffusivity);
                    }
                }
            }

            if (diff <= 0.0 || !std::isfinite(diff))
            {
                throw std::runtime_error(
                    "TimeStep - invalid timestep diffusivity at element "
                    + std::to_string(ie));
            }

            return diff;
        }


        // Calculates the local element stability time step:
        //
        //     dt_adv  = h / |u|                    fluid elements
        //
        //     dt_diff = h^2 / (2 D_e)
        //
        //     dt_e = csafm min(dt_adv, dt_diff)
        //
        // Solid elements have no advective restriction. When ilots = 2, the
        // result is additionally limited by deltr.
        Real local_stability_dt(
            const CBSStateSI& s,
            const Int ie,
            const Real h,
            const Real speed)
        {
            const Real diffusivity = element_diffusivity_for_timestep(s, ie);

            const Real advective_dt =
                (s.mat_elem(ie) == 0)
                ? h / safe_speed(speed)
                : std::numeric_limits<Real>::max();

            const Real diffusive_dt = h * h * (0.5 / diffusivity);

            Real dt = s.cfg.csafm * std::min(advective_dt, diffusive_dt);

            if (s.cfg.ilots == 2)
            {
                // Legacy ilots==2 keeps the local timestep structure but caps it
                // by the real-time/diffusion control deltr when provided.
                if (s.cfg.deltr > 0.0 && std::isfinite(s.cfg.deltr))
                {
                    dt = std::min(dt, s.cfg.deltr);
                }
            }

            if (dt <= 0.0 || !std::isfinite(dt))
            {
                throw std::runtime_error(
                    "TimeStep - invalid local stability timestep at element "
                    + std::to_string(ie));
            }

            return dt;
        }

        // Calculates the characteristic length at every node from the
        // minimum altitude of the connected tetrahedra.
        //
        // For a tetrahedron:
        //
        //     V_e = A_i h_i / 3
        //
        // therefore the altitude opposite local node i is:
        //
        //     h_i = 3 V_e / A_i
        //
        // Since det(J_e) = 6 V_e, the implemented three-dimensional form is:
        //
        //     h_i = 0.5 det(J_e) / A_i
        //
        // The nodal length is the minimum altitude among all elements sharing
        // that node.
        void compute_nodal_lengths_from_elements(CBSStateSI& s)
        {
            // Use the minimum tetrahedral altitude connected to each node.
            // In three dimensions:
            //
            //     h_i = 3 V_e / A_i = 0.5 det(J_e) / A_i
            const Real lfact = (s.cfg.ndim == 2) ? 1.0 : 0.5;
            const Real big = std::numeric_limits<Real>::max();

            s.alen.fill(big);

            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                if (s.detJ(ie) <= 0.0 || !std::isfinite(s.detJ(ie)))
                {
                    throw std::runtime_error(
                        "TimeStep - invalid detJ in nodal length; run Preprocess first");
                }

                const Real volum = lfact * s.detJ(ie);

                for (Int in = 1; in <= s.cfg.nep; ++in)
                {
                    const Int ip = s.intma(in, ie);
                    const Real face_area = s.annxf(s.cfg.ndim1, in, ie);

                    if (face_area <= 0.0 || !std::isfinite(face_area))
                    {
                        throw std::runtime_error(
                            "TimeStep - invalid face area in nodal length; run Preprocess::getNormals first");
                    }

                    const Real hite = volum / face_area;
                    s.alen(ip) = std::min(s.alen(ip), hite);
                }
            }
        }

        // Stores the nodal velocity magnitude and unit direction:
        //
        //     velcp_i = |u_i|
        //
        //     e_u,i = u_i / |u_i|
        //
        // The unit vector is set to zero at stationary nodes.
        void compute_velocity_helpers(CBSStateSI& s)
        {
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                const Real speed = node_speed(s, ip);
                s.velcp(ip) = speed;

                if (speed > 1.0e-14)
                {
                    s.unkno_unit_vec(1, ip) = s.unkno(1, ip) / speed;
                    s.unkno_unit_vec(2, ip) = s.unkno(2, ip) / speed;
                    s.unkno_unit_vec(3, ip) = s.unkno(3, ip) / speed;
                }
                else
                {
                    s.unkno_unit_vec(1, ip) = 0.0;
                    s.unkno_unit_vec(2, ip) = 0.0;
                    s.unkno_unit_vec(3, ip) = 0.0;
                }
            }
        }

        // Identifies whether each node is connected to at least one fluid
        // element and/or at least one solid element.
        //
        // Interface nodes therefore have both flags equal to one.
        void build_node_material_touch_flags(
            const CBSStateSI& s,
            std::vector<char>& touches_fluid,
            std::vector<char>& touches_solid)
        {
            touches_fluid.assign(static_cast<std::size_t>(s.cfg.npoin + 1), 0);
            touches_solid.assign(static_cast<std::size_t>(s.cfg.npoin + 1), 0);

            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                const bool fluid_element = (s.mat_elem(ie) == 0);

                for (Int in = 1; in <= s.cfg.nep; ++in)
                {
                    const Int ip = s.intma(in, ie);

                    if (fluid_element)
                    {
                        touches_fluid[static_cast<std::size_t>(ip)] = 1;
                    }
                    else
                    {
                        touches_solid[static_cast<std::size_t>(ip)] = 1;
                    }
                }
            }
        }

        // Checks the problem sizes and time-step controls required by all
        // public routines in this module.
        void validate_common_inputs(const CBSStateSI& s)
        {
            if (s.cfg.nelem < 1 || s.cfg.npoin < 1)
            {
                throw std::runtime_error(
                    "TimeStep - problem sizes are not initialised");
            }

            if (s.cfg.csafm <= 0.0 || !std::isfinite(s.cfg.csafm))
            {
                throw std::runtime_error(
                    "TimeStep - csafm must be positive");
            }

            if (s.cfg.dtfix <= 0.0 || !std::isfinite(s.cfg.dtfix))
            {
                throw std::runtime_error(
                    "TimeStep - dtfix must be positive");
            }

            if (s.cfg.htype < 1 || s.cfg.htype > 3)
            {
                throw std::runtime_error(
                    "TimeStep - htype must be 1, 2, or 3");
            }
        }
    }

    //=========================================================================
    // Calculates the element and nodal time steps for one CBS iteration.
    //
    // Fixed-step mode assigns:
    //
    //     delte(e) = deltp(i) = deltp2(i) = dtfix
    //
    // Otherwise, each element receives:
    //
    //     dt_e = csafm min(h_e / U_e, h_e^2 / (2 D_e))
    //
    // and the global time step is:
    //
    //     dtreal = min_e(dt_e)
    //
    // Local values are retained only for ilots = 1 or ilots = 2. Transient
    // calculations and the pressure-solver modes requiring a common operator
    // use the global minimum.
    //
    // Inputs:
    //     velocity field, element lengths, material diffusivities and controls
    //
    // Outputs:
    //     delte, deltp, deltp1, deltp2, beta1 and dtreal
    //=========================================================================
    void TimeStep::computeTimeStep(CBSStateSI& s, const Int iitime)
    {
        validate_common_inputs(s);

        if (s.cfg.cbs_scheme != 1 && s.cfg.cbs_scheme != 0)
        {
            throw std::runtime_error(
                "TimeStep::computeTimeStep - cbs_scheme must be 1 semi-implicit or 0 explicit");
        }

        if (s.cfg.cbs_scheme == 0)
        {
            // The timestep infrastructure accepts the explicit-mode flag, but
            // the explicit CBS3D step kernels must be ported before this path
            // can be used safely by the solver driver.
            throw std::runtime_error(
                "Fully explicit CBS3D mode was requested, but explicit CBS3D step kernels are not yet ported.");
        }

        if (s.cfg.ilots <= -1 || s.cfg.dtfixed > 0 ||
            (s.cfg.dtfix_end > 0 && iitime <= s.cfg.dtfix_end))
        {
            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                s.delte(ie) = s.cfg.dtfix;
            }

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                s.deltp(ip) = s.cfg.dtfix;
                s.deltp2(ip) = s.cfg.dtfix;
                s.beta1(ip) = 1.0;
            }

            s.cfg.dtreal = s.cfg.dtfix;
            return;
        }

        compute_velocity_helpers(s);

        // Nodal characteristic lengths are needed by the htype=2 path and by
        // legacy-style diagnostics even when element timesteps are ultimately
        // forced to a global minimum.
        compute_nodal_lengths_from_elements(s);

        computeNodalLocalTimeStep(s);

        Real dt_min = std::numeric_limits<Real>::max();

        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            const Real h = timestep_length(s, ie);
            const Real speed = timestep_speed(s, ie);

            s.delte(ie) = local_stability_dt(s, ie, h, speed);
            dt_min = std::min(dt_min, s.delte(ie));
        }

        if (dt_min <= 0.0 || !std::isfinite(dt_min))
        {
            throw std::runtime_error(
                "TimeStep::computeTimeStep - invalid global minimum timestep");
        }

        // Correct legacy ilots semantics:
        //
        //   ilots == 1 : retain local element/nodal timestep
        //   ilots == 2 : retain local timestep with extra cap already applied
        //   otherwise  : force global minimum timestep
        //
        // Transient real-time runs also use one global dt.
        if (s.cfg.transient_on > 0 ||
            s.cfg.solver_opt > 1 ||
            s.cfg.dtfixed == -1 ||
            (s.cfg.ilots != 1 && s.cfg.ilots != 2))
        {
            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                s.delte(ie) = dt_min;
            }

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                s.deltp(ip) = dt_min;
                s.deltp2(ip) = dt_min;
            }
        }

        if (s.cfg.rem_deltp > 0)
        {
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                if (s.deltp1(ip) > 0.0)
                {
                    s.deltp(ip) = std::min(s.deltp(ip), s.deltp1(ip));
                    s.deltp2(ip) = std::min(s.deltp2(ip), s.deltp1(ip));
                }

                s.deltp1(ip) = s.deltp(ip);
            }
        }
        else
        {
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                s.deltp1(ip) = s.deltp(ip);
            }
        }

        s.cfg.dtreal = dt_min;
    }

    //=========================================================================
    // Calculates a local stability time step at every node.
    //
    //     dt_adv,i  = h_i / |u_i|
    //
    //     dt_diff,i = h_i^2 / (2 ani)
    //
    //     dt_i = csafm min(dt_adv,i, dt_diff,i)
    //
    // When beta_opt is enabled, the nodal blending parameter is:
    //
    //     beta_i = dt_diff,i / (dt_adv,i + dt_diff,i)
    //
    // limited to the interval [0,1].
    //=========================================================================
    void TimeStep::computeNodalLocalTimeStep(CBSStateSI& s)
    {
        validate_common_inputs(s);

        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            const Real h =
                (s.alen(ip) > 0.0 && std::isfinite(s.alen(ip)))
                    ? s.alen(ip)
                    : 0.0;

            const Real speed = safe_speed(node_speed(s, ip));

            if (h <= 0.0)
            {
                s.deltp(ip) = s.cfg.dtfix;
                s.deltp2(ip) = s.cfg.dtfix;
                s.beta1(ip) = 1.0;
                s.vvis(ip) = 0.0;
                continue;
            }

            const Real advective_dt = h / speed;
            const Real diffusive_dt = h * h * (0.5 / s.cfg.ani);

            s.vvis(ip) = diffusive_dt;

            Real dt = s.cfg.csafm * std::min(advective_dt, diffusive_dt);

            if (s.cfg.ilots == 2 &&
                s.cfg.deltr > 0.0 &&
                std::isfinite(s.cfg.deltr))
            {
                dt = std::min(dt, s.cfg.deltr);
            }

            if (dt <= 0.0 || !std::isfinite(dt))
            {
                throw std::runtime_error(
                    "TimeStep::computeNodalLocalTimeStep - invalid nodal timestep at node "
                    + std::to_string(ip));
            }

            s.deltp(ip) = dt;
            s.deltp2(ip) = dt;

            if (s.cfg.beta_opt > 0)
            {
                const Real denom = advective_dt + diffusive_dt + 1.0e-30;
                s.beta1(ip) = std::clamp(diffusive_dt / denom, 0.0, 1.0);
            }
            else
            {
                s.beta1(ip) = 1.0;
            }
        }
    }

    //=========================================================================
    // Forces one global real-time step after calculating the local stability
    // limits:
    //
    //     dt_global = min_e(dt_e)
    //
    // The same value is assigned to every element and node.
    //=========================================================================
    void TimeStep::computeGlobalRealTimeStep(CBSStateSI& s, const Int iitime)
    {
        computeTimeStep(s, iitime);

        Real dt_min = std::numeric_limits<Real>::max();

        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            dt_min = std::min(dt_min, s.delte(ie));
        }

        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            s.delte(ie) = dt_min;
        }

        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            s.deltp(ip) = dt_min;
            s.deltp2(ip) = dt_min;
        }

        s.cfg.dtreal = dt_min;
    }

    //=========================================================================
    // Applies the optional pressure-based time-step restriction before CBS
    // Step 2.
    //
    // At every element node, the candidate value is:
    //
    //     dt_p = csafm2 h (|p_old| + epsilon1) / |u|
    //
    // and the element time step is reduced to the minimum admissible value.
    //
    // This routine only decreases delte; it never increases the stability
    // time step calculated previously.
    //=========================================================================
    void TimeStep::applyStep2PressureTimeStepCorrection(CBSStateSI& s)
    {
        if (s.cfg.step2_check < 1)
        {
            return;
        }

        if (s.cfg.csafm2 <= 0.0 || !std::isfinite(s.cfg.csafm2))
        {
            throw std::runtime_error(
                "TimeStep::applyStep2PressureTimeStepCorrection - csafm2 must be positive");
        }

        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            Real dt_face_min = s.delte(ie);

            for (Int in = 1; in <= s.cfg.nep; ++in)
            {
                const Int ip = s.intma(in, ie);

                const Real pscale =
                    std::abs(s.pres1(ip)) + s.cfg.epsilon1;

                const Real speed = node_speed(s, ip);

                if (speed > 1.0e-14 && pscale > 0.0)
                {
                    const Real h =
                        (s.alen(ip) > 0.0)
                            ? s.alen(ip)
                            : standard_element_length(s, ie);

                    const Real candidate =
                        s.cfg.csafm2 * h * pscale / speed;

                    if (candidate > 0.0 && std::isfinite(candidate))
                    {
                        dt_face_min = std::min(dt_face_min, candidate);
                    }
                }
            }

            s.delte(ie) = dt_face_min;
        }
    }

    //=========================================================================
    // Forms the time-dependent left-hand-side coefficients used by the
    // momentum, energy and pressure equations.
    //
    // Momentum diagonal, fluid domain only:
    //
    //     M_u,i / dt = sum_e m_i^(e) / dt_e
    //
    //     elcoe2(i) = [M_u,i / dt]^-1
    //
    // Thermal capacitance diagonal, fluid and solid domains:
    //
    //     C_T,i / dt = sum_e rho_e cp_e m_i^(e) / dt_e
    //
    //     elcoe2p(i) = [C_T,i / dt]^-1
    //
    // Pressure operator, fluid elements only:
    //
    //     A_p^(e) = dt_e H^(e)
    //
    //     gstif^(e) = dt_e gstifE^(e)
    //
    //     pdiag_i = sum_e dt_e pdiagE_i^(e)
    //
    // Prescribed-pressure nodes are converted to unit diagonal equations.
    // Solid-only nodes remain inactive in momentum and pressure.
    //=========================================================================
    void TimeStep::updateLhsDiagonal(CBSStateSI& s)
    {
        s.elcoe2.fill(0.0);
        s.elcoe2p.fill(0.0);
        s.pdiag.fill(0.0);

        std::vector<char> touches_fluid;
        std::vector<char> touches_solid;
        build_node_material_touch_flags(s, touches_fluid, touches_solid);

        // Momentum and thermal mass/capacitance are different fields in CHT.
        //
        //   elcoe2  : inverse momentum mass/time diagonal, fluid domain only.
        //   elcoe2p : inverse thermal capacitance/time diagonal, fluid + solid.
        //
        // Solid elements have no velocity equation and must not contribute to
        // elcoe2.  They do contribute to elcoe2p because Step 4 solves thermal
        // diffusion in both domains.
        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            const Real dt = s.delte(ie);

            if (dt <= 0.0 || !std::isfinite(dt))
            {
                throw std::runtime_error(
                    "TimeStep::updateLhsDiagonal - non-positive local timestep at element "
                    + std::to_string(ie));
            }

            const bool fluid_element = (s.mat_elem(ie) == 0);

            for (Int in = 1; in <= s.cfg.nep; ++in)
            {
                const Int ip = s.intma(in, ie);
                const Real mass = s.elcoe_e(element_node_index(s, ie, in));

                if (mass <= 0.0 || !std::isfinite(mass))
                {
                    throw std::runtime_error(
                        "TimeStep::updateLhsDiagonal - invalid element-node mass");
                }

                if (fluid_element)
                {
                    s.elcoe2(ip) += mass / dt;
                }

                s.elcoe2p(ip) += s.rho_cp_e(ie) * mass / dt;
            }
        }

        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            const bool fluid =
                touches_fluid[static_cast<std::size_t>(ip)] != 0;
            const bool solid =
                touches_solid[static_cast<std::size_t>(ip)] != 0;

            if (!fluid && !solid)
            {
                throw std::runtime_error(
                    "TimeStep::updateLhsDiagonal - orphan node detected at node "
                    + std::to_string(ip));
            }

            if (fluid)
            {
                if (s.elcoe2(ip) <= 0.0 || !std::isfinite(s.elcoe2(ip)))
                {
                    throw std::runtime_error(
                        "TimeStep::updateLhsDiagonal - non-positive fluid momentum mass/time diagonal at node "
                        + std::to_string(ip));
                }

                s.elcoe2(ip) = 1.0 / s.elcoe2(ip);
            }
            else
            {
                // Solid-only node: no velocity/momentum equation.
                s.elcoe2(ip) = 0.0;
            }

            if (s.elcoe2p(ip) <= 0.0 || !std::isfinite(s.elcoe2p(ip)))
            {
                throw std::runtime_error(
                    "TimeStep::updateLhsDiagonal - non-positive capacitance/time diagonal at node "
                    + std::to_string(ip));
            }

            s.elcoe2p(ip) = 1.0 / s.elcoe2p(ip);
        }

        // Pressure stiffness and pressure Jacobi diagonal are fluid-domain
        // operators.  Solid elements are explicitly skipped rather than relying
        // on their assembled pressure coefficients being zero.
        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            const Real dt = s.delte(ie);

            if (s.mat_elem(ie) != 0)
            {
                for (Int ig = 1; ig <= s.cfg.gsdim; ++ig)
                {
                    const Int idx = pressure_offdiag_index(s, ie, ig);
                    s.gstif(idx) = 0.0;
                }

                continue;
            }

            for (Int ig = 1; ig <= s.cfg.gsdim; ++ig)
            {
                const Int idx = pressure_offdiag_index(s, ie, ig);
                s.gstif(idx) = s.gstifE(idx) * dt;
            }

            for (Int in = 1; in <= s.cfg.nep; ++in)
            {
                const Int idx = element_node_index(s, ie, in);
                const Int ip = s.intma(in, ie);

                s.pdiag(ip) += s.pdiagE(idx) * dt;
            }
        }

        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            const bool fluid =
                touches_fluid[static_cast<std::size_t>(ip)] != 0;
            const bool solid =
                touches_solid[static_cast<std::size_t>(ip)] != 0;

            if (fluid)
            {
                if (std::abs(s.pdiag(ip)) <= 1.0e-14 || !std::isfinite(s.pdiag(ip)))
                {
                    throw std::runtime_error(
                        "TimeStep::updateLhsDiagonal - zero pressure diagonal on fluid-connected node "
                        + std::to_string(ip));
                }
            }
            else if (solid)
            {
                // Solid-only node: inactive pressure location, not a pressure BC.
                // Do not add it to bc_list; CG/matvec logic must exclude it from
                // pressure dot products and pressure updates.
                s.pdiag(ip) = 0.0;
                s.pres(ip) = 0.0;
                s.rhs1(ip) = 0.0;
            }
        }

        // Real pressure constraints only: outlet nodes or one reference node for
        // closed fluid cavities.  Solid-only inactive nodes must not be inserted
        // into this list.
        for (Int i = 1; i <= s.cfg.bc_fixed; ++i)
        {
            const Int ip = s.bc_list(i);

            if (ip < 1 || ip > s.cfg.npoin)
            {
                throw std::runtime_error(
                    "TimeStep::updateLhsDiagonal - fixed pressure node out of range");
            }

            if (touches_fluid[static_cast<std::size_t>(ip)] == 0)
            {
                throw std::runtime_error(
                    "TimeStep::updateLhsDiagonal - fixed pressure node is not connected to the fluid pressure domain, node "
                    + std::to_string(ip));
            }

            s.pdiag(ip) = 1.0;
            s.pres(ip) = s.bc_values(i);
            s.rhs1(ip) = s.bc_values(i);
        }
    }

    //=========================================================================
    // Reserved hook for additional real-time contributions.
    //
    // The current implementation intentionally performs no operation.
    //=========================================================================
    void TimeStep::updateRealTimeTerms(CBSStateSI& s)
    {
        (void)s;
    }
}
