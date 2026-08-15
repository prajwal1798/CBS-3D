//=============================================================================
// CBS3D++_SI
//
// Convergence measures for the CBS iteration loop.
//
// For a nodal variable phi:
//
//     delta_phi_i = phi_i^(n+1) - phi_i^n
//
// Relative residual:
//
//     R_phi =
//         sqrt[ sum_i delta_phi_i^2
//               /
//               (sum_i (phi_i^(n+1))^2 + epsilon) ]
//
// Current-field L2 norm:
//
//     ||phi^(n+1)||_2 = sqrt[ sum_i (phi_i^(n+1))^2 ]
//
// Velocity and temperature also report a residual reconstructed on the
// assembled RHS scale:
//
//     r_u,i = delta_u_i / elcoe2(i)
//
//     r_T,i = delta_T_i / elcoe2p(i)
//
// The convergence norms are evaluated only over the physically active nodes of
// each field.
//=============================================================================

#include "cbs/solver/Convergence.hpp"

#ifdef CBS3D_USE_MPI
#include <mpi.h>
#endif

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace cbs
{
    namespace
    {
        // Small positive number used to prevent division by zero when the
        // current-field L2 norm is zero.
        constexpr Real eps_norm = 1.0e-30;


        // Checks the dimensions required by the three-dimensional convergence
        // calculations.
        void validate_convergence_state(const CBSStateSI& s)
        {
            if (s.cfg.ndim != 3)
            {
                throw std::runtime_error(
                    "Convergence - CBS3D convergence requires ndim=3");
            }

            if (s.cfg.npoin < 1)
            {
                throw std::runtime_error(
                    "Convergence - npoin is not initialised");
            }
        }


        // Converts an inverse diagonal back to its original diagonal weight.
        //
        // For a nodal update:
        //
        //     delta_phi_i = D_i^(-1) r_i
        //
        // the corresponding RHS-scale residual is:
        //
        //     r_i = delta_phi_i / D_i^(-1)
        //
        // Invalid or effectively zero inverse weights are replaced by one to
        // keep the diagnostic calculation finite.
        Real safe_inverse_weight(const Real inv_mass)
        {
            if (std::abs(inv_mass) <= 1.0e-300 || !std::isfinite(inv_mass))
            {
                return 1.0;
            }

            return 1.0 / inv_mass;
        }


        // Material-connectivity masks used to define the physical domain of
        // each convergence norm.
        struct MaterialNodeMasks
        {
            std::vector<char> touches_fluid;
            std::vector<char> touches_solid;
            std::vector<char> velocity_active;
            std::vector<char> pressure_active;
            std::vector<char> thermal_active;
        };


        // Builds the material-domain masks from element connectivity.
        //
        // A node may be:
        //
        //     fluid-only       touches fluid but not solid
        //     solid-only       touches solid but not fluid
        //     interface        touches both fluid and solid
        //
        // The active convergence spaces are:
        //
        //     velocity_active = fluid AND NOT solid
        //
        //     pressure_active = fluid
        //
        //     thermal_active  = fluid OR solid
        MaterialNodeMasks build_material_node_masks(const CBSStateSI& s)
        {
            MaterialNodeMasks masks;
            masks.touches_fluid.assign(static_cast<std::size_t>(s.cfg.npoin + 1), 0);
            masks.touches_solid.assign(static_cast<std::size_t>(s.cfg.npoin + 1), 0);
            masks.velocity_active.assign(static_cast<std::size_t>(s.cfg.npoin + 1), 0);
            masks.pressure_active.assign(static_cast<std::size_t>(s.cfg.npoin + 1), 0);
            masks.thermal_active.assign(static_cast<std::size_t>(s.cfg.npoin + 1), 0);

            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                const bool fluid_element = (s.mat_elem(ie) == 0);

                for (Int in = 1; in <= s.cfg.nep; ++in)
                {
                    const Int ip = s.intma(in, ie);

                    if (ip < 1 || ip > s.cfg.npoin)
                    {
                        throw std::runtime_error(
                            "Convergence - element connectivity node out of range");
                    }

                    if (fluid_element)
                    {
                        masks.touches_fluid[static_cast<std::size_t>(ip)] = 1;
                    }
                    else
                    {
                        masks.touches_solid[static_cast<std::size_t>(ip)] = 1;
                    }
                }
            }

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                const bool fluid =
                    masks.touches_fluid[static_cast<std::size_t>(ip)] != 0;
                const bool solid =
                    masks.touches_solid[static_cast<std::size_t>(ip)] != 0;

                if (!fluid && !solid)
                {
                    throw std::runtime_error(
                        "Convergence - orphan node not touched by any material element");
                }

                // Momentum convergence is measured only where velocity is a
                // genuine flow unknown. Shared fluid-solid nodes are no-slip
                // interface nodes in the conformal CHT mesh and must not dilute
                // or mask the velocity residual.
                masks.velocity_active[static_cast<std::size_t>(ip)] =
                    static_cast<char>(fluid && !solid);

                // Pressure is a fluid-domain field. Interface nodes are still
                // pressure-active because they belong to fluid tetrahedra; only
                // solid-only nodes are excluded from the pressure norm.
                masks.pressure_active[static_cast<std::size_t>(ip)] =
                    static_cast<char>(fluid);

                // Temperature is solved in both fluid and solid. Conformal
                // shared nodes naturally carry the continuous interface value.
                masks.thermal_active[static_cast<std::size_t>(ip)] =
                    static_cast<char>(fluid || solid);
            }

            return masks;
        }


        // Tests whether node ip belongs to the selected active field space.
        bool active_at(const std::vector<char>& mask, const Int ip)
        {
            return mask[static_cast<std::size_t>(ip)] != 0;
        }


        // Accumulates one Cartesian velocity-component residual.
        //
        // For component q:
        //
        //     delta_q_i = q_i^(n+1) - q_i^n
        //
        // The three stored sums are:
        //
        //     sum_i delta_q_i^2
        //
        //     sum_i (q_i^(n+1))^2
        //
        //     sum_i [delta_q_i / elcoe2(i)]^2
        void accumulate_velocity_component(
            CBSStateSI& s,
            const MaterialNodeMasks& masks,
            const Int ip,
            const Int idim,
            const Int base)
        {
            if (!active_at(masks.velocity_active, ip))
            {
                return;
            }

            const Real old_value = s.unkn1(idim, ip);
            const Real new_value = s.unkno(idim, ip);
            const Real delta = new_value - old_value;

            s.hb[base + 0] += delta * delta;
            s.hb[base + 1] += new_value * new_value;

            // Step 1/3 velocity updates are multiplied by elcoe2. Dividing by
            // elcoe2 recovers the RHS-scale absolute residual for active flow
            // nodes only. Solid-only and conformal interface nodes are skipped.
            const Real rhs_scale = delta * safe_inverse_weight(s.elcoe2(ip));
            s.hb[base + 2] += rhs_scale * rhs_scale;
        }


        // Converts the three accumulated sums for one field component into:
        //
        //     relative residual
        //     current-field L2 norm
        //     RHS-scale residual norm
        void finalise_component(CBSStateSI& s, const Int base)
        {
            s.hb[base + 0] =
                std::sqrt(s.hb[base + 0] / (s.hb[base + 1] + eps_norm));

            s.hb[base + 1] =
                std::sqrt(s.hb[base + 1]);

            s.hb[base + 2] =
                std::sqrt(s.hb[base + 2]);
        }
    }


    //=========================================================================
    // Calculates all convergence measures after one complete CBS iteration.
    //
    // Velocity:
    //
    //     evaluated separately for u, v and w on fluid-only flow nodes
    //
    // Pressure:
    //
    //     evaluated on every fluid-connected node
    //
    // Temperature:
    //
    //     evaluated over the complete fluid-solid thermal domain
    //
    // Inputs:
    //     unkn1, pres1, temperature1   previous CBS iteration
    //     unkno, pres, temperature     current CBS iteration
    //     elcoe2                       inverse momentum diagonal
    //     elcoe2p                      inverse thermal diagonal
    //
    // Output:
    //     hb[0:14]                     convergence measures
    //=========================================================================
    void Convergence::evaluate(CBSStateSI& s)
    {
        validate_convergence_state(s);

        const MaterialNodeMasks masks = build_material_node_masks(s);

        for (Real& value : s.hb)
        {
            value = 0.0;
        }

        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            accumulate_velocity_component(s, masks, ip, 1, 0);
            accumulate_velocity_component(s, masks, ip, 2, 3);
            accumulate_velocity_component(s, masks, ip, 3, 6);

            // Pressure convergence is measured only on fluid-connected nodes.
            // Solid-only nodes carry no pressure degree of freedom and must not
            // enter either the pressure update norm or its denominator.
            if (active_at(masks.pressure_active, ip))
            {
                const Real delta = s.pres(ip) - s.pres1(ip);

                s.hb[9] += delta * delta;
                s.hb[10] += s.pres(ip) * s.pres(ip);
                s.hb[11] += delta * delta;
            }

            // Temperature convergence is measured over the whole thermal domain:
            // fluid, solid, and shared conformal interface nodes.
            if (s.cfg.temp_calc > 0 && active_at(masks.thermal_active, ip))
            {
                const Real delta = s.temperature(ip) - s.temperature1(ip);

                s.hb[12] += delta * delta;
                s.hb[13] += s.temperature(ip) * s.temperature(ip);

                const Real rhs_scale = delta * safe_inverse_weight(s.elcoe2p(ip));
                s.hb[14] += rhs_scale * rhs_scale;
            }
        }

        finalise_component(s, 0);
        finalise_component(s, 3);
        finalise_component(s, 6);

        s.hb[9] =
            std::sqrt(s.hb[9] / (s.hb[10] + eps_norm));

        s.hb[10] =
            std::sqrt(s.hb[10]);

        s.hb[11] =
            std::sqrt(s.hb[11]);

        if (s.cfg.temp_calc > 0)
        {
            s.hb[12] =
                std::sqrt(s.hb[12] / (s.hb[13] + eps_norm));

            s.hb[13] =
                std::sqrt(s.hb[13]);

            s.hb[14] =
                std::sqrt(s.hb[14]);
        }
        else
        {
            s.hb[12] = 0.0;
            s.hb[13] = 0.0;
            s.hb[14] = 0.0;
        }
    }


    //=========================================================================
    // Returns the largest relative velocity-component residual:
    //
    //     R_u = max(R_u_x, R_u_y, R_u_z)
    //=========================================================================
    Real Convergence::velocityResidual(const CBSStateSI& s)
    {
        return std::max({ s.hb[0], s.hb[3], s.hb[6] });
    }


    //=========================================================================
    // Returns the relative pressure residual.
    //=========================================================================
    Real Convergence::pressureResidual(const CBSStateSI& s)
    {
        return s.hb[9];
    }


    //=========================================================================
    // Returns the relative temperature residual.
    //=========================================================================
    Real Convergence::temperatureResidual(const CBSStateSI& s)
    {
        return s.hb[12];
    }


    //=========================================================================
    // Returns the global relative L2 residual of the SA working variable.
    //=========================================================================
    Real Convergence::turbulenceResidual(const CBSStateSI& s)
    {
        if (s.cfg.turbulence_on < 1)
        {
            return 0.0;
        }

        Real delta_sum = 0.0;
        Real value_sum = 0.0;

        // Owned nodes only.  owned_nodes is empty in a serial run, in which case
        // every local node is owned by definition.
        if (!s.owned_nodes.empty())
        {
            for (Size i = 0; i < s.owned_nodes.size(); ++i)
            {
                const Int ip = s.owned_nodes[i];

                if (s.sa_active_node(ip) == 0)
                {
                    continue;
                }

                const Real delta = s.sa_residual(ip);
                const Real value = s.nu_tilde(ip);

                delta_sum += delta * delta;
                value_sum += value * value;
            }
        }
        else
        {
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                if (s.sa_active_node(ip) == 0)
                {
                    continue;
                }

                const Real delta = s.sa_residual(ip);
                const Real value = s.nu_tilde(ip);

                delta_sum += delta * delta;
                value_sum += value * value;
            }
        }

#ifdef CBS3D_USE_MPI
        if (s.mpi_enabled && s.mpi_size > 1)
        {
            Real local[2] = { delta_sum, value_sum };
            Real global[2] = { 0.0, 0.0 };

            MPI_Allreduce(
                local,
                global,
                2,
                MPI_DOUBLE,
                MPI_SUM,
                MPI_COMM_WORLD);

            delta_sum = global[0];
            value_sum = global[1];
        }
#endif

        const Real denominator = std::sqrt(value_sum);

        if (!(denominator > 1.0e-30) || !std::isfinite(denominator))
        {
            return 0.0;
        }

        return std::sqrt(delta_sum) / denominator;
    }


    //=========================================================================
    // Collects globally reduced SA extrema for the iteration monitor.
    //=========================================================================
    Convergence::TurbulenceDiagnostics
    Convergence::turbulenceDiagnostics(const CBSStateSI& s)
    {
        TurbulenceDiagnostics d;

        if (s.cfg.turbulence_on < 1)
        {
            return d;
        }

        Real nu_tilde_min = 1.0e300;
        Real nu_tilde_max = -1.0e300;
        Real nu_t_min = 1.0e300;
        Real nu_t_max = -1.0e300;

        const bool have_owned = !s.owned_nodes.empty();
        const Size node_count = have_owned
            ? s.owned_nodes.size()
            : static_cast<Size>(s.cfg.npoin);

        for (Size i = 0; i < node_count; ++i)
        {
            const Int ip = have_owned
                ? s.owned_nodes[i]
                : static_cast<Int>(i) + 1;

            if (s.sa_active_node(ip) == 0)
            {
                continue;
            }

            nu_tilde_min = std::min(nu_tilde_min, s.nu_tilde(ip));
            nu_tilde_max = std::max(nu_tilde_max, s.nu_tilde(ip));
            nu_t_min = std::min(nu_t_min, s.nu_t(ip));
            nu_t_max = std::max(nu_t_max, s.nu_t(ip));
        }

        Real mu_t_max = -1.0e300;
        Real mu_eff_max = -1.0e300;

        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            mu_t_max = std::max(mu_t_max, s.mu_t_e(ie));
            mu_eff_max = std::max(mu_eff_max, s.mu_eff_e(ie));
        }

        // Elements are owned by exactly one rank, so no ownership filter is
        // needed for the element extrema.

#ifdef CBS3D_USE_MPI
        if (s.mpi_enabled && s.mpi_size > 1)
        {
            Real local_min[2] = { nu_tilde_min, nu_t_min };
            Real global_min[2] = { 0.0, 0.0 };

            MPI_Allreduce(
                local_min, global_min, 2, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

            nu_tilde_min = global_min[0];
            nu_t_min = global_min[1];

            Real local_max[4] =
                { nu_tilde_max, nu_t_max, mu_t_max, mu_eff_max };
            Real global_max[4] = { 0.0, 0.0, 0.0, 0.0 };

            MPI_Allreduce(
                local_max, global_max, 4, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

            nu_tilde_max = global_max[0];
            nu_t_max = global_max[1];
            mu_t_max = global_max[2];
            mu_eff_max = global_max[3];
        }
#endif

        d.nu_tilde_min = nu_tilde_min;
        d.nu_tilde_max = nu_tilde_max;
        d.nu_t_min = nu_t_min;
        d.nu_t_max = nu_t_max;
        d.mu_t_max = mu_t_max;
        d.mu_eff_max = mu_eff_max;
        d.residual = turbulenceResidual(s);

        return d;
    }


    //=========================================================================
    // Checks the enabled steady-state stopping criteria.
    //
    // Velocity criterion:
    //
    //     max(R_u, R_v, R_w) < l2norm_vel_tolerance
    //
    // Temperature criterion:
    //
    //     R_T < l2norm_temp_tolerance
    //
    // Pressure criterion:
    //
    //     R_p < l2norm_pres_tolerance
    //
    // A criterion is ignored when its corresponding check is disabled. The
    // pressure criterion is also ignored when its tolerance is non-positive.
    //=========================================================================
    bool Convergence::steadyStateReached(const CBSStateSI& s)
    {
        const bool velocity_ok =
            (s.cfg.vel_check < 1) ||
            (velocityResidual(s) < s.cfg.l2norm_vel_tolerance);

        const bool temperature_ok =
            (s.cfg.temp_check < 1) ||
            (s.cfg.temp_calc < 1) ||
            (temperatureResidual(s) < s.cfg.l2norm_temp_tolerance);

        // Pressure convergence is reported, but the existing CBS2D++_SI stop
        // rule is velocity/temperature based. Preserve that behaviour unless
        // l2norm_pres_tolerance is explicitly set positive.
        const bool pressure_ok =
            (s.cfg.l2norm_pres_tolerance <= 0.0) ||
            (pressureResidual(s) < s.cfg.l2norm_pres_tolerance);

        return velocity_ok && temperature_ok && pressure_ok;
    }
}
