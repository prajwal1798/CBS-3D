#include "cbs/assembly/SpalartAllmarasAssembly.hpp"

#include "cbs/turbulence/SpalartAllmaras.hpp"

#include <algorithm>
#include <vector>

namespace cbs
{
    namespace
    {
        bool is_fluid_element(const CBSStateSI& s, Int ie)
        {
            return s.mat_elem(ie) == 0;
        }
    }

    //=========================================================================
    // Resets all turbulence-derived properties to their laminar molecular values.
    //
    // This routine is intentionally called before any eddy-viscosity update.
    // It guarantees that the solver has a safe laminar state when turbulence is
    // disabled or when no active fluid elements are present.
    //
    // Element quantities after this routine:
    //
    //     nu_tilde_e = 0
    //     nu_t_e     = 0
    //     mu_t_e     = 0
    //     mu_eff_e   = mu_e
    //     k_eff_e    = k_e
    //
    // Nodal quantities after this routine:
    //
    //     nu_t = 0
    //     mu_t = 0
    //=========================================================================
    void SpalartAllmarasAssembly::resetEffectiveProperties(CBSStateSI& s)
    {
        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            s.nu_tilde_e(ie) = 0.0;
            s.nu_t_e(ie) = 0.0;
            s.mu_t_e(ie) = 0.0;
            s.mu_eff_e(ie) = s.mu_e(ie);
            s.k_eff_e(ie) = s.k_e(ie);
        }

        s.nu_t.fill(0.0);
        s.mu_t.fill(0.0);
    }

    //=========================================================================
    // Converts the transported SA variable into eddy viscosity and effective
    // material properties.
    //
    // Inputs:
    //
    //     nu_tilde(node)    transported SA working variable
    //     rho_e(element)    density
    //     mu_e(element)     molecular dynamic viscosity
    //     k_e(element)      molecular thermal conductivity
    //     rho_cp_e(element) volumetric heat capacity
    //
    // Outputs:
    //
    //     nu_tilde_e(element) element-averaged SA working variable
    //     nu_t_e(element)     turbulent kinematic viscosity
    //     mu_t_e(element)     turbulent dynamic viscosity
    //     mu_eff_e(element)   molecular plus turbulent dynamic viscosity
    //     k_eff_e(element)    molecular plus turbulent thermal conductivity
    //     nu_t(node)          averaged nodal turbulent kinematic viscosity
    //     mu_t(node)          averaged nodal turbulent dynamic viscosity
    //
    // Important restrictions:
    //
    //     1. Only fluid elements receive turbulent viscosity.
    //     2. Solid elements keep mu_eff_e = mu_e and k_eff_e = k_e.
    //     3. The molecular fields mu_e and k_e are never overwritten.
    //     4. If turbulence_on = 0, the routine returns after the reset.
    //
    // This routine does not solve the SA transport equation.  It only evaluates
    // properties from the current value of nu_tilde.
    //=========================================================================
    void SpalartAllmarasAssembly::updateEddyViscosity(CBSStateSI& s)
    {
        resetEffectiveProperties(s);

        if (s.cfg.turbulence_on < 1)
        {
            return;
        }

        turbulence::SpalartAllmarasConstants constants;

        std::vector<Real> nodal_nu_t_sum(
            static_cast<std::size_t>(s.cfg.npoin + 1),
            0.0);

        std::vector<Real> nodal_mu_t_sum(
            static_cast<std::size_t>(s.cfg.npoin + 1),
            0.0);

        std::vector<Int> nodal_count(
            static_cast<std::size_t>(s.cfg.npoin + 1),
            0);

        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            if (!is_fluid_element(s, ie))
            {
                continue;
            }

            Real nu_tilde_avg = 0.0;

            for (Int in = 1; in <= s.cfg.nep; ++in)
            {
                const Int ip = s.intma(in, ie);
                const Real limited_value = std::max(
                    s.cfg.sa_nu_tilde_floor,
                    s.nu_tilde(ip));

                nu_tilde_avg += limited_value;
            }

            nu_tilde_avg /= static_cast<Real>(s.cfg.nep);

            Real molecular_nu = 0.0;
            if (s.rho_e(ie) > 0.0)
            {
                molecular_nu = s.mu_e(ie) / s.rho_e(ie);
            }

            Real nu_t = 0.0;
            if (molecular_nu > 0.0)
            {
                nu_t = turbulence::eddyKinematicViscosity(
                    nu_tilde_avg,
                    molecular_nu,
                    constants);
            }

            const Real mu_t = s.rho_e(ie) * nu_t;

            s.nu_tilde_e(ie) = nu_tilde_avg;
            s.nu_t_e(ie) = nu_t;
            s.mu_t_e(ie) = mu_t;
            s.mu_eff_e(ie) = s.mu_e(ie) + mu_t;

            if (s.cfg.turbulent_thermal_diffusivity_on > 0)
            {
                s.k_eff_e(ie) = s.k_e(ie)
                    + s.rho_cp_e(ie) * nu_t / s.cfg.sa_prandtl_t;
            }

            for (Int in = 1; in <= s.cfg.nep; ++in)
            {
                const Int ip = s.intma(in, ie);

                nodal_nu_t_sum[static_cast<std::size_t>(ip)] += nu_t;
                nodal_mu_t_sum[static_cast<std::size_t>(ip)] += mu_t;
                ++nodal_count[static_cast<std::size_t>(ip)];
            }
        }

        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            const Int count = nodal_count[static_cast<std::size_t>(ip)];

            if (count > 0)
            {
                s.nu_t(ip) = nodal_nu_t_sum[static_cast<std::size_t>(ip)]
                    / static_cast<Real>(count);

                s.mu_t(ip) = nodal_mu_t_sum[static_cast<std::size_t>(ip)]
                    / static_cast<Real>(count);
            }
        }
    }
}
