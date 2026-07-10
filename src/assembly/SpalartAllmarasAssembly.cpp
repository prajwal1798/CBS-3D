#include "cbs/assembly/SpalartAllmarasAssembly.hpp"

#include "cbs/turbulence/SpalartAllmaras.hpp"

#include <algorithm>
#include <vector>

namespace cbs
{
    namespace
    {
        [[nodiscard]] bool is_fluid_element(const CBSStateSI& s, const Int ie)
        {
            return s.mat_elem(ie) == 0;
        }
    }

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

    void SpalartAllmarasAssembly::updateEddyViscosity(CBSStateSI& s)
    {
        resetEffectiveProperties(s);

        if (s.cfg.turbulence_on < 1)
        {
            return;
        }

        std::vector<Real> nodal_nu_t_sum(static_cast<std::size_t>(s.cfg.npoin + 1), 0.0);
        std::vector<Real> nodal_mu_t_sum(static_cast<std::size_t>(s.cfg.npoin + 1), 0.0);
        std::vector<Int> nodal_count(static_cast<std::size_t>(s.cfg.npoin + 1), 0);

        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            if (!is_fluid_element(s, ie))
            {
                continue;
            }

            Real nu_tilde_avg = 0.0;
            for (Int in = 1; in <= s.cfg.nep; ++in)
            {
                nu_tilde_avg += std::max(s.cfg.sa_nu_tilde_floor, s.nu_tilde(s.intma(in, ie)));
            }
            nu_tilde_avg /= static_cast<Real>(s.cfg.nep);

            const Real molecular_nu = s.rho_e(ie) > 0.0
                ? s.mu_e(ie) / s.rho_e(ie)
                : 0.0;

            const Real nu_t = molecular_nu > 0.0
                ? turbulence::eddyKinematicViscosity(nu_tilde_avg, molecular_nu)
                : 0.0;

            const Real mu_t = s.rho_e(ie) * nu_t;

            s.nu_tilde_e(ie) = nu_tilde_avg;
            s.nu_t_e(ie) = nu_t;
            s.mu_t_e(ie) = mu_t;
            s.mu_eff_e(ie) = s.mu_e(ie) + mu_t;

            if (s.cfg.turbulent_thermal_diffusivity_on > 0)
            {
                s.k_eff_e(ie) = s.k_e(ie) + s.rho_cp_e(ie) * nu_t / s.cfg.sa_prandtl_t;
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
            const Int n = nodal_count[static_cast<std::size_t>(ip)];
            if (n > 0)
            {
                s.nu_t(ip) = nodal_nu_t_sum[static_cast<std::size_t>(ip)] / static_cast<Real>(n);
                s.mu_t(ip) = nodal_mu_t_sum[static_cast<std::size_t>(ip)] / static_cast<Real>(n);
            }
        }
    }
}
