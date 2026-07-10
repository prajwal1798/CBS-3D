#include "cbs/boundary/TurbulenceBoundary.hpp"

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

        [[nodiscard]] bool is_wall_bc(const CBSStateSI& s, const Int bc)
        {
            return bc == s.cfg.bc_noslip_adiabatic_wall
                || bc == s.cfg.bc_noslip_heatflux_wall
                || bc == s.cfg.bc_cht_interface;
        }

        [[nodiscard]] bool is_inlet_bc(const CBSStateSI& s, const Int bc)
        {
            return bc == s.cfg.bc_velocity_temperature_inlet
                || bc == s.cfg.bc_massflow_temperature_inlet;
        }

        [[nodiscard]] Real nodal_molecular_nu(
            const CBSStateSI& s,
            const Int ip,
            const std::vector<Int>& count)
        {
            Real sum = 0.0;

            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                if (!is_fluid_element(s, ie))
                {
                    continue;
                }

                bool touches_node = false;
                for (Int in = 1; in <= s.cfg.nep; ++in)
                {
                    touches_node = touches_node || (s.intma(in, ie) == ip);
                }

                if (touches_node && s.rho_e(ie) > 0.0)
                {
                    sum += s.mu_e(ie) / s.rho_e(ie);
                }
            }

            const Int n = count[static_cast<std::size_t>(ip)];
            return n > 0 ? sum / static_cast<Real>(n) : 0.0;
        }
    }

    void TurbulenceBoundary::classifyNodes(CBSStateSI& s)
    {
        s.sa_active_node.fill(0);
        s.sa_wall_node.fill(0);
        s.sa_inlet_node.fill(0);

        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            if (!is_fluid_element(s, ie))
            {
                continue;
            }

            for (Int in = 1; in <= s.cfg.nep; ++in)
            {
                s.sa_active_node(s.intma(in, ie)) = 1;
            }
        }

        for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
        {
            const Int bc = s.iside(s.cfg.bsid, ib);

            if (!is_wall_bc(s, bc) && !is_inlet_bc(s, bc))
            {
                continue;
            }

            for (Int i = 1; i <= s.cfg.nsidp; ++i)
            {
                const Int ip = s.iside(i, ib);

                if (is_wall_bc(s, bc))
                {
                    s.sa_wall_node(ip) = 1;
                }

                if (is_inlet_bc(s, bc))
                {
                    s.sa_inlet_node(ip) = 1;
                }
            }
        }
    }

    void TurbulenceBoundary::initialiseNuTilde(CBSStateSI& s)
    {
        classifyNodes(s);

        s.nu_tilde.fill(0.0);
        s.nu_tilde1.fill(0.0);
        s.nu_t.fill(0.0);
        s.mu_t.fill(0.0);

        if (s.cfg.turbulence_on < 1)
        {
            return;
        }

        std::vector<Int> fluid_touch_count(static_cast<std::size_t>(s.cfg.npoin + 1), 0);

        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            if (!is_fluid_element(s, ie))
            {
                continue;
            }

            for (Int in = 1; in <= s.cfg.nep; ++in)
            {
                ++fluid_touch_count[static_cast<std::size_t>(s.intma(in, ie))];
            }
        }

        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            if (s.sa_active_node(ip) == 0 || s.sa_wall_node(ip) != 0)
            {
                continue;
            }

            const Real nu = nodal_molecular_nu(s, ip, fluid_touch_count);
            s.nu_tilde(ip) = std::max(0.0, s.cfg.sa_inlet_ratio * nu);
            s.nu_tilde1(ip) = s.nu_tilde(ip);
        }

        applyWallValues(s);
        applyInletValues(s);
    }

    void TurbulenceBoundary::applyWallValues(CBSStateSI& s)
    {
        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            if (s.sa_wall_node(ip) != 0)
            {
                s.nu_tilde(ip) = 0.0;
                s.nu_tilde1(ip) = 0.0;
                s.nu_t(ip) = 0.0;
                s.mu_t(ip) = 0.0;
            }
        }
    }

    void TurbulenceBoundary::applyInletValues(CBSStateSI& s)
    {
        if (s.cfg.turbulence_on < 1)
        {
            return;
        }

        std::vector<Int> fluid_touch_count(static_cast<std::size_t>(s.cfg.npoin + 1), 0);

        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            if (!is_fluid_element(s, ie))
            {
                continue;
            }

            for (Int in = 1; in <= s.cfg.nep; ++in)
            {
                ++fluid_touch_count[static_cast<std::size_t>(s.intma(in, ie))];
            }
        }

        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            if (s.sa_inlet_node(ip) != 0 && s.sa_wall_node(ip) == 0)
            {
                const Real nu = nodal_molecular_nu(s, ip, fluid_touch_count);
                s.nu_tilde(ip) = std::max(0.0, s.cfg.sa_inlet_ratio * nu);
                s.nu_tilde1(ip) = s.nu_tilde(ip);
            }
        }
    }
}
