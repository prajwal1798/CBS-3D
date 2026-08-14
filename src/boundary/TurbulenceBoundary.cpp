#include "cbs/boundary/TurbulenceBoundary.hpp"

#include "cbs/parallel/HaloExchange.hpp"

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

        bool is_wall_bc(const CBSStateSI& s, Int bc)
        {
            return bc == s.cfg.bc_noslip_adiabatic_wall
                || bc == s.cfg.bc_noslip_heatflux_wall
                || bc == s.cfg.bc_cht_interface;
        }

        bool is_inlet_bc(const CBSStateSI& s, Int bc)
        {
            return bc == s.cfg.bc_velocity_temperature_inlet
                || bc == s.cfg.bc_massflow_temperature_inlet;
        }

        //=====================================================================
        // Returns an averaged molecular kinematic viscosity at one node.
        //
        // The SA inlet value is prescribed as
        //
        //     nu_tilde_inlet = sa_inlet_ratio * nu
        //
        // where nu is the molecular kinematic viscosity.  Material data are
        // stored per element, while the SA working variable is nodal.  Therefore
        // the nodal molecular viscosity is obtained by averaging all fluid
        // elements touching the node.
        //
        // This is a simple serial implementation.  It is acceptable for the
        // initial turbulence scaffold and will be replaced by a precomputed nodal
        // material field when the full SA equation is coupled to the solver.
        //=====================================================================
        Real nodal_molecular_nu(
            const CBSStateSI& s,
            Int ip,
            const std::vector<Int>& fluid_touch_count)
        {
            Real sum = 0.0;

            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                if (!is_fluid_element(s, ie))
                {
                    continue;
                }

                bool element_touches_node = false;

                for (Int in = 1; in <= s.cfg.nep; ++in)
                {
                    if (s.intma(in, ie) == ip)
                    {
                        element_touches_node = true;
                    }
                }

                if (element_touches_node && s.rho_e(ie) > 0.0)
                {
                    sum += s.mu_e(ie) / s.rho_e(ie);
                }
            }

            const Int count = fluid_touch_count[static_cast<std::size_t>(ip)];

            if (count > 0)
            {
                return sum / static_cast<Real>(count);
            }

            return 0.0;
        }
    }

    //=========================================================================
    // Classifies nodes for the Spalart-Allmaras boundary treatment.
    //
    // A node is SA-active when it belongs to at least one fluid element.  The SA
    // equation is not solved in pure solid regions.
    //
    // Wall nodes are nodes on no-slip physical walls and on the conformal CHT
    // interface.  These nodes receive
    //
    //     nu_tilde = 0
    //
    // because the SA working variable vanishes at no-slip walls.
    //
    // Inlet nodes receive a prescribed fully turbulent inlet value.  Partition
    // interfaces are not physical boundaries and are not classified here.
    //=========================================================================
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
                const Int ip = s.intma(in, ie);
                s.sa_active_node(ip) = 1;
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

    //=========================================================================
    // Reduces the three SA node flags over every copy of a shared node.
    //
    // Each flag is written into a real buffer, summed onto the owning rank and
    // broadcast back to the ghost copies.  A non-zero total means at least one
    // rank classified the node, so the reduction is a logical OR.  Summing and
    // then thresholding is used because the halo exchange is defined for real
    // arrays, and an OR is what is wanted here rather than a sum.
    //=========================================================================
    void TurbulenceBoundary::synchroniseClassification(CBSStateSI& s)
    {
#ifdef CBS3D_USE_MPI
        if (!s.mpi_enabled || s.mpi_size <= 1)
        {
            return;
        }

        Array1D<Int>* flags[3] =
        {
            &s.sa_active_node,
            &s.sa_wall_node,
            &s.sa_inlet_node
        };

        Array1D<Real> buffer(s.cfg.npoin);

        for (Int f = 0; f < 3; ++f)
        {
            Array1D<Int>& flag = *flags[f];

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                buffer(ip) = flag(ip) != 0 ? 1.0 : 0.0;
            }

            HaloExchange::sumGhostContributionsToOwners(
                buffer,
                s.partition_metadata,
                MPI_COMM_WORLD);

            HaloExchange::broadcastOwnedToGhosts(
                buffer,
                s.partition_metadata,
                MPI_COMM_WORLD);

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                flag(ip) = buffer(ip) > 0.5 ? 1 : 0;
            }
        }
#else
        (void)s;
#endif
    }


    //=========================================================================
    // Initialises the transported SA working variable.
    //
    // The first value is deliberately simple:
    //
    //     nu_tilde = sa_inlet_ratio * nu
    //
    // for active non-wall fluid nodes.  Wall and solid-only nodes are kept at
    // zero.  Inlet nodes are then overwritten by the same prescribed inlet rule
    // to make the boundary condition explicit.
    //=========================================================================
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

        std::vector<Int> fluid_touch_count(
            static_cast<std::size_t>(s.cfg.npoin + 1),
            0);

        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            if (!is_fluid_element(s, ie))
            {
                continue;
            }

            for (Int in = 1; in <= s.cfg.nep; ++in)
            {
                const Int ip = s.intma(in, ie);
                ++fluid_touch_count[static_cast<std::size_t>(ip)];
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

    //=========================================================================
    // Applies the wall condition for the SA working variable.
    //
    // At all no-slip turbulence walls:
    //
    //     nu_tilde = 0
    //     nu_t     = 0
    //     mu_t     = 0
    //
    // The CHT interface is included here for turbulence because the fluid sees a
    // no-slip solid wall.  This does not alter the thermal continuity treatment.
    //=========================================================================
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

    //=========================================================================
    // Applies the inlet condition for the SA working variable.
    //
    // The prescribed value is
    //
    //     nu_tilde_inlet = sa_inlet_ratio * nu
    //
    // Wall nodes are not overwritten by the inlet rule.  If a node is both an
    // inlet node and a wall node due to boundary-face adjacency, the wall value
    // remains dominant.
    //=========================================================================
    void TurbulenceBoundary::applyInletValues(CBSStateSI& s)
    {
        if (s.cfg.turbulence_on < 1)
        {
            return;
        }

        std::vector<Int> fluid_touch_count(
            static_cast<std::size_t>(s.cfg.npoin + 1),
            0);

        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            if (!is_fluid_element(s, ie))
            {
                continue;
            }

            for (Int in = 1; in <= s.cfg.nep; ++in)
            {
                const Int ip = s.intma(in, ie);
                ++fluid_touch_count[static_cast<std::size_t>(ip)];
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
