//=============================================================================
// CBS3D++_SI
//
// Boundary-condition application for the three-dimensional semi-implicit
// Characteristic-Based Split finite-element solver.
//
// Strong Dirichlet conditions are applied directly to nodal velocity,
// pressure and temperature values. Symmetry and outlet conditions are applied
// by projecting the velocity onto or away from the boundary normal.
//
// The module also enforces zero velocity at solid and conformal fluid-solid
// interface nodes in conjugate heat-transfer calculations.
//=============================================================================

#include "cbs/boundary/Boundary.hpp"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace cbs
{
    namespace
    {
        //-------------------------------------------------------------------------
        // Returns true when the boundary identifier is implemented by this
        // module. The actual numerical action associated with each identifier
        // is applied in the public boundary routines below.
        //-------------------------------------------------------------------------
        bool supported_bc(const CBSStateSI& s, Int bc)
        {
            return
                bc == s.cfg.bc_adiabatic_prescribed_velocity ||
                bc == s.cfg.bc_temperature_one_noslip ||
                bc == s.cfg.bc_temperature_zero_noslip ||
                bc == s.cfg.bc_temperature_zero_prescribed_velocity ||
                bc == s.cfg.bc_pressure ||
                bc == s.cfg.bc_symmetry_no_flux ||
                bc == s.cfg.bc_bfs_parabolic_inlet ||
                bc == s.cfg.bc_parabolic_inlet ||
                bc == s.cfg.bc_velocity_temperature_inlet ||
                bc == s.cfg.bc_massflow_temperature_inlet ||
                bc == s.cfg.bc_pressure_outlet ||
                bc == s.cfg.bc_noslip_adiabatic_wall ||
                bc == s.cfg.bc_noslip_heatflux_wall ||
                bc == s.cfg.bc_cht_interface ||
                bc == s.cfg.bc_heatflux_marker;
        }

        //-------------------------------------------------------------------------
        // Checks every mapped boundary face before applying any condition.
        // This prevents an unknown boundary identifier from being silently
        // ignored during the numerical calculation.
        //-------------------------------------------------------------------------
        void validate_boundary_flags(const CBSStateSI& s)
        {
            for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
            {
                const Int bc = s.iside(s.cfg.bsid, ib);

                if (!supported_bc(s, bc))
                {
                    throw std::runtime_error(
                        "Boundary - unsupported BC_ID "
                        + std::to_string(bc)
                        + " at boundary face "
                        + std::to_string(ib)
                        + ". This CBS3D boundary module supports the CBS2D++_SI BC family adapted to 3D.");
                }
            }
        }

        //-------------------------------------------------------------------------
        // Verifies the one-based global node index used by the solver arrays.
        //-------------------------------------------------------------------------
        void check_node_range(const CBSStateSI& s, Int ip, const char* context)
        {
            if (ip < 1 || ip > s.cfg.npoin)
            {
                throw std::runtime_error(
                    std::string(context) + " - node index out of range");
            }
        }

        //-------------------------------------------------------------------------
        // Recovers the outward unit normal from the area-weighted face normal:
        //
        //     n = (A_f n) / A_f
        //
        // where face_norm(1:3,ib) stores A_f n and face_norm(4,ib)
        // stores the triangular face area A_f.
        //-------------------------------------------------------------------------
        void unit_outward_normal(
            const CBSStateSI& s,
            Int ib,
            Real& nx,
            Real& ny,
            Real& nz)
        {
            const Real area = s.face_norm(4, ib);

            if (area <= 0.0 || !std::isfinite(area))
            {
                throw std::runtime_error(
                    "Boundary - invalid boundary-face area in face_norm(4,ib). "
                    "Preprocess::getNormals must be called before boundary application.");
            }

            nx = s.face_norm(1, ib) / area;
            ny = s.face_norm(2, ib) / area;
            nz = s.face_norm(3, ib) / area;
        }

        //-------------------------------------------------------------------------
        // Assigns the three Cartesian velocity components at one mesh node.
        //-------------------------------------------------------------------------
        void set_velocity(
            CBSStateSI& s,
            Int ip,
            Real u,
            Real v,
            Real w)
        {
            check_node_range(s, ip, "Boundary::set_velocity");

            s.unkno(1, ip) = u;
            s.unkno(2, ip) = v;
            s.unkno(3, ip) = w;
        }

        //-------------------------------------------------------------------------
        // Applies the no-slip condition:
        //
        //     u = v = w = 0
        //-------------------------------------------------------------------------
        void set_zero_velocity(CBSStateSI& s, Int ip)
        {
            set_velocity(s, ip, 0.0, 0.0, 0.0);
        }


        //-------------------------------------------------------------------------
        // Identifies pressure boundaries on which no strong velocity value is
        // imposed by applyVelocity().
        //-------------------------------------------------------------------------
        bool is_pressure_only_bc(const CBSStateSI& s, Int bc)
        {
            return
                bc == s.cfg.bc_pressure ||
                bc == s.cfg.bc_pressure_outlet;
        }

        //-------------------------------------------------------------------------
        // Legacy backward-facing-step inlet profile:
        //
        //     u(y) = 0.6624 y^6 - 7.5547 y^5 + 33.9 y^4
        //          - 75.283 y^3 + 83.368 y^2 - 37.793 y + 2.6959
        //-------------------------------------------------------------------------
        Real bfs_profile_u(Real y)
        {
            return
                0.6624 * std::pow(y, 6)
                - 7.5547 * std::pow(y, 5)
                + 33.9 * std::pow(y, 4)
                - 75.283 * std::pow(y, 3)
                + 83.368 * std::pow(y, 2)
                - 37.793 * y
                + 2.6959;
        }

        //-------------------------------------------------------------------------
        // Parabolic inlet profile for the rectangular benchmark:
        //
        //     u(y) = 6 y (1 - y)
        //
        // The profile is zero at y = 0 and y = 1 and has unit mean value over
        // the interval 0 <= y <= 1.
        //-------------------------------------------------------------------------
        Real rectangular_profile_u(Real y)
        {
            return 6.0 * y * (1.0 - y);
        }

        //-------------------------------------------------------------------------
        // Determines whether each global node belongs to at least one fluid
        // element and/or at least one solid element.
        //
        // A conformal CHT interface node is characterised by:
        //
        //     touches_fluid(ip) = true
        //     touches_solid(ip) = true
        //-------------------------------------------------------------------------
        void build_node_material_touch_masks(
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
                    check_node_range(s, ip, "Boundary::build_node_material_touch_masks");

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

        //-------------------------------------------------------------------------
        // Removes every node touched by a solid element from the free velocity
        // space:
        //
        //     u(ip) = 0
        //
        // This includes both solid-only nodes and shared fluid-solid interface
        // nodes. Temperature remains continuous and is not modified here.
        //-------------------------------------------------------------------------
        void enforce_zero_velocity_on_material_solid_touch_nodes(CBSStateSI& s)
        {
            std::vector<char> touches_fluid;
            std::vector<char> touches_solid;
            build_node_material_touch_masks(s, touches_fluid, touches_solid);

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                const bool solid =
                    touches_solid[static_cast<std::size_t>(ip)] != 0;

                if (solid)
                {
                    // In conformal CHT meshes, every node touched by a solid
                    // element is outside the free velocity space.  This covers
                    // both solid-only nodes and shared fluid-solid interface
                    // nodes.  Temperature remains unconstrained here.
                    set_zero_velocity(s, ip);
                }
            }
        }
    }

    //=========================================================================
    // Applies the impermeability condition on symmetry or slip boundaries.
    //
    // For velocity u and outward unit normal n, the normal component is:
    //
    //     u_n = u . n
    //
    // The corrected velocity is:
    //
    //     u_new = u - u_n n
    //
    // Hence:
    //
    //     u_new . n = 0
    //
    // while both tangential velocity components remain unchanged.
    //
    // Modified array:
    //     s.unkno(1:3,ip)
    //=========================================================================
    void Boundary::applySymmetry(CBSStateSI& s)
    {
        validate_boundary_flags(s);

        for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
        {
            const Int bc = s.iside(s.cfg.bsid, ib);

            if (bc != s.cfg.bc_symmetry_no_flux)
            {
                continue;
            }

            Real nx = 0.0;
            Real ny = 0.0;
            Real nz = 0.0;
            unit_outward_normal(s, ib, nx, ny, nz);

            for (Int in = 1; in <= s.cfg.nsidp; ++in)
            {
                const Int ip = s.iside(in, ib);
                check_node_range(s, ip, "Boundary::applySymmetry");

                const Real un =
                    s.unkno(1, ip) * nx +
                    s.unkno(2, ip) * ny +
                    s.unkno(3, ip) * nz;

                // Symmetry/slip boundary: remove only the normal velocity
                // component. Tangential components remain untouched.
                s.unkno(1, ip) -= un * nx;
                s.unkno(2, ip) -= un * ny;
                s.unkno(3, ip) -= un * nz;
            }
        }
    }

    //=========================================================================
    // Applies all strong nodal velocity boundary conditions.
    //
    // Prescribed velocity:
    //
    //     u = u_bc
    //
    // No-slip wall, solid and conformal CHT interface:
    //
    //     u = 0
    //
    // For a mass-flow inlet, the scalar inlet speed is calculated during
    // preprocessing and directed into the domain:
    //
    //     u_in = -U_in n
    //
    // where n is the outward unit normal. Pressure and symmetry boundaries do
    // not receive a strong velocity value in this routine.
    //
    // The final priority is:
    //
    //     ordinary face conditions
    //         -> wall-node no-slip
    //         -> moving-wall BC 500
    //         -> material-domain no-slip
    //
    // Modified array:
    //     s.unkno(1:3,ip)
    //=========================================================================
    void Boundary::applyVelocity(CBSStateSI& s)
    {
        validate_boundary_flags(s);

        for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
        {
            const Int bc = s.iside(s.cfg.bsid, ib);

            Real nx = 0.0;
            Real ny = 0.0;
            Real nz = 0.0;

            if (bc == s.cfg.bc_massflow_temperature_inlet &&
                s.cfg.mass_flow_inlet_enabled > 0)
            {
                // Mass-flow inlet in 3D uses the inward normal direction.
                // face_norm is outward from the domain, so inflow is -n.
                unit_outward_normal(s, ib, nx, ny, nz);
            }

            for (Int in = 1; in <= s.cfg.nsidp; ++in)
            {
                const Int ip = s.iside(in, ib);
                check_node_range(s, ip, "Boundary::applyVelocity");

                // -------------------------------------------------------------
                // 500 : adiabatic with prescribed velocity.
                // -------------------------------------------------------------
                if (bc == s.cfg.bc_adiabatic_prescribed_velocity)
                {
                    set_velocity(s, ip, s.cfg.inlet_u, s.cfg.inlet_v, s.cfg.inlet_w);
                }

                // -------------------------------------------------------------
                // 501 / 502 : constant-temperature no-slip walls.
                // -------------------------------------------------------------
                else if (bc == s.cfg.bc_temperature_one_noslip ||
                         bc == s.cfg.bc_temperature_zero_noslip)
                {
                    set_zero_velocity(s, ip);
                }

                // -------------------------------------------------------------
                // 503 : legacy prescribed x-velocity inlet.
                // -------------------------------------------------------------
                else if (bc == s.cfg.bc_temperature_zero_prescribed_velocity)
                {
                    set_velocity(s, ip, 1.0, 0.0, 0.0);
                }

                // -------------------------------------------------------------
                // 504 / 520 : pressure boundary / pressure outlet.
                // No strong velocity is imposed.
                // -------------------------------------------------------------
                else if (is_pressure_only_bc(s, bc))
                {
                    continue;
                }

                // -------------------------------------------------------------
                // 506 : symmetry/slip. The velocity projection is handled in
                // applySymmetry().
                // -------------------------------------------------------------
                else if (bc == s.cfg.bc_symmetry_no_flux)
                {
                    continue;
                }

                // -------------------------------------------------------------
                // 507 : legacy backward-facing-step parabolic inlet.
                // The original profile is x-directed and varies with y.
                // -------------------------------------------------------------
                else if (bc == s.cfg.bc_bfs_parabolic_inlet)
                {
                    set_velocity(s, ip, bfs_profile_u(s.coord(2, ip)), 0.0, 0.0);
                }

                // -------------------------------------------------------------
                // 508 : rectangular-channel parabolic inlet.
                // The original benchmark profile is x-directed and varies with y.
                // -------------------------------------------------------------
                else if (bc == s.cfg.bc_parabolic_inlet)
                {
                    set_velocity(s, ip, rectangular_profile_u(s.coord(2, ip)), 0.0, 0.0);
                }

                // -------------------------------------------------------------
                // 510 : prescribed 3D inlet velocity + prescribed temperature.
                // -------------------------------------------------------------
                else if (bc == s.cfg.bc_velocity_temperature_inlet)
                {
                    set_velocity(s, ip, s.cfg.inlet_u, s.cfg.inlet_v, s.cfg.inlet_w);
                }

                // -------------------------------------------------------------
                // 511 : mass-flow inlet + prescribed temperature.
                //
                // The scalar speed is computed from the inlet area in
                // Preprocess::computeMassFlowInletVelocity().
                // -------------------------------------------------------------
                else if (bc == s.cfg.bc_massflow_temperature_inlet)
                {
                    if (s.cfg.mass_flow_inlet_enabled > 0)
                    {
                        set_velocity(
                            s,
                            ip,
                            -s.cfg.inlet_u_from_massflow * nx,
                            -s.cfg.inlet_u_from_massflow * ny,
                            -s.cfg.inlet_u_from_massflow * nz);
                    }
                    else
                    {
                        // Flow-only benchmarks, such as the turbulent flat plate,
                        // often use BC 511 as a uniform velocity inlet without a
                        // prescribed mass-flow rate.  In that case the velocity
                        // components are read directly from the .par file.
                        set_velocity(
                            s,
                            ip,
                            s.cfg.inlet_u,
                            s.cfg.inlet_v,
                            s.cfg.inlet_w);
                    }
                }

                // -------------------------------------------------------------
                // 530 / 532 / 901 : no-slip wall or CHT interface.
                // -------------------------------------------------------------
                else if (bc == s.cfg.bc_noslip_adiabatic_wall ||
                         bc == s.cfg.bc_noslip_heatflux_wall ||
                         bc == s.cfg.bc_cht_interface)
                {
                    set_zero_velocity(s, ip);
                }

                // -------------------------------------------------------------
                // 902 : heat-flux marker only.
                // It is not a velocity boundary and not a CHT interface marker.
                // -------------------------------------------------------------
                else if (bc == s.cfg.bc_heatflux_marker)
                {
                    continue;
                }
            }
        }

        // Wall-node list is applied after ordinary boundary faces so strict
        // no-slip walls/interfaces remain enforced.  BC 500 is intentionally
        // not part of wall_node_list.
        for (Int iw = 1; iw <= s.cfg.npoin_wall; ++iw)
        {
            const Int ip = s.wall_node_list(iw);

            if (ip >= 1 && ip <= s.cfg.npoin)
            {
                set_zero_velocity(s, ip);
            }
        }

        // Re-apply BC 500 after the no-slip wall-node pass.  This gives the
        // moving lid a clear priority in lid-driven-cavity cases where lid
        // nodes also lie on geometric wall edges.
        for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
        {
            const Int bc = s.iside(s.cfg.bsid, ib);

            if (bc != s.cfg.bc_adiabatic_prescribed_velocity)
            {
                continue;
            }

            for (Int in = 1; in <= s.cfg.nsidp; ++in)
            {
                const Int ip = s.iside(in, ib);
                set_velocity(s, ip, s.cfg.inlet_u, s.cfg.inlet_v, s.cfg.inlet_w);
            }
        }

        // Material-domain enforcement is the final velocity priority for CHT.
        // The conformal fluid-solid interface is usually internal and therefore
        // absent from the .plt boundary-face block.  No-slip at that interface
        // must therefore be enforced from material adjacency, not from BC 901.
        enforce_zero_velocity_on_material_solid_touch_nodes(s);
    }

    //=========================================================================
    // Applies strong Dirichlet temperature boundary conditions.
    //
    // A prescribed temperature is imposed directly as:
    //
    //     T(ip) = T_bc
    //
    // Adiabatic walls, pressure boundaries, heat-flux walls and conformal CHT
    // interfaces are not modified here. Their thermal behaviour is handled by
    // the energy-equation assembly:
    //
    //     -k grad(T) . n = q''
    //
    // For a conformal interface, shared mesh nodes provide temperature
    // continuity without imposing a separate nodal value.
    //
    // Modified array:
    //     s.temperature(ip)
    //=========================================================================
    void Boundary::applyTemperature(CBSStateSI& s)
    {
        validate_boundary_flags(s);

        for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
        {
            const Int bc = s.iside(s.cfg.bsid, ib);

            for (Int in = 1; in <= s.cfg.nsidp; ++in)
            {
                const Int ip = s.iside(in, ib);
                check_node_range(s, ip, "Boundary::applyTemperature");

                // -------------------------------------------------------------
                // 501 : constant nondimensional hot wall.
                // -------------------------------------------------------------
                if (bc == s.cfg.bc_temperature_one_noslip)
                {
                    s.temperature(ip) = 1.0;
                }

                // -------------------------------------------------------------
                // 502 : constant nondimensional cold wall.
                // -------------------------------------------------------------
                else if (bc == s.cfg.bc_temperature_zero_noslip)
                {
                    s.temperature(ip) = 0.0;
                }

                // -------------------------------------------------------------
                // 503 : legacy prescribed-velocity inlet with T = 0.
                // -------------------------------------------------------------
                else if (bc == s.cfg.bc_temperature_zero_prescribed_velocity)
                {
                    s.temperature(ip) = 0.0;
                }

                // -------------------------------------------------------------
                // 508 : rectangular parabolic benchmark inlet with T = 0.
                // -------------------------------------------------------------
                else if (bc == s.cfg.bc_parabolic_inlet)
                {
                    s.temperature(ip) = 0.0;
                }

                // -------------------------------------------------------------
                // 510 / 511 : prescribed-temperature inlet.
                // -------------------------------------------------------------
                else if (bc == s.cfg.bc_velocity_temperature_inlet ||
                         bc == s.cfg.bc_massflow_temperature_inlet)
                {
                    s.temperature(ip) = s.cfg.inlet_temperature;
                }

                // -------------------------------------------------------------
                // 500, 504, 506, 507, 520, 530:
                // no strong thermal value.
                // -------------------------------------------------------------
                else if (bc == s.cfg.bc_adiabatic_prescribed_velocity ||
                         bc == s.cfg.bc_pressure ||
                         bc == s.cfg.bc_symmetry_no_flux ||
                         bc == s.cfg.bc_bfs_parabolic_inlet ||
                         bc == s.cfg.bc_pressure_outlet ||
                         bc == s.cfg.bc_noslip_adiabatic_wall)
                {
                    continue;
                }

                // -------------------------------------------------------------
                // 532 : heat-flux wall.
                // Neumann heat flux must be assembled in EnergyAssembly.
                // -------------------------------------------------------------
                else if (bc == s.cfg.bc_noslip_heatflux_wall)
                {
                    continue;
                }

                // -------------------------------------------------------------
                // 901 : conformal CHT interface.
                // Do not impose temperature strongly; shared nodes enforce
                // temperature continuity and material properties define flux.
                // -------------------------------------------------------------
                else if (bc == s.cfg.bc_cht_interface)
                {
                    continue;
                }

                // -------------------------------------------------------------
                // 902 : heat-flux marker only.
                // It must not be treated as a CHT interface or fixed temperature.
                // -------------------------------------------------------------
                else if (bc == s.cfg.bc_heatflux_marker)
                {
                    continue;
                }
            }
        }
    }

    //=========================================================================
    // Applies prescribed nodal pressure values.
    //
    // For every fixed-pressure node stored during preprocessing:
    //
    //     p(ip) = p_bc
    //
    // Inputs:
    //     s.bc_list(1:bc_fixed)
    //     s.bc_values(1:bc_fixed)
    //
    // Modified array:
    //     s.pres(ip)
    //=========================================================================
    void Boundary::applyPressure(CBSStateSI& s)
    {
        for (Int i = 1; i <= s.cfg.bc_fixed; ++i)
        {
            const Int ip = s.bc_list(i);

            check_node_range(s, ip, "Boundary::applyPressure");

            s.pres(ip) = s.bc_values(i);
        }
    }

    //=========================================================================
    // Prevents inward flow through a pressure-outlet boundary.
    //
    // The normal velocity is:
    //
    //     u_n = u . n
    //
    // Positive u_n denotes flow leaving the domain. When u_n < 0, the inward
    // normal component is removed:
    //
    //     u_new = u - u_n n
    //
    // Tangential velocity is retained. Nodes touched by solid elements remain
    // no-slip and are excluded from the outlet correction.
    //
    // Modified array:
    //     s.unkno(1:3,ip)
    //=========================================================================
    void Boundary::applyOutletBackflowControl(CBSStateSI& s)
    {
        validate_boundary_flags(s);

        std::vector<char> touches_fluid;
        std::vector<char> touches_solid;
        build_node_material_touch_masks(s, touches_fluid, touches_solid);

        for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
        {
            const Int bc = s.iside(s.cfg.bsid, ib);

            if (bc != s.cfg.bc_pressure_outlet)
            {
                continue;
            }

            Real nx = 0.0;
            Real ny = 0.0;
            Real nz = 0.0;
            unit_outward_normal(s, ib, nx, ny, nz);

            for (Int in = 1; in <= s.cfg.nsidp; ++in)
            {
                const Int ip = s.iside(in, ib);
                check_node_range(s, ip, "Boundary::applyOutletBackflowControl");

                if (touches_solid[static_cast<std::size_t>(ip)] != 0)
                {
                    // Outlet backflow treatment belongs only to the free fluid
                    // velocity space.  Interface/solid nodes remain no-slip.
                    set_zero_velocity(s, ip);
                    continue;
                }

                const Real un =
                    s.unkno(1, ip) * nx +
                    s.unkno(2, ip) * ny +
                    s.unkno(3, ip) * nz;

                // At a pressure outlet, positive normal velocity leaves the
                // domain. Negative normal velocity represents backflow through
                // the outlet. Remove only that inward normal component.
                if (un < 0.0)
                {
                    s.unkno(1, ip) -= un * nx;
                    s.unkno(2, ip) -= un * ny;
                    s.unkno(3, ip) -= un * nz;
                }
            }
        }
    }
}
