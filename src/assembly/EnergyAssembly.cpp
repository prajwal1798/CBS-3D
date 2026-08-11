//=============================================================================
// CBS3D++_SI
//
// Energy-equation residual assembly for three-dimensional conjugate heat
// transfer.
//
// Governing equation:
//
//     rho cp [ dT/dt + u . grad(T) ]
//         = div(k grad(T)) + Q
//
// The assembled thermal residual is:
//
//     r_T = r_conv + r_stab + r_diff + r_source + r_flux
//=============================================================================

#include "cbs/assembly/EnergyAssembly.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace cbs
{
    namespace
    {
        Int dNkdx_index(
            const CBSStateSI& s,
            Int ie,
            Int dim,
            Int local_node)
        {
            return (ie - 1) * s.cfg.ndim * s.cfg.nep
                + (dim - 1) * s.cfg.nep
                + local_node;
        }

        bool is_fluid_element(
            const CBSStateSI& s,
            Int ie)
        {
            return s.mat_elem(ie) == 0;
        }

        void validate_energy_dimensions(const CBSStateSI& s)
        {
            if (s.cfg.ndim != 3 ||
                s.cfg.nep != 4 ||
                s.cfg.nsid != 4 ||
                s.cfg.nsidp != 3 ||
                s.cfg.ndim1 != 4)
            {
                throw std::runtime_error(
                    "EnergyAssembly - CBS3D energy assembly requires ndim=3, nep=4, nsid=4, nsidp=3, ndim1=4");
            }
        }

        Real grad(
            const CBSStateSI& s,
            Int ie,
            Int dim,
            Int local_node)
        {
            return s.dNkdx(dNkdx_index(s, ie, dim, local_node));
        }

        void compute_temperature_gradient(
            const CBSStateSI& s,
            Int ie,
            Real& dTdx,
            Real& dTdy,
            Real& dTdz)
        {
            dTdx = 0.0;
            dTdy = 0.0;
            dTdz = 0.0;

            for (Int a = 1; a <= s.cfg.nep; ++a)
            {
                const Int ip = s.intma(a, ie);
                const Real T = s.temperature1(ip);

                dTdx += grad(s, ie, 1, a) * T;
                dTdy += grad(s, ie, 2, a) * T;
                dTdz += grad(s, ie, 3, a) * T;
            }
        }

        void validate_material_properties(
            const CBSStateSI& s,
            Int ie)
        {
            if (s.rho_cp_e(ie) <= 0.0 || !std::isfinite(s.rho_cp_e(ie)))
            {
                throw std::runtime_error(
                    "EnergyAssembly::assembleStep4Rhs - invalid rho*Cp at element "
                    + std::to_string(ie));
            }

            if (s.k_e(ie) <= 0.0 || !std::isfinite(s.k_e(ie)))
            {
                throw std::runtime_error(
                    "EnergyAssembly::assembleStep4Rhs - invalid thermal conductivity at element "
                    + std::to_string(ie));
            }

            if (!std::isfinite(s.Qvol_e(ie)))
            {
                throw std::runtime_error(
                    "EnergyAssembly::assembleStep4Rhs - invalid volumetric source at element "
                    + std::to_string(ie));
            }
        }

        void add_fluid_convection(
            const CBSStateSI& s,
            Int ie,
            Real dTdx,
            Real dTdy,
            Real dTdz,
            Real lrhs[5])
        {
            const Real adv_factor = s.rho_cp_e(ie) * s.detJ(ie) / 120.0;

            Real u_sum = 0.0;
            Real v_sum = 0.0;
            Real w_sum = 0.0;

            for (Int b = 1; b <= s.cfg.nep; ++b)
            {
                const Int ip = s.intma(b, ie);
                u_sum += s.unkno(1, ip);
                v_sum += s.unkno(2, ip);
                w_sum += s.unkno(3, ip);
            }

            for (Int a = 1; a <= s.cfg.nep; ++a)
            {
                const Int ip = s.intma(a, ie);

                const Real u_weight = u_sum + s.unkno(1, ip);
                const Real v_weight = v_sum + s.unkno(2, ip);
                const Real w_weight = w_sum + s.unkno(3, ip);

                const Real advT =
                    u_weight * dTdx +
                    v_weight * dTdy +
                    w_weight * dTdz;

                lrhs[a] -= adv_factor * advT;
            }
        }

        void add_fluid_convection_stabilisation(
            const CBSStateSI& s,
            Int ie,
            Real dTdx,
            Real dTdy,
            Real dTdz,
            Real lrhs[5])
        {
            const Real dt = s.delte(ie);

            if (dt <= 0.0 || !std::isfinite(dt))
            {
                throw std::runtime_error(
                    "EnergyAssembly::assembleStep4Rhs - element timestep must be positive for thermal stabilisation at element "
                    + std::to_string(ie));
            }

            Real ubar = 0.0;
            Real vbar = 0.0;
            Real wbar = 0.0;

            for (Int a = 1; a <= s.cfg.nep; ++a)
            {
                const Int ip = s.intma(a, ie);
                ubar += s.unkno(1, ip);
                vbar += s.unkno(2, ip);
                wbar += s.unkno(3, ip);
            }

            ubar /= static_cast<Real>(s.cfg.nep);
            vbar /= static_cast<Real>(s.cfg.nep);
            wbar /= static_cast<Real>(s.cfg.nep);

            const Real advT =
                ubar * dTdx +
                vbar * dTdy +
                wbar * dTdz;

            const Real volume = s.detJ(ie) / 6.0;
            const Real tau_factor = 0.5 * dt * s.rho_cp_e(ie) * volume;

            for (Int a = 1; a <= s.cfg.nep; ++a)
            {
                const Real u_gradNa =
                    ubar * grad(s, ie, 1, a) +
                    vbar * grad(s, ie, 2, a) +
                    wbar * grad(s, ie, 3, a);

                lrhs[a] += tau_factor * u_gradNa * advT;
            }
        }

        void add_thermal_diffusion(
            const CBSStateSI& s,
            Int ie,
            Real dTdx,
            Real dTdy,
            Real dTdz,
            Real lrhs[5])
        {
            const Real volume = s.detJ(ie) / 6.0;
            Real k = s.k_e(ie);

            if (s.cfg.turbulence_on > 0 &&
                s.cfg.turbulent_thermal_diffusivity_on > 0 &&
                is_fluid_element(s, ie))
            {
                k = s.k_eff_e(ie);
            }

            if (k <= 0.0 || !std::isfinite(k))
            {
                throw std::runtime_error(
                    "EnergyAssembly::assembleStep4Rhs - invalid effective thermal conductivity at element "
                    + std::to_string(ie));
            }

            for (Int a = 1; a <= s.cfg.nep; ++a)
            {
                const Real gradNi_dot_gradT =
                    grad(s, ie, 1, a) * dTdx +
                    grad(s, ie, 2, a) * dTdy +
                    grad(s, ie, 3, a) * dTdz;

                lrhs[a] -= k * volume * gradNi_dot_gradT;
            }
        }

        // Returns the effective volumetric source for one element.
        //
        // Fluid (material ID 0): source comes only from .matprop Qvol.
        // Solid (material ID != 0): non-zero .par source_solid is applied to
        // every solid element. If source_solid is zero, .matprop Qvol remains
        // available for backward compatibility/material-specific heating.
        // Supplying both for the same solid is rejected to prevent double load.
        Real element_volumetric_source(
            const CBSStateSI& s,
            Int ie)
        {
            const Real material_source = s.Qvol_e(ie);

            if (!std::isfinite(material_source))
            {
                throw std::runtime_error(
                    "EnergyAssembly::assembleStep4Rhs - non-finite .matprop volumetric source at element "
                    + std::to_string(ie));
            }

            if (is_fluid_element(s, ie))
            {
                return material_source;
            }

            const Real parameter_source = s.cfg.source_solid;

            if (!std::isfinite(parameter_source))
            {
                throw std::runtime_error(
                    "EnergyAssembly::assembleStep4Rhs - non-finite source_solid in .par");
            }

            if (parameter_source != 0.0 && material_source != 0.0)
            {
                throw std::runtime_error(
                    "EnergyAssembly::assembleStep4Rhs - solid volumetric heat source is specified in both .par source_solid and .matprop Qvol; use exactly one source definition");
            }

            return parameter_source != 0.0
                ? parameter_source
                : material_source;
        }

        // P1 tetrahedron source:
        //     int_V N_a Q dV = Q V/4 = Q det(J)/24.
        void add_volumetric_source(
            const CBSStateSI& s,
            Int ie,
            Real lrhs[5])
        {
            const Real qvol = element_volumetric_source(s, ie);

            if (qvol == 0.0)
            {
                return;
            }

            const Real source =
                qvol * s.detJ(ie) * s.cfg.mass_factor;

            for (Int a = 1; a <= s.cfg.nep; ++a)
            {
                lrhs[a] += source;
            }
        }

        void add_prescribed_heat_flux(CBSStateSI& s)
        {
            if (s.cfg.heat_flux_bc == 0.0)
            {
                return;
            }

            for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
            {
                const Int bc = s.iside(s.cfg.bsid, ib);

                if (bc != s.cfg.bc_noslip_heatflux_wall)
                {
                    continue;
                }

                const Real area = s.face_norm(4, ib);

                if (area <= 0.0 || !std::isfinite(area))
                {
                    throw std::runtime_error(
                        "EnergyAssembly::assembleStep4Rhs - invalid heat-flux boundary face area");
                }

                const Real contribution =
                    s.cfg.heat_flux_bc * area /
                    static_cast<Real>(s.cfg.nsidp);

                for (Int in = 1; in <= s.cfg.nsidp; ++in)
                {
                    const Int ip = s.iside(in, ib);

                    if (ip < 1 || ip > s.cfg.npoin)
                    {
                        throw std::runtime_error(
                            "EnergyAssembly::assembleStep4Rhs - heat-flux boundary node out of range");
                    }

                    s.rhs1(ip) += contribution;
                }
            }
        }
    }

    void EnergyAssembly::assembleStep4Rhs(CBSStateSI& s)
    {
        validate_energy_dimensions(s);
        s.rhs1.fill(0.0);

        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            if (s.detJ(ie) <= 0.0 || !std::isfinite(s.detJ(ie)))
            {
                throw std::runtime_error(
                    "EnergyAssembly::assembleStep4Rhs - invalid detJ at element "
                    + std::to_string(ie));
            }

            validate_material_properties(s, ie);

            Real dTdx = 0.0;
            Real dTdy = 0.0;
            Real dTdz = 0.0;
            compute_temperature_gradient(s, ie, dTdx, dTdy, dTdz);

            Real lrhs[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };

            if (is_fluid_element(s, ie))
            {
                add_fluid_convection(s, ie, dTdx, dTdy, dTdz, lrhs);
                add_fluid_convection_stabilisation(s, ie, dTdx, dTdy, dTdz, lrhs);
            }

            add_thermal_diffusion(s, ie, dTdx, dTdy, dTdz, lrhs);
            add_volumetric_source(s, ie, lrhs);

            for (Int a = 1; a <= s.cfg.nep; ++a)
            {
                const Int ip = s.intma(a, ie);
                s.rhs1(ip) += lrhs[a];
            }
        }

        add_prescribed_heat_flux(s);
    }

    void EnergyAssembly::applyRealTimeEnergyTerm(CBSStateSI& s)
    {
        (void)s;
    }
}
