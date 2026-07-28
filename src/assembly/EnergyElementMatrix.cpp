//=============================================================================
// CBS3D++_SI
//
// Element transport matrix used by the bounded thermal AFC update.
//
// The matrix is defined so that the temperature-dependent element residual is:
//
//     r_e(T) = -A_e T_e.
//
// Volumetric heat generation and prescribed boundary heat flux are independent
// loads and are not included here. The coefficients reproduce the convection,
// characteristic and diffusion terms assembled in EnergyAssembly.cpp.
//=============================================================================

#include "cbs/assembly/EnergyAssembly.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace cbs
{
    namespace
    {
        Int gradient_index(
            const CBSStateSI& s,
            const Int element,
            const Int dimension,
            const Int local_node)
        {
            return (element - 1) * s.cfg.ndim * s.cfg.nep
                 + (dimension - 1) * s.cfg.nep
                 + local_node;
        }

        Real gradient(
            const CBSStateSI& s,
            const Int element,
            const Int dimension,
            const Int local_node)
        {
            return s.dNkdx(
                gradient_index(s, element, dimension, local_node));
        }

        bool fluid_element(const CBSStateSI& s, const Int element)
        {
            return s.mat_elem(element) == 0;
        }

        Real effective_conductivity(
            const CBSStateSI& s,
            const Int element)
        {
            Real conductivity = s.k_e(element);

            if (s.cfg.turbulence_on > 0
                && s.cfg.turbulent_thermal_diffusivity_on > 0
                && fluid_element(s, element))
            {
                conductivity = s.k_eff_e(element);
            }

            return conductivity;
        }
    }

    void EnergyAssembly::buildElementTransportMatrix(
        const CBSStateSI& s,
        const Int ie,
        Real matrix[5][5])
    {
        if (s.cfg.ndim != 3 || s.cfg.nep != 4)
        {
            throw std::runtime_error(
                "EnergyAssembly::buildElementTransportMatrix requires a "
                "four-node three-dimensional tetrahedron");
        }

        if (ie < 1 || ie > s.cfg.nelem)
        {
            throw std::runtime_error(
                "EnergyAssembly::buildElementTransportMatrix - element index "
                "out of range");
        }

        const Real determinant = s.detJ(ie);
        const Real rho_cp = s.rho_cp_e(ie);
        const Real conductivity = effective_conductivity(s, ie);

        if (determinant <= 0.0 || !std::isfinite(determinant))
        {
            throw std::runtime_error(
                "EnergyAssembly::buildElementTransportMatrix - invalid detJ at "
                "element " + std::to_string(ie));
        }

        if (rho_cp <= 0.0 || !std::isfinite(rho_cp))
        {
            throw std::runtime_error(
                "EnergyAssembly::buildElementTransportMatrix - invalid rho*Cp "
                "at element " + std::to_string(ie));
        }

        if (conductivity <= 0.0 || !std::isfinite(conductivity))
        {
            throw std::runtime_error(
                "EnergyAssembly::buildElementTransportMatrix - invalid thermal "
                "conductivity at element " + std::to_string(ie));
        }

        for (Int a = 0; a <= 4; ++a)
        {
            for (Int b = 0; b <= 4; ++b)
            {
                matrix[a][b] = 0.0;
            }
        }

        const Real volume = determinant / 6.0;

        // Diffusion:
        //
        //     r_a^diff = -sum_b k V grad(N_a).grad(N_b) T_b.
        for (Int a = 1; a <= 4; ++a)
        {
            for (Int b = 1; b <= 4; ++b)
            {
                const Real gradient_dot =
                    gradient(s, ie, 1, a) * gradient(s, ie, 1, b)
                  + gradient(s, ie, 2, a) * gradient(s, ie, 2, b)
                  + gradient(s, ie, 3, a) * gradient(s, ie, 3, b);

                matrix[a][b] += conductivity * volume * gradient_dot;
            }
        }

        if (!fluid_element(s, ie))
        {
            return;
        }

        Real velocity_sum[4] = {0.0, 0.0, 0.0, 0.0};

        for (Int local_node = 1; local_node <= 4; ++local_node)
        {
            const Int node = s.intma(local_node, ie);
            velocity_sum[1] += s.unkno(1, node);
            velocity_sum[2] += s.unkno(2, node);
            velocity_sum[3] += s.unkno(3, node);
        }

        // Exact P1 Galerkin convection:
        //
        //     A_ab^conv = rhoCp detJ/120
        //                   (sum_c u_c + u_a).grad(N_b).
        const Real convection_factor = rho_cp * determinant / 120.0;

        for (Int a = 1; a <= 4; ++a)
        {
            const Int node_a = s.intma(a, ie);

            const Real weighted_u =
                velocity_sum[1] + s.unkno(1, node_a);
            const Real weighted_v =
                velocity_sum[2] + s.unkno(2, node_a);
            const Real weighted_w =
                velocity_sum[3] + s.unkno(3, node_a);

            for (Int b = 1; b <= 4; ++b)
            {
                matrix[a][b] += convection_factor *
                    (weighted_u * gradient(s, ie, 1, b)
                   + weighted_v * gradient(s, ie, 2, b)
                   + weighted_w * gradient(s, ie, 3, b));
            }
        }

        const Real dt = s.delte(ie);
        if (dt <= 0.0 || !std::isfinite(dt))
        {
            throw std::runtime_error(
                "EnergyAssembly::buildElementTransportMatrix - invalid element "
                "timestep at element " + std::to_string(ie));
        }

        const Real ubar = velocity_sum[1] / 4.0;
        const Real vbar = velocity_sum[2] / 4.0;
        const Real wbar = velocity_sum[3] / 4.0;
        const Real characteristic_factor = 0.5 * dt * rho_cp * volume;

        // Existing CBS characteristic contribution:
        //
        //     r_a^char = +tau rhoCp V (u.gradN_a)(u.gradT)
        //              = -sum_b A_ab^char T_b,
        //
        // hence A_ab^char carries the negative sign below.
        for (Int a = 1; a <= 4; ++a)
        {
            const Real u_grad_a =
                ubar * gradient(s, ie, 1, a)
              + vbar * gradient(s, ie, 2, a)
              + wbar * gradient(s, ie, 3, a);

            for (Int b = 1; b <= 4; ++b)
            {
                const Real u_grad_b =
                    ubar * gradient(s, ie, 1, b)
                  + vbar * gradient(s, ie, 2, b)
                  + wbar * gradient(s, ie, 3, b);

                matrix[a][b] -=
                    characteristic_factor * u_grad_a * u_grad_b;
            }
        }
    }
}
