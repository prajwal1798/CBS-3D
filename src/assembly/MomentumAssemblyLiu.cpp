//=============================================================================
// CBS3D++_SI
//
// CBS Step-1 momentum assembly following Chun-Bin Liu (2005), Chapter 4 and
// Appendix B for the three-dimensional P1/TET4 RANS CBS formulation.
//
// Source equations used here:
//
//   Eq. (4.20)  characteristic Step-1 intermediate momentum equation
//   Eq. (B.3)   conservative momentum convection matrix C_u
//   Eq. (B.5)   deviatoric molecular+turbulent diffusion matrix K_tau
//   Eq. (B.9)   viscous/turbulent traction vector
//   Eq. (B.13)  CBS stabilisation matrix K_u
//
// Two details are intentionally different from the historical C++ port:
//
//   1. convection is assembled as the strong conservative operator
//        integral N_a div(u u_i) dV,
//      rather than by linearly interpolating the nodal products u_j u_i;
//
//   2. diffusion uses the full deviatoric RANS stress
//        nu_eff [grad(u)+grad(u)^T-(2/3)div(u)I],
//      rather than three independent scalar Laplacians.
//
// The CBS characteristic term is the Appendix-B matrix
//
//   -dt/2 integral div(u N_a) div(u N_b) dV,
//
// applied to each momentum component.  No separate advective/characteristic
// boundary correction is added: Appendix B represents those terms by the
// volume matrices C_u and K_u, while the explicit boundary vector in this part
// of the formulation is the stress traction f_u.
//=============================================================================

#include "cbs/assembly/MomentumAssembly.hpp"
#include "cbs/assembly/LiuNithiarasuKernels.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace cbs
{
    namespace
    {
        Int dNkdx_index(
            const CBSStateSI& s,
            const Int ie,
            const Int dim,
            const Int local_node)
        {
            return (ie - 1) * s.cfg.ndim * s.cfg.nep
                + (dim - 1) * s.cfg.nep
                + local_node;
        }

        Real grad(
            const CBSStateSI& s,
            const Int ie,
            const Int dim,
            const Int local_node)
        {
            return s.dNkdx(dNkdx_index(s, ie, dim, local_node));
        }

        bool is_fluid_element(const CBSStateSI& s, const Int ie)
        {
            return s.mat_elem(ie) == 0;
        }

        void validate_step1_dimensions(const CBSStateSI& s)
        {
            if (s.cfg.ndim != 3 || s.cfg.nep != 4 ||
                s.cfg.nsid != 4 || s.cfg.nsidp != 3 || s.cfg.gdim != 13)
            {
                throw std::runtime_error(
                    "MomentumAssembly - Liu/Nithiarasu Step 1 requires "
                    "ndim=3, nep=4, nsid=4, nsidp=3, gdim=13");
            }
        }

        Real momentum_diffusivity(const CBSStateSI& s, const Int ie)
        {
            if (s.cfg.dimensional_mode > 0 &&
                s.cfg.material_properties_enabled > 0)
            {
                const Real rho = s.rho_e(ie);
                if (!(rho > 0.0) || !std::isfinite(rho))
                {
                    throw std::runtime_error(
                        "MomentumAssembly - invalid fluid density at element "
                        + std::to_string(ie));
                }

                Real mu = s.mu_e(ie);
                if (s.cfg.turbulence_on > 0)
                {
                    mu = s.mu_eff_e(ie);
                }

                if (mu < 0.0 || !std::isfinite(mu))
                {
                    throw std::runtime_error(
                        "MomentumAssembly - invalid effective viscosity at element "
                        + std::to_string(ie));
                }

                return mu / rho;
            }

            if (!(s.cfg.ani >= 0.0) || !std::isfinite(s.cfg.ani))
            {
                throw std::runtime_error(
                    "MomentumAssembly - invalid non-dimensional viscosity");
            }
            return s.cfg.ani;
        }

        void velocity_gradient_and_divergence(
            const CBSStateSI& s,
            const Int ie,
            const Real velocity[5][4],
            Real velocity_gradient[4][4],
            Real& div_u)
        {
            for (Int i = 1; i <= 3; ++i)
            {
                for (Int a = 1; a <= 4; ++a)
                {
                    for (Int j = 1; j <= 3; ++j)
                    {
                        velocity_gradient[i][j] +=
                            velocity[a][i] * grad(s, ie, j, a);
                    }
                }
            }

            div_u =
                velocity_gradient[1][1]
                + velocity_gradient[2][2]
                + velocity_gradient[3][3];
        }

        // Liu Appendix B Eq. (B.3):
        //
        //   C_u = integral N_u^T div(u N_u) dV.
        //
        // For each momentum component q=u_i,
        //
        //   div(u q) = u.grad(q) + q div(u).
        //
        // With P1 u and q this is a P1 scalar, so the V/20 consistent-load
        // identity integrates N_a div(u q) exactly.
        void step1_conservative_convection(
            const CBSStateSI& s,
            const Int ie,
            const Real velocity[5][4],
            const Real velocity_gradient[4][4],
            const Real div_u,
            Real lrhs[4][5])
        {
            const Real volume = s.detJ(ie) / 6.0;

            for (Int i = 1; i <= 3; ++i)
            {
                Real q[5] = {};
                Real grad_q[4] = {};
                Real div_uq[5] = {};

                for (Int a = 1; a <= 4; ++a)
                {
                    q[a] = velocity[a][i];
                }
                for (Int j = 1; j <= 3; ++j)
                {
                    grad_q[j] = velocity_gradient[i][j];
                }

                liu_nithiarasu::conservative_scalar_divergence_nodes(
                    velocity,
                    q,
                    grad_q,
                    div_u,
                    div_uq);

                for (Int a = 1; a <= 4; ++a)
                {
                    lrhs[i][a] -=
                        liu_nithiarasu::consistent_linear_load(
                            volume,
                            div_uq,
                            a);
                }
            }
        }

        // Liu Appendix B Eq. (B.13):
        //
        //   K_u = -1/2 integral div(u N)^T div(u N) dV.
        //
        // Multiplication by the old nodal momentum gives, component by
        // component,
        //
        //   -dt/2 integral div(u N_a) div(u u_i) dV.
        //
        // Both factors are P1 scalars, so the product integral is exact.
        void step1_characteristic_stabilisation(
            const CBSStateSI& s,
            const Int ie,
            const Real velocity[5][4],
            const Real velocity_gradient[4][4],
            const Real div_u,
            Real lrhs[4][5])
        {
            const Real dt = s.delte(ie);
            if (!(dt > 0.0) || !std::isfinite(dt))
            {
                throw std::runtime_error(
                    "MomentumAssembly - invalid element timestep at element "
                    + std::to_string(ie));
            }

            const Real volume = s.detJ(ie) / 6.0;
            Real div_momentum[4][5] = {};

            for (Int i = 1; i <= 3; ++i)
            {
                Real q[5] = {};
                Real grad_q[4] = {};

                for (Int b = 1; b <= 4; ++b)
                {
                    q[b] = velocity[b][i];
                }
                for (Int j = 1; j <= 3; ++j)
                {
                    grad_q[j] = velocity_gradient[i][j];
                }

                liu_nithiarasu::conservative_scalar_divergence_nodes(
                    velocity,
                    q,
                    grad_q,
                    div_u,
                    div_momentum[i]);
            }

            for (Int a = 1; a <= 4; ++a)
            {
                Real grad_na[4] = {};
                Real div_u_na[5] = {};
                for (Int j = 1; j <= 3; ++j)
                {
                    grad_na[j] = grad(s, ie, j, a);
                }

                liu_nithiarasu::test_function_divergence_nodes(
                    velocity,
                    grad_na,
                    div_u,
                    a,
                    div_u_na);

                for (Int i = 1; i <= 3; ++i)
                {
                    lrhs[i][a] -= 0.5 * dt *
                        liu_nithiarasu::p1_product_integral(
                            volume,
                            div_u_na,
                            div_momentum[i]);
                }
            }
        }

        // Liu Appendix B Eqs. (B.5), (B.9)-(B.12).
        //
        //   tau_ij = nu_eff [du_i/dx_j + du_j/dx_i
        //                    - (2/3) div(u) delta_ij].
        //
        // Since a P1/TET4 velocity has a constant element gradient, tau is
        // constant within the element when nu_eff is represented elementwise.
        void step1_deviatoric_diffusion(
            const CBSStateSI& s,
            const Int ie,
            const Real velocity_gradient[4][4],
            const Real nu_eff,
            Real lrhs[4][5])
        {
            Real tau[4][4] = {};
            liu_nithiarasu::deviatoric_stress(
                velocity_gradient,
                nu_eff,
                tau);

            const Real volume = s.detJ(ie) / 6.0;

            // Volume weak term: -integral grad(N_a)_j tau_ij dV.
            for (Int a = 1; a <= 4; ++a)
            {
                for (Int i = 1; i <= 3; ++i)
                {
                    Real value = 0.0;
                    for (Int j = 1; j <= 3; ++j)
                    {
                        value += grad(s, ie, j, a) * tau[i][j];
                    }
                    lrhs[i][a] -= volume * value;
                }
            }

            // Natural traction vector f_u.  annxf is area-weighted outward
            // normal and fdif[2]=1/3 integrates a constant traction against a
            // TRI3 face shape function.
            for (Int is = 1; is <= 4; ++is)
            {
                if (s.fedge(is, ie) == 0)
                {
                    continue;
                }

                Real traction_area[4] = {};
                for (Int i = 1; i <= 3; ++i)
                {
                    for (Int j = 1; j <= 3; ++j)
                    {
                        traction_area[i] +=
                            tau[i][j] * s.annxf(j, is, ie);
                    }
                }

                for (Int face_node = 1; face_node <= 3; ++face_node)
                {
                    const Int a = s.ippn1(is, face_node);
                    for (Int i = 1; i <= 3; ++i)
                    {
                        lrhs[i][a] +=
                            traction_area[i] * s.cfg.fdif[2];
                    }
                }
            }
        }
    }

    void MomentumAssembly::assembleStep1Rhs(CBSStateSI& s)
    {
        validate_step1_dimensions(s);
        s.rhs.fill(0.0);

        if (!(s.cfg.dtreal > 0.0) || !std::isfinite(s.cfg.dtreal))
        {
            throw std::runtime_error(
                "MomentumAssembly - dtreal must be positive before Step 1");
        }

        int bad_detj = 0;

        for (Int colour = 0; colour < s.ncolor; ++colour)
        {
            const Int cbeg = s.color_ptr[static_cast<Size>(colour)];
            const Int cend = s.color_ptr[static_cast<Size>(colour) + 1];

#ifdef CBS3D_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (Int k = cbeg; k < cend; ++k)
            {
                const Int ie = s.color_elem[static_cast<Size>(k)];
                if (!is_fluid_element(s, ie))
                {
                    continue;
                }

                if (!(s.detJ(ie) > 0.0) || !std::isfinite(s.detJ(ie)))
                {
#ifdef CBS3D_USE_OPENMP
#pragma omp atomic write
#endif
                    bad_detj = 1;
                    continue;
                }

                Real velocity[5][4] = {};
                Real velocity_gradient[4][4] = {};
                Real lrhs[4][5] = {};
                Int ipn[5] = {};

                for (Int a = 1; a <= 4; ++a)
                {
                    const Int ip = s.intma(a, ie);
                    ipn[a] = ip;
                    for (Int i = 1; i <= 3; ++i)
                    {
                        velocity[a][i] = s.unkn1(i, ip);
                    }
                }

                Real div_u = 0.0;
                velocity_gradient_and_divergence(
                    s,
                    ie,
                    velocity,
                    velocity_gradient,
                    div_u);

                step1_conservative_convection(
                    s,
                    ie,
                    velocity,
                    velocity_gradient,
                    div_u,
                    lrhs);

                step1_characteristic_stabilisation(
                    s,
                    ie,
                    velocity,
                    velocity_gradient,
                    div_u,
                    lrhs);

                step1_deviatoric_diffusion(
                    s,
                    ie,
                    velocity_gradient,
                    momentum_diffusivity(s, ie),
                    lrhs);

                for (Int a = 1; a <= 4; ++a)
                {
                    const Int ip = ipn[a];
                    for (Int i = 1; i <= 3; ++i)
                    {
                        s.rhs(i, ip) += lrhs[i][a];
                    }
                }
            }
        }

        if (bad_detj != 0)
        {
            throw std::runtime_error(
                "MomentumAssembly - invalid detJ at one or more fluid elements");
        }
    }

    void MomentumAssembly::applyRealTimeMomentumTerm(CBSStateSI& s)
    {
        (void)s;
    }
}
