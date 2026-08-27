//=============================================================================
// CBS3D++_SI
//
// CBS Step-1 momentum assembly aligned with the 3-D RANS matrices in
// Chun-Bin Liu (2005), Appendix B.
//
// The crucial distinction from the historical C++ port is the diffusion/stress
// operator.  Liu Eq. B.5 uses
//
//   B^T (I_o - 2/3 m m^T) B
//
// multiplied by the molecular+turbulent viscosity.  In tensor notation the
// kinematic stress is
//
//   tau_ij = nu_eff (du_i/dx_j + du_j/dx_i
//                    - 2/3 div(u) delta_ij).
//
// This file therefore does NOT apply a scalar Laplacian independently to u, v
// and w.  The cross derivatives and deviatoric trace subtraction are retained.
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

        // Galerkin conservative momentum convection.  This is retained from the
        // established CBS implementation.
        void step1_convective_galerkin(
            const CBSStateSI& s,
            const Int ie,
            const Real lunk[4][4][5],
            Real lrhs[4][5])
        {
            Real lunksum[4][4] = {};
            Real lunksum_face[4][4] = {};
            const Real volume_factor = s.detJ(ie) * s.cfg.fcon[1];

            for (Int i = 1; i <= 3; ++i)
            {
                for (Int j = 1; j <= 3; ++j)
                {
                    for (Int a = 1; a <= 4; ++a)
                    {
                        lunksum[i][j] += lunk[i][j][a];
                    }
                }
            }

            for (Int a = 1; a <= 4; ++a)
            {
                for (Int i = 1; i <= 3; ++i)
                {
                    for (Int j = 1; j <= 3; ++j)
                    {
                        lrhs[i][a] +=
                            grad(s, ie, j, a) * volume_factor * lunksum[i][j];
                    }
                }
            }

            for (Int is = 1; is <= 4; ++is)
            {
                if (s.fedge(is, ie) == 0)
                {
                    continue;
                }

                for (Int i = 1; i <= 3; ++i)
                {
                    for (Int j = 1; j <= 3; ++j)
                    {
                        lunksum_face[i][j] =
                            lunksum[i][j] - lunk[i][j][is];
                    }
                }

                for (Int j = 1; j <= 3; ++j)
                {
                    const Real face_factor =
                        s.annxf(j, is, ie) * s.cfg.fcon[2];

                    for (Int face_node = 1; face_node <= 3; ++face_node)
                    {
                        const Int a = s.ippn1(is, face_node);
                        for (Int i = 1; i <= 3; ++i)
                        {
                            lrhs[i][a] -=
                                (lunksum_face[i][j] + lunk[i][j][a])
                                * face_factor;
                        }
                    }
                }
            }
        }

        // Existing transferred-CBS momentum stabilisation.  The RANS correction
        // in this source concerns the stress operator; this characteristic term
        // is kept identical to the established laminar CBS path so the change is
        // isolated and testable.
        void step1_characteristic_correction(
            const CBSStateSI& s,
            const Int ie,
            const Real lunk[4][4][5],
            const Real lunkno[4][5],
            Real lrhs[4][5])
        {
            Real lunksum[4][4] = {};
            Real lunksumk[4][4] = {};
            Real umean[4] = {};

            for (Int i = 1; i <= 3; ++i)
            {
                for (Int j = 1; j <= 3; ++j)
                {
                    for (Int a = 1; a <= 4; ++a)
                    {
                        lunksum[i][j] +=
                            grad(s, ie, j, a) * lunk[i][j][a];
                    }
                }
            }

            for (Int i = 1; i <= 3; ++i)
            {
                for (Int a = 1; a <= 4; ++a)
                {
                    umean[i] += lunkno[i][a];
                }
                umean[i] *= 0.25;
            }

            for (Int i = 1; i <= 3; ++i)
            {
                Real div_flux = 0.0;
                for (Int j = 1; j <= 3; ++j)
                {
                    div_flux += lunksum[i][j];
                }
                for (Int k = 1; k <= 3; ++k)
                {
                    lunksumk[i][k] = umean[k] * div_flux;
                }
            }

            const Real dt = s.delte(ie);
            if (!(dt > 0.0) || !std::isfinite(dt))
            {
                throw std::runtime_error(
                    "MomentumAssembly - invalid element timestep at element "
                    + std::to_string(ie));
            }

            const Real volume_factor =
                0.5 * dt * s.detJ(ie) * s.cfg.fcon[3];
            const Real face_factor_base =
                0.5 * dt * s.cfg.fcon[4];

            for (Int a = 1; a <= 4; ++a)
            {
                for (Int k = 1; k <= 3; ++k)
                {
                    for (Int i = 1; i <= 3; ++i)
                    {
                        lrhs[i][a] -=
                            grad(s, ie, k, a)
                            * lunksumk[i][k]
                            * volume_factor;
                    }
                }
            }

            for (Int is = 1; is <= 4; ++is)
            {
                if (s.fedge(is, ie) == 0)
                {
                    continue;
                }

                for (Int k = 1; k <= 3; ++k)
                {
                    const Real face_factor =
                        s.annxf(k, is, ie) * face_factor_base;
                    for (Int face_node = 1; face_node <= 3; ++face_node)
                    {
                        const Int a = s.ippn1(is, face_node);
                        for (Int i = 1; i <= 3; ++i)
                        {
                            lrhs[i][a] += lunksumk[i][k] * face_factor;
                        }
                    }
                }
            }
        }

        // Liu Appendix B Eq. B.5.  This is the formulation-level correction.
        void step1_deviatoric_diffusion(
            const CBSStateSI& s,
            const Int ie,
            const Real lunkno[4][5],
            const Real nu_eff,
            Real lrhs[4][5])
        {
            Real velocity_gradient[4][4] = {};
            for (Int i = 1; i <= 3; ++i)
            {
                for (Int a = 1; a <= 4; ++a)
                {
                    for (Int j = 1; j <= 3; ++j)
                    {
                        velocity_gradient[i][j] +=
                            lunkno[i][a] * grad(s, ie, j, a);
                    }
                }
            }

            Real tau[4][4] = {};
            liu_nithiarasu::deviatoric_stress(
                velocity_gradient,
                nu_eff,
                tau);

            const Real volume = s.detJ(ie) / 6.0;

            // Weak volume term: -int grad(N_a)_j tau_ij dV.
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

            // Natural traction term: +int N_a tau_ij n_j dGamma.
            // annxf stores the area-weighted outward normal and fdif[2] is the
            // established TRI3 integral factor (1/3).
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

#pragma omp parallel for schedule(static)
            for (Int k = cbeg; k < cend; ++k)
            {
                const Int ie = s.color_elem[static_cast<Size>(k)];
                if (!is_fluid_element(s, ie))
                {
                    continue;
                }

                if (!(s.detJ(ie) > 0.0) || !std::isfinite(s.detJ(ie)))
                {
#pragma omp atomic write
                    bad_detj = 1;
                    continue;
                }

                Real lunkno[4][5] = {};
                Real lunk[4][4][5] = {};
                Real lrhs[4][5] = {};
                Int ipn[5] = {};

                for (Int a = 1; a <= 4; ++a)
                {
                    const Int ip = s.intma(a, ie);
                    ipn[a] = ip;
                    for (Int i = 1; i <= 3; ++i)
                    {
                        lunkno[i][a] = s.unkn1(i, ip);
                    }
                }

                for (Int a = 1; a <= 4; ++a)
                {
                    for (Int j = 1; j <= 3; ++j)
                    {
                        for (Int i = 1; i <= 3; ++i)
                        {
                            lunk[i][j][a] =
                                lunkno[j][a] * lunkno[i][a];
                        }
                    }
                }

                step1_convective_galerkin(s, ie, lunk, lrhs);
                step1_characteristic_correction(s, ie, lunk, lunkno, lrhs);
                step1_deviatoric_diffusion(
                    s,
                    ie,
                    lunkno,
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
