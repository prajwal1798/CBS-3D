//=============================================================================
// CBS3D++_SI
//
// Momentum residual assembly for CBS Step 1.
//
// The momentum predictor has the nodal form:
//
//     M_L / dt * (u* - u^n) = r_m
//
// The residual assembled in this file is:
//
//     r_m = r_conv + r_char + r_diff
//
// where:
//
//     r_conv   Galerkin convective contribution
//     r_char   CBS characteristic correction
//     r_diff   viscous diffusion contribution
//
// The subsequent nodal update is performed in Steps::step1SemiImplicit():
//
//     u* = u^n + elcoe2 * rhs
//
// Only fluid elements contribute to the momentum equation.
//=============================================================================

#include "cbs/assembly/MomentumAssembly.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace cbs
{
    namespace
    {
        // Returns the one-dimensional storage position of:
        //
        //     dN_local_node / dx_dim
        //
        // for tetrahedral element ie.
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


        // Momentum is assembled only over fluid elements.
        bool is_fluid_element(
            const CBSStateSI& s,
            Int ie)
        {
            return s.mat_elem(ie) == 0;
        }


        // Checks the fixed dimensions required by the current
        // three-dimensional P1 tetrahedral CBS formulation.
        void validate_step1_dimensions(const CBSStateSI& s)
        {
            if (s.cfg.ndim != 3 ||
                s.cfg.nep != 4 ||
                s.cfg.nsid != 4 ||
                s.cfg.nsidp != 3 ||
                s.cfg.gdim != 13)
            {
                throw std::runtime_error(
                    "MomentumAssembly - CBS3D Step 1 requires ndim=3, nep=4, nsid=4, nsidp=3, gdim=13");
            }
        }


        // Returns one Cartesian derivative of a tetrahedral shape function:
        //
        //     grad(N_a)_dim = dN_a / dx_dim
        Real grad(
            const CBSStateSI& s,
            Int ie,
            Int dim,
            Int local_node)
        {
            return s.dNkdx(dNkdx_index(s, ie, dim, local_node));
        }


        // Returns the kinematic viscosity used in the Step-1 viscous term.
        //
        // Laminar dimensional mode:
        //
        //     nu_e = mu_e / rho_e
        //
        // Spalart-Allmaras dimensional mode:
        //
        //     nu_eff,e = mu_eff_e / rho_e
        //
        // where:
        //
        //     mu_eff_e = mu_e + mu_t_e
        //
        // Non-dimensional laminar mode continues to use ani.  The first SA
        // implementation is intended for dimensional material cases, where rho
        // and mu are available element by element.
        Real momentum_diffusivity(
            const CBSStateSI& s,
            Int ie)
        {
            if (s.cfg.dimensional_mode > 0 && s.cfg.material_properties_enabled > 0)
            {
                if (s.rho_e(ie) <= 0.0 || !std::isfinite(s.rho_e(ie)))
                {
                    throw std::runtime_error(
                        "MomentumAssembly::assembleStep1Rhs - invalid fluid density at element "
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
                        "MomentumAssembly::assembleStep1Rhs - invalid effective viscosity at element "
                        + std::to_string(ie));
                }

                return mu / s.rho_e(ie);
            }

            return s.cfg.ani;
        }


        //=====================================================================
        // Assembles the Galerkin convective contribution for one element.
        //
        // The strong convective term is:
        //
        //     div(u_j u_i)
        //
        // After integration by parts, the contribution to the predictor RHS is:
        //
        //     r_conv,i,a^(e)
        //       = integral(V_e) grad(N_a)_j (u_j u_i) dV
        //       - integral(Gamma_e) N_a (u_j u_i) n_j dGamma
        //
        // For a P1 tetrahedron, the implemented volume integration uses:
        //
        //     det(J_e) * fcon[1] = det(J_e) / 24 = V_e / 4
        //
        // Inputs:
        //     lunk[i][j][a] = u_j(a) u_i(a)
        //
        // Output:
        //     lrhs           local momentum residual
        //=====================================================================
        void step1_convective_galerkin(
            const CBSStateSI& s,
            Int ie,
            const Real lunk[4][4][5],
            Real lrhs[4][5])
        {
            Real lunksum[4][4] = {};
            Real lunksum_face[4][4] = {};

            const Real volume_factor = s.detJ(ie) * s.cfg.fcon[1];

            // Sum the nodal momentum-flux products:
            //
            //     lunksum(i,j) = sum_a u_j(a) u_i(a)
            for (Int i = 1; i <= s.cfg.ndim; ++i)
            {
                for (Int j = 1; j <= s.cfg.ndim; ++j)
                {
                    for (Int in = 1; in <= s.cfg.nep; ++in)
                    {
                        lunksum[i][j] += lunk[i][j][in];
                    }
                }
            }

            // Add the element-volume convective contribution.
            for (Int a = 1; a <= s.cfg.nep; ++a)
            {
                for (Int i = 1; i <= s.cfg.ndim; ++i)
                {
                    for (Int j = 1; j <= s.cfg.ndim; ++j)
                    {
                        lrhs[i][a] +=
                            grad(s, ie, j, a) * volume_factor * lunksum[i][j];
                    }
                }
            }

            // Add the convective boundary contribution on marked faces.
            //
            // Local face is is opposite local node is. Therefore:
            //
            //     lunksum_face = lunksum - lunk at the opposite node
            //
            // gives the sum over the three nodes lying on the triangular face.
            for (Int is = 1; is <= s.cfg.nsid; ++is)
            {
                if (s.fedge(is, ie) == 0)
                {
                    continue;
                }

                for (Int i = 1; i <= s.cfg.ndim; ++i)
                {
                    for (Int j = 1; j <= s.cfg.ndim; ++j)
                    {
                        lunksum_face[i][j] = lunksum[i][j] - lunk[i][j][is];
                    }
                }

                for (Int j = 1; j <= s.cfg.ndim; ++j)
                {
                    const Real face_factor = s.annxf(j, is, ie) * s.cfg.fcon[2];

                    for (Int face_node = 1; face_node <= s.cfg.nsidp; ++face_node)
                    {
                        const Int a = s.ippn1(is, face_node);

                        for (Int i = 1; i <= s.cfg.ndim; ++i)
                        {
                            lrhs[i][a] -=
                                (lunksum_face[i][j] + lunk[i][j][a]) * face_factor;
                        }
                    }
                }
            }
        }


        //=====================================================================
        // Assembles the CBS characteristic correction for one element.
        //
        // The correction has the form:
        //
        //     r_char,i,a^(e)
        //       = -dt_e/2 integral(V_e)
        //           dN_a/dx_k u_bar_k d(u_j u_i)/dx_j dV
        //
        //         +dt_e/2 integral(Gamma_e)
        //           N_a u_bar_k d(u_j u_i)/dx_j n_k dGamma
        //
        // where the element-average velocity is:
        //
        //     u_bar_k = (1 / 4) sum_a u_k(a)
        //
        // This term provides the characteristic stabilisation used by the CBS
        // momentum predictor.
        //=====================================================================
        void step1_characteristic_correction(
            const CBSStateSI& s,
            Int ie,
            const Real lunk[4][4][5],
            const Real lunkno[4][5],
            Real lrhs[4][5])
        {
            Real lunksum[4][4] = {};
            Real lunksumk[4][4] = {};
            Real umean[4] = {};

            // Calculate:
            //
            //     lunksum(i,j) =
            //         sum_a dN_a/dx_j [u_j(a) u_i(a)]
            for (Int i = 1; i <= s.cfg.ndim; ++i)
            {
                for (Int j = 1; j <= s.cfg.ndim; ++j)
                {
                    for (Int a = 1; a <= s.cfg.nep; ++a)
                    {
                        lunksum[i][j] += grad(s, ie, j, a) * lunk[i][j][a];
                    }
                }
            }

            // Calculate the element-average velocity.
            for (Int i = 1; i <= s.cfg.ndim; ++i)
            {
                for (Int a = 1; a <= s.cfg.nep; ++a)
                {
                    umean[i] += lunkno[i][a];
                }

                umean[i] /= static_cast<Real>(s.cfg.nep);
            }

            // Form the characteristic directional product used by the
            // transferred CBS formulation.
            for (Int i = 1; i <= s.cfg.ndim; ++i)
            {
                Real sum_j = 0.0;

                for (Int j = 1; j <= s.cfg.ndim; ++j)
                {
                    sum_j += lunksum[i][j];
                }

                for (Int k = 1; k <= s.cfg.ndim; ++k)
                {
                    lunksumk[i][k] = umean[k] * sum_j;
                }
            }

            const Real dt = s.delte(ie);
            if (dt <= 0.0 || !std::isfinite(dt))
            {
                throw std::runtime_error(
                    "MomentumAssembly::assembleStep1Rhs - invalid element timestep at element "
                    + std::to_string(ie));
            }

            const Real ldelte =
                0.5 * dt * s.detJ(ie) * s.cfg.fcon[3];

            const Real ldelte2 =
                0.5 * dt * s.cfg.fcon[4];

            // Add the characteristic volume contribution.
            for (Int a = 1; a <= s.cfg.nep; ++a)
            {
                for (Int k = 1; k <= s.cfg.ndim; ++k)
                {
                    for (Int i = 1; i <= s.cfg.ndim; ++i)
                    {
                        lrhs[i][a] -= grad(s, ie, k, a) * lunksumk[i][k] * ldelte;
                    }
                }
            }

            // Add the characteristic boundary contribution.
            for (Int is = 1; is <= s.cfg.nsid; ++is)
            {
                if (s.fedge(is, ie) == 0)
                {
                    continue;
                }

                for (Int k = 1; k <= s.cfg.ndim; ++k)
                {
                    const Real face_factor = s.annxf(k, is, ie) * ldelte2;

                    for (Int face_node = 1; face_node <= s.cfg.nsidp; ++face_node)
                    {
                        const Int a = s.ippn1(is, face_node);

                        for (Int i = 1; i <= s.cfg.ndim; ++i)
                        {
                            lrhs[i][a] += lunksumk[i][k] * face_factor;
                        }
                    }
                }
            }
        }


        //=====================================================================
        // Assembles the viscous diffusion contribution for one element.
        //
        // The velocity gradient is:
        //
        //     du_i/dx_j = sum_a u_i(a) dN_a/dx_j
        //
        // The weak diffusion contribution is:
        //
        //     r_diff,i,a^(e)
        //       = -nu_e integral(V_e)
        //           grad(N_a) . grad(u_i) dV
        //
        //         +nu_e integral(Gamma_e)
        //           N_a grad(u_i) . n dGamma
        //
        // where:
        //
        //     nu_e = mu_e / rho_e      dimensional material mode
        //     nu_e = ani               non-dimensional mode
        //=====================================================================
        void step1_diffusion(
            const CBSStateSI& s,
            Int ie,
            const Real lunkno[4][5],
            Real ani,
            Real lrhs[4][5])
        {
            Real dNuidxj[4][4] = {};

            // Calculate the constant element velocity gradient:
            //
            //     dNuidxj(i,j) = du_i / dx_j
            for (Int i = 1; i <= s.cfg.ndim; ++i)
            {
                for (Int a = 1; a <= s.cfg.nep; ++a)
                {
                    for (Int j = 1; j <= s.cfg.ndim; ++j)
                    {
                        dNuidxj[i][j] += grad(s, ie, j, a) * lunkno[i][a];
                    }
                }
            }

            const Real volume_factor = ani * s.detJ(ie) * s.cfg.fdif[1];

            // Add the element-volume diffusion contribution.
            for (Int i = 1; i <= s.cfg.ndim; ++i)
            {
                for (Int a = 1; a <= s.cfg.nep; ++a)
                {
                    for (Int j = 1; j <= s.cfg.ndim; ++j)
                    {
                        lrhs[i][a] -=
                            grad(s, ie, j, a) * dNuidxj[i][j] * volume_factor;
                    }
                }
            }

            const Real face_diffusion_factor = ani * s.cfg.fdif[2];

            // Add the natural boundary diffusion contribution.
            for (Int is = 1; is <= s.cfg.nsid; ++is)
            {
                if (s.fedge(is, ie) == 0)
                {
                    continue;
                }

                for (Int j = 1; j <= s.cfg.ndim; ++j)
                {
                    const Real face_factor =
                        s.annxf(j, is, ie) * face_diffusion_factor;

                    for (Int i = 1; i <= s.cfg.ndim; ++i)
                    {
                        for (Int face_node = 1; face_node <= s.cfg.nsidp; ++face_node)
                        {
                            const Int a = s.ippn1(is, face_node);
                            lrhs[i][a] += dNuidxj[i][j] * face_factor;
                        }
                    }
                }
            }
        }
    }


    //=========================================================================
    // Assembles the complete CBS Step 1 momentum residual.
    //
    // For each fluid tetrahedron:
    //
    //     1. Gather the previous nodal velocity u^n.
    //     2. Form the nodal products u_j u_i.
    //     3. Add Galerkin convection.
    //     4. Add the CBS characteristic correction.
    //     5. Add viscous diffusion.
    //     6. Scatter the local residual into rhs.
    //
    // Elements are processed colour by colour. Within one colour, no two
    // tetrahedra share a node, so the nodal scatter is race-free under OpenMP.
    //
    // Inputs:
    //     unkn1       velocity at the previous CBS iteration
    //     dNkdx       shape-function gradients
    //     detJ        tetrahedral Jacobian determinant
    //     delte       local element time step
    //     annxf       area-weighted boundary normals
    //     mat_elem    material identifier
    //
    // Output:
    //     rhs         global nodal momentum residual
    //=========================================================================
    void MomentumAssembly::assembleStep1Rhs(CBSStateSI& s)
    {
        validate_step1_dimensions(s);

        s.rhs.fill(0.0);

        if (s.cfg.dtreal <= 0.0 || !std::isfinite(s.cfg.dtreal))
        {
            throw std::runtime_error(
                "MomentumAssembly::assembleStep1Rhs - dtreal must be positive before Step 1 assembly");
        }

        bool bad_detj = false;

        for (Int c = 0; c < s.ncolor; ++c)
        {
            const Int cbeg = s.color_ptr[static_cast<std::size_t>(c)];
            const Int cend = s.color_ptr[static_cast<std::size_t>(c) + 1];
#pragma omp parallel for schedule(static)
            for (Int k = cbeg; k < cend; ++k)
            {
                const Int ie = s.color_elem[static_cast<std::size_t>(k)];

                if (!is_fluid_element(s, ie))
                {
                    continue;
                }

                if (s.detJ(ie) <= 0.0 || !std::isfinite(s.detJ(ie)))
                {
                    bad_detj = true;
                    continue;
                }

                Real lunkno[4][5] = {};
                Real lunk[4][4][5] = {};
                Real lrhs[4][5] = {};

                Int ipn[5] = {};

                // Gather the previous nodal velocity for the current element.
                for (Int a = 1; a <= s.cfg.nep; ++a)
                {
                    const Int ip = s.intma(a, ie);
                    ipn[a] = ip;
                    for (Int i = 1; i <= s.cfg.ndim; ++i)
                    {
                        lunkno[i][a] = s.unkn1(i, ip);
                    }
                }

                // Form the nodal momentum-flux products:
                //
                //     lunk(i,j,a) = u_j(a) u_i(a)
                for (Int a = 1; a <= s.cfg.nep; ++a)
                {
                    for (Int j = 1; j <= s.cfg.ndim; ++j)
                    {
                        for (Int i = 1; i <= s.cfg.ndim; ++i)
                        {
                            lunk[i][j][a] = lunkno[j][a] * lunkno[i][a];
                        }
                    }
                }

                step1_convective_galerkin(s, ie, lunk, lrhs);
                step1_characteristic_correction(s, ie, lunk, lunkno, lrhs);
                step1_diffusion(s, ie, lunkno, momentum_diffusivity(s, ie), lrhs);

                // Scatter the local residual into the global nodal array.
                // Element colouring prevents shared-node writes within the
                // current OpenMP loop.
                for (Int a = 1; a <= s.cfg.nep; ++a)
                {
                    const Int ip = ipn[a];
                    for (Int i = 1; i <= s.cfg.ndim; ++i)
                    {
                        s.rhs(i, ip) += lrhs[i][a];
                    }
                }
            }
        }

        if (bad_detj)
        {
            throw std::runtime_error(
                "MomentumAssembly::assembleStep1Rhs - invalid detJ at one or more elements");
        }

        // Artificial diffusion is not assembled in the current core Step 1
        // implementation. It requires the additional stabilisation arrays and
        // should be introduced only as a separately validated numerical term.
    }


    //=========================================================================
    // Reserved hook for additional real-time momentum terms.
    //
    // The current implementation intentionally performs no operation. The
    // Step 1 nodal update is carried out in Steps::step1SemiImplicit() using
    // the inverse diagonal elcoe2.
    //=========================================================================
    void MomentumAssembly::applyRealTimeMomentumTerm(CBSStateSI& s)
    {
        // BDF and consistent-mass correction terms are handled in the legacy
        // step1ggMPI driver after RHS assembly.  The first CBS3D++_SI serial
        // port keeps Step 1 assembly separate from the nodal update, matching
        // CBS2D++_SI.  Solver/Steps code will apply the update using elcoe2.
        (void)s;
    }
}
