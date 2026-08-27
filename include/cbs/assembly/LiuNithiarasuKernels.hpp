#pragma once

//=============================================================================
// CBS3D++_SI
//
// Small algebraic kernels that encode the 3-D P1/TET4 matrices documented in
// Chun-Bin Liu (2005), Appendix B, for the CBS RANS formulation developed at
// Swansea under P. Nithiarasu.
//
// These routines are deliberately mesh-agnostic.  They make the two formulation
// details that are easiest to lose during a C++ rewrite explicit and testable:
//
//   1. momentum diffusion is the full deviatoric symmetric stress, not three
//      unrelated scalar Laplacians;
//   2. scalar convection/stabilisation uses div(u q) and div(u N_a), exactly as
//      in the C_nu and K_{u nu} matrices of Appendix B.
//=============================================================================

#include "cbs/core/Types.hpp"

namespace cbs
{
    namespace liu_nithiarasu
    {
        // Integral of N_a q_h over a four-node tetrahedron.  For P1 TET4:
        //
        //   int N_a q_h dV = V/20 * (sum_b q_b + q_a).
        inline Real consistent_scalar_load(
            const Real volume,
            const Real q[5],
            const Int a)
        {
            Real sum = 0.0;
            for (Int b = 1; b <= 4; ++b)
            {
                sum += q[b];
            }
            return volume * (sum + q[a]) / 20.0;
        }

        // Integral of N_a f_h for a linear nodal scalar f_h.
        inline Real consistent_linear_load(
            const Real volume,
            const Real f[5],
            const Int a)
        {
            return consistent_scalar_load(volume, f, a);
        }

        // Exact integral of the product of two P1 scalars whose nodal values are
        // f_a and g_a:
        //
        //   int f_h g_h dV = V/20 * [(sum f)(sum g) + sum(f_a g_a)].
        inline Real p1_product_integral(
            const Real volume,
            const Real f[5],
            const Real g[5])
        {
            Real sf = 0.0;
            Real sg = 0.0;
            Real sfg = 0.0;
            for (Int a = 1; a <= 4; ++a)
            {
                sf += f[a];
                sg += g[a];
                sfg += f[a] * g[a];
            }
            return volume * (sf * sg + sfg) / 20.0;
        }

        // Forms the kinematic deviatoric stress used by Liu Appendix B Eq. B.5:
        //
        // tau_ij = nu_eff [du_i/dx_j + du_j/dx_i
        //                   - 2/3 div(u) delta_ij].
        //
        // The gradient arrays use indices 1..3.
        inline void deviatoric_stress(
            const Real velocity_gradient[4][4],
            const Real nu_eff,
            Real tau[4][4])
        {
            const Real div_u =
                velocity_gradient[1][1]
                + velocity_gradient[2][2]
                + velocity_gradient[3][3];

            for (Int i = 1; i <= 3; ++i)
            {
                for (Int j = 1; j <= 3; ++j)
                {
                    tau[i][j] = nu_eff *
                        (velocity_gradient[i][j]
                         + velocity_gradient[j][i]
                         - (i == j ? (2.0 / 3.0) * div_u : 0.0));
                }
            }
        }

        // Nodal values of div(u q) for P1 u and q.  The divergence itself is
        // linear because u_h q_h is quadratic:
        //
        //   div(u q) = u . grad(q) + q div(u).
        inline void conservative_scalar_divergence_nodes(
            const Real velocity[5][4],
            const Real q[5],
            const Real grad_q[4],
            const Real div_u,
            Real div_uq[5])
        {
            for (Int a = 1; a <= 4; ++a)
            {
                div_uq[a] =
                    velocity[a][1] * grad_q[1]
                    + velocity[a][2] * grad_q[2]
                    + velocity[a][3] * grad_q[3]
                    + q[a] * div_u;
            }
        }

        // Nodal values of div(u N_a), used by the Liu CBS stabilisation matrix
        // K_{u nu} = -1/2 int div(u N)^T div(u N) dV (Appendix B Eq. B.60).
        inline void test_function_divergence_nodes(
            const Real velocity[5][4],
            const Real grad_na[4],
            const Real div_u,
            const Int test_node,
            Real div_u_na[5])
        {
            for (Int b = 1; b <= 4; ++b)
            {
                div_u_na[b] =
                    velocity[b][1] * grad_na[1]
                    + velocity[b][2] * grad_na[2]
                    + velocity[b][3] * grad_na[3]
                    + (b == test_node ? div_u : 0.0);
            }
        }
    }
}
