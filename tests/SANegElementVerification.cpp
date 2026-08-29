//=============================================================================
// Deterministic one-element verification of the NASA/TMR SA-neg negative branch
// inside the production Liu/Nithiarasu P1/TET4 scalar transport assembly.
//
// The test deliberately uses a skew tetrahedron, asymmetric negative nu_tilde
// values and nonzero velocity gradients.  Expected element residuals are
// evaluated independently from the production SA helper routines.
//=============================================================================

#include "cbs/assembly/SpalartAllmarasAssembly.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace
{
    using cbs::CBSStateSI;
    using cbs::Int;
    using cbs::Real;

    constexpr Real cb1 = 0.1355;
    constexpr Real sigma = 2.0 / 3.0;
    constexpr Real cb2 = 0.622;
    constexpr Real kappa = 0.41;
    constexpr Real cw2 = 0.3;
    constexpr Real cw3 = 2.0;
    constexpr Real ct3 = 1.2;
    constexpr Real cn1 = 16.0;

    struct Geometry
    {
        Real det_j = 0.0;
        Real volume = 0.0;
        Real grad[5][4] = {};
    };

    Real determinant3(const Real a[4][4])
    {
        return
            a[1][1] * (a[2][2] * a[3][3] - a[2][3] * a[3][2])
          - a[1][2] * (a[2][1] * a[3][3] - a[2][3] * a[3][1])
          + a[1][3] * (a[2][1] * a[3][2] - a[2][2] * a[3][1]);
    }

    Geometry build_geometry(const Real xyz[5][4])
    {
        Real j[4][4] = {};
        for (Int row = 1; row <= 3; ++row)
        {
            j[row][1] = xyz[2][row] - xyz[1][row];
            j[row][2] = xyz[3][row] - xyz[1][row];
            j[row][3] = xyz[4][row] - xyz[1][row];
        }

        Geometry g;
        g.det_j = determinant3(j);
        g.volume = g.det_j / 6.0;

        Real inv[4][4] = {};
        inv[1][1] =  (j[2][2] * j[3][3] - j[2][3] * j[3][2]) / g.det_j;
        inv[1][2] = -(j[1][2] * j[3][3] - j[1][3] * j[3][2]) / g.det_j;
        inv[1][3] =  (j[1][2] * j[2][3] - j[1][3] * j[2][2]) / g.det_j;
        inv[2][1] = -(j[2][1] * j[3][3] - j[2][3] * j[3][1]) / g.det_j;
        inv[2][2] =  (j[1][1] * j[3][3] - j[1][3] * j[3][1]) / g.det_j;
        inv[2][3] = -(j[1][1] * j[2][3] - j[1][3] * j[2][1]) / g.det_j;
        inv[3][1] =  (j[2][1] * j[3][2] - j[2][2] * j[3][1]) / g.det_j;
        inv[3][2] = -(j[1][1] * j[3][2] - j[1][2] * j[3][1]) / g.det_j;
        inv[3][3] =  (j[1][1] * j[2][2] - j[1][2] * j[2][1]) / g.det_j;

        for (Int dim = 1; dim <= 3; ++dim)
        {
            g.grad[2][dim] = inv[1][dim];
            g.grad[3][dim] = inv[2][dim];
            g.grad[4][dim] = inv[3][dim];
            g.grad[1][dim] =
                -g.grad[2][dim] - g.grad[3][dim] - g.grad[4][dim];
        }
        return g;
    }

    Real consistent_load(const Real volume, const Real f[5], const Int a)
    {
        Real sum = 0.0;
        for (Int b = 1; b <= 4; ++b)
        {
            sum += f[b];
        }
        return volume * (sum + f[a]) / 20.0;
    }

    Real product_integral(const Real volume, const Real f[5], const Real h[5])
    {
        Real sf = 0.0;
        Real sh = 0.0;
        Real sfh = 0.0;
        for (Int a = 1; a <= 4; ++a)
        {
            sf += f[a];
            sh += h[a];
            sfh += f[a] * h[a];
        }
        return volume * (sf * sh + sfh) / 20.0;
    }

    bool close_small(const Real got, const Real expected)
    {
        const Real scale = std::max({Real(1.0e-18), std::abs(got), std::abs(expected)});
        return std::abs(got - expected) <= 8.0e-11 * scale + 2.0e-22;
    }

    CBSStateSI make_state(
        const Real xyz[5][4],
        const Geometry& g,
        const Real velocity[5][4],
        const Real q[5],
        const Real distance[5],
        const Real nu,
        const Real dt)
    {
        CBSStateSI s;
        s.cfg.turbulence_on = 1;
        s.cfg.turbulence_model = 1;
        s.cfg.dimensional_mode = 0;
        s.cfg.material_properties_enabled = 0;
        s.cfg.ani = nu;
        s.cfg.sa_nu_tilde_floor = 0.0;
        s.cfg.sa_min_wall_distance = 1.0e-14;
        s.cfg.sa_min_stilde = 1.0e-14;

        s.initialise_local_topology();
        s.set_problem_sizes(1, 4, 0, 0);
        s.ncolor = 1;
        s.color_ptr = {0, 1};
        s.color_elem = {1};

        s.mat_elem(1) = 0;
        s.detJ(1) = g.det_j;
        s.delte(1) = dt;

        for (Int a = 1; a <= 4; ++a)
        {
            s.intma(a, 1) = a;
            s.sa_active_node(a) = 1;
            s.sa_wall_node(a) = 0;
            s.sa_inlet_node(a) = 0;
            s.wall_distance(a) = distance[a];
            s.nu_tilde1(a) = q[a];
            s.nu_tilde(a) = q[a];
            s.elcoe2(a) = 1.0;

            for (Int dim = 1; dim <= 3; ++dim)
            {
                s.coord(dim, a) = xyz[a][dim];
                s.unkno(dim, a) = velocity[a][dim];
                s.unkn1(dim, a) = velocity[a][dim];
                s.dNkdx((dim - 1) * 4 + a) = g.grad[a][dim];
            }
        }

        return s;
    }
}

int main()
{
    const Real xyz[5][4] =
    {
        {},
        {0.0, 0.10, -0.20, 0.30},
        {0.0, 1.20,  0.10, 0.20},
        {0.0, 0.20,  1.40, 0.40},
        {0.0, 0.30,  0.20, 1.70}
    };

    const Real velocity[5][4] =
    {
        {},
        {0.0,  0.37, -0.11,  0.23},
        {0.0,  1.12,  0.29, -0.07},
        {0.0,  0.08,  0.91,  0.34},
        {0.0, -0.16,  0.27,  1.05}
    };

    const Real q[5] =
    {
        0.0,
        -1.0e-6,
        -2.2e-6,
        -4.1e-6,
        -3.0e-6
    };

    const Real distance[5] =
    {
        0.0,
        4.0e-4,
        1.3e-3,
        2.2e-3,
        3.7e-3
    };

    const Real nu = 1.7e-5;
    const Real dt = 2.3e-4;
    const Geometry g = build_geometry(xyz);

    if (!(g.det_j > 0.0))
    {
        std::printf("FAIL SA-neg geometry detJ=%.17e\n", g.det_j);
        return 2;
    }

    CBSStateSI s = make_state(xyz, g, velocity, q, distance, nu, dt);

    Real du[4][4] = {};
    Real grad_q[4] = {};
    Real q_sum = 0.0;
    Real d_sum = 0.0;

    for (Int a = 1; a <= 4; ++a)
    {
        q_sum += q[a];
        d_sum += distance[a];
        for (Int i = 1; i <= 3; ++i)
        {
            for (Int j = 1; j <= 3; ++j)
            {
                du[i][j] += velocity[a][i] * g.grad[a][j];
            }
            grad_q[i] += q[a] * g.grad[a][i];
        }
    }

    const Real qbar = 0.25 * q_sum;
    const Real dbar = 0.25 * d_sum;
    if (!(qbar < 0.0))
    {
        std::printf("FAIL SA-neg test state is not negative\n");
        return 3;
    }

    const Real div_u = du[1][1] + du[2][2] + du[3][3];
    const Real wx = du[3][2] - du[2][3];
    const Real wy = du[1][3] - du[3][1];
    const Real wz = du[2][1] - du[1][2];
    const Real omega = std::sqrt(wx * wx + wy * wy + wz * wz);

    const Real chi = qbar / nu;
    const Real chi3 = chi * chi * chi;
    const Real fn = (cn1 + chi3) / (cn1 - chi3);
    const Real cw1 = cb1 / (kappa * kappa) + (1.0 + cb2) / sigma;
    const Real production_coefficient = cb1 * (1.0 - ct3) * omega;
    const Real recovery_coefficient = cw1 / (dbar * dbar);
    const Real diffusion_coefficient = (nu + qbar * fn) / sigma;

    Real grad_q_sq = 0.0;
    for (Int j = 1; j <= 3; ++j)
    {
        grad_q_sq += grad_q[j] * grad_q[j];
    }
    const Real nonlinear_source = cb2 * grad_q_sq / sigma;

    Real div_uq[5] = {};
    for (Int b = 1; b <= 4; ++b)
    {
        Real advective = 0.0;
        for (Int j = 1; j <= 3; ++j)
        {
            advective += velocity[b][j] * grad_q[j];
        }
        div_uq[b] = advective + q[b] * div_u;
    }

    Real expected_rhs[5] = {};
    Real expected_source[5] = {};

    for (Int a = 1; a <= 4; ++a)
    {
        const Real advection = -consistent_load(g.volume, div_uq, a);

        Real grad_dot = 0.0;
        for (Int j = 1; j <= 3; ++j)
        {
            grad_dot += g.grad[a][j] * grad_q[j];
        }
        const Real diffusion = -g.volume * diffusion_coefficient * grad_dot;
        const Real nonlinear = 0.25 * g.volume * nonlinear_source;

        Real div_u_na[5] = {};
        for (Int b = 1; b <= 4; ++b)
        {
            Real u_dot_grad_na = 0.0;
            for (Int j = 1; j <= 3; ++j)
            {
                u_dot_grad_na += velocity[b][j] * g.grad[a][j];
            }
            div_u_na[b] = u_dot_grad_na + (b == a ? div_u : 0.0);
        }
        const Real characteristic =
            -0.5 * dt * product_integral(g.volume, div_u_na, div_uq);

        const Real q_load = consistent_load(g.volume, q, a);
        const Real production = production_coefficient * q_load;
        const Real recovery = recovery_coefficient * qbar * q_load;
        expected_source[a] = production + recovery;
        expected_rhs[a] =
            advection + characteristic + diffusion + nonlinear
            + expected_source[a];
    }

    cbs::SpalartAllmarasAssembly::assembleTransportRhs(s);

    for (Int a = 1; a <= 4; ++a)
    {
        if (!close_small(s.sa_rhs(a), expected_rhs[a]))
        {
            std::printf(
                "FAIL SA-neg RHS a=%d got=%.17e expected=%.17e\n",
                a, s.sa_rhs(a), expected_rhs[a]);
            return 10 + a;
        }
        if (!close_small(s.sa_source(a), expected_source[a]))
        {
            std::printf(
                "FAIL SA-neg source a=%d got=%.17e expected=%.17e\n",
                a, s.sa_source(a), expected_source[a]);
            return 20 + a;
        }
    }

    // The defining numerical behavior: negative finite states are transported,
    // not hard-clipped to zero.  Force a known update independent of assembly.
    for (Int a = 1; a <= 4; ++a)
    {
        s.nu_tilde1(a) = -2.0e-6;
        s.nu_tilde(a) = -2.0e-6;
        s.sa_rhs(a) = -3.0e-7;
        s.elcoe2(a) = 1.0;
    }

    cbs::SpalartAllmarasAssembly::updateNuTilde(s);

    for (Int a = 1; a <= 4; ++a)
    {
        const Real expected = -2.3e-6;
        if (!close_small(s.nu_tilde(a), expected))
        {
            std::printf(
                "FAIL SA-neg update clipped/altered node=%d got=%.17e expected=%.17e\n",
                a, s.nu_tilde(a), expected);
            return 30 + a;
        }
    }

    // Negative SA working variable must not create negative eddy viscosity.
    cbs::SpalartAllmarasAssembly::updateEddyViscosity(s);
    if (!(s.mu_t_e(1) == 0.0 && s.nu_t_e(1) == 0.0))
    {
        std::printf(
            "FAIL SA-neg negative state produced eddy viscosity: nu_t=%.17e mu_t=%.17e\n",
            s.nu_t_e(1), s.mu_t_e(1));
        return 40;
    }

    std::printf("PASS NASA/TMR SA-neg negative-branch element residual\n");
    std::printf("PASS SA-neg preserves negative transported states\n");
    std::printf("PASS SA-neg negative state supplies zero eddy viscosity to momentum\n");
    return 0;
}
