//=============================================================================
// Deterministic one-element verification of the Liu/Nithiarasu 3-D P1/TET4
// RANS-CBS and Spalart-Allmaras production assemblies.
//
// Reference: Chun-Bin Liu (2005), Appendix B.
//
// IMPORTANT: the expected residuals below are evaluated independently from the
// production helper kernels.  The test deliberately uses a skew tetrahedron and
// asymmetric nodal fields so that index swaps, sign mistakes and missing factors
// cannot be hidden by symmetry.
//=============================================================================

#include "cbs/assembly/MomentumAssembly.hpp"
#include "cbs/assembly/SpalartAllmarasAssembly.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>

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
    constexpr Real cv1 = 7.1;

    struct Geometry
    {
        Real det_j = 0.0;
        Real volume = 0.0;
        Real grad[5][4] = {};
    };

    Real det3(const Real a[4][4])
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
        g.det_j = det3(j);
        if (!(g.det_j > 0.0))
        {
            std::printf("FAIL geometry: detJ=%.17e\n", g.det_j);
            std::exit(2);
        }
        g.volume = g.det_j / 6.0;

        // J^{-1}: rows are gradients of xi, eta, zeta with respect to x,y,z.
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

        // N2=xi, N3=eta, N4=zeta; N1=1-xi-eta-zeta.
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

    bool close(const Real got, const Real expected)
    {
        const Real scale = std::max({Real(1.0), std::abs(got), std::abs(expected)});
        return std::abs(got - expected) <= 2.0e-11 * scale;
    }

    bool close_small(const Real got, const Real expected)
    {
        const Real scale = std::max({Real(1.0e-16), std::abs(got), std::abs(expected)});
        return std::abs(got - expected) <= 5.0e-11 * scale + 1.0e-22;
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

    Real p1_product(const Real volume, const Real f[5], const Real h[5])
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

    Real fv1_reference(const Real chi)
    {
        if (!(chi > 0.0))
        {
            return 0.0;
        }
        const Real chi3 = chi * chi * chi;
        const Real cv13 = cv1 * cv1 * cv1;
        return chi3 / (chi3 + cv13);
    }

    Real fv2_reference(const Real chi, const Real fv1)
    {
        return 1.0 - chi / (1.0 + chi * fv1);
    }

    Real fw_reference(Real r)
    {
        r = std::clamp(r, Real(0.0), Real(10.0));
        const Real r6 = std::pow(r, 6.0);
        const Real gg = r + cw2 * (r6 - r);
        const Real gg6 = std::pow(gg, 6.0);
        const Real cw36 = std::pow(cw3, 6.0);
        return gg * std::pow((1.0 + cw36) / (gg6 + cw36), 1.0 / 6.0);
    }

    CBSStateSI make_state(
        const Real xyz[5][4],
        const Geometry& geometry,
        const Real velocity[5][4],
        const Real q[5],
        const Real distance[5],
        const Real molecular_nu,
        const Real dt)
    {
        CBSStateSI s;
        s.cfg.turbulence_on = 1;
        s.cfg.turbulence_model = 0;
        s.cfg.dimensional_mode = 0;
        s.cfg.material_properties_enabled = 0;
        s.cfg.ani = molecular_nu;
        s.cfg.dtreal = dt;
        s.cfg.sa_nu_tilde_floor = 0.0;
        s.cfg.sa_min_wall_distance = 1.0e-14;
        s.cfg.sa_min_stilde = 1.0e-14;

        s.initialise_local_topology();
        s.set_problem_sizes(1, 4, 0, 0);

        s.ncolor = 1;
        s.color_ptr = {0, 1};
        s.color_elem = {1};

        for (Int a = 1; a <= 4; ++a)
        {
            s.intma(a, 1) = a;
            for (Int dim = 1; dim <= 3; ++dim)
            {
                s.coord(dim, a) = xyz[a][dim];
                s.unkn1(dim, a) = velocity[a][dim];
                s.unkno(dim, a) = velocity[a][dim];
            }

            s.nu_tilde1(a) = q[a];
            s.nu_tilde(a) = q[a];
            s.wall_distance(a) = distance[a];
            s.sa_active_node(a) = 1;
            s.sa_wall_node(a) = 0;
            s.sa_inlet_node(a) = 0;
        }

        s.mat_elem(1) = 0;
        s.detJ(1) = geometry.det_j;
        s.delte(1) = dt;
        s.fedge.fill(0);
        s.annxf.fill(0.0);

        for (Int dim = 1; dim <= 3; ++dim)
        {
            for (Int a = 1; a <= 4; ++a)
            {
                const Int index = (dim - 1) * 4 + a;
                s.dNkdx(index) = geometry.grad[a][dim];
            }
        }
        return s;
    }

    int verify_momentum(
        CBSStateSI& s,
        const Geometry& geometry,
        const Real velocity[5][4],
        const Real molecular_nu,
        const Real dt)
    {
        Real du[4][4] = {};
        for (Int i = 1; i <= 3; ++i)
        {
            for (Int j = 1; j <= 3; ++j)
            {
                for (Int a = 1; a <= 4; ++a)
                {
                    du[i][j] += velocity[a][i] * geometry.grad[a][j];
                }
            }
        }
        const Real div_u = du[1][1] + du[2][2] + du[3][3];

        Real tau[4][4] = {};
        for (Int i = 1; i <= 3; ++i)
        {
            for (Int j = 1; j <= 3; ++j)
            {
                tau[i][j] = molecular_nu *
                    (du[i][j] + du[j][i]
                     - (i == j ? (2.0 / 3.0) * div_u : 0.0));
            }
        }

        Real expected[4][5] = {};
        for (Int i = 1; i <= 3; ++i)
        {
            Real div_momentum[5] = {};
            for (Int b = 1; b <= 4; ++b)
            {
                Real advective = 0.0;
                for (Int j = 1; j <= 3; ++j)
                {
                    advective += velocity[b][j] * du[i][j];
                }
                div_momentum[b] = advective + velocity[b][i] * div_u;
            }

            for (Int a = 1; a <= 4; ++a)
            {
                const Real convection =
                    -consistent_load(geometry.volume, div_momentum, a);

                Real div_u_na[5] = {};
                for (Int b = 1; b <= 4; ++b)
                {
                    Real u_dot_grad_na = 0.0;
                    for (Int j = 1; j <= 3; ++j)
                    {
                        u_dot_grad_na += velocity[b][j] * geometry.grad[a][j];
                    }
                    div_u_na[b] =
                        u_dot_grad_na + (b == a ? div_u : 0.0);
                }

                const Real characteristic =
                    -0.5 * dt * p1_product(
                        geometry.volume,
                        div_u_na,
                        div_momentum);

                Real stress_contraction = 0.0;
                for (Int j = 1; j <= 3; ++j)
                {
                    stress_contraction += geometry.grad[a][j] * tau[i][j];
                }
                const Real diffusion =
                    -geometry.volume * stress_contraction;

                expected[i][a] = convection + characteristic + diffusion;
            }
        }

        cbs::MomentumAssembly::assembleStep1Rhs(s);

        for (Int i = 1; i <= 3; ++i)
        {
            for (Int a = 1; a <= 4; ++a)
            {
                if (!close(s.rhs(i, a), expected[i][a]))
                {
                    std::printf(
                        "FAIL Liu momentum i=%d a=%d got=%.17e expected=%.17e\n",
                        i, a, s.rhs(i, a), expected[i][a]);
                    return 10 + 4 * (i - 1) + a;
                }
            }
        }

        std::printf("PASS Liu momentum B.3/B.5/B.13 one-TET residual\n");
        return 0;
    }

    int verify_sa(
        CBSStateSI& s,
        const Geometry& geometry,
        const Real velocity[5][4],
        const Real q[5],
        const Real distance[5],
        const Real molecular_nu,
        const Real dt)
    {
        Real du[4][4] = {};
        for (Int i = 1; i <= 3; ++i)
        {
            for (Int j = 1; j <= 3; ++j)
            {
                for (Int a = 1; a <= 4; ++a)
                {
                    du[i][j] += velocity[a][i] * geometry.grad[a][j];
                }
            }
        }
        const Real div_u = du[1][1] + du[2][2] + du[3][3];

        Real grad_q[4] = {};
        Real q_sum = 0.0;
        Real d_sum = 0.0;
        for (Int a = 1; a <= 4; ++a)
        {
            q_sum += q[a];
            d_sum += distance[a];
            for (Int j = 1; j <= 3; ++j)
            {
                grad_q[j] += q[a] * geometry.grad[a][j];
            }
        }
        const Real qbar = 0.25 * q_sum;
        const Real dbar = 0.25 * d_sum;

        const Real wx = du[3][2] - du[2][3];
        const Real wy = du[1][3] - du[3][1];
        const Real wz = du[2][1] - du[1][2];
        const Real omega = std::sqrt(wx * wx + wy * wy + wz * wz);

        const Real chi = qbar / molecular_nu;
        const Real fv1 = fv1_reference(chi);
        const Real fv2 = fv2_reference(chi, fv1);
        const Real sbar = qbar * fv2 / (kappa * kappa * dbar * dbar);
        const Real stilde = std::max(omega + sbar, s.cfg.sa_min_stilde);

        Real r = qbar / (stilde * kappa * kappa * dbar * dbar);
        r = std::clamp(r, Real(0.0), Real(10.0));
        const Real fw = fw_reference(r);
        const Real cw1 = cb1 / (kappa * kappa) + (1.0 + cb2) / sigma;

        const Real production_coefficient = cb1 * stilde; // Liu/SA-noft2
        const Real destruction_coefficient = cw1 * fw / (dbar * dbar);
        const Real diffusion_coefficient = (molecular_nu + qbar) / sigma;

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
        Real expected_prod[5] = {};
        Real expected_dest[5] = {};
        Real expected_diff[5] = {};

        for (Int a = 1; a <= 4; ++a)
        {
            const Real advection =
                -consistent_load(geometry.volume, div_uq, a);

            Real grad_dot = 0.0;
            for (Int j = 1; j <= 3; ++j)
            {
                grad_dot += geometry.grad[a][j] * grad_q[j];
            }
            const Real diffusion =
                -geometry.volume * diffusion_coefficient * grad_dot;
            const Real nonlinear =
                geometry.volume * 0.25 * nonlinear_source;

            Real div_u_na[5] = {};
            for (Int b = 1; b <= 4; ++b)
            {
                Real u_dot_grad_na = 0.0;
                for (Int j = 1; j <= 3; ++j)
                {
                    u_dot_grad_na += velocity[b][j] * geometry.grad[a][j];
                }
                div_u_na[b] =
                    u_dot_grad_na + (b == a ? div_u : 0.0);
            }
            const Real characteristic =
                -0.5 * dt * p1_product(
                    geometry.volume,
                    div_u_na,
                    div_uq);

            const Real q_load = consistent_load(geometry.volume, q, a);
            const Real production = production_coefficient * q_load;
            const Real destruction = destruction_coefficient * qbar * q_load;

            expected_prod[a] = production;
            expected_dest[a] = destruction;
            expected_diff[a] = diffusion + nonlinear;
            expected_rhs[a] =
                advection + characteristic + diffusion + nonlinear
                + production - destruction;
        }

        cbs::SpalartAllmarasAssembly::assembleTransportRhs(s);

        for (Int a = 1; a <= 4; ++a)
        {
            if (!close_small(s.sa_rhs(a), expected_rhs[a]))
            {
                std::printf(
                    "FAIL Liu SA RHS a=%d got=%.17e expected=%.17e\n",
                    a, s.sa_rhs(a), expected_rhs[a]);
                return 30 + a;
            }
            if (!close_small(s.sa_production(a), expected_prod[a]))
            {
                std::printf(
                    "FAIL Liu SA production a=%d got=%.17e expected=%.17e\n",
                    a, s.sa_production(a), expected_prod[a]);
                return 40 + a;
            }
            if (!close_small(s.sa_destruction(a), expected_dest[a]))
            {
                std::printf(
                    "FAIL Liu SA destruction a=%d got=%.17e expected=%.17e\n",
                    a, s.sa_destruction(a), expected_dest[a]);
                return 50 + a;
            }
            if (!close_small(s.sa_diffusion(a), expected_diff[a]))
            {
                std::printf(
                    "FAIL Liu SA diffusion a=%d got=%.17e expected=%.17e\n",
                    a, s.sa_diffusion(a), expected_diff[a]);
                return 60 + a;
            }
        }

        // Prevent a vacuous source test: unequal q_a must produce unequal
        // consistent source loads for the frozen element closure coefficient.
        if (close_small(expected_prod[1], expected_prod[4]))
        {
            std::printf("FAIL Liu SA test state is insufficiently asymmetric\n");
            return 70;
        }

        std::printf("PASS Liu SA B.50/B.53/B.56/B.60 one-TET residual\n");
        std::printf(
            "INFO SA vorticity convention checked as |curl(u)| = %.17e\n",
            omega);
        return 0;
    }
}

int main()
{
    // Deliberately skew, positively oriented TET4.
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
        2.1e-5,
        4.7e-5,
        1.3e-5,
        6.2e-5
    };

    const Real distance[5] =
    {
        0.0,
        4.0e-4,
        1.3e-3,
        2.2e-3,
        3.7e-3
    };

    const Real molecular_nu = 1.7e-5;
    const Real dt = 2.3e-4;

    const Geometry geometry = build_geometry(xyz);
    CBSStateSI state = make_state(
        xyz, geometry, velocity, q, distance, molecular_nu, dt);

    const int momentum_status = verify_momentum(
        state, geometry, velocity, molecular_nu, dt);
    if (momentum_status != 0)
    {
        return momentum_status;
    }

    // Momentum assembly overwrites only rhs.  The same deterministic state can
    // therefore be used for the independent SA residual check.
    const int sa_status = verify_sa(
        state, geometry, velocity, q, distance, molecular_nu, dt);
    if (sa_status != 0)
    {
        return sa_status;
    }

    std::printf("PASS Liu/Nithiarasu deterministic element verification\n");
    return 0;
}
