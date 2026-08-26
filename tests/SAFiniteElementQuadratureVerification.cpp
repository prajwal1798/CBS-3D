//=============================================================================
// SA P1/TET4 nonlinear-source quadrature regression.
//
// The historical implementation evaluated production/destruction once at the
// element centroid and distributed V/4 to every node.  That is not the finite-
// element integral of N_a P(q,d) or N_a D(q,d) when q and wall distance vary
// through a near-wall tetrahedron.  This test independently reconstructs the
// four-point positive tetrahedron quadrature used by the verified assembly and
// checks the assembled nodal production/destruction exactly.
//=============================================================================

#include "cbs/assembly/SpalartAllmarasAssembly.hpp"
#include "cbs/turbulence/SpalartAllmaras.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>

namespace
{
    using cbs::CBSStateSI;
    using cbs::Int;
    using cbs::Real;

    constexpr Real qa = 0.58541019662496845446;
    constexpr Real qb = 0.13819660112501051518;
    constexpr Real qw = 0.25;

    bool close(const Real a, const Real b)
    {
        const Real scale = std::max({Real(1.0e-30), std::abs(a), std::abs(b)});
        return std::abs(a - b) <= 5.0e-12 * scale + 1.0e-24;
    }

    std::array<Real, 4> shape(const Int qp)
    {
        std::array<Real, 4> n = {qb, qb, qb, qb};
        n[static_cast<std::size_t>(qp)] = qa;
        return n;
    }

    CBSStateSI make_state()
    {
        CBSStateSI s;
        s.cfg.turbulence_on = 1;
        s.cfg.turbulence_model = 0;
        s.cfg.dimensional_mode = 0;
        s.cfg.material_properties_enabled = 0;
        s.cfg.ani = 1.0e-5;
        s.cfg.sa_nu_tilde_floor = 0.0;
        s.cfg.sa_min_wall_distance = 1.0e-14;
        s.cfg.sa_min_stilde = 1.0e-14;
        s.cfg.sa_use_stilde_limiter = 1;
        s.cfg.sa_implicit_destruction = 0;

        s.initialise_local_topology();
        s.set_problem_sizes(1, 4, 0, 0);

        s.ncolor = 1;
        s.color_ptr = {0, 1};
        s.color_elem = {1};

        for (Int a = 1; a <= 4; ++a)
        {
            s.intma(a, 1) = a;
        }
        s.mat_elem(1) = 0;
        s.detJ(1) = 1.0; // V = 1/6
        s.delte(1) = 1.0e-4;

        // Unit reference tetrahedron:
        // N1=1-x-y-z, N2=x, N3=y, N4=z.
        const Real g[3][4] =
        {
            {-1.0, 1.0, 0.0, 0.0},
            {-1.0, 0.0, 1.0, 0.0},
            {-1.0, 0.0, 0.0, 1.0}
        };
        for (Int j = 1; j <= 3; ++j)
        {
            for (Int a = 1; a <= 4; ++a)
            {
                const Int index = (j - 1) * 4 + a;
                s.dNkdx(index) = g[j - 1][a - 1];
            }
        }

        s.fedge.fill(0);
        s.annxf.fill(0.0);
        s.unkno.fill(0.0); // omega = 0, isolates nonlinear wall/source algebra.

        const Real q[4] = {0.8e-5, 1.7e-5, 3.1e-5, 5.2e-5};
        const Real d[4] = {2.0e-4, 8.0e-4, 2.5e-3, 7.0e-3};

        for (Int a = 1; a <= 4; ++a)
        {
            s.nu_tilde1(a) = q[a - 1];
            s.nu_tilde(a) = q[a - 1];
            s.wall_distance(a) = d[a - 1];
            s.sa_active_node(a) = 1;
            s.sa_wall_node(a) = 0;
            s.sa_inlet_node(a) = 0;
        }

        return s;
    }
}

int main()
{
    CBSStateSI s = make_state();
    cbs::SpalartAllmarasAssembly::assembleTransportRhs(s);

    const Real volume = 1.0 / 6.0;
    const Real molecular_nu = s.cfg.ani;
    const Real qnode[4] = {0.8e-5, 1.7e-5, 3.1e-5, 5.2e-5};
    const Real dnode[4] = {2.0e-4, 8.0e-4, 2.5e-3, 7.0e-3};

    cbs::turbulence::SpalartAllmarasConstants c;
    Real expected_p[4] = {};
    Real expected_d[4] = {};

    for (Int qp = 0; qp < 4; ++qp)
    {
        const auto n = shape(qp);
        Real q = 0.0;
        Real d = 0.0;
        for (Int a = 0; a < 4; ++a)
        {
            q += n[static_cast<std::size_t>(a)] * qnode[a];
            d += n[static_cast<std::size_t>(a)] * dnode[a];
        }

        const Real sbar = cbs::turbulence::sBar(q, molecular_nu, d, c);
        Real stilde = cbs::turbulence::limitedSTilde(0.0, sbar, c);
        stilde = std::max(stilde, s.cfg.sa_min_stilde);
        const Real chi = cbs::turbulence::chi(q, molecular_nu);
        const Real ft2 = cbs::turbulence::ft2(chi, c);
        const Real r = cbs::turbulence::rFunction(q, stilde, d, c);
        const Real fw = cbs::turbulence::fw(r, c);
        const Real p = cbs::turbulence::productionTerm(q, stilde, ft2, c);
        const Real destr = cbs::turbulence::destructionTerm(q, d, fw, ft2, c);

        for (Int a = 0; a < 4; ++a)
        {
            const Real weight = volume * qw * n[static_cast<std::size_t>(a)];
            expected_p[a] += weight * p;
            expected_d[a] += weight * destr;
        }
    }

    // Also reconstruct the retired centroid approximation and require the test
    // state to distinguish it.  Otherwise this regression would be ceremonial.
    Real qbar = 0.0;
    Real dbar = 0.0;
    for (Int a = 0; a < 4; ++a)
    {
        qbar += 0.25 * qnode[a];
        dbar += 0.25 * dnode[a];
    }
    const Real sbar_c = cbs::turbulence::sBar(qbar, molecular_nu, dbar, c);
    Real stilde_c = cbs::turbulence::limitedSTilde(0.0, sbar_c, c);
    stilde_c = std::max(stilde_c, s.cfg.sa_min_stilde);
    const Real chi_c = cbs::turbulence::chi(qbar, molecular_nu);
    const Real ft2_c = cbs::turbulence::ft2(chi_c, c);
    const Real r_c = cbs::turbulence::rFunction(qbar, stilde_c, dbar, c);
    const Real fw_c = cbs::turbulence::fw(r_c, c);
    const Real p_c = cbs::turbulence::productionTerm(qbar, stilde_c, ft2_c, c);
    const Real d_c = cbs::turbulence::destructionTerm(qbar, dbar, fw_c, ft2_c, c);

    Real centroid_difference = 0.0;
    for (Int a = 1; a <= 4; ++a)
    {
        if (!close(s.sa_production(a), expected_p[a - 1]))
        {
            std::printf("FAIL production node %d: got %.17e expected %.17e\n",
                        a, s.sa_production(a), expected_p[a - 1]);
            return 1;
        }
        if (!close(s.sa_destruction(a), expected_d[a - 1]))
        {
            std::printf("FAIL destruction node %d: got %.17e expected %.17e\n",
                        a, s.sa_destruction(a), expected_d[a - 1]);
            return 2;
        }

        centroid_difference +=
            std::abs(expected_p[a - 1] - volume * 0.25 * p_c) +
            std::abs(expected_d[a - 1] - volume * 0.25 * d_c);
    }

    if (!(centroid_difference > 1.0e-14))
    {
        std::printf("FAIL: regression state does not distinguish centroid quadrature\n");
        return 3;
    }

    std::printf("PASS: SA nonlinear source quadrature and nodal test-function weighting\n");
    return 0;
}
