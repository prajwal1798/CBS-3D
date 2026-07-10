#include "cbs/turbulence/WallDistance.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cbs::turbulence
{
    namespace
    {
        [[nodiscard]] std::array<Real, 3> node_point(const CBSStateSI& s, const Int ip)
        {
            return {s.coord(1, ip), s.coord(2, ip), s.coord(3, ip)};
        }

        [[nodiscard]] bool is_wall_bc(const CBSStateSI& s, const Int bc)
        {
            return bc == s.cfg.bc_noslip_adiabatic_wall
                || bc == s.cfg.bc_noslip_heatflux_wall
                || bc == s.cfg.bc_cht_interface;
        }

        [[nodiscard]] std::array<Real, 3> sub(
            const std::array<Real, 3>& a,
            const std::array<Real, 3>& b)
        {
            return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
        }

        [[nodiscard]] std::array<Real, 3> add_scaled(
            const std::array<Real, 3>& a,
            const std::array<Real, 3>& b,
            const Real scale)
        {
            return {a[0] + scale * b[0], a[1] + scale * b[1], a[2] + scale * b[2]};
        }

        [[nodiscard]] Real dot(
            const std::array<Real, 3>& a,
            const std::array<Real, 3>& b)
        {
            return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        }

        [[nodiscard]] Real distance_squared(
            const std::array<Real, 3>& a,
            const std::array<Real, 3>& b)
        {
            const auto d = sub(a, b);
            return dot(d, d);
        }
    }

    std::vector<WallTriangle> WallDistance::collectPhysicalWallTriangles(const CBSStateSI& s)
    {
        std::vector<WallTriangle> triangles;
        triangles.reserve(static_cast<std::size_t>(s.cfg.nboun));

        for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
        {
            const Int bc = s.iside(s.cfg.bsid, ib);

            if (!is_wall_bc(s, bc))
            {
                continue;
            }

            WallTriangle tri;
            tri.a = node_point(s, s.iside(1, ib));
            tri.b = node_point(s, s.iside(2, ib));
            tri.c = node_point(s, s.iside(3, ib));
            triangles.push_back(tri);
        }

        return triangles;
    }

    Real WallDistance::pointTriangleDistance(
        const std::array<Real, 3>& p,
        const WallTriangle& tri)
    {
        // Real-Time Collision Detection, Christer Ericson, point-triangle
        // closest-point test.  The formula is purely geometric and does not
        // assume structured wall-normal mesh lines.
        const auto ab = sub(tri.b, tri.a);
        const auto ac = sub(tri.c, tri.a);
        const auto ap = sub(p, tri.a);

        const Real d1 = dot(ab, ap);
        const Real d2 = dot(ac, ap);
        if (d1 <= 0.0 && d2 <= 0.0)
        {
            return std::sqrt(distance_squared(p, tri.a));
        }

        const auto bp = sub(p, tri.b);
        const Real d3 = dot(ab, bp);
        const Real d4 = dot(ac, bp);
        if (d3 >= 0.0 && d4 <= d3)
        {
            return std::sqrt(distance_squared(p, tri.b));
        }

        const Real vc = d1 * d4 - d3 * d2;
        if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
        {
            const Real v = d1 / (d1 - d3);
            const auto closest = add_scaled(tri.a, ab, v);
            return std::sqrt(distance_squared(p, closest));
        }

        const auto cp = sub(p, tri.c);
        const Real d5 = dot(ab, cp);
        const Real d6 = dot(ac, cp);
        if (d6 >= 0.0 && d5 <= d6)
        {
            return std::sqrt(distance_squared(p, tri.c));
        }

        const Real vb = d5 * d2 - d1 * d6;
        if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
        {
            const Real w = d2 / (d2 - d6);
            const auto closest = add_scaled(tri.a, ac, w);
            return std::sqrt(distance_squared(p, closest));
        }

        const Real va = d3 * d6 - d5 * d4;
        if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
        {
            const Real w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
            const auto bc = sub(tri.c, tri.b);
            const auto closest = add_scaled(tri.b, bc, w);
            return std::sqrt(distance_squared(p, closest));
        }

        const Real denom = 1.0 / (va + vb + vc);
        const Real v = vb * denom;
        const Real w = vc * denom;
        const auto closest = add_scaled(add_scaled(tri.a, ab, v), ac, w);
        return std::sqrt(distance_squared(p, closest));
    }

    void WallDistance::computeSerial(CBSStateSI& s)
    {
        const auto wall_triangles = collectPhysicalWallTriangles(s);

        if (wall_triangles.empty())
        {
            throw std::runtime_error(
                "WallDistance::computeSerial - no physical no-slip wall triangles found");
        }

        s.wall_distance.fill(std::numeric_limits<Real>::max());

        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            if (s.sa_active_node(ip) == 0)
            {
                continue;
            }

            if (s.sa_wall_node(ip) != 0)
            {
                s.wall_distance(ip) = s.cfg.sa_min_wall_distance;
                continue;
            }

            const auto p = node_point(s, ip);
            Real d_min = std::numeric_limits<Real>::max();

            for (const WallTriangle& tri : wall_triangles)
            {
                d_min = std::min(d_min, pointTriangleDistance(p, tri));
            }

            s.wall_distance(ip) = std::max(d_min, s.cfg.sa_min_wall_distance);
        }
    }
}
