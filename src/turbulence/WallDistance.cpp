#include "cbs/turbulence/WallDistance.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cbs
{
    namespace turbulence
    {
        namespace
        {
            std::array<Real, 3> node_point(const CBSStateSI& s, Int ip)
            {
                std::array<Real, 3> p;

                p[0] = s.coord(1, ip);
                p[1] = s.coord(2, ip);
                p[2] = s.coord(3, ip);

                return p;
            }

            bool is_wall_bc(const CBSStateSI& s, Int bc)
            {
                return bc == s.cfg.bc_noslip_adiabatic_wall
                    || bc == s.cfg.bc_noslip_heatflux_wall
                    || bc == s.cfg.bc_cht_interface;
            }

            std::array<Real, 3> subtract_points(
                const std::array<Real, 3>& a,
                const std::array<Real, 3>& b)
            {
                std::array<Real, 3> result;

                result[0] = a[0] - b[0];
                result[1] = a[1] - b[1];
                result[2] = a[2] - b[2];

                return result;
            }

            std::array<Real, 3> add_scaled_vector(
                const std::array<Real, 3>& a,
                const std::array<Real, 3>& b,
                Real scale)
            {
                std::array<Real, 3> result;

                result[0] = a[0] + scale * b[0];
                result[1] = a[1] + scale * b[1];
                result[2] = a[2] + scale * b[2];

                return result;
            }

            Real dot_product(
                const std::array<Real, 3>& a,
                const std::array<Real, 3>& b)
            {
                return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
            }

            Real distance_squared(
                const std::array<Real, 3>& a,
                const std::array<Real, 3>& b)
            {
                const std::array<Real, 3> d = subtract_points(a, b);
                return dot_product(d, d);
            }

            void check_wall_triangle_list(
                const std::vector<WallTriangle>& wall_triangles,
                const char* routine_name)
            {
                if (wall_triangles.empty())
                {
                    throw std::runtime_error(
                        std::string(routine_name)
                        + " - no physical no-slip wall triangles found");
                }
            }

            void compute_one_node_distance(
                CBSStateSI& s,
                Int ip,
                const std::vector<WallTriangle>& wall_triangles)
            {
                if (s.sa_active_node(ip) == 0)
                {
                    return;
                }

                if (s.sa_wall_node(ip) != 0)
                {
                    s.wall_distance(ip) = s.cfg.sa_min_wall_distance;
                    return;
                }

                const std::array<Real, 3> p = node_point(s, ip);
                Real d_min = std::numeric_limits<Real>::max();

                for (Size i = 0; i < wall_triangles.size(); ++i)
                {
                    const Real distance =
                        WallDistance::pointTriangleDistance(p, wall_triangles[i]);

                    d_min = std::min(d_min, distance);
                }

                s.wall_distance(ip) = std::max(d_min, s.cfg.sa_min_wall_distance);
            }
        }

        //=====================================================================
        // Collects physical wall triangles from the boundary-face list.
        //
        // Only real no-slip wall boundaries are used.  Artificial MPI partition
        // interfaces must never appear here because they are communication
        // surfaces, not physical walls.  The conformal CHT interface is included
        // because the fluid side sees it as a no-slip wall for turbulence.
        //=====================================================================
        std::vector<WallTriangle> WallDistance::collectPhysicalWallTriangles(
            const CBSStateSI& s)
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

        //=====================================================================
        // Computes the exact distance from one point to one triangle.
        //
        // The procedure evaluates which region of the triangle contains the
        // closest point:
        //
        //     1. vertex A region,
        //     2. vertex B region,
        //     3. edge AB region,
        //     4. vertex C region,
        //     5. edge AC region,
        //     6. edge BC region,
        //     7. interior face region.
        //
        // It is a purely geometric test and does not rely on structured mesh
        // lines.  This is important for unstructured tetrahedral meshes.
        //=====================================================================
        Real WallDistance::pointTriangleDistance(
            const std::array<Real, 3>& p,
            const WallTriangle& tri)
        {
            const std::array<Real, 3> ab = subtract_points(tri.b, tri.a);
            const std::array<Real, 3> ac = subtract_points(tri.c, tri.a);
            const std::array<Real, 3> ap = subtract_points(p, tri.a);

            const Real d1 = dot_product(ab, ap);
            const Real d2 = dot_product(ac, ap);

            if (d1 <= 0.0 && d2 <= 0.0)
            {
                return std::sqrt(distance_squared(p, tri.a));
            }

            const std::array<Real, 3> bp = subtract_points(p, tri.b);
            const Real d3 = dot_product(ab, bp);
            const Real d4 = dot_product(ac, bp);

            if (d3 >= 0.0 && d4 <= d3)
            {
                return std::sqrt(distance_squared(p, tri.b));
            }

            const Real vc = d1 * d4 - d3 * d2;

            if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
            {
                const Real v = d1 / (d1 - d3);
                const std::array<Real, 3> closest = add_scaled_vector(tri.a, ab, v);
                return std::sqrt(distance_squared(p, closest));
            }

            const std::array<Real, 3> cp = subtract_points(p, tri.c);
            const Real d5 = dot_product(ab, cp);
            const Real d6 = dot_product(ac, cp);

            if (d6 >= 0.0 && d5 <= d6)
            {
                return std::sqrt(distance_squared(p, tri.c));
            }

            const Real vb = d5 * d2 - d1 * d6;

            if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
            {
                const Real w = d2 / (d2 - d6);
                const std::array<Real, 3> closest = add_scaled_vector(tri.a, ac, w);
                return std::sqrt(distance_squared(p, closest));
            }

            const Real va = d3 * d6 - d5 * d4;

            if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
            {
                const Real w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
                const std::array<Real, 3> bc = subtract_points(tri.c, tri.b);
                const std::array<Real, 3> closest = add_scaled_vector(tri.b, bc, w);
                return std::sqrt(distance_squared(p, closest));
            }

            const Real denominator = 1.0 / (va + vb + vc);
            const Real v = vb * denominator;
            const Real w = vc * denominator;

            const std::array<Real, 3> point_on_ab = add_scaled_vector(tri.a, ab, v);
            const std::array<Real, 3> closest = add_scaled_vector(point_on_ab, ac, w);

            return std::sqrt(distance_squared(p, closest));
        }

        //=====================================================================
        // Computes wall distance for every SA-active node using OpenMP.
        //
        // The input wall triangles are fixed after boundary preprocessing.  Each
        // nodal distance is independent of every other nodal distance.  Therefore
        // the outer node loop is safe to parallelise: every thread writes to a
        // different wall_distance(ip) entry, and the wall-triangle list is only
        // read.
        //
        // Wall nodes receive the small positive lower bound specified by
        // sa_min_wall_distance.  This avoids division by zero in the SA
        // destruction term while still representing the wall boundary value.
        //
        // Non-wall active fluid nodes receive the minimum distance to all
        // physical wall triangles.  This implementation is O(Nnode Nwall).  It
        // is deliberately simple and verifiable before a spatial search tree is
        // introduced.
        //=====================================================================
        void WallDistance::compute(CBSStateSI& s)
        {
            const std::vector<WallTriangle> wall_triangles =
                collectPhysicalWallTriangles(s);

            check_wall_triangle_list(wall_triangles, "WallDistance::compute");

            s.wall_distance.fill(std::numeric_limits<Real>::max());

#ifdef CBS3D_USE_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                compute_one_node_distance(s, ip, wall_triangles);
            }
        }

        //=====================================================================
        // Serial reference version of the wall-distance calculation.
        //
        // This routine is retained to allow direct comparison with the OpenMP
        // version during debugging.  It performs exactly the same geometric
        // calculation but uses one thread only.
        //=====================================================================
        void WallDistance::computeSerial(CBSStateSI& s)
        {
            const std::vector<WallTriangle> wall_triangles =
                collectPhysicalWallTriangles(s);

            check_wall_triangle_list(wall_triangles, "WallDistance::computeSerial");

            s.wall_distance.fill(std::numeric_limits<Real>::max());

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                compute_one_node_distance(s, ip, wall_triangles);
            }
        }
    }
}
