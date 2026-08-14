#include "cbs/turbulence/WallDistance.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#ifdef CBS3D_USE_OPENMP
#include <omp.h>
#endif

#ifdef CBS3D_USE_MPI
#include <mpi.h>
#endif

namespace cbs
{
    namespace turbulence
    {
        namespace
        {
            const Real huge_real = 1.0e300;

            Real seconds_since(
                const std::chrono::steady_clock::time_point& start)
            {
                const std::chrono::duration<Real> elapsed =
                    std::chrono::steady_clock::now() - start;

                return elapsed.count();
            }

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

            Real dot3(const Real a[3], const Real b[3])
            {
                return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
            }

            //-----------------------------------------------------------------
            // Squared distance from a point to an axis-aligned box.  Returns
            // zero when the point is inside the box.  This is the pruning test
            // used during the BVH descent.
            //-----------------------------------------------------------------
            Real point_box_distance_squared(
                const std::array<Real, 3>& p,
                const Real bmin[3],
                const Real bmax[3])
            {
                Real total = 0.0;

                for (Int k = 0; k < 3; ++k)
                {
                    const Real value = p[static_cast<Size>(k)];

                    if (value < bmin[k])
                    {
                        const Real delta = bmin[k] - value;
                        total += delta * delta;
                    }
                    else if (value > bmax[k])
                    {
                        const Real delta = value - bmax[k];
                        total += delta * delta;
                    }
                }

                return total;
            }

            void expand_box(Real bmin[3], Real bmax[3], const std::array<Real, 3>& p)
            {
                for (Int k = 0; k < 3; ++k)
                {
                    const Size index = static_cast<Size>(k);

                    bmin[k] = std::min(bmin[k], p[index]);
                    bmax[k] = std::max(bmax[k], p[index]);
                }
            }

            void reset_box(Real bmin[3], Real bmax[3])
            {
                for (Int k = 0; k < 3; ++k)
                {
                    bmin[k] = huge_real;
                    bmax[k] = -huge_real;
                }
            }

            std::array<Real, 3> triangle_centroid(const WallTriangle& tri)
            {
                std::array<Real, 3> c;

                for (Size k = 0; k < 3; ++k)
                {
                    c[k] = (tri.a[k] + tri.b[k] + tri.c[k]) / 3.0;
                }

                return c;
            }

            //-----------------------------------------------------------------
            // Interleaves the low 21 bits of three integers into a 63-bit
            // Morton key.  Used to order query points along a Z-curve so that
            // consecutive queries touch the same BVH nodes.
            //-----------------------------------------------------------------
            std::uint64_t split_by_3(std::uint64_t value)
            {
                value &= 0x1fffffULL;
                value = (value | (value << 32)) & 0x1f00000000ffffULL;
                value = (value | (value << 16)) & 0x1f0000ff0000ffULL;
                value = (value | (value << 8)) & 0x100f00f00f00f00fULL;
                value = (value | (value << 4)) & 0x10c30c30c30c30c3ULL;
                value = (value | (value << 2)) & 0x1249249249249249ULL;

                return value;
            }

            std::uint64_t morton_key(
                const std::array<Real, 3>& p,
                const Real bmin[3],
                const Real inverse_extent[3])
            {
                std::uint64_t part[3];

                for (Int k = 0; k < 3; ++k)
                {
                    Real normalised =
                        (p[static_cast<Size>(k)] - bmin[k]) * inverse_extent[k];

                    if (!(normalised > 0.0))
                    {
                        normalised = 0.0;
                    }

                    if (normalised > 1.0)
                    {
                        normalised = 1.0;
                    }

                    const Real scaled = normalised * 2097151.0;
                    part[k] = static_cast<std::uint64_t>(scaled);
                }

                return split_by_3(part[0])
                    | (split_by_3(part[1]) << 1)
                    | (split_by_3(part[2]) << 2);
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
        }

        //=====================================================================
        // Exact squared distance from a point to a triangle.
        //
        // The procedure identifies which Voronoi region of the triangle
        // contains the closest point:
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
        // lines, which matters for unstructured tetrahedral meshes.
        //
        // The final barycentric branch divides by (va + vb + vc), which is
        // proportional to twice the squared triangle area.  A sliver or
        // duplicated wall face makes that quantity zero, so the denominator is
        // tested.  A degenerate triangle is handled by returning the smallest of
        // its three edge distances, which is the correct limit.
        //=====================================================================
        Real WallDistance::pointTriangleDistanceSquared(
            const std::array<Real, 3>& p,
            const WallTriangle& tri)
        {
            Real ab[3];
            Real ac[3];
            Real ap[3];

            for (Size k = 0; k < 3; ++k)
            {
                ab[k] = tri.b[k] - tri.a[k];
                ac[k] = tri.c[k] - tri.a[k];
                ap[k] = p[k] - tri.a[k];
            }

            const Real d1 = dot3(ab, ap);
            const Real d2 = dot3(ac, ap);

            if (d1 <= 0.0 && d2 <= 0.0)
            {
                return dot3(ap, ap);
            }

            Real bp[3];

            for (Size k = 0; k < 3; ++k)
            {
                bp[k] = p[k] - tri.b[k];
            }

            const Real d3 = dot3(ab, bp);
            const Real d4 = dot3(ac, bp);

            if (d3 >= 0.0 && d4 <= d3)
            {
                return dot3(bp, bp);
            }

            const Real vc = d1 * d4 - d3 * d2;

            if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
            {
                const Real denominator = d1 - d3;
                const Real v = denominator != 0.0 ? d1 / denominator : 0.0;
                Real diff[3];

                for (Size k = 0; k < 3; ++k)
                {
                    diff[k] = ap[k] - v * ab[k];
                }

                return dot3(diff, diff);
            }

            Real cp[3];

            for (Size k = 0; k < 3; ++k)
            {
                cp[k] = p[k] - tri.c[k];
            }

            const Real d5 = dot3(ab, cp);
            const Real d6 = dot3(ac, cp);

            if (d6 >= 0.0 && d5 <= d6)
            {
                return dot3(cp, cp);
            }

            const Real vb = d5 * d2 - d1 * d6;

            if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
            {
                const Real denominator = d2 - d6;
                const Real w = denominator != 0.0 ? d2 / denominator : 0.0;
                Real diff[3];

                for (Size k = 0; k < 3; ++k)
                {
                    diff[k] = ap[k] - w * ac[k];
                }

                return dot3(diff, diff);
            }

            const Real va = d3 * d6 - d5 * d4;

            if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
            {
                const Real numerator = d4 - d3;
                const Real denominator = numerator + (d5 - d6);
                const Real w = denominator != 0.0 ? numerator / denominator : 0.0;
                Real diff[3];

                for (Size k = 0; k < 3; ++k)
                {
                    diff[k] = bp[k] - w * (tri.c[k] - tri.b[k]);
                }

                return dot3(diff, diff);
            }

            const Real sum = va + vb + vc;

            // Degenerate face: the barycentric interior branch is not defined.
            // Fall back to the minimum of the three edge distances.
            if (!(std::abs(sum) > 0.0) || !std::isfinite(sum))
            {
                Real best = dot3(ap, ap);
                best = std::min(best, dot3(bp, bp));
                best = std::min(best, dot3(cp, cp));

                return best;
            }

            const Real inverse_sum = 1.0 / sum;
            const Real v = vb * inverse_sum;
            const Real w = vc * inverse_sum;

            Real diff[3];

            for (Size k = 0; k < 3; ++k)
            {
                diff[k] = ap[k] - (v * ab[k] + w * ac[k]);
            }

            return dot3(diff, diff);
        }


        Real WallDistance::pointTriangleDistance(
            const std::array<Real, 3>& p,
            const WallTriangle& tri)
        {
            return std::sqrt(pointTriangleDistanceSquared(p, tri));
        }


        //=====================================================================
        // BVH construction.
        //
        // Median splitting on the longest axis of the centroid bounding box is
        // used.  std::nth_element performs an O(n) partition per level, giving
        // O(n log n) overall, and produces a balanced tree without the cost or
        // the implementation risk of a full surface-area-heuristic sweep.  For a
        // wall surface of a few tens of thousands of triangles the build takes a
        // few milliseconds; the queries, not the build, are the bottleneck.
        //=====================================================================
        void WallTriangleBVH::build(
            std::vector<WallTriangle>& triangles,
            Int leaf_size)
        {
            nodes_.clear();
            triangles_.clear();
            depth_ = 0;

            if (triangles.empty())
            {
                return;
            }

            if (leaf_size < 1)
            {
                leaf_size = 1;
            }

            triangles_.swap(triangles);

            const Int count = static_cast<Int>(triangles_.size());

            std::vector<Int> index(static_cast<Size>(count));

            for (Int i = 0; i < count; ++i)
            {
                index[static_cast<Size>(i)] = i;
            }

            // A median-split tree with leaf_size >= 1 contains at most
            // 2 * count - 1 nodes.  Reserving up front keeps the recursive
            // push_back calls free of reallocation.
            nodes_.reserve(static_cast<Size>(2 * count));

            buildRecursive(index, 0, count, leaf_size, 1);

            // Reorder triangle storage so that every leaf reads one contiguous
            // block of memory.  index now holds the leaf-ordered permutation and
            // each leaf covers the half-open range [first_triangle, +count).
            std::vector<WallTriangle> reordered;
            reordered.resize(static_cast<Size>(count));

            for (Int i = 0; i < count; ++i)
            {
                reordered[static_cast<Size>(i)] =
                    triangles_[static_cast<Size>(index[static_cast<Size>(i)])];
            }

            triangles_.swap(reordered);
        }


        Int WallTriangleBVH::buildRecursive(
            std::vector<Int>& index,
            Int begin,
            Int end,
            Int leaf_size,
            Int level)
        {
            const Int node_index = static_cast<Int>(nodes_.size());

            nodes_.push_back(Node());

            depth_ = std::max(depth_, level);

            Real bmin[3];
            Real bmax[3];
            reset_box(bmin, bmax);

            Real cmin[3];
            Real cmax[3];
            reset_box(cmin, cmax);

            for (Int i = begin; i < end; ++i)
            {
                const WallTriangle& tri =
                    triangles_[static_cast<Size>(index[static_cast<Size>(i)])];

                expand_box(bmin, bmax, tri.a);
                expand_box(bmin, bmax, tri.b);
                expand_box(bmin, bmax, tri.c);

                expand_box(cmin, cmax, triangle_centroid(tri));
            }

            for (Int k = 0; k < 3; ++k)
            {
                nodes_[static_cast<Size>(node_index)].bmin[k] = bmin[k];
                nodes_[static_cast<Size>(node_index)].bmax[k] = bmax[k];
            }

            const Int count = end - begin;

            if (count <= leaf_size)
            {
                nodes_[static_cast<Size>(node_index)].first_triangle = begin;
                nodes_[static_cast<Size>(node_index)].triangle_count = count;
                nodes_[static_cast<Size>(node_index)].left_child = -1;
                nodes_[static_cast<Size>(node_index)].right_child = -1;

                return node_index;
            }

            Int axis = 0;
            Real best_extent = cmax[0] - cmin[0];

            for (Int k = 1; k < 3; ++k)
            {
                const Real extent = cmax[k] - cmin[k];

                if (extent > best_extent)
                {
                    best_extent = extent;
                    axis = k;
                }
            }

            const Int middle = begin + count / 2;

            // Coincident centroids give a zero-extent centroid box on every
            // axis.  A positional median split is still well defined, so the
            // tree stays balanced instead of degenerating into a linked list.
            const std::vector<WallTriangle>& triangle_storage = triangles_;

            std::nth_element(
                index.begin() + begin,
                index.begin() + middle,
                index.begin() + end,
                [&triangle_storage, axis](Int lhs, Int rhs)
                {
                    const Real lhs_value = triangle_centroid(
                        triangle_storage[static_cast<Size>(lhs)])
                            [static_cast<Size>(axis)];

                    const Real rhs_value = triangle_centroid(
                        triangle_storage[static_cast<Size>(rhs)])
                            [static_cast<Size>(axis)];

                    return lhs_value < rhs_value;
                });

            const Int left_index =
                buildRecursive(index, begin, middle, leaf_size, level + 1);

            const Int right_index =
                buildRecursive(index, middle, end, leaf_size, level + 1);

            nodes_[static_cast<Size>(node_index)].left_child = left_index;
            nodes_[static_cast<Size>(node_index)].right_child = right_index;
            nodes_[static_cast<Size>(node_index)].first_triangle = -1;
            nodes_[static_cast<Size>(node_index)].triangle_count = 0;

            return node_index;
        }


        bool WallTriangleBVH::empty() const
        {
            return triangles_.empty();
        }


        Size WallTriangleBVH::triangleCount() const
        {
            return triangles_.size();
        }


        Int WallTriangleBVH::nodeCount() const
        {
            return static_cast<Int>(nodes_.size());
        }


        Int WallTriangleBVH::depth() const
        {
            return depth_;
        }


        //=====================================================================
        // Nearest-triangle query.
        //
        // Depth-first descent with an explicit stack.  A child is pushed only if
        // its bounding box is closer than the best squared distance found so
        // far, and the nearer child is pushed last so that it is popped first
        // and tightens the bound as early as possible.
        //
        // upper_bound_squared is a pruning hint only.  Provided it is greater
        // than or equal to the true squared distance, the returned value is the
        // exact nearest squared distance: any triangle achieving the true
        // minimum has squared distance strictly below the seed, so neither it
        // nor any box containing it is ever pruned.
        //
        // The traversal stack is a fixed local array.  A median-split tree over
        // n triangles has depth ceil(log2(n)) + 1, so 128 entries covers any
        // wall surface that fits in memory.  The bound is checked rather than
        // assumed.
        //=====================================================================
        Real WallTriangleBVH::nearestDistanceSquared(
            const std::array<Real, 3>& p,
            Real upper_bound_squared) const
        {
            if (nodes_.empty())
            {
                return huge_real;
            }

            Real best = upper_bound_squared;

            if (!(best > 0.0) || !std::isfinite(best))
            {
                best = huge_real;
            }

            const Int stack_capacity = 128;
            Int stack[stack_capacity];
            Int stack_size = 0;

            stack[stack_size++] = 0;

            while (stack_size > 0)
            {
                const Int node_index = stack[--stack_size];
                const Node& node = nodes_[static_cast<Size>(node_index)];

                if (point_box_distance_squared(p, node.bmin, node.bmax) >= best)
                {
                    continue;
                }

                if (node.triangle_count > 0)
                {
                    const Int first = node.first_triangle;
                    const Int last = first + node.triangle_count;

                    for (Int i = first; i < last; ++i)
                    {
                        const Real distance =
                            WallDistance::pointTriangleDistanceSquared(
                                p,
                                triangles_[static_cast<Size>(i)]);

                        if (distance < best)
                        {
                            best = distance;
                        }
                    }

                    continue;
                }

                const Int left_child = node.left_child;
                const Int right_child = node.right_child;

                const Real left_distance = point_box_distance_squared(
                    p,
                    nodes_[static_cast<Size>(left_child)].bmin,
                    nodes_[static_cast<Size>(left_child)].bmax);

                const Real right_distance = point_box_distance_squared(
                    p,
                    nodes_[static_cast<Size>(right_child)].bmin,
                    nodes_[static_cast<Size>(right_child)].bmax);

                if (stack_size + 2 > stack_capacity)
                {
                    throw std::runtime_error(
                        "WallTriangleBVH::nearestDistanceSquared"
                        " - traversal stack overflow");
                }

                if (left_distance <= right_distance)
                {
                    if (right_distance < best)
                    {
                        stack[stack_size++] = right_child;
                    }

                    if (left_distance < best)
                    {
                        stack[stack_size++] = left_child;
                    }
                }
                else
                {
                    if (left_distance < best)
                    {
                        stack[stack_size++] = left_child;
                    }

                    if (right_distance < best)
                    {
                        stack[stack_size++] = right_child;
                    }
                }
            }

            return best;
        }


        //=====================================================================
        // Collects physical wall triangles owned by this rank.
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
            triangles.reserve(static_cast<Size>(s.cfg.nboun));

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
        // Collects the complete physical wall surface visible to this rank.
        //
        // Under domain decomposition the nearest wall triangle to a local node
        // frequently belongs to another rank, and a rank lying entirely in the
        // freestream owns no wall faces at all.  Computing the distance from
        // locally owned faces alone is therefore not an approximation, it is
        // wrong: it silently inflates d wherever the true nearest wall is off
        // rank, which suppresses the SA destruction term and corrupts the
        // boundary layer.
        //
        // The wall surface is a two-dimensional subset of a three-dimensional
        // mesh, so its size grows like npoin^(2/3).  Gathering it in full is
        // therefore cheap relative to the mesh itself, and it removes all
        // communication from the query phase: after the gather, every rank
        // performs a purely local BVH search.
        //
        // Duplicate faces are harmless for correctness, because the minimum over
        // a multiset equals the minimum over the underlying set.  Deduplication
        // is applied only when global node numbering is available, and purely as
        // a cost reduction.
        //=====================================================================
        std::vector<WallTriangle> WallDistance::collectGlobalWallTriangles(
            const CBSStateSI& s,
            WallDistanceStats& stats)
        {
            const std::chrono::steady_clock::time_point gather_start =
                std::chrono::steady_clock::now();

            std::vector<WallTriangle> local = collectPhysicalWallTriangles(s);

            stats.local_wall_triangles = local.size();

#ifdef CBS3D_USE_MPI
            if (s.mpi_enabled && s.mpi_size > 1)
            {
                const bool have_global_ids =
                    s.local_to_global_node.size()
                        >= static_cast<Size>(s.cfg.npoin) + 1U;

                // Pack the local wall faces as nine coordinates each, and, when
                // global numbering exists, the sorted global node triple used as
                // the deduplication key.
                std::vector<double> local_coordinates;
                std::vector<long long> local_keys;

                local_coordinates.reserve(local.size() * 9U);

                if (have_global_ids)
                {
                    local_keys.reserve(local.size() * 3U);
                }

                Size face_counter = 0;

                for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
                {
                    const Int bc = s.iside(s.cfg.bsid, ib);

                    if (!is_wall_bc(s, bc))
                    {
                        continue;
                    }

                    const WallTriangle& tri = local[face_counter];
                    ++face_counter;

                    for (Size k = 0; k < 3; ++k)
                    {
                        local_coordinates.push_back(tri.a[k]);
                    }

                    for (Size k = 0; k < 3; ++k)
                    {
                        local_coordinates.push_back(tri.b[k]);
                    }

                    for (Size k = 0; k < 3; ++k)
                    {
                        local_coordinates.push_back(tri.c[k]);
                    }

                    if (have_global_ids)
                    {
                        long long key[3];

                        for (Int corner = 0; corner < 3; ++corner)
                        {
                            const Int local_node = s.iside(corner + 1, ib);

                            key[corner] = static_cast<long long>(
                                s.local_to_global_node[
                                    static_cast<Size>(local_node)]);
                        }

                        std::sort(key, key + 3);

                        local_keys.push_back(key[0]);
                        local_keys.push_back(key[1]);
                        local_keys.push_back(key[2]);
                    }
                }

                const Int rank_count = s.mpi_size;

                std::vector<int> coordinate_counts(
                    static_cast<Size>(rank_count),
                    0);

                int local_coordinate_count =
                    static_cast<int>(local_coordinates.size());

                MPI_Allgather(
                    &local_coordinate_count,
                    1,
                    MPI_INT,
                    coordinate_counts.data(),
                    1,
                    MPI_INT,
                    MPI_COMM_WORLD);

                std::vector<int> coordinate_offsets(
                    static_cast<Size>(rank_count),
                    0);

                long long total_coordinates = 0;

                for (Int r = 0; r < rank_count; ++r)
                {
                    coordinate_offsets[static_cast<Size>(r)] =
                        static_cast<int>(total_coordinates);

                    total_coordinates +=
                        coordinate_counts[static_cast<Size>(r)];
                }

                if (total_coordinates > static_cast<long long>(
                        std::numeric_limits<int>::max()))
                {
                    throw std::runtime_error(
                        "WallDistance::collectGlobalWallTriangles - gathered wall"
                        " surface exceeds the MPI_Allgatherv element limit;"
                        " enable a pruned wall-distance exchange for this mesh");
                }

                std::vector<double> all_coordinates(
                    static_cast<Size>(total_coordinates));

                MPI_Allgatherv(
                    local_coordinates.data(),
                    local_coordinate_count,
                    MPI_DOUBLE,
                    all_coordinates.data(),
                    coordinate_counts.data(),
                    coordinate_offsets.data(),
                    MPI_DOUBLE,
                    MPI_COMM_WORLD);

                const Size gathered_faces =
                    static_cast<Size>(total_coordinates) / 9U;

                std::vector<long long> all_keys;

                if (have_global_ids)
                {
                    std::vector<int> key_counts(
                        static_cast<Size>(rank_count),
                        0);

                    std::vector<int> key_offsets(
                        static_cast<Size>(rank_count),
                        0);

                    for (Int r = 0; r < rank_count; ++r)
                    {
                        key_counts[static_cast<Size>(r)] =
                            coordinate_counts[static_cast<Size>(r)] / 3;
                    }

                    int running = 0;

                    for (Int r = 0; r < rank_count; ++r)
                    {
                        key_offsets[static_cast<Size>(r)] = running;
                        running += key_counts[static_cast<Size>(r)];
                    }

                    all_keys.resize(static_cast<Size>(running));

                    MPI_Allgatherv(
                        local_keys.data(),
                        static_cast<int>(local_keys.size()),
                        MPI_LONG_LONG,
                        all_keys.data(),
                        key_counts.data(),
                        key_offsets.data(),
                        MPI_LONG_LONG,
                        MPI_COMM_WORLD);
                }

                std::vector<WallTriangle> gathered;
                gathered.reserve(gathered_faces);

                std::vector<std::array<long long, 3> > seen;

                if (have_global_ids)
                {
                    seen.reserve(gathered_faces);

                    for (Size f = 0; f < gathered_faces; ++f)
                    {
                        std::array<long long, 3> key;

                        key[0] = all_keys[f * 3U + 0U];
                        key[1] = all_keys[f * 3U + 1U];
                        key[2] = all_keys[f * 3U + 2U];

                        seen.push_back(key);
                    }
                }

                std::vector<Size> keep;
                keep.reserve(gathered_faces);

                if (have_global_ids)
                {
                    std::vector<Size> order(gathered_faces);

                    for (Size f = 0; f < gathered_faces; ++f)
                    {
                        order[f] = f;
                    }

                    std::sort(
                        order.begin(),
                        order.end(),
                        [&seen](Size lhs, Size rhs)
                        {
                            return seen[lhs] < seen[rhs];
                        });

                    for (Size i = 0; i < order.size(); ++i)
                    {
                        if (i == 0 || seen[order[i]] != seen[order[i - 1U]])
                        {
                            keep.push_back(order[i]);
                        }
                    }
                }
                else
                {
                    for (Size f = 0; f < gathered_faces; ++f)
                    {
                        keep.push_back(f);
                    }
                }

                for (Size i = 0; i < keep.size(); ++i)
                {
                    const Size base = keep[i] * 9U;

                    WallTriangle tri;

                    for (Size k = 0; k < 3; ++k)
                    {
                        tri.a[k] = all_coordinates[base + 0U + k];
                        tri.b[k] = all_coordinates[base + 3U + k];
                        tri.c[k] = all_coordinates[base + 6U + k];
                    }

                    gathered.push_back(tri);
                }

                stats.global_wall_triangles = gathered.size();
                stats.gather_seconds = seconds_since(gather_start);

                return gathered;
            }
#endif

            stats.global_wall_triangles = local.size();
            stats.gather_seconds = seconds_since(gather_start);

            return local;
        }


        namespace
        {
            //-----------------------------------------------------------------
            // Discards gathered wall triangles that provably cannot be nearest
            // to any local query point.
            //
            // Let B be the axis-aligned box containing every local query point
            // and let v be any single wall-surface vertex.  For every p in B,
            //
            //     d(p) <= |p - v| <= R,      R = max over corners c of B of |c - v|
            //
            // so R is a valid upper bound on the wall distance of every local
            // query point.  A triangle T can only be nearest to some p in B if
            //
            //     minDist(B, AABB(T)) <= R
            //
            // and every triangle failing that test is discarded without changing
            // any computed distance.  Choosing v as the surface vertex nearest
            // the centre of B makes R as tight as this argument allows.
            //
            // On one rank B typically spans the whole domain and nothing is
            // discarded, which costs one linear scan.  Under decomposition B is
            // a small sub-block, so the far field of the wall surface is
            // discarded and both the BVH and the queries shrink accordingly.
            //-----------------------------------------------------------------
            void prune_unreachable_triangles(
                std::vector<WallTriangle>& triangles,
                const Real query_min[3],
                const Real query_max[3],
                bool have_query_box)
            {
                if (!have_query_box || triangles.empty())
                {
                    return;
                }

                Real centre[3];

                for (Int k = 0; k < 3; ++k)
                {
                    centre[k] = 0.5 * (query_min[k] + query_max[k]);
                }

                const std::array<Real, 3> centre_point =
                    { centre[0], centre[1], centre[2] };

                std::array<Real, 3> anchor = triangles[0].a;
                Real anchor_distance = huge_real;

                for (Size i = 0; i < triangles.size(); ++i)
                {
                    const std::array<Real, 3>* corner[3] =
                        { &triangles[i].a, &triangles[i].b, &triangles[i].c };

                    for (Int c = 0; c < 3; ++c)
                    {
                        Real total = 0.0;

                        for (Size k = 0; k < 3; ++k)
                        {
                            const Real delta =
                                (*corner[c])[k] - centre_point[k];

                            total += delta * delta;
                        }

                        if (total < anchor_distance)
                        {
                            anchor_distance = total;
                            anchor = *corner[c];
                        }
                    }
                }

                // R^2 = max over the eight corners of B of |corner - anchor|^2.
                Real radius_squared = 0.0;

                for (Int mask = 0; mask < 8; ++mask)
                {
                    Real total = 0.0;

                    for (Int k = 0; k < 3; ++k)
                    {
                        const Real value = (mask & (1 << k)) != 0
                            ? query_max[k]
                            : query_min[k];

                        const Real delta = value - anchor[static_cast<Size>(k)];

                        total += delta * delta;
                    }

                    radius_squared = std::max(radius_squared, total);
                }

                std::vector<WallTriangle> kept;
                kept.reserve(triangles.size());

                for (Size i = 0; i < triangles.size(); ++i)
                {
                    Real tmin[3];
                    Real tmax[3];
                    reset_box(tmin, tmax);

                    expand_box(tmin, tmax, triangles[i].a);
                    expand_box(tmin, tmax, triangles[i].b);
                    expand_box(tmin, tmax, triangles[i].c);

                    // Squared separation between the query box and the triangle
                    // box, zero when they overlap.
                    Real separation = 0.0;

                    for (Int k = 0; k < 3; ++k)
                    {
                        if (tmin[k] > query_max[k])
                        {
                            const Real delta = tmin[k] - query_max[k];
                            separation += delta * delta;
                        }
                        else if (query_min[k] > tmax[k])
                        {
                            const Real delta = query_min[k] - tmax[k];
                            separation += delta * delta;
                        }
                    }

                    if (separation <= radius_squared)
                    {
                        kept.push_back(triangles[i]);
                    }
                }

                triangles.swap(kept);
            }
        }


        //=====================================================================
        // Computes wall_distance(node) for every mesh node on this rank.
        //
        // The distance is evaluated for all npoin nodes rather than for SA-active
        // nodes only.  The extra work is negligible once the search is a BVH
        // query, and it guarantees that no consumer of the field can read an
        // uninitialised entry: element-averaged wall distances, y-plus
        // post-processing and VTU output are all well defined everywhere.
        //
        // Wall nodes are set to sa_min_wall_distance, a small positive floor that
        // keeps the SA destruction denominator finite while still representing
        // the wall boundary value.
        //
        // Query points are visited in Morton order.  Each OpenMP thread takes a
        // contiguous chunk of that curve, so its successive query points are
        // spatially adjacent: they descend the same BVH branches, reuse the same
        // cache lines, and the previous point's distance supplies a valid
        // starting bound through the 1-Lipschitz property
        //
        //     d(p) <= d(q) + |p - q|.
        //
        // The seed is inflated by a relative tolerance so that it stays a strict
        // upper bound under rounding, which keeps the returned distance exact.
        //=====================================================================
        void WallDistance::compute(CBSStateSI& s)
        {
            WallDistanceStats stats;
            compute(s, stats);
        }


        void WallDistance::compute(CBSStateSI& s, WallDistanceStats& stats)
        {
            stats = WallDistanceStats();

            std::vector<WallTriangle> wall_triangles =
                collectGlobalWallTriangles(s, stats);

            check_wall_triangle_list(wall_triangles, "WallDistance::compute");

            // Nodes needing a search: everything except the wall nodes, which
            // take the floor value directly.
            std::vector<Int> query_nodes;
            query_nodes.reserve(static_cast<Size>(s.cfg.npoin));

            Real query_min[3];
            Real query_max[3];
            reset_box(query_min, query_max);

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                if (s.sa_wall_node(ip) != 0)
                {
                    s.wall_distance(ip) = s.cfg.sa_min_wall_distance;
                    continue;
                }

                query_nodes.push_back(ip);
                expand_box(query_min, query_max, node_point(s, ip));
            }

            const bool have_query_box = !query_nodes.empty();

            prune_unreachable_triangles(
                wall_triangles,
                query_min,
                query_max,
                have_query_box);

            check_wall_triangle_list(
                wall_triangles,
                "WallDistance::compute after pruning");

            const std::chrono::steady_clock::time_point build_start =
                std::chrono::steady_clock::now();

            stats.searched_wall_triangles = wall_triangles.size();

            WallTriangleBVH bvh;
            bvh.build(wall_triangles, 8);

            stats.bvh_nodes = bvh.nodeCount();
            stats.bvh_depth = bvh.depth();
            stats.build_seconds = seconds_since(build_start);

            if (!have_query_box)
            {
                stats.queried_nodes = 0;
                return;
            }

            // Morton ordering of the query nodes.
            Real inverse_extent[3];

            for (Int k = 0; k < 3; ++k)
            {
                const Real extent = query_max[k] - query_min[k];

                inverse_extent[k] = extent > 0.0 ? 1.0 / extent : 0.0;
            }

            std::vector<std::uint64_t> keys(query_nodes.size());

            for (Size i = 0; i < query_nodes.size(); ++i)
            {
                keys[i] = morton_key(
                    node_point(s, query_nodes[i]),
                    query_min,
                    inverse_extent);
            }

            std::vector<Size> order(query_nodes.size());

            for (Size i = 0; i < order.size(); ++i)
            {
                order[i] = i;
            }

            std::sort(
                order.begin(),
                order.end(),
                [&keys](Size lhs, Size rhs)
                {
                    return keys[lhs] < keys[rhs];
                });

            const std::chrono::steady_clock::time_point query_start =
                std::chrono::steady_clock::now();

            const Real floor_distance = s.cfg.sa_min_wall_distance;
            const Int query_count = static_cast<Int>(order.size());

            // Relative inflation of the Lipschitz seed.  Without it a point
            // whose true distance exactly equals the bound could be pruned at
            // the root; with it the seed is strictly larger than the true
            // distance, so the nearest triangle is always visited.
            const Real seed_tolerance = 1.0 + 1.0e-10;

#ifdef CBS3D_USE_OPENMP
#pragma omp parallel
            {
                // Per-thread Lipschitz state.  Each thread seeds only from its
                // own previous query point, so the bound is always valid for the
                // point being evaluated and no thread reads another's state.
                bool have_previous = false;
                std::array<Real, 3> previous_point = { 0.0, 0.0, 0.0 };
                Real previous_distance = 0.0;

#pragma omp for schedule(static)
                for (Int i = 0; i < query_count; ++i)
                {
                    const Int ip = query_nodes[order[static_cast<Size>(i)]];
                    const std::array<Real, 3> p = node_point(s, ip);

                    Real seed = huge_real;

                    if (have_previous)
                    {
                        Real shift = 0.0;

                        for (Size k = 0; k < 3; ++k)
                        {
                            const Real delta = p[k] - previous_point[k];
                            shift += delta * delta;
                        }

                        const Real bound =
                            (previous_distance + std::sqrt(shift))
                                * seed_tolerance;

                        seed = bound * bound;
                    }

                    const Real distance_squared =
                        bvh.nearestDistanceSquared(p, seed);

                    const Real distance = std::sqrt(distance_squared);

                    s.wall_distance(ip) = std::max(distance, floor_distance);

                    have_previous = true;
                    previous_point = p;
                    previous_distance = distance;
                }
            }
#else
            {
                bool have_previous = false;
                std::array<Real, 3> previous_point = { 0.0, 0.0, 0.0 };
                Real previous_distance = 0.0;

                for (Int i = 0; i < query_count; ++i)
                {
                    const Int ip = query_nodes[order[static_cast<Size>(i)]];
                    const std::array<Real, 3> p = node_point(s, ip);

                    Real seed = huge_real;

                    if (have_previous)
                    {
                        Real shift = 0.0;

                        for (Size k = 0; k < 3; ++k)
                        {
                            const Real delta = p[k] - previous_point[k];
                            shift += delta * delta;
                        }

                        const Real bound =
                            (previous_distance + std::sqrt(shift))
                                * seed_tolerance;

                        seed = bound * bound;
                    }

                    const Real distance_squared =
                        bvh.nearestDistanceSquared(p, seed);

                    const Real distance = std::sqrt(distance_squared);

                    s.wall_distance(ip) = std::max(distance, floor_distance);

                    have_previous = true;
                    previous_point = p;
                    previous_distance = distance;
                }
            }
#endif

            stats.query_seconds = seconds_since(query_start);
            stats.queried_nodes = query_count;

            Real minimum = huge_real;
            Real maximum = 0.0;

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                minimum = std::min(minimum, s.wall_distance(ip));
                maximum = std::max(maximum, s.wall_distance(ip));
            }

            stats.min_distance = minimum;
            stats.max_distance = maximum;
        }


        //=====================================================================
        // Brute-force reference version.
        //
        // Retained so that the BVH result can be verified node by node against
        // an exhaustive search.  This is O(npoin * nwall) and must not be used on
        // production meshes.
        //
        // It deliberately uses only locally collected wall triangles, matching
        // the serial case for which it is a valid reference.
        //=====================================================================
        void WallDistance::computeSerial(CBSStateSI& s)
        {
            const std::vector<WallTriangle> wall_triangles =
                collectPhysicalWallTriangles(s);

            check_wall_triangle_list(
                wall_triangles,
                "WallDistance::computeSerial");

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                if (s.sa_wall_node(ip) != 0)
                {
                    s.wall_distance(ip) = s.cfg.sa_min_wall_distance;
                    continue;
                }

                const std::array<Real, 3> p = node_point(s, ip);
                Real best = huge_real;

                for (Size i = 0; i < wall_triangles.size(); ++i)
                {
                    const Real distance = pointTriangleDistanceSquared(
                        p,
                        wall_triangles[i]);

                    best = std::min(best, distance);
                }

                s.wall_distance(ip) = std::max(
                    std::sqrt(best),
                    s.cfg.sa_min_wall_distance);
            }
        }
    }
}
