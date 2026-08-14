#pragma once

//=============================================================================
// CBS3D++_SI
//
// Wall-distance computation for the Spalart-Allmaras turbulence model.
//
// The SA destruction term depends on
//
//     d = distance to nearest no-slip wall
//
// This distance must be the true minimum Euclidean distance from a fluid node to
// a physical wall triangle.  It must not be approximated by nearest wall node,
// grid-line search, or artificial partition-boundary distance.
//
//-----------------------------------------------------------------------------
// Algorithm
//-----------------------------------------------------------------------------
//
// The previous implementation compared every fluid node against every wall
// triangle, giving O(npoin * nwall) work.  For a flat-plate mesh with 1e6 nodes
// and 2e4 wall triangles that is 2e10 point-triangle tests, which dominates
// preprocessing.  This module replaces that with a three-part scheme:
//
//     1. Bounding volume hierarchy (BVH).
//        The wall triangles are stored in an axis-aligned bounding-box tree
//        built by median splitting on the longest axis of the centroid box.
//        A nearest-triangle query descends the tree nearest-child-first and
//        prunes any subtree whose box is further away than the best distance
//        found so far.  Expected cost per query is O(log nwall).
//
//     2. Squared-distance arithmetic.
//        Every comparison is done on squared distances.  Exactly one square
//        root is taken per node instead of one per point-triangle test.
//
//     3. Morton-ordered queries with Lipschitz bound seeding.
//        Query points are visited along a Morton (Z-order) curve so that
//        consecutive queries are spatially adjacent and reuse the same BVH
//        nodes and cache lines.  Because the distance function is 1-Lipschitz,
//
//            d(p) <= d(q) + |p - q|
//
//        the distance of the previous point in the curve provides a valid
//        starting upper bound for the next point, which prunes the descent
//        immediately instead of starting from infinity.
//
//-----------------------------------------------------------------------------
// Parallelism
//-----------------------------------------------------------------------------
//
// OpenMP: the query loop is parallelised over the Morton-ordered node list with
// a static schedule.  Each thread owns a contiguous chunk of the curve, writes
// only to its own wall_distance(ip) entries, and maintains its own private
// Lipschitz seed taken from its own previous point.  The BVH is read-only.
//
// MPI: a node near a partition interface may have its nearest wall triangle on
// a neighbouring rank, and a rank in the freestream may hold no wall faces at
// all.  Computing the distance from local wall faces only is therefore wrong,
// not merely inexact.  The wall surface is a two-dimensional subset of a
// three-dimensional mesh, so gathering it in full is cheap: this module
// all-gathers the physical wall triangles, deduplicates them by sorted global
// node id, and then prunes the gathered set with a provably safe local-radius
// test before building the BVH.  See collectGlobalWallTriangles.
//
// This is a preprocessing operation.  The distance field is computed once after
// mesh and boundary preprocessing, stored in wall_distance(node), and then
// reused by every SA transport-equation assembly.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"

#include <array>
#include <cstdint>
#include <vector>

namespace cbs
{
    namespace turbulence
    {
        struct WallTriangle
        {
            std::array<Real, 3> a;
            std::array<Real, 3> b;
            std::array<Real, 3> c;
        };

        //=====================================================================
        // Statistics reported by the wall-distance preprocessing step.
        //
        // These are printed by TurbulencePreprocess so that the cost of the
        // wall-distance search is visible in the solver log, and so that the
        // effect of the BVH can be confirmed on a production mesh instead of
        // assumed.
        //=====================================================================
        struct WallDistanceStats
        {
            Size local_wall_triangles = 0;
            Size global_wall_triangles = 0;
            Size searched_wall_triangles = 0;
            Int bvh_nodes = 0;
            Int bvh_depth = 0;
            Int queried_nodes = 0;
            Real build_seconds = 0.0;
            Real gather_seconds = 0.0;
            Real query_seconds = 0.0;
            Real min_distance = 0.0;
            Real max_distance = 0.0;
        };

        //=====================================================================
        // Bounding volume hierarchy over a fixed set of wall triangles.
        //
        // The tree is immutable once built.  All query methods are const and
        // therefore safe to call concurrently from multiple OpenMP threads.
        //=====================================================================
        class WallTriangleBVH
        {
        public:
            // Takes ownership of the triangle list and reorders it so that the
            // triangles of each leaf are contiguous in memory.
            void build(std::vector<WallTriangle>& triangles, Int leaf_size = 8);

            bool empty() const;
            Size triangleCount() const;
            Int nodeCount() const;
            Int depth() const;

            // Returns the squared distance from p to the nearest triangle.
            //
            // upper_bound_squared is a pruning hint.  It must be greater than or
            // equal to the true squared distance; pass a large value when no
            // bound is known.  The returned value is always the exact nearest
            // squared distance regardless of the hint.
            Real nearestDistanceSquared(
                const std::array<Real, 3>& p,
                Real upper_bound_squared) const;

        private:
            struct Node
            {
                Real bmin[3];
                Real bmax[3];

                // A leaf has triangle_count > 0; first_triangle is then the
                // offset of its first triangle in the reordered triangle array
                // and both child indices are -1.
                //
                // An internal node has triangle_count == 0 and stores both
                // child indices explicitly.  The right child is not adjacent to
                // the left child in the depth-first build order, so it cannot be
                // derived and must be recorded.
                Int left_child = -1;
                Int right_child = -1;
                Int first_triangle = -1;
                Int triangle_count = 0;
            };

            Int buildRecursive(
                std::vector<Int>& index,
                Int begin,
                Int end,
                Int leaf_size,
                Int level);

            std::vector<Node> nodes_;
            std::vector<WallTriangle> triangles_;
            Int depth_ = 0;
        };

        class WallDistance
        {
        public:
            // Physical no-slip wall triangles owned by this rank.
            static std::vector<WallTriangle> collectPhysicalWallTriangles(
                const CBSStateSI& s);

            // Complete physical wall surface visible to this rank.
            //
            // In a serial run this is identical to collectPhysicalWallTriangles.
            // Under MPI the local lists are all-gathered and deduplicated so
            // that every rank can compute a true global nearest-wall distance.
            static std::vector<WallTriangle> collectGlobalWallTriangles(
                const CBSStateSI& s,
                WallDistanceStats& stats);

            // Exact squared distance from a point to one triangle.  Degenerate
            // triangles fall back to the nearest-edge result instead of
            // dividing by a zero barycentric denominator.
            static Real pointTriangleDistanceSquared(
                const std::array<Real, 3>& p,
                const WallTriangle& tri);

            // Retained for interface compatibility.  Returns the square root of
            // pointTriangleDistanceSquared.
            static Real pointTriangleDistance(
                const std::array<Real, 3>& p,
                const WallTriangle& tri);

            // Computes wall_distance(node) for every SA-active node using the
            // BVH search, Morton query ordering and OpenMP.
            static void compute(CBSStateSI& s);

            static void compute(CBSStateSI& s, WallDistanceStats& stats);

            // Brute-force O(npoin * nwall) reference version.  Retained for
            // verification of the BVH result and for regression testing only.
            // It is not intended for production meshes.
            static void computeSerial(CBSStateSI& s);
        };
    }
}
