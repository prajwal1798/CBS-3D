//=============================================================================
// CBS3D++_SI
//
// Verification and performance test for the Spalart-Allmaras wall-distance
// module.
//
// The test builds a flat-plate style structured node distribution with
// wall-normal grading, extracts the y = 0 wall as a triangulated surface, and
// then checks the BVH search against the exhaustive point-to-triangle search.
//
// Two properties are asserted:
//
//     1. Exactness.  The BVH distance must equal the brute-force distance for
//        every node to within floating-point round-off.  The BVH is a pruning
//        structure, not an approximation, so any difference above round-off is a
//        defect.
//
//     2. Geometric consistency.  For this configuration the analytic distance to
//        the plate is known away from the plate edges, so the interior nodes are
//        additionally checked against y.
//
// The measured timings are reported for both searches.
//
// Build:
//     g++ -std=c++17 -O2 -fopenmp -DCBS3D_USE_OPENMP -Iinclude \
//         tests/WallDistanceVerification.cpp src/turbulence/WallDistance.cpp \
//         -o build/wall_distance_verification
//=============================================================================

#include "cbs/turbulence/WallDistance.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace
{
    using cbs::Int;
    using cbs::Real;
    using cbs::Size;

    Real seconds_since(const std::chrono::steady_clock::time_point& start)
    {
        const std::chrono::duration<Real> elapsed =
            std::chrono::steady_clock::now() - start;

        return elapsed.count();
    }

    struct PlateMesh
    {
        Int nx = 0;
        Int ny = 0;
        Int nz = 0;
        Real length = 1.0;
        Real height = 0.1;
        Real width = 0.05;
        Real first_cell = 1.0e-5;
    };

    // Geometric wall-normal grading, similar to a real flat-plate boundary-layer
    // mesh: the first cell is very thin and the spacing expands away from the
    // plate.
    std::vector<Real> graded_coordinates(Int count, Real height, Real first_cell)
    {
        std::vector<Real> y(static_cast<Size>(count), 0.0);

        if (count < 2)
        {
            return y;
        }

        // Bisect for the expansion ratio giving the requested total height.
        const Real n = static_cast<Real>(count - 1);

        const auto total_height = [&](Real r)
        {
            if (std::fabs(r - 1.0) < 1.0e-12)
            {
                return first_cell * n;
            }

            return first_cell * (std::pow(r, n) - 1.0) / (r - 1.0);
        };

        Real low = 1.0;
        Real high = 2.0;

        while (total_height(high) < height && high < 100.0)
        {
            high *= 1.5;
        }

        for (Int iteration = 0; iteration < 200; ++iteration)
        {
            const Real mid = 0.5 * (low + high);

            if (total_height(mid) < height)
            {
                low = mid;
            }
            else
            {
                high = mid;
            }
        }

        const Real ratio = 0.5 * (low + high);

        Real spacing = first_cell;

        for (Int j = 1; j < count; ++j)
        {
            y[static_cast<Size>(j)] = y[static_cast<Size>(j - 1)] + spacing;
            spacing *= ratio;
        }

        return y;
    }

    void build_state(cbs::CBSStateSI& s, const PlateMesh& mesh)
    {
        const Int npoin = mesh.nx * mesh.ny * mesh.nz;

        // Two triangles per structured quad on the y = 0 plate.
        const Int nboun = 2 * (mesh.nx - 1) * (mesh.nz - 1);

        s.initialise_local_topology();
        s.set_problem_sizes(1, npoin, nboun, 0);

        const std::vector<Real> y = graded_coordinates(
            mesh.ny,
            mesh.height,
            mesh.first_cell);

        // One-based node numbering, i fastest.
        const auto node_id = [&mesh](Int i, Int j, Int k)
        {
            return 1 + i + mesh.nx * (j + mesh.ny * k);
        };

        for (Int k = 0; k < mesh.nz; ++k)
        {
            for (Int j = 0; j < mesh.ny; ++j)
            {
                for (Int i = 0; i < mesh.nx; ++i)
                {
                    const Int ip = node_id(i, j, k);

                    s.coord(1, ip) = mesh.length
                        * static_cast<Real>(i) / static_cast<Real>(mesh.nx - 1);

                    s.coord(2, ip) = y[static_cast<Size>(j)];

                    s.coord(3, ip) = mesh.width
                        * static_cast<Real>(k) / static_cast<Real>(mesh.nz - 1);
                }
            }
        }

        Int ib = 0;

        for (Int k = 0; k < mesh.nz - 1; ++k)
        {
            for (Int i = 0; i < mesh.nx - 1; ++i)
            {
                const Int n00 = node_id(i, 0, k);
                const Int n10 = node_id(i + 1, 0, k);
                const Int n01 = node_id(i, 0, k + 1);
                const Int n11 = node_id(i + 1, 0, k + 1);

                ++ib;
                s.iside(1, ib) = n00;
                s.iside(2, ib) = n10;
                s.iside(3, ib) = n11;
                s.iside(s.cfg.bsid, ib) = s.cfg.bc_noslip_adiabatic_wall;

                ++ib;
                s.iside(1, ib) = n00;
                s.iside(2, ib) = n11;
                s.iside(3, ib) = n01;
                s.iside(s.cfg.bsid, ib) = s.cfg.bc_noslip_adiabatic_wall;
            }
        }

        // Mark the plate nodes as SA wall nodes and everything else as active.
        s.sa_active_node.fill(1);
        s.sa_wall_node.fill(0);

        for (Int k = 0; k < mesh.nz; ++k)
        {
            for (Int i = 0; i < mesh.nx; ++i)
            {
                s.sa_wall_node(node_id(i, 0, k)) = 1;
            }
        }

        s.cfg.sa_min_wall_distance = 1.0e-14;
    }
}


int main(int argc, char** argv)
{
    PlateMesh mesh;

    mesh.nx = argc > 1 ? std::atoi(argv[1]) : 120;
    mesh.ny = argc > 2 ? std::atoi(argv[2]) : 80;
    mesh.nz = argc > 3 ? std::atoi(argv[3]) : 40;

    cbs::CBSStateSI state;
    build_state(state, mesh);

    const Int npoin = state.cfg.npoin;
    const Int nboun = state.cfg.nboun;

    std::printf("Flat-plate wall-distance verification\n");
    std::printf("  nodes           %d\n", npoin);
    std::printf("  wall triangles  %d\n", nboun);
    std::printf("  brute-force pair tests %.3e\n",
        static_cast<Real>(npoin) * static_cast<Real>(nboun));

    // Reference: exhaustive search.
    const std::chrono::steady_clock::time_point brute_start =
        std::chrono::steady_clock::now();

    cbs::turbulence::WallDistance::computeSerial(state);

    const Real brute_seconds = seconds_since(brute_start);

    std::vector<Real> reference(static_cast<Size>(npoin) + 1U, 0.0);

    for (Int ip = 1; ip <= npoin; ++ip)
    {
        reference[static_cast<Size>(ip)] = state.wall_distance(ip);
    }

    state.wall_distance.fill(-1.0);

    // Fast path: BVH plus Morton-ordered queries.
    cbs::turbulence::WallDistanceStats stats;

    const std::chrono::steady_clock::time_point bvh_start =
        std::chrono::steady_clock::now();

    cbs::turbulence::WallDistance::compute(state, stats);

    const Real bvh_seconds = seconds_since(bvh_start);

    // Exactness check.
    Real max_absolute = 0.0;
    Real max_relative = 0.0;
    Int worst_node = 0;

    for (Int ip = 1; ip <= npoin; ++ip)
    {
        const Real expected = reference[static_cast<Size>(ip)];
        const Real actual = state.wall_distance(ip);
        const Real absolute = std::fabs(actual - expected);

        if (absolute > max_absolute)
        {
            max_absolute = absolute;
            worst_node = ip;
        }

        if (expected > 0.0)
        {
            max_relative = std::max(max_relative, absolute / expected);
        }
    }

    // Analytic check on nodes away from the plate edges, where the nearest point
    // of an infinite plate is directly below and d = y.
    const Real edge_margin = 0.02;
    Real max_analytic = 0.0;

    for (Int ip = 1; ip <= npoin; ++ip)
    {
        const Real x = state.coord(1, ip);
        const Real y = state.coord(2, ip);
        const Real z = state.coord(3, ip);

        const bool interior =
            x > edge_margin && x < mesh.length - edge_margin &&
            z > edge_margin * 0.5 && z < mesh.width - edge_margin * 0.5;

        if (!interior || state.sa_wall_node(ip) != 0)
        {
            continue;
        }

        max_analytic = std::max(
            max_analytic,
            std::fabs(state.wall_distance(ip) - y));
    }

    std::printf("\nTiming\n");
    std::printf("  brute force        %8.3f s\n", brute_seconds);
    std::printf("  BVH total          %8.3f s\n", bvh_seconds);
    std::printf("    gather           %8.3f s\n", stats.gather_seconds);
    std::printf("    build            %8.3f s\n", stats.build_seconds);
    std::printf("    query            %8.3f s\n", stats.query_seconds);
    std::printf("  speedup            %8.1fx\n", brute_seconds / bvh_seconds);

    std::printf("\nBVH\n");
    std::printf("  searched triangles %zu\n", stats.searched_wall_triangles);
    std::printf("  nodes              %d\n", stats.bvh_nodes);
    std::printf("  depth              %d\n", stats.bvh_depth);
    std::printf("  queried nodes      %d\n", stats.queried_nodes);
    std::printf("  min distance       %.6e\n", stats.min_distance);
    std::printf("  max distance       %.6e\n", stats.max_distance);

    std::printf("\nAccuracy\n");
    std::printf("  max |BVH - brute|      %.3e (node %d)\n",
        max_absolute, worst_node);
    std::printf("  max relative error     %.3e\n", max_relative);
    std::printf("  max |d - y| interior   %.3e\n", max_analytic);

    const Real tolerance = 1.0e-12;
    bool failed = false;

    if (max_relative > tolerance)
    {
        std::printf("\nFAIL: BVH result differs from the exhaustive search\n");
        failed = true;
    }

    if (max_analytic > 1.0e-12)
    {
        std::printf("\nFAIL: interior distance does not match the analytic value\n");
        failed = true;
    }

    if (!failed)
    {
        std::printf("\nPASS\n");
    }

    return failed ? 1 : 0;
}
