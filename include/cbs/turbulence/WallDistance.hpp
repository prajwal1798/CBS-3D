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
// The production CBS solver uses unstructured tetrahedral meshes.  Therefore the
// wall distance is computed by a point-to-triangle search over the physical wall
// boundary faces.  The threaded routine is still deliberately simple: each
// active fluid node is independent, so OpenMP parallelises the outer node loop.
//
// This is a preprocessing operation.  The distance field is computed once after
// mesh and boundary preprocessing, stored in wall_distance(node), and then reused
// by every SA transport-equation assembly.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"

#include <array>
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

        class WallDistance
        {
        public:
            static std::vector<WallTriangle> collectPhysicalWallTriangles(
                const CBSStateSI& s);

            static Real pointTriangleDistance(
                const std::array<Real, 3>& p,
                const WallTriangle& tri);

            // Computes wall_distance(node) using OpenMP when the solver was
            // built with CBS3D_ENABLE_OPENMP=ON.  If OpenMP is not available,
            // the same routine automatically falls back to the serial loop.
            static void compute(CBSStateSI& s);

            // Serial reference version kept for debugging and verification.
            static void computeSerial(CBSStateSI& s);
        };
    }
}
