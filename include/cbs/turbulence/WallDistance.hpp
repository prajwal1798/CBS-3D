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
// The first implementation is a serial reference implementation.  It loops over
// all active fluid nodes and all physical wall triangles.  This is intentionally
// simple and easy to verify before replacing it with an accelerated search tree.
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

            static void computeSerial(CBSStateSI& s);
        };
    }
}
