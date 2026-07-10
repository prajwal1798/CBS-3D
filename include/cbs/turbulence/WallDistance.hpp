#pragma once

//=============================================================================
// CBS3D++_SI
//
// Wall-distance computation for the Spalart-Allmaras model.  The required
// quantity is the true minimum Euclidean distance from each active fluid node to
// the nearest physical no-slip wall triangle.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"

#include <array>
#include <vector>

namespace cbs::turbulence
{
    struct WallTriangle
    {
        std::array<Real, 3> a{};
        std::array<Real, 3> b{};
        std::array<Real, 3> c{};
    };

    class WallDistance
    {
    public:
        static std::vector<WallTriangle> collectPhysicalWallTriangles(const CBSStateSI& s);
        static Real pointTriangleDistance(const std::array<Real, 3>& p, const WallTriangle& tri);
        static void computeSerial(CBSStateSI& s);
    };
}
