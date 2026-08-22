//=============================================================================
// Serial regression for the material-aware SA turbulence-wall topology.
//=============================================================================

#include "cbs/turbulence/TurbulenceWallTopology.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <exception>
#include <vector>

namespace
{
    using cbs::CBSStateSI;
    using cbs::Int;
    using cbs::Real;
    using cbs::turbulence::TurbulenceWallTopology;
    using cbs::turbulence::WallTriangle;

    void set_standard_tet_face_map(CBSStateSI& s)
    {
        s.ippn1.resize(4, 3);

        // Face opposite local node 1.
        s.ippn1(1, 1) = 2;
        s.ippn1(1, 2) = 3;
        s.ippn1(1, 3) = 4;

        // Face opposite local node 2.
        s.ippn1(2, 1) = 1;
        s.ippn1(2, 2) = 4;
        s.ippn1(2, 3) = 3;

        // Face opposite local node 3.
        s.ippn1(3, 1) = 1;
        s.ippn1(3, 2) = 2;
        s.ippn1(3, 3) = 4;

        // Face opposite local node 4.
        s.ippn1(4, 1) = 1;
        s.ippn1(4, 2) = 3;
        s.ippn1(4, 3) = 2;
    }

    bool contains_point(
        const WallTriangle& triangle,
        const std::array<Real, 3>& point)
    {
        const std::array<std::array<Real, 3>, 3> vertices =
        {
            triangle.a,
            triangle.b,
            triangle.c
        };

        for (const auto& vertex : vertices)
        {
            Real squared = 0.0;

            for (std::size_t k = 0; k < 3; ++k)
            {
                const Real delta = vertex[k] - point[k];
                squared += delta * delta;
            }

            if (squared < 1.0e-28)
            {
                return true;
            }
        }

        return false;
    }

    bool nearly_equal(
        const Real a,
        const Real b,
        const Real relative_tolerance,
        const Real absolute_tolerance = 0.0)
    {
        return std::fabs(a - b) <=
            absolute_tolerance +
            relative_tolerance * std::max(std::fabs(a), std::fabs(b));
    }
}

int main()
{
    try
    {
        CBSStateSI s;

        s.cfg.npoin = 6;
        s.cfg.nelem = 2;
        s.cfg.nboun = 2;
        s.cfg.ndim = 3;
        s.cfg.nep = 4;
        s.cfg.nsid = 4;
        s.cfg.nsidp = 3;
        s.cfg.bsid = 6;
        s.cfg.nsidpe = 5;

        s.intma.resize(4, 2);
        s.iside.resize(6, 2);
        s.coord.resize(3, 6);
        s.mat_elem.resize(2);
        s.node_material_mask.resize(6);
        s.sa_active_node.resize(6);
        s.sa_wall_node.resize(6);
        s.sa_inlet_node.resize(6);
        s.wall_distance.resize(6);

        set_standard_tet_face_map(s);

        // Fluid tetrahedron: nodes 1,2,3,4.
        s.intma(1, 1) = 1;
        s.intma(2, 1) = 2;
        s.intma(3, 1) = 3;
        s.intma(4, 1) = 4;

        // Solid tetrahedron: nodes 5,2,3,4. The shared internal face is
        // exactly {2,3,4}; there is deliberately no BC 901 record for it.
        s.intma(1, 2) = 5;
        s.intma(2, 2) = 2;
        s.intma(3, 2) = 3;
        s.intma(4, 2) = 4;

        s.mat_elem(1) = 0;
        s.mat_elem(2) = 1;

        // Node 6 is an interior fluid query point used only to verify the
        // corrected nearest-wall distance. It is not part of the synthetic
        // topology because the two tetrahedra above already define the exact
        // wall surface under test.
        const std::array<std::array<Real, 3>, 6> coordinates =
        {{
            {{0.0, 0.0, 0.0}},
            {{1.0, 0.0, 0.0}},
            {{0.0, 1.0, 0.0}},
            {{0.0, 0.0, 1.0}},
            {{1.0, 1.0, 1.0}},
            {{0.2, 0.2, 0.2}}
        }};

        for (Int ip = 1; ip <= 6; ++ip)
        {
            for (Int dim = 1; dim <= 3; ++dim)
            {
                s.coord(dim, ip) =
                    coordinates[static_cast<std::size_t>(ip - 1)]
                               [static_cast<std::size_t>(dim - 1)];
            }
        }

        s.node_material_mask(1) = CBSStateSI::node_touches_fluid;
        s.node_material_mask(2) =
            CBSStateSI::node_touches_fluid | CBSStateSI::node_touches_solid;
        s.node_material_mask(3) = s.node_material_mask(2);
        s.node_material_mask(4) = s.node_material_mask(2);
        s.node_material_mask(5) = CBSStateSI::node_touches_solid;
        s.node_material_mask(6) = CBSStateSI::node_touches_fluid;

        // Fluid tetrahedron nodes plus the interior query point are SA active.
        for (Int ip = 1; ip <= 6; ++ip)
        {
            s.sa_active_node(ip) = (ip <= 4 || ip == 6) ? 1 : 0;
            s.sa_wall_node(ip) = 0;
            s.sa_inlet_node(ip) = 0;
        }

        // Explicit fluid BC 530: face {1,4,3}, parent fluid tetrahedron 1.
        s.iside(1, 1) = 1;
        s.iside(2, 1) = 4;
        s.iside(3, 1) = 3;
        s.iside(4, 1) = 2;
        s.iside(5, 1) = 1;
        s.iside(6, 1) = s.cfg.bc_noslip_adiabatic_wall;

        // Explicit solid BC 532: face {5,4,3}, parent solid tetrahedron 2.
        // This is analogous to a plasma-facing heat-flux surface and must not
        // become a helium SA wall merely because its BC number is 532.
        s.iside(1, 2) = 5;
        s.iside(2, 2) = 4;
        s.iside(3, 2) = 3;
        s.iside(4, 2) = 2;
        s.iside(5, 2) = 2;
        s.iside(6, 2) = s.cfg.bc_noslip_heatflux_wall;

        TurbulenceWallTopology::reconcileWallNodeClassification(s);

        // Node 1 is on the explicit fluid wall. Nodes 2,3,4 are reconstructed
        // material-interface wall nodes. Node 5 is solid-only and node 6 is an
        // interior fluid point; neither may be classified as an SA wall.
        if (s.sa_wall_node(1) != 1 ||
            s.sa_wall_node(2) != 1 ||
            s.sa_wall_node(3) != 1 ||
            s.sa_wall_node(4) != 1 ||
            s.sa_wall_node(5) != 0 ||
            s.sa_wall_node(6) != 0)
        {
            std::printf("FAIL: serial wall-node classification\n");
            return 1;
        }

        const std::vector<WallTriangle> walls =
            TurbulenceWallTopology::collectFluidWallTriangles(s);

        // Expected surface:
        //   1. explicit fluid BC 530;
        //   2. internal fluid-solid face {2,3,4} reconstructed from mat_elem.
        // The solid-parent BC 532 must be absent.
        if (walls.size() != 2U)
        {
            std::printf(
                "FAIL: expected 2 turbulence-wall triangles, found %zu\n",
                walls.size());
            return 1;
        }

        const std::array<Real, 3> solid_only_point = coordinates[4];

        for (const WallTriangle& wall : walls)
        {
            if (contains_point(wall, solid_only_point))
            {
                std::printf("FAIL: solid-parent BC 532 entered SA wall geometry\n");
                return 1;
            }
        }

        cbs::turbulence::WallDistanceStats stats;
        TurbulenceWallTopology::computeWallDistance(s, stats);

        if (stats.global_wall_triangles != 2U)
        {
            std::printf("FAIL: corrected wall-distance inventory count\n");
            return 1;
        }

        for (Int ip = 1; ip <= 4; ++ip)
        {
            if (s.wall_distance(ip) != 0.0)
            {
                std::printf("FAIL: wall node received non-zero wall distance\n");
                return 1;
            }
        }

        // The nearest accepted wall to (0.2,0.2,0.2) is the explicit fluid
        // triangle in plane x=0, so the exact Euclidean distance is 0.2. The
        // reconstructed fluid-solid interface lies on x+y+z=1 and is farther
        // away: 0.4/sqrt(3) approximately 0.23094.
        if (!nearly_equal(s.wall_distance(6), 0.2, 2.0e-13, 1.0e-14))
        {
            std::printf(
                "FAIL: interior nearest-wall distance=% .17e expected=2.0e-1\n",
                s.wall_distance(6));
            return 1;
        }

        std::printf("PASS: serial CHT turbulence-wall topology\n");
        return 0;
    }
    catch (const std::exception& error)
    {
        std::printf("FAIL: unexpected exception: %s\n", error.what());
        return 1;
    }
}
