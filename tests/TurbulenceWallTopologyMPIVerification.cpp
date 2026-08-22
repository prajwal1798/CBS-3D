//=============================================================================
// Two-rank regression for CHT interface reconstruction across an MPI cut.
//=============================================================================

#include "cbs/turbulence/TurbulenceWallTopology.hpp"

#ifdef CBS3D_USE_MPI
#include <mpi.h>
#endif

#include <array>
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

        s.ippn1(1, 1) = 2;
        s.ippn1(1, 2) = 3;
        s.ippn1(1, 3) = 4;

        s.ippn1(2, 1) = 1;
        s.ippn1(2, 2) = 4;
        s.ippn1(2, 3) = 3;

        s.ippn1(3, 1) = 1;
        s.ippn1(3, 2) = 2;
        s.ippn1(3, 3) = 4;

        s.ippn1(4, 1) = 1;
        s.ippn1(4, 2) = 3;
        s.ippn1(4, 3) = 2;
    }
}

int main(int argc, char** argv)
{
#ifndef CBS3D_USE_MPI
    (void)argc;
    (void)argv;
    std::printf("FAIL: MPI topology test built without CBS3D_USE_MPI\n");
    return 1;
#else
    MPI_Init(&argc, &argv);

    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int local_failure = 0;

    try
    {
        if (size != 2)
        {
            throw std::runtime_error(
                "TurbulenceWallTopologyMPIVerification requires exactly 2 ranks");
        }

        CBSStateSI s;
        s.mpi_enabled = true;
        s.mpi_rank = rank;
        s.mpi_size = size;

        s.cfg.npoin = 4;
        s.cfg.nelem = 1;
        s.cfg.nboun = 0;
        s.cfg.ndim = 3;
        s.cfg.nep = 4;
        s.cfg.nsid = 4;
        s.cfg.nsidp = 3;
        s.cfg.bsid = 6;
        s.cfg.nsidpe = 5;

        s.intma.resize(4, 1);
        s.iside.resize(6, 0);
        s.coord.resize(3, 4);
        s.mat_elem.resize(1);
        s.node_material_mask.resize(4);
        s.sa_active_node.resize(4);
        s.sa_wall_node.resize(4);
        s.sa_inlet_node.resize(4);
        s.wall_distance.resize(4);

        set_standard_tet_face_map(s);

        // Each rank owns one tetrahedron. Local nodes 2,3,4 represent the same
        // global face {10,11,12}. Rank 0 owns the fluid tetrahedron and rank 1
        // owns the solid tetrahedron, so the interface exists only after the
        // partition-face records are matched communicator-wide.
        s.intma(1, 1) = 1;
        s.intma(2, 1) = 2;
        s.intma(3, 1) = 3;
        s.intma(4, 1) = 4;
        s.mat_elem(1) = rank == 0 ? 0 : 1;

        const std::array<std::array<Real, 3>, 3> shared =
        {{
            {{1.0, 0.0, 0.0}},
            {{0.0, 1.0, 0.0}},
            {{0.0, 0.0, 1.0}}
        }};

        // Unique node lies on opposite sides of the shared face.
        const std::array<Real, 3> unique =
            rank == 0
                ? std::array<Real, 3>{{0.0, 0.0, 0.0}}
                : std::array<Real, 3>{{1.0, 1.0, 1.0}};

        for (Int dim = 1; dim <= 3; ++dim)
        {
            s.coord(dim, 1) = unique[static_cast<std::size_t>(dim - 1)];
            s.coord(dim, 2) = shared[0][static_cast<std::size_t>(dim - 1)];
            s.coord(dim, 3) = shared[1][static_cast<std::size_t>(dim - 1)];
            s.coord(dim, 4) = shared[2][static_cast<std::size_t>(dim - 1)];
        }

        s.local_to_global_node.assign(5U, 0);
        s.local_to_global_node[1] = rank == 0 ? 1 : 2;
        s.local_to_global_node[2] = 10;
        s.local_to_global_node[3] = 11;
        s.local_to_global_node[4] = 12;

        // The shared-node material masks represent the already-reconciled state
        // produced by Preprocess::buildMaterialNodeMasks().
        const Int interface_mask =
            CBSStateSI::node_touches_fluid |
            CBSStateSI::node_touches_solid;

        s.node_material_mask(1) =
            rank == 0
                ? CBSStateSI::node_touches_fluid
                : CBSStateSI::node_touches_solid;
        s.node_material_mask(2) = interface_mask;
        s.node_material_mask(3) = interface_mask;
        s.node_material_mask(4) = interface_mask;

        // Active flags are also represented after their owner/ghost OR. The
        // solid rank can hold interface-node copies that are active because the
        // same global nodes touch the fluid tetrahedron on rank 0.
        s.sa_active_node(1) = rank == 0 ? 1 : 0;
        s.sa_active_node(2) = 1;
        s.sa_active_node(3) = 1;
        s.sa_active_node(4) = 1;
        s.sa_wall_node.fill(0);
        s.sa_inlet_node.fill(0);

        TurbulenceWallTopology::reconcileWallNodeClassification(s);

        const std::vector<WallTriangle> walls =
            TurbulenceWallTopology::collectFluidWallTriangles(s);

        if (walls.size() != 1U)
        {
            std::printf(
                "rank %d FAIL: expected one reconstructed cross-rank interface, found %zu\n",
                rank,
                walls.size());
            local_failure = 1;
        }

        if (s.sa_wall_node(2) != 1 ||
            s.sa_wall_node(3) != 1 ||
            s.sa_wall_node(4) != 1)
        {
            std::printf(
                "rank %d FAIL: interface nodes were not classified as SA walls\n",
                rank);
            local_failure = 1;
        }
    }
    catch (const std::exception& error)
    {
        std::printf("rank %d FAIL: %s\n", rank, error.what());
        local_failure = 1;
    }

    int global_failure = 0;
    MPI_Allreduce(
        &local_failure,
        &global_failure,
        1,
        MPI_INT,
        MPI_MAX,
        MPI_COMM_WORLD);

    if (rank == 0 && global_failure == 0)
    {
        std::printf("PASS: MPI CHT turbulence-wall topology\n");
    }

    MPI_Finalize();
    return global_failure;
#endif
}
