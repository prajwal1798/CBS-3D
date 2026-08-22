//=============================================================================
// MPI regression for rank-local wall-model topology.
//
// A distributed flat plate does NOT require every rank to touch the plate.
// Rank 0 owns one valid BC530 TRI3 wall face; rank 1 owns no model-wall face.
// Both ranks must initialise the production wall treatment, agree that the
// global wall-face count is one, and preserve their correct local wall-node
// inventory without throwing.
//=============================================================================

#include "cbs/turbulence/WallModelCoupling.hpp"

#include <mpi.h>

#include <array>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <vector>

namespace
{
    using cbs::CBSStateSI;
    using cbs::Int;
    using cbs::Real;
    using cbs::turbulence::WallModelCoupling;

    CBSStateSI make_state(const int rank, const int size)
    {
        const Int nboun = rank == 0 ? 1 : 0;

        CBSStateSI s;
        s.initialise_local_topology();
        s.set_problem_sizes(1, 4, nboun, 0);

        s.cfg.turbulence_on = 1;
        s.cfg.turbulence_model = 0;
        s.cfg.dimensional_mode = 0;
        s.cfg.material_properties_enabled = 0;
        s.cfg.ani = 1.0e-5;

        s.mpi_rank = rank;
        s.mpi_size = size;
        s.mpi_enabled = size > 1;

        // No halo neighbours are required for this regression: the two local
        // tetrahedra are intentionally independent.  The wall-face count is a
        // communicator-wide property and is audited with MPI collectives.
        s.partition_metadata.neighbours.clear();

        s.local_to_global_node.assign(5, -1);
        for (Int ip = 1; ip <= 4; ++ip)
        {
            s.local_to_global_node[static_cast<std::size_t>(ip)] =
                4 * rank + ip;
        }

        s.intma(1, 1) = 1;
        s.intma(2, 1) = 2;
        s.intma(3, 1) = 3;
        s.intma(4, 1) = 4;
        s.mat_elem(1) = 0;

        // Unit-base tetrahedron with height 0.2, detJ = 0.2.
        const std::array<std::array<Real, 3>, 4> xyz =
        {{
            {{0.0, 0.0, 0.0}},
            {{1.0, 0.0, 0.0}},
            {{0.0, 1.0, 0.0}},
            {{0.0, 0.0, 0.2}}
        }};

        for (Int ip = 1; ip <= 4; ++ip)
        {
            for (Int dim = 1; dim <= 3; ++dim)
            {
                s.coord(dim, ip) =
                    xyz[static_cast<std::size_t>(ip - 1)]
                       [static_cast<std::size_t>(dim - 1)];
            }
        }

        s.detJ(1) = 0.2;
        s.fedge.fill(0);
        s.node_velocity_bc_type.fill(CBSStateSI::velocity_bc_free);
        s.unkno.fill(0.0);
        s.unkn1.fill(0.0);
        s.rhs.fill(0.0);

        if (rank == 0)
        {
            // Face 4 is opposite local node 4.  The fluid occupies z>0, so
            // the outward unit normal on z=0 is -ez and A=0.5.
            s.iside(1, 1) = 1;
            s.iside(2, 1) = 2;
            s.iside(3, 1) = 3;
            s.iside(s.cfg.nsidpl, 1) = 4;
            s.iside(s.cfg.nsidpe, 1) = 1;
            s.iside(s.cfg.bsid, 1) = s.cfg.bc_noslip_adiabatic_wall;

            s.face_norm(1, 1) = 0.0;
            s.face_norm(2, 1) = 0.0;
            s.face_norm(3, 1) = -0.5;
            s.face_norm(4, 1) = 0.5;
            s.fedge(4, 1) = 1;

            s.node_velocity_bc_type(1) = CBSStateSI::velocity_bc_noslip;
            s.node_velocity_bc_type(2) = CBSStateSI::velocity_bc_noslip;
            s.node_velocity_bc_type(3) = CBSStateSI::velocity_bc_noslip;
        }

        return s;
    }
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size != 2)
    {
        if (rank == 0)
        {
            std::fprintf(stderr, "FAIL: regression requires exactly two MPI ranks\n");
        }
        MPI_Finalize();
        return 2;
    }

    setenv("CBS3D_SA_WALL_TREATMENT", "1", 1);

    int local_ok = 1;

    try
    {
        CBSStateSI s = make_state(rank, size);

        const long long global_faces =
            WallModelCoupling::globalWallFaceCount(s);

        if (global_faces != 1)
        {
            std::fprintf(
                stderr,
                "rank %d FAIL: global wall-face count=%lld, expected 1\n",
                rank,
                global_faces);
            local_ok = 0;
        }

        const auto captured = WallModelCoupling::captureVelocity(s, false);

        if (rank == 0)
        {
            if (captured.nodes.size() != 3 ||
                !WallModelCoupling::isModelWallNode(s, 1) ||
                !WallModelCoupling::isModelWallNode(s, 2) ||
                !WallModelCoupling::isModelWallNode(s, 3) ||
                WallModelCoupling::isModelWallNode(s, 4))
            {
                std::fprintf(stderr, "rank 0 FAIL: incorrect local wall-node inventory\n");
                local_ok = 0;
            }
        }
        else
        {
            if (!captured.nodes.empty())
            {
                std::fprintf(
                    stderr,
                    "rank 1 FAIL: interior rank invented %zu wall nodes\n",
                    captured.nodes.size());
                local_ok = 0;
            }

            for (Int ip = 1; ip <= 4; ++ip)
            {
                if (WallModelCoupling::isModelWallNode(s, ip))
                {
                    std::fprintf(
                        stderr,
                        "rank 1 FAIL: local node %d incorrectly classified as wall\n",
                        ip);
                    local_ok = 0;
                }
            }
        }
    }
    catch (const std::exception& error)
    {
        std::fprintf(
            stderr,
            "rank %d FAIL: unexpected wall-model exception: %s\n",
            rank,
            error.what());
        local_ok = 0;
    }

    int global_ok = 0;
    MPI_Allreduce(
        &local_ok,
        &global_ok,
        1,
        MPI_INT,
        MPI_MIN,
        MPI_COMM_WORLD);

    if (rank == 0 && global_ok)
    {
        std::printf(
            "PASS: MPI wall-model permits ranks with zero local wall faces/nodes\n");
    }

    MPI_Finalize();
    return global_ok ? 0 : 1;
}
