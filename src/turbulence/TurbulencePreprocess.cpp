#include "cbs/turbulence/TurbulencePreprocess.hpp"

#include "cbs/assembly/SpalartAllmarasAssembly.hpp"
#include "cbs/boundary/TurbulenceBoundary.hpp"
#include "cbs/parallel/HaloExchange.hpp"
#include "cbs/turbulence/WallDistance.hpp"

#include <iostream>
#include <stdexcept>

namespace cbs
{
    //=========================================================================
    // Prepares the Spalart-Allmaras turbulence model before time advancement.
    //
    // This routine performs only preprocessing and initialisation.  It does not
    // assemble the SA transport equation and it does not modify the CBS Step 1,
    // Step 2, Step 3 or Step 4 algorithms.
    //
    // Operations performed when turbulence_on is enabled:
    //
    //     1. classify SA-active, SA-wall and SA-inlet nodes;
    //     2. reduce that classification over partition interfaces;
    //     3. compute wall_distance(node) with the BVH nearest-wall search;
    //     4. initialise the transported SA working variable nu_tilde;
    //     5. compute the first eddy-viscosity and effective-property fields.
    //
    // The wall-distance search is MPI aware.  It gathers the complete physical
    // wall surface before searching, because under domain decomposition the
    // nearest wall triangle to a local node is frequently owned by another rank
    // and a rank lying in the freestream owns no wall faces at all.  See
    // WallDistance::collectGlobalWallTriangles.
    //=========================================================================
    void TurbulencePreprocess::prepareSpalartAllmaras(CBSStateSI& s)
    {
        SpalartAllmarasAssembly::resetEffectiveProperties(s);

        if (s.cfg.turbulence_on < 1)
        {
            return;
        }

        TurbulenceBoundary::classifyNodes(s);
        TurbulenceBoundary::synchroniseClassification(s);

        turbulence::WallDistanceStats stats;
        turbulence::WallDistance::compute(s, stats);

        if (s.mpi_rank == 0)
        {
            std::cout
                << "  Spalart-Allmaras wall distance\n"
                << "    local wall triangles     "
                << stats.local_wall_triangles << "\n"
                << "    global wall triangles    "
                << stats.global_wall_triangles << "\n"
                << "    searched after pruning   "
                << stats.searched_wall_triangles << "\n"
                << "    BVH nodes / depth        "
                << stats.bvh_nodes << " / " << stats.bvh_depth << "\n"
                << "    queried mesh nodes       "
                << stats.queried_nodes << "\n"
                << "    gather / build / query   "
                << stats.gather_seconds << " / "
                << stats.build_seconds << " / "
                << stats.query_seconds << " s\n"
                << "    min / max distance       "
                << stats.min_distance << " / "
                << stats.max_distance << "\n";
        }

        TurbulenceBoundary::initialiseNuTilde(s);

        // initialiseNuTilde re-runs the local classification, so the interface
        // reduction has to be repeated before the field is exchanged.
        TurbulenceBoundary::synchroniseClassification(s);

#ifdef CBS3D_USE_MPI
        if (s.mpi_enabled && s.mpi_size > 1)
        {
            HaloExchange::broadcastOwnedToGhosts(
                s.nu_tilde,
                s.partition_metadata,
                MPI_COMM_WORLD);
        }
#endif

        // updateEddyViscosity reduces its own nodal averages across interfaces.
        SpalartAllmarasAssembly::updateEddyViscosity(s);
    }
}
