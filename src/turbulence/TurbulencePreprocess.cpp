#include "cbs/turbulence/TurbulencePreprocess.hpp"

#include "cbs/assembly/SpalartAllmarasAssembly.hpp"
#include "cbs/boundary/TurbulenceBoundary.hpp"
#include "cbs/parallel/HaloExchange.hpp"
#include "cbs/turbulence/TurbulenceWallTopology.hpp"

#include <cctype>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace cbs
{
    namespace
    {
        bool legacy_vtu_restart_requested()
        {
            const char* value = std::getenv("CBS3D_RESTART_FORMAT");

            if (value == nullptr || value[0] == '\0')
            {
                return false;
            }

            std::string format(value);

            for (char& character : format)
            {
                character = static_cast<char>(
                    std::tolower(static_cast<unsigned char>(character)));
            }

            return
                format == "legacy_vtu" ||
                format == "vtu" ||
                format == "legacy";
        }
    }

    //=========================================================================
    // Prepares the Spalart-Allmaras turbulence model before time advancement.
    //=========================================================================
    void TurbulencePreprocess::prepareSpalartAllmaras(CBSStateSI& s)
    {
        SpalartAllmarasAssembly::resetEffectiveProperties(s);

        if (s.cfg.turbulence_on < 1)
        {
            return;
        }

        // Build the ordinary SA active/inlet flags first and reconcile them
        // across the partition. The material node mask was already reconciled
        // by preprocessing, so the authoritative wall topology can then replace
        // the old boundary-ID-only wall classification deterministically.
        TurbulenceBoundary::classifyNodes(s);
        TurbulenceBoundary::synchroniseClassification(s);
        turbulence::TurbulenceWallTopology::reconcileWallNodeClassification(s);
        TurbulenceBoundary::synchroniseClassification(s);

        turbulence::WallDistanceStats stats;
        turbulence::TurbulenceWallTopology::computeWallDistance(s, stats);

        if (s.mpi_rank == 0)
        {
            std::cout
                << "  Spalart-Allmaras wall distance\n"
                << "    local wall triangles     "
                << stats.local_wall_triangles << "\n"
                << "    global wall triangles    "
                << stats.global_wall_triangles << "\n"
                << "    searched wall triangles  "
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

        // Native restart safety.
        //
        // Native .cbsrst checkpoints contain nu_tilde and therefore preserve it.
        // The legacy VTU importer deliberately does not contain an SA history;
        // RestartIO sets nu_tilde to zero on that path, so treating every
        // istart>1 state as a native restart would silently start SA from an
        // invalid zero field. Distinguish the explicitly selected legacy
        // importer and initialise SA from the configured freestream value there.
        const bool loaded_restart_state = s.cfg.istart > 1;
        const bool legacy_vtu_restart =
            loaded_restart_state && legacy_vtu_restart_requested();
        const bool native_restart =
            loaded_restart_state && !legacy_vtu_restart;

        if (native_restart)
        {
            // Keep the checkpoint field and re-enforce only strong SA BCs.
            TurbulenceBoundary::applyWallValues(s);
            TurbulenceBoundary::applyInletValues(s);
        }
        else
        {
            // Fresh calculation or velocity/pressure-only legacy VTU import.
            // initialiseNuTilde currently re-runs the legacy local classifier,
            // so the authoritative topology is restored immediately below
            // before the field is exposed to the momentum equation.
            TurbulenceBoundary::initialiseNuTilde(s);
        }

        // Restore globally consistent active/inlet flags, then replace the
        // legacy wall flag with the material-aware topology again. This is
        // required after initialiseNuTilde(), which deliberately performs its
        // own local classification. Reapplying wall/inlet values guarantees
        // nu_tilde=0 on every reconstructed fluid-solid interface node.
        TurbulenceBoundary::synchroniseClassification(s);
        turbulence::TurbulenceWallTopology::reconcileWallNodeClassification(s);
        TurbulenceBoundary::synchroniseClassification(s);
        TurbulenceBoundary::applyWallValues(s);
        TurbulenceBoundary::applyInletValues(s);

#ifdef CBS3D_USE_MPI
        if (s.mpi_enabled && s.mpi_size > 1)
        {
            HaloExchange::broadcastOwnedToGhosts(
                s.nu_tilde,
                s.partition_metadata,
                MPI_COMM_WORLD);
        }
#endif

        // The first continuation residual must be measured against the restored
        // state, not against the original freestream initialisation.
        s.nu_tilde1 = s.nu_tilde;

        if (s.mpi_rank == 0)
        {
            if (native_restart)
            {
                std::cout
                    << "  Spalart-Allmaras: nu_tilde restored from native checkpoint\n";
            }
            else if (legacy_vtu_restart)
            {
                std::cout
                    << "  Spalart-Allmaras: legacy VTU imported; nu_tilde reinitialised from freestream\n";
            }
            else
            {
                std::cout
                    << "  Spalart-Allmaras: nu_tilde initialised from freestream\n";
            }
        }

        // Rebuild all derived turbulent/effective material properties from the
        // authoritative nu_tilde field before the first momentum timestep.
        SpalartAllmarasAssembly::updateEddyViscosity(s);
    }
}
