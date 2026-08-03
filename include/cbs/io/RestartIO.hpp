#pragma once

//=============================================================================
// CBS3D++_SI
//
// Distributed restart/checkpoint support.
//
// The native restart format stores the complete rank-local solution state for
// the current MPI partition.  A one-time legacy importer is also provided for
// the existing ASCII VTU pieces written by DistributedPost.
//
// The first implementation intentionally requires the same mesh,
// local/global numbering, MPI process count and domain decomposition.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"
#include "cbs/core/Types.hpp"

#include <algorithm>
#include <string>

namespace cbs
{
    class RestartIO
    {
    public:
        struct LoadResult
        {
            bool loaded = false;
            bool imported_from_legacy_vtu = false;
            Int completed_iteration = 0;
            Real physical_time = 0.0;
            Real time_step = 0.0;
        };

        // Loads a restart when either RESTART_OPT=1 in the .par file or the
        // CBS3D_RESTART environment flag is enabled.
        //
        // Environment controls:
        //
        //   CBS3D_RESTART=1
        //   CBS3D_RESTART_FORMAT=native | legacy_vtu
        //   CBS3D_RESTART_ROOT=<checkpoint directory or legacy VTU directory>
        //   CBS3D_RESTART_ITERATION=<required for legacy_vtu>
        [[nodiscard]] static LoadResult loadIfRequested(
            CBSStateSI& s,
            const std::string& rank_local_case_name);

        // Writes one failure-safe native MPI checkpoint.  Each rank writes one
        // binary file; rank zero commits the manifest only after every rank
        // file has been written and renamed successfully.
        static void writeCheckpoint(
            const CBSStateSI& s,
            const std::string& rank_local_case_name,
            Int completed_iteration);

        // Returns CBS3D_CHECKPOINT_EVERY.  Zero disables periodic checkpoints.
        [[nodiscard]] static Int checkpointInterval();

        // Returns the configured checkpoint root.  When
        // CBS3D_CHECKPOINT_ROOT is not set, the default is:
        //
        //   output/<distributed-case>/restart
        [[nodiscard]] static std::string checkpointRoot(
            const std::string& rank_local_case_name);
    };
}
