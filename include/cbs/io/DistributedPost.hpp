#pragma once

//=============================================================================
// CBS3D++_SI
//
// Distributed-memory residual and ParaView output.
//
// Every MPI rank writes one VTU piece containing its owned tetrahedra and
// rank-local node set. Rank zero writes the PVTU descriptor, the PVD time-series
// collection and the globally reduced residual CSV.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"
#include "cbs/core/Types.hpp"

#include <string>

namespace cbs
{
    class DistributedPost
    {
    public:
        // Derives the common distributed case name from a rank-local case path.
        [[nodiscard]] static std::string distributedCaseName(
            const std::string& rank_local_case_name);

        // Creates the distributed output directory, residual CSV and optional
        // iteration-zero PVTU/PVD solution.
        static void initialise(
            CBSStateSI& s,
            const std::string& rank_local_case_name);

        // Returns true when the current iteration satisfies an iteration-based
        // or physical-time-based distributed output trigger.
        [[nodiscard]] static bool shouldWriteSolution(
            CBSStateSI& s,
            Int iteration);

        // Writes one VTU piece per rank, then a rank-zero PVTU descriptor and
        // updated PVD collection.
        static void writeSolution(
            CBSStateSI& s,
            const std::string& rank_local_case_name,
            Int iteration);

        // Appends one globally reduced convergence row on rank zero.
        static void writeResidualRow(
            const CBSStateSI& s,
            const std::string& rank_local_case_name,
            Int iteration,
            Real continuity_rms,
            Real continuity_max,
            Real maximum_velocity,
            Real maximum_velocity_correction,
            Real iteration_wall_seconds);

        // Tests whether the requested iteration has already been recorded in
        // the distributed PVD history.
        [[nodiscard]] static bool solutionAlreadyWritten(
            const CBSStateSI& s,
            Int iteration);
    };
}
