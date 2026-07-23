//=============================================================================
// CBS3D++_SI
//
// Distributed-memory execution of CBS Step 3: pressure-gradient velocity
// correction followed by a discrete continuity-residual measurement.
//=============================================================================

#include "cbs/solver/Solver.hpp"

#include "cbs/assembly/PressureAssembly.hpp"
#include "cbs/assembly/VelocityCorrectionAssembly.hpp"
#include "cbs/boundary/Boundary.hpp"
#include "cbs/parallel/HaloExchange.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#ifdef CBS3D_USE_MPI
#include <mpi.h>
#endif

namespace cbs
{
#ifdef CBS3D_USE_MPI
    namespace
    {
        void check_step3_mpi(
            const int error_code,
            const char* operation)
        {
            if (error_code != MPI_SUCCESS)
            {
                throw std::runtime_error(
                    std::string(
                        "Solver::runDistributedStep3 - MPI failure in ")
                    + operation);
            }
        }


        bool pressure_active(
            const CBSStateSI& s,
            const Int ip)
        {
            return
                (s.node_material_mask(ip) &
                 CBSStateSI::node_touches_fluid) != 0;
        }


        struct ContinuityMetrics
        {
            Real l2_squared = 0.0;
            Real max_abs = 0.0;
            long long count = 0;
        };


        ContinuityMetrics global_free_continuity_metrics(
            const CBSStateSI& s)
        {
            ContinuityMetrics local;

            for (const Int ip : s.owned_nodes)
            {
                if (!pressure_active(s, ip) ||
                    s.node_pressure_fixed(ip) != 0)
                {
                    continue;
                }

                const Real residual = s.rhs1(ip);

                if (!std::isfinite(residual))
                {
                    throw std::runtime_error(
                        "Solver::runDistributedStep3 found a non-finite "
                        "continuity residual");
                }

                local.l2_squared += residual * residual;
                local.max_abs =
                    std::max(local.max_abs, std::abs(residual));
                ++local.count;
            }

            ContinuityMetrics global;

            check_step3_mpi(
                MPI_Allreduce(
                    &local.l2_squared,
                    &global.l2_squared,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM,
                    MPI_COMM_WORLD),
                "MPI_Allreduce continuity L2");

            check_step3_mpi(
                MPI_Allreduce(
                    &local.max_abs,
                    &global.max_abs,
                    1,
                    MPI_DOUBLE,
                    MPI_MAX,
                    MPI_COMM_WORLD),
                "MPI_Allreduce continuity maximum");

            check_step3_mpi(
                MPI_Allreduce(
                    &local.count,
                    &global.count,
                    1,
                    MPI_LONG_LONG,
                    MPI_SUM,
                    MPI_COMM_WORLD),
                "MPI_Allreduce continuity row count");

            if (global.count <= 0)
            {
                throw std::runtime_error(
                    "Solver::runDistributedStep3 found no free pressure rows");
            }

            return global;
        }
    }
#endif


    //=========================================================================
    // Executes distributed CBS Steps 1, 2 and 3.
    //
    // Step 3 follows the standard owner/ghost finite-element sequence:
    //
    //     1. Broadcast owner pressure to ghost copies.
    //     2. Assemble pressure-gradient contributions on owned elements.
    //     3. Reverse-sum shared-node vectors to node owners.
    //     4. Correct velocity on owners only.
    //     5. Reapply strong velocity constraints on owners.
    //     6. Broadcast corrected owner velocity to ghosts.
    //
    // The weak-divergence residual is then reassembled with the corrected
    // velocity. Norms are evaluated only on unconstrained pressure rows because
    // prescribed-pressure equations replace the original continuity rows.
    //=========================================================================
    void Solver::runDistributedStep3()
    {
#if defined(CBS3D_USE_MPI) && defined(CBS3D_USE_PETSC)
        if (!s_.mpi_enabled)
        {
            throw std::runtime_error(
                "Solver::runDistributedStep3 requires more than one MPI rank");
        }

        // Step 2 leaves u*, pressure and the assembled Step-2 RHS available.
        runDistributedStep2();

        const ContinuityMetrics before =
            global_free_continuity_metrics(s_);

        // Every rank-owned tetrahedron must read current pressure values at all
        // four local nodes, including pressure values owned by neighbours.
        HaloExchange::broadcastOwnedToGhosts(
            s_.pres,
            s_.partition_metadata,
            MPI_COMM_WORLD);

        // Rank-local element integration followed by distributed nodal assembly.
        VelocityCorrectionAssembly::assembleStep3Rhs(s_);

        HaloExchange::sumGhostContributionsToOwners(
            s_.rhs,
            s_.partition_metadata,
            MPI_COMM_WORLD);

        Real local_correction_l2_squared = 0.0;
        Real local_correction_max = 0.0;

        for (const Int ip : s_.owned_nodes)
        {
            if (!pressure_active(s_, ip))
            {
                continue;
            }

            for (Int dim = 1;
                 dim <= s_.cfg.ndim;
                 ++dim)
            {
                const Real correction =
                    s_.rhs(dim, ip) * s_.elcoe2(ip);

                if (!std::isfinite(correction))
                {
                    throw std::runtime_error(
                        "Solver::runDistributedStep3 produced a non-finite "
                        "velocity correction");
                }

                s_.unkno(dim, ip) += correction;
                local_correction_l2_squared +=
                    correction * correction;
                local_correction_max =
                    std::max(
                        local_correction_max,
                        std::abs(correction));
            }
        }

        // The current blanket case contains no symmetry nodes. Stop explicitly
        // rather than applying an unreconciled rank-local symmetry projection.
        for (const Int ip : s_.owned_nodes)
        {
            if (s_.node_symmetry(ip) != 0)
            {
                throw std::runtime_error(
                    "Solver::runDistributedStep3 requires a persistent "
                    "distributed symmetry-normal operator");
            }
        }

        Boundary::applyOwnedVelocityConstraints(s_);

        HaloExchange::broadcastOwnedToGhosts(
            s_.unkno,
            s_.partition_metadata,
            MPI_COMM_WORLD);

        Real local_corrected_velocity_max = 0.0;

        for (const Int ip : s_.owned_nodes)
        {
            const Real u = s_.unkno(1, ip);
            const Real v = s_.unkno(2, ip);
            const Real w = s_.unkno(3, ip);

            if (!std::isfinite(u) ||
                !std::isfinite(v) ||
                !std::isfinite(w))
            {
                throw std::runtime_error(
                    "Solver::runDistributedStep3 produced non-finite "
                    "corrected velocity");
            }

            local_corrected_velocity_max =
                std::max(
                    local_corrected_velocity_max,
                    std::sqrt(u * u + v * v + w * w));
        }

        Real global_correction_l2_squared = 0.0;
        Real global_correction_max = 0.0;
        Real global_corrected_velocity_max = 0.0;

        check_step3_mpi(
            MPI_Allreduce(
                &local_correction_l2_squared,
                &global_correction_l2_squared,
                1,
                MPI_DOUBLE,
                MPI_SUM,
                MPI_COMM_WORLD),
            "MPI_Allreduce velocity-correction L2");

        check_step3_mpi(
            MPI_Allreduce(
                &local_correction_max,
                &global_correction_max,
                1,
                MPI_DOUBLE,
                MPI_MAX,
                MPI_COMM_WORLD),
            "MPI_Allreduce velocity-correction maximum");

        check_step3_mpi(
            MPI_Allreduce(
                &local_corrected_velocity_max,
                &global_corrected_velocity_max,
                1,
                MPI_DOUBLE,
                MPI_MAX,
                MPI_COMM_WORLD),
            "MPI_Allreduce corrected-velocity maximum");

        // Reassemble the same weak-divergence operator with u^(n+1).
        PressureAssembly::assembleStep2Rhs(s_);

        HaloExchange::sumGhostContributionsToOwners(
            s_.rhs1,
            s_.partition_metadata,
            MPI_COMM_WORLD);

        const ContinuityMetrics after =
            global_free_continuity_metrics(s_);

        const Real before_l2 =
            std::sqrt(before.l2_squared);

        const Real after_l2 =
            std::sqrt(after.l2_squared);

        const Real before_rms =
            before_l2 /
            std::sqrt(static_cast<Real>(before.count));

        const Real after_rms =
            after_l2 /
            std::sqrt(static_cast<Real>(after.count));

        const Real after_over_before =
            before_l2 > 1.0e-300
                ? after_l2 / before_l2
                : 0.0;

        if (!std::isfinite(after_over_before))
        {
            throw std::runtime_error(
                "Solver::runDistributedStep3 produced an invalid "
                "continuity-reduction factor");
        }

        if (s_.mpi_rank == 0)
        {
            std::cout
                << "============================================================\n"
                << "CBS3D DISTRIBUTED STEP 3 COMPLETE\n"
                << "============================================================\n"
                << "MPI ranks                    : " << s_.mpi_size << "\n"
                << "free pressure rows           : " << before.count << "\n"
                << "continuity L2 before         : " << before_l2 << "\n"
                << "continuity L2 after          : " << after_l2 << "\n"
                << "continuity RMS before        : " << before_rms << "\n"
                << "continuity RMS after         : " << after_rms << "\n"
                << "continuity max before        : " << before.max_abs << "\n"
                << "continuity max after         : " << after.max_abs << "\n"
                << "after / before L2            : " << after_over_before << "\n"
                << "velocity-correction L2       : "
                << std::sqrt(global_correction_l2_squared) << "\n"
                << "maximum |velocity correction|: "
                << global_correction_max << "\n"
                << "maximum corrected velocity   : "
                << global_corrected_velocity_max << "\n"
                << "continuity status             : "
                << (after_l2 < before_l2 ? "REDUCED" : "NOT REDUCED")
                << "\n"
                << "next solver stage            : distributed energy/loop\n"
                << "============================================================\n";
        }
#else
        throw std::runtime_error(
            "Solver::runDistributedStep3 requires MPI and PETSc support");
#endif
    }
}
