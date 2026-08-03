//=============================================================================
// CBS3D++_SI
//
// Restart-aware wrapper around the validated distributed production loop.
//
// The numerical loop remains in SolverDistributedProductionLoop.cpp.  This
// translation unit compiles that implementation once with three controlled
// substitutions:
//
//   1. preprocessing loads a native or legacy-VTU restart after the mesh and
//      partition metadata are reconstructed;
//   2. output iteration numbers are offset by the completed restart iteration;
//   3. timestep logic receives the global continuation iteration.
//
// This keeps the validated CBS assembly sequence unchanged.
//=============================================================================

#include "cbs/io/DistributedPost.hpp"
#include "cbs/io/RestartIO.hpp"
#include "cbs/solver/Solver.hpp"
#include "cbs/timestep/TimeStep.hpp"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <string>

namespace cbs
{
    namespace
    {
        bool restart_environment_flag_enabled(const char* name)
        {
            const char* value = std::getenv(name);

            return value != nullptr &&
                   value[0] != '\0' &&
                   std::string(value) != "0";
        }


        class RestartTimeStep
        {
        public:
            static void computeTimeStep(
                CBSStateSI& s,
                const Int local_iteration)
            {
                const Int global_iteration =
                    s.cfg.iiter_total + local_iteration;

                TimeStep::computeTimeStep(s, global_iteration);
            }
        };


        class RestartDistributedPost
        {
        public:
            static void reset()
            {
                last_completed_iteration_ = 0;
                last_checkpoint_iteration_ = -1;
            }

            [[nodiscard]] static std::string distributedCaseName(
                const std::string& rank_local_case_name)
            {
                return DistributedPost::distributedCaseName(
                    rank_local_case_name);
            }

            static void initialise(
                CBSStateSI& s,
                const std::string& rank_local_case_name)
            {
                last_completed_iteration_ = s.cfg.iiter_total;
                last_checkpoint_iteration_ = -1;

                if (s.cfg.iiter_total <= 0)
                {
                    DistributedPost::initialise(
                        s,
                        rank_local_case_name);
                    return;
                }

                // Initialise a new continuation output directory and residual
                // CSV without creating a misleading iteration-zero field.
                const Int saved_vtu_output = s.cfg.vtu_output_enabled;
                s.cfg.vtu_output_enabled = 0;

                DistributedPost::initialise(
                    s,
                    rank_local_case_name);

                s.cfg.vtu_output_enabled = saved_vtu_output;

                // Record the imported/restarted state as the first field in the
                // new continuation PVD using its true global iteration number.
                if (saved_vtu_output > 0)
                {
                    DistributedPost::writeSolution(
                        s,
                        rank_local_case_name,
                        s.cfg.iiter_total);
                }
            }

            [[nodiscard]] static bool shouldWriteSolution(
                CBSStateSI& s,
                const Int local_iteration)
            {
                return DistributedPost::shouldWriteSolution(
                    s,
                    globalIteration(s, local_iteration));
            }

            static void writeSolution(
                CBSStateSI& s,
                const std::string& rank_local_case_name,
                const Int local_iteration)
            {
                DistributedPost::writeSolution(
                    s,
                    rank_local_case_name,
                    globalIteration(s, local_iteration));
            }

            static void writeResidualRow(
                const CBSStateSI& s,
                const std::string& rank_local_case_name,
                const Int local_iteration,
                const Real continuity_rms,
                const Real continuity_max,
                const Real maximum_velocity,
                const Real maximum_velocity_correction,
                const Real iteration_wall_seconds)
            {
                const Int global_iteration =
                    globalIteration(s, local_iteration);

                last_completed_iteration_ = global_iteration;

                DistributedPost::writeResidualRow(
                    s,
                    rank_local_case_name,
                    global_iteration,
                    continuity_rms,
                    continuity_max,
                    maximum_velocity,
                    maximum_velocity_correction,
                    iteration_wall_seconds);

                const Int checkpoint_interval =
                    RestartIO::checkpointInterval();

                if (checkpoint_interval > 0 &&
                    (global_iteration % checkpoint_interval) == 0)
                {
                    RestartIO::writeCheckpoint(
                        s,
                        rank_local_case_name,
                        global_iteration);

                    last_checkpoint_iteration_ = global_iteration;
                }
            }

            [[nodiscard]] static bool solutionAlreadyWritten(
                const CBSStateSI& s,
                const Int local_iteration)
            {
                return DistributedPost::solutionAlreadyWritten(
                    s,
                    globalIteration(s, local_iteration));
            }

            [[nodiscard]] static Int lastCompletedIteration()
            {
                return last_completed_iteration_;
            }

            [[nodiscard]] static Int lastCheckpointIteration()
            {
                return last_checkpoint_iteration_;
            }

            static void markCheckpoint(const Int iteration)
            {
                last_checkpoint_iteration_ = iteration;
            }

        private:
            [[nodiscard]] static Int globalIteration(
                const CBSStateSI& s,
                const Int local_iteration)
            {
                return s.cfg.iiter_total + local_iteration;
            }

            inline static Int last_completed_iteration_ = 0;
            inline static Int last_checkpoint_iteration_ = -1;
        };
    }
}

// Compile the validated production implementation under private method and
// helper names.  The public wrapper below performs the restart-specific work.
#define DistributedPost RestartDistributedPost
#define TimeStep RestartTimeStep
#define runDistributedProductionLoop runDistributedProductionLoopImpl
#define runDistributedPreprocessing runDistributedPreprocessingWithRestart
#include "SolverDistributedProductionLoop.cpp"
#undef runDistributedPreprocessing
#undef runDistributedProductionLoop
#undef TimeStep
#undef DistributedPost

namespace cbs
{
    void Solver::runDistributedPreprocessingWithRestart()
    {
        // Reconstruct the exact rank-local mesh, ownership maps, geometry,
        // material masks and boundary state before any solution is imported.
        runDistributedPreprocessing();

        const RestartIO::LoadResult restart =
            RestartIO::loadIfRequested(s_, case_name_);

        if (!restart.loaded)
        {
            return;
        }

        if (s_.mpi_rank == 0)
        {
            std::cout
                << "Restart loaded\n"
                << "  source      : "
                << (restart.imported_from_legacy_vtu
                    ? "legacy distributed VTU"
                    : "native CBS3D checkpoint")
                << "\n"
                << "  completed   : "
                << restart.completed_iteration << " iterations\n"
                << "  continuation: "
                << restart.completed_iteration + 1
                << " to "
                << restart.completed_iteration + s_.cfg.ntime
                << "\n";
        }

        if (restart.imported_from_legacy_vtu)
        {
            // Immediately convert the existing 100,000-iteration VTU state to
            // the native format.  Subsequent jobs must restart from this native
            // checkpoint rather than reparsing ParaView files.
            RestartIO::writeCheckpoint(
                s_,
                case_name_,
                restart.completed_iteration);

            RestartDistributedPost::markCheckpoint(
                restart.completed_iteration);
        }
    }


    void Solver::runDistributedProductionLoop()
    {
        RestartDistributedPost::reset();

        runDistributedProductionLoopImpl();

        const Int checkpoint_interval =
            RestartIO::checkpointInterval();

        const Int last_completed_iteration =
            RestartDistributedPost::lastCompletedIteration();

        if (checkpoint_interval > 0 &&
            last_completed_iteration > 0 &&
            last_completed_iteration !=
                RestartDistributedPost::lastCheckpointIteration())
        {
            // Preserve a final checkpoint even when a short validation job or
            // an early convergence stop does not coincide with the interval.
            RestartIO::writeCheckpoint(
                s_,
                case_name_,
                last_completed_iteration);

            RestartDistributedPost::markCheckpoint(
                last_completed_iteration);
        }

        if (s_.mpi_rank == 0 && s_.cfg.iiter_total > 0)
        {
            std::cout
                << "Restart note: console iteration counters above are local "
                << "to this submitted segment; residual and solution files "
                << "use global continuation iterations.\n";
        }
    }
}
