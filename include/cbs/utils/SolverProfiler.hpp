#pragma once

//=============================================================================
// CBS3D++_SI
//
// Per-iteration wall-clock profiler for the main CBS solver stages.
//
// The profiler measures the elapsed time of selected parts of one CBS
// iteration:
//
//     time-step calculation
//     left-hand-side diagonal construction
//     CBS Steps 1 to 4
//     convergence evaluation
//     post-processing and solution output
//
// A timed scope is created in the solver as:
//
//     auto timer = profiler.time(section);
//
// The timer records the starting wall-clock time when it is created. When the
// program leaves the surrounding scope, its destructor calculates:
//
//     elapsed time = finish time - start time
//
// and adds the result to the selected section.
//
// For one reported iteration:
//
//     total measured time = sum_s time_s
//
//     percentage_s = 100 time_s / total measured time
//
// The profiler is intentionally lightweight. It is used to identify the
// dominant numerical cost before changing algorithms, OpenMP assembly, PETSc
// configuration or MPI decomposition.
//
// The profiler measures wall-clock time through std::chrono::steady_clock.
//=============================================================================

#include "cbs/core/Types.hpp"

#include <array>
#include <chrono>
#include <iosfwd>
#include <string_view>

namespace cbs
{
    class SolverProfiler
    {
    public:
        // Timed sections listed in the same general order as one CBS
        // iteration. Count is used only to determine the array size.
        enum class Section : std::size_t
        {
            TimeStepCompute = 0,
            LhsDiagonal,
            Step1Momentum,
            Step2TimeStepCorrection,
            Step2PressureRhs,
            Step2PressureSolve,
            Step3VelocityCorrection,
            Step4Energy,
            VelocityMagnitude,
            Convergence,
            PostOutput,
            Count
        };


        // Records the wall-clock time spent inside one C++ scope.
        //
        // The object owns no solver data. It stores only:
        //
        //     pointer to the parent profiler
        //     selected timing section
        //     starting wall-clock time
        //     active/inactive state
        //
        // Copying is disabled so that one timed scope cannot be counted twice.
        class ScopedTimer
        {
        public:
            ScopedTimer(
                SolverProfiler& profiler,
                Section section) noexcept;

            ScopedTimer(const ScopedTimer&) = delete;
            ScopedTimer& operator=(const ScopedTimer&) = delete;

            ScopedTimer(ScopedTimer&& other) noexcept;
            ScopedTimer& operator=(ScopedTimer&& other) = delete;

            // Adds the elapsed wall-clock time to the parent profiler.
            ~ScopedTimer();

        private:
            SolverProfiler* profiler_;
            Section section_;
            std::chrono::steady_clock::time_point start_;
            bool active_;
        };


        // Creates an empty profiler and resets all per-iteration values.
        SolverProfiler();

        // Clears all section times, call counts and Pressure CG information
        // before a new CBS iteration begins.
        void resetIteration();

        // Starts timing one solver section.
        ScopedTimer time(Section section) noexcept;

        // Adds one measured duration and increments the corresponding call
        // count.
        void addElapsed(
            Section section,
            double seconds) noexcept;

        // Stores Pressure CG information for the optional profiler summary.
        void setCgIterations(Int iterations) noexcept;
        void setCgInitialResidual(Real value) noexcept;
        void setCgFinalResidual(Real value) noexcept;
        void setCgRelativeResidual(Real value) noexcept;

        // Prints the section times for the selected CBS iteration interval.
        void printIteration(
            std::ostream& os,
            Int iteration,
            Int print_every) const;

        // Returns:
        //
        //     sum_s time_s
        double totalMeasuredSeconds() const noexcept;

        // Returns a readable label for one profiling section.
        static std::string_view sectionName(Section section) noexcept;

    private:
        static constexpr std::size_t nsection_ =
            static_cast<std::size_t>(Section::Count);

        // Accumulated wall-clock seconds and number of timed calls for the
        // current CBS iteration.
        std::array<double, nsection_> seconds_;
        std::array<Int, nsection_> calls_;

        // Pressure CG data printed with the iteration profile when available.
        Int cg_iterations_;
        Real cg_initial_residual_;
        Real cg_final_residual_;
        Real cg_relative_residual_;
    };
}
