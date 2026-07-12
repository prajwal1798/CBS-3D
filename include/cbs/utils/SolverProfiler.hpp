#pragma once

//=============================================================================
// CBS3D++_SI
//
// Per-iteration wall-clock profiler for the main solver stages.
//
// The profiler measures the elapsed time of selected parts of one solver
// iteration:
//
//     time-step calculation
//     left-hand-side diagonal construction
//     CBS Steps 1 to 4
//     optional Spalart-Allmaras turbulence step
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
        // Timed sections listed in the same general order as one solver
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
            StepSpalartAllmaras,
            Step4Energy,
            VelocityMagnitude,
            Convergence,
            PostOutput,
            Count
        };


        // Records the wall-clock time spent inside one C++ scope.
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

            ~ScopedTimer();

        private:
            SolverProfiler* profiler_;
            Section section_;
            std::chrono::steady_clock::time_point start_;
            bool active_;
        };


        SolverProfiler();

        void resetIteration();

        ScopedTimer time(Section section) noexcept;

        void addElapsed(
            Section section,
            double seconds) noexcept;

        void setCgIterations(Int iterations) noexcept;
        void setCgInitialResidual(Real value) noexcept;
        void setCgFinalResidual(Real value) noexcept;
        void setCgRelativeResidual(Real value) noexcept;

        void printIteration(
            std::ostream& os,
            Int iteration,
            Int print_every) const;

        double totalMeasuredSeconds() const noexcept;

        static std::string_view sectionName(Section section) noexcept;

    private:
        static constexpr std::size_t nsection_ =
            static_cast<std::size_t>(Section::Count);

        std::array<double, nsection_> seconds_;
        std::array<Int, nsection_> calls_;

        Int cg_iterations_;
        Real cg_initial_residual_;
        Real cg_final_residual_;
        Real cg_relative_residual_;
    };
}
