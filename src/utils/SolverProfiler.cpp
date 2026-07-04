//=============================================================================
// CBS3D++_SI
//
// Per-iteration wall-clock profiler for the main CBS solver stages.
//
// Each ScopedTimer measures:
//
//     elapsed time = finish time - start time
//
// The measured duration is accumulated in the selected solver section. The
// final report prints:
//
//     total measured time = sum_s time_s
//
//     percentage_s = 100 time_s / total measured time
//
// together with the number of times each section was entered.
//=============================================================================

#include "cbs/utils/SolverProfiler.hpp"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <ostream>

namespace cbs
{
    //=========================================================================
    // Starts timing one solver section.
    //=========================================================================
    SolverProfiler::ScopedTimer::ScopedTimer(
        SolverProfiler& profiler,
        const Section section) noexcept
        : profiler_(&profiler),
          section_(section),
          start_(std::chrono::steady_clock::now()),
          active_(true)
    {
    }


    //=========================================================================
    // Transfers ownership of an active timer.
    //
    // The moved-from timer is disabled so that the same elapsed interval is
    // not added twice.
    //=========================================================================
    SolverProfiler::ScopedTimer::ScopedTimer(ScopedTimer&& other) noexcept
        : profiler_(other.profiler_),
          section_(other.section_),
          start_(other.start_),
          active_(other.active_)
    {
        other.profiler_ = nullptr;
        other.active_ = false;
    }


    //=========================================================================
    // Stops the timer and adds its elapsed wall-clock time to the selected
    // profiler section.
    //=========================================================================
    SolverProfiler::ScopedTimer::~ScopedTimer()
    {
        if (!active_ || profiler_ == nullptr)
        {
            return;
        }

        const auto finish = std::chrono::steady_clock::now();
        const std::chrono::duration<double> elapsed = finish - start_;
        profiler_->addElapsed(section_, elapsed.count());
    }


    //=========================================================================
    // Creates a profiler with all per-iteration values set to zero.
    //=========================================================================
    SolverProfiler::SolverProfiler()
        : seconds_{},
          calls_{},
          cg_iterations_(0),
          cg_initial_residual_(0.0),
          cg_final_residual_(0.0),
          cg_relative_residual_(0.0)
    {
        resetIteration();
    }


    //=========================================================================
    // Clears all timing and Pressure CG information before the next CBS
    // iteration begins.
    //=========================================================================
    void SolverProfiler::resetIteration()
    {
        seconds_.fill(0.0);
        calls_.fill(0);

        cg_iterations_ = 0;
        cg_initial_residual_ = 0.0;
        cg_final_residual_ = 0.0;
        cg_relative_residual_ = 0.0;
    }


    //=========================================================================
    // Creates one scope timer for the requested solver section.
    //=========================================================================
    SolverProfiler::ScopedTimer SolverProfiler::time(const Section section) noexcept
    {
        return ScopedTimer(*this, section);
    }


    //=========================================================================
    // Adds one measured duration to a profiling section.
    //
    // For section s:
    //
    //     seconds_s <- seconds_s + elapsed
    //
    //     calls_s   <- calls_s + 1
    //=========================================================================
    void SolverProfiler::addElapsed(
        const Section section,
        const double seconds) noexcept
    {
        const auto idx = static_cast<std::size_t>(section);

        if (idx >= nsection_)
        {
            return;
        }

        seconds_[idx] += seconds;
        ++calls_[idx];
    }


    //=========================================================================
    // Stores the number of Pressure CG iterations completed during the current
    // CBS iteration.
    //=========================================================================
    void SolverProfiler::setCgIterations(const Int iterations) noexcept
    {
        cg_iterations_ = iterations;
    }


    //=========================================================================
    // Stores the initial Pressure CG L2 residual.
    //=========================================================================
    void SolverProfiler::setCgInitialResidual(const Real value) noexcept
    {
        cg_initial_residual_ = value;
    }


    //=========================================================================
    // Stores the final Pressure CG L2 residual.
    //=========================================================================
    void SolverProfiler::setCgFinalResidual(const Real value) noexcept
    {
        cg_final_residual_ = value;
    }


    //=========================================================================
    // Stores the final relative Pressure CG residual.
    //=========================================================================
    void SolverProfiler::setCgRelativeResidual(const Real value) noexcept
    {
        cg_relative_residual_ = value;
    }


    //=========================================================================
    // Returns the total measured time for the current CBS iteration:
    //
    //     total = sum_s seconds_s
    //
    // Only explicitly timed sections contribute to this total.
    //=========================================================================
    double SolverProfiler::totalMeasuredSeconds() const noexcept
    {
        double total = 0.0;

        for (const double value : seconds_)
        {
            total += value;
        }

        return total;
    }


    //=========================================================================
    // Returns the terminal label associated with one profiler section.
    //=========================================================================
    std::string_view SolverProfiler::sectionName(const Section section) noexcept
    {
        switch (section)
        {
        case Section::TimeStepCompute:
            return "TimeStep::computeTimeStep";
        case Section::LhsDiagonal:
            return "TimeStep::updateLhsDiagonal";
        case Section::Step1Momentum:
            return "Steps::step1 momentum";
        case Section::Step2TimeStepCorrection:
            return "Step2 pressure dt correction";
        case Section::Step2PressureRhs:
            return "Steps::step2 pressure RHS/CG";
        case Section::Step2PressureSolve:
            return "Steps::step2 pressure solve";
        case Section::Step3VelocityCorrection:
            return "Steps::step3 velocity correction";
        case Section::Step4Energy:
            return "Steps::step4 energy";
        case Section::VelocityMagnitude:
            return "Solver::updateVelocityMagnitude";
        case Section::Convergence:
            return "Convergence::evaluate";
        case Section::PostOutput:
            return "Post::writeSolution";
        case Section::Count:
            return "invalid";
        }

        return "unknown";
    }


    //=========================================================================
    // Prints the timing profile at the requested iteration interval.
    //
    // A report is produced only when:
    //
    //     print_every > 0
    //
    // and:
    //
    //     iteration mod print_every = 0
    //
    // For each active section:
    //
    //     percentage_s =
    //         100 seconds_s / total measured seconds
    //
    // Sections with zero calls are omitted.
    //=========================================================================
    void SolverProfiler::printIteration(
        std::ostream& os,
        const Int iteration,
        const Int print_every) const
    {
        if (print_every <= 0 || iteration <= 0 || (iteration % print_every) != 0)
        {
            return;
        }

        const double total = totalMeasuredSeconds();

        os << "\n[Profiler] iteration " << iteration
           << "  total_measured=" << std::fixed << std::setprecision(6)
           << total << " s\n";

        os << "  "
           << std::left << std::setw(36) << "section"
           << std::right << std::setw(14) << "seconds"
           << std::setw(12) << "percent"
           << std::setw(10) << "calls"
           << "\n";

        for (std::size_t i = 0; i < nsection_; ++i)
        {
            if (calls_[i] == 0)
            {
                continue;
            }

            const auto section = static_cast<Section>(i);
            const double pct = (total > 0.0) ? (100.0 * seconds_[i] / total) : 0.0;

            os << "  "
               << std::left << std::setw(36) << sectionName(section)
               << std::right << std::setw(14) << std::fixed << std::setprecision(6)
               << seconds_[i]
               << std::setw(11) << std::fixed << std::setprecision(2)
               << pct << "%"
               << std::setw(10) << calls_[i]
               << "\n";
        }

        if (cg_iterations_ > 0)
        {
            os << "  CG summary:"
               << " iterations=" << cg_iterations_
               << " initial_l2=" << std::scientific << std::setprecision(6)
               << cg_initial_residual_
               << " final_l2=" << cg_final_residual_
               << " relative=" << cg_relative_residual_
               << std::fixed << "\n";
        }

        os << std::flush;
    }
}
