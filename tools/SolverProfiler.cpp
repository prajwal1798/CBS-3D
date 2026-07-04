#include "cbs/utils/SolverProfiler.hpp"

#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace cbs
{
    SolverProfiler::ScopedTimer::ScopedTimer(SolverProfiler& profiler, const Section section)
        : profiler_(profiler),
          section_(section),
          start_(std::chrono::steady_clock::now())
    {
    }

    SolverProfiler::ScopedTimer::~ScopedTimer()
    {
        const auto stop = std::chrono::steady_clock::now();
        const std::chrono::duration<double> elapsed = stop - start_;
        profiler_.addTime(section_, elapsed.count());
    }

    SolverProfiler::SolverProfiler()
        : entries_(static_cast<std::size_t>(Section::Count))
    {
        for (std::size_t i = 0; i < entries_.size(); ++i)
        {
            entries_[i].name = sectionName(static_cast<Section>(i));
        }
    }

    void SolverProfiler::resetIteration()
    {
        for (Entry& e : entries_)
        {
            e.seconds = 0.0;
            e.count = 0;
        }

        cg_iterations_ = 0;
        cg_initial_residual_ = 0.0;
        cg_final_residual_ = 0.0;
        cg_relative_residual_ = 0.0;
    }

    void SolverProfiler::addTime(const Section section, const double seconds)
    {
        const auto i = static_cast<std::size_t>(section);
        if (i >= entries_.size())
        {
            throw std::runtime_error("SolverProfiler::addTime - invalid profiler section");
        }

        entries_[i].seconds += seconds;
        ++entries_[i].count;
    }

    void SolverProfiler::addCount(const Section section, const std::size_t count)
    {
        const auto i = static_cast<std::size_t>(section);
        if (i >= entries_.size())
        {
            throw std::runtime_error("SolverProfiler::addCount - invalid profiler section");
        }

        entries_[i].count += count;
    }

    SolverProfiler::ScopedTimer SolverProfiler::time(const Section section)
    {
        return ScopedTimer(*this, section);
    }

    void SolverProfiler::setCgIterations(const int iterations)
    {
        cg_iterations_ = iterations;
    }

    void SolverProfiler::setCgInitialResidual(const double value)
    {
        cg_initial_residual_ = value;
    }

    void SolverProfiler::setCgFinalResidual(const double value)
    {
        cg_final_residual_ = value;
    }

    void SolverProfiler::setCgRelativeResidual(const double value)
    {
        cg_relative_residual_ = value;
    }

    double SolverProfiler::totalIterationSeconds() const
    {
        double total = 0.0;
        for (const Entry& e : entries_)
        {
            total += e.seconds;
        }
        return total;
    }

    void SolverProfiler::printIteration(std::ostream& os, const int iteration, const int printEvery) const
    {
        if (printEvery <= 0 || iteration % printEvery != 0)
        {
            return;
        }

        const double total = totalIterationSeconds();

        os << "\n[Profiler] iteration " << iteration
           << "  total_measured=" << std::fixed << std::setprecision(6) << total << " s\n";

        os << "  " << std::left << std::setw(30) << "section"
           << std::right << std::setw(14) << "seconds"
           << std::setw(12) << "percent"
           << std::setw(12) << "calls" << "\n";

        for (const Entry& e : entries_)
        {
            if (e.seconds <= 0.0 && e.count == 0)
            {
                continue;
            }

            const double pct = (total > 0.0) ? (100.0 * e.seconds / total) : 0.0;
            os << "  " << std::left << std::setw(30) << e.name
               << std::right << std::setw(14) << std::fixed << std::setprecision(6) << e.seconds
               << std::setw(11) << std::fixed << std::setprecision(2) << pct << "%"
               << std::setw(12) << e.count << "\n";
        }

        os << "  CG summary: iterations=" << cg_iterations_
           << " initial=" << std::scientific << cg_initial_residual_
           << " final=" << cg_final_residual_
           << " relative=" << cg_relative_residual_ << std::defaultfloat << "\n";
    }

    const char* SolverProfiler::sectionName(const Section section)
    {
        switch (section)
        {
        case Section::TimeStepCompute: return "TimeStep::compute";
        case Section::LhsDiagonal: return "TimeStep::updateLHS";
        case Section::Step1Momentum: return "Step1 momentum";
        case Section::BoundaryAfterStep1: return "Boundary after Step1";
        case Section::Step2PressureRhs: return "Step2 pressure RHS";
        case Section::PressureCg: return "Pressure CG total";
        case Section::PressureCgMatvec: return "CG pressure matvec";
        case Section::PressureCgDot: return "CG dot products";
        case Section::PressureCgPreconditioner: return "CG preconditioner";
        case Section::PressureCgUpdate: return "CG vector updates";
        case Section::Step3VelocityCorrection: return "Step3 velocity correction";
        case Section::BoundaryAfterStep3: return "Boundary after Step3";
        case Section::Step4Energy: return "Step4 energy";
        case Section::VelocityMagnitude: return "Velocity magnitude";
        case Section::Convergence: return "Convergence";
        case Section::PostOutput: return "Post output";
        case Section::Count: return "Count";
        }

        return "Unknown";
    }
}
