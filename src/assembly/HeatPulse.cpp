//=============================================================================
// CBS3D++_SI
//
// FILE      : HeatPulse.cpp
// PURPOSE   : PIECEWISE-LINEAR PLASMA-POWER FACTOR f(t), SEE HeatPulse.hpp.
//=============================================================================

#include "cbs/assembly/HeatPulse.hpp"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#ifdef CBS3D_USE_MPI
#include <mpi.h>
#endif

namespace cbs
{
    std::vector<Real> HeatPulse::times_;
    std::vector<Real> HeatPulse::values_;
    Real HeatPulse::period_ = 0.0;
    bool HeatPulse::active_ = false;

    namespace
    {
        // ================================
        // READ "t f" PAIRS; THROWS ON A MALFORMED OR NON-INCREASING TABLE
        // ================================
        void read_table(const std::string& path, std::vector<Real>& t, std::vector<Real>& f)
        {
            std::ifstream in(path);

            if (!in)
            {
                throw std::runtime_error("HeatPulse - cannot open CBS3D_HEAT_PULSE_FILE " + path);
            }

            std::string line;
            int line_no = 0;

            while (std::getline(in, line))
            {
                ++line_no;
                const std::size_t first = line.find_first_not_of(" \t\r");

                if (first == std::string::npos || line[first] == '#')
                {
                    continue;
                }

                std::istringstream row(line);
                Real a = 0.0;
                Real b = 0.0;

                if (!(row >> a >> b) || !std::isfinite(a) || !std::isfinite(b))
                {
                    throw std::runtime_error(
                        "HeatPulse - " + path + " line " + std::to_string(line_no) +
                        ": expected two numbers 't f'");
                }

                if (b < 0.0)
                {
                    throw std::runtime_error(
                        "HeatPulse - " + path + " line " + std::to_string(line_no) +
                        ": the factor f must be >= 0");
                }

                if (!t.empty() && !(a > t.back()))
                {
                    throw std::runtime_error(
                        "HeatPulse - " + path + " line " + std::to_string(line_no) +
                        ": times must increase strictly");
                }

                t.push_back(a);
                f.push_back(b);
            }

            if (t.empty())
            {
                throw std::runtime_error("HeatPulse - " + path + " has no 't f' rows");
            }
        }
    }

    bool HeatPulse::active()
    {
        return active_;
    }

    Real HeatPulse::factor(Real t)
    {
        if (!active_)
        {
            return 1.0;
        }

        if (period_ > 0.0)
        {
            t = std::fmod(t, period_);

            if (t < 0.0)
            {
                t += period_;
            }
        }

        if (t <= times_.front())
        {
            return values_.front();
        }

        if (t >= times_.back())
        {
            return values_.back();
        }

        std::size_t hi = 1;

        while (times_[hi] < t)
        {
            ++hi;
        }

        const std::size_t lo = hi - 1;
        const Real w = (t - times_[lo]) / (times_[hi] - times_[lo]);
        return (1.0 - w) * values_[lo] + w * values_[hi];
    }

    void HeatPulse::update(CBSStateSI& s)
    {
        s.heat_source_scale = factor(s.cfg.rtime);
    }

    void HeatPulse::initialise(CBSStateSI& s)
    {
        times_.clear();
        values_.clear();
        active_ = false;
        period_ = 0.0;

        const char* file = std::getenv("CBS3D_HEAT_PULSE_FILE");

        if (file == nullptr || std::string(file).empty())
        {
            s.heat_source_scale = 1.0;
            return;
        }

        // ================================
        // RANK 0 READS, EVERY RANK RECEIVES THE SAME TABLE
        // ================================
        std::vector<double> packed;
        std::string error;

        if (s.mpi_rank == 0)
        {
            try
            {
                std::vector<Real> t;
                std::vector<Real> f;
                read_table(file, t, f);

                for (std::size_t i = 0; i < t.size(); ++i)
                {
                    packed.push_back(t[i]);
                    packed.push_back(f[i]);
                }
            }
            catch (const std::exception& e)
            {
                error = e.what();
            }
        }

        int count = static_cast<int>(packed.size());

#ifdef CBS3D_USE_MPI
        int failed = error.empty() ? 0 : 1;

        if (s.mpi_enabled)
        {
            MPI_Bcast(&failed, 1, MPI_INT, 0, MPI_COMM_WORLD);
        }

        if (failed)
        {
            throw std::runtime_error(error.empty() ? "HeatPulse - table rejected on rank 0" : error);
        }

        if (s.mpi_enabled)
        {
            MPI_Bcast(&count, 1, MPI_INT, 0, MPI_COMM_WORLD);
            packed.resize(static_cast<std::size_t>(count));
            MPI_Bcast(packed.data(), count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
#else
        if (!error.empty())
        {
            throw std::runtime_error(error);
        }
#endif

        for (int i = 0; i + 1 < count; i += 2)
        {
            times_.push_back(packed[static_cast<std::size_t>(i)]);
            values_.push_back(packed[static_cast<std::size_t>(i + 1)]);
        }

        const char* period = std::getenv("CBS3D_HEAT_PULSE_PERIOD");

        if (period != nullptr && std::string(period).size() > 0)
        {
            period_ = std::atof(period);

            if (!(period_ > 0.0))
            {
                throw std::runtime_error("HeatPulse - CBS3D_HEAT_PULSE_PERIOD must be > 0");
            }
        }

        active_ = true;
        update(s);

        if (s.mpi_rank == 0)
        {
            std::cout << "HEAT PULSE  : " << file << " (" << times_.size() << " points";

            if (period_ > 0.0)
            {
                std::cout << ", period " << period_ << " s";
            }

            std::cout << ")\n";

            for (std::size_t i = 0; i < times_.size(); ++i)
            {
                std::cout << "  t = " << times_[i] << " s   f = " << values_[i] << "\n";
            }

            std::cout << "  f(t = " << s.cfg.rtime << " s) = " << s.heat_source_scale
                      << "  (SCALES BC 532 FLUX AND ALL VOLUMETRIC SOURCES)\n";
        }
    }
}
