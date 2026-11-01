//=============================================================================
// CBS3D++_SI
//
// FILE      : TimeAverage.cpp
// PURPOSE   : RUNNING TIME AVERAGES OF A TRANSIENT RUN (SEE TimeAverage.hpp).
//=============================================================================

#include "cbs/io/TimeAverage.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;

namespace cbs
{
    namespace
    {
        constexpr char average_magic[8] = {'C', 'B', 'S', '3', 'D', 'A', 'V', 'G'};
        // VERSION 2 ADDS sum_ut_ (THE u_i T COVARIANCE ACCUMULATOR).
        // A VERSION 1 RANK FILE CANNOT SUPPLY <u'T'> FOR THE WHOLE WINDOW
        // AND IS REJECTED; START A NEW WINDOW WITH THIS BUILD.
        constexpr std::uint32_t average_version = 2U;

        std::string env_text(const char* name)
        {
            const char* value = std::getenv(name);
            return (value == nullptr) ? std::string() : std::string(value);
        }

        bool env_flag(const char* name)
        {
            const std::string v = env_text(name);
            return !v.empty() && v != "0" && v != "false" && v != "FALSE" &&
                   v != "off" && v != "OFF" && v != "no" && v != "NO";
        }

        Real parse_real(const char* name, const std::string& text)
        {
            std::size_t used = 0;
            Real value = 0.0;
            try
            {
                value = std::stod(text, &used);
            }
            catch (const std::exception&)
            {
                used = 0;
            }
            if (used != text.size() || !std::isfinite(value) || value < 0.0)
            {
                throw std::runtime_error(std::string(name) +
                                         " MUST BE A NON-NEGATIVE NUMBER");
            }
            return value;
        }

        Int parse_int(const char* name, const std::string& text)
        {
            std::size_t used = 0;
            long long value = -1;
            try
            {
                value = std::stoll(text, &used, 10);
            }
            catch (const std::exception&)
            {
                used = 0;
            }
            if (used != text.size() || value < 0 ||
                value > std::numeric_limits<Int>::max())
            {
                throw std::runtime_error(std::string(name) +
                                         " MUST BE A NON-NEGATIVE INTEGER");
            }
            return static_cast<Int>(value);
        }

        bool local_time_stepping(const CBSStateSI& s)
        {
            return (s.cfg.ilots == 1 || s.cfg.ilots == 2) &&
                   s.cfg.transient_on == 0 && s.cfg.dtfixed == 0;
        }

        template <typename T>
        void put(std::ofstream& out, const T& value)
        {
            out.write(reinterpret_cast<const char*>(&value), sizeof(T));
        }

        template <typename T>
        T get(std::ifstream& in)
        {
            T value{};
            in.read(reinterpret_cast<char*>(&value), sizeof(T));
            if (!in)
            {
                throw std::runtime_error("TimeAverage: TRUNCATED RANK FILE");
            }
            return value;
        }

        void put_raw(std::ofstream& out, const std::vector<Real>& data)
        {
            // POSITION ZERO IS UNUSED BY Array1D/Array2D; IT IS WRITTEN TOO
            out.write(reinterpret_cast<const char*>(data.data()),
                      static_cast<std::streamsize>(data.size() * sizeof(Real)));
        }

        void get_raw(std::ifstream& in, std::vector<Real>& data)
        {
            in.read(reinterpret_cast<char*>(data.data()),
                    static_cast<std::streamsize>(data.size() * sizeof(Real)));
            if (!in)
            {
                throw std::runtime_error("TimeAverage: TRUNCATED RANK FILE");
            }
        }
    }

    //=========================================================================
    // CONTROLS, ALLOCATION, RESTART
    //=========================================================================
    void TimeAverage::initialise(CBSStateSI& s)
    {
        const std::string t0 = env_text("CBS3D_AVERAGE_START_TIME");
        const std::string n0 = env_text("CBS3D_AVERAGE_START_ITERATION");

        enabled_ = false;
        window_ = 0.0;
        samples_ = 0;

        if (t0.empty() && n0.empty())
        {
            return;
        }
        if (!t0.empty() && !n0.empty())
        {
            throw std::runtime_error(
                "TimeAverage: SET CBS3D_AVERAGE_START_TIME OR "
                "CBS3D_AVERAGE_START_ITERATION, NOT BOTH");
        }
        if (local_time_stepping(s))
        {
            throw std::runtime_error(
                "TimeAverage: LOCAL TIME STEPPING (ilots = 1/2, transient_on = 0) "
                "HAS NO COMMON TIME; USE transient_on = 1 FOR A TIME AVERAGE");
        }

        by_time_ = !t0.empty();
        if (by_time_)
        {
            if (s.cfg.transient_on < 1)
            {
                throw std::runtime_error(
                    "TimeAverage: CBS3D_AVERAGE_START_TIME NEEDS transient_on = 1 "
                    "(THE PHYSICAL TIME IS ONLY ADVANCED THEN)");
            }
            start_time_ = parse_real("CBS3D_AVERAGE_START_TIME", t0);
        }
        else
        {
            start_iteration_ = parse_int("CBS3D_AVERAGE_START_ITERATION", n0);
        }

        npoin_ = s.cfg.npoin;
        sum_u_.resize(3, npoin_);
        sum_uu_.resize(3, npoin_);
        sum_ut_.resize(3, npoin_);
        sum_p_.resize(npoin_);
        sum_t_.resize(npoin_);
        sum_tt_.resize(npoin_);
        max_t_.resize(npoin_);
        ref_t_.resize(npoin_);
        sum_u_.fill(0.0);
        sum_uu_.fill(0.0);
        sum_ut_.fill(0.0);
        sum_p_.fill(0.0);
        sum_t_.fill(0.0);
        sum_tt_.fill(0.0);
        max_t_.fill(-std::numeric_limits<Real>::infinity());
        ref_t_.fill(0.0);
        enabled_ = true;

        // RESTART: CONTINUE THE WINDOW OF THE CHECKPOINT WHEN IT HAS ONE
        bool loaded = false;
        const bool restart = s.cfg.restart_opt > 0 || env_flag("CBS3D_RESTART");
        const std::string root = env_text("CBS3D_RESTART_ROOT");
        if (restart && !root.empty())
        {
            const fs::path file = rankFilePath(fs::path(root), s.mpi_rank);
            std::error_code ec;
            if (fs::exists(file, ec))
            {
                std::ifstream in(file, std::ios::binary);
                char magic[8];
                in.read(magic, 8);
                if (!in || std::memcmp(magic, average_magic, 8) != 0)
                {
                    throw std::runtime_error("TimeAverage: BAD RANK FILE " + file.string());
                }
                const auto version = get<std::uint32_t>(in);
                const auto npoin = get<Int>(in);
                if (version != average_version)
                {
                    throw std::runtime_error(
                        "TimeAverage: RANK FILE " + file.string() +
                        " HAS AVERAGE VERSION " + std::to_string(version) +
                        "; THIS BUILD WRITES VERSION " +
                        std::to_string(average_version) +
                        " AND NEEDS <u'T'> FOR THE WHOLE WINDOW. "
                        "START A NEW WINDOW (CBS3D_AVERAGE_START_TIME)");
                }
                if (npoin != npoin_)
                {
                    throw std::runtime_error(
                        "TimeAverage: RANK FILE " + file.string() +
                        " DOES NOT MATCH THIS PARTITION");
                }
                window_ = get<Real>(in);
                samples_ = get<Int>(in);
                get_raw(in, sum_u_.raw());
                get_raw(in, sum_uu_.raw());
                get_raw(in, sum_p_.raw());
                get_raw(in, sum_t_.raw());
                get_raw(in, sum_tt_.raw());
                get_raw(in, max_t_.raw());
                get_raw(in, ref_t_.raw());
                get_raw(in, sum_ut_.raw());
                loaded = true;
            }
        }

        if (s.mpi_rank == 0)
        {
            std::cout << "Time averaging\n"
                      << "  start      : "
                      << (by_time_ ? "physical time " + std::to_string(start_time_) + " s"
                                   : "global iteration " + std::to_string(start_iteration_))
                      << "\n"
                      << "  restart    : "
                      << (loaded ? "window continued from the checkpoint ("
                                       + std::to_string(window_) + " s, "
                                       + std::to_string(samples_) + " steps)"
                                 : "new window")
                      << "\n"
                      << "  output     : mean_velocity rms_velocity mean_pressure "
                         "mean_temperature rms_temperature max_temperature "
                         "covariance_velocity_temperature\n";
        }
    }

    //=========================================================================
    // ACCUMULATION AFTER ONE COMPLETED STEP
    //=========================================================================
    void TimeAverage::accumulate(const CBSStateSI& s, const Int global_iteration,
                                 const Real dt)
    {
        if (!enabled_)
        {
            return;
        }
        if (by_time_ ? (s.cfg.rtime < start_time_)
                     : (global_iteration < start_iteration_))
        {
            return;
        }
        if (!(dt > 0.0) || !std::isfinite(dt))
        {
            throw std::runtime_error("TimeAverage: INVALID TIME STEP");
        }

        const bool energy = s.cfg.temp_calc > 0;
        if (energy && samples_ == 0)
        {
            // TEMPERATURE IS ACCUMULATED RELATIVE TO ITS FIRST SAMPLE, SO THE
            // VARIANCE IS NOT A SMALL DIFFERENCE OF TWO NUMBERS ~ (300 K)^2
            for (Int ip = 1; ip <= npoin_; ++ip)
            {
                ref_t_(ip) = s.temperature(ip);
            }
        }
        for (Int ip = 1; ip <= npoin_; ++ip)
        {
            for (Int d = 1; d <= 3; ++d)
            {
                const Real u = s.unkno(d, ip);
                sum_u_(d, ip) += u * dt;
                sum_uu_(d, ip) += u * u * dt;
            }
            sum_p_(ip) += s.pres(ip) * dt;
            if (energy)
            {
                const Real t = s.temperature(ip);
                const Real dtheta = t - ref_t_(ip);
                sum_t_(ip) += dtheta * dt;
                sum_tt_(ip) += dtheta * dtheta * dt;
                for (Int d = 1; d <= 3; ++d)
                {
                    // <u_i T> - <u_i><T> = <u_i theta> - <u_i><theta>
                    sum_ut_(d, ip) += s.unkno(d, ip) * dtheta * dt;
                }
                max_t_(ip) = std::max(max_t_(ip), t);
            }
        }
        window_ += dt;
        ++samples_;
    }

    //=========================================================================
    // NODAL STATISTICS
    //=========================================================================
    Real TimeAverage::meanVelocity(const Int dim, const Int ip)
    {
        return sum_u_(dim, ip) / window_;
    }

    Real TimeAverage::rmsVelocity(const Int dim, const Int ip)
    {
        const Real m = sum_u_(dim, ip) / window_;
        return std::sqrt(std::max(0.0, sum_uu_(dim, ip) / window_ - m * m));
    }

    Real TimeAverage::meanPressure(const Int ip)
    {
        return sum_p_(ip) / window_;
    }

    Real TimeAverage::meanTemperature(const Int ip)
    {
        return ref_t_(ip) + sum_t_(ip) / window_;
    }

    Real TimeAverage::rmsTemperature(const Int ip)
    {
        // sum_t_ AND sum_tt_ HOLD T - T_first, SO NO (300 K)^2 CANCELLATION
        const Real m = sum_t_(ip) / window_;
        return std::sqrt(std::max(0.0, sum_tt_(ip) / window_ - m * m));
    }

    Real TimeAverage::maxTemperature(const Int ip)
    {
        return max_t_(ip);
    }

    Real TimeAverage::covarianceVelocityTemperature(const Int dim, const Int ip)
    {
        // THE SHIFT CANCELS: <u theta> - <u><theta> = <u T> - <u><T>,
        // WITH theta = T - T_FIRST.
        const Real mu = sum_u_(dim, ip) / window_;
        const Real mt = sum_t_(ip) / window_;
        return sum_ut_(dim, ip) / window_ - mu * mt;
    }

    //=========================================================================
    // CHECKPOINT RANK FILE
    //=========================================================================
    fs::path TimeAverage::rankFilePath(const fs::path& directory, const int rank)
    {
        char name[64];
        std::snprintf(name, sizeof(name), "time_average_rank_%04d.bin", rank);
        return directory / name;
    }

    void TimeAverage::writeRankFile(const CBSStateSI& s, const fs::path& directory)
    {
        if (!enabled_)
        {
            return;
        }
        const fs::path final_path = rankFilePath(directory, s.mpi_rank);
        const fs::path tmp_path = final_path.string() + ".tmp";
        {
            std::ofstream out(tmp_path, std::ios::binary | std::ios::trunc);
            if (!out)
            {
                throw std::runtime_error("TimeAverage: CANNOT CREATE " + tmp_path.string());
            }
            out.write(average_magic, 8);
            put(out, average_version);
            put(out, npoin_);
            put(out, window_);
            put(out, samples_);
            put_raw(out, sum_u_.raw());
            put_raw(out, sum_uu_.raw());
            put_raw(out, sum_p_.raw());
            put_raw(out, sum_t_.raw());
            put_raw(out, sum_tt_.raw());
            put_raw(out, max_t_.raw());
            put_raw(out, ref_t_.raw());
            put_raw(out, sum_ut_.raw());
            out.flush();
            if (!out)
            {
                throw std::runtime_error("TimeAverage: WRITE FAILED " + tmp_path.string());
            }
        }
        std::error_code ec;
        fs::remove(final_path, ec);
        fs::rename(tmp_path, final_path);
    }
}
