#pragma once

//=============================================================================
// CBS3D++_SI
//
// FILE      : TimeAverage.hpp
// PURPOSE   : RUNNING TIME AVERAGES OF A TRANSIENT (GLOBAL TIME STEP) MPI RUN.
//             FOR EVERY LOCAL NODE THE TIME INTEGRALS
//
//                 INT u dt,  INT u_i^2 dt,  INT p dt,  INT T' dt,  INT T'^2 dt
//
//             (T' = T - T AT THE FIRST SAMPLE, AGAINST ROUND-OFF IN THE VARIANCE)
//
//             AND THE WINDOW MAXIMUM OF T ARE ACCUMULATED AFTER EACH STEP
//             (RECTANGLE RULE: f^{n+1} dt^n). THE VTU THEN CARRIES
//
//                 mean_velocity, rms_velocity (PER COMPONENT), mean_pressure,
//                 mean_temperature, rms_temperature, max_temperature
//
//             WITH mean = INT f dt / T_w AND rms = sqrt(INT f^2 dt / T_w - mean^2).
//
// CONTROLS  : ENVIRONMENT (AVERAGING IS OFF WHEN NEITHER IS SET)
//   CBS3D_AVERAGE_START_TIME=t0       START AT PHYSICAL TIME t0 [s]
//                                     (NEEDS transient_on = 1)
//   CBS3D_AVERAGE_START_ITERATION=n0  START AT GLOBAL ITERATION n0
//
//             A LOCAL-TIME-STEPPING RUN (ilots = 1/2, transient_on = 0) HAS
//             NO COMMON TIME, SO AVERAGING IS REFUSED THERE.
//
// RESTART   : EVERY NATIVE CHECKPOINT ALSO HOLDS
//             time_average_rank_XXXX.bin; A RESTART FROM THAT CHECKPOINT
//             (CBS3D_RESTART=1, SAME RANK COUNT) CONTINUES THE SAME WINDOW.
//=============================================================================

#include "cbs/core/Array.hpp"
#include "cbs/core/CBSStateSI.hpp"

#include <filesystem>
#include <string>

namespace cbs
{
    class TimeAverage
    {
    public:
        // READS THE CONTROLS, ALLOCATES THE ACCUMULATORS AND, ON A RESTART,
        // LOADS THE ACCUMULATORS OF THE CHECKPOINT WHEN THEY EXIST.
        static void initialise(CBSStateSI& s);

        // ADDS THE CURRENT FIELDS WITH WEIGHT dt ONCE THE WINDOW HAS STARTED.
        static void accumulate(const CBSStateSI& s, Int global_iteration, Real dt);

        [[nodiscard]] static bool enabled() { return enabled_; }

        // TRUE WHEN AT LEAST ONE STEP HAS BEEN ACCUMULATED (THE VTU ARRAYS
        // ARE WRITTEN ONLY THEN; PIECE AND .pvtu USE THE SAME TEST).
        [[nodiscard]] static bool hasSamples() { return enabled_ && window_ > 0.0; }

        [[nodiscard]] static Real window() { return window_; }
        [[nodiscard]] static Int samples() { return samples_; }

        // MEAN AND RMS AT ONE NODE (dim = 1..3 FOR VELOCITY).
        [[nodiscard]] static Real meanVelocity(Int dim, Int ip);
        [[nodiscard]] static Real rmsVelocity(Int dim, Int ip);
        [[nodiscard]] static Real meanPressure(Int ip);
        [[nodiscard]] static Real meanTemperature(Int ip);
        [[nodiscard]] static Real rmsTemperature(Int ip);
        [[nodiscard]] static Real maxTemperature(Int ip);

        // RANK FILE IN A CHECKPOINT DIRECTORY (WRITTEN THROUGH A .tmp + RENAME).
        static void writeRankFile(const CBSStateSI& s,
                                  const std::filesystem::path& directory);

        [[nodiscard]] static std::filesystem::path rankFilePath(
            const std::filesystem::path& directory, int rank);

    private:
        inline static bool enabled_ = false;
        inline static bool by_time_ = false;
        inline static Real start_time_ = 0.0;
        inline static Int start_iteration_ = 0;
        inline static Real window_ = 0.0;
        inline static Int samples_ = 0;
        inline static Int npoin_ = 0;

        inline static Array2D<Real> sum_u_;
        inline static Array2D<Real> sum_uu_;
        inline static Array1D<Real> sum_p_;
        inline static Array1D<Real> sum_t_;
        inline static Array1D<Real> sum_tt_;
        inline static Array1D<Real> max_t_;
        inline static Array1D<Real> ref_t_;   // FIRST SAMPLE OF T (SHIFT)
    };
}
