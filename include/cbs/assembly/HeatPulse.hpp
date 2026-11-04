#pragma once

//=============================================================================
// CBS3D++_SI
//
// FILE      : HeatPulse.hpp
// PURPOSE   : TIME-DEPENDENT PLASMA POWER FOR TRANSIENT RUNS.
//
//             ALL PRESCRIBED HEATING IS MULTIPLIED BY ONE FACTOR f(t):
//
//                 q''_wall(x, t) = f(t) q''_wall(x)      BC 532 SURFACE FLUX
//                 Q'''(x, t)     = f(t) Q'''(x)          VOLUMETRIC SOURCES
//                                                        (NUCLEAR PROFILES,
//                                                         .matprop Qvol,
//                                                         source_solid)
//
//             f(t) IS PIECEWISE LINEAR BETWEEN THE POINTS OF A TWO-COLUMN
//             TABLE "t  f" (SIMULATION TIME [s], FACTOR [-]); CONSTANT
//             BEFORE THE FIRST AND AFTER THE LAST POINT. LINES STARTING
//             WITH # ARE COMMENTS.
//
//             t IS THE SOLVER'S SIMULATION TIME rtime: 0 AT THE START OF A
//             RUN THAT IMPORTS A VTU AS ITS INITIAL CONDITION, AND RESTORED
//             FROM A NATIVE CHECKPOINT, SO A RESUMED RUN CONTINUES THE SAME
//             PULSE.
//
// USAGE     : CBS3D_HEAT_PULSE_FILE=path/to/pulse.txt
//             CBS3D_HEAT_PULSE_PERIOD=T   (OPTIONAL: REPEAT THE TABLE WITH
//                                          PERIOD T, USING t MOD T)
//             WITHOUT CBS3D_HEAT_PULSE_FILE, f = 1 (UNCHANGED BEHAVIOUR).
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"

#include <vector>

namespace cbs
{
    class HeatPulse
    {
    public:
        // READS THE TABLE (RANK 0) AND BROADCASTS IT; SETS THE SCALE FOR THE
        // CURRENT TIME. PRINTS THE TABLE ON RANK 0.
        static void initialise(CBSStateSI& s);

        // s.heat_source_scale = f(s.cfg.rtime)
        static void update(CBSStateSI& s);

        static bool active();

        static Real factor(Real t);

    private:
        static std::vector<Real> times_;
        static std::vector<Real> values_;
        static Real period_;
        static bool active_;
    };
}
