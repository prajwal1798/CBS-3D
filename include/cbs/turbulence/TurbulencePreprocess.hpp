#pragma once

//=============================================================================
// CBS3D++_SI
//
// Preprocessing entry point for the Spalart-Allmaras turbulence model.
//
// This routine is called after the usual mesh, material and boundary
// preprocessing has finished and before the CBS time loop begins.  It prepares
// all geometry-dependent SA quantities that do not change during the simulation.
//
// The most important quantity prepared here is the wall distance d.  The wall
// distance is not a time-dependent unknown.  It is a fixed geometric field and
// must therefore be computed once, stored, and reused during SA assembly.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"

namespace cbs
{
    class TurbulencePreprocess
    {
    public:
        static void prepareSpalartAllmaras(CBSStateSI& s);
    };
}
