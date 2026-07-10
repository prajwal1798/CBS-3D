#pragma once

//=============================================================================
// CBS3D++_SI
//
// Boundary and initial-value handling for the Spalart-Allmaras working variable.
// Partition boundaries are never treated as turbulence walls.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"

namespace cbs
{
    class TurbulenceBoundary
    {
    public:
        static void classifyNodes(CBSStateSI& s);
        static void initialiseNuTilde(CBSStateSI& s);
        static void applyWallValues(CBSStateSI& s);
        static void applyInletValues(CBSStateSI& s);
    };
}
