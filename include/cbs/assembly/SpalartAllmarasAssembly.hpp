#pragma once

//=============================================================================
// CBS3D++_SI
//
// Spalart-Allmaras turbulence assembly and effective-property update entry
// points.  Milestone 2 intentionally provides the state-safe scaffolding and
// algebraic eddy-viscosity update, while the full transport equation is added in
// a later milestone.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"

namespace cbs
{
    class SpalartAllmarasAssembly
    {
    public:
        static void resetEffectiveProperties(CBSStateSI& s);
        static void updateEddyViscosity(CBSStateSI& s);
    };
}
