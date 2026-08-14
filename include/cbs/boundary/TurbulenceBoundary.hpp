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

        // Makes the SA node classification identical on every rank that holds a
        // copy of a shared node.
        //
        // classifyNodes derives the flags from the locally visible element and
        // boundary-face lists.  A node lying on a physical wall whose wall face
        // is owned by a neighbouring rank is therefore not marked as a wall node
        // here, so the same node would be a wall node on one rank and an
        // ordinary fluid node on another.  The two ranks would then disagree on
        // whether to hold nu_tilde at zero and on whether to floor its wall
        // distance, and the disagreement would enter the halo exchange every
        // step.  This routine performs a logical-OR reduction of the three flags
        // over all copies of each shared node.  It is a no-op in serial.
        static void synchroniseClassification(CBSStateSI& s);
        static void initialiseNuTilde(CBSStateSI& s);
        static void applyWallValues(CBSStateSI& s);
        static void applyInletValues(CBSStateSI& s);
    };
}
