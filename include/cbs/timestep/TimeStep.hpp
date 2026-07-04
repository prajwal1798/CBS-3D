#pragma once

//=============================================================================
// CBS3D++_SI
//
// Time-step calculation for the three-dimensional semi-implicit CBS solver.
//
// The module supports fixed, local and global pseudo-time stepping. The local
// stability limit is obtained from the advective and diffusive restrictions:
//
//     dt_adv  = h / |u|
//
//     dt_diff = h^2 / (2 D)
//
//     dt_e = C_safety min(dt_adv, dt_diff)
//
// where h is the selected element length, |u| is the characteristic velocity
// and D is the controlling diffusivity.
//
// The parameter ilots selects the time-step treatment:
//
//     ilots <= -1  fixed global time step
//     ilots == 1  retain local nodal and element time steps
//     ilots == 2  retain local time steps with an additional deltr cap
//     otherwise   use the global minimum time step
//
// The parameter htype selects the element length and velocity estimate:
//
//     htype == 1  standard geometric element length
//     htype == 2  nodal/SUPG-style characteristic length
//     htype == 3  geometric length with centroidal velocity
//
// The module also forms the inverse momentum and thermal time diagonals and
// applies the local time-step scaling to the pressure operator.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"

namespace cbs
{
    class TimeStep
    {
    public:
        // Calculates the element and nodal time steps for one CBS iteration.
        static void computeTimeStep(CBSStateSI& s, Int iitime);

        // Calculates the local stability time step at every mesh node.
        static void computeNodalLocalTimeStep(CBSStateSI& s);

        // Calculates one global real-time step from the minimum local value.
        static void computeGlobalRealTimeStep(CBSStateSI& s, Int iitime);

        // Applies the optional pressure-based Step 2 time-step restriction.
        static void applyStep2PressureTimeStepCorrection(CBSStateSI& s);

        // Forms the momentum, thermal and pressure left-hand-side diagonals.
        static void updateLhsDiagonal(CBSStateSI& s);

        // Reserved hook for additional real-time terms.
        static void updateRealTimeTerms(CBSStateSI& s);
    };
}
