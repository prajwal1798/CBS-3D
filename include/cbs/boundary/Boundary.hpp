#pragma once

//=============================================================================
// CBS3D++_SI
//
// Boundary conditions for the three-dimensional semi-implicit CBS solver.
//
// The module applies strong nodal conditions for velocity, pressure and
// temperature. It also applies symmetry and outlet-backflow corrections.
//
// Boundary-face storage:
//
//     iside(1:3, ib) = nodes of triangular boundary face ib
//     iside(4, ib)   = local face number in the parent tetrahedron
//     iside(5, ib)   = parent tetrahedron
//     iside(6, ib)   = solver boundary-condition identifier
//
// The outward area-weighted normal is stored as:
//
//     face_norm(1:3, ib) = A_f n
//     face_norm(4, ib)   = A_f
//
// where A_f is the face area and n is the outward unit normal.
//
// Supported boundary identifiers:
//
//     500  prescribed velocity, adiabatic
//     501  T = 1, no-slip
//     502  T = 0, no-slip
//     503  prescribed x-velocity, T = 0
//     504  pressure boundary
//     506  symmetry / slip
//     507  backward-facing-step parabolic inlet
//     508  rectangular-channel parabolic inlet
//     510  prescribed velocity and prescribed temperature
//     511  mass-flow inlet and prescribed temperature
//     520  pressure outlet
//     530  no-slip adiabatic wall
//     532  no-slip prescribed-heat-flux wall
//     901  conformal fluid-solid interface
//     902  legacy heat-flux marker
//
// BC 902 is not a fluid-solid interface condition.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"

namespace cbs
{
    class Boundary
    {
    public:
        // Applies prescribed nodal temperature values.
        static void applyTemperature(CBSStateSI& s);

        // Applies prescribed nodal velocity values and no-slip conditions.
        static void applyVelocity(CBSStateSI& s);

        // Applies the persistent strong velocity state only at MPI-owned nodes.
        // The owner values must subsequently be broadcast to ghost copies.
        static void applyOwnedVelocityConstraints(CBSStateSI& s);

        // Applies prescribed nodal pressure values.
        static void applyPressure(CBSStateSI& s);

        // Removes the normal velocity component on symmetry boundaries.
        static void applySymmetry(CBSStateSI& s);

        // Removes inward normal velocity on pressure-outlet boundaries.
        static void applyOutletBackflowControl(CBSStateSI& s);
    };
}