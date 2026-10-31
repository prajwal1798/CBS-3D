#pragma once

//=============================================================================
// CBS3D++_SI
//
// Depth-dependent volumetric heat source (nuclear heating in a blanket).
//
// Neutron and gamma heating in a fusion blanket is largest just behind the
// first wall and decays roughly exponentially with depth. For every material
// m listed in the optional .par block the element source is
//
//     Q_e = A_m exp(-d_e / lambda_m)                        [W/m^3]
//
//     d_e = max(0, (x_c,e - x_0) . n)                       [m]
//
// where x_c,e is the element centroid, x_0 any point on the first-wall plane
// and n the unit normal from that plane into the blanket. The value is stored
// in Qvol_e once, before the time loop, so the Step-4 assembly and its cost
// are unchanged.
//
// Element sizes in the blanket meshes are 0.75-12.5 mm and lambda is about
// 150 mm, so one centroid value per P1 element is accurate to a fraction of a
// percent where the heating is strongest.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"

namespace cbs
{
    class NuclearHeating
    {
    public:
        // Fill Qvol_e from the profiles (called after the .matprop is read).
        // Rejects a profile on a material that also has a constant Qvol or,
        // for solids, a non-zero .par source_solid.
        static void applyProfiles(CBSStateSI& s);

        // Print the total heat input of the whole (all-rank) domain:
        // volumetric source per material and prescribed BC 532 flux.
        // Needs detJ and face_norm, i.e. call after shape derivatives and
        // boundary normals. Collective in MPI runs.
        static void reportHeatInput(const CBSStateSI& s);

        // Q = A exp(-d / lambda) at point x for one profile.
        static Real profileValue(
            const RunConfig& cfg,
            const NuclearHeatingProfile& profile,
            Real x,
            Real y,
            Real z);
    };
}
