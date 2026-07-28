#pragma once

//=============================================================================
// CBS3D++_SI
//
// Bounded algebraic-flux-correction update for CBS Step 4.
//
// The existing Galerkin/characteristic residual remains the high-order target.
// A conservative element-edge diffusion operator first creates a monotone
// low-order update. A local-extremum-preserving Zalesak limiter then restores
// as much of the removed high-order transport as the nodal bounds permit.
//
// The class operates on an already assembled high-order rhs1. This keeps the
// established EnergyAssembly source, heat-flux and material treatment intact.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"

#ifdef CBS3D_USE_MPI
#include <mpi.h>
#endif

namespace cbs
{
    class ThermalAfc
    {
    public:
        // Applies the bounded AFC update to a serial/OpenMP calculation.
        // EnergyAssembly::assembleStep4Rhs() must have populated rhs1 first.
        static void applySerial(CBSStateSI& s);

#ifdef CBS3D_USE_MPI
        // Applies the same AFC update to a partitioned calculation.
        //
        // rhs1 contains rank-local owned-element contributions. The method
        // performs all required reverse/forward halo operations, adds the
        // owner-assembled external heat-flux load, advances owner nodes and
        // broadcasts the corrected temperature to ghosts.
        static void applyDistributed(
            CBSStateSI& s,
            const Array1D<Int>& fixed_temperature,
            const Array1D<Real>& fixed_value,
            const Array1D<Real>& external_heat_flux_load,
            MPI_Comm communicator = MPI_COMM_WORLD);
#endif
    };
}
