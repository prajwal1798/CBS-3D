#pragma once

//=============================================================================
// CBS3D++_SI
//
// Distributed algebraic flux correction (AFC) for the explicit Step-4
// temperature update.
//
// The high-order thermal residual is first augmented by a symmetric graph-
// viscosity operator that makes every elemental off-diagonal transport
// coefficient non-negative.  The resulting low-order update is then corrected
// by limited pairwise antidiffusive fluxes using a Zalesak-style nodal limiter.
//
// This module operates on the already assembled rank-local high-order residual
// in rhs1.  It performs all owner/ghost reductions required by the limiter and
// returns a globally assembled limited residual on owner nodes.  Ghost rhs1
// entries are zeroed so the established outer reverse-add operation is a no-op
// for the AFC path.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"

namespace cbs
{
    class TemperatureAFC
    {
    public:
        struct Diagnostics
        {
            long long elemental_edges = 0;
            long long limited_edges = 0;
            Real minimum_alpha = 1.0;
            Real minimum_low_order_temperature = 0.0;
            Real minimum_limited_temperature = 0.0;
            Real maximum_limited_temperature = 0.0;
            Real global_correction_balance = 0.0;
        };

        // Environment control:
        //
        //     CBS3D_TEMPERATURE_AFC=1
        //
        // Optional strict physical lower bound:
        //
        //     CBS3D_AFC_MIN_TEMPERATURE=<Kelvin>
        //
        // Optional diagnostic cadence, default 1000 calls:
        //
        //     CBS3D_AFC_LOG_EVERY=<positive integer>
        static bool enabled();

        // Converts the rank-local high-order residual in rhs1 into a
        // distributed AFC-limited residual.  The caller must have already
        // broadcast temperature1 and the corrected velocity to ghosts.
        static Diagnostics limitDistributedResidual(CBSStateSI& s);
    };
}
