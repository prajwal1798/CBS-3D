//=============================================================================
// CBS3D++_SI
//
// Restart-aware distributed production wrapper with optional temperature AFC
// and opt-in SA wall treatment.
//
// The validated production loop is compiled unchanged.  Narrow wrapper classes
// intercept only the established extension points:
//
//   Step 1 momentum RHS  -> replace natural wall flux by Spalding traction
//   strong wall velocity -> restore tangent-space velocity and keep u.n = 0
//   Step 4 energy RHS    -> optional AFC limiter
//=============================================================================

#include "cbs/assembly/EnergyAssembly.hpp"
#include "cbs/assembly/MomentumAssembly.hpp"
#include "cbs/assembly/TemperatureAFC.hpp"
#include "cbs/boundary/Boundary.hpp"
#include "cbs/turbulence/WallModelCoupling.hpp"

namespace cbs
{
    class RestartAfcEnergyAssembly
    {
    public:
        static void assembleStep4Rhs(CBSStateSI& s)
        {
            EnergyAssembly::assembleStep4Rhs(s);

            if (TemperatureAFC::enabled())
            {
                TemperatureAFC::limitDistributedResidual(s);
            }
        }

        static void applyRealTimeEnergyTerm(CBSStateSI& s)
        {
            EnergyAssembly::applyRealTimeEnergyTerm(s);
        }
    };

    class RestartWallModelMomentumAssembly
    {
    public:
        static void assembleStep1Rhs(CBSStateSI& s)
        {
            MomentumAssembly::assembleStep1Rhs(s);
            (void)turbulence::WallModelCoupling::replaceMomentumWallFlux(s);
        }

        static void applyRealTimeMomentumTerm(CBSStateSI& s)
        {
            MomentumAssembly::applyRealTimeMomentumTerm(s);
        }
    };

    class RestartWallModelBoundary
    {
    public:
        static void applyOwnedVelocityConstraints(CBSStateSI& s)
        {
            const auto captured =
                turbulence::WallModelCoupling::captureVelocity(s, false);

            Boundary::applyOwnedVelocityConstraints(s);

            turbulence::WallModelCoupling::restoreTangentialAndEnforceImpermeability(
                s,
                captured);
        }
    };
}

#define EnergyAssembly RestartAfcEnergyAssembly
#define MomentumAssembly RestartWallModelMomentumAssembly
#define Boundary RestartWallModelBoundary
#include "SolverDistributedProductionLoopRestart.cpp"
#undef Boundary
#undef MomentumAssembly
#undef EnergyAssembly
