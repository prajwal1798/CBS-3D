//=============================================================================
// CBS3D++_SI
//
// Restart-aware distributed production wrapper with optional temperature AFC.
//
// The validated production loop is compiled unchanged.  Only the Step-4
// EnergyAssembly class name is substituted so that the established high-order
// residual can be converted to an AFC-limited residual when
// CBS3D_TEMPERATURE_AFC=1 is present in the environment.
//=============================================================================

#include "cbs/assembly/EnergyAssembly.hpp"
#include "cbs/assembly/TemperatureAFC.hpp"

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
}

#define EnergyAssembly RestartAfcEnergyAssembly
#include "SolverDistributedProductionLoopRestart.cpp"
#undef EnergyAssembly
