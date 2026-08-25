//=============================================================================
// CBS3D++_SI
//
// Restart-aware distributed production wrapper with optional temperature AFC,
// pure-fluid Spalding treatment and conformal CHT Spalding/Kader treatment.
//
// The validated production loop remains unchanged.  Narrow wrappers intercept
// only the established extension points:
//
//   SA preprocessing       -> prepare CHT thermal stability conductivity
//   Step 1 momentum RHS    -> add the selected wall-model momentum coupling
//   strong wall velocity   -> restore tangent-space velocity, keep u.n = 0
//   Step 4 energy RHS      -> CHT Kader normal-diffusion replacement OR AFC
//
// Temperature AFC is intentionally not combined with the CHT wall operator yet.
// AFC reconstructs an isotropic elemental operator from k_eff_e; the CHT wall
// model replaces that operator by an anisotropic rank-one normal correction.
// Running both without deriving the same operator inside AFC would make the
// limiter inconsistent with the residual it is limiting.
//=============================================================================

#include "cbs/assembly/EnergyAssembly.hpp"
#include "cbs/assembly/MomentumAssembly.hpp"
#include "cbs/assembly/TemperatureAFC.hpp"
#include "cbs/boundary/Boundary.hpp"
#include "cbs/turbulence/CHTWallModelCoupling.hpp"
#include "cbs/turbulence/TurbulencePreprocess.hpp"
#include "cbs/turbulence/WallModelCoupling.hpp"

#include <stdexcept>

namespace cbs
{
    class RestartCHTTurbulencePreprocess
    {
    public:
        static void prepareSpalartAllmaras(CBSStateSI& s)
        {
            TurbulencePreprocess::prepareSpalartAllmaras(s);

            // Fail before the first pressure/timestep setup if an unsupported
            // combination was requested.  Silently applying AFC to a different
            // thermal operator would be numerically worse than stopping here.
            if (turbulence::CHTWallModelCoupling::enabled(s) &&
                TemperatureAFC::enabled())
            {
                throw std::runtime_error(
                    "CBS3D CHT wall treatment cannot currently be combined with CBS3D_TEMPERATURE_AFC: AFC reconstructs an isotropic k_eff operator while the Kader wall treatment uses anisotropic wall-normal conductivity");
            }

            // This is the critical first-step stability hook.  SA preprocessing
            // has just rebuilt nu_t_e/k_eff_e, so wall-adjacent conductivity can
            // now be raised to a conservative Kader normal-diffusion bound before
            // TimeStep::computeTimeStep is called for the first time.
            turbulence::CHTWallModelCoupling::prepareThermalStabilityConductivity(s);
        }
    };

    class RestartAfcEnergyAssembly
    {
    public:
        static void assembleStep4Rhs(CBSStateSI& s)
        {
            // SA has updated nu_t_e immediately before Step 4. Rebuild the exact
            // ordinary bulk turbulent conductivity and install the wall-normal
            // stability bound for both this residual and the next timestep.
            turbulence::CHTWallModelCoupling::prepareThermalStabilityConductivity(s);

            EnergyAssembly::assembleStep4Rhs(s);

            // EnergyAssembly saw the conservative isotropic upper bound. Replace
            // it by the intended tensor: ordinary SA k in tangent directions and
            // Kader k_n in each incident CHT wall-normal direction.
            (void)turbulence::CHTWallModelCoupling::correctThermalWallDiffusion(s);

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

            // The two wall modes are intentionally independent.  The CHT module
            // rejects simultaneous activation so an accidental double wall load
            // cannot pass silently.
            (void)turbulence::WallModelCoupling::replaceMomentumWallFlux(s);
            (void)turbulence::CHTWallModelCoupling::addMomentumWallFlux(s);
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
            const auto captured_pure =
                turbulence::WallModelCoupling::captureVelocity(s, false);
            const auto captured_cht =
                turbulence::CHTWallModelCoupling::captureVelocity(s, false);

            Boundary::applyOwnedVelocityConstraints(s);

            turbulence::WallModelCoupling::restoreTangentialAndEnforceImpermeability(
                s,
                captured_pure);
            turbulence::CHTWallModelCoupling::restoreTangentialAndEnforceImpermeability(
                s,
                captured_cht);
        }
    };
}

#define TurbulencePreprocess RestartCHTTurbulencePreprocess
#define EnergyAssembly RestartAfcEnergyAssembly
#define MomentumAssembly RestartWallModelMomentumAssembly
#define Boundary RestartWallModelBoundary
#include "SolverDistributedProductionLoopRestart.cpp"
#undef Boundary
#undef MomentumAssembly
#undef EnergyAssembly
#undef TurbulencePreprocess