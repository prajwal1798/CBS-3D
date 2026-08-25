//=============================================================================
// CBS3D++_SI
//
// Compile the validated Steps implementation with narrow wall-model wrappers.
// The underlying Steps.cpp remains unchanged.  Pure-fluid explicit-wall
// Spalding coupling and conformal CHT Spalding/Kader coupling are opt-in and
// remain separate modules.
//=============================================================================

#include "cbs/assembly/EnergyAssembly.hpp"
#include "cbs/assembly/MomentumAssembly.hpp"
#include "cbs/boundary/Boundary.hpp"
#include "cbs/turbulence/CHTWallModelCoupling.hpp"
#include "cbs/turbulence/WallModelCoupling.hpp"

namespace cbs
{
    class WallModelMomentumAssembly
    {
    public:
        static void assembleStep1Rhs(CBSStateSI& s)
        {
            MomentumAssembly::assembleStep1Rhs(s);
            (void)turbulence::WallModelCoupling::replaceMomentumWallFlux(s);
            (void)turbulence::CHTWallModelCoupling::addMomentumWallFlux(s);
        }

        static void applyRealTimeMomentumTerm(CBSStateSI& s)
        {
            MomentumAssembly::applyRealTimeMomentumTerm(s);
        }
    };

    class WallModelEnergyAssembly
    {
    public:
        static void assembleStep4Rhs(CBSStateSI& s)
        {
            // The preparation reconstructs the ordinary SA bulk conductivity
            // and installs a conservative isotropic upper bound for explicit
            // stability.  The subsequent correction removes that artificial
            // tangential inflation and installs the intended Kader wall-normal
            // conductivity tensor.
            turbulence::CHTWallModelCoupling::prepareThermalStabilityConductivity(s);
            EnergyAssembly::assembleStep4Rhs(s);
            (void)turbulence::CHTWallModelCoupling::correctThermalWallDiffusion(s);
        }
    };

    class WallModelBoundary
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

        static void applySymmetry(CBSStateSI& s)
        {
            Boundary::applySymmetry(s);
        }

        static void applyVelocity(CBSStateSI& s)
        {
            const auto captured_pure =
                turbulence::WallModelCoupling::captureVelocity(s, false);
            const auto captured_cht =
                turbulence::CHTWallModelCoupling::captureVelocity(s, false);

            Boundary::applyVelocity(s);

            turbulence::WallModelCoupling::restoreTangentialAndEnforceImpermeability(
                s,
                captured_pure);
            turbulence::CHTWallModelCoupling::restoreTangentialAndEnforceImpermeability(
                s,
                captured_cht);
        }

        static void applyOutletBackflowControl(CBSStateSI& s)
        {
            Boundary::applyOutletBackflowControl(s);
            turbulence::WallModelCoupling::enforceImpermeability(s, false);
            turbulence::CHTWallModelCoupling::enforceImpermeability(s, false);
        }

        static void applyPressure(CBSStateSI& s)
        {
            Boundary::applyPressure(s);
        }

        static void applyTemperature(CBSStateSI& s)
        {
            Boundary::applyTemperature(s);
        }
    };
}

#define MomentumAssembly WallModelMomentumAssembly
#define EnergyAssembly WallModelEnergyAssembly
#define Boundary WallModelBoundary
#include "Steps.cpp"
#undef Boundary
#undef EnergyAssembly
#undef MomentumAssembly
