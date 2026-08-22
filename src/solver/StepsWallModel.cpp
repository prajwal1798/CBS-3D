//=============================================================================
// CBS3D++_SI
//
// Compile the validated Steps implementation with narrow wall-model wrappers.
// The underlying Steps.cpp remains unchanged: only Step-1 wall flux replacement
// and post-Dirichlet normal-only wall projection are intercepted here.
//=============================================================================

#include "cbs/assembly/MomentumAssembly.hpp"
#include "cbs/boundary/Boundary.hpp"
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
        }

        static void applyRealTimeMomentumTerm(CBSStateSI& s)
        {
            MomentumAssembly::applyRealTimeMomentumTerm(s);
        }
    };

    class WallModelBoundary
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

        static void applySymmetry(CBSStateSI& s)
        {
            Boundary::applySymmetry(s);
        }

        static void applyVelocity(CBSStateSI& s)
        {
            const auto captured =
                turbulence::WallModelCoupling::captureVelocity(s, false);

            Boundary::applyVelocity(s);

            turbulence::WallModelCoupling::restoreTangentialAndEnforceImpermeability(
                s,
                captured);
        }

        static void applyOutletBackflowControl(CBSStateSI& s)
        {
            Boundary::applyOutletBackflowControl(s);
            turbulence::WallModelCoupling::enforceImpermeability(s, false);
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
#define Boundary WallModelBoundary
#include "Steps.cpp"
#undef Boundary
#undef MomentumAssembly
