#pragma once

//=============================================================================
// CBS3D++_SI
//
// Production coupling between the verified Spalding wall law and the CBS
// momentum equation.  The coupling is opt-in through
//
//     CBS3D_SA_WALL_TREATMENT=1
//
// and is deliberately restricted to pure-fluid SA calculations until the
// thermal/CHT wall treatment has been derived and verified.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"

#include <array>
#include <vector>

namespace cbs
{
    namespace turbulence
    {
        struct CapturedWallVelocity
        {
            std::vector<Int> nodes;
            std::vector<std::array<Real, 3>> values;
        };

        struct WallModelMomentumDiagnostics
        {
            long long local_faces = 0;
            Real local_area = 0.0;
            Real minimum_y_plus = 0.0;
            Real maximum_y_plus = 0.0;
            Real modeled_wall_work = 0.0;

            // Kinematic force-like integrals [m^4/s^2] for the residual form:
            //     integral_Gamma tau_w/rho dA.
            std::array<Real, 3> modeled_surface_load = {0.0, 0.0, 0.0};
            std::array<Real, 3> assembled_nodal_load = {0.0, 0.0, 0.0};
        };

        class WallModelCoupling
        {
        public:
            // True only when the production wall treatment was explicitly
            // requested.  A requested but invalid configuration throws rather
            // than silently falling back to wall-resolved behaviour.
            static bool enabled(const CBSStateSI& s);

            // Capture the unconstrained velocity produced by Step 1 or Step 3
            // before the legacy no-slip package zeros wall nodes.
            static CapturedWallVelocity captureVelocity(
                const CBSStateSI& s,
                bool owned_only);

            // Restore only the tangent-space component of the captured velocity.
            // The normal-space projector is constructed from the span of all
            // incident modelled-wall normals, so edges and corners are handled
            // without averaging normals.
            static void restoreTangentialAndEnforceImpermeability(
                CBSStateSI& s,
                const CapturedWallVelocity& captured);

            // Re-project the current velocity without restoring a previous value.
            static void enforceImpermeability(
                CBSStateSI& s,
                bool owned_only);

            // Replace the natural viscous boundary term already assembled by
            // MomentumAssembly on modelled wall faces:
            //
            //   R <- R - integral N (nu_eff grad(u).n) dA
            //          + integral N (tau_w/rho) dA.
            //
            // This must be called before MPI ghost-to-owner residual summation.
            static WallModelMomentumDiagnostics replaceMomentumWallFlux(
                CBSStateSI& s);

            // Test/audit queries.  They trigger the same validated cache build
            // used by production coupling.
            static long long globalWallFaceCount(const CBSStateSI& s);
            static bool isModelWallNode(const CBSStateSI& s, Int ip);
        };
    }
}
