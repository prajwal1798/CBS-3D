#pragma once

//=============================================================================
// CBS3D++_SI
//
// Coarse-wall treatment for a conformal fluid/solid CHT interface.
//
// Runtime activation:
//
//     CBS3D_CHT_WALL_TREATMENT=1
//
// Momentum uses the existing continuous Spalding wall law on the fluid side of
// every reconstructed fluid/solid tetrahedral face.  Strong material-interface
// no-slip is converted to strong impermeability plus weak tangential traction.
//
// Energy uses a continuous Kader temperature law.  Since fluid and solid share
// the same interface temperature DOF, the unresolved fluid-side thermal layer
// is represented by replacing only the first fluid tetrahedron's wall-normal
// conductivity.  Solid conduction and tangential fluid diffusion are retained.
// This avoids adding an artificial interface heat source or double-counting the
// ordinary first-cell diffusion operator.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"

#include <array>
#include <vector>

namespace cbs
{
    namespace turbulence
    {
        struct CHTCapturedWallVelocity
        {
            std::vector<Int> nodes;
            std::vector<std::array<Real, 3>> values;
        };

        struct CHTWallMomentumDiagnostics
        {
            long long local_faces = 0;
            Real local_area = 0.0;
            Real minimum_y_plus = 0.0;
            Real maximum_y_plus = 0.0;
            Real modeled_wall_work = 0.0;
            std::array<Real, 3> modeled_surface_load = {0.0, 0.0, 0.0};
            std::array<Real, 3> assembled_nodal_load = {0.0, 0.0, 0.0};
        };

        struct CHTWallThermalDiagnostics
        {
            long long local_faces = 0;
            Real local_area = 0.0;
            Real minimum_y_plus = 0.0;
            Real maximum_y_plus = 0.0;
            Real minimum_prandtl = 0.0;
            Real maximum_prandtl = 0.0;
            Real minimum_kn_over_k = 0.0;
            Real maximum_kn_over_k = 0.0;
            Real residual_correction_sum = 0.0;
        };

        class CHTWallModelCoupling
        {
        public:
            static bool enabled(const CBSStateSI& s);
            static long long globalWallFaceCount(const CBSStateSI& s);
            static bool isModelWallNode(const CBSStateSI& s, Int ip);

            static CHTCapturedWallVelocity captureVelocity(
                const CBSStateSI& s,
                bool owned_only);

            static void restoreTangentialAndEnforceImpermeability(
                CBSStateSI& s,
                const CHTCapturedWallVelocity& captured);

            static void enforceImpermeability(
                CBSStateSI& s,
                bool owned_only);

            static CHTWallMomentumDiagnostics addMomentumWallFlux(
                CBSStateSI& s);

            // Reconstruct ordinary SA turbulent conductivity and place a
            // conservative wall-normal stability bound in k_eff_e before the
            // explicit timestep/energy stage.
            static void prepareThermalStabilityConductivity(
                CBSStateSI& s);

            // Replace the temporary isotropic stability conductivity by the
            // intended anisotropic Kader wall-normal correction in rhs1.
            static CHTWallThermalDiagnostics correctThermalWallDiffusion(
                CBSStateSI& s);
        };
    }
}
