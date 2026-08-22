#pragma once

//=============================================================================
// CBS3D++_SI
//
// Algebraic foundation for wall-modelled momentum boundary conditions.
//
// This module deliberately has no mesh, MPI or CBS-stage dependency. It turns
// one off-wall velocity sample into a wall shear traction. Keeping this kernel
// independent lets the law-of-the-wall mathematics be verified before any
// existing wall boundary semantics are changed.
//=============================================================================

#include "cbs/core/Types.hpp"

#include <array>

namespace cbs
{
    namespace turbulence
    {
        struct WallTreatmentOptions
        {
            // Smooth-wall Spalding constants.
            Real kappa = 0.41;
            Real log_intercept = 5.2;

            // Safeguarded Newton/bisection convergence controls for u+.
            Real root_relative_tolerance = 1.0e-13;
            Int root_max_iterations = 80;

            // Below this tangential speed the wall shear is exactly zero.
            Real tangential_speed_floor = 1.0e-14;
        };

        struct WallTreatmentResult
        {
            std::array<Real, 3> tangential_velocity = { 0.0, 0.0, 0.0 };
            std::array<Real, 3> unit_tangent = { 0.0, 0.0, 0.0 };

            // Dynamic traction exerted by the wall on the fluid [Pa].
            // It opposes the sampled tangential velocity.
            std::array<Real, 3> wall_shear_on_fluid = { 0.0, 0.0, 0.0 };

            // wall_shear_on_fluid / rho [m^2/s^2]. This has the same units as
            // nu * du/dn in the current kinematic Step-1 diffusion weak form.
            std::array<Real, 3> kinematic_wall_traction = { 0.0, 0.0, 0.0 };

            Real tangential_speed = 0.0;
            Real friction_velocity = 0.0;
            Real u_plus = 0.0;
            Real y_plus = 0.0;
            Real wall_shear_magnitude = 0.0;
            Int root_iterations = 0;
        };

        class WallTreatment
        {
        public:
            // Quartic-subtraction Spalding smooth-wall relation:
            //
            // y+ = u+ + exp(-kappa B)
            //      [exp(kappa u+) - sum_{n=0}^4 (kappa u+)^n/n!].
            //
            // The implementation evaluates the exponential remainder with a
            // cancellation-safe series near the viscous sublayer.
            static Real spaldingYPlus(
                Real u_plus,
                const WallTreatmentOptions& options = WallTreatmentOptions());

            // Removes the velocity component normal to the supplied wall normal.
            static std::array<Real, 3> tangentialVelocity(
                const std::array<Real, 3>& velocity,
                const std::array<Real, 3>& wall_normal);

            // Evaluates the continuous wall law from one off-wall sample.
            //
            // Re_y = |u_t| y / nu = u+ y+
            //
            // is used to solve a monotone scalar equation for u+. The returned
            // traction is purely tangential and dissipative:
            //
            // wall_shear_on_fluid . u_t <= 0.
            static WallTreatmentResult evaluateSpalding(
                const std::array<Real, 3>& sample_velocity,
                const std::array<Real, 3>& wall_normal,
                Real wall_distance,
                Real density,
                Real molecular_nu,
                const WallTreatmentOptions& options = WallTreatmentOptions());
        };
    }
}
