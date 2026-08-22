//=============================================================================
// CBS3D++_SI
// Algebraic verification of the continuous wall-treatment kernel.
//=============================================================================

#include "cbs/turbulence/WallTreatment.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <exception>

namespace
{
    using cbs::Real;
    using cbs::turbulence::WallTreatment;

    Real dot(const std::array<Real, 3>& a, const std::array<Real, 3>& b)
    {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    Real magnitude(const std::array<Real, 3>& value)
    {
        return std::sqrt(dot(value, value));
    }

    bool nearly_equal(
        Real a,
        Real b,
        Real relative_tolerance,
        Real absolute_tolerance = 0.0)
    {
        return std::fabs(a - b) <=
            absolute_tolerance +
            relative_tolerance * std::max(std::fabs(a), std::fabs(b));
    }
}

int main()
{
    bool failed = false;

    try
    {
        const std::array<Real, 3> tangent_positive =
            WallTreatment::tangentialVelocity(
                { 3.0, 4.0, 5.0 }, { 0.0, 0.0, 2.0 });
        const std::array<Real, 3> tangent_negative =
            WallTreatment::tangentialVelocity(
                { 3.0, 4.0, 5.0 }, { 0.0, 0.0, -7.0 });

        if (!nearly_equal(tangent_positive[0], 3.0, 1.0e-14) ||
            !nearly_equal(tangent_positive[1], 4.0, 1.0e-14) ||
            !nearly_equal(tangent_positive[2], 0.0, 1.0e-14) ||
            magnitude({
                tangent_positive[0] - tangent_negative[0],
                tangent_positive[1] - tangent_negative[1],
                tangent_positive[2] - tangent_negative[2] }) > 1.0e-14)
        {
            std::printf("FAIL: tangential projection\n");
            failed = true;
        }

        // Viscous-sublayer asymptote: tau_w -> rho nu U_t / y.
        const Real rho = 4.0;
        const Real nu = 1.0e-5;
        const Real y = 1.0e-6;
        const Real velocity = 1.0e-2;
        const auto viscous = WallTreatment::evaluateSpalding(
            { velocity, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, y, rho, nu);
        const Real analytic_viscous_shear = rho * nu * velocity / y;

        if (!nearly_equal(
                viscous.wall_shear_magnitude,
                analytic_viscous_shear,
                2.0e-6))
        {
            std::printf(
                "FAIL: viscous asymptote tau=% .17e expected=% .17e\n",
                viscous.wall_shear_magnitude,
                analytic_viscous_shear);
            failed = true;
        }

        // Representative wall-function range and exact Re_y identity.
        const auto log_layer = WallTreatment::evaluateSpalding(
            { 10.5, 0.0, 0.0 },
            { 0.0, 1.0, 0.0 },
            2.0e-4,
            1.0,
            3.0e-6);
        const Real wall_reynolds =
            log_layer.tangential_speed * 2.0e-4 / 3.0e-6;
        const Real reconstructed_wall_reynolds =
            log_layer.u_plus * log_layer.y_plus;

        if (!nearly_equal(
                reconstructed_wall_reynolds,
                wall_reynolds,
                2.0e-12) ||
            log_layer.y_plus < 30.0 ||
            log_layer.y_plus > 100.0)
        {
            std::printf(
                "FAIL: wall-unit solve Re_y=% .17e reconstructed=% .17e y+=% .17e\n",
                wall_reynolds,
                reconstructed_wall_reynolds,
                log_layer.y_plus);
            failed = true;
        }

        // Wall shear must be tangent and dissipative.
        const std::array<Real, 3> normal = { 0.0, 1.0, 0.0 };

        if (std::fabs(dot(log_layer.wall_shear_on_fluid, normal)) >
            1.0e-13 * std::max(1.0, log_layer.wall_shear_magnitude))
        {
            std::printf("FAIL: wall shear has a normal component\n");
            failed = true;
        }

        if (!(dot(
                log_layer.wall_shear_on_fluid,
                log_layer.tangential_velocity) < 0.0))
        {
            std::printf("FAIL: wall shear is not dissipative\n");
            failed = true;
        }

        // Geometric normal orientation must not change shear.
        const auto flipped = WallTreatment::evaluateSpalding(
            { 10.5, 0.0, 0.0 },
            { 0.0, -3.0, 0.0 },
            2.0e-4,
            1.0,
            3.0e-6);

        if (magnitude({
                log_layer.wall_shear_on_fluid[0] - flipped.wall_shear_on_fluid[0],
                log_layer.wall_shear_on_fluid[1] - flipped.wall_shear_on_fluid[1],
                log_layer.wall_shear_on_fluid[2] - flipped.wall_shear_on_fluid[2] }) > 1.0e-12)
        {
            std::printf("FAIL: normal-orientation invariance\n");
            failed = true;
        }

        // Purely normal sample has zero modeled wall shear.
        const auto zero_tangent = WallTreatment::evaluateSpalding(
            { 0.0, 0.0, 3.0 },
            { 0.0, 0.0, 1.0 },
            1.0e-3,
            1.0,
            1.0e-5);

        if (zero_tangent.wall_shear_magnitude != 0.0 ||
            zero_tangent.friction_velocity != 0.0)
        {
            std::printf("FAIL: zero tangential velocity\n");
            failed = true;
        }

        bool invalid_normal_rejected = false;
        try
        {
            (void)WallTreatment::evaluateSpalding(
                { 1.0, 0.0, 0.0 },
                { 0.0, 0.0, 0.0 },
                1.0e-3,
                1.0,
                1.0e-5);
        }
        catch (const std::exception&)
        {
            invalid_normal_rejected = true;
        }

        if (!invalid_normal_rejected)
        {
            std::printf("FAIL: zero wall normal was accepted\n");
            failed = true;
        }
    }
    catch (const std::exception& error)
    {
        std::printf("FAIL: unexpected exception: %s\n", error.what());
        failed = true;
    }

    if (!failed)
    {
        std::printf("PASS: wall-treatment algebraic verification\n");
    }

    return failed ? 1 : 0;
}
