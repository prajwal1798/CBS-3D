//=============================================================================
// CBS3D++_SI
//
// Wide-range numerical robustness sweep for the Spalding wall-law root solve.
//=============================================================================

#include "cbs/turbulence/WallTreatment.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>

int main()
{
    using cbs::Real;
    using cbs::turbulence::WallTreatment;

    const std::array<Real, 12> wall_reynolds_values =
    {
        1.0e-12,
        1.0e-9,
        1.0e-6,
        1.0e-3,
        1.0,
        10.0,
        100.0,
        1.0e3,
        1.0e4,
        1.0e6,
        1.0e8,
        1.0e10
    };

    Real maximum_relative_identity_error = 0.0;
    cbs::Int maximum_iterations = 0;

    for (const Real wall_reynolds : wall_reynolds_values)
    {
        // With y=1 and nu=1, U_t is numerically equal to Re_y.
        const auto result = WallTreatment::evaluateSpalding(
            { wall_reynolds, 0.0, 0.0 },
            { 0.0, 1.0, 0.0 },
            1.0,
            1.0,
            1.0);

        const Real reconstructed =
            result.u_plus * result.y_plus;

        const Real relative_error =
            std::fabs(reconstructed - wall_reynolds) /
            wall_reynolds;

        maximum_relative_identity_error = std::max(
            maximum_relative_identity_error,
            relative_error);

        maximum_iterations = std::max(
            maximum_iterations,
            result.root_iterations);

        if (!std::isfinite(result.u_plus) ||
            !std::isfinite(result.y_plus) ||
            !std::isfinite(result.wall_shear_magnitude) ||
            relative_error > 2.0e-12)
        {
            std::printf(
                "FAIL: Re_y=% .6e u+=% .6e y+=% .6e iterations=%d relerr=% .6e\n",
                wall_reynolds,
                result.u_plus,
                result.y_plus,
                result.root_iterations,
                relative_error);

            return 1;
        }
    }

    std::printf(
        "PASS: Re_y sweep max_identity_error=% .6e max_iterations=%d\n",
        maximum_relative_identity_error,
        maximum_iterations);

    return 0;
}
