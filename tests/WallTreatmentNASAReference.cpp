//=============================================================================
// CBS3D++_SI
//
// Independent reference-point check for the quartic Spalding wall law used by
// this foundation. NASA/TP-2006-214497 reproduces the Spalding velocity profile
// as
//
//   y+ = u+ + (1/E) [ exp(kappa u+) - 1 - kappa u+
//                     - (kappa u+)^2/2!
//                     - (kappa u+)^3/3!
//                     - (kappa u+)^4/4! ]
//
// and records the original values kappa=0.40, B=5.5 and 1/E=0.1108.
// This test deliberately uses a hard-coded value evaluated independently from
// that published formula. It therefore detects an accidental change to the
// cubic-truncation variant used by some other software packages.
//=============================================================================

#include "cbs/turbulence/WallTreatment.hpp"

#include <cmath>
#include <cstdio>

int main()
{
    cbs::turbulence::WallTreatmentOptions options;
    options.kappa = 0.40;
    options.log_intercept = 5.5;

    const cbs::Real u_plus = 10.0;
    const cbs::Real expected_y_plus = 12.245405693972817;

    const cbs::Real actual_y_plus =
        cbs::turbulence::WallTreatment::spaldingYPlus(u_plus, options);

    const cbs::Real error =
        std::fabs(actual_y_plus - expected_y_plus);

    if (error > 2.0e-13)
    {
        std::printf(
            "FAIL: NASA Spalding reference u+=%.17e y+=%.17e expected=%.17e error=%.3e\n",
            u_plus,
            actual_y_plus,
            expected_y_plus,
            error);

        return 1;
    }

    std::printf(
        "PASS: NASA Spalding quartic reference u+=%.6f y+=%.15f\n",
        u_plus,
        actual_y_plus);

    return 0;
}
