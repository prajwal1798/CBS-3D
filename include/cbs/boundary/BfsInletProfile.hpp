#pragma once

//=============================================================================
// CBS3D++_SI
//
// FILE      : BfsInletProfile.hpp
// PURPOSE   : STREAMWISE INLET VELOCITY u(y) OF BC 507 (BACKWARD-FACING STEP).
//
//   DEFAULT   THE LEGACY POLYNOMIAL OF CBS2D++_SI (UNCHANGED BEHAVIOUR):
//
//               u(y) = 0.6624 y^6 - 7.5547 y^5 + 33.9 y^4 - 75.283 y^3
//                    + 83.368 y^2 - 37.793 y + 2.6959
//
//   TABLE     CBS3D_BFS_INLET_PROFILE=path/to/profile.txt
//               LINES "y u" (# COMMENTS), y STRICTLY INCREASING, AT LEAST
//               TWO ROWS. u IS LINEAR BETWEEN THE ROWS AND 0 OUTSIDE
//               [y_first, y_last]. E.G. THE MEASURED PROFILE OF DENHAM,
//               BRIARD & PATRICK (1975) FOR THE Re = 3025 STEP:
//               cases/backward_facing_step/reference/inlet_denham_station0.txt
//
//             v = w = 0 IN BOTH CASES (SET BY THE CALLER).
//=============================================================================

#include "cbs/core/Types.hpp"

#include <string>
#include <vector>

namespace cbs
{
    class BfsInletProfile
    {
    public:
        // u AT HEIGHT y. THE TABLE (IF ANY) IS READ ON THE FIRST CALL; A
        // MALFORMED TABLE THROWS std::runtime_error.
        static Real u(Real y);

        // TRUE WHEN CBS3D_BFS_INLET_PROFILE IS SET (THE TABLE IS IN USE).
        static bool tabulated();

        // ONE LOG LINE: "TABLE path (n points, y = a .. b, INT u dy = q)"
        // OR "LEGACY POLYNOMIAL".
        static std::string summary();

        // THE LEGACY POLYNOMIAL, ALSO USED WHEN NO TABLE IS SET.
        static Real legacyPolynomial(Real y);

        // LINEAR INTERPOLATION IN A TABLE; 0 OUTSIDE [ys.front(), ys.back()].
        static Real interpolate(const std::vector<Real>& ys, const std::vector<Real>& us, Real y);

        // READS "y u" ROWS; THROWS ON A MISSING FILE OR A BAD ROW.
        static void readTable(const std::string& path, std::vector<Real>& ys, std::vector<Real>& us);
    };
}
