//=============================================================================
// CBS3D++_SI
//
// FILE      : BfsInletProfileVerification.cpp
// PURPOSE   : BC 507 INLET PROFILE (BfsInletProfile):
//               1. WITHOUT CBS3D_BFS_INLET_PROFILE, u(y) IS THE LEGACY
//                  CBS2D++_SI POLYNOMIAL (UNCHANGED BEHAVIOUR)
//               2. A "y u" TABLE IS READ AND INTERPOLATED LINEARLY, 0 OUTSIDE
//               3. A NON-INCREASING OR MALFORMED TABLE IS REJECTED
//=============================================================================

#include "cbs/boundary/BfsInletProfile.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using cbs::BfsInletProfile;
using cbs::Real;

namespace
{
    int failures = 0;

    void check(bool ok, const std::string& what)
    {
        std::cout << (ok ? "PASS  " : "FAIL  ") << what << "\n";

        if (!ok)
        {
            ++failures;
        }
    }

    bool close(Real a, Real b, Real tol)
    {
        return std::abs(a - b) <= tol;
    }

    std::string write_file(const std::string& name, const std::string& text)
    {
        std::ofstream out(name);
        out << text;
        return name;
    }
}

int main()
{
    // ================================
    // 1. DEFAULT: LEGACY POLYNOMIAL
    // ================================
#ifdef _WIN32
    _putenv_s("CBS3D_BFS_INLET_PROFILE", "");
#else
    unsetenv("CBS3D_BFS_INLET_PROFILE");
#endif
    check(!BfsInletProfile::tabulated(), "no table without CBS3D_BFS_INLET_PROFILE");

    for (const Real y : {1.0, 1.5, 1.7679, 2.0, 2.5, 3.0})
    {
        const Real direct =
            0.6624 * std::pow(y, 6) - 7.5547 * std::pow(y, 5) + 33.9 * std::pow(y, 4)
            - 75.283 * std::pow(y, 3) + 83.368 * std::pow(y, 2) - 37.793 * y + 2.6959;
        check(BfsInletProfile::u(y) == direct, "u(" + std::to_string(y) + ") = legacy polynomial");
    }

    check(close(BfsInletProfile::legacyPolynomial(1.0), 0.0, 1.0e-2) &&
          close(BfsInletProfile::legacyPolynomial(3.0), 0.0, 2.0e-2),
          "legacy polynomial ~ 0 at both walls y = 1, 3");

    // ================================
    // 2. TABLE: LINEAR INTERPOLATION, 0 OUTSIDE
    // ================================
    const std::string good = write_file(
        "bfs_inlet_profile_good.txt",
        "# y u\n1.0 0.0\n1.5 1.0\n\n2.5 1.2\n3.0 0.0\n");

    std::vector<Real> ys;
    std::vector<Real> us;
    BfsInletProfile::readTable(good, ys, us);
    check(ys.size() == 4 && us.size() == 4, "table: 4 rows read, comment and blank line skipped");
    check(close(BfsInletProfile::interpolate(ys, us, 1.25), 0.5, 1.0e-14), "table: u(1.25) = 0.5");
    check(close(BfsInletProfile::interpolate(ys, us, 2.0), 1.1, 1.0e-14), "table: u(2.0) = 1.1");
    check(close(BfsInletProfile::interpolate(ys, us, 1.5), 1.0, 1.0e-14), "table: u at a row = row value");
    check(BfsInletProfile::interpolate(ys, us, 0.5) == 0.0 &&
          BfsInletProfile::interpolate(ys, us, 3.5) == 0.0, "table: u = 0 outside [y_first, y_last]");
    check(close(BfsInletProfile::interpolate(ys, us, 3.0), 0.0, 1.0e-14), "table: u(y_last) = last value");

    // ================================
    // 3. REJECTED TABLES
    // ================================
    const auto rejects = [](const std::string& file)
    {
        std::vector<Real> a;
        std::vector<Real> b;

        try
        {
            BfsInletProfile::readTable(file, a, b);
        }
        catch (const std::runtime_error&)
        {
            return true;
        }

        return false;
    };

    check(rejects(write_file("bfs_inlet_profile_order.txt", "1.0 0.0\n2.0 1.0\n2.0 1.0\n")),
          "table: repeated y rejected");
    check(rejects(write_file("bfs_inlet_profile_text.txt", "1.0 0.0\n2.0 abc\n")),
          "table: malformed row rejected");
    check(rejects(write_file("bfs_inlet_profile_short.txt", "# only one row\n1.0 0.0\n")),
          "table: fewer than two rows rejected");
    check(rejects("bfs_inlet_profile_missing_file.txt"), "table: missing file rejected");

    for (const char* f : {"bfs_inlet_profile_good.txt", "bfs_inlet_profile_order.txt",
                          "bfs_inlet_profile_text.txt", "bfs_inlet_profile_short.txt"})
    {
        std::remove(f);
    }

    std::cout << (failures == 0 ? "ALL PASSED\n" : "FAILURES: " + std::to_string(failures) + "\n");
    return failures == 0 ? 0 : 1;
}
