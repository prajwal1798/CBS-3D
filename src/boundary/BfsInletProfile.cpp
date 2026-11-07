//=============================================================================
// CBS3D++_SI
//
// FILE      : BfsInletProfile.cpp
// PURPOSE   : BC 507 INLET PROFILE: LEGACY POLYNOMIAL OR A "y u" TABLE,
//             SEE BfsInletProfile.hpp.
//=============================================================================

#include "cbs/boundary/BfsInletProfile.hpp"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace cbs
{
    namespace
    {
        // ================================
        // THE TABLE OF CBS3D_BFS_INLET_PROFILE, READ ONCE (THREAD-SAFE
        // FUNCTION-LOCAL STATIC). EVERY RANK READS THE SAME SMALL FILE.
        // ================================
        struct Table
        {
            bool active = false;
            std::string path;
            std::vector<Real> ys;
            std::vector<Real> us;
        };

        const Table& table()
        {
            static const Table t = []()
            {
                Table r;
                const char* file = std::getenv("CBS3D_BFS_INLET_PROFILE");

                if (file == nullptr || std::string(file).empty())
                {
                    return r;
                }

                r.path = file;
                BfsInletProfile::readTable(r.path, r.ys, r.us);
                r.active = true;
                return r;
            }();

            return t;
        }
    }

    Real BfsInletProfile::legacyPolynomial(Real y)
    {
        return
            0.6624 * std::pow(y, 6)
            - 7.5547 * std::pow(y, 5)
            + 33.9 * std::pow(y, 4)
            - 75.283 * std::pow(y, 3)
            + 83.368 * std::pow(y, 2)
            - 37.793 * y
            + 2.6959;
    }

    Real BfsInletProfile::interpolate(const std::vector<Real>& ys, const std::vector<Real>& us, Real y)
    {
        if (ys.size() < 2 || y < ys.front() || y > ys.back())
        {
            return 0.0;
        }

        std::size_t hi = 1;

        while (hi + 1 < ys.size() && ys[hi] < y)
        {
            ++hi;
        }

        const std::size_t lo = hi - 1;
        const Real w = (y - ys[lo]) / (ys[hi] - ys[lo]);
        return (1.0 - w) * us[lo] + w * us[hi];
    }

    void BfsInletProfile::readTable(const std::string& path, std::vector<Real>& ys, std::vector<Real>& us)
    {
        ys.clear();
        us.clear();

        std::ifstream in(path);

        if (!in)
        {
            throw std::runtime_error("BfsInletProfile - cannot open CBS3D_BFS_INLET_PROFILE " + path);
        }

        std::string line;
        int line_no = 0;

        while (std::getline(in, line))
        {
            ++line_no;
            const std::size_t first = line.find_first_not_of(" \t\r");

            if (first == std::string::npos || line[first] == '#')
            {
                continue;
            }

            std::istringstream row(line);
            Real y = 0.0;
            Real u = 0.0;

            if (!(row >> y >> u) || !std::isfinite(y) || !std::isfinite(u))
            {
                throw std::runtime_error(
                    "BfsInletProfile - " + path + " line " + std::to_string(line_no) +
                    ": expected two numbers 'y u'");
            }

            if (!ys.empty() && !(y > ys.back()))
            {
                throw std::runtime_error(
                    "BfsInletProfile - " + path + " line " + std::to_string(line_no) +
                    ": y must increase strictly");
            }

            ys.push_back(y);
            us.push_back(u);
        }

        if (ys.size() < 2)
        {
            throw std::runtime_error("BfsInletProfile - " + path + " needs at least two 'y u' rows");
        }
    }

    bool BfsInletProfile::tabulated()
    {
        return table().active;
    }

    std::string BfsInletProfile::summary()
    {
        const Table& t = table();

        if (!t.active)
        {
            return "LEGACY POLYNOMIAL";
        }

        Real flux = 0.0;

        for (std::size_t i = 1; i < t.ys.size(); ++i)
        {
            flux += 0.5 * (t.us[i] + t.us[i - 1]) * (t.ys[i] - t.ys[i - 1]);
        }

        std::ostringstream line;
        line << "TABLE " << t.path << " (" << t.ys.size() << " points, y = "
             << t.ys.front() << " .. " << t.ys.back() << ", INT u dy = " << flux << ")";
        return line.str();
    }

    Real BfsInletProfile::u(Real y)
    {
        const Table& t = table();
        return t.active ? interpolate(t.ys, t.us, y) : legacyPolynomial(y);
    }
}
