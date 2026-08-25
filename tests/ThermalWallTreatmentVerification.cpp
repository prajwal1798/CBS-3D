#include "cbs/turbulence/ThermalWallTreatment.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>

namespace
{
    using cbs::Real;
    using cbs::turbulence::ThermalWallTreatment;

    bool close(
        const Real a,
        const Real b,
        const Real relative = 2.0e-11,
        const Real absolute = 2.0e-13)
    {
        return std::abs(a - b) <=
            absolute + relative * std::max(std::abs(a), std::abs(b));
    }
}

int main()
{
    try
    {
        // Independent Kader reference point evaluated directly from the
        // published formula.
        const Real tplus = ThermalWallTreatment::kaderTPlus(74.0, 0.66);

        if (!close(tplus, 12.459500372920848, 2.0e-12, 2.0e-13))
        {
            std::printf("FAIL: Kader reference T+ = %.17e\n", tplus);
            return 1;
        }

        // Conductive sublayer limit: k_n must tend to molecular k.
        const Real rho = 6.7;
        const Real cp = 5200.0;
        const Real mu = 3.1e-5;
        const Real k = 0.24;
        const Real y = 1.0e-7;
        const Real u_tau = 0.02;

        const auto conductive = ThermalWallTreatment::evaluateKader(
            u_tau, y, rho, cp, mu, k);

        if (!close(
                conductive.wall_normal_conductivity,
                k,
                2.0e-7,
                2.0e-10))
        {
            std::printf(
                "FAIL: conductive limit k_n=%.17e k=%.17e\n",
                conductive.wall_normal_conductivity,
                k);
            return 1;
        }

        // Exact zero-shear startup/stagnation limit.
        const auto zero = ThermalWallTreatment::evaluateKader(
            0.0, 1.5e-4, rho, cp, mu, k);

        if (!close(zero.wall_normal_conductivity, k) ||
            !close(zero.y_plus, 0.0) ||
            !close(zero.t_plus, 0.0))
        {
            std::printf("FAIL: zero-shear thermal limit\n");
            return 1;
        }

        // Helium-like wall-function state must remain finite, positive and
        // enhance the wall-normal transport above molecular conduction.
        const auto helium = ThermalWallTreatment::evaluateKader(
            2.3, 1.5e-4, rho, cp, mu, k);

        if (!(helium.y_plus > 0.0) ||
            !(helium.t_plus > 0.0) ||
            !(helium.wall_normal_conductivity > k) ||
            !std::isfinite(helium.wall_normal_conductivity))
        {
            std::printf("FAIL: helium-like Kader state\n");
            return 1;
        }
    }
    catch (const std::exception& error)
    {
        std::printf("FAIL: %s\n", error.what());
        return 1;
    }

    std::printf("PASS: Kader thermal wall treatment\n");
    return 0;
}
