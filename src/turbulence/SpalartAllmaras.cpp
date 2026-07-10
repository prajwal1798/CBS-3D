#include "cbs/turbulence/SpalartAllmaras.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cbs::turbulence
{
    namespace
    {
        constexpr Real tiny = 1.0e-30;

        [[nodiscard]] Real safe_positive(const Real value, const char* name)
        {
            if (!std::isfinite(value) || value <= 0.0)
            {
                throw std::runtime_error(
                    std::string("SpalartAllmaras - non-positive or non-finite ") + name);
            }
            return value;
        }
    }

    Real SpalartAllmarasConstants::cw1() const noexcept
    {
        return cb1 / (kappa * kappa) + (1.0 + cb2) / sigma;
    }

    Real chi(const Real nu_tilde, const Real molecular_nu)
    {
        return nu_tilde / safe_positive(molecular_nu, "molecular_nu");
    }

    Real fv1(const Real chi_value, const SpalartAllmarasConstants& c)
    {
        const Real chi3 = chi_value * chi_value * chi_value;
        const Real cv13 = c.cv1 * c.cv1 * c.cv1;
        const Real denominator = chi3 + cv13;

        if (std::abs(denominator) <= tiny)
        {
            return 0.0;
        }

        return chi3 / denominator;
    }

    Real fv2(const Real chi_value, const Real fv1_value)
    {
        const Real denominator = 1.0 + chi_value * fv1_value;

        if (std::abs(denominator) <= tiny)
        {
            return 0.0;
        }

        return 1.0 - chi_value / denominator;
    }

    Real ft2(const Real chi_value, const SpalartAllmarasConstants& c)
    {
        return c.ct3 * std::exp(-c.ct4 * chi_value * chi_value);
    }

    Real sBar(
        const Real nu_tilde,
        const Real molecular_nu,
        const Real wall_distance,
        const SpalartAllmarasConstants& c)
    {
        const Real d = safe_positive(wall_distance, "wall_distance");
        const Real chi_value = chi(nu_tilde, molecular_nu);
        const Real fv1_value = fv1(chi_value, c);
        const Real fv2_value = fv2(chi_value, fv1_value);

        return nu_tilde * fv2_value / (c.kappa * c.kappa * d * d);
    }

    Real limitedSTilde(
        const Real omega,
        const Real s_bar,
        const SpalartAllmarasConstants& c)
    {
        if (!std::isfinite(omega) || omega < 0.0)
        {
            throw std::runtime_error("SpalartAllmaras - omega must be non-negative and finite");
        }

        if (!std::isfinite(s_bar))
        {
            throw std::runtime_error("SpalartAllmaras - s_bar must be finite");
        }

        if (s_bar >= -c.c2 * omega)
        {
            return std::max(omega + s_bar, 0.0);
        }

        const Real denominator = (c.c3 - 2.0 * c.c2) * omega - s_bar;

        if (std::abs(denominator) <= tiny)
        {
            return std::max(omega, 0.0);
        }

        const Real limited = omega
            + omega * (c.c2 * c.c2 * omega + c.c3 * s_bar) / denominator;

        return std::max(limited, 0.0);
    }

    Real rFunction(
        const Real nu_tilde,
        const Real s_tilde,
        const Real wall_distance,
        const SpalartAllmarasConstants& c)
    {
        const Real d = safe_positive(wall_distance, "wall_distance");

        if (!std::isfinite(s_tilde) || s_tilde <= tiny)
        {
            return 10.0;
        }

        const Real denominator = s_tilde * c.kappa * c.kappa * d * d;

        if (denominator <= tiny)
        {
            return 10.0;
        }

        return std::clamp(nu_tilde / denominator, 0.0, 10.0);
    }

    Real fw(const Real r_value, const SpalartAllmarasConstants& c)
    {
        const Real r = std::clamp(r_value, 0.0, 10.0);
        const Real r6 = std::pow(r, 6.0);
        const Real g = r + c.cw2 * (r6 - r);
        const Real g6 = std::pow(g, 6.0);
        const Real cw36 = std::pow(c.cw3, 6.0);

        return g * std::pow((1.0 + cw36) / (g6 + cw36), 1.0 / 6.0);
    }

    Real eddyKinematicViscosity(
        const Real nu_tilde,
        const Real molecular_nu,
        const SpalartAllmarasConstants& c)
    {
        if (nu_tilde <= 0.0)
        {
            return 0.0;
        }

        const Real chi_value = chi(nu_tilde, molecular_nu);
        return nu_tilde * fv1(chi_value, c);
    }

    Real negativeBranchFn(
        const Real chi_value,
        const SpalartAllmarasConstants& c)
    {
        const Real chi3 = chi_value * chi_value * chi_value;
        const Real denominator = c.cn1 - chi3;

        if (std::abs(denominator) <= tiny)
        {
            throw std::runtime_error("SpalartAllmaras - singular SA-neg fn denominator");
        }

        return (c.cn1 + chi3) / denominator;
    }

    Real productionTerm(
        const Real nu_tilde,
        const Real s_tilde,
        const Real ft2_value,
        const SpalartAllmarasConstants& c)
    {
        return c.cb1 * (1.0 - ft2_value) * s_tilde * nu_tilde;
    }

    Real destructionCoefficient(
        const Real wall_distance,
        const Real fw_value,
        const Real ft2_value,
        const SpalartAllmarasConstants& c)
    {
        const Real d = safe_positive(wall_distance, "wall_distance");
        return (c.cw1() * fw_value - c.cb1 * ft2_value / (c.kappa * c.kappa)) / (d * d);
    }

    Real destructionTerm(
        const Real nu_tilde,
        const Real wall_distance,
        const Real fw_value,
        const Real ft2_value,
        const SpalartAllmarasConstants& c)
    {
        return destructionCoefficient(wall_distance, fw_value, ft2_value, c)
            * nu_tilde * nu_tilde;
    }
}
