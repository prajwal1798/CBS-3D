#include "cbs/turbulence/SpalartAllmaras.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace cbs
{
    namespace turbulence
    {
        namespace
        {
            const Real small_number = 1.0e-30;

            Real require_positive_value(Real value, const char* name)
            {
                if (!std::isfinite(value) || value <= 0.0)
                {
                    throw std::runtime_error(
                        std::string("SpalartAllmaras - non-positive or non-finite ")
                        + name);
                }

                return value;
            }

            Real clamp_real(Real value, Real lower, Real upper)
            {
                if (value < lower)
                {
                    return lower;
                }

                if (value > upper)
                {
                    return upper;
                }

                return value;
            }
        }

        //=========================================================================
        // Returns the derived destruction constant
        //
        //     cw1 = cb1/kappa^2 + (1 + cb2)/sigma
        //
        // This value is not an independent model constant.  Computing it here
        // prevents the code from silently using an inconsistent set of constants.
        //=========================================================================
        Real SpalartAllmarasConstants::cw1() const
        {
            return cb1 / (kappa * kappa) + (1.0 + cb2) / sigma;
        }

        //=========================================================================
        // Calculates
        //
        //     chi = nu_tilde / nu
        //
        // where nu is the molecular kinematic viscosity.  The molecular viscosity
        // must be positive; otherwise the SA non-dimensional functions are not
        // mathematically defined.
        //=========================================================================
        Real chi(Real nu_tilde, Real molecular_nu)
        {
            return nu_tilde / require_positive_value(molecular_nu, "molecular_nu");
        }

        //=========================================================================
        // Calculates the viscous damping function
        //
        //     fv1 = chi^3 / (chi^3 + cv1^3)
        //
        // This function converts the transported SA variable into eddy viscosity.
        //=========================================================================
        Real fv1(Real chi_value, const SpalartAllmarasConstants& c)
        {
            if (!(chi_value > 0.0) || !std::isfinite(chi_value))
            {
                return 0.0;
            }

            const Real cv13 = c.cv1 * c.cv1 * c.cv1;

            // For large chi the direct form overflows: chi^3 exceeds the double
            // range once chi passes roughly 5.6e102, and the expression then
            // evaluates as inf/inf = NaN, which propagates silently into nu_t,
            // mu_t and mu_eff.  The algebraically identical reciprocal form
            //
            //     fv1 = 1 / (1 + (cv1/chi)^3)
            //
            // is exact for the same inputs and cannot overflow, because
            // cv1/chi < 1 in this branch.  The crossover is placed well below
            // the overflow limit so that both branches are evaluated in a range
            // where each is accurate.
            //
            // This is a robustness guard, not a stability fix.  A chi large
            // enough to reach either branch already indicates a diverging
            // solution; the value of this change is that the divergence is
            // reported by the SA ceiling check rather than surfacing several
            // steps later as a non-finite viscosity in the momentum assembly.
            if (chi_value > 1.0e6)
            {
                const Real ratio = c.cv1 / chi_value;
                const Real ratio3 = ratio * ratio * ratio;

                return 1.0 / (1.0 + ratio3);
            }

            const Real chi3 = chi_value * chi_value * chi_value;
            const Real denominator = chi3 + cv13;

            if (std::abs(denominator) <= small_number)
            {
                return 0.0;
            }

            return chi3 / denominator;
        }

        //=========================================================================
        // Calculates
        //
        //     fv2 = 1 - chi/(1 + chi fv1)
        //
        // This term appears in the modified vorticity S_tilde.
        //=========================================================================
        Real fv2(Real chi_value, Real fv1_value)
        {
            const Real denominator = 1.0 + chi_value * fv1_value;

            if (std::abs(denominator) <= small_number)
            {
                return 0.0;
            }

            return 1.0 - chi_value / denominator;
        }

        //=========================================================================
        // Calculates the trip-suppression function used by the standard SA
        // equation:
        //
        //     ft2 = ct3 exp(-ct4 chi^2)
        //=========================================================================
        Real ft2(Real chi_value, const SpalartAllmarasConstants& c)
        {
            return c.ct3 * std::exp(-c.ct4 * chi_value * chi_value);
        }

        //=========================================================================
        // Calculates the correction term
        //
        //     S_bar = nu_tilde fv2 / (kappa^2 d^2)
        //
        // where d is the true distance to the nearest no-slip wall.  The wall
        // distance must already have been computed by the wall-distance module.
        //=========================================================================
        Real sBar(
            Real nu_tilde,
            Real molecular_nu,
            Real wall_distance,
            const SpalartAllmarasConstants& c)
        {
            const Real d = require_positive_value(wall_distance, "wall_distance");
            const Real chi_value = chi(nu_tilde, molecular_nu);
            const Real fv1_value = fv1(chi_value, c);
            const Real fv2_value = fv2(chi_value, fv1_value);

            return nu_tilde * fv2_value / (c.kappa * c.kappa * d * d);
        }

        //=========================================================================
        // Calculates the limited modified vorticity S_tilde.
        //
        // Without protection, S_tilde can become zero or negative when S_bar is
        // sufficiently negative.  That would make the r-function undefined or
        // numerically unstable.  The limiter used here follows the usual
        // Allmaras-Johnson-Spalart protected form.
        //=========================================================================
        Real limitedSTilde(
            Real omega,
            Real s_bar,
            const SpalartAllmarasConstants& c)
        {
            if (!std::isfinite(omega) || omega < 0.0)
            {
                throw std::runtime_error(
                    "SpalartAllmaras - omega must be non-negative and finite");
            }

            if (!std::isfinite(s_bar))
            {
                throw std::runtime_error(
                    "SpalartAllmaras - s_bar must be finite");
            }

            if (s_bar >= -c.c2 * omega)
            {
                const Real value = omega + s_bar;
                return value > 0.0 ? value : 0.0;
            }

            const Real denominator = (c.c3 - 2.0 * c.c2) * omega - s_bar;

            if (std::abs(denominator) <= small_number)
            {
                return omega > 0.0 ? omega : 0.0;
            }

            const Real value = omega
                + omega * (c.c2 * c.c2 * omega + c.c3 * s_bar) / denominator;

            return value > 0.0 ? value : 0.0;
        }

        //=========================================================================
        // Calculates the wall-destruction argument
        //
        //     r = nu_tilde / (S_tilde kappa^2 d^2)
        //
        // and limits it to [0, 10].  If S_tilde is zero or extremely small, the
        // limiting value 10 is returned.
        //=========================================================================
        Real rFunction(
            Real nu_tilde,
            Real s_tilde,
            Real wall_distance,
            const SpalartAllmarasConstants& c)
        {
            const Real d = require_positive_value(wall_distance, "wall_distance");

            if (!std::isfinite(s_tilde) || s_tilde <= small_number)
            {
                return 10.0;
            }

            const Real denominator = s_tilde * c.kappa * c.kappa * d * d;

            if (denominator <= small_number)
            {
                return 10.0;
            }

            return clamp_real(nu_tilde / denominator, 0.0, 10.0);
        }

        //=========================================================================
        // Calculates the wall function
        //
        //     fw = g [(1 + cw3^6)/(g^6 + cw3^6)]^(1/6)
        //
        // where
        //
        //     g = r + cw2 (r^6 - r)
        //=========================================================================
        Real fw(Real r_value, const SpalartAllmarasConstants& c)
        {
            const Real r = clamp_real(r_value, 0.0, 10.0);
            const Real r6 = std::pow(r, 6.0);
            const Real g = r + c.cw2 * (r6 - r);
            const Real g6 = std::pow(g, 6.0);
            const Real cw36 = std::pow(c.cw3, 6.0);

            return g * std::pow((1.0 + cw36) / (g6 + cw36), 1.0 / 6.0);
        }

        //=========================================================================
        // Converts the transported SA variable into turbulent kinematic viscosity:
        //
        //     nu_t = nu_tilde fv1
        //
        // The first implementation is the non-negative standard SA path.  If a
        // numerical update gives a non-positive nu_tilde, it is not allowed to
        // produce negative turbulent viscosity.
        //=========================================================================
        Real eddyKinematicViscosity(
            Real nu_tilde,
            Real molecular_nu,
            const SpalartAllmarasConstants& c)
        {
            if (nu_tilde <= 0.0)
            {
                return 0.0;
            }

            const Real chi_value = chi(nu_tilde, molecular_nu);
            return nu_tilde * fv1(chi_value, c);
        }

        //=========================================================================
        // SA-neg auxiliary function:
        //
        //     fn = (cn1 + chi^3)/(cn1 - chi^3)
        //
        // This function is not used in the standard SA branch.  It is kept here so
        // that the later robust negative branch can be added without changing the
        // public turbulence interface.
        //=========================================================================
        Real negativeBranchFn(
            Real chi_value,
            const SpalartAllmarasConstants& c)
        {
            const Real chi3 = chi_value * chi_value * chi_value;
            const Real denominator = c.cn1 - chi3;

            if (std::abs(denominator) <= small_number)
            {
                throw std::runtime_error(
                    "SpalartAllmaras - singular SA-neg fn denominator");
            }

            return (c.cn1 + chi3) / denominator;
        }

        //=========================================================================
        // Calculates the pointwise production term
        //
        //     P = cb1 (1 - ft2) S_tilde nu_tilde
        //=========================================================================
        Real productionTerm(
            Real nu_tilde,
            Real s_tilde,
            Real ft2_value,
            const SpalartAllmarasConstants& c)
        {
            return c.cb1 * (1.0 - ft2_value) * s_tilde * nu_tilde;
        }

        //=========================================================================
        // Calculates the coefficient multiplying nu_tilde^2 in the destruction
        // term:
        //
        //     D = C_D nu_tilde^2
        //
        // where
        //
        //     C_D = [cw1 fw - cb1 ft2/kappa^2]/d^2
        //
        // This coefficient is useful when treating destruction semi-implicitly.
        //=========================================================================
        Real destructionCoefficient(
            Real wall_distance,
            Real fw_value,
            Real ft2_value,
            const SpalartAllmarasConstants& c)
        {
            const Real d = require_positive_value(wall_distance, "wall_distance");
            return (c.cw1() * fw_value - c.cb1 * ft2_value / (c.kappa * c.kappa))
                / (d * d);
        }

        //=========================================================================
        // Calculates the pointwise destruction term
        //
        //     D = [cw1 fw - cb1 ft2/kappa^2] (nu_tilde/d)^2
        //=========================================================================
        Real destructionTerm(
            Real nu_tilde,
            Real wall_distance,
            Real fw_value,
            Real ft2_value,
            const SpalartAllmarasConstants& c)
        {
            return destructionCoefficient(wall_distance, fw_value, ft2_value, c)
                * nu_tilde * nu_tilde;
        }
    }
}
