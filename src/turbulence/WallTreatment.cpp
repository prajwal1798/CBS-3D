//=============================================================================
// CBS3D++_SI
//
// Continuous Spalding wall-law kernel.
//=============================================================================

#include "cbs/turbulence/WallTreatment.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cbs
{
    namespace turbulence
    {
        namespace
        {
            Real vector_magnitude(const std::array<Real, 3>& value)
            {
                return std::sqrt(
                    value[0] * value[0] +
                    value[1] * value[1] +
                    value[2] * value[2]);
            }

            std::array<Real, 3> normalized_wall_normal(
                const std::array<Real, 3>& normal)
            {
                const Real magnitude = vector_magnitude(normal);

                if (!(magnitude > 0.0) || !std::isfinite(magnitude))
                {
                    throw std::runtime_error(
                        "WallTreatment - wall normal must have finite non-zero magnitude");
                }

                return
                {
                    normal[0] / magnitude,
                    normal[1] / magnitude,
                    normal[2] / magnitude
                };
            }

            void validate_options(const WallTreatmentOptions& options)
            {
                if (!(options.kappa > 0.0) || !std::isfinite(options.kappa))
                {
                    throw std::runtime_error(
                        "WallTreatment - kappa must be positive and finite");
                }

                if (!std::isfinite(options.log_intercept))
                {
                    throw std::runtime_error(
                        "WallTreatment - log_intercept must be finite");
                }

                const Real exponent =
                    -options.kappa * options.log_intercept;

                if (exponent >= std::log(std::numeric_limits<Real>::max()))
                {
                    throw std::runtime_error(
                        "WallTreatment - exp(-kappa*B) would overflow");
                }

                if (!(options.root_relative_tolerance > 0.0) ||
                    !std::isfinite(options.root_relative_tolerance))
                {
                    throw std::runtime_error(
                        "WallTreatment - root tolerance must be positive and finite");
                }

                if (options.root_max_iterations <= 0)
                {
                    throw std::runtime_error(
                        "WallTreatment - root_max_iterations must be positive");
                }

                if (options.tangential_speed_floor < 0.0 ||
                    !std::isfinite(options.tangential_speed_floor))
                {
                    throw std::runtime_error(
                        "WallTreatment - tangential_speed_floor must be non-negative and finite");
                }
            }

            // Returns exp(x) - sum_{n=0}^order x^n/n! without subtractive
            // cancellation when x is small. Only orders 3 and 4 are required.
            Real exponential_remainder(Real x, Int order)
            {
                if (x < 0.0 || !std::isfinite(x))
                {
                    throw std::runtime_error(
                        "WallTreatment - invalid exponential remainder argument");
                }

                if (order != 3 && order != 4)
                {
                    throw std::runtime_error(
                        "WallTreatment - unsupported exponential remainder order");
                }

                if (x < 1.0)
                {
                    const Int first_power = order + 1;
                    Real term = 1.0;

                    for (Int n = 1; n <= first_power; ++n)
                    {
                        term *= x / static_cast<Real>(n);
                    }

                    Real sum = term;

                    for (Int n = first_power + 1; n <= 100; ++n)
                    {
                        term *= x / static_cast<Real>(n);
                        sum += term;

                        if (std::fabs(term) <=
                            1.0e-16 * std::max(1.0, std::fabs(sum)))
                        {
                            break;
                        }
                    }

                    return sum;
                }

                const Real maximum_log =
                    std::log(std::numeric_limits<Real>::max());

                if (x >= maximum_log)
                {
                    return std::numeric_limits<Real>::infinity();
                }

                const Real x2 = x * x;
                const Real x3 = x2 * x;
                Real polynomial =
                    1.0 + x + 0.5 * x2 + x3 / 6.0;

                if (order == 4)
                {
                    polynomial += x3 * x / 24.0;
                }

                return std::exp(x) - polynomial;
            }

            void spalding_value_and_derivative(
                Real u_plus,
                const WallTreatmentOptions& options,
                Real& y_plus,
                Real& dyplus_duplus)
            {
                const Real x = options.kappa * u_plus;
                const Real coefficient =
                    std::exp(-options.kappa * options.log_intercept);

                const Real remainder4 = exponential_remainder(x, 4);
                const Real remainder3 = exponential_remainder(x, 3);

                y_plus = u_plus + coefficient * remainder4;
                dyplus_duplus =
                    1.0 + coefficient * options.kappa * remainder3;
            }
        }

        Real WallTreatment::spaldingYPlus(
            Real u_plus,
            const WallTreatmentOptions& options)
        {
            validate_options(options);

            if (u_plus < 0.0 || !std::isfinite(u_plus))
            {
                throw std::runtime_error(
                    "WallTreatment::spaldingYPlus - u_plus must be non-negative and finite");
            }

            Real y_plus = 0.0;
            Real derivative = 0.0;

            spalding_value_and_derivative(
                u_plus,
                options,
                y_plus,
                derivative);

            if (!(y_plus >= 0.0) || std::isnan(y_plus) ||
                !(derivative > 0.0) || std::isnan(derivative))
            {
                throw std::runtime_error(
                    "WallTreatment::spaldingYPlus - invalid Spalding evaluation");
            }

            return y_plus;
        }

        std::array<Real, 3> WallTreatment::tangentialVelocity(
            const std::array<Real, 3>& velocity,
            const std::array<Real, 3>& wall_normal)
        {
            for (const Real component : velocity)
            {
                if (!std::isfinite(component))
                {
                    throw std::runtime_error(
                        "WallTreatment::tangentialVelocity - velocity must be finite");
                }
            }

            const std::array<Real, 3> normal =
                normalized_wall_normal(wall_normal);

            const Real normal_velocity =
                velocity[0] * normal[0] +
                velocity[1] * normal[1] +
                velocity[2] * normal[2];

            return
            {
                velocity[0] - normal_velocity * normal[0],
                velocity[1] - normal_velocity * normal[1],
                velocity[2] - normal_velocity * normal[2]
            };
        }

        WallTreatmentResult WallTreatment::evaluateSpalding(
            const std::array<Real, 3>& sample_velocity,
            const std::array<Real, 3>& wall_normal,
            Real wall_distance,
            Real density,
            Real molecular_nu,
            const WallTreatmentOptions& options)
        {
            validate_options(options);

            if (!(wall_distance > 0.0) || !std::isfinite(wall_distance))
            {
                throw std::runtime_error(
                    "WallTreatment::evaluateSpalding - wall_distance must be positive and finite");
            }

            if (!(density > 0.0) || !std::isfinite(density))
            {
                throw std::runtime_error(
                    "WallTreatment::evaluateSpalding - density must be positive and finite");
            }

            if (!(molecular_nu > 0.0) || !std::isfinite(molecular_nu))
            {
                throw std::runtime_error(
                    "WallTreatment::evaluateSpalding - molecular_nu must be positive and finite");
            }

            WallTreatmentResult result;
            result.tangential_velocity =
                tangentialVelocity(sample_velocity, wall_normal);

            result.tangential_speed =
                vector_magnitude(result.tangential_velocity);

            if (!(result.tangential_speed > options.tangential_speed_floor))
            {
                return result;
            }

            for (Int i = 0; i < 3; ++i)
            {
                result.unit_tangent[static_cast<Size>(i)] =
                    result.tangential_velocity[static_cast<Size>(i)] /
                    result.tangential_speed;
            }

            // Product identity:
            //
            //     Re_y = U_t y / nu = u+ y+
            //
            // The quartic Spalding law satisfies y+(u+) >= u+ for u+ >= 0,
            // hence F(u+) = u+ y+(u+) - Re_y is strictly monotone and the root
            // is guaranteed to lie in [0, sqrt(Re_y)].
            const Real wall_reynolds =
                result.tangential_speed * wall_distance / molecular_nu;

            if (!(wall_reynolds > 0.0) || !std::isfinite(wall_reynolds))
            {
                throw std::runtime_error(
                    "WallTreatment::evaluateSpalding - invalid wall Reynolds number");
            }

            Real lower = 0.0;
            Real upper = std::sqrt(wall_reynolds);

            const Real logarithmic_guess =
                std::log(std::max(1.0, wall_reynolds)) / options.kappa +
                options.log_intercept;

            Real u_plus = logarithmic_guess;

            if (!(u_plus > lower && u_plus < upper) ||
                !std::isfinite(u_plus))
            {
                u_plus = 0.5 * (lower + upper);
            }

            bool converged = false;

            for (Int iteration = 1;
                 iteration <= options.root_max_iterations;
                 ++iteration)
            {
                Real y_plus = 0.0;
                Real dyplus_duplus = 0.0;

                spalding_value_and_derivative(
                    u_plus,
                    options,
                    y_plus,
                    dyplus_duplus);

                const Real residual =
                    u_plus * y_plus - wall_reynolds;

                result.root_iterations = iteration;

                if (std::fabs(residual) <=
                    options.root_relative_tolerance * wall_reynolds)
                {
                    converged = true;
                    break;
                }

                if (residual < 0.0)
                {
                    lower = u_plus;
                }
                else
                {
                    upper = u_plus;
                }

                const Real bracket_width = upper - lower;
                const Real root_scale = std::max(
                    u_plus,
                    std::sqrt(std::numeric_limits<Real>::epsilon()));

                if (bracket_width <=
                    options.root_relative_tolerance * root_scale)
                {
                    converged = true;
                    u_plus = 0.5 * (lower + upper);
                    break;
                }

                const Real derivative =
                    y_plus + u_plus * dyplus_duplus;

                Real candidate =
                    u_plus - residual / derivative;

                if (!std::isfinite(candidate) ||
                    !(candidate > lower && candidate < upper))
                {
                    candidate = 0.5 * (lower + upper);
                }

                u_plus = candidate;
            }

            if (!converged)
            {
                throw std::runtime_error(
                    "WallTreatment::evaluateSpalding - safeguarded root solve did not converge");
            }

            result.u_plus = u_plus;
            result.y_plus = spaldingYPlus(result.u_plus, options);
            result.friction_velocity =
                result.tangential_speed / result.u_plus;

            result.wall_shear_magnitude =
                density *
                result.friction_velocity *
                result.friction_velocity;

            if (!(result.wall_shear_magnitude > 0.0) ||
                !std::isfinite(result.wall_shear_magnitude))
            {
                throw std::runtime_error(
                    "WallTreatment::evaluateSpalding - invalid wall shear magnitude");
            }

            for (Int i = 0; i < 3; ++i)
            {
                const Size index = static_cast<Size>(i);

                result.wall_shear_on_fluid[index] =
                    -result.wall_shear_magnitude * result.unit_tangent[index];

                result.kinematic_wall_traction[index] =
                    result.wall_shear_on_fluid[index] / density;
            }

            return result;
        }
    }
}
