#pragma once

//=============================================================================
// CBS3D++_SI
//
// Continuous thermal law-of-the-wall algebra used by the CHT wall treatment.
//
// The implementation uses the Kader temperature wall function
//
//   T+ = Pr y+ exp(-Gamma)
//      + [2.12 ln(1+y+) + beta] exp(-1/Gamma)
//
//   beta  = (3.85 Pr^(1/3) - 1.3)^2 + 2.12 ln(Pr)
//   Gamma = 0.01 (Pr y+)^4 / [1 + 5 Pr^3 y+]
//
// and converts it to a wall-normal effective conductivity for a first fluid
// cell of height y:
//
//   q_w = rho cp u_tau (T_w - T_p) / T+
//   k_n = rho cp u_tau y / T+.
//
// In the viscous/conductive limit y+ -> 0, T+ -> Pr y+ and therefore
// k_n -> k_molecular exactly.  No empirical clipping is used.
//=============================================================================

#include "cbs/core/Types.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cbs
{
    namespace turbulence
    {
        struct ThermalWallTreatmentResult
        {
            Real molecular_prandtl = 0.0;
            Real y_plus = 0.0;
            Real t_plus = 0.0;
            Real wall_normal_conductivity = 0.0;
        };

        class ThermalWallTreatment
        {
        private:
            static Real checked_prandtl(
                const Real density,
                const Real specific_heat,
                const Real dynamic_viscosity,
                const Real molecular_conductivity)
            {
                if (!(density > 0.0) || !std::isfinite(density) ||
                    !(specific_heat > 0.0) || !std::isfinite(specific_heat) ||
                    !(dynamic_viscosity > 0.0) || !std::isfinite(dynamic_viscosity) ||
                    !(molecular_conductivity > 0.0) || !std::isfinite(molecular_conductivity))
                {
                    throw std::runtime_error(
                        "ThermalWallTreatment - molecular material properties must be positive and finite");
                }

                const Real pr =
                    dynamic_viscosity * specific_heat / molecular_conductivity;

                if (!(pr > 0.0) || !std::isfinite(pr))
                {
                    throw std::runtime_error(
                        "ThermalWallTreatment - invalid molecular Prandtl number");
                }

                return pr;
            }

        public:
            // Kader (1981) continuous temperature law of the wall.
            static Real kaderTPlus(
                const Real y_plus,
                const Real molecular_prandtl)
            {
                if (y_plus < 0.0 || !std::isfinite(y_plus))
                {
                    throw std::runtime_error(
                        "ThermalWallTreatment::kaderTPlus - y_plus must be non-negative and finite");
                }

                if (!(molecular_prandtl > 0.0) ||
                    !std::isfinite(molecular_prandtl))
                {
                    throw std::runtime_error(
                        "ThermalWallTreatment::kaderTPlus - Pr must be positive and finite");
                }

                if (y_plus == 0.0)
                {
                    return 0.0;
                }

                const Real pr13 = std::cbrt(molecular_prandtl);
                const Real beta_core = 3.85 * pr13 - 1.3;
                const Real beta =
                    beta_core * beta_core +
                    2.12 * std::log(molecular_prandtl);

                // Evaluate Gamma through logarithms.  This avoids an artificial
                // overflow in (Pr*y+)^4 on pathological input while retaining
                // the exact Kader expression over the practical wall-function
                // range.
                const Real log_pr = std::log(molecular_prandtl);
                const Real log_y = std::log(y_plus);

                const Real pr3_y =
                    std::exp(std::min(
                        std::log(std::numeric_limits<Real>::max()),
                        3.0 * log_pr + log_y));

                const Real log_denominator =
                    std::log1p(5.0 * pr3_y);

                const Real log_gamma =
                    std::log(0.01) +
                    4.0 * (log_pr + log_y) -
                    log_denominator;

                Real gamma = 0.0;

                if (log_gamma >= std::log(std::numeric_limits<Real>::max()))
                {
                    gamma = std::numeric_limits<Real>::infinity();
                }
                else if (log_gamma > std::log(std::numeric_limits<Real>::min()))
                {
                    gamma = std::exp(log_gamma);
                }

                const Real conductive_weight =
                    std::isinf(gamma) ? 0.0 : std::exp(-gamma);

                const Real logarithmic_weight =
                    gamma > 0.0
                        ? std::exp(-1.0 / gamma)
                        : 0.0;

                const Real t_plus =
                    molecular_prandtl * y_plus * conductive_weight +
                    (2.12 * std::log1p(y_plus) + beta) * logarithmic_weight;

                if (!(t_plus > 0.0) || !std::isfinite(t_plus))
                {
                    throw std::runtime_error(
                        "ThermalWallTreatment::kaderTPlus - invalid T+ evaluation");
                }

                return t_plus;
            }

            static ThermalWallTreatmentResult evaluateKader(
                const Real friction_velocity,
                const Real wall_distance,
                const Real density,
                const Real specific_heat,
                const Real dynamic_viscosity,
                const Real molecular_conductivity)
            {
                if (friction_velocity < 0.0 ||
                    !std::isfinite(friction_velocity))
                {
                    throw std::runtime_error(
                        "ThermalWallTreatment::evaluateKader - u_tau must be non-negative and finite");
                }

                if (!(wall_distance > 0.0) || !std::isfinite(wall_distance))
                {
                    throw std::runtime_error(
                        "ThermalWallTreatment::evaluateKader - wall distance must be positive and finite");
                }

                ThermalWallTreatmentResult result;
                result.molecular_prandtl =
                    checked_prandtl(
                        density,
                        specific_heat,
                        dynamic_viscosity,
                        molecular_conductivity);

                // Exact zero-shear limit.  It is important to retain molecular
                // conduction rather than manufacture a singular thermal wall
                // resistance at stagnation or during startup.
                if (friction_velocity <=
                    16.0 * std::numeric_limits<Real>::epsilon())
                {
                    result.wall_normal_conductivity =
                        molecular_conductivity;
                    return result;
                }

                const Real molecular_nu = dynamic_viscosity / density;
                result.y_plus =
                    friction_velocity * wall_distance / molecular_nu;

                if (!(result.y_plus > 0.0) || !std::isfinite(result.y_plus))
                {
                    throw std::runtime_error(
                        "ThermalWallTreatment::evaluateKader - invalid y+");
                }

                result.t_plus =
                    kaderTPlus(
                        result.y_plus,
                        result.molecular_prandtl);

                result.wall_normal_conductivity =
                    density * specific_heat *
                    friction_velocity * wall_distance /
                    result.t_plus;

                if (!(result.wall_normal_conductivity > 0.0) ||
                    !std::isfinite(result.wall_normal_conductivity))
                {
                    throw std::runtime_error(
                        "ThermalWallTreatment::evaluateKader - invalid wall-normal conductivity");
                }

                return result;
            }
        };
    }
}
