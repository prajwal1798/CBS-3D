#pragma once

//=============================================================================
// CBS3D++_SI
//
// Pure algebraic functions and constants for the Spalart-Allmaras turbulence
// model.  These routines do not access mesh storage and are therefore suitable
// for deterministic unit testing before the transport equation is coupled into
// the CBS solver.
//=============================================================================

#include "cbs/core/Types.hpp"

namespace cbs::turbulence
{
    struct SpalartAllmarasConstants
    {
        Real cb1 = 0.1355;
        Real sigma = 2.0 / 3.0;
        Real cb2 = 0.622;
        Real kappa = 0.41;
        Real cw2 = 0.3;
        Real cw3 = 2.0;
        Real cv1 = 7.1;
        Real ct3 = 1.2;
        Real ct4 = 0.5;

        // Allmaras-Johnson-Spalart S_tilde limiter constants.
        Real c2 = 0.7;
        Real c3 = 0.9;

        // SA-neg reserved constant.  The negative branch is not activated in
        // Milestone 2 but the value is kept here so the public API is stable.
        Real cn1 = 16.0;

        [[nodiscard]] Real cw1() const noexcept;
    };

    [[nodiscard]] Real chi(Real nu_tilde, Real molecular_nu);
    [[nodiscard]] Real fv1(Real chi_value, const SpalartAllmarasConstants& c = {});
    [[nodiscard]] Real fv2(Real chi_value, Real fv1_value);
    [[nodiscard]] Real ft2(Real chi_value, const SpalartAllmarasConstants& c = {});

    [[nodiscard]] Real sBar(
        Real nu_tilde,
        Real molecular_nu,
        Real wall_distance,
        const SpalartAllmarasConstants& c = {});

    [[nodiscard]] Real limitedSTilde(
        Real omega,
        Real s_bar,
        const SpalartAllmarasConstants& c = {});

    [[nodiscard]] Real rFunction(
        Real nu_tilde,
        Real s_tilde,
        Real wall_distance,
        const SpalartAllmarasConstants& c = {});

    [[nodiscard]] Real fw(Real r_value, const SpalartAllmarasConstants& c = {});

    [[nodiscard]] Real eddyKinematicViscosity(
        Real nu_tilde,
        Real molecular_nu,
        const SpalartAllmarasConstants& c = {});

    [[nodiscard]] Real negativeBranchFn(
        Real chi_value,
        const SpalartAllmarasConstants& c = {});

    [[nodiscard]] Real productionTerm(
        Real nu_tilde,
        Real s_tilde,
        Real ft2_value,
        const SpalartAllmarasConstants& c = {});

    [[nodiscard]] Real destructionCoefficient(
        Real wall_distance,
        Real fw_value,
        Real ft2_value,
        const SpalartAllmarasConstants& c = {});

    [[nodiscard]] Real destructionTerm(
        Real nu_tilde,
        Real wall_distance,
        Real fw_value,
        Real ft2_value,
        const SpalartAllmarasConstants& c = {});
}
