#pragma once

//=============================================================================
// CBS3D++_SI
//
// Algebraic support routines for the Spalart-Allmaras turbulence model.
//
// This file deliberately contains only small scalar formulae.  It does not know
// anything about the mesh, elements, nodes, MPI partitions, or the CBS time
// integration loop.  Keeping these formulae separate makes the turbulence model
// easier to check before it is coupled to the finite-element assembly.
//
// Coding-style rule for this public solver:
//
//     * use plain C++ syntax;
//     * avoid decorative attributes such as [[nodiscard]];
//     * avoid compact namespace syntax;
//     * keep every formula close to the notation used in the turbulence model.
//
// The full SA transport equation is not assembled here.  These routines are used
// by the turbulence assembly to evaluate fv1, fv2, fw, production, destruction,
// and the eddy-viscosity relation.
//=============================================================================

#include "cbs/core/Types.hpp"

namespace cbs
{
    namespace turbulence
    {
        // Standard Spalart-Allmaras constants.
        //
        // cb1, sigma, cb2, kappa, cw2, cw3, cv1, ct3 and ct4 are the usual
        // constants of the standard one-equation SA model.  The constant cw1 is
        // not stored directly because it is derived from cb1, cb2, sigma and
        // kappa.
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

            // Constants used in the Allmaras-Johnson-Spalart limiting of
            // S_tilde.  The limiter avoids an invalid denominator when the
            // unmodified correction S_bar is too negative.
            Real c2 = 0.7;
            Real c3 = 0.9;

            // SA-neg constant.  The negative branch is not activated by the
            // first SA implementation, but the constant is kept here so the
            // interface remains stable when SA-neg is added later.
            Real cn1 = 16.0;

            Real cw1() const;
        };

        // Ratio between the transported SA working variable and the molecular
        // kinematic viscosity.
        Real chi(Real nu_tilde, Real molecular_nu);

        // Viscosity damping functions.
        Real fv1(Real chi_value, const SpalartAllmarasConstants& c);
        Real fv2(Real chi_value, Real fv1_value);
        Real ft2(Real chi_value, const SpalartAllmarasConstants& c);

        // Computes the correction S_bar used in S_tilde.
        Real sBar(
            Real nu_tilde,
            Real molecular_nu,
            Real wall_distance,
            const SpalartAllmarasConstants& c);

        // Applies the protected S_tilde formula.  The input omega is the
        // magnitude of vorticity used by the SA model.
        Real limitedSTilde(
            Real omega,
            Real s_bar,
            const SpalartAllmarasConstants& c);

        // Computes the wall-destruction argument r, limited to the interval
        // [0, 10] for the standard non-negative SA branch.
        Real rFunction(
            Real nu_tilde,
            Real s_tilde,
            Real wall_distance,
            const SpalartAllmarasConstants& c);

        // Wall-destruction function.
        Real fw(Real r_value, const SpalartAllmarasConstants& c);

        // Eddy-viscosity relation:
        //
        //     nu_t = nu_tilde * fv1
        //
        // For the first non-negative SA implementation, negative nu_tilde does
        // not create negative eddy viscosity.  The returned nu_t is zero when
        // nu_tilde is not positive.
        Real eddyKinematicViscosity(
            Real nu_tilde,
            Real molecular_nu,
            const SpalartAllmarasConstants& c);

        // SA-neg auxiliary function.  It is provided for a later robust branch;
        // the standard SA path does not call it during the first implementation.
        Real negativeBranchFn(
            Real chi_value,
            const SpalartAllmarasConstants& c);

        // Pointwise production and destruction terms of the standard SA model.
        Real productionTerm(
            Real nu_tilde,
            Real s_tilde,
            Real ft2_value,
            const SpalartAllmarasConstants& c);

        Real destructionCoefficient(
            Real wall_distance,
            Real fw_value,
            Real ft2_value,
            const SpalartAllmarasConstants& c);

        Real destructionTerm(
            Real nu_tilde,
            Real wall_distance,
            Real fw_value,
            Real ft2_value,
            const SpalartAllmarasConstants& c);
    }
}
