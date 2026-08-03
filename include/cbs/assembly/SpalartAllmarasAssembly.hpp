#pragma once

//=============================================================================
// CBS3D++_SI
//
// Spalart-Allmaras turbulence assembly interface.
//
// This class owns the CBS-level routines for the transported SA working
// variable.  The variable advanced by the model is
//
//     nu_tilde
//
// not the eddy viscosity itself.  The eddy viscosity is calculated afterwards
// from
//
//     nu_t = nu_tilde fv1
//     mu_t = rho nu_t
//
// The class therefore contains two separate groups of operations:
//
//     1. SA transport equation assembly and nodal update.
//     2. Conversion of nu_tilde into nu_t, mu_t, mu_eff and k_eff.
//
// The laminar path must remain untouched.  When turbulence_on = 0, the
// effective properties are reset to their molecular values:
//
//     mu_eff_e = mu_e
//     k_eff_e  = k_e
//
// No pressure equation is assembled here.  The pressure correction remains the
// original CBS Step 2 operator.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"

namespace cbs
{
    class SpalartAllmarasAssembly
    {
    public:
        // Restores the turbulence-derived element properties to laminar values.
        static void resetEffectiveProperties(CBSStateSI& s);

        // Converts the current nu_tilde field into nu_t, mu_t, mu_eff and k_eff.
        static void updateEddyViscosity(CBSStateSI& s);

        // Assembles the finite-element residual of the SA transport equation.
        static void assembleTransportRhs(CBSStateSI& s);

        // Applies the nodal SA update using the assembled residual.
        static void updateNuTilde(CBSStateSI& s);
    };
}
