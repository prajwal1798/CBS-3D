#pragma once

//=============================================================================
// CBS3D++_SI
//
// Spalart-Allmaras turbulence assembly interface.
//
// This class owns the CBS-level entry points for turbulence quantities.  The
// first milestone does not yet assemble the transported SA equation.  It only
// performs the safe conversion
//
//     nu_tilde  ->  nu_t  ->  mu_t  ->  mu_eff
//
// and, when requested, the turbulent thermal conductivity correction
//
//     k_eff = k + rho cp nu_t / Pr_t
//
// The laminar path must remain untouched.  Therefore, when turbulence_on = 0,
// the effective properties are reset to their molecular values:
//
//     mu_eff_e = mu_e
//     k_eff_e  = k_e
//
// No pressure equation, momentum equation, or temperature equation is advanced
// in this class.  Those solver steps will read mu_eff_e and k_eff_e only after
// the turbulence coupling is deliberately enabled in a later milestone.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"

namespace cbs
{
    class SpalartAllmarasAssembly
    {
    public:
        static void resetEffectiveProperties(CBSStateSI& s);
        static void updateEddyViscosity(CBSStateSI& s);
    };
}
