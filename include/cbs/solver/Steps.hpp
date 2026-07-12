#pragma once

//=============================================================================
// CBS3D++_SI
//
// High-level execution of the semi-implicit CBS steps and the optional
// Spalart-Allmaras turbulence transport step.
//
// The numerical sequence is:
//
//     Step 1     Momentum predictor
//     Step 2     Pressure solve
//     Step 3     Pressure-gradient velocity correction
//     Step SA    Spalart-Allmaras turbulence transport, when enabled
//     Step 4     Energy equation
//
// Step SA is deliberately not called Step 5.  It is not a fifth CBS splitting
// step.  It is an additional transported scalar equation inserted after the
// corrected velocity is available from Step 3 and before the energy equation uses
// the updated turbulent thermal diffusivity.
//
// For the semi-implicit CBS formulation, the standard nodal updates have the
// general form:
//
//     u*       = u^n + D_u^(-1) r_m
//
//     A_p p    = b_p
//
//     u^(n+1)  = u* + D_u^(-1) r_p
//
//     nu_tilde^(n+1) = nu_tilde^n + D_u^(-1) r_SA
//
//     T^(n+1)  = T^n + D_T^(-1) r_T
//
// where:
//
//     D_u^(-1) = elcoe2   inverse fluid lumped mass/time diagonal
//     D_T^(-1) = elcoe2p  inverse thermal-capacitance/time diagonal
//
// The detailed finite-element residuals are assembled in MomentumAssembly,
// PressureAssembly, SpalartAllmarasAssembly and EnergyAssembly.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"

namespace cbs
{
    class Steps
    {
    public:
        // Executes CBS Step 1: momentum predictor.
        static void step1(CBSStateSI& s);

        // Executes CBS Step 2: pressure-system assembly and solution.
        static void step2(CBSStateSI& s);

        // Executes CBS Step 3: velocity correction from the pressure gradient.
        static void step3(CBSStateSI& s);

        // Executes the optional Spalart-Allmaras turbulence transport step.
        static void stepSpalartAllmaras(CBSStateSI& s);

        // Executes CBS Step 4: thermal or CHT energy update.
        static void step4(CBSStateSI& s);

    private:
        static void step1SemiImplicit(CBSStateSI& s);
        static void step2SemiImplicit(CBSStateSI& s);
        static void step3SemiImplicit(CBSStateSI& s);
        static void step4Energy(CBSStateSI& s);

        // Reports that the explicit CBS kernels are not yet available.
        static void rejectExplicitMode(const char* step_name);
    };
}
