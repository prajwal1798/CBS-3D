#pragma once

//=============================================================================
// CBS3D++_SI
//
// High-level execution of the four semi-implicit
// Characteristic-Based Split steps.
//
// The numerical sequence is:
//
//     Step 1  Momentum predictor
//     Step 2  Pressure solve
//     Step 3  Pressure-gradient velocity correction
//     Step 4  Energy equation
//
// For the semi-implicit CBS formulation, the nodal updates have the general
// form:
//
//     u*       = u^n + D_u^(-1) r_m
//
//     A_p p    = b_p
//
//     u^(n+1)  = u* + D_u^(-1) r_p
//
//     T^(n+1)  = T^n + D_T^(-1) r_T
//
// where:
//
//     D_u^(-1) = elcoe2   inverse momentum mass/time diagonal
//     D_T^(-1) = elcoe2p  inverse thermal-capacitance/time diagonal
//
// The detailed finite-element residuals are assembled in MomentumAssembly,
// PressureAssembly and EnergyAssembly.
//
// cbs_scheme:
//
//     1  semi-implicit CBS
//     0  explicit CBS, not yet implemented in the current solver
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
