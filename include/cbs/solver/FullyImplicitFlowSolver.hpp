#pragma once

//=============================================================================
// CBS3D++_SI
//
// Fully implicit steady incompressible flow solver.
//
// The nonlinear CBS fixed point is linearised by Picard/Oseen iteration.
// Predictor velocity is retained as an auxiliary unknown so the monolithic
// solve reproduces the existing Step-1/Step-2/Step-3 algebra:
//
//     M_dt (u* - u) = r_m(u)
//     tau H p        = B u*
//     M_dt (u - u*) + G p = 0
//
// Per fluid node the unknowns are [u,v,w,u*,v*,w*,p]. The pressure Laplacian
// tau H is the existing CBS equal-order stabilisation; no unstabilised P1/P1
// saddle-point discretisation is introduced.
//
// Strong velocity and pressure conditions are imposed algebraically. PETSc
// FGMRES + PCFIELDSPLIT/Schur is the default linear solution strategy.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"

namespace cbs
{
    class FullyImplicitFlowSolver
    {
    public:
        struct Result
        {
            bool converged = false;
            Int nonlinear_iterations = 0;
            Int total_krylov_iterations = 0;
            Real relative_velocity_change = 0.0;
            Real relative_pressure_change = 0.0;
            Real divergence_rms = 0.0;
            Real maximum_velocity = 0.0;
            Real maximum_velocity_change = 0.0;
            Real linear_relative_residual = 0.0;
            double seconds = 0.0;
        };

        // Steady nonlinear CBS fixed-point solve (pseudo-time only).
        static Result solve(CBSStateSI& s);

        // One physical backward-Euler step. s.unkn1 must contain u^n and
        // s.unkno is used as the Picard initial guess for u^(n+1).
        static Result solveBackwardEulerStep(CBSStateSI& s, Real physical_dt);

    private:
        static Result solveImpl(
            CBSStateSI& s,
            bool backward_euler,
            Real physical_dt);
    };
}
