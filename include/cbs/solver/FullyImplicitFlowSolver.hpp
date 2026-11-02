#pragma once

//=============================================================================
// CBS3D++_SI
//
// Fully implicit steady incompressible flow solver.
//
// The nonlinear Navier-Stokes problem is linearised by Picard/Oseen iteration.
// Each nonlinear iteration solves velocity and kinematic pressure in one
// distributed PETSc block system:
//
//     [ A_Oseen(u^k)   G ] [ u^(k+1) ] = [ 0 ]
//     [   -D         tau H] [ p^(k+1) ]   [ 0 ]
//
// The pressure-Laplacian block tau H is the existing CBS equal-order pressure
// stabilisation; this deliberately retains the solver's P1/P1 spatial
// formulation rather than introducing an unstabilised saddle-point method.
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
            Real linear_relative_residual = 0.0;
            double seconds = 0.0;
        };

        static Result solve(CBSStateSI& s);
    };
}
