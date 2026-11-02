//=============================================================================
// CBS3D++_SI
//
// Monolithic CBS flow solver for both:
//   1. steady pseudo-time fixed points, and
//   2. one physical backward-Euler timestep.
//
// Per fluid node the unknowns are [u,v,w,u*,v*,w*,p].  The same production
// P1/TET4 Step-1, Step-2 and Step-3 operators are retained.  Picard/Oseen
// iteration freezes only the advecting velocity in the nonlinear Step-1 flux.
//
// Physical BE mode solves:
//
//   A_u(w^k) u^(n+1) + M/dt u*^(n+1) = M/dt u^n
//   dt H p^(n+1) - B u*^(n+1)         = boundary_flux(u^n)
//   M/dt (u^(n+1)-u*^(n+1)) + Gp      = 0
//
// so the first and third rows combine to the backward-Euler momentum balance
//
//   M/dt (u^(n+1)-u^n) + A_u(w^k)u^(n+1) + Gp^(n+1) = 0.
//
// The pressure Laplacian dt H is the existing CBS equal-order stabilisation.
// PETSc FGMRES + PCFIELDSPLIT/Schur is the default linear strategy.
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
            Real maximum_velocity_correction = 0.0;
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
