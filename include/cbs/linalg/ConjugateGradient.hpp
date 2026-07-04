#pragma once

//=============================================================================
// CBS3D++_SI
//
// Serial Conjugate Gradient solver for the CBS pressure equation.
//
// The pressure system is:
//
//     A p = b
//
// where A is the symmetric pressure operator assembled over the fluid domain.
// Prescribed-pressure nodes are imposed as fixed rows, while solid-only nodes
// are excluded from the pressure space.
//
// The solver supports:
//
//     cg_preconditioner = 0   no preconditioner
//     cg_preconditioner = 1   Jacobi preconditioner
//
// For Jacobi preconditioning:
//
//     z_i = r_i / A_ii
//
// The Result structure stores the convergence measures and the time spent in
// the main parts of one pressure CG solve.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"

namespace cbs
{
    class ConjugateGradient
    {
    public:
        struct Result
        {
            // Number of completed Pressure CG iterations.
            Int iterations = 0;

            // Initial and final Euclidean residual norms.
            Real initial_l2 = 0.0;
            Real final_l2 = 0.0;

            // Relative residual:
            //
            //     ||r_k||_2 / ||r_0||_2
            Real final_relative_l2 = 0.0;

            // Maximum absolute residual over free pressure nodes.
            Real final_max_abs = 0.0;

            bool converged = false;

            // Timing information for one pressure solve.
            double total_seconds = 0.0;
            double setup_seconds = 0.0;
            double matvec_seconds = 0.0;
            double dot_seconds = 0.0;
            double preconditioner_seconds = 0.0;
            double vector_update_seconds = 0.0;
            double constraint_seconds = 0.0;

            // Number of major algebraic operations performed.
            Int matvec_calls = 0;
            Int dot_calls = 0;
            Int preconditioner_calls = 0;
        };

        // Solves the constrained fluid-pressure system using Conjugate
        // Gradient iteration.
        static Result solvePressure(CBSStateSI& s);
    };
}
