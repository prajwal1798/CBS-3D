#pragma once

//=============================================================================
// CBS3D++_SI
//
// PETSc pressure solver for semi-implicit CBS Step 2.
//
// The pressure equation is:
//
//     A p = b
//
// where the matrix A is assembled from the existing CBS pressure arrays:
//
//     pdiag(i)   assembled nodal diagonal coefficient
//     gstif      compact tetrahedral off-diagonal coefficients
//
// Only nodes connected to fluid elements are included in the PETSc pressure
// space. Solid-only CHT nodes are excluded.
//
// Prescribed-pressure nodes are imposed by symmetric Dirichlet elimination.
// The free pressure block therefore remains symmetric and suitable for:
//
//     KSPCG + BoomerAMG
//
// when PETSc is built with HYPRE. PETSc GAMG is used as the fallback
// preconditioner when HYPRE is unavailable.
//
// The current implementation creates sequential PETSc objects on:
//
//     PETSC_COMM_SELF
//
// It is therefore a serial PETSc pressure solver, not yet a distributed MPI
// pressure solve.
//
// The PETSc matrix, vectors, KSP object and AMG hierarchy are cached and reused
// between CBS iterations while the pressure-system structure is unchanged.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"
#include "cbs/linalg/ConjugateGradient.hpp"

namespace cbs
{
    class PetscPressureSolver
    {
    public:
        // Solves the constrained fluid-pressure system with PETSc CG and AMG.
        static ConjugateGradient::Result solvePressure(CBSStateSI& s);

        // Destroys the persistent PETSc matrix, vectors, KSP object and
        // preconditioner hierarchy.
        static void shutdown();
    };
}
