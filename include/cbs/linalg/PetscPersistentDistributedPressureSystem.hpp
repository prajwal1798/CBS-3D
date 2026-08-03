#pragma once

//=============================================================================
// CBS3D++_SI
//
// Persistent distributed PETSc pressure system used by the MPI CBS loop.
//
// The geometry-only pressure matrix is assembled once on MPI_COMM_WORLD:
//
//     K_ij = integral grad(N_i) . grad(N_j) dV
//
// For the global CBS time step dt, the usual equation
//
//     dt K p = b
//
// is solved in the equivalent persistent form
//
//     K p = b / dt.
//
// This keeps the matrix sparsity, coefficients, KSP object and AMG hierarchy
// unchanged between CBS iterations while preserving the original algebra.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"
#include "cbs/linalg/ConjugateGradient.hpp"

namespace cbs
{
    class PetscPersistentDistributedPressureSystem
    {
    public:
        PetscPersistentDistributedPressureSystem();
        ~PetscPersistentDistributedPressureSystem();

        PetscPersistentDistributedPressureSystem(
            const PetscPersistentDistributedPressureSystem&) = delete;

        PetscPersistentDistributedPressureSystem& operator=(
            const PetscPersistentDistributedPressureSystem&) = delete;

        // Builds the distributed pressure DOF map, matrix, vectors, KSP and AMG
        // hierarchy. This routine must be called collectively exactly once.
        void initialise(CBSStateSI& s);

        // Updates only the pressure RHS and initial guess, then solves using the
        // already-built matrix and preconditioner hierarchy.
        ConjugateGradient::Result solve(CBSStateSI& s);

        // Collectively releases all PETSc resources owned by this object.
        void shutdown() noexcept;

        [[nodiscard]] bool ready() const noexcept;

    private:
        struct Impl;
        Impl* impl_ = nullptr;
    };
}
