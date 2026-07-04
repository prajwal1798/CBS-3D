//=============================================================================
// CBS3D++_SI
//
// Persistent PETSc pressure solver for semi-implicit CBS Step 2.
//
// The constrained pressure equation is:
//
//     A p = b
//
// The solver uses:
//
//     KSP type       : Conjugate Gradient
//     preconditioner : HYPRE BoomerAMG when available
//                      PETSc GAMG otherwise
//
// Fluid-connected nodes form the PETSc pressure space. Solid-only CHT nodes
// are excluded, and prescribed-pressure nodes are imposed by symmetric
// Dirichlet elimination.
//
// The current PETSc objects are created with PETSC_COMM_SELF and sequential
// matrix/vector constructors. This is therefore a serial PETSc pressure solve.
//
// The matrix, vectors, KSP object and AMG hierarchy are cached and reused
// between CBS iterations.
//=============================================================================

#include "cbs/linalg/PetscPressureSolver.hpp"

#include <stdexcept>

#ifdef CBS3D_USE_PETSC

#include <petscksp.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

namespace cbs
{
    namespace
    {
        // Returns true when a global mesh node belongs to the fluid
        // pressure space.
        bool is_pressure_active(
            const std::vector<char>& active,
            const Int ip)
        {
            return active[static_cast<std::size_t>(ip)] != 0;
        }

        // Returns true when pressure is prescribed at the node.
        bool is_fixed_pressure_node(
            const std::vector<char>& fixed,
            const Int ip)
        {
            return fixed[static_cast<std::size_t>(ip)] != 0;
        }

        // Converts a non-zero PETSc error code into a C++ runtime error.
        void check_petsc(const PetscErrorCode ierr, const char* context)
        {
            if (ierr != 0)
            {
                throw std::runtime_error(
                    std::string("PetscPressureSolver::solvePressure - PETSc call failed in ")
                    + context
                    + ", ierr="
                    + std::to_string(static_cast<int>(ierr)));
            }
        }

        bool g_petsc_started_by_cbs3d = false;
        bool g_petsc_finalize_registered = false;

        // Finalises PETSc only when this solver performed the initialisation.
        // PETSc is not finalised here when it was initialised externally.
        void finalise_petsc_if_owned()
        {
            PetscBool is_finalised = PETSC_FALSE;
            PetscFinalized(&is_finalised);

            if (is_finalised)
            {
                return;
            }

            PetscBool is_initialised = PETSC_FALSE;
            PetscInitialized(&is_initialised);

            if (is_initialised && g_petsc_started_by_cbs3d)
            {
                PetscFinalize();
            }
        }

        // Initialises PETSc lazily on the first pressure solve and registers
        // a process-exit cleanup routine.
        void ensure_petsc_initialised()
        {
            PetscBool is_initialised = PETSC_FALSE;
            check_petsc(PetscInitialized(&is_initialised), "PetscInitialized");

            if (!is_initialised)
            {
                check_petsc(PetscInitialize(nullptr, nullptr, nullptr, nullptr),
                            "PetscInitialize");

                g_petsc_started_by_cbs3d = true;

                if (!g_petsc_finalize_registered)
                {
                    std::atexit(finalise_petsc_if_owned);
                    g_petsc_finalize_registered = true;
                }
            }
        }

        // Returns the compact storage position of one tetrahedral
        // off-diagonal pressure coefficient.
        Int compact_offdiag_index(
            const CBSStateSI& s,
            const Int ie,
            const Int pair_id)
        {
            return (ie - 1) * s.cfg.gsdim + pair_id;
        }

        // Checks the tetrahedral pressure dimensions and PETSc solver
        // tolerances required by the current implementation.
        void validate_pressure_dimensions(const CBSStateSI& s)
        {
            if (s.cfg.ndim != 3 ||
                s.cfg.nep != 4 ||
                s.cfg.gsdim != 6 ||
                s.cfg.gdim != 13)
            {
                throw std::runtime_error(
                    "PetscPressureSolver::solvePressure - CBS3D pressure solve requires ndim=3, nep=4, gsdim=6, gdim=13");
            }

            if (s.cfg.cg_max_iter < 1)
            {
                throw std::runtime_error(
                    "PetscPressureSolver::solvePressure - cg_max_iter must be positive");
            }

            if (s.cfg.relToler <= 0.0 || s.cfg.absToler <= 0.0)
            {
                throw std::runtime_error(
                    "PetscPressureSolver::solvePressure - relToler and absToler must be positive");
            }
        }

        // Marks every node connected to at least one fluid element:
        //
        //     mat_elem(e) = 0
        //
        // Only these nodes become PETSc pressure degrees of freedom.
        void build_pressure_active_mask(
            const CBSStateSI& s,
            std::vector<char>& active)
        {
            active.assign(static_cast<std::size_t>(s.cfg.npoin + 1), 0);

            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                if (s.mat_elem(ie) != 0)
                {
                    continue;
                }

                for (Int in = 1; in <= s.cfg.nep; ++in)
                {
                    const Int ip = s.intma(in, ie);

                    if (ip < 1 || ip > s.cfg.npoin)
                    {
                        throw std::runtime_error(
                            "PetscPressureSolver::solvePressure - element connectivity node out of range");
                    }

                    active[static_cast<std::size_t>(ip)] = 1;
                }
            }

            Int active_count = 0;
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                if (is_pressure_active(active, ip))
                {
                    ++active_count;
                }
            }

            if (active_count == 0)
            {
                throw std::runtime_error(
                    "PetscPressureSolver::solvePressure - no fluid-connected pressure nodes found");
            }
        }

        // Marks the prescribed/reference pressure nodes and stores their
        // imposed pressure values.
        void build_fixed_pressure_mask(
            const CBSStateSI& s,
            const std::vector<char>& active,
            std::vector<char>& fixed,
            std::vector<Real>& fixed_values)
        {
            fixed.assign(static_cast<std::size_t>(s.cfg.npoin + 1), 0);
            fixed_values.assign(static_cast<std::size_t>(s.cfg.npoin + 1), 0.0);

            for (Int i = 1; i <= s.cfg.bc_fixed; ++i)
            {
                const Int ip = s.bc_list(i);

                if (ip < 1 || ip > s.cfg.npoin)
                {
                    throw std::runtime_error(
                        "PetscPressureSolver::solvePressure - fixed pressure node out of range");
                }

                if (!is_pressure_active(active, ip))
                {
                    throw std::runtime_error(
                        "PetscPressureSolver::solvePressure - fixed pressure node is not fluid-connected: node "
                        + std::to_string(ip));
                }

                fixed[static_cast<std::size_t>(ip)] = 1;
                fixed_values[static_cast<std::size_t>(ip)] = s.bc_values(i);
            }
        }

        // Creates the mapping:
        //
        //     global mesh node -> contiguous PETSc pressure DOF
        //
        // Solid-only nodes receive the value -1 and are excluded.
        PetscInt build_pressure_dof_map(
            const CBSStateSI& s,
            const std::vector<char>& active,
            std::vector<PetscInt>& pressure_dof)
        {
            pressure_dof.assign(static_cast<std::size_t>(s.cfg.npoin + 1), -1);

            PetscInt ndof = 0;
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                if (is_pressure_active(active, ip))
                {
                    pressure_dof[static_cast<std::size_t>(ip)] = ndof;
                    ++ndof;
                }
            }

            if (ndof <= 0)
            {
                throw std::runtime_error(
                    "PetscPressureSolver::solvePressure - no PETSc pressure DOFs created");
            }

            return ndof;
        }

        // Adds one scalar contribution to the PETSc right-hand side.
        void set_rhs_value(Vec b, const PetscInt row, const Real value)
        {
            const PetscScalar scalar_value = static_cast<PetscScalar>(value);
            check_petsc(VecSetValue(b, row, scalar_value, ADD_VALUES),
                        "VecSetValue(ADD_VALUES)");
        }

        // Adds one additional scalar contribution to the PETSc RHS.
        void add_rhs_value(Vec b, const PetscInt row, const Real value)
        {
            const PetscScalar scalar_value = static_cast<PetscScalar>(value);
            check_petsc(VecSetValue(b, row, scalar_value, ADD_VALUES),
                        "VecSetValue(ADD_VALUES)");
        }

        // Adds one coefficient to the PETSc sparse pressure matrix.
        void add_matrix_value(
            Mat A,
            const PetscInt row,
            const PetscInt col,
            const Real value)
        {
            const PetscScalar scalar_value = static_cast<PetscScalar>(value);
            check_petsc(MatSetValue(A, row, col, scalar_value, ADD_VALUES),
                        "MatSetValue");
        }

        // Calculates the Euclidean norm:
        //
        //     ||v||_2
        Real vector_l2_norm(Vec v)
        {
            PetscReal norm = 0.0;
            check_petsc(VecNorm(v, NORM_2, &norm), "VecNorm(NORM_2)");
            return static_cast<Real>(norm);
        }

        // Calculates the infinity norm:
        //
        //     ||v||_infinity = max_i |v_i|
        Real vector_max_abs(Vec v)
        {
            PetscReal norm = 0.0;
            check_petsc(VecNorm(v, NORM_INFINITY, &norm), "VecNorm(NORM_INFINITY)");
            return static_cast<Real>(norm);
        }

        // Calculates the true algebraic residual:
        //
        //     r = b - A x
        void compute_residual(Mat A, Vec x, Vec b, Vec residual)
        {
            check_petsc(MatMult(A, x, residual), "MatMult(residual)");
            check_petsc(VecAYPX(residual, -1.0, b), "VecAYPX(residual)");
        }

        // Configures the PETSc linear solver.
        //
        //     KSPCG                Conjugate Gradient iteration
        //     PCHYPRE/boomeramg    AMG preconditioner when HYPRE is available
        //     PCGAMG               PETSc AMG fallback
        //
        // KSPSetFromOptions is called last so command-line PETSc options can
        // override the defaults used here.
        void configure_ksp(const CBSStateSI& s, KSP ksp)
        {
            check_petsc(KSPSetType(ksp, KSPCG), "KSPSetType(KSPCG)");

            PC pc = nullptr;
            check_petsc(KSPGetPC(ksp, &pc), "KSPGetPC");

#if defined(PETSC_HAVE_HYPRE)
            check_petsc(PCSetType(pc, PCHYPRE), "PCSetType(PCHYPRE)");
            check_petsc(PCHYPRESetType(pc, "boomeramg"),
                        "PCHYPRESetType(boomeramg)");
#else
            check_petsc(PCSetType(pc, PCGAMG), "PCSetType(PCGAMG)");
#endif

            check_petsc(KSPSetTolerances(ksp,
                                         static_cast<PetscReal>(s.cfg.relToler),
                                         static_cast<PetscReal>(s.cfg.absToler),
                                         PETSC_DEFAULT,
                                         static_cast<PetscInt>(s.cfg.cg_max_iter)),
                        "KSPSetTolerances");

            check_petsc(KSPSetInitialGuessNonzero(ksp, PETSC_TRUE),
                        "KSPSetInitialGuessNonzero");

            check_petsc(KSPSetFromOptions(ksp), "KSPSetFromOptions");
        }

        struct PersistentPressureSystem
        {
            bool ready = false;
            Int npoin = 0;
            Int nelem = 0;
            Int bc_fixed = 0;
            PetscInt ndof = 0;

            std::vector<char> active;
            std::vector<char> fixed;
            std::vector<Real> fixed_values;
            std::vector<PetscInt> pressure_dof;

            Mat A = nullptr;
            Vec b = nullptr;
            Vec x = nullptr;
            Vec r = nullptr;
            KSP ksp = nullptr;

            // Destroys all cached PETSc objects before PETSc finalisation.
            void destroy()
            {
                PetscBool is_finalised = PETSC_FALSE;
                PetscFinalized(&is_finalised);

                if (is_finalised)
                {
                    A = nullptr;
                    b = nullptr;
                    x = nullptr;
                    r = nullptr;
                    ksp = nullptr;
                    ready = false;
                    return;
                }

                if (ksp != nullptr)
                {
                    check_petsc(KSPDestroy(&ksp), "KSPDestroy(cache)");
                }
                if (r != nullptr)
                {
                    check_petsc(VecDestroy(&r), "VecDestroy(r cache)");
                }
                if (x != nullptr)
                {
                    check_petsc(VecDestroy(&x), "VecDestroy(x cache)");
                }
                if (b != nullptr)
                {
                    check_petsc(VecDestroy(&b), "VecDestroy(b cache)");
                }
                if (A != nullptr)
                {
                    check_petsc(MatDestroy(&A), "MatDestroy(A cache)");
                }

                active.clear();
                fixed.clear();
                fixed_values.clear();
                pressure_dof.clear();

                ready = false;
                npoin = 0;
                nelem = 0;
                bc_fixed = 0;
                ndof = 0;
            }

            // Checks whether the cached pressure-system dimensions match the
            // current case. This check does not compare matrix coefficients.
            bool matches_current_problem(const CBSStateSI& s) const
            {
                return ready &&
                       npoin == s.cfg.npoin &&
                       nelem == s.cfg.nelem &&
                       bc_fixed == s.cfg.bc_fixed;
            }

            // Rebuilds the complete persistent PETSc pressure system:
            //
            //     1. active/fixed pressure masks
            //     2. global-node to PETSc-DOF map
            //     3. sequential matrix and vectors
            //     4. pressure matrix
            //     5. KSP and AMG hierarchy
            void initialise(const CBSStateSI& s)
            {
                destroy();

                npoin = s.cfg.npoin;
                nelem = s.cfg.nelem;
                bc_fixed = s.cfg.bc_fixed;

                build_pressure_active_mask(s, active);
                build_fixed_pressure_mask(s, active, fixed, fixed_values);
                ndof = build_pressure_dof_map(s, active, pressure_dof);

                create_petsc_objects();
                assemble_matrix_once(s);
                configure_solver_once(s);

                ready = true;
            }

            // Creates sequential PETSc objects on PETSC_COMM_SELF.
            //
            // This is a serial pressure solve even when the application is
            // launched under MPI.
            void create_petsc_objects()
            {
                check_petsc(MatCreateSeqAIJ(PETSC_COMM_SELF,
                                            ndof,
                                            ndof,
                                            80,
                                            nullptr,
                                            &A),
                            "MatCreateSeqAIJ(cache)");

                check_petsc(MatSetOption(A,
                                         MAT_NEW_NONZERO_ALLOCATION_ERR,
                                         PETSC_FALSE),
                            "MatSetOption(MAT_NEW_NONZERO_ALLOCATION_ERR)");

                check_petsc(VecCreateSeq(PETSC_COMM_SELF, ndof, &b),
                            "VecCreateSeq(b cache)");
                check_petsc(VecDuplicate(b, &x), "VecDuplicate(x cache)");
                check_petsc(VecDuplicate(b, &r), "VecDuplicate(r cache)");

                check_petsc(KSPCreate(PETSC_COMM_SELF, &ksp),
                            "KSPCreate(cache)");
            }

            // Assembles the constrained PETSc pressure matrix.
            //
            // Free-node diagonal:
            //
            //     A_ii = pdiag(i)
            //
            // Free-free tetrahedral couplings are inserted symmetrically:
            //
            //     A_ij = A_ji = gstif_ij
            //
            // Fixed-pressure rows are inserted as identity rows:
            //
            //     A_ii = 1
            //
            // Couplings involving a fixed node are omitted here and transferred
            // to the free-node RHS during every CBS iteration.
            void assemble_matrix_once(const CBSStateSI& s)
            {
                // The matrix is assembled once and reused by the current cache.
                // This is valid only while pdiag and gstif remain unchanged.
                // The per-iteration pressure RHS and initial guess are updated
                // separately.
                for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
                {
                    if (!is_pressure_active(active, ip))
                    {
                        continue;
                    }

                    const PetscInt row = pressure_dof[static_cast<std::size_t>(ip)];

                    if (is_fixed_pressure_node(fixed, ip))
                    {
                        add_matrix_value(A, row, row, 1.0);
                        continue;
                    }

                    const Real diag = s.pdiag(ip);

                    if (!std::isfinite(diag) || diag <= 1.0e-30)
                    {
                        throw std::runtime_error(
                            "PetscPressureSolver::solvePressure - invalid pressure diagonal at active free node "
                            + std::to_string(ip)
                            + ", pdiag="
                            + std::to_string(diag));
                    }

                    add_matrix_value(A, row, row, diag);
                }

                const Int pair_a[7] = {0, 1, 1, 1, 2, 2, 3};
                const Int pair_b[7] = {0, 2, 3, 4, 3, 4, 4};

                for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
                {
                    if (s.mat_elem(ie) != 0)
                    {
                        continue;
                    }

                    for (Int pair_id = 1; pair_id <= s.cfg.gsdim; ++pair_id)
                    {
                        const Int local_a = pair_a[pair_id];
                        const Int local_b = pair_b[pair_id];

                        const Int ipa = s.intma(local_a, ie);
                        const Int ipb = s.intma(local_b, ie);

                        if (!is_pressure_active(active, ipa) ||
                            !is_pressure_active(active, ipb))
                        {
                            throw std::runtime_error(
                                "PetscPressureSolver::solvePressure - fluid element contains inactive pressure node");
                        }

                        const bool fixed_a = is_fixed_pressure_node(fixed, ipa);
                        const bool fixed_b = is_fixed_pressure_node(fixed, ipb);

                        if (fixed_a || fixed_b)
                        {
                            continue;
                        }

                        const Real kij = s.gstif(compact_offdiag_index(s, ie, pair_id));

                        if (std::abs(kij) <= 0.0)
                        {
                            continue;
                        }

                        if (!std::isfinite(kij))
                        {
                            throw std::runtime_error(
                                "PetscPressureSolver::solvePressure - invalid compact pressure off-diagonal at element "
                                + std::to_string(ie));
                        }

                        const PetscInt row_a = pressure_dof[static_cast<std::size_t>(ipa)];
                        const PetscInt row_b = pressure_dof[static_cast<std::size_t>(ipb)];

                        add_matrix_value(A, row_a, row_b, kij);
                        add_matrix_value(A, row_b, row_a, kij);
                    }
                }

                check_petsc(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY),
                            "MatAssemblyBegin(cache)");
                check_petsc(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY),
                            "MatAssemblyEnd(cache)");
            }

            // Attaches the pressure matrix to KSP and explicitly builds the
            // AMG hierarchy once through KSPSetUp().
            void configure_solver_once(const CBSStateSI& s)
            {
                check_petsc(KSPSetOperators(ksp, A, A),
                            "KSPSetOperators(cache)");
                configure_ksp(s, ksp);

                // This explicitly builds the AMG hierarchy once.  Later calls only
                // update b/x and reuse the same KSP/PC objects.
                check_petsc(KSPSetUp(ksp), "KSPSetUp(cache)");
            }

            // Updates the per-iteration pressure RHS and non-zero initial
            // pressure guess while reusing the cached matrix and AMG hierarchy.
            void update_rhs_and_initial_guess(const CBSStateSI& s)
            {
                check_petsc(VecSet(b, 0.0), "VecSet(b cache)");
                check_petsc(VecSet(x, 0.0), "VecSet(x cache)");

                for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
                {
                    if (!is_pressure_active(active, ip))
                    {
                        continue;
                    }

                    const PetscInt row = pressure_dof[static_cast<std::size_t>(ip)];

                    if (is_fixed_pressure_node(fixed, ip))
                    {
                        const Real fixed_value = fixed_values[static_cast<std::size_t>(ip)];
                        set_rhs_value(b, row, fixed_value);
                        check_petsc(VecSetValue(x,
                                                 row,
                                                 static_cast<PetscScalar>(fixed_value),
                                                 INSERT_VALUES),
                                    "VecSetValue(x fixed cache)");
                        continue;
                    }

                    set_rhs_value(b, row, s.rhs1(ip));
                    check_petsc(VecSetValue(x,
                                             row,
                                             static_cast<PetscScalar>(s.pres(ip)),
                                             INSERT_VALUES),
                                "VecSetValue(x initial cache)");
                }

                add_fixed_pressure_elimination_to_rhs(s);

                check_petsc(VecAssemblyBegin(b), "VecAssemblyBegin(b cache)");
                check_petsc(VecAssemblyEnd(b), "VecAssemblyEnd(b cache)");

                check_petsc(VecAssemblyBegin(x), "VecAssemblyBegin(x cache)");
                check_petsc(VecAssemblyEnd(x), "VecAssemblyEnd(x cache)");
            }

            // Applies symmetric Dirichlet elimination.
            //
            // Partitioning the pressure system into free and prescribed DOFs:
            //
            //     A_ff p_f + A_fd p_d = b_f
            //
            // gives the reduced free-node equation:
            //
            //     A_ff p_f = b_f - A_fd p_d
            //
            // Each free-fixed tetrahedral coupling therefore contributes:
            //
            //     b_f <- b_f - A_fd p_d
            void add_fixed_pressure_elimination_to_rhs(const CBSStateSI& s)
            {
                const Int pair_a[7] = {0, 1, 1, 1, 2, 2, 3};
                const Int pair_b[7] = {0, 2, 3, 4, 3, 4, 4};

                for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
                {
                    if (s.mat_elem(ie) != 0)
                    {
                        continue;
                    }

                    for (Int pair_id = 1; pair_id <= s.cfg.gsdim; ++pair_id)
                    {
                        const Int ipa = s.intma(pair_a[pair_id], ie);
                        const Int ipb = s.intma(pair_b[pair_id], ie);

                        const bool fixed_a = is_fixed_pressure_node(fixed, ipa);
                        const bool fixed_b = is_fixed_pressure_node(fixed, ipb);

                        if (fixed_a == fixed_b)
                        {
                            continue;
                        }

                        const Real kij = s.gstif(compact_offdiag_index(s, ie, pair_id));

                        if (std::abs(kij) <= 0.0)
                        {
                            continue;
                        }

                        if (!std::isfinite(kij))
                        {
                            throw std::runtime_error(
                                "PetscPressureSolver::solvePressure - invalid fixed-node pressure coupling at element "
                                + std::to_string(ie));
                        }

                        if (!fixed_a && fixed_b)
                        {
                            const PetscInt row_a = pressure_dof[static_cast<std::size_t>(ipa)];
                            add_rhs_value(b,
                                          row_a,
                                          -kij * fixed_values[static_cast<std::size_t>(ipb)]);
                        }
                        else if (fixed_a && !fixed_b)
                        {
                            const PetscInt row_b = pressure_dof[static_cast<std::size_t>(ipb)];
                            add_rhs_value(b,
                                          row_b,
                                          -kij * fixed_values[static_cast<std::size_t>(ipa)]);
                        }
                    }
                }
            }

            // Copies the PETSc pressure solution back to the global solver
            // state. Solid-only nodes are reset to zero and prescribed-pressure
            // nodes are restored exactly.
            void copy_solution_to_state(CBSStateSI& s)
            {
                const PetscScalar* values = nullptr;
                check_petsc(VecGetArrayRead(x, &values), "VecGetArrayRead(cache)");

                for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
                {
                    if (!is_pressure_active(active, ip))
                    {
                        s.pres(ip) = 0.0;
                        s.rhs1(ip) = 0.0;
                        continue;
                    }

                    if (is_fixed_pressure_node(fixed, ip))
                    {
                        s.pres(ip) = fixed_values[static_cast<std::size_t>(ip)];
                        continue;
                    }

                    const PetscInt dof = pressure_dof[static_cast<std::size_t>(ip)];
                    s.pres(ip) = static_cast<Real>(PetscRealPart(values[dof]));
                }

                check_petsc(VecRestoreArrayRead(x, &values),
                            "VecRestoreArrayRead(cache)");
            }
        };

        PersistentPressureSystem* g_pressure_system = nullptr;
        bool g_pressure_system_cleanup_registered = false;

        // Destroys the persistent pressure cache and finalises PETSc when
        // PETSc was initialised by this solver.
        void cleanup_pressure_system()
        {
            if (g_pressure_system != nullptr)
            {
                g_pressure_system->destroy();
                delete g_pressure_system;
                g_pressure_system = nullptr;
            }

            // The persistent PETSc cache owns PETSc in the current CBS3D++_SI
            // executable path because PETSc is initialised lazily from this solver.
            // Finalise here after all Mat/Vec/KSP objects have been destroyed so
            // mpirun sees a clean MPI/PETSc shutdown even when the process exits
            // immediately after the solver dashboard is printed.
            finalise_petsc_if_owned();
        }

        // Returns the process-wide persistent PETSc pressure-system cache.
        PersistentPressureSystem& pressure_system()
        {
            if (g_pressure_system == nullptr)
            {
                g_pressure_system = new PersistentPressureSystem();

                if (!g_pressure_system_cleanup_registered)
                {
                    std::atexit(cleanup_pressure_system);
                    g_pressure_system_cleanup_registered = true;
                }
            }

            return *g_pressure_system;
        }
    }

    //=========================================================================
    // Explicitly releases the persistent PETSc pressure-system cache.
    //
    // This destroys the matrix, vectors, KSP object and AMG hierarchy before
    // PETSc is finalised.
    //=========================================================================
    void PetscPressureSolver::shutdown()
    {
        cleanup_pressure_system();
    }

    //=========================================================================
    // Solves the constrained fluid-pressure system with PETSc.
    //
    // The algebraic problem is:
    //
    //     A p = b
    //
    // The initial pressure stored in s.pres is used as a non-zero initial guess.
    //
    // Initial residual:
    //
    //     r_0 = b - A p_0
    //
    // PETSc then performs preconditioned Conjugate Gradient iterations using
    // the cached AMG preconditioner.
    //
    // Final convergence measures:
    //
    //     relative L2 = ||r_k||_2 / ||r_0||_2
    //
    //     maximum residual = ||r_k||_infinity
    //
    // Inputs:
    //     pdiag, gstif      pressure matrix coefficients
    //     rhs1              pressure right-hand side
    //     pres              initial pressure guess
    //     bc_list/values    prescribed pressure constraints
    //
    // Output:
    //     pres              converged pressure field
    //=========================================================================
    ConjugateGradient::Result PetscPressureSolver::solvePressure(CBSStateSI& s)
    {
        validate_pressure_dimensions(s);
        ensure_petsc_initialised();

        PersistentPressureSystem& system = pressure_system();

        if (!system.matches_current_problem(s))
        {
            system.initialise(s);
        }

        ConjugateGradient::Result result;

        system.update_rhs_and_initial_guess(s);
        compute_residual(system.A, system.x, system.b, system.r);

        result.initial_l2 = vector_l2_norm(system.r);
        result.final_l2 = result.initial_l2;
        result.final_relative_l2 = (result.initial_l2 > 1.0e-300) ? 1.0 : 0.0;
        result.final_max_abs = vector_max_abs(system.r);

        if (result.initial_l2 <= 1.0e-300)
        {
            result.converged = true;
            system.copy_solution_to_state(s);
            return result;
        }

        check_petsc(KSPSolve(system.ksp, system.b, system.x), "KSPSolve(cache)");

        PetscInt iterations = 0;
        KSPConvergedReason reason = KSP_CONVERGED_ITERATING;

        check_petsc(KSPGetIterationNumber(system.ksp, &iterations),
                    "KSPGetIterationNumber(cache)");
        check_petsc(KSPGetConvergedReason(system.ksp, &reason),
                    "KSPGetConvergedReason(cache)");

        compute_residual(system.A, system.x, system.b, system.r);

        result.iterations = static_cast<Int>(iterations);
        result.final_l2 = vector_l2_norm(system.r);
        result.final_relative_l2 =
            result.initial_l2 > 1.0e-300 ? result.final_l2 / result.initial_l2 : 0.0;
        result.final_max_abs = vector_max_abs(system.r);
        result.converged = (reason > 0);

        system.copy_solution_to_state(s);

        return result;
    }
}

#else

namespace cbs
{
    ConjugateGradient::Result PetscPressureSolver::solvePressure(CBSStateSI&)
    {
        throw std::runtime_error(
            "PetscPressureSolver::solvePressure - CBS3D_USE_PETSC is not enabled in this build");
    }
}

#endif
