#include "cbs/linalg/PetscPressureSolver.hpp"

#ifdef CBS3D_USE_PETSC

#include "cbs/parallel/HaloExchange.hpp"
#include <petscksp.h>

#ifdef CBS3D_USE_MPI
#include <mpi.h>
#endif

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace cbs
{
    namespace
    {
        void petsc_check(PetscErrorCode ierr, const char* where)
        {
            if (ierr != 0)
            {
                throw std::runtime_error(
                    std::string("Distributed PETSc pressure failure in ") + where
                    + ", ierr=" + std::to_string(static_cast<int>(ierr)));
            }
        }

#ifdef CBS3D_USE_MPI
        void mpi_check(int ierr, const char* where)
        {
            if (ierr != MPI_SUCCESS)
            {
                throw std::runtime_error(
                    std::string("Distributed pressure MPI failure in ") + where);
            }
        }
#endif

        bool pressure_active(const CBSStateSI& s, Int ip)
        {
            return
                (s.node_material_mask(ip) &
                 CBSStateSI::node_touches_fluid) != 0;
        }

        Int element_node_index(const CBSStateSI& s, Int ie, Int a)
        {
            return (ie - 1) * s.cfg.nep + a;
        }

        Int offdiag_index(const CBSStateSI& s, Int ie, Int pair)
        {
            return (ie - 1) * s.cfg.gsdim + pair;
        }

        void assemble(Vec v)
        {
            petsc_check(VecAssemblyBegin(v), "VecAssemblyBegin");
            petsc_check(VecAssemblyEnd(v), "VecAssemblyEnd");
        }

        Real norm2(Vec v)
        {
            PetscReal value = 0.0;
            petsc_check(VecNorm(v, NORM_2, &value), "VecNorm(NORM_2)");
            return static_cast<Real>(value);
        }

        Real norm_inf(Vec v)
        {
            PetscReal value = 0.0;
            petsc_check(
                VecNorm(v, NORM_INFINITY, &value),
                "VecNorm(NORM_INFINITY)");
            return static_cast<Real>(value);
        }

        void true_residual(Mat A, Vec x, Vec b, Vec r)
        {
            petsc_check(MatMult(A, x, r), "MatMult(residual)");
            petsc_check(VecAYPX(r, -1.0, b), "VecAYPX(residual)");
        }

        struct PetscResources
        {
            Mat A = nullptr;
            Vec b = nullptr;
            Vec x = nullptr;
            Vec r = nullptr;
            Vec fixed_values = nullptr;
            Vec shift = nullptr;
            IS fixed_rows = nullptr;
            KSP ksp = nullptr;

            ~PetscResources()
            {
                PetscBool finalised = PETSC_FALSE;
                PetscFinalized(&finalised);

                if (finalised)
                {
                    return;
                }

                if (ksp)
                {
                    KSPDestroy(&ksp);
                }
                if (fixed_rows)
                {
                    ISDestroy(&fixed_rows);
                }
                if (shift)
                {
                    VecDestroy(&shift);
                }
                if (fixed_values)
                {
                    VecDestroy(&fixed_values);
                }
                if (r)
                {
                    VecDestroy(&r);
                }
                if (x)
                {
                    VecDestroy(&x);
                }
                if (b)
                {
                    VecDestroy(&b);
                }
                if (A)
                {
                    MatDestroy(&A);
                }
            }
        };
    }

    ConjugateGradient::Result
    PetscPressureSolver::solveDistributedPressure(CBSStateSI& s)
    {
#ifdef CBS3D_USE_MPI
        if (!s.mpi_enabled)
        {
            throw std::runtime_error(
                "solveDistributedPressure requires more than one MPI rank");
        }

        if (s.cfg.ndim != 3 ||
            s.cfg.nep != 4 ||
            s.cfg.gsdim != 6)
        {
            throw std::runtime_error(
                "solveDistributedPressure requires P1 tetrahedra");
        }

        PetscBool petsc_initialised = PETSC_FALSE;
        petsc_check(PetscInitialized(&petsc_initialised), "PetscInitialized");

        const bool started_petsc_here =
            petsc_initialised == PETSC_FALSE;

        if (started_petsc_here)
        {
            petsc_check(
                PetscInitializeNoArguments(),
                "PetscInitializeNoArguments");
        }

        ConjugateGradient::Result result;

        {
            PetscResources p;

            std::vector<Int> active_owned;
            active_owned.reserve(s.owned_nodes.size());

            for (const Int ip : s.owned_nodes)
            {
                if (pressure_active(s, ip))
                {
                    active_owned.push_back(ip);
                }
            }

            const long long nlocal64 =
                static_cast<long long>(active_owned.size());

            long long nglobal64 = 0;
            long long offset64 = 0;

            mpi_check(
                MPI_Allreduce(
                    &nlocal64,
                    &nglobal64,
                    1,
                    MPI_LONG_LONG,
                    MPI_SUM,
                    MPI_COMM_WORLD),
                "MPI_Allreduce pressure DOFs");

            mpi_check(
                MPI_Exscan(
                    &nlocal64,
                    &offset64,
                    1,
                    MPI_LONG_LONG,
                    MPI_SUM,
                    MPI_COMM_WORLD),
                "MPI_Exscan pressure offset");

            if (s.mpi_rank == 0)
            {
                offset64 = 0;
            }

            if (nglobal64 <= 0 ||
                nglobal64 >
                    static_cast<long long>(
                        std::numeric_limits<Int>::max()))
            {
                throw std::runtime_error(
                    "Invalid distributed pressure DOF count");
            }

            Array1D<Int> pressure_dof(s.cfg.npoin);
            pressure_dof.fill(-1);

            for (std::size_t k = 0;
                 k < active_owned.size();
                 ++k)
            {
                pressure_dof(active_owned[k]) =
                    static_cast<Int>(
                        offset64 + static_cast<long long>(k));
            }

            HaloExchange::broadcastOwnedToGhosts(
                pressure_dof,
                s.partition_metadata,
                MPI_COMM_WORLD);

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                if (pressure_active(s, ip) &&
                    pressure_dof(ip) < 0)
                {
                    throw std::runtime_error(
                        "Active pressure node has no PETSc DOF");
                }
            }

            const PetscInt nlocal =
                static_cast<PetscInt>(nlocal64);

            const PetscInt nglobal =
                static_cast<PetscInt>(nglobal64);

            petsc_check(
                MatCreateAIJ(
                    MPI_COMM_WORLD,
                    nlocal,
                    nlocal,
                    nglobal,
                    nglobal,
                    32,
                    nullptr,
                    32,
                    nullptr,
                    &p.A),
                "MatCreateAIJ");

            petsc_check(
                MatSetOption(
                    p.A,
                    MAT_NEW_NONZERO_ALLOCATION_ERR,
                    PETSC_FALSE),
                "MatSetOption(allocation)");

            petsc_check(
                MatSetOption(
                    p.A,
                    MAT_IGNORE_ZERO_ENTRIES,
                    PETSC_TRUE),
                "MatSetOption(ignore zero)");

            constexpr Int pair_a[6] = {1, 1, 1, 2, 2, 3};
            constexpr Int pair_b[6] = {2, 3, 4, 3, 4, 4};

            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                if (s.mat_elem(ie) != 0)
                {
                    continue;
                }

                const Real dt = s.delte(ie);

                if (dt <= 0.0 || !std::isfinite(dt))
                {
                    throw std::runtime_error(
                        "Invalid pressure element timestep");
                }

                PetscInt rows[4] = {};
                PetscScalar Ke[16] = {};

                for (Int a = 1; a <= 4; ++a)
                {
                    const Int ip = s.intma(a, ie);
                    rows[a - 1] =
                        static_cast<PetscInt>(pressure_dof(ip));

                    const Real value =
                        dt *
                        s.pdiagE(
                            element_node_index(s, ie, a));

                    if (!std::isfinite(value) ||
                        value <= 0.0)
                    {
                        throw std::runtime_error(
                            "Invalid element pressure diagonal");
                    }

                    Ke[(a - 1) * 4 + (a - 1)] =
                        static_cast<PetscScalar>(value);
                }

                for (Int pair = 0; pair < 6; ++pair)
                {
                    const Real value =
                        dt *
                        s.gstifE(
                            offdiag_index(s, ie, pair + 1));

                    if (!std::isfinite(value))
                    {
                        throw std::runtime_error(
                            "Invalid element pressure coupling");
                    }

                    const Int a = pair_a[pair] - 1;
                    const Int b = pair_b[pair] - 1;

                    Ke[a * 4 + b] =
                        static_cast<PetscScalar>(value);

                    Ke[b * 4 + a] =
                        static_cast<PetscScalar>(value);
                }

                petsc_check(
                    MatSetValues(
                        p.A,
                        4,
                        rows,
                        4,
                        rows,
                        Ke,
                        ADD_VALUES),
                    "MatSetValues(element)");
            }

            petsc_check(
                MatAssemblyBegin(p.A, MAT_FINAL_ASSEMBLY),
                "MatAssemblyBegin");

            petsc_check(
                MatAssemblyEnd(p.A, MAT_FINAL_ASSEMBLY),
                "MatAssemblyEnd");

            petsc_check(
                MatCreateVecs(p.A, &p.x, &p.b),
                "MatCreateVecs");

            petsc_check(
                VecDuplicate(p.b, &p.r),
                "VecDuplicate(residual)");

            petsc_check(
                VecDuplicate(p.b, &p.fixed_values),
                "VecDuplicate(fixed values)");

            petsc_check(
                VecDuplicate(p.b, &p.shift),
                "VecDuplicate(shift)");

            PetscInt begin = 0;
            PetscInt end = 0;

            petsc_check(
                VecGetOwnershipRange(p.x, &begin, &end),
                "VecGetOwnershipRange");

            if (begin != static_cast<PetscInt>(offset64) ||
                end - begin != nlocal)
            {
                throw std::runtime_error(
                    "PETSc and CBS pressure ownership disagree");
            }

            petsc_check(
                VecSet(p.fixed_values, 0.0),
                "VecSet(fixed values)");

            std::vector<PetscInt> fixed_dofs;

            for (const Int ip : active_owned)
            {
                if (s.node_pressure_fixed(ip) == 0)
                {
                    continue;
                }

                const PetscInt dof =
                    static_cast<PetscInt>(pressure_dof(ip));

                fixed_dofs.push_back(dof);

                petsc_check(
                    VecSetValue(
                        p.fixed_values,
                        dof,
                        static_cast<PetscScalar>(
                            s.cfg.outlet_pressure_gauge),
                        INSERT_VALUES),
                    "VecSetValue(fixed value)");
            }

            assemble(p.fixed_values);

            petsc_check(
                MatMult(p.A, p.fixed_values, p.shift),
                "MatMult(Dirichlet shift)");

            petsc_check(
                ISCreateGeneral(
                    MPI_COMM_WORLD,
                    static_cast<PetscInt>(fixed_dofs.size()),
                    fixed_dofs.empty()
                        ? nullptr
                        : fixed_dofs.data(),
                    PETSC_COPY_VALUES,
                    &p.fixed_rows),
                "ISCreateGeneral(fixed rows)");

            petsc_check(
                MatZeroRowsColumnsIS(
                    p.A,
                    p.fixed_rows,
                    1.0,
                    nullptr,
                    nullptr),
                "MatZeroRowsColumnsIS");

            petsc_check(
                MatSetOption(
                    p.A,
                    MAT_SYMMETRIC,
                    PETSC_TRUE),
                "MatSetOption(symmetric)");

            petsc_check(
                MatSetOption(
                    p.A,
                    MAT_SPD,
                    PETSC_TRUE),
                "MatSetOption(SPD)");

            petsc_check(VecSet(p.b, 0.0), "VecSet(rhs)");
            petsc_check(VecSet(p.x, 0.0), "VecSet(initial pressure)");

            for (const Int ip : active_owned)
            {
                const PetscInt dof =
                    static_cast<PetscInt>(pressure_dof(ip));

                if (!std::isfinite(s.rhs1(ip)) ||
                    !std::isfinite(s.pres(ip)))
                {
                    throw std::runtime_error(
                        "Non-finite owned pressure data");
                }

                petsc_check(
                    VecSetValue(
                        p.b,
                        dof,
                        static_cast<PetscScalar>(s.rhs1(ip)),
                        INSERT_VALUES),
                    "VecSetValue(rhs)");

                petsc_check(
                    VecSetValue(
                        p.x,
                        dof,
                        static_cast<PetscScalar>(s.pres(ip)),
                        INSERT_VALUES),
                    "VecSetValue(initial pressure)");
            }

            assemble(p.b);
            assemble(p.x);

            petsc_check(
                VecAXPY(p.b, -1.0, p.shift),
                "VecAXPY(Dirichlet shift)");

            for (const Int ip : active_owned)
            {
                if (s.node_pressure_fixed(ip) == 0)
                {
                    continue;
                }

                const PetscInt dof =
                    static_cast<PetscInt>(pressure_dof(ip));

                const PetscScalar value =
                    static_cast<PetscScalar>(
                        s.cfg.outlet_pressure_gauge);

                petsc_check(
                    VecSetValue(
                        p.b,
                        dof,
                        value,
                        INSERT_VALUES),
                    "VecSetValue(fixed rhs)");

                petsc_check(
                    VecSetValue(
                        p.x,
                        dof,
                        value,
                        INSERT_VALUES),
                    "VecSetValue(fixed initial pressure)");
            }

            assemble(p.b);
            assemble(p.x);

            petsc_check(
                KSPCreate(MPI_COMM_WORLD, &p.ksp),
                "KSPCreate");

            petsc_check(
                KSPSetOperators(p.ksp, p.A, p.A),
                "KSPSetOperators");

            petsc_check(
                KSPSetType(p.ksp, KSPCG),
                "KSPSetType(KSPCG)");

            PC pc = nullptr;

            petsc_check(
                KSPGetPC(p.ksp, &pc),
                "KSPGetPC");

#if defined(PETSC_HAVE_HYPRE)
            petsc_check(
                PCSetType(pc, PCHYPRE),
                "PCSetType(PCHYPRE)");

            petsc_check(
                PCHYPRESetType(pc, "boomeramg"),
                "PCHYPRESetType(boomeramg)");
#else
            petsc_check(
                PCSetType(pc, PCGAMG),
                "PCSetType(PCGAMG)");
#endif

            petsc_check(
                KSPSetTolerances(
                    p.ksp,
                    static_cast<PetscReal>(s.cfg.relToler),
                    static_cast<PetscReal>(s.cfg.absToler),
                    PETSC_DEFAULT,
                    static_cast<PetscInt>(s.cfg.cg_max_iter)),
                "KSPSetTolerances");

            petsc_check(
                KSPSetInitialGuessNonzero(
                    p.ksp,
                    PETSC_TRUE),
                "KSPSetInitialGuessNonzero");

            petsc_check(
                KSPSetErrorIfNotConverged(
                    p.ksp,
                    PETSC_FALSE),
                "KSPSetErrorIfNotConverged");

            petsc_check(
                KSPSetFromOptions(p.ksp),
                "KSPSetFromOptions");

            true_residual(p.A, p.x, p.b, p.r);

            result.initial_l2 = norm2(p.r);
            result.final_l2 = result.initial_l2;
            result.final_relative_l2 =
                result.initial_l2 > 1.0e-300
                    ? 1.0
                    : 0.0;
            result.final_max_abs = norm_inf(p.r);

            if (result.initial_l2 > 1.0e-300)
            {
                petsc_check(
                    KSPSolve(p.ksp, p.b, p.x),
                    "KSPSolve");

                PetscInt iterations = 0;
                KSPConvergedReason reason =
                    KSP_CONVERGED_ITERATING;

                petsc_check(
                    KSPGetIterationNumber(
                        p.ksp,
                        &iterations),
                    "KSPGetIterationNumber");

                petsc_check(
                    KSPGetConvergedReason(
                        p.ksp,
                        &reason),
                    "KSPGetConvergedReason");

                true_residual(p.A, p.x, p.b, p.r);

                result.iterations =
                    static_cast<Int>(iterations);

                result.final_l2 = norm2(p.r);
                result.final_relative_l2 =
                    result.final_l2 / result.initial_l2;
                result.final_max_abs = norm_inf(p.r);
                result.converged = reason > 0;
            }
            else
            {
                result.converged = true;
            }

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                if (!pressure_active(s, ip))
                {
                    s.pres(ip) = 0.0;
                    s.rhs1(ip) = 0.0;
                }
            }

            const PetscScalar* values = nullptr;

            petsc_check(
                VecGetArrayRead(p.x, &values),
                "VecGetArrayRead");

            for (std::size_t k = 0;
                 k < active_owned.size();
                 ++k)
            {
                const Int ip = active_owned[k];

                s.pres(ip) =
                    static_cast<Real>(
                        PetscRealPart(values[k]));

                if (s.node_pressure_fixed(ip) != 0)
                {
                    s.pres(ip) =
                        s.cfg.outlet_pressure_gauge;
                }
            }

            petsc_check(
                VecRestoreArrayRead(p.x, &values),
                "VecRestoreArrayRead");
        }

        if (started_petsc_here)
        {
            petsc_check(PetscFinalize(), "PetscFinalize");
        }

        return result;
#else
        (void)s;
        throw std::runtime_error(
            "solveDistributedPressure requires an MPI-enabled build");
#endif
    }
}

#endif
