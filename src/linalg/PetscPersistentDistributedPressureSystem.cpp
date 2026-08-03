#include "cbs/linalg/PetscPersistentDistributedPressureSystem.hpp"

#include <stdexcept>

#if defined(CBS3D_USE_MPI) && defined(CBS3D_USE_PETSC)

#include "cbs/parallel/HaloExchange.hpp"

#include <mpi.h>
#include <petscksp.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace cbs
{
    namespace
    {
        void check_mpi(const int error_code, const char* operation)
        {
            if (error_code != MPI_SUCCESS)
            {
                throw std::runtime_error(
                    std::string(
                        "PetscPersistentDistributedPressureSystem MPI failure in ")
                    + operation);
            }
        }

        void check_petsc(
            const PetscErrorCode error_code,
            const char* operation)
        {
            if (error_code != 0)
            {
                throw std::runtime_error(
                    std::string(
                        "PetscPersistentDistributedPressureSystem PETSc failure in ")
                    + operation
                    + ", ierr="
                    + std::to_string(static_cast<int>(error_code)));
            }
        }

        bool pressure_active(
            const CBSStateSI& s,
            const Int ip)
        {
            return
                (s.node_material_mask(ip) &
                 CBSStateSI::node_touches_fluid) != 0;
        }

        Int element_node_index(
            const CBSStateSI& s,
            const Int ie,
            const Int local_node)
        {
            return (ie - 1) * s.cfg.nep + local_node;
        }

        Int offdiag_index(
            const CBSStateSI& s,
            const Int ie,
            const Int pair_id)
        {
            return (ie - 1) * s.cfg.gsdim + pair_id;
        }

        void assemble_vector(Vec vector)
        {
            check_petsc(
                VecAssemblyBegin(vector),
                "VecAssemblyBegin");

            check_petsc(
                VecAssemblyEnd(vector),
                "VecAssemblyEnd");
        }

        Real vector_l2_norm(Vec vector)
        {
            PetscReal norm = 0.0;

            check_petsc(
                VecNorm(vector, NORM_2, &norm),
                "VecNorm(NORM_2)");

            return static_cast<Real>(norm);
        }

        Real vector_max_norm(Vec vector)
        {
            PetscReal norm = 0.0;

            check_petsc(
                VecNorm(vector, NORM_INFINITY, &norm),
                "VecNorm(NORM_INFINITY)");

            return static_cast<Real>(norm);
        }

        void compute_true_residual(
            Mat matrix,
            Vec solution,
            Vec rhs,
            Vec residual)
        {
            check_petsc(
                MatMult(matrix, solution, residual),
                "MatMult(residual)");

            check_petsc(
                VecAYPX(residual, -1.0, rhs),
                "VecAYPX(residual)");
        }

        Real prescribed_pressure_value(
            const CBSStateSI& s,
            const Int ip)
        {
            for (Int i = 1; i <= s.cfg.bc_fixed; ++i)
            {
                if (s.bc_list(i) == ip)
                {
                    const Real value = s.bc_values(i);

                    if (!std::isfinite(value))
                    {
                        throw std::runtime_error(
                            "Non-finite prescribed pressure value");
                    }

                    return value;
                }
            }

            throw std::runtime_error(
                "Fixed pressure node is missing from the local pressure list");
        }
    }

    struct PetscPersistentDistributedPressureSystem::Impl
    {
        bool ready = false;
        bool started_petsc_here = false;

        Int npoin = 0;
        Int nelem = 0;

        PetscInt local_dofs = 0;
        PetscInt global_dofs = 0;
        PetscInt ownership_begin = 0;

        std::vector<Int> active_owned_nodes;
        std::vector<Real> active_owned_fixed_values;
        Array1D<Int> pressure_dof;

        Mat matrix = nullptr;
        Vec rhs = nullptr;
        Vec solution = nullptr;
        Vec residual = nullptr;
        Vec fixed_values = nullptr;
        Vec fixed_shift = nullptr;
        IS fixed_rows = nullptr;
        KSP ksp = nullptr;

        void destroy_objects() noexcept
        {
            PetscBool finalised = PETSC_FALSE;
            PetscFinalized(&finalised);

            if (finalised)
            {
                matrix = nullptr;
                rhs = nullptr;
                solution = nullptr;
                residual = nullptr;
                fixed_values = nullptr;
                fixed_shift = nullptr;
                fixed_rows = nullptr;
                ksp = nullptr;
                ready = false;
                return;
            }

            if (ksp != nullptr)
            {
                KSPDestroy(&ksp);
            }
            if (fixed_rows != nullptr)
            {
                ISDestroy(&fixed_rows);
            }
            if (fixed_shift != nullptr)
            {
                VecDestroy(&fixed_shift);
            }
            if (fixed_values != nullptr)
            {
                VecDestroy(&fixed_values);
            }
            if (residual != nullptr)
            {
                VecDestroy(&residual);
            }
            if (solution != nullptr)
            {
                VecDestroy(&solution);
            }
            if (rhs != nullptr)
            {
                VecDestroy(&rhs);
            }
            if (matrix != nullptr)
            {
                MatDestroy(&matrix);
            }

            ready = false;
        }
    };

    PetscPersistentDistributedPressureSystem::
    PetscPersistentDistributedPressureSystem()
        : impl_(new Impl())
    {
    }

    PetscPersistentDistributedPressureSystem::
    ~PetscPersistentDistributedPressureSystem()
    {
        shutdown();
        delete impl_;
        impl_ = nullptr;
    }

    bool PetscPersistentDistributedPressureSystem::ready() const noexcept
    {
        return impl_ != nullptr && impl_->ready;
    }

    void PetscPersistentDistributedPressureSystem::initialise(CBSStateSI& s)
    {
        if (impl_ == nullptr)
        {
            throw std::runtime_error(
                "Persistent distributed pressure implementation is unavailable");
        }

        if (impl_->ready)
        {
            if (impl_->npoin == s.cfg.npoin &&
                impl_->nelem == s.cfg.nelem)
            {
                return;
            }

            throw std::runtime_error(
                "Persistent pressure system already represents another mesh");
        }

        if (!s.mpi_enabled)
        {
            throw std::runtime_error(
                "Persistent distributed pressure requires more than one MPI rank");
        }

        if (s.cfg.ndim != 3 ||
            s.cfg.nep != 4 ||
            s.cfg.gsdim != 6)
        {
            throw std::runtime_error(
                "Persistent distributed pressure requires P1 tetrahedra");
        }

        PetscBool petsc_initialised = PETSC_FALSE;

        check_petsc(
            PetscInitialized(&petsc_initialised),
            "PetscInitialized");

        if (petsc_initialised == PETSC_FALSE)
        {
            check_petsc(
                PetscInitializeNoArguments(),
                "PetscInitializeNoArguments");

            impl_->started_petsc_here = true;
        }

        impl_->npoin = s.cfg.npoin;
        impl_->nelem = s.cfg.nelem;

        impl_->active_owned_nodes.clear();
        impl_->active_owned_nodes.reserve(s.owned_nodes.size());

        for (const Int ip : s.owned_nodes)
        {
            if (pressure_active(s, ip))
            {
                impl_->active_owned_nodes.push_back(ip);
            }
        }

        const long long local_count =
            static_cast<long long>(impl_->active_owned_nodes.size());

        long long global_count = 0;
        long long offset = 0;

        check_mpi(
            MPI_Allreduce(
                &local_count,
                &global_count,
                1,
                MPI_LONG_LONG,
                MPI_SUM,
                MPI_COMM_WORLD),
            "MPI_Allreduce pressure DOF count");

        check_mpi(
            MPI_Exscan(
                &local_count,
                &offset,
                1,
                MPI_LONG_LONG,
                MPI_SUM,
                MPI_COMM_WORLD),
            "MPI_Exscan pressure ownership offset");

        if (s.mpi_rank == 0)
        {
            offset = 0;
        }

        if (global_count <= 0 ||
            global_count >
                static_cast<long long>(
                    std::numeric_limits<Int>::max()) ||
            global_count >
                static_cast<long long>(
                    std::numeric_limits<PetscInt>::max()))
        {
            throw std::runtime_error(
                "Invalid persistent distributed pressure DOF count");
        }

        impl_->local_dofs = static_cast<PetscInt>(local_count);
        impl_->global_dofs = static_cast<PetscInt>(global_count);
        impl_->ownership_begin = static_cast<PetscInt>(offset);

        impl_->pressure_dof.resize(s.cfg.npoin);
        impl_->pressure_dof.fill(-1);

        for (std::size_t k = 0;
             k < impl_->active_owned_nodes.size();
             ++k)
        {
            impl_->pressure_dof(
                impl_->active_owned_nodes[k]) =
                static_cast<Int>(
                    offset + static_cast<long long>(k));
        }

        HaloExchange::broadcastOwnedToGhosts(
            impl_->pressure_dof,
            s.partition_metadata,
            MPI_COMM_WORLD);

        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            if (pressure_active(s, ip) &&
                impl_->pressure_dof(ip) < 0)
            {
                throw std::runtime_error(
                    "Active pressure node has no persistent PETSc DOF");
            }
        }

        check_petsc(
            MatCreateAIJ(
                MPI_COMM_WORLD,
                impl_->local_dofs,
                impl_->local_dofs,
                impl_->global_dofs,
                impl_->global_dofs,
                48,
                nullptr,
                48,
                nullptr,
                &impl_->matrix),
            "MatCreateAIJ");

        check_petsc(
            MatSetOption(
                impl_->matrix,
                MAT_NEW_NONZERO_ALLOCATION_ERR,
                PETSC_FALSE),
            "MatSetOption(allocation)");

        check_petsc(
            MatSetOption(
                impl_->matrix,
                MAT_IGNORE_ZERO_ENTRIES,
                PETSC_TRUE),
            "MatSetOption(ignore zero entries)");

        constexpr Int pair_a[6] = {1, 1, 1, 2, 2, 3};
        constexpr Int pair_b[6] = {2, 3, 4, 3, 4, 4};

        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            if (s.mat_elem(ie) != 0)
            {
                continue;
            }

            PetscInt rows[4] = {};
            PetscScalar element_matrix[16] = {};

            for (Int a = 1; a <= s.cfg.nep; ++a)
            {
                const Int ip = s.intma(a, ie);
                const Int dof = impl_->pressure_dof(ip);

                if (dof < 0)
                {
                    throw std::runtime_error(
                        "Fluid element contains an inactive pressure node");
                }

                rows[a - 1] = static_cast<PetscInt>(dof);

                const Real diagonal =
                    s.pdiagE(element_node_index(s, ie, a));

                if (diagonal <= 0.0 ||
                    !std::isfinite(diagonal))
                {
                    throw std::runtime_error(
                        "Invalid geometry-only pressure diagonal");
                }

                element_matrix[(a - 1) * 4 + (a - 1)] =
                    static_cast<PetscScalar>(diagonal);
            }

            for (Int pair = 0; pair < s.cfg.gsdim; ++pair)
            {
                const Real coupling =
                    s.gstifE(
                        offdiag_index(s, ie, pair + 1));

                if (!std::isfinite(coupling))
                {
                    throw std::runtime_error(
                        "Invalid geometry-only pressure coupling");
                }

                const Int a = pair_a[pair] - 1;
                const Int b = pair_b[pair] - 1;

                element_matrix[a * 4 + b] =
                    static_cast<PetscScalar>(coupling);

                element_matrix[b * 4 + a] =
                    static_cast<PetscScalar>(coupling);
            }

            check_petsc(
                MatSetValues(
                    impl_->matrix,
                    4,
                    rows,
                    4,
                    rows,
                    element_matrix,
                    ADD_VALUES),
                "MatSetValues(element pressure matrix)");
        }

        check_petsc(
            MatAssemblyBegin(
                impl_->matrix,
                MAT_FINAL_ASSEMBLY),
            "MatAssemblyBegin");

        check_petsc(
            MatAssemblyEnd(
                impl_->matrix,
                MAT_FINAL_ASSEMBLY),
            "MatAssemblyEnd");

        check_petsc(
            MatCreateVecs(
                impl_->matrix,
                &impl_->solution,
                &impl_->rhs),
            "MatCreateVecs");

        check_petsc(
            VecDuplicate(
                impl_->rhs,
                &impl_->residual),
            "VecDuplicate(residual)");

        check_petsc(
            VecDuplicate(
                impl_->rhs,
                &impl_->fixed_values),
            "VecDuplicate(fixed values)");

        check_petsc(
            VecDuplicate(
                impl_->rhs,
                &impl_->fixed_shift),
            "VecDuplicate(fixed shift)");

        PetscInt ownership_begin = 0;
        PetscInt ownership_end = 0;

        check_petsc(
            VecGetOwnershipRange(
                impl_->solution,
                &ownership_begin,
                &ownership_end),
            "VecGetOwnershipRange");

        if (ownership_begin != impl_->ownership_begin ||
            ownership_end - ownership_begin != impl_->local_dofs)
        {
            throw std::runtime_error(
                "PETSc and CBS persistent pressure ownership disagree");
        }

        check_petsc(
            VecSet(impl_->fixed_values, 0.0),
            "VecSet(fixed values)");

        std::vector<PetscInt> fixed_dofs;
        impl_->active_owned_fixed_values.assign(
            impl_->active_owned_nodes.size(),
            0.0);

        for (std::size_t k = 0;
             k < impl_->active_owned_nodes.size();
             ++k)
        {
            const Int ip = impl_->active_owned_nodes[k];

            if (s.node_pressure_fixed(ip) == 0)
            {
                continue;
            }

            const PetscInt dof =
                static_cast<PetscInt>(
                    impl_->pressure_dof(ip));

            const Real value =
                prescribed_pressure_value(s, ip);

            fixed_dofs.push_back(dof);
            impl_->active_owned_fixed_values[k] = value;

            check_petsc(
                VecSetValue(
                    impl_->fixed_values,
                    dof,
                    static_cast<PetscScalar>(value),
                    INSERT_VALUES),
                "VecSetValue(fixed pressure value)");
        }

        assemble_vector(impl_->fixed_values);

        check_petsc(
            MatMult(
                impl_->matrix,
                impl_->fixed_values,
                impl_->fixed_shift),
            "MatMult(fixed pressure shift)");

        check_petsc(
            ISCreateGeneral(
                MPI_COMM_WORLD,
                static_cast<PetscInt>(fixed_dofs.size()),
                fixed_dofs.empty()
                    ? nullptr
                    : fixed_dofs.data(),
                PETSC_COPY_VALUES,
                &impl_->fixed_rows),
            "ISCreateGeneral(fixed pressure rows)");

        check_petsc(
            MatZeroRowsColumnsIS(
                impl_->matrix,
                impl_->fixed_rows,
                1.0,
                nullptr,
                nullptr),
            "MatZeroRowsColumnsIS");

        check_petsc(
            MatSetOption(
                impl_->matrix,
                MAT_SYMMETRIC,
                PETSC_TRUE),
            "MatSetOption(MAT_SYMMETRIC)");

        check_petsc(
            MatSetOption(
                impl_->matrix,
                MAT_SPD,
                PETSC_TRUE),
            "MatSetOption(MAT_SPD)");

        check_petsc(
            KSPCreate(
                MPI_COMM_WORLD,
                &impl_->ksp),
            "KSPCreate");

        check_petsc(
            KSPSetOperators(
                impl_->ksp,
                impl_->matrix,
                impl_->matrix),
            "KSPSetOperators");

        check_petsc(
            KSPSetType(
                impl_->ksp,
                KSPCG),
            "KSPSetType(KSPCG)");

        PC preconditioner = nullptr;

        check_petsc(
            KSPGetPC(
                impl_->ksp,
                &preconditioner),
            "KSPGetPC");

#if defined(PETSC_HAVE_HYPRE)
        check_petsc(
            PCSetType(
                preconditioner,
                PCHYPRE),
            "PCSetType(PCHYPRE)");

        check_petsc(
            PCHYPRESetType(
                preconditioner,
                "boomeramg"),
            "PCHYPRESetType(boomeramg)");
#else
        check_petsc(
            PCSetType(
                preconditioner,
                PCGAMG),
            "PCSetType(PCGAMG)");
#endif

        check_petsc(
            KSPSetInitialGuessNonzero(
                impl_->ksp,
                PETSC_TRUE),
            "KSPSetInitialGuessNonzero");

        check_petsc(
            KSPSetErrorIfNotConverged(
                impl_->ksp,
                PETSC_FALSE),
            "KSPSetErrorIfNotConverged");

        check_petsc(
            KSPSetFromOptions(impl_->ksp),
            "KSPSetFromOptions");

        check_petsc(
            KSPSetUp(impl_->ksp),
            "KSPSetUp(persistent pressure)");

        impl_->ready = true;
    }

    ConjugateGradient::Result
    PetscPersistentDistributedPressureSystem::solve(CBSStateSI& s)
    {
        if (impl_ == nullptr || !impl_->ready)
        {
            throw std::runtime_error(
                "Persistent distributed pressure system is not initialised");
        }

        if (impl_->npoin != s.cfg.npoin ||
            impl_->nelem != s.cfg.nelem)
        {
            throw std::runtime_error(
                "Persistent distributed pressure mesh changed after setup");
        }

        const Real local_dt = s.cfg.dtreal;

        if (local_dt <= 0.0 ||
            !std::isfinite(local_dt))
        {
            throw std::runtime_error(
                "Persistent pressure solve received an invalid global timestep");
        }

        Real minimum_dt = 0.0;
        Real maximum_dt = 0.0;

        check_mpi(
            MPI_Allreduce(
                &local_dt,
                &minimum_dt,
                1,
                MPI_DOUBLE,
                MPI_MIN,
                MPI_COMM_WORLD),
            "MPI_Allreduce minimum pressure timestep");

        check_mpi(
            MPI_Allreduce(
                &local_dt,
                &maximum_dt,
                1,
                MPI_DOUBLE,
                MPI_MAX,
                MPI_COMM_WORLD),
            "MPI_Allreduce maximum pressure timestep");

        const Real timestep_scale =
            std::max(1.0, std::abs(local_dt));

        if (std::abs(maximum_dt - minimum_dt) >
            1.0e-13 * timestep_scale)
        {
            throw std::runtime_error(
                "Persistent pressure system requires one communicator-wide timestep");
        }

        check_petsc(
            VecSet(impl_->rhs, 0.0),
            "VecSet(rhs)");

        check_petsc(
            VecSet(impl_->solution, 0.0),
            "VecSet(initial pressure)");

        for (std::size_t k = 0;
             k < impl_->active_owned_nodes.size();
             ++k)
        {
            const Int ip = impl_->active_owned_nodes[k];
            const PetscInt dof =
                static_cast<PetscInt>(
                    impl_->pressure_dof(ip));

            if (!std::isfinite(s.rhs1(ip)) ||
                !std::isfinite(s.pres(ip)))
            {
                throw std::runtime_error(
                    "Non-finite owned pressure data in persistent solve");
            }

            if (s.node_pressure_fixed(ip) != 0)
            {
                const Real fixed_value =
                    impl_->active_owned_fixed_values[k];

                check_petsc(
                    VecSetValue(
                        impl_->rhs,
                        dof,
                        static_cast<PetscScalar>(fixed_value),
                        INSERT_VALUES),
                    "VecSetValue(fixed rhs)");

                check_petsc(
                    VecSetValue(
                        impl_->solution,
                        dof,
                        static_cast<PetscScalar>(fixed_value),
                        INSERT_VALUES),
                    "VecSetValue(fixed initial pressure)");
            }
            else
            {
                const Real scaled_rhs =
                    s.rhs1(ip) / local_dt;

                check_petsc(
                    VecSetValue(
                        impl_->rhs,
                        dof,
                        static_cast<PetscScalar>(scaled_rhs),
                        INSERT_VALUES),
                    "VecSetValue(free rhs)");

                check_petsc(
                    VecSetValue(
                        impl_->solution,
                        dof,
                        static_cast<PetscScalar>(s.pres(ip)),
                        INSERT_VALUES),
                    "VecSetValue(initial pressure)");
            }
        }

        assemble_vector(impl_->rhs);
        assemble_vector(impl_->solution);

        check_petsc(
            VecAXPY(
                impl_->rhs,
                -1.0,
                impl_->fixed_shift),
            "VecAXPY(fixed pressure shift)");

        for (std::size_t k = 0;
             k < impl_->active_owned_nodes.size();
             ++k)
        {
            const Int ip = impl_->active_owned_nodes[k];

            if (s.node_pressure_fixed(ip) == 0)
            {
                continue;
            }

            const PetscInt dof =
                static_cast<PetscInt>(
                    impl_->pressure_dof(ip));

            const PetscScalar fixed_value =
                static_cast<PetscScalar>(
                    impl_->active_owned_fixed_values[k]);

            check_petsc(
                VecSetValue(
                    impl_->rhs,
                    dof,
                    fixed_value,
                    INSERT_VALUES),
                "VecSetValue(restored fixed rhs)");

            check_petsc(
                VecSetValue(
                    impl_->solution,
                    dof,
                    fixed_value,
                    INSERT_VALUES),
                "VecSetValue(restored fixed initial pressure)");
        }

        assemble_vector(impl_->rhs);
        assemble_vector(impl_->solution);

        check_petsc(
            KSPSetTolerances(
                impl_->ksp,
                static_cast<PetscReal>(s.cfg.relToler),
                static_cast<PetscReal>(s.cfg.absToler / local_dt),
                PETSC_DEFAULT,
                static_cast<PetscInt>(s.cfg.cg_max_iter)),
            "KSPSetTolerances");

        ConjugateGradient::Result result;

        compute_true_residual(
            impl_->matrix,
            impl_->solution,
            impl_->rhs,
            impl_->residual);

        const Real scaled_initial_l2 =
            vector_l2_norm(impl_->residual);

        result.initial_l2 =
            local_dt * scaled_initial_l2;

        result.final_l2 = result.initial_l2;
        result.final_relative_l2 =
            scaled_initial_l2 > 1.0e-300
                ? 1.0
                : 0.0;

        result.final_max_abs =
            local_dt * vector_max_norm(impl_->residual);

        if (scaled_initial_l2 > 1.0e-300)
        {
            check_petsc(
                KSPSolve(
                    impl_->ksp,
                    impl_->rhs,
                    impl_->solution),
                "KSPSolve(persistent pressure)");

            PetscInt iterations = 0;
            KSPConvergedReason reason =
                KSP_CONVERGED_ITERATING;

            check_petsc(
                KSPGetIterationNumber(
                    impl_->ksp,
                    &iterations),
                "KSPGetIterationNumber");

            check_petsc(
                KSPGetConvergedReason(
                    impl_->ksp,
                    &reason),
                "KSPGetConvergedReason");

            compute_true_residual(
                impl_->matrix,
                impl_->solution,
                impl_->rhs,
                impl_->residual);

            const Real scaled_final_l2 =
                vector_l2_norm(impl_->residual);

            result.iterations =
                static_cast<Int>(iterations);

            result.final_l2 =
                local_dt * scaled_final_l2;

            result.final_relative_l2 =
                scaled_final_l2 / scaled_initial_l2;

            result.final_max_abs =
                local_dt * vector_max_norm(impl_->residual);

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

        const PetscScalar* local_values = nullptr;

        check_petsc(
            VecGetArrayRead(
                impl_->solution,
                &local_values),
            "VecGetArrayRead(solution)");

        for (std::size_t k = 0;
             k < impl_->active_owned_nodes.size();
             ++k)
        {
            const Int ip = impl_->active_owned_nodes[k];

            s.pres(ip) =
                static_cast<Real>(
                    PetscRealPart(local_values[k]));

            if (s.node_pressure_fixed(ip) != 0)
            {
                s.pres(ip) =
                    impl_->active_owned_fixed_values[k];
            }
        }

        check_petsc(
            VecRestoreArrayRead(
                impl_->solution,
                &local_values),
            "VecRestoreArrayRead(solution)");

        return result;
    }

    void PetscPersistentDistributedPressureSystem::shutdown() noexcept
    {
        if (impl_ == nullptr)
        {
            return;
        }

        const bool finalise_petsc =
            impl_->started_petsc_here;

        impl_->destroy_objects();

        if (finalise_petsc)
        {
            PetscBool finalised = PETSC_FALSE;
            PetscFinalized(&finalised);

            if (!finalised)
            {
                PetscFinalize();
            }
        }

        impl_->started_petsc_here = false;
        impl_->active_owned_nodes.clear();
        impl_->active_owned_fixed_values.clear();
        impl_->pressure_dof.resize(0);
        impl_->npoin = 0;
        impl_->nelem = 0;
        impl_->local_dofs = 0;
        impl_->global_dofs = 0;
        impl_->ownership_begin = 0;
    }
}

#else

namespace cbs
{
    struct PetscPersistentDistributedPressureSystem::Impl
    {
        bool ready = false;
    };

    PetscPersistentDistributedPressureSystem::
    PetscPersistentDistributedPressureSystem()
        : impl_(new Impl())
    {
    }

    PetscPersistentDistributedPressureSystem::
    ~PetscPersistentDistributedPressureSystem()
    {
        delete impl_;
        impl_ = nullptr;
    }

    bool PetscPersistentDistributedPressureSystem::ready() const noexcept
    {
        return false;
    }

    void PetscPersistentDistributedPressureSystem::initialise(CBSStateSI&)
    {
        throw std::runtime_error(
            "Persistent distributed pressure requires MPI and PETSc support");
    }

    ConjugateGradient::Result
    PetscPersistentDistributedPressureSystem::solve(CBSStateSI&)
    {
        throw std::runtime_error(
            "Persistent distributed pressure requires MPI and PETSc support");
    }

    void PetscPersistentDistributedPressureSystem::shutdown() noexcept
    {
    }
}

#endif
