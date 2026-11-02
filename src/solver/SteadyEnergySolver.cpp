//=============================================================================
// CBS3D++_SI
//
// FILE      : SteadyEnergySolver.cpp
// PURPOSE   : IMPLICIT STEADY CONJUGATE ENERGY SOLVE (SEE THE HEADER).
//=============================================================================

#include "cbs/solver/SteadyEnergySolver.hpp"

#include <stdexcept>

#if defined(CBS3D_USE_MPI) && defined(CBS3D_USE_PETSC)

#include "cbs/assembly/EnergyAssembly.hpp"
#include "cbs/assembly/TemperatureAFC.hpp"
#include "cbs/parallel/HaloExchange.hpp"

#include <mpi.h>
#include "cbs/linalg/PetscHeaders.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace cbs
{
    namespace
    {
        constexpr Int pair_a[6] = {1, 1, 1, 2, 2, 3};
        constexpr Int pair_b[6] = {2, 3, 4, 3, 4, 4};

        void check_mpi(const int error_code, const char* operation)
        {
            if (error_code != MPI_SUCCESS)
            {
                throw std::runtime_error(
                    std::string("SteadyEnergySolver MPI failure in ") + operation);
            }
        }

        void check_petsc(const PetscErrorCode error_code, const char* operation)
        {
            if (error_code != 0)
            {
                throw std::runtime_error(
                    std::string("SteadyEnergySolver PETSc failure in ") + operation +
                    ", ierr=" + std::to_string(static_cast<int>(error_code)));
            }
        }

        bool env_flag(const char* name)
        {
            const char* v = std::getenv(name);
            return v != nullptr && v[0] != '\0' && std::string(v) != "0";
        }

        Real env_real(const char* name, const Real fallback)
        {
            const char* v = std::getenv(name);
            if (v == nullptr || v[0] == '\0')
            {
                return fallback;
            }
            return std::stod(v);
        }

        Real global_max(Real local)
        {
            Real global = 0.0;
            check_mpi(MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX,
                                    MPI_COMM_WORLD), "Allreduce max");
            return global;
        }

        Real global_sum(Real local)
        {
            Real global = 0.0;
            check_mpi(MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM,
                                    MPI_COMM_WORLD), "Allreduce sum");
            return global;
        }

        // ================================
        // FULL EXPLICIT RESIDUAL R(T) = M T + c + q ON THE OWNERS, ASSEMBLED
        // WITH THE ORIGINAL STEP-4 CODE (HEAT FLUX ADDED FROM THE NODAL LOAD,
        // AS IN THE PRODUCTION LOOP). T IS TAKEN FROM s.temperature.
        // ================================
        void explicit_residual(CBSStateSI& s, const Array1D<Real>& heat_flux_load,
                               const bool apply_afc = false)
        {
            s.temperature1 = s.temperature;

            const Real saved_flux = s.cfg.heat_flux_bc;
            s.cfg.heat_flux_bc = 0.0;
            try
            {
                EnergyAssembly::assembleStep4Rhs(s);
            }
            catch (...)
            {
                s.cfg.heat_flux_bc = saved_flux;
                throw;
            }
            s.cfg.heat_flux_bc = saved_flux;

            if (apply_afc)
            {
                // Production Step 4 limits the volume residual first. The
                // prescribed face load must be added AFTER the limiter or its
                // low-order predictor/bounds differ from explicit AFC.
                TemperatureAFC::limitDistributedResidual(s);
            }
            else
            {
                HaloExchange::sumGhostContributionsToOwners(
                    s.rhs1, s.partition_metadata, MPI_COMM_WORLD);
            }

            for (const Int ip : s.owned_nodes)
            {
                s.rhs1(ip) += heat_flux_load(ip);
            }
        }

        // ================================
        // OWNER VECTOR <-> NODAL ARRAY
        // ================================
        void vec_from_nodes(Vec v, const CBSStateSI& s, const Array1D<Int>& dof,
                            const Array1D<Real>& nodal)
        {
            for (const Int ip : s.owned_nodes)
            {
                check_petsc(VecSetValue(v, static_cast<PetscInt>(dof(ip)),
                                        static_cast<PetscScalar>(nodal(ip)),
                                        INSERT_VALUES), "VecSetValue");
            }
            check_petsc(VecAssemblyBegin(v), "VecAssemblyBegin");
            check_petsc(VecAssemblyEnd(v), "VecAssemblyEnd");
        }

        void nodes_from_vec(Array1D<Real>& nodal, CBSStateSI& s,
                            const Array1D<Int>& dof, Vec v)
        {
            const PetscScalar* xa = nullptr;
            PetscInt begin = 0;
            PetscInt end = 0;
            check_petsc(VecGetOwnershipRange(v, &begin, &end), "VecGetOwnershipRange");
            check_petsc(VecGetArrayRead(v, &xa), "VecGetArrayRead");
            for (const Int ip : s.owned_nodes)
            {
                nodal(ip) = static_cast<Real>(xa[dof(ip) - begin]);
            }
            check_petsc(VecRestoreArrayRead(v, &xa), "VecRestoreArrayRead");
            HaloExchange::broadcastOwnedToGhosts(nodal, s.partition_metadata,
                                                 MPI_COMM_WORLD);
        }
    }

    SteadyEnergySolver::Result SteadyEnergySolver::solve(
        CBSStateSI& s,
        const Array1D<Int>& fixed,
        const Array1D<Real>& value,
        const Array1D<Real>& heat_flux_load)
    {
        if (!s.mpi_enabled || s.cfg.nep != 4 || s.cfg.ndim != 3)
        {
            throw std::runtime_error(
                "SteadyEnergySolver needs an MPI run on four-node tetrahedra");
        }

        const double start = MPI_Wtime();
        Result result;

        // AFC: solve the production explicit Step-4 limited residual R(T)=0
        // using a frozen low-order PETSc matrix as defect-correction
        // preconditioner. The limiter itself is NOT independently reimplemented.
        const bool afc = env_flag("CBS3D_TEMPERATURE_AFC");
        const bool low_order_only =
            env_flag("CBS3D_STEADY_ENERGY_LOW_ORDER_ONLY");
        // Genuine *physical-time* lumped-capacitance backward Euler for the
        // same frozen-velocity, graph-viscosity-stabilised CHT operator.
        // This is NOT the local pseudo-time relaxation used in explicit CBS.
        const Real be_dt =
            env_real("CBS3D_IMPLICIT_ENERGY_BE_DT", 0.0);
        const bool transient_be = be_dt > 0.0;
        const Real be_steps_input =
            env_real("CBS3D_IMPLICIT_ENERGY_BE_STEPS", 1.0);
        if (!std::isfinite(be_dt) || be_dt < 0.0 ||
            !std::isfinite(be_steps_input) || be_steps_input < 1.0 ||
            be_steps_input > 100000.0 ||
            std::floor(be_steps_input) != be_steps_input)
        {
            throw std::runtime_error(
                "Implicit physical-time BE: dt must be finite and >=0; "
                "steps must be an integer in [1,100000]");
        }
        const Int be_steps = static_cast<Int>(be_steps_input);
        const bool be_check_floor =
            std::getenv("CBS3D_IMPLICIT_ENERGY_BE_MIN_TEMPERATURE") != nullptr;
        const Real be_floor = be_check_floor
            ? env_real("CBS3D_IMPLICIT_ENERGY_BE_MIN_TEMPERATURE", 0.0)
            : 0.0;
        if (be_check_floor && !std::isfinite(be_floor))
        {
            throw std::runtime_error(
                "Implicit BE minimum temperature check must be finite");
        }
        if (transient_be && !low_order_only)
        {
            throw std::runtime_error(
                "Physical-time BE currently requires "
                "CBS3D_STEADY_ENERGY_LOW_ORDER_ONLY=1. "
                "High-order bounded transient AFC has NOT been implemented.");
        }

        // EnergyAssembly's characteristic streamline diffusion uses the
        // CBS pseudo-time delte_e in tau_e = 0.5*delte_e*rhoCp*V.
        // That pseudo-time has NO physical interpretation in BE. Remove
        // that contribution *only during BE energy assembly*: the low-order
        // graph viscosity remains the chosen spatial stabilization.
        // Restore the original element pseudo-timesteps on every exit.
        struct PhysicalBEOperatorGuard
        {
            CBSStateSI& state;
            bool active;
            Array1D<Real> original;
            PhysicalBEOperatorGuard(CBSStateSI& s, bool enabled)
                : state(s), active(enabled)
            {
                if (active)
                {
                    original = state.delte;
                    state.delte.fill(0.0);
                }
            }
            ~PhysicalBEOperatorGuard()
            {
                if (active)
                {
                    state.delte = original;
                }
            }
        };
        PhysicalBEOperatorGuard physical_be_operator(s, transient_be);
        if (low_order_only && !afc)
        {
            throw std::runtime_error(
                "CBS3D_STEADY_ENERGY_LOW_ORDER_ONLY requires "
                "CBS3D_TEMPERATURE_AFC=1 to construct the fluid graph viscosity");
        }
        const Real afc_tol = env_real("CBS3D_STEADY_ENERGY_TOL", 1.0e-6);
        const Real residual_tol = env_real("CBS3D_STEADY_ENERGY_RESIDUAL_TOL", 1.0e-6);
        if (!(residual_tol > 0.0) || !std::isfinite(residual_tol))
        {
            throw std::runtime_error("CBS3D_STEADY_ENERGY_RESIDUAL_TOL must be positive and finite");
        }
        const Int afc_maxit = static_cast<Int>(env_real("CBS3D_STEADY_ENERGY_MAXIT", 300.0));
        // OPTIONAL ANDERSON ACCELERATION OF THE AFC DEFECT-CORRECTION MAP.
        // ZERO RETAINS THE PREVIOUS DEFECT-CORRECTION ALGORITHM.
        const Real anderson_input =
            env_real("CBS3D_STEADY_ENERGY_ANDERSON_M", 0.0);
        if (!std::isfinite(anderson_input) || anderson_input < 0.0 ||
            anderson_input > 5.0 || std::floor(anderson_input) != anderson_input)
        {
            throw std::runtime_error(
                "CBS3D_STEADY_ENERGY_ANDERSON_M must be an integer in [0,5]");
        }
        const Int anderson_m = afc ? static_cast<Int>(anderson_input) : 0;
        // UNDER-RELAXATION OF THE DEFECT CORRECTION (THE ZALESAK LIMITER IS
        // NOT SMOOTH; 0.5 CONVERGED WHERE 1.0 CYCLED). NO AFC: ONE EXACT SOLVE.
        const Real afc_relax = afc ? env_real("CBS3D_STEADY_ENERGY_RELAX", 0.5) : 1.0;
        // Positive diagonal shift represents transient heat-capacity / delta t:
        // (A_low + gamma * C/dt) delta T = R_AFC(T). gamma=0 recovers the
        // unshifted steady preconditioner that amplified the U-channel defect.
        const Real pseudotime_shift =
            afc ? env_real("CBS3D_STEADY_ENERGY_PSEUDOTIME_SHIFT", 1.0) : 0.0;
        if (transient_be && pseudotime_shift != 0.0)
        {
            throw std::runtime_error(
                "Physical-time BE requires CBS3D_STEADY_ENERGY_PSEUDOTIME_SHIFT=0; "
                "the physical C/dt diagonal is assembled independently.");
        }
        if (!std::isfinite(pseudotime_shift) || pseudotime_shift < 0.0)
        {
            throw std::runtime_error(
                "CBS3D_STEADY_ENERGY_PSEUDOTIME_SHIFT must be finite and nonnegative");
        }
        const Real armijo_c1 =
            afc ? env_real("CBS3D_STEADY_ENERGY_ARMIJO_C1", 1.0e-5) : 1.0e-4;
        if (!std::isfinite(armijo_c1) || armijo_c1 < 0.0 || armijo_c1 >= 1.0)
        {
            throw std::runtime_error(
                "CBS3D_STEADY_ENERGY_ARMIJO_C1 must be finite and in [0,1)");
        }

        // ================================
        // 1. DOF MAP: ONE ROW PER GLOBAL NODE, OWNED NODES CONTIGUOUS
        // ================================
        long long local_count = static_cast<long long>(s.owned_nodes.size());
        long long global_count = 0;
        long long offset = 0;

        check_mpi(MPI_Allreduce(&local_count, &global_count, 1, MPI_LONG_LONG,
                                MPI_SUM, MPI_COMM_WORLD), "DOF count");
        check_mpi(MPI_Exscan(&local_count, &offset, 1, MPI_LONG_LONG,
                             MPI_SUM, MPI_COMM_WORLD), "DOF offset");
        if (s.mpi_rank == 0)
        {
            offset = 0;
        }

        Array1D<Int> dof(s.cfg.npoin);
        dof.fill(-1);
        {
            long long k = 0;
            for (const Int ip : s.owned_nodes)
            {
                dof(ip) = static_cast<Int>(offset + k);
                ++k;
            }
        }
        HaloExchange::broadcastOwnedToGhosts(dof, s.partition_metadata, MPI_COMM_WORLD);

        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            if (dof(ip) < 0)
            {
                throw std::runtime_error("SteadyEnergySolver: NODE WITHOUT A DOF");
            }
        }

        // ================================
        // 2. CONSTANT PART  c + q  =  R(T = 0)
        // ================================
        const Array1D<Real> saved_temperature = s.temperature;
        s.temperature.fill(0.0);
        explicit_residual(s, heat_flux_load);
        const Array1D<Real> constant_part = s.rhs1;
        s.temperature = saved_temperature;

        // ================================
        // 3. MATRIX  A = -(M + D) + gamma * diag(1/elcoe2p) (FREE ROWS),
        //    IDENTITY (DIRICHLET ROWS). gamma=0 is the original steady form.
        //    Positive gamma adds implicit pseudo-transient heat capacity.
        //    D = GRAPH VISCOSITY OF THE FLUID ELEMENTS WHEN AFC IS ON
        //    (d_ab = max(0, -M_ab, -M_ba), AS IN TemperatureAFC.cpp), WHICH
        //    MAKES THE LOW-ORDER OFF-DIAGONALS l_ab = M_ab + d_ab >= 0.
        //    The matrix is used only as a fixed defect-correction preconditioner.
        // ================================
        Mat A = nullptr;
        check_petsc(MatCreateAIJ(MPI_COMM_WORLD,
                                 static_cast<PetscInt>(local_count),
                                 static_cast<PetscInt>(local_count),
                                 static_cast<PetscInt>(global_count),
                                 static_cast<PetscInt>(global_count),
                                 48, nullptr, 48, nullptr, &A), "MatCreateAIJ");
        check_petsc(MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE),
                    "MatSetOption");

        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            Real M[5][5];
            EnergyAssembly::elementOperator(s, ie, M);

            if (afc && s.mat_elem(ie) == 0)
            {
                for (Int k = 0; k < 6; ++k)
                {
                    const Int a = pair_a[k];
                    const Int b = pair_b[k];
                    const Real d = std::max({0.0, -M[a][b], -M[b][a]});
                    M[a][b] += d;
                    M[b][a] += d;
                    M[a][a] -= d;
                    M[b][b] -= d;
                }
            }

            PetscInt rows[4];
            PetscInt cols[4];
            PetscScalar values[16];
            PetscInt nrow = 0;

            for (Int b = 1; b <= 4; ++b)
            {
                cols[b - 1] = static_cast<PetscInt>(dof(s.intma(b, ie)));
            }

            for (Int a = 1; a <= 4; ++a)
            {
                const Int ip = s.intma(a, ie);
                if (fixed(ip) != 0)
                {
                    continue;
                }
                rows[nrow] = static_cast<PetscInt>(dof(ip));
                for (Int b = 1; b <= 4; ++b)
                {
                    const Real entry = -M[a][b];
                    if (!std::isfinite(entry))
                    {
                        throw std::runtime_error(
                            "SteadyEnergySolver: NON-FINITE ELEMENT OPERATOR AT ELEMENT "
                            + std::to_string(ie));
                    }
                    values[nrow * 4 + (b - 1)] = static_cast<PetscScalar>(entry);
                }
                ++nrow;
            }

            if (nrow > 0)
            {
                check_petsc(MatSetValues(A, nrow, rows, 4, cols, values, ADD_VALUES),
                            "MatSetValues(element)");
            }
        }

        for (const Int ip : s.owned_nodes)
        {
            const PetscInt row = static_cast<PetscInt>(dof(ip));
            if (fixed(ip) != 0)
            {
                check_petsc(MatSetValue(A, row, row, 1.0, ADD_VALUES),
                            "MatSetValue(Dirichlet)");
            }
            else if (pseudotime_shift > 0.0)
            {
                const Real inv_ctimestep = 1.0 / s.elcoe2p(ip);
                if (!(inv_ctimestep > 0.0) || !std::isfinite(inv_ctimestep))
                {
                    throw std::runtime_error(
                        "SteadyEnergySolver: invalid pseudo-time thermal diagonal");
                }
                const Real mass_shift = pseudotime_shift * inv_ctimestep;
                if (!std::isfinite(mass_shift))
                {
                    throw std::runtime_error(
                        "SteadyEnergySolver: non-finite pseudo-time shift");
                }
                check_petsc(
                    MatSetValue(A, row, row, mass_shift, ADD_VALUES),
                    "MatSetValue(pseudo-time thermal diagonal)");
            }
        }

        check_petsc(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY), "MatAssemblyBegin");
        check_petsc(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY), "MatAssemblyEnd");

        // Experimental nonlinear AFC initialisation: obtain the stationary
        // monotone *low-order* solution once, rather than starting the AFC
        // defect iteration from a nearly frozen explicit pseudo-time state.
        // This is strictly opt-in and leaves production defaults unchanged.
        // Remove gamma*C/dt only from the temporary bootstrap matrix.
        const bool low_order_bootstrap =
            afc && !transient_be && (low_order_only ||
                env_flag("CBS3D_STEADY_ENERGY_LOW_ORDER_BOOTSTRAP"));
        Mat A_bootstrap = nullptr;
        if (low_order_bootstrap)
        {
            check_petsc(MatDuplicate(A, MAT_COPY_VALUES, &A_bootstrap),
                        "MatDuplicate(low-order stationary bootstrap)");
            if (pseudotime_shift > 0.0)
            {
                for (const Int ip : s.owned_nodes)
                {
                    if (fixed(ip) == 0)
                    {
                        const PetscInt row = static_cast<PetscInt>(dof(ip));
                        const Real shift = pseudotime_shift / s.elcoe2p(ip);
                        check_petsc(
                            MatSetValue(A_bootstrap, row, row, -shift, ADD_VALUES),
                            "MatSetValue(remove pseudo-time shift)");
                    }
                }
                check_petsc(MatAssemblyBegin(A_bootstrap, MAT_FINAL_ASSEMBLY),
                            "MatAssemblyBegin(bootstrap)");
                check_petsc(MatAssemblyEnd(A_bootstrap, MAT_FINAL_ASSEMBLY),
                            "MatAssemblyEnd(bootstrap)");
            }
        }

        // ================================
        // 4. KRYLOV SOLVER: GMRES + ASM/ILU(1) BY DEFAULT (ROBUST FOR THE
        //    NON-SYMMETRIC ADVECTION OPERATOR); PREFIX energy_ OVERRIDES,
        //    E.G. -energy_pc_type hypre. THE MATRIX IS FIXED, SO THE
        //    PRECONDITIONER IS BUILT ONCE FOR ALL AFC ITERATIONS.
        // ================================
        Vec x = nullptr;
        Vec b = nullptr;
        Vec r = nullptr;
        KSP ksp = nullptr;
        check_petsc(MatCreateVecs(A, &x, &b), "MatCreateVecs");
        check_petsc(VecDuplicate(b, &r), "VecDuplicate");

        check_petsc(KSPCreate(MPI_COMM_WORLD, &ksp), "KSPCreate");
        check_petsc(KSPSetOptionsPrefix(ksp, "energy_"), "KSPSetOptionsPrefix");
        check_petsc(KSPSetOperators(ksp, A, A), "KSPSetOperators");
        check_petsc(KSPSetType(ksp, KSPGMRES), "KSPSetType");
        check_petsc(KSPGMRESSetRestart(ksp, 200), "KSPGMRESSetRestart");
        check_petsc(KSPSetTolerances(ksp, 1.0e-12, 0.0, PETSC_DEFAULT, 5000),
                    "KSPSetTolerances");
        check_petsc(KSPSetInitialGuessNonzero(ksp, PETSC_TRUE), "KSPSetInitialGuessNonzero");
        {
            PC pc = nullptr;
            check_petsc(KSPGetPC(ksp, &pc), "KSPGetPC");
            check_petsc(PCSetType(pc, PCASM), "PCSetType(PCASM)");
            check_petsc(PCASMSetOverlap(pc, 1), "PCASMSetOverlap");
        }
        {
            // DEFAULT SUBDOMAIN SOLVER ILU(1), UNLESS GIVEN IN PETSC_OPTIONS
            PetscBool given = PETSC_FALSE;
            check_petsc(PetscOptionsHasName(nullptr, nullptr, "-energy_sub_pc_type", &given),
                        "PetscOptionsHasName");
            if (!given)
            {
                check_petsc(PetscOptionsSetValue(nullptr, "-energy_sub_pc_type", "ilu"),
                            "PetscOptionsSetValue");
                check_petsc(PetscOptionsSetValue(nullptr, "-energy_sub_pc_factor_levels", "1"),
                            "PetscOptionsSetValue");
            }
        }
        check_petsc(KSPSetFromOptions(ksp), "KSPSetFromOptions");

        // UNLIMITED GALERKIN RHS: c + q (FREE), T_D (DIRICHLET).
        // AFC INSTEAD USES THE SHARED EXPLICIT LIMITED RESIDUAL.
        Array1D<Real> base_rhs(s.cfg.npoin);
        base_rhs.fill(0.0);
        Real local_rhs_max = 0.0;
        for (const Int ip : s.owned_nodes)
        {
            if (fixed(ip) != 0)
            {
                base_rhs(ip) = value(ip);
            }
            else
            {
                base_rhs(ip) = constant_part(ip);
                local_rhs_max = std::max(local_rhs_max, std::abs(constant_part(ip)));
            }
        }
        const Real rhs_max = global_max(local_rhs_max);

        // INITIAL GUESS: THE CURRENT TEMPERATURE WITH THE DIRICHLET VALUES
        Array1D<Real> temperature = s.temperature;
        for (const Int ip : s.owned_nodes)
        {
            if (fixed(ip) != 0)
            {
                temperature(ip) = value(ip);
            }
        }
        HaloExchange::broadcastOwnedToGhosts(temperature, s.partition_metadata,
                                             MPI_COMM_WORLD);

        Array1D<Real> rhs_nodal(s.cfg.npoin);
        Array1D<Real> solved(s.cfg.npoin);

        Int total_krylov = 0;

        if (transient_be)
        {
            // P1 tetrahedral lumped physical thermal capacitance:
            // C_i = sum_{e contains i} (rho*cp)_e * Volume_e / 4.
            // In contrast, elcoe2p is the inverse *pseudo-time* diagonal
            // sum_e C_i,e / delte_e, and MUST NOT represent physical C/dt
            // when the CBS solver uses local timesteps.
            Array1D<Real> capacity(s.cfg.npoin);
            capacity.fill(0.0);
            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                const Real contribution =
                    s.rho_cp_e(ie) * s.detJ(ie) / 24.0;
                if (!(contribution > 0.0) ||
                    !std::isfinite(contribution))
                {
                    throw std::runtime_error(
                        "Implicit BE invalid physical thermal element capacity");
                }
                for (Int a = 1; a <= s.cfg.nep; ++a)
                {
                    capacity(s.intma(a, ie)) += contribution;
                }
            }
            HaloExchange::sumGhostContributionsToOwners(
                capacity, s.partition_metadata, MPI_COMM_WORLD);

            for (const Int ip : s.owned_nodes)
            {
                const Real ci = capacity(ip);
                if (!(ci > 0.0) || !std::isfinite(ci))
                {
                    throw std::runtime_error(
                        "Implicit BE invalid assembled physical nodal capacitance");
                }
                if (fixed(ip) == 0)
                {
                    const PetscInt row = static_cast<PetscInt>(dof(ip));
                    check_petsc(MatSetValue(A, row, row, ci / be_dt, ADD_VALUES),
                                "MatSetValue(physical C/dt)");
                }
            }
            check_petsc(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY),
                        "MatAssemblyBegin(physical BE)");
            check_petsc(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY),
                        "MatAssemblyEnd(physical BE)");
            check_petsc(KSPSetOperators(ksp, A, A),
                        "KSPSetOperators(physical BE)");

            Real max_linear_residual = 0.0;
            Real last_storage_rate = 0.0;
            for (Int step = 1; step <= be_steps; ++step)
            {
                Real local_storage_rate = 0.0;
                for (const Int ip : s.owned_nodes)
                {
                    rhs_nodal(ip) = fixed(ip) != 0 ? value(ip)
                        : base_rhs(ip) + capacity(ip) / be_dt * temperature(ip);
                }
                vec_from_nodes(b, s, dof, rhs_nodal);
                vec_from_nodes(x, s, dof, temperature);
                check_petsc(KSPSolve(ksp, b, x), "KSPSolve(physical BE)");
                KSPConvergedReason reason;
                PetscInt iterations = 0;
                check_petsc(KSPGetConvergedReason(ksp, &reason),
                            "KSPGetConvergedReason(physical BE)");
                check_petsc(KSPGetIterationNumber(ksp, &iterations),
                            "KSPGetIterationNumber(physical BE)");
                if (reason <= 0)
                {
                    throw std::runtime_error(
                        "Implicit physical BE PETSc KSP failed at step " +
                        std::to_string(step) + " reason " +
                        std::to_string(static_cast<int>(reason)));
                }
                total_krylov += static_cast<Int>(iterations);

                // Independently evaluate the full BE algebraic equation.
                check_petsc(MatMult(A, x, r), "MatMult(BE audit)");
                check_petsc(VecAYPX(r, -1.0, b), "VecAYPX(BE audit)");
                PetscReal rn = 0.0;
                PetscReal bn = 0.0;
                check_petsc(VecNorm(r, NORM_INFINITY, &rn),
                            "VecNorm(BE residual)");
                check_petsc(VecNorm(b, NORM_INFINITY, &bn),
                            "VecNorm(BE right-hand side)");
                const Real relative = bn > 0.0 ? rn / bn : rn;
                max_linear_residual = std::max(max_linear_residual, relative);
                if (!std::isfinite(relative) || relative >= residual_tol)
                {
                    throw std::runtime_error(
                        "Implicit BE algebraic residual failed at step " +
                        std::to_string(step));
                }

                nodes_from_vec(solved, s, dof, x);
                Real local_minimum = std::numeric_limits<Real>::max();
                Real local_maximum = -std::numeric_limits<Real>::max();
                for (const Int ip : s.owned_nodes)
                {
                    if (!std::isfinite(solved(ip)))
                    {
                        throw std::runtime_error(
                            "Implicit BE produced nonfinite temperature");
                    }
                    if (fixed(ip) == 0)
                    {
                        local_storage_rate += capacity(ip) *
                            (solved(ip) - temperature(ip)) / be_dt;
                    }
                    local_minimum = std::min(local_minimum, solved(ip));
                    local_maximum = std::max(local_maximum, solved(ip));
                }
                temperature = solved;
                result.t_min = -global_max(-local_minimum);
                result.t_max = global_max(local_maximum);
                last_storage_rate = global_sum(local_storage_rate);
                if (be_check_floor &&
                    result.t_min < be_floor -
                        1.0e-10 * std::max(1.0, std::abs(be_floor)))
                {
                    throw std::runtime_error(
                        "Implicit BE temperature fell below the requested "
                        "physical verification floor at step " +
                        std::to_string(step));
                }
                if (s.mpi_rank == 0 &&
                    (step <= 3 || step == be_steps || step % 10 == 0))
                {
                    std::cout
                        << "  [PHYSICAL-BE] step=" << step
                        << " time=" << be_dt * step
                        << " s Krylov=" << iterations
                        << " residual_inf_rel=" << relative
                        << " Tmin=" << result.t_min
                        << " Tmax=" << result.t_max
                        << " K storage_rate_free=" << last_storage_rate
                        << " W\n" << std::flush;
                }
            }
            // No post-hoc clipping of temperatures: the linear discrete
            // solution is retained exactly, including all edge flux balances.
            result.iterations = total_krylov;
            result.low_order = true;
            result.transient_be = true;
            result.be_steps = be_steps;
            result.be_dt = be_dt;
            result.residual_check = max_linear_residual;
            result.ksp_relative_residual = max_linear_residual;
            result.converged = std::isfinite(max_linear_residual) &&
                               max_linear_residual < residual_tol;
            s.temperature = temperature;
            s.temperature1 = temperature;
            result.seconds = MPI_Wtime() - start;
            check_petsc(KSPDestroy(&ksp), "KSPDestroy(physical BE)");
            check_petsc(VecDestroy(&x), "VecDestroy(x physical BE)");
            check_petsc(VecDestroy(&b), "VecDestroy(b physical BE)");
            check_petsc(VecDestroy(&r), "VecDestroy(r physical BE)");
            check_petsc(MatDestroy(&A), "MatDestroy(A physical BE)");
            return result;
        }

        if (low_order_bootstrap)
        {
            // Use the SAME stationary source vector and Dirichlet constraints
            // as the unlimited Galerkin solve, but on A_low = -(M+D).
            // The subsequent nonlinear AFC solve still targets the EXACT
            // shared, single-application production AFC residual.
            check_petsc(KSPSetOperators(ksp, A_bootstrap, A_bootstrap),
                        "KSPSetOperators(bootstrap)");
            vec_from_nodes(b, s, dof, base_rhs);
            vec_from_nodes(x, s, dof, temperature);
            check_petsc(KSPSolve(ksp, b, x), "KSPSolve(low-order bootstrap)");

            KSPConvergedReason bootstrap_reason;
            PetscInt bootstrap_iterations = 0;
            check_petsc(KSPGetConvergedReason(ksp, &bootstrap_reason),
                        "KSPGetConvergedReason(bootstrap)");
            check_petsc(KSPGetIterationNumber(ksp, &bootstrap_iterations),
                        "KSPGetIterationNumber(bootstrap)");
            if (bootstrap_reason <= 0)
            {
                throw std::runtime_error(
                    "Low-order steady AFC bootstrap KSP did not converge; reason "
                    + std::to_string(static_cast<int>(bootstrap_reason)));
            }
            total_krylov += static_cast<Int>(bootstrap_iterations);
            nodes_from_vec(temperature, s, dof, x);

            Real local_min = std::numeric_limits<Real>::max();
            Real local_max = -std::numeric_limits<Real>::max();
            for (const Int ip : s.owned_nodes)
            {
                if (!std::isfinite(temperature(ip)))
                {
                    throw std::runtime_error(
                        "Low-order steady AFC bootstrap returned nonfinite T");
                }
                local_min = std::min(local_min, temperature(ip));
                local_max = std::max(local_max, temperature(ip));
            }
            const Real minimum = -global_max(-local_min);
            const Real maximum = global_max(local_max);
            if (s.mpi_rank == 0)
            {
                std::cout << "  [AFC-BOOTSTRAP] low-order stationary KSP="
                          << bootstrap_iterations
                          << " Tmin=" << minimum
                          << " Tmax=" << maximum
                          << (low_order_only
                                 ? " K (stationary low-order result)"
                                 : " K (subsequent AFC residual still nonlinear)")
                          << "\n" << std::flush;
            }
            if (low_order_only)
            {
                // This solves a DIFFERENT discretisation from high-order AFC:
                // stationary graph-viscosity-stabilised P1 FEM. It is an
                // experimental bounded/robust baseline, NOT high-fidelity AFC.
                // Check the actual assembled low-order linear equation rather
                // than evaluating a different high-order nonlinear residual.
                check_petsc(MatMult(A_bootstrap, x, r),
                            "MatMult(low-order independent residual)");
                check_petsc(VecAYPX(r, -1.0, b),
                            "VecAYPX(b - A_low*T)");
                PetscReal residual_inf = 0.0;
                PetscReal rhs_inf = 0.0;
                check_petsc(VecNorm(r, NORM_INFINITY, &residual_inf),
                            "VecNorm(low-order residual infinity)");
                check_petsc(VecNorm(b, NORM_INFINITY, &rhs_inf),
                            "VecNorm(low-order RHS infinity)");
                result.residual_check = rhs_inf > 0.0
                    ? residual_inf / rhs_inf : residual_inf;
                result.ksp_relative_residual = result.residual_check;
                result.t_min = minimum;
                result.t_max = maximum;
                result.iterations = total_krylov;
                result.low_order = true;
                result.afc = false;
                result.converged =
                    std::isfinite(result.residual_check) &&
                    result.residual_check < residual_tol;
                s.temperature = temperature;
                s.temperature1 = temperature;
                result.seconds = MPI_Wtime() - start;
                if (s.mpi_rank == 0)
                {
                    std::cout << "  [LOW-ORDER] implicit stationary solve"
                              << " relative Linf residual="
                              << result.residual_check
                              << " accepted=" << result.converged
                              << "\n" << std::flush;
                }
                check_petsc(KSPDestroy(&ksp), "KSPDestroy(low-order)");
                check_petsc(VecDestroy(&x), "VecDestroy(x low-order)");
                check_petsc(VecDestroy(&b), "VecDestroy(b low-order)");
                check_petsc(VecDestroy(&r), "VecDestroy(r low-order)");
                check_petsc(MatDestroy(&A_bootstrap), "MatDestroy(A_low)");
                check_petsc(MatDestroy(&A), "MatDestroy(A)");
                return result;
            }

            // Restore the intended shifted defect preconditioner before
            // nonlinear iterations; PETSc rebuilds this preconditioner.
            check_petsc(KSPSetOperators(ksp, A, A),
                        "KSPSetOperators(restore defect)");
            check_petsc(MatDestroy(&A_bootstrap),
                        "MatDestroy(low-order bootstrap)");
        }

        Int afc_iterations = 0;
        Real last_change = 0.0;
        bool all_converged = true;
        const Int outer_max = afc ? afc_maxit : 1;

        std::vector<Array1D<Real>> anderson_x;
        std::vector<Array1D<Real>> anderson_f;
        bool linear_system_solved = false;
        Real accepted_rate = 0.0;
        Real accepted_merit = 0.0;
        Real accepted_relative_residual = 0.0;
        Real accepted_scale = 0.0;
        std::string accepted_method = "none";

        for (Int it = 1; it <= outer_max; ++it)
        {
            Real initial_rate = 0.0;
            Real initial_merit = 0.0;
            Real initial_relative_residual = 0.0;
            Real initial_raw_residual = 0.0;
            Array1D<Real> initial_residual(s.cfg.npoin);
            if (afc)
            {
                // Compute EXACTLY the same nonlinear residual and Zalesak
                // limiter as explicit Step 4, including its low-order predictor.
                // This changes the nonlinear target, not the PETSc matrix:
                // A_low * deltaT = R_explicit_AFC(T), with Dirichlet deltaT=0.
                s.temperature = temperature;
                explicit_residual(s, heat_flux_load, true);
                Real local_rate = 0.0;
                Real local_merit_squared = 0.0;
                Real local_raw = 0.0;
                for (const Int ip : s.owned_nodes)
                {
                    rhs_nodal(ip) = fixed(ip) != 0 ? 0.0 : s.rhs1(ip);
                    if (fixed(ip) == 0)
                    {
                        const Real scaled = s.rhs1(ip) * s.elcoe2p(ip);
                        local_rate = std::max(local_rate, std::abs(scaled));
                        local_merit_squared += scaled * scaled;
                        local_raw = std::max(local_raw, std::abs(s.rhs1(ip)));
                    }
                }
                initial_rate = global_max(local_rate);
                initial_merit = std::sqrt(global_sum(local_merit_squared));
                initial_raw_residual = global_max(local_raw);
                initial_residual = s.rhs1;

                // A tiny physical pseudo-time update AND a small raw defect
                // certify the loaded iterate without a potentially unstable
                // PETSc correction. Do not confuse a small update with an
                // acceptable nonlinear residual by itself.
                const Real normalised_raw =
                    rhs_max > 0.0 ? initial_raw_residual / rhs_max
                                  : initial_raw_residual;
                initial_relative_residual = normalised_raw;
                if (initial_rate < afc_tol &&
                    normalised_raw < residual_tol)
                {
                    last_change = 0.0;
                    afc_iterations = it - 1;
                    accepted_rate = initial_rate;
                    accepted_merit = initial_merit;
                    accepted_relative_residual = normalised_raw;
                    accepted_scale = 0.0;
                    accepted_method = "already converged";
                    if (s.mpi_rank == 0)
                    {
                        std::cout << "  [AFC] initial stationary residual pass"
                                  << "  update_est=" << initial_rate
                                  << " K, raw_rel=" << normalised_raw
                                  << "\n" << std::flush;
                    }
                    break;
                }

                vec_from_nodes(b, s, dof, rhs_nodal);
                check_petsc(VecSet(x, 0.0), "VecSet(delta T)");
            }
            else
            {
                // Original direct, un-limited Galerkin solve is unchanged.
                for (const Int ip : s.owned_nodes)
                {
                    rhs_nodal(ip) = base_rhs(ip);
                }
                vec_from_nodes(b, s, dof, rhs_nodal);
                vec_from_nodes(x, s, dof, temperature);
            }
            check_petsc(KSPSolve(ksp, b, x), "KSPSolve");
            linear_system_solved = true;

            KSPConvergedReason reason;
            PetscInt iterations = 0;
            check_petsc(KSPGetConvergedReason(ksp, &reason), "KSPGetConvergedReason");
            check_petsc(KSPGetIterationNumber(ksp, &iterations), "KSPGetIterationNumber");
            total_krylov += static_cast<Int>(iterations);
            if (reason <= 0)
            {
                throw std::runtime_error(
                    "Implicit energy PETSc KSP failed at AFC iteration " +
                    std::to_string(it) + ", reason " +
                    std::to_string(static_cast<int>(reason)));
            }
            all_converged = all_converged && reason > 0;

            nodes_from_vec(solved, s, dof, x);

            // ================================
            // OPTIONAL ANDERSON MIXING OF THE AFC DEFECT-CORRECTION MAP.
            //
            // F_K = G(T_K) - T_K = A_LOW^-1 R_EXPLICIT_AFC(T_K).
            // MINIMISE ||F_K - DELTA_F * GAMMA||_2 AND UPDATE
            // T_NEXT = T_K + BETA F_K - (DELTA_T + BETA DELTA_F) GAMMA.
            //
            // THE OLD UNDER-RELAXED UPDATE IS THE FALLBACK IF THE SMALL
            // LEAST-SQUARES SYSTEM IS DEGENERATE OR EXTRAPOLATION IS EXCESSIVE.
            // ================================
            Array1D<Real> fixed_point(s.cfg.npoin);
            fixed_point.fill(0.0);
            for (const Int ip : s.owned_nodes)
            {
                if (fixed(ip) == 0)
                {
                    fixed_point(ip) = afc ? solved(ip) : solved(ip) - temperature(ip);
                }
            }

            if (anderson_m > 0)
            {
                anderson_x.push_back(temperature);
                anderson_f.push_back(fixed_point);
                if (static_cast<Int>(anderson_x.size()) > anderson_m + 1)
                {
                    anderson_x.erase(anderson_x.begin());
                    anderson_f.erase(anderson_f.begin());
                }
            }

            Int n_mix = 0;
            std::vector<Real> gamma;
            if (anderson_m > 0 && anderson_x.size() > 1U)
            {
                n_mix = static_cast<Int>(anderson_x.size()) - 1;
                const Int count = n_mix * n_mix + n_mix;
                std::vector<Real> local(static_cast<std::size_t>(count), 0.0);
                std::vector<Real> global(static_cast<std::size_t>(count), 0.0);

                for (const Int ip : s.owned_nodes)
                {
                    if (fixed(ip) != 0)
                    {
                        continue;
                    }
                    Real df[5];
                    for (Int j = 0; j < n_mix; ++j)
                    {
                        df[j] = anderson_f[static_cast<std::size_t>(j + 1)](ip)
                              - anderson_f[static_cast<std::size_t>(j)](ip);
                        local[static_cast<std::size_t>(n_mix * n_mix + j)] +=
                            df[j] * fixed_point(ip);
                    }
                    for (Int i = 0; i < n_mix; ++i)
                    {
                        for (Int j = 0; j < n_mix; ++j)
                        {
                            local[static_cast<std::size_t>(i * n_mix + j)] +=
                                df[i] * df[j];
                        }
                    }
                }
                check_mpi(MPI_Allreduce(local.data(), global.data(), count,
                                        MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD),
                          "Anderson global least-squares products");

                // TIKHONOV REGULARISATION KEEPS DEPENDENT HISTORY STABLE.
                Real trace = 0.0;
                for (Int i = 0; i < n_mix; ++i)
                {
                    trace += global[static_cast<std::size_t>(i * n_mix + i)];
                }
                bool valid = std::isfinite(trace) && trace > 0.0;
                gamma.assign(static_cast<std::size_t>(n_mix), 0.0);
                if (valid)
                {
                    const Real shift = 1.0e-10 * trace / n_mix;
                    for (Int i = 0; i < n_mix; ++i)
                    {
                        global[static_cast<std::size_t>(i * n_mix + i)] += shift;
                    }
                    // SMALL DENSE GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING.
                    for (Int i = 0; i < n_mix && valid; ++i)
                    {
                        Int pivot = i;
                        for (Int j = i + 1; j < n_mix; ++j)
                        {
                            if (std::abs(global[static_cast<std::size_t>(j * n_mix + i)]) >
                                std::abs(global[static_cast<std::size_t>(pivot * n_mix + i)]))
                            {
                                pivot = j;
                            }
                        }
                        if (!(std::abs(global[static_cast<std::size_t>(pivot * n_mix + i)]) >
                              1.0e-14 * trace))
                        {
                            valid = false;
                            break;
                        }
                        if (pivot != i)
                        {
                            for (Int j = 0; j < n_mix; ++j)
                            {
                                std::swap(global[static_cast<std::size_t>(i * n_mix + j)],
                                          global[static_cast<std::size_t>(pivot * n_mix + j)]);
                            }
                            std::swap(global[static_cast<std::size_t>(n_mix * n_mix + i)],
                                      global[static_cast<std::size_t>(n_mix * n_mix + pivot)]);
                        }
                        for (Int j = i + 1; j < n_mix; ++j)
                        {
                            const Real scale =
                                global[static_cast<std::size_t>(j * n_mix + i)] /
                                global[static_cast<std::size_t>(i * n_mix + i)];
                            for (Int k = i; k < n_mix; ++k)
                            {
                                global[static_cast<std::size_t>(j * n_mix + k)] -=
                                    scale * global[static_cast<std::size_t>(i * n_mix + k)];
                            }
                            global[static_cast<std::size_t>(n_mix * n_mix + j)] -=
                                scale * global[static_cast<std::size_t>(n_mix * n_mix + i)];
                        }
                    }
                    for (Int i = n_mix - 1; i >= 0 && valid; --i)
                    {
                        Real rhs = global[static_cast<std::size_t>(n_mix * n_mix + i)];
                        for (Int j = i + 1; j < n_mix; ++j)
                        {
                            rhs -= global[static_cast<std::size_t>(i * n_mix + j)] *
                                   gamma[static_cast<std::size_t>(j)];
                        }
                        gamma[static_cast<std::size_t>(i)] =
                            rhs / global[static_cast<std::size_t>(i * n_mix + i)];
                        valid = valid && std::isfinite(gamma[static_cast<std::size_t>(i)]);
                    }
                }
                if (!valid)
                {
                    n_mix = 0;
                }
            }

            Array1D<Real> candidate = temperature;
            Real local_change = 0.0;
            Real local_baseline = 0.0;
            for (const Int ip : s.owned_nodes)
            {
                if (fixed(ip) != 0)
                {
                    continue;
                }
                const Real step = afc_relax * fixed_point(ip);
                local_baseline = std::max(local_baseline, std::abs(step));
                Real accelerated = step;
                for (Int j = 0; j < n_mix; ++j)
                {
                    const Real dx =
                        anderson_x[static_cast<std::size_t>(j + 1)](ip) -
                        anderson_x[static_cast<std::size_t>(j)](ip);
                    const Real df =
                        anderson_f[static_cast<std::size_t>(j + 1)](ip) -
                        anderson_f[static_cast<std::size_t>(j)](ip);
                    accelerated -= gamma[static_cast<std::size_t>(j)] *
                                   (dx + afc_relax * df);
                }
                candidate(ip) += accelerated;
                local_change = std::max(local_change, std::abs(accelerated));
            }

            Real local_values[3] = {local_change, local_baseline, 0.0};
            for (const Int ip : s.owned_nodes)
            {
                if (!std::isfinite(candidate(ip)))
                {
                    local_values[2] = 1.0;
                }
            }
            Real global_values[3] = {};
            check_mpi(MPI_Allreduce(local_values, global_values, 3, MPI_DOUBLE,
                                    MPI_MAX, MPI_COMM_WORLD),
                      "Anderson update safeguards");

            if (n_mix > 0 &&
                (global_values[2] != 0.0 ||
                 global_values[0] > 100.0 * std::max(global_values[1], 1.0e-12)))
            {
                // RESTART THE HISTORY AND USE THE ORIGINAL BOUNDED STEP.
                anderson_x.clear();
                anderson_f.clear();
                n_mix = 0;
            }

            if (afc)
            {
                // Globalize using the L2 norm of the EXACT production AFC
                // pseudo-time residual C_i^-1 R_i(T). Max-norm line searches
                // stall at limiter-switching nodes: a small step may reduce the
                // largest component without making useful global progress.
                //
                // The PETSc defect direction is only a preconditioned trial.
                // Evaluate it AND the stable explicit pseudo-time direction,
                // then retain the strongest residual reduction. Never accept
                // an infinitesimal PETSc step merely because it was tried first.
                const Real sufficient_decrease = armijo_c1;
                const Real max_rate_growth = 1.25;
                Array1D<Real> trial = temperature;
                Array1D<Real> best_trial = temperature;
                Real best_merit = initial_merit;
                Real best_rate = initial_rate;
                Real best_relative_residual = initial_relative_residual;
                Real best_scale = 0.0;
                std::string best_method = "none";
                bool step_accepted = false;

                auto search_direction = [&](const Array1D<Real>& direction,
                                            const char* direction_name,
                                            const Int search_depth)
                {
                    Real alpha = 1.0;
                    for (Int ls = 0; ls < search_depth; ++ls, alpha *= 0.5)
                    {
                        bool local_finite = true;
                        for (const Int ip : s.owned_nodes)
                        {
                            if (fixed(ip) == 0)
                            {
                                trial(ip) = temperature(ip) + alpha * direction(ip);
                                local_finite =
                                    local_finite && std::isfinite(trial(ip));
                            }
                        }

                        if (global_max(local_finite ? 0.0 : 1.0) != 0.0)
                        {
                            continue;
                        }

                        HaloExchange::broadcastOwnedToGhosts(
                            trial, s.partition_metadata, MPI_COMM_WORLD);
                        s.temperature = trial;
                        explicit_residual(s, heat_flux_load, true);

                        Real local_merit_squared = 0.0;
                        Real local_rate = 0.0;
                        Real local_raw = 0.0;
                        for (const Int ip : s.owned_nodes)
                        {
                            if (fixed(ip) == 0)
                            {
                                const Real scaled =
                                    s.rhs1(ip) * s.elcoe2p(ip);
                                local_merit_squared += scaled * scaled;
                                local_rate = std::max(
                                    local_rate, std::abs(scaled));
                                local_raw = std::max(
                                    local_raw, std::abs(s.rhs1(ip)));
                            }
                        }

                        const Real trial_merit =
                            std::sqrt(global_sum(local_merit_squared));
                        const Real trial_rate = global_max(local_rate);
                        const Real trial_raw = global_max(local_raw);
                        const Real trial_relative_residual =
                            rhs_max > 0.0 ? trial_raw / rhs_max : trial_raw;

                        // Audit-only: log the actual nonlinear trial residuals
                        // BEFORE imposing the L2 Armijo and max-rate guards.
                        // This distinguishes (1) an intrinsically non-descent
                        // direction from (2) an overly restrictive safeguard.
                        // The probe does not change the accepted iterates.
                        if (it == 1 && ls < 8 && s.mpi_rank == 0)
                        {
                            const bool finite =
                                std::isfinite(trial_merit) &&
                                std::isfinite(trial_rate) &&
                                std::isfinite(trial_relative_residual);
                            const bool max_rate_ok =
                                trial_rate <= max_rate_growth * initial_rate;
                            const bool merit_ok =
                                trial_merit <
                                    initial_merit *
                                        (1.0 - sufficient_decrease * alpha);
                            std::cout
                                << "  [AFC-PROBE] direction=" << direction_name
                                << " alpha=" << alpha
                                << " L2=" << trial_merit
                                << " L2_initial=" << initial_merit
                                << " max_rate=" << trial_rate
                                << " max_rate_initial=" << initial_rate
                                << " raw_rel=" << trial_relative_residual
                                << " finite=" << finite
                                << " Armijo=" << merit_ok
                                << " rate_guard=" << max_rate_ok
                                << "\n" << std::flush;
                        }

                        if (!std::isfinite(trial_merit) ||
                            !std::isfinite(trial_rate) ||
                            !std::isfinite(trial_relative_residual) ||
                            trial_rate > max_rate_growth * initial_rate)
                        {
                            continue;
                        }

                        if (trial_merit < best_merit &&
                            trial_merit <
                                initial_merit * (1.0 - sufficient_decrease * alpha))
                        {
                            best_trial = trial;
                            best_merit = trial_merit;
                            best_rate = trial_rate;
                            best_relative_residual = trial_relative_residual;
                            best_scale = alpha;
                            best_method = direction_name;
                            step_accepted = true;
                        }
                    }
                };

                Array1D<Real> direction = candidate;
                for (const Int ip : s.owned_nodes)
                {
                    if (fixed(ip) == 0)
                    {
                        direction(ip) -= temperature(ip);
                    }
                }
                search_direction(direction,
                                 n_mix > 0 ? "Anderson" : "PETSc defect", 14);

                if (n_mix > 0)
                {
                    direction.fill(0.0);
                    for (const Int ip : s.owned_nodes)
                    {
                        if (fixed(ip) == 0)
                        {
                            direction(ip) = afc_relax * fixed_point(ip);
                        }
                    }
                    search_direction(direction, "PETSc fallback", 14);
                }

                direction.fill(0.0);
                for (const Int ip : s.owned_nodes)
                {
                    if (fixed(ip) == 0)
                    {
                        direction(ip) =
                            initial_residual(ip) * s.elcoe2p(ip);
                    }
                }
                search_direction(direction, "explicit fallback", 8);

                if (!step_accepted)
                {
                    throw std::runtime_error(
                        "SteadyEnergySolver: neither PETSc nor explicit AFC "
                        "could reduce the L2 pseudo-time residual. "
                        "Initial merit=" + std::to_string(initial_merit)
                        + ", max update=" + std::to_string(initial_rate)
                        + ", relative raw residual="
                        + std::to_string(initial_relative_residual));
                }

                if (n_mix > 0 && best_method != "Anderson")
                {
                    anderson_x.clear();
                    anderson_f.clear();
                }

                local_change = 0.0;
                for (const Int ip : s.owned_nodes)
                {
                    if (fixed(ip) == 0)
                    {
                        local_change = std::max(
                            local_change,
                            std::abs(best_trial(ip) - temperature(ip)));
                    }
                }
                temperature = best_trial;
                last_change = global_max(local_change);
                accepted_rate = best_rate;
                accepted_merit = best_merit;
                accepted_relative_residual = best_relative_residual;
                accepted_scale = best_scale;
                accepted_method = best_method;
            }
            else
            {
                local_change = 0.0;
                for (const Int ip : s.owned_nodes)
                {
                    if (fixed(ip) != 0)
                    {
                        continue;
                    }
                    const Real updated = n_mix > 0 ? candidate(ip) :
                        temperature(ip) + afc_relax * fixed_point(ip);
                    local_change = std::max(
                        local_change, std::abs(updated - temperature(ip)));
                    temperature(ip) = updated;
                }
                HaloExchange::broadcastOwnedToGhosts(
                    temperature, s.partition_metadata, MPI_COMM_WORLD);
                last_change = global_max(local_change);
            }
            afc_iterations = it;

            if (afc && s.mpi_rank == 0 && (it <= 5 || it % 10 == 0))
            {
                std::cout << "  [AFC] iteration " << it << ": max |dT| = " << last_change
                          << " K, Krylov " << iterations;
                if (anderson_m > 0)
                {
                    std::cout << ", Anderson history " << n_mix;
                }
                std::cout
                    << ", rate=" << accepted_rate
                    << " K, L2=" << accepted_merit
                    << " K, raw_rel=" << accepted_relative_residual
                    << ", alpha=" << accepted_scale
                    << ", method=" << accepted_method
                    << "\n" << std::flush;
            }
            // A small step can result solely from aggressive backtracking.
            // Accept stationary convergence only when BOTH independent
            // residual checks pass. Do not mistake line-search stagnation
            // for convergence.
            if (afc &&
                accepted_rate < afc_tol &&
                accepted_relative_residual < residual_tol)
            {
                last_change = 0.0;
                break;
            }
        }

        // TRUE RELATIVE RESIDUAL OF THE LAST LINEAR SYSTEM.
        // An already-converged restart needs no PETSc solve.
        if (linear_system_solved)
        {
            // x holds the last PETSc solution: deltaT for AFC, T for Galerkin.
            check_petsc(MatMult(A, x, r), "MatMult");
            check_petsc(VecAYPX(r, -1.0, b), "VecAYPX");
            PetscReal rn = 0.0;
            PetscReal bn = 0.0;
            check_petsc(VecNorm(r, NORM_2, &rn), "VecNorm(r)");
            check_petsc(VecNorm(b, NORM_2, &bn), "VecNorm(b)");
            result.ksp_relative_residual = bn > 0.0 ? rn / bn : rn;
        }

        check_petsc(KSPDestroy(&ksp), "KSPDestroy");
        check_petsc(VecDestroy(&x), "VecDestroy(x)");
        check_petsc(VecDestroy(&b), "VecDestroy(b)");
        check_petsc(VecDestroy(&r), "VecDestroy(r)");
        check_petsc(MatDestroy(&A), "MatDestroy");

        result.iterations = total_krylov;
        result.afc_iterations = afc_iterations;
        result.afc_change = last_change;
        result.afc = afc;
        s.temperature = temperature;

        // Independent final nonlinear check: reevaluate the production
        // Step-4 residual at the final candidate, not at the previous iterate.
        // The explicit and steady solvers now share the AFC evaluator.
        explicit_residual(s, heat_flux_load, afc);

        Real local_res = 0.0;
        Real local_tmin = std::numeric_limits<Real>::max();
        Real local_tmax = -std::numeric_limits<Real>::max();
        for (const Int ip : s.owned_nodes)
        {
            if (fixed(ip) == 0)
            {
                local_res = std::max(local_res, std::abs(s.rhs1(ip)));
            }
            local_tmin = std::min(local_tmin, s.temperature(ip));
            local_tmax = std::max(local_tmax, s.temperature(ip));
        }
        check_mpi(MPI_Allreduce(&local_tmin, &result.t_min, 1, MPI_DOUBLE, MPI_MIN,
                                MPI_COMM_WORLD), "T min");
        check_mpi(MPI_Allreduce(&local_tmax, &result.t_max, 1, MPI_DOUBLE, MPI_MAX,
                                MPI_COMM_WORLD), "T max");
        const Real global_res = global_max(local_res);
        result.residual_check = rhs_max > 0.0 ? global_res / rhs_max : global_res;
        result.converged =
            all_converged &&
            (!afc || last_change < afc_tol) &&
            std::isfinite(result.residual_check) &&
            result.residual_check < residual_tol;

        s.temperature1 = s.temperature;
        result.seconds = MPI_Wtime() - start;
        return result;
    }
}

#else

namespace cbs
{
    SteadyEnergySolver::Result SteadyEnergySolver::solve(
        CBSStateSI&, const Array1D<Int>&, const Array1D<Real>&, const Array1D<Real>&)
    {
        throw std::runtime_error("SteadyEnergySolver needs MPI and PETSc");
    }
}

#endif
