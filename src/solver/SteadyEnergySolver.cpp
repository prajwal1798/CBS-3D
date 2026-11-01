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
                // TemperatureAFC performs its own reverse-add, flux limiting and
                // owner residual assembly. Add the face load once beforehand.
                for (const Int ip : s.owned_nodes)
                {
                    s.rhs1(ip) += heat_flux_load(ip);
                }
                TemperatureAFC::limitDistributedResidual(s);
            }
            else
            {
                HaloExchange::sumGhostContributionsToOwners(
                    s.rhs1, s.partition_metadata, MPI_COMM_WORLD);
                for (const Int ip : s.owned_nodes)
                {
                    s.rhs1(ip) += heat_flux_load(ip);
                }
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
        // 3. MATRIX  A = -(M + D)  (FREE ROWS), IDENTITY (DIRICHLET ROWS).
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
            if (fixed(ip) != 0)
            {
                const PetscInt row = static_cast<PetscInt>(dof(ip));
                check_petsc(MatSetValue(A, row, row, 1.0, ADD_VALUES),
                            "MatSetValue(Dirichlet)");
            }
        }

        check_petsc(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY), "MatAssemblyBegin");
        check_petsc(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY), "MatAssemblyEnd");

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
        Int afc_iterations = 0;
        Real last_change = 0.0;
        bool all_converged = true;
        const Int outer_max = afc ? afc_maxit : 1;

        std::vector<Array1D<Real>> anderson_x;
        std::vector<Array1D<Real>> anderson_f;

        for (Int it = 1; it <= outer_max; ++it)
        {
            if (afc)
            {
                // Compute EXACTLY the same nonlinear residual and Zalesak
                // limiter as explicit Step 4, including its low-order predictor.
                // This changes the nonlinear target, not the PETSc matrix:
                // A_low * deltaT = R_explicit_AFC(T), with Dirichlet deltaT=0.
                s.temperature = temperature;
                explicit_residual(s, heat_flux_load, true);
                for (const Int ip : s.owned_nodes)
                {
                    rhs_nodal(ip) = fixed(ip) != 0 ? 0.0 : s.rhs1(ip);
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

            local_change = 0.0;
            for (const Int ip : s.owned_nodes)
            {
                if (fixed(ip) != 0)
                {
                    continue;
                }
                const Real updated = n_mix > 0 ? candidate(ip) :
                    temperature(ip) + afc_relax * fixed_point(ip);
                local_change = std::max(local_change, std::abs(updated - temperature(ip)));
                temperature(ip) = updated;
            }
            HaloExchange::broadcastOwnedToGhosts(temperature, s.partition_metadata,
                                                 MPI_COMM_WORLD);
            last_change = global_max(local_change);
            afc_iterations = it;

            if (afc && s.mpi_rank == 0 && (it <= 5 || it % 10 == 0))
            {
                std::cout << "  [AFC] iteration " << it << ": max |dT| = " << last_change
                          << " K, Krylov " << iterations;
                if (anderson_m > 0)
                {
                    std::cout << ", Anderson history " << n_mix;
                }
                std::cout << "\n" << std::flush;
            }
            if (afc && last_change < afc_tol)
            {
                break;
            }
        }

        // TRUE RELATIVE RESIDUAL OF THE LAST LINEAR SYSTEM
        {
            // x holds the last PETSc solution: deltaT for AFC, T for Galerkin.
            // Do not overwrite deltaT with T before checking the linear residual.
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
