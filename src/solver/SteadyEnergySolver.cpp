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
        void explicit_residual(CBSStateSI& s, const Array1D<Real>& heat_flux_load)
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

            HaloExchange::sumGhostContributionsToOwners(
                s.rhs1, s.partition_metadata, MPI_COMM_WORLD);

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

        // ALGEBRAIC FLUX CORRECTION (SAME SWITCH AS THE EXPLICIT SOLVER):
        // LOW-ORDER MATRIX + LIMITED ANTIDIFFUSION BY DEFECT CORRECTION
        const bool afc = env_flag("CBS3D_TEMPERATURE_AFC");
        const Real afc_tol = env_real("CBS3D_STEADY_ENERGY_TOL", 1.0e-6);
        const Real residual_tol = env_real("CBS3D_STEADY_ENERGY_RESIDUAL_TOL", 1.0e-6);
        if (!(residual_tol > 0.0) || !std::isfinite(residual_tol))
        {
            throw std::runtime_error("CBS3D_STEADY_ENERGY_RESIDUAL_TOL must be positive and finite");
        }
        const Int afc_maxit = static_cast<Int>(env_real("CBS3D_STEADY_ENERGY_MAXIT", 300.0));
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
        //    q_i = SUM_b max(0, l_ib) IS KEPT FOR THE LIMITER.
        // ================================
        std::vector<Real> graph(afc ? static_cast<std::size_t>(s.cfg.nelem) * 6U : 0U, 0.0);
        Array1D<Real> q_low(s.cfg.npoin);
        q_low.fill(0.0);

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
                    graph[static_cast<std::size_t>(ie - 1) * 6U + static_cast<std::size_t>(k)] = d;
                    M[a][b] += d;
                    M[b][a] += d;
                    M[a][a] -= d;
                    M[b][b] -= d;
                }
            }

            for (Int a = 1; a <= 4; ++a)
            {
                for (Int b = 1; b <= 4; ++b)
                {
                    if (a != b)
                    {
                        q_low(s.intma(a, ie)) += std::max(0.0, M[a][b]);
                    }
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

        HaloExchange::sumGhostContributionsToOwners(q_low, s.partition_metadata,
                                                    MPI_COMM_WORLD);
        HaloExchange::broadcastOwnedToGhosts(q_low, s.partition_metadata,
                                             MPI_COMM_WORLD);

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

        // RIGHT-HAND SIDE WITHOUT ANTIDIFFUSION: c + q (FREE), T_D (DIRICHLET)
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
        Array1D<Real> t_lo(s.cfg.npoin);
        Array1D<Real> t_hi(s.cfg.npoin);
        Array1D<Real> p_plus(s.cfg.npoin);
        Array1D<Real> p_minus(s.cfg.npoin);
        Array1D<Real> r_plus(s.cfg.npoin);
        Array1D<Real> r_minus(s.cfg.npoin);
        Array1D<Real> correction(s.cfg.npoin);

        Int total_krylov = 0;
        Int afc_iterations = 0;
        Real last_change = 0.0;
        bool all_converged = true;
        const Int outer_max = afc ? afc_maxit : 1;

        // USE THE SAME LIMITER FOR EACH DEFECT CORRECTION AND FOR THE
        // INDEPENDENT NONLINEAR RESIDUAL CHECK AT THE FINAL TEMPERATURE.
        auto assemble_afc_correction = [&]()
        {
            correction.fill(0.0);
            if (afc)
            {
                t_lo = temperature;
                t_hi = temperature;
                for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
                {
                    Real emin = std::numeric_limits<Real>::max();
                    Real emax = -std::numeric_limits<Real>::max();
                    for (Int a = 1; a <= 4; ++a)
                    {
                        const Real t = temperature(s.intma(a, ie));
                        emin = std::min(emin, t);
                        emax = std::max(emax, t);
                    }
                    for (Int a = 1; a <= 4; ++a)
                    {
                        const Int ip = s.intma(a, ie);
                        t_lo(ip) = std::min(t_lo(ip), emin);
                        t_hi(ip) = std::max(t_hi(ip), emax);
                    }
                }
                // MAXIMUM AS MINIMUM OF -T
                for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
                {
                    t_hi(ip) = -t_hi(ip);
                }
                HaloExchange::minGhostContributionsToOwners(t_lo, s.partition_metadata,
                                                            MPI_COMM_WORLD);
                HaloExchange::minGhostContributionsToOwners(t_hi, s.partition_metadata,
                                                            MPI_COMM_WORLD);
                HaloExchange::broadcastOwnedToGhosts(t_lo, s.partition_metadata,
                                                     MPI_COMM_WORLD);
                HaloExchange::broadcastOwnedToGhosts(t_hi, s.partition_metadata,
                                                     MPI_COMM_WORLD);
                for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
                {
                    t_hi(ip) = -t_hi(ip);
                }

                p_plus.fill(0.0);
                p_minus.fill(0.0);
                for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
                {
                    if (s.mat_elem(ie) != 0)
                    {
                        continue;
                    }
                    for (Int k = 0; k < 6; ++k)
                    {
                        const Real d = graph[static_cast<std::size_t>(ie - 1) * 6U + static_cast<std::size_t>(k)];
                        if (d <= 0.0)
                        {
                            continue;
                        }
                        const Int i = s.intma(pair_a[k], ie);
                        const Int j = s.intma(pair_b[k], ie);
                        const Real f = d * (temperature(i) - temperature(j));
                        if (f >= 0.0)
                        {
                            p_plus(i) += f;
                            p_minus(j) -= f;
                        }
                        else
                        {
                            p_minus(i) += f;
                            p_plus(j) -= f;
                        }
                    }
                }
                HaloExchange::sumGhostContributionsToOwners(p_plus, s.partition_metadata,
                                                            MPI_COMM_WORLD);
                HaloExchange::sumGhostContributionsToOwners(p_minus, s.partition_metadata,
                                                            MPI_COMM_WORLD);

                for (const Int ip : s.owned_nodes)
                {
                    const Real qp = q_low(ip) * (t_hi(ip) - temperature(ip));
                    const Real qm = q_low(ip) * (t_lo(ip) - temperature(ip));
                    r_plus(ip) = p_plus(ip) > 0.0 ? std::min(1.0, qp / p_plus(ip)) : 1.0;
                    r_minus(ip) = p_minus(ip) < 0.0 ? std::min(1.0, qm / p_minus(ip)) : 1.0;
                    if (fixed(ip) != 0)
                    {
                        r_plus(ip) = 1.0;
                        r_minus(ip) = 1.0;
                    }
                }
                HaloExchange::broadcastOwnedToGhosts(r_plus, s.partition_metadata,
                                                     MPI_COMM_WORLD);
                HaloExchange::broadcastOwnedToGhosts(r_minus, s.partition_metadata,
                                                     MPI_COMM_WORLD);

                for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
                {
                    if (s.mat_elem(ie) != 0)
                    {
                        continue;
                    }
                    for (Int k = 0; k < 6; ++k)
                    {
                        const Real d = graph[static_cast<std::size_t>(ie - 1) * 6U + static_cast<std::size_t>(k)];
                        if (d <= 0.0)
                        {
                            continue;
                        }
                        const Int i = s.intma(pair_a[k], ie);
                        const Int j = s.intma(pair_b[k], ie);
                        const Real f = d * (temperature(i) - temperature(j));
                        const Real alpha = f >= 0.0
                            ? std::min(r_plus(i), r_minus(j))
                            : std::min(r_minus(i), r_plus(j));
                        correction(i) += alpha * f;
                        correction(j) -= alpha * f;
                    }
                }
                HaloExchange::sumGhostContributionsToOwners(correction, s.partition_metadata,
                                                            MPI_COMM_WORLD);
            }

        };

        for (Int it = 1; it <= outer_max; ++it)
        {
            // ================================
            // 5. LIMITED ANTIDIFFUSIVE FLUXES FROM THE CURRENT ITERATE
            //    f_ab = d_ab (T_a - T_b) INTO a; ZALESAK LIMITER WITH THE
            //    STEADY BOUNDS Q_i = q_i (T_i^max/min - T_i)
            // ================================
            assemble_afc_correction();

            // ================================
            // 6. SOLVE  A T* = b + g
            // ================================
            for (const Int ip : s.owned_nodes)
            {
                rhs_nodal(ip) = base_rhs(ip) + (fixed(ip) != 0 ? 0.0 : correction(ip));
            }
            vec_from_nodes(b, s, dof, rhs_nodal);
            vec_from_nodes(x, s, dof, temperature);
            check_petsc(KSPSolve(ksp, b, x), "KSPSolve");

            KSPConvergedReason reason;
            PetscInt iterations = 0;
            check_petsc(KSPGetConvergedReason(ksp, &reason), "KSPGetConvergedReason");
            check_petsc(KSPGetIterationNumber(ksp, &iterations), "KSPGetIterationNumber");
            total_krylov += static_cast<Int>(iterations);
            all_converged = all_converged && reason > 0;

            nodes_from_vec(solved, s, dof, x);

            Real local_change = 0.0;
            for (const Int ip : s.owned_nodes)
            {
                const Real updated = temperature(ip) + afc_relax * (solved(ip) - temperature(ip));
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
                          << " K, Krylov " << iterations << "\n" << std::flush;
            }
            if (afc && last_change < afc_tol)
            {
                break;
            }
        }

        // TRUE RELATIVE RESIDUAL OF THE LAST LINEAR SYSTEM
        {
            vec_from_nodes(x, s, dof, temperature);
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

        // ================================
        // 7. INDEPENDENT CHECK AT THE FINAL TEMPERATURE.
        //    HIGH-ORDER RESIDUAL: R_HIGH(T).
        //    FULL GRAPH FLUX:     F_i = SUM_j d_ij (T_i - T_j).
        //    LIMITED GRAPH FLUX:  G_i = SUM_j alpha_ij F_ij.
        //    LOW-ORDER AFC EQUATION: R_HIGH(T) - F(T) + G(T) = 0.
        //    THE FINAL LIMITER MUST BE RECOMPUTED: THE LAST LINEAR RHS
        //    USED THE PREVIOUS (POSSIBLY UNDER-RELAXED) ITERATE.
        // ================================
        Array1D<Real> full_antidiffusion(s.cfg.npoin);
        full_antidiffusion.fill(0.0);
        if (afc)
        {
            assemble_afc_correction();
            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                if (s.mat_elem(ie) != 0)
                {
                    continue;
                }
                for (Int k = 0; k < 6; ++k)
                {
                    const Real d = graph[static_cast<std::size_t>(ie - 1) * 6U +
                                         static_cast<std::size_t>(k)];
                    if (d <= 0.0)
                    {
                        continue;
                    }
                    const Int i = s.intma(pair_a[k], ie);
                    const Int j = s.intma(pair_b[k], ie);
                    const Real f = d * (temperature(i) - temperature(j));
                    full_antidiffusion(i) += f;
                    full_antidiffusion(j) -= f;
                }
            }
            HaloExchange::sumGhostContributionsToOwners(
                full_antidiffusion, s.partition_metadata, MPI_COMM_WORLD);
        }

        explicit_residual(s, heat_flux_load);

        Real local_res = 0.0;
        Real local_tmin = std::numeric_limits<Real>::max();
        Real local_tmax = -std::numeric_limits<Real>::max();
        for (const Int ip : s.owned_nodes)
        {
            if (fixed(ip) == 0)
            {
                Real residual = s.rhs1(ip);
                if (afc)
                {
                    residual -= full_antidiffusion(ip) - correction(ip);
                }
                local_res = std::max(local_res, std::abs(residual));
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
