//=============================================================================
// CBS3D++_SI
//
// Monolithic CBS flow solve for steady fixed points and physical backward Euler.
//
// Unknown ordering per active fluid node:
//
//   [u, v, w, u*, v*, w*, p].
//
// The solver retains the production P1/TET4 CBS Step-1/2/3 discretisation.
// Picard iteration freezes the advecting velocity in the nonlinear Step-1
// fluxes. In physical-time mode the old accepted velocity u^n appears only in
// the BE mass RHS and in the existing Step-2 boundary flux convention:
//
//   A_u(w^k) u^(n+1) + M/dt u*^(n+1) = M/dt u^n
//   dt H p^(n+1) - B_vol u*^(n+1)     = -B_bnd u^n
//   M/dt (u^(n+1)-u*^(n+1)) + Gp      = 0.
//
// Temperature is not part of this PETSc block because the present blanket case
// has constant properties and no T -> flow feedback. The transient driver
// advances implicit BE energy after the converged flow step at the same dt.
//=============================================================================

#include "cbs/solver/FullyImplicitFlowSolver.hpp"

#include <stdexcept>

#if defined(CBS3D_USE_MPI) && defined(CBS3D_USE_PETSC)

#include "cbs/boundary/Boundary.hpp"
#include "cbs/linalg/PetscHeaders.hpp"
#include "cbs/parallel/HaloExchange.hpp"
#include "cbs/timestep/TimeStep.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace cbs
{
    namespace
    {
        constexpr Int corrected_offset = 0;
        constexpr Int predictor_offset = 3;
        constexpr Int pressure_offset = 6;
        constexpr Int dofs_per_active_node = 7;

        void check_mpi(const int code, const char* what)
        {
            if (code != MPI_SUCCESS)
            {
                throw std::runtime_error(std::string("MPI failure in ") + what);
            }
        }

        void check_petsc(const PetscErrorCode code, const char* what)
        {
            if (code != PETSC_SUCCESS)
            {
                // Return both the numeric code and PETSc's message. In
                // particular, PETSC_ERR_MPI (98) is NOT a matrix preallocation
                // error, and requires checking MPI/PETSc runtime compatibility.
                const char* error_text = nullptr;
                const char* detail = nullptr;
                (void)PetscErrorMessage(code, &error_text, &detail);
                std::ostringstream message;
                message << "PETSc failure in " << what
                        << " (code " << static_cast<int>(code) << ")";
                if (error_text != nullptr)
                {
                    message << ": " << error_text;
                }
                if (detail != nullptr)
                {
                    message << " [" << detail << "]";
                }
                throw std::runtime_error(message.str());
            }
        }

        Real env_real(const char* name, const Real fallback)
        {
            const char* text = std::getenv(name);
            if (text == nullptr || *text == '\0')
            {
                return fallback;
            }

            char* end = nullptr;
            const double value = std::strtod(text, &end);
            if (end == text || *end != '\0' || !std::isfinite(value))
            {
                throw std::runtime_error(
                    std::string("Invalid real value for ") + name);
            }
            return static_cast<Real>(value);
        }

        Int env_int(const char* name, const Int fallback)
        {
            const Real value = env_real(name, static_cast<Real>(fallback));
            if (value < 1.0 || std::floor(value) != value)
            {
                throw std::runtime_error(
                    std::string("Invalid positive integer value for ") + name);
            }
            return static_cast<Int>(value);
        }

        bool touches_fluid(const CBSStateSI& s, const Int ip)
        {
            return
                (s.node_material_mask(ip) & CBSStateSI::node_touches_fluid) != 0;
        }

        bool velocity_fixed(const CBSStateSI& s, const Int ip)
        {
            return s.node_velocity_bc_type(ip) != CBSStateSI::velocity_bc_free;
        }

        Int derivative_index(
            const CBSStateSI& s,
            const Int ie,
            const Int dim,
            const Int local_node)
        {
            return (ie - 1) * s.cfg.ndim * s.cfg.nep
                + (dim - 1) * s.cfg.nep
                + local_node;
        }

        Real grad(
            const CBSStateSI& s,
            const Int ie,
            const Int dim,
            const Int local_node)
        {
            return s.dNkdx(derivative_index(s, ie, dim, local_node));
        }

        Real fluid_kinematic_viscosity(const CBSStateSI& s, const Int ie)
        {
            if (s.cfg.dimensional_mode > 0 &&
                s.cfg.material_properties_enabled > 0)
            {
                if (!(s.rho_e(ie) > 0.0) || !std::isfinite(s.rho_e(ie)))
                {
                    throw std::runtime_error(
                        "Fully implicit CBS: invalid fluid density");
                }

                const Real mu =
                    s.cfg.turbulence_on > 0 ? s.mu_eff_e(ie) : s.mu_e(ie);

                if (!(mu >= 0.0) || !std::isfinite(mu))
                {
                    throw std::runtime_error(
                        "Fully implicit CBS: invalid fluid viscosity");
                }
                return mu / s.rho_e(ie);
            }

            if (!(s.cfg.ani >= 0.0) || !std::isfinite(s.cfg.ani))
            {
                throw std::runtime_error(
                    "Fully implicit CBS: invalid nondimensional viscosity");
            }
            return s.cfg.ani;
        }

        void set_default_option(const char* name, const char* value)
        {
            PetscBool given = PETSC_FALSE;
            check_petsc(
                PetscOptionsHasName(nullptr, nullptr, name, &given),
                "PetscOptionsHasName");
            if (!given)
            {
                check_petsc(
                    PetscOptionsSetValue(nullptr, name, value),
                    "PetscOptionsSetValue");
            }
        }

        void install_default_petsc_options()
        {
            set_default_option("-flow_ksp_type", "fgmres");
            set_default_option("-flow_ksp_rtol", "1e-8");
            set_default_option("-flow_ksp_max_it", "500");
            set_default_option("-flow_pc_type", "fieldsplit");
            set_default_option("-flow_pc_fieldsplit_type", "schur");
            set_default_option(
                "-flow_pc_fieldsplit_schur_fact_type", "full");
            set_default_option(
                "-flow_pc_fieldsplit_schur_precondition", "a11");

            // The six-component velocity/predictor block is nonsymmetric.
            set_default_option(
                "-flow_fieldsplit_velocity_ksp_type", "gmres");
            set_default_option(
                "-flow_fieldsplit_velocity_ksp_rtol", "1e-2");
            set_default_option(
                "-flow_fieldsplit_velocity_ksp_max_it", "40");
            set_default_option(
                "-flow_fieldsplit_velocity_pc_type", "asm");
            set_default_option(
                "-flow_fieldsplit_velocity_pc_asm_overlap", "1");
            set_default_option(
                "-flow_fieldsplit_velocity_sub_pc_type", "ilu");
            set_default_option(
                "-flow_fieldsplit_velocity_sub_pc_factor_levels", "1");

            // The pressure preconditioning block is tau times a Laplacian.
            set_default_option(
                "-flow_fieldsplit_pressure_ksp_type", "gmres");
            set_default_option(
                "-flow_fieldsplit_pressure_ksp_rtol", "1e-2");
            set_default_option(
                "-flow_fieldsplit_pressure_ksp_max_it", "40");
            set_default_option(
                "-flow_fieldsplit_pressure_pc_type", "hypre");
            set_default_option(
                "-flow_fieldsplit_pressure_pc_hypre_type", "boomeramg");
        }

        struct Metrics
        {
            Real rel_u = 0.0;
            Real rel_p = 0.0;
            Real div_rms = 0.0;
            Real umax = 0.0;
            Real correction_max = 0.0;
        };

        Metrics update_state_from_solution(
            CBSStateSI& s,
            Vec x,
            const Real relaxation)
        {
            PetscInt local_start = 0;
            PetscInt local_end = 0;
            check_petsc(
                VecGetOwnershipRange(x, &local_start, &local_end),
                "VecGetOwnershipRange");

            const PetscScalar* values = nullptr;
            check_petsc(VecGetArrayRead(x, &values), "VecGetArrayRead");

            long long cursor = 0;
            Real local_du2 = 0.0;
            Real local_u2 = 0.0;
            Real local_dp2 = 0.0;
            Real local_p2 = 0.0;
            Real local_umax = 0.0;
            Real local_correction_max = 0.0;

            for (const Int ip : s.owned_nodes)
            {
                if (!touches_fluid(s, ip))
                {
                    continue;
                }

                for (Int dim = 1; dim <= 3; ++dim)
                {
                    const Real old_value = s.unkno(dim, ip);
                    const Real linear_value =
                        static_cast<Real>(
                            values[cursor + corrected_offset + (dim - 1)]);
                    const Real new_value =
                        velocity_fixed(s, ip)
                            ? linear_value
                            : old_value +
                                relaxation * (linear_value - old_value);

                    const Real delta = new_value - old_value;
                    if (!velocity_fixed(s, ip))
                    {
                        local_du2 += delta * delta;
                        local_u2 += new_value * new_value;
                    }
                    s.unkno(dim, ip) = new_value;

                    const Real predictor_value =
                        static_cast<Real>(
                            values[cursor + predictor_offset + (dim - 1)]);
                    local_correction_max =
                        std::max(
                            local_correction_max,
                            std::abs(linear_value - predictor_value));
                }

                const Real old_pressure = s.pres(ip);
                const Real linear_pressure =
                    static_cast<Real>(values[cursor + pressure_offset]);
                const Real new_pressure =
                    s.node_pressure_fixed(ip) != 0
                        ? linear_pressure
                        : old_pressure +
                            relaxation *
                                (linear_pressure - old_pressure);

                const Real dp = new_pressure - old_pressure;
                if (s.node_pressure_fixed(ip) == 0)
                {
                    local_dp2 += dp * dp;
                    local_p2 += new_pressure * new_pressure;
                }
                s.pres(ip) = new_pressure;

                const Real ux = s.unkno(1, ip);
                const Real uy = s.unkno(2, ip);
                const Real uz = s.unkno(3, ip);
                local_umax = std::max(
                    local_umax,
                    std::sqrt(ux * ux + uy * uy + uz * uz));

                cursor += dofs_per_active_node;
            }

            check_petsc(VecRestoreArrayRead(x, &values), "VecRestoreArrayRead");

            if (local_start + static_cast<PetscInt>(cursor) != local_end)
            {
                throw std::runtime_error(
                    "Fully implicit CBS: local vector/DOF ownership mismatch");
            }

            Boundary::applyOwnedVelocityConstraints(s);

            HaloExchange::broadcastOwnedToGhosts(
                s.unkno, s.partition_metadata, MPI_COMM_WORLD);
            HaloExchange::broadcastOwnedToGhosts(
                s.pres, s.partition_metadata, MPI_COMM_WORLD);

            Real local_div2 = 0.0;
            Real local_volume = 0.0;
            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                if (s.mat_elem(ie) != 0)
                {
                    continue;
                }

                Real div_u = 0.0;
                for (Int a = 1; a <= 4; ++a)
                {
                    const Int ip = s.intma(a, ie);
                    for (Int dim = 1; dim <= 3; ++dim)
                    {
                        div_u +=
                            s.unkno(dim, ip) * grad(s, ie, dim, a);
                    }
                }

                const Real volume = s.detJ(ie) / 6.0;
                local_div2 += volume * div_u * div_u;
                local_volume += volume;
            }

            Real global_du2 = 0.0;
            Real global_u2 = 0.0;
            Real global_dp2 = 0.0;
            Real global_p2 = 0.0;
            Real global_div2 = 0.0;
            Real global_volume = 0.0;
            Real global_umax = 0.0;
            Real global_correction_max = 0.0;

            check_mpi(MPI_Allreduce(
                &local_du2, &global_du2, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD), "flow du norm");
            check_mpi(MPI_Allreduce(
                &local_u2, &global_u2, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD), "flow u norm");
            check_mpi(MPI_Allreduce(
                &local_dp2, &global_dp2, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD), "flow dp norm");
            check_mpi(MPI_Allreduce(
                &local_p2, &global_p2, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD), "flow p norm");
            check_mpi(MPI_Allreduce(
                &local_div2, &global_div2, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD), "flow divergence norm");
            check_mpi(MPI_Allreduce(
                &local_volume, &global_volume, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD), "flow fluid volume");
            check_mpi(MPI_Allreduce(
                &local_umax, &global_umax, 1, MPI_DOUBLE, MPI_MAX,
                MPI_COMM_WORLD), "flow Umax");
            check_mpi(MPI_Allreduce(
                &local_correction_max, &global_correction_max,
                1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD),
                "flow velocity correction max");

            Metrics m;
            m.rel_u =
                std::sqrt(global_du2) /
                std::max(
                    std::sqrt(global_u2),
                    static_cast<Real>(1.0e-30));
            m.rel_p =
                std::sqrt(global_dp2) /
                std::max(
                    std::sqrt(global_p2),
                    static_cast<Real>(1.0));
            m.div_rms =
                global_volume > 0.0
                    ? std::sqrt(global_div2 / global_volume)
                    : 0.0;
            m.umax = global_umax;
            m.correction_max = global_correction_max;
            return m;
        }
    }

    FullyImplicitFlowSolver::Result
    FullyImplicitFlowSolver::solveImpl(
        CBSStateSI& s,
        const bool backward_euler,
        const Real physical_dt)
    {
        if (!s.mpi_enabled || s.mpi_size < 2)
        {
            throw std::runtime_error(
                "Fully implicit CBS requires distributed MPI execution");
        }
        if (s.cfg.ndim != 3 || s.cfg.nep != 4)
        {
            throw std::runtime_error(
                "Fully implicit CBS currently requires 3D P1/TET4");
        }
        if (s.cfg.turbulence_on > 0)
        {
            throw std::runtime_error(
                "Fully implicit CBS baseline is laminar only");
        }

        for (const Int ip : s.owned_nodes)
        {
            if (s.node_symmetry(ip) != 0)
            {
                throw std::runtime_error(
                    "Fully implicit CBS does not yet support distributed "
                    "symmetry constraints");
            }
        }

        const double start = MPI_Wtime();
        Result result;

        if (backward_euler &&
            (!(physical_dt > 0.0) || !std::isfinite(physical_dt)))
        {
            throw std::runtime_error(
                "Fully implicit CBS backward Euler needs a positive finite dt");
        }

        const Real nonlinear_tol =
            env_real("CBS3D_IMPLICIT_FLOW_TOL", 1.0e-7);
        const Real pressure_tol =
            env_real("CBS3D_IMPLICIT_FLOW_PRESSURE_TOL", 1.0e-6);
        const Real divergence_tol =
            env_real("CBS3D_IMPLICIT_FLOW_DIVERGENCE_TOL", 1.0e-7);
        const Int max_nonlinear =
            env_int("CBS3D_IMPLICIT_FLOW_MAXIT", 100);
        const Real requested_relaxation =
            env_real("CBS3D_IMPLICIT_FLOW_RELAX", 0.7);

        if (!(nonlinear_tol > 0.0) ||
            !(pressure_tol > 0.0) ||
            !(divergence_tol > 0.0) ||
            !(requested_relaxation > 0.0) ||
            requested_relaxation > 1.0)
        {
            throw std::runtime_error(
                "Invalid fully implicit CBS nonlinear controls");
        }

        long long local_active = 0;
        for (const Int ip : s.owned_nodes)
        {
            if (touches_fluid(s, ip))
            {
                ++local_active;
            }
        }

        long long global_active = 0;
        long long active_offset = 0;
        check_mpi(MPI_Allreduce(
            &local_active, &global_active, 1, MPI_LONG_LONG, MPI_SUM,
            MPI_COMM_WORLD), "flow active node count");
        check_mpi(MPI_Exscan(
            &local_active, &active_offset, 1, MPI_LONG_LONG, MPI_SUM,
            MPI_COMM_WORLD), "flow active node offset");
        if (s.mpi_rank == 0)
        {
            active_offset = 0;
        }

        const long long global_dofs =
            dofs_per_active_node * global_active;
        if (global_dofs > std::numeric_limits<Int>::max())
        {
            throw std::runtime_error(
                "Fully implicit CBS global DOF count exceeds CBS Int range");
        }

        Array1D<Int> base_dof(s.cfg.npoin);
        base_dof.fill(-1);

        long long local_active_index = 0;
        for (const Int ip : s.owned_nodes)
        {
            if (!touches_fluid(s, ip))
            {
                continue;
            }

            base_dof(ip) = static_cast<Int>(
                dofs_per_active_node *
                (active_offset + local_active_index));
            ++local_active_index;
        }

        HaloExchange::broadcastOwnedToGhosts(
            base_dof, s.partition_metadata, MPI_COMM_WORLD);

        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            if (s.mat_elem(ie) != 0)
            {
                continue;
            }
            for (Int a = 1; a <= 4; ++a)
            {
                if (base_dof(s.intma(a, ie)) < 0)
                {
                    throw std::runtime_error(
                        "Fully implicit CBS: fluid node has no global DOF");
                }
            }
        }

        const PetscInt local_dofs =
            static_cast<PetscInt>(
                dofs_per_active_node * local_active);
        const PetscInt n_global =
            static_cast<PetscInt>(global_dofs);

        Mat A = nullptr;
        Vec x = nullptr;
        Vec b = nullptr;
        Vec residual = nullptr;
        KSP ksp = nullptr;
        IS velocity_is = nullptr;
        IS pressure_is = nullptr;
        IS fixed_is = nullptr;

        check_petsc(MatCreateAIJ(
            MPI_COMM_WORLD,
            local_dofs, local_dofs,
            n_global, n_global,
            96, nullptr, 96, nullptr,
            &A), "MatCreateAIJ(fully implicit CBS)");
        check_petsc(MatSetOption(
            A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE),
            "MatSetOption(flow)");

        check_petsc(MatCreateVecs(A, &x, &b), "MatCreateVecs(flow)");
        check_petsc(VecDuplicate(b, &residual), "VecDuplicate(flow)");

        std::vector<PetscInt> velocity_indices;
        std::vector<PetscInt> pressure_indices;
        std::vector<PetscInt> fixed_indices;

        velocity_indices.reserve(
            static_cast<std::size_t>(6 * local_active));
        pressure_indices.reserve(
            static_cast<std::size_t>(local_active));

        for (const Int ip : s.owned_nodes)
        {
            if (!touches_fluid(s, ip))
            {
                continue;
            }

            const PetscInt base =
                static_cast<PetscInt>(base_dof(ip));

            for (Int q = 0; q < 6; ++q)
            {
                velocity_indices.push_back(base + q);
                if (velocity_fixed(s, ip))
                {
                    fixed_indices.push_back(base + q);
                }
            }

            pressure_indices.push_back(base + pressure_offset);
            if (s.node_pressure_fixed(ip) != 0)
            {
                fixed_indices.push_back(base + pressure_offset);
            }
        }

        check_petsc(ISCreateGeneral(
            MPI_COMM_WORLD,
            static_cast<PetscInt>(velocity_indices.size()),
            velocity_indices.data(),
            PETSC_COPY_VALUES,
            &velocity_is), "ISCreateGeneral(flow velocity)");
        check_petsc(ISCreateGeneral(
            MPI_COMM_WORLD,
            static_cast<PetscInt>(pressure_indices.size()),
            pressure_indices.data(),
            PETSC_COPY_VALUES,
            &pressure_is), "ISCreateGeneral(flow pressure)");
        check_petsc(ISCreateGeneral(
            MPI_COMM_WORLD,
            static_cast<PetscInt>(fixed_indices.size()),
            fixed_indices.data(),
            PETSC_COPY_VALUES,
            &fixed_is), "ISCreateGeneral(flow fixed)");

        check_petsc(KSPCreate(MPI_COMM_WORLD, &ksp), "KSPCreate(flow)");
        check_petsc(KSPSetOptionsPrefix(
            ksp, "flow_"), "KSPSetOptionsPrefix(flow)");
        check_petsc(KSPSetType(
            ksp, KSPFGMRES), "KSPSetType(flow FGMRES)");
        check_petsc(KSPSetTolerances(
            ksp, 1.0e-8, 0.0, PETSC_DEFAULT, 500),
            "KSPSetTolerances(flow)");
        check_petsc(KSPSetInitialGuessNonzero(
            ksp, PETSC_TRUE), "KSPSetInitialGuessNonzero(flow)");

        PC pc = nullptr;
        check_petsc(KSPGetPC(ksp, &pc), "KSPGetPC(flow)");
        check_petsc(PCSetType(
            pc, PCFIELDSPLIT), "PCSetType(flow fieldsplit)");
        check_petsc(PCFieldSplitSetIS(
            pc, "velocity", velocity_is),
            "PCFieldSplitSetIS(velocity)");
        check_petsc(PCFieldSplitSetIS(
            pc, "pressure", pressure_is),
            "PCFieldSplitSetIS(pressure)");
        check_petsc(PCFieldSplitSetType(
            pc, PC_COMPOSITE_SCHUR),
            "PCFieldSplitSetType(Schur)");
        check_petsc(PCFieldSplitSetSchurFactType(
            pc, PC_FIELDSPLIT_SCHUR_FACT_FULL),
            "PCFieldSplitSetSchurFactType");
        check_petsc(PCFieldSplitSetSchurPre(
            pc, PC_FIELDSPLIT_SCHUR_PRE_A11, nullptr),
            "PCFieldSplitSetSchurPre");

        install_default_petsc_options();
        check_petsc(KSPSetFromOptions(
            ksp), "KSPSetFromOptions(flow)");

        if (s.mpi_rank == 0)
        {
            std::cout
                << (backward_euler
                        ? "Fully implicit backward-Euler CBS flow step\n"
                        : "Fully implicit steady CBS flow\n")
                << "  active fluid nodes : " << global_active << "\n"
                << "  coupled DOFs       : " << global_dofs
                << " (u,v,w,u*,v*,w*,p)\n"
                << "  nonlinear method   : Picard/Oseen\n"
                << "  linear method      : PETSc FGMRES + "
                   "PCFIELDSPLIT Schur\n"
                << "  nonlinear tol      : " << nonlinear_tol << "\n"
                << "  pressure tol       : " << pressure_tol << "\n"
                << "  divergence tol     : "
                << divergence_tol << " 1/s\n"
                << "  max Picard iters   : " << max_nonlinear << "\n"
                << "  initial relaxation : "
                << requested_relaxation << "\n";
            if (backward_euler)
            {
                std::cout
                    << "  physical BE dt     : "
                    << physical_dt << " s\n";
            }
        }

        Real relaxation = requested_relaxation;
        Real previous_rel_u = std::numeric_limits<Real>::max();
        Array1D<Real> mass_over_tau(s.cfg.npoin);

        for (Int nonlinear = 1;
             nonlinear <= max_nonlinear;
             ++nonlinear)
        {
            result.nonlinear_iterations = nonlinear;

            HaloExchange::broadcastOwnedToGhosts(
                s.unkno, s.partition_metadata, MPI_COMM_WORLD);
            HaloExchange::broadcastOwnedToGhosts(
                s.pres, s.partition_metadata, MPI_COMM_WORLD);

            mass_over_tau.fill(0.0);

            Real global_tau = physical_dt;

            if (!backward_euler)
            {
                // Steady solve: reproduce the production pseudo-time convention.
                // The pseudo-time is a nonlinear continuation/stabilisation
                // parameter and is not physical elapsed time.
                TimeStep::computeTimeStep(s, nonlinear);

                Real local_tau = std::numeric_limits<Real>::max();
                for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
                {
                    local_tau = std::min(local_tau, s.delte(ie));
                }

                check_mpi(MPI_Allreduce(
                    &local_tau, &global_tau, 1, MPI_DOUBLE, MPI_MIN,
                    MPI_COMM_WORLD), "flow global CBS timestep");

                const Real tau_scale =
                    env_real("CBS3D_IMPLICIT_FLOW_TAU_SCALE", 1.0);
                if (!(tau_scale > 0.0) || !std::isfinite(tau_scale))
                {
                    throw std::runtime_error(
                        "CBS3D_IMPLICIT_FLOW_TAU_SCALE must be positive and finite");
                }
                global_tau *= tau_scale;
            }

            if (!(global_tau > 0.0) || !std::isfinite(global_tau))
            {
                throw std::runtime_error(
                    backward_euler
                        ? "Fully implicit CBS obtained invalid physical BE timestep"
                        : "Fully implicit CBS obtained invalid pseudo-timestep");
            }

            // In BE mode the physical step is used consistently by all three
            // CBS blocks. In steady mode this remains the global pseudo-time.
            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                s.delte(ie) = global_tau;
            }
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                s.deltp(ip) = global_tau;
                s.deltp1(ip) = global_tau;
                s.deltp2(ip) = global_tau;
            }
            s.cfg.dtreal = global_tau;

            // Assemble M_dt = sum_e (V_e/4)/dt on the owner nodes.
            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                if (s.mat_elem(ie) != 0)
                {
                    continue;
                }

                const Real volume = s.detJ(ie) / 6.0;
                const Real contribution = (volume * 0.25) / global_tau;
                for (Int a = 1; a <= 4; ++a)
                {
                    mass_over_tau(s.intma(a, ie)) += contribution;
                }
            }

            HaloExchange::sumGhostContributionsToOwners(
                mass_over_tau, s.partition_metadata, MPI_COMM_WORLD);
            HaloExchange::broadcastOwnedToGhosts(
                mass_over_tau, s.partition_metadata, MPI_COMM_WORLD);

            check_petsc(MatZeroEntries(A), "MatZeroEntries(flow)");
            check_petsc(VecSet(b, 0.0), "VecSet(flow b)");

            // Element operators: A_u, G, B and tau H.
            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                if (s.mat_elem(ie) != 0)
                {
                    continue;
                }

                const Real volume = s.detJ(ie) / 6.0;
                if (!(volume > 0.0) || !std::isfinite(volume))
                {
                    throw std::runtime_error(
                        "Fully implicit CBS: invalid fluid element volume");
                }

                const Real nu = fluid_kinematic_viscosity(s, ie);
                const Real tau = global_tau;

                // Picard field w freezes the advecting velocity in the
                // EXACT discrete Step-1 residual used by MomentumAssembly:
                // nodal conservative flux products, CBS characteristic term,
                // and the componentwise Laplacian diffusion.
                Real w[5][4] = {};
                Real wmean[4] = {};
                Real q_coeff[5] = {};

                for (Int a = 1; a <= 4; ++a)
                {
                    const Int ip = s.intma(a, ie);
                    for (Int dim = 1; dim <= 3; ++dim)
                    {
                        w[a][dim] = s.unkno(dim, ip);
                        wmean[dim] += 0.25 * w[a][dim];
                    }
                }

                // q_coeff[b] is the coefficient multiplying transported
                // component u_i(b) in div(I_h(w_j u_i)):
                //
                //   q_b = sum_j dN_b/dx_j * w_j(b).
                for (Int bnode = 1; bnode <= 4; ++bnode)
                {
                    for (Int j = 1; j <= 3; ++j)
                    {
                        q_coeff[bnode] +=
                            grad(s, ie, j, bnode) * w[bnode][j];
                    }
                }

                PetscInt corrected_rows[12] = {};
                PetscInt corrected_cols[12] = {};
                PetscInt predictor_rows[12] = {};
                PetscInt predictor_cols[12] = {};
                PetscInt pressure_rows[4] = {};
                PetscInt pressure_cols[4] = {};

                PetscScalar a_uu[12 * 12] = {};
                PetscScalar b_volume[4 * 12] = {};
                PetscScalar h_pressure[4 * 4] = {};
                PetscScalar g_pressure[12 * 4] = {};
                PetscScalar b_boundary[4 * 12] = {};
                bool has_step2_boundary_flux = false;

                for (Int a = 1; a <= 4; ++a)
                {
                    const PetscInt base =
                        static_cast<PetscInt>(
                            base_dof(s.intma(a, ie)));

                    pressure_rows[a - 1] = base + pressure_offset;
                    pressure_cols[a - 1] = base + pressure_offset;

                    for (Int i = 1; i <= 3; ++i)
                    {
                        const Int q = (a - 1) * 3 + (i - 1);
                        corrected_rows[q] =
                            base + corrected_offset + (i - 1);
                        corrected_cols[q] =
                            base + corrected_offset + (i - 1);
                        predictor_rows[q] =
                            base + predictor_offset + (i - 1);
                        predictor_cols[q] =
                            base + predictor_offset + (i - 1);
                    }
                }

                for (Int a = 1; a <= 4; ++a)
                {
                    Real test_direction = 0.0;
                    for (Int j = 1; j <= 3; ++j)
                    {
                        test_direction +=
                            grad(s, ie, j, a) * wmean[j];
                    }

                    for (Int bnode = 1; bnode <= 4; ++bnode)
                    {
                        Real grad_dot_ab = 0.0;
                        Real row_dot_advecting_node = 0.0;

                        for (Int j = 1; j <= 3; ++j)
                        {
                            grad_dot_ab +=
                                grad(s, ie, j, a) *
                                grad(s, ie, j, bnode);
                            row_dot_advecting_node +=
                                grad(s, ie, j, a) *
                                w[bnode][j];
                        }

                        // A_u = - d(r_m)/d(u_i) under Picard freezing:
                        //
                        // r_conv volume =
                        //   +(V/4) grad(N_a).sum_b[w(b) u_i(b)]
                        //
                        // r_char volume =
                        //   -(dt/2) V (wbar.grad N_a)
                        //                sum_b[q_b u_i(b)]
                        //
                        // r_diff volume =
                        //   -nu V grad(N_a).grad(u_i).
                        const Real convection =
                            -volume * 0.25 *
                            row_dot_advecting_node;
                        const Real characteristic =
                            0.5 * tau * volume *
                            test_direction * q_coeff[bnode];
                        const Real diffusion =
                            nu * volume * grad_dot_ab;

                        for (Int i = 1; i <= 3; ++i)
                        {
                            const Int row =
                                (a - 1) * 3 + (i - 1);
                            const Int col =
                                (bnode - 1) * 3 + (i - 1);

                            a_uu[row * 12 + col] +=
                                static_cast<PetscScalar>(
                                    convection +
                                    characteristic +
                                    diffusion);
                        }

                        // Step 2: tau H p - B_volume u* = ...
                        for (Int i = 1; i <= 3; ++i)
                        {
                            const Int col =
                                (bnode - 1) * 3 + (i - 1);

                            b_volume[(a - 1) * 12 + col] -=
                                static_cast<PetscScalar>(
                                    volume * 0.25 *
                                    grad(s, ie, i, a));
                        }

                        h_pressure[(a - 1) * 4 + (bnode - 1)] +=
                            static_cast<PetscScalar>(
                                tau * volume * grad_dot_ab);

                        // Step 3: M_dt(u-u*) + Gp = 0.
                        for (Int i = 1; i <= 3; ++i)
                        {
                            const Int row =
                                (a - 1) * 3 + (i - 1);

                            g_pressure[row * 4 + (bnode - 1)] +=
                                static_cast<PetscScalar>(
                                    volume * 0.25 *
                                    grad(s, ie, i, bnode));
                        }
                    }
                }

                // Exact Picard linearisation of production Step-1 boundary
                // terms. MomentumAssembly uses the nodal conservative face
                // flux, the CBS characteristic boundary correction and the
                // componentwise natural diffusion flux.
                for (Int face = 1; face <= 4; ++face)
                {
                    if (s.fedge(face, ie) == 0)
                    {
                        continue;
                    }

                    Real wmean_dot_normal = 0.0;
                    for (Int j = 1; j <= 3; ++j)
                    {
                        wmean_dot_normal +=
                            wmean[j] * s.annxf(j, face, ie);
                    }

                    for (Int fa = 1; fa <= 3; ++fa)
                    {
                        const Int a = s.ippn1(face, fa);

                        // Diffusion and characteristic boundary terms involve
                        // the gradient/divergence reconstructed from all four
                        // element nodes.
                        for (Int bnode = 1; bnode <= 4; ++bnode)
                        {
                            Real grad_dot_normal = 0.0;
                            for (Int j = 1; j <= 3; ++j)
                            {
                                grad_dot_normal +=
                                    grad(s, ie, j, bnode) *
                                    s.annxf(j, face, ie);
                            }

                            const Real diffusion_boundary =
                                nu * s.cfg.fdif[2] *
                                grad_dot_normal;
                            const Real characteristic_boundary =
                                0.5 * tau * s.cfg.fcon[4] *
                                wmean_dot_normal * q_coeff[bnode];

                            for (Int i = 1; i <= 3; ++i)
                            {
                                const Int row =
                                    (a - 1) * 3 + (i - 1);
                                const Int col =
                                    (bnode - 1) * 3 + (i - 1);

                                // Both terms are positive in r_m at the face;
                                // A_u = -r_m.
                                a_uu[row * 12 + col] -=
                                    static_cast<PetscScalar>(
                                        diffusion_boundary +
                                        characteristic_boundary);
                            }
                        }

                        // Conservative convective boundary term:
                        // r_conv has - integral N_a (w_j u_i) n_j.
                        for (Int fb = 1; fb <= 3; ++fb)
                        {
                            const Int bnode =
                                s.ippn1(face, fb);

                            Real advective_flux = 0.0;
                            for (Int j = 1; j <= 3; ++j)
                            {
                                advective_flux +=
                                    w[bnode][j] *
                                    s.annxf(j, face, ie);
                            }

                            const Real face_mass_factor =
                                s.cfg.fcon[2] *
                                (a == bnode ? 2.0 : 1.0);

                            for (Int i = 1; i <= 3; ++i)
                            {
                                const Int row =
                                    (a - 1) * 3 + (i - 1);
                                const Int col =
                                    (bnode - 1) * 3 + (i - 1);

                                a_uu[row * 12 + col] +=
                                    static_cast<PetscScalar>(
                                        face_mass_factor *
                                        advective_flux);
                            }
                        }
                    }

                    // Production Step 2 volume divergence uses u*, while its
                    // boundary flux uses the previous physical velocity in a
                    // transient step. At a steady fixed point u^n == u.
                    if (s.fedge(face, ie) == 1)
                    {
                        has_step2_boundary_flux = true;

                        for (Int fa = 1; fa <= 3; ++fa)
                        {
                            const Int a = s.ippn1(face, fa);

                            for (Int fb = 1; fb <= 3; ++fb)
                            {
                                const Int bnode =
                                    s.ippn1(face, fb);

                                const Real face_mass_factor =
                                    s.cfg.fcon[2] *
                                    (a == bnode ? 2.0 : 1.0);

                                for (Int i = 1; i <= 3; ++i)
                                {
                                    const Int col =
                                        (bnode - 1) * 3 + (i - 1);

                                    b_boundary[(a - 1) * 12 + col] +=
                                        static_cast<PetscScalar>(
                                            face_mass_factor *
                                            s.annxf(i, face, ie));
                                }
                            }
                        }
                    }
                }

                // Insert only the mathematical nonzero blocks. Avoiding the
                // former dense 28x28 structural insertion is essential on the
                // multi-million-tetrahedron production blanket.
                check_petsc(MatSetValues(
                    A, 12, corrected_rows,
                    12, corrected_cols,
                    a_uu, ADD_VALUES),
                    "MatSetValues(flow Auu)");

                check_petsc(MatSetValues(
                    A, 4, pressure_rows,
                    12, predictor_cols,
                    b_volume, ADD_VALUES),
                    "MatSetValues(flow Bvolume)");

                check_petsc(MatSetValues(
                    A, 4, pressure_rows,
                    4, pressure_cols,
                    h_pressure, ADD_VALUES),
                    "MatSetValues(flow Hpressure)");

                check_petsc(MatSetValues(
                    A, 12, predictor_rows,
                    4, pressure_cols,
                    g_pressure, ADD_VALUES),
                    "MatSetValues(flow Gpressure)");

                if (has_step2_boundary_flux)
                {
                    if (backward_euler)
                    {
                        PetscScalar pressure_boundary_rhs[4] = {};
                        for (Int a = 0; a < 4; ++a)
                        {
                            for (Int bnode = 1; bnode <= 4; ++bnode)
                            {
                                const Int ip = s.intma(bnode, ie);
                                for (Int dim = 1; dim <= 3; ++dim)
                                {
                                    const Int col =
                                        (bnode - 1) * 3 + (dim - 1);
                                    pressure_boundary_rhs[a] -=
                                        b_boundary[a * 12 + col] *
                                        static_cast<PetscScalar>(
                                            s.unkn1(dim, ip));
                                }
                            }
                        }

                        check_petsc(VecSetValues(
                            b, 4, pressure_rows,
                            pressure_boundary_rhs, ADD_VALUES),
                            "VecSetValues(flow BE pressure boundary RHS)");
                    }
                    else
                    {
                        check_petsc(MatSetValues(
                            A, 4, pressure_rows,
                            12, corrected_cols,
                            b_boundary, ADD_VALUES),
                            "MatSetValues(flow Bboundary)");
                    }
                }
            }

            // Add the assembled lumped M_dt coupling exactly once per owned
            // active node. It appears in Step 1 and Step 3 with opposite signs.
            for (const Int ip : s.owned_nodes)
            {
                if (!touches_fluid(s, ip))
                {
                    continue;
                }

                const Real mdt = mass_over_tau(ip);
                if (!(mdt > 0.0) || !std::isfinite(mdt))
                {
                    throw std::runtime_error(
                        "Fully implicit CBS: invalid assembled M_dt");
                }

                const PetscInt base =
                    static_cast<PetscInt>(base_dof(ip));

                for (Int dim = 0; dim < 3; ++dim)
                {
                    const PetscInt step1_row =
                        base + corrected_offset + dim;
                    const PetscInt step3_row =
                        base + predictor_offset + dim;
                    const PetscInt u_col =
                        base + corrected_offset + dim;
                    const PetscInt predictor_col =
                        base + predictor_offset + dim;

                    if (!backward_euler)
                    {
                        check_petsc(MatSetValue(
                            A, step1_row, u_col,
                            -mdt, ADD_VALUES),
                            "MatSetValue(Step1 steady -Mdt)");
                    }

                    check_petsc(MatSetValue(
                        A, step1_row, predictor_col,
                        mdt, ADD_VALUES),
                        "MatSetValue(Step1 +Mdt)");

                    if (backward_euler)
                    {
                        const PetscScalar old_mass_rhs =
                            static_cast<PetscScalar>(
                                mdt * s.unkn1(dim + 1, ip));
                        check_petsc(VecSetValue(
                            b, step1_row, old_mass_rhs, ADD_VALUES),
                            "VecSetValue(Step1 BE old-velocity RHS)");
                    }

                    check_petsc(MatSetValue(
                        A, step3_row, u_col,
                        mdt, ADD_VALUES),
                        "MatSetValue(Step3 +Mdt)");
                    check_petsc(MatSetValue(
                        A, step3_row, predictor_col,
                        -mdt, ADD_VALUES),
                        "MatSetValue(Step3 -Mdt)");
                }
            }

            check_petsc(MatAssemblyBegin(
                A, MAT_FINAL_ASSEMBLY),
                "MatAssemblyBegin(flow)");
            check_petsc(MatAssemblyEnd(
                A, MAT_FINAL_ASSEMBLY),
                "MatAssemblyEnd(flow)");

            // Current corrected velocity is also a good initial guess for u*.
            check_petsc(VecSet(x, 0.0), "VecSet(flow x)");

            PetscScalar* x_values = nullptr;
            check_petsc(VecGetArray(
                x, &x_values), "VecGetArray(flow x)");

            long long cursor = 0;
            for (const Int ip : s.owned_nodes)
            {
                if (!touches_fluid(s, ip))
                {
                    continue;
                }

                for (Int dim = 0; dim < 3; ++dim)
                {
                    const Real value = s.unkno(dim + 1, ip);
                    x_values[
                        cursor + corrected_offset + dim] = value;
                    x_values[
                        cursor + predictor_offset + dim] = value;
                }
                x_values[cursor + pressure_offset] = s.pres(ip);

                cursor += dofs_per_active_node;
            }

            check_petsc(VecRestoreArray(
                x, &x_values), "VecRestoreArray(flow x)");
            check_petsc(VecAssemblyBegin(
                x), "VecAssemblyBegin(flow x)");
            check_petsc(VecAssemblyEnd(
                x), "VecAssemblyEnd(flow x)");
            check_petsc(VecAssemblyBegin(
                b), "VecAssemblyBegin(flow b)");
            check_petsc(VecAssemblyEnd(
                b), "VecAssemblyEnd(flow b)");

            check_petsc(MatZeroRowsColumnsIS(
                A, fixed_is, 1.0, x, b),
                "MatZeroRowsColumnsIS(flow BC)");

            check_petsc(KSPSetOperators(
                ksp, A, A), "KSPSetOperators(flow)");
            check_petsc(KSPSolve(
                ksp, b, x), "KSPSolve(flow)");

            KSPConvergedReason reason;
            PetscInt linear_iterations = 0;

            check_petsc(KSPGetConvergedReason(
                ksp, &reason), "KSPGetConvergedReason(flow)");
            check_petsc(KSPGetIterationNumber(
                ksp, &linear_iterations),
                "KSPGetIterationNumber(flow)");

            if (reason <= 0)
            {
                std::ostringstream message;
                message
                    << "Fully implicit CBS linear solve failed at Picard "
                    << nonlinear
                    << ": KSP reason=" << static_cast<int>(reason)
                    << ", iterations=" << linear_iterations;
                throw std::runtime_error(message.str());
            }

            result.total_krylov_iterations +=
                static_cast<Int>(linear_iterations);

            check_petsc(MatMult(
                A, x, residual), "MatMult(flow audit)");
            check_petsc(VecAYPX(
                residual, -1.0, b), "VecAYPX(flow audit)");

            PetscReal rn = 0.0;
            PetscReal bn = 0.0;

            check_petsc(VecNorm(
                residual, NORM_2, &rn), "VecNorm(flow residual)");
            check_petsc(VecNorm(
                b, NORM_2, &bn), "VecNorm(flow rhs)");

            result.linear_relative_residual =
                bn > 0.0 ? rn / bn : rn;

            if (!std::isfinite(result.linear_relative_residual))
            {
                throw std::runtime_error(
                    "Fully implicit CBS produced nonfinite linear residual");
            }

            const Metrics metrics =
                update_state_from_solution(s, x, relaxation);

            result.relative_velocity_change = metrics.rel_u;
            result.relative_pressure_change = metrics.rel_p;
            result.divergence_rms = metrics.div_rms;
            result.maximum_velocity = metrics.umax;
            result.maximum_velocity_correction = metrics.correction_max;

            const Real tau_min = global_tau;
            const Real tau_max = global_tau;

            if (s.mpi_rank == 0)
            {
                std::cout
                    << std::setprecision(8)
                    << "[IMPLICIT-FLOW] Picard " << nonlinear
                    << "  KSP=" << linear_iterations
                    << "  KSPRel=" << result.linear_relative_residual
                    << "  RelU=" << metrics.rel_u
                    << "  RelP=" << metrics.rel_p
                    << "  DivRMS=" << metrics.div_rms
                    << "  Umax=" << metrics.umax
                    << "  dUmax=" << metrics.correction_max
                    << "  omega=" << relaxation
                    << "  tau=[" << tau_min
                    << "," << tau_max << "]\n"
                    << std::flush;
            }

            if (metrics.rel_u < nonlinear_tol &&
                metrics.rel_p < pressure_tol &&
                metrics.div_rms < divergence_tol)
            {
                result.converged = true;
                break;
            }

            if (nonlinear > 2 &&
                std::isfinite(previous_rel_u) &&
                metrics.rel_u > 1.25 * previous_rel_u)
            {
                relaxation = std::max(
                    static_cast<Real>(0.2),
                    static_cast<Real>(0.5) * relaxation);
            }

            previous_rel_u = metrics.rel_u;
        }

        result.seconds = MPI_Wtime() - start;

        check_petsc(ISDestroy(
            &fixed_is), "ISDestroy(flow fixed)");
        check_petsc(ISDestroy(
            &pressure_is), "ISDestroy(flow pressure)");
        check_petsc(ISDestroy(
            &velocity_is), "ISDestroy(flow velocity)");
        check_petsc(KSPDestroy(
            &ksp), "KSPDestroy(flow)");
        check_petsc(VecDestroy(
            &residual), "VecDestroy(flow residual)");
        check_petsc(VecDestroy(
            &b), "VecDestroy(flow b)");
        check_petsc(VecDestroy(
            &x), "VecDestroy(flow x)");
        check_petsc(MatDestroy(
            &A), "MatDestroy(flow)");

        return result;
    }

    FullyImplicitFlowSolver::Result
    FullyImplicitFlowSolver::solve(CBSStateSI& s)
    {
        return solveImpl(s, false, 0.0);
    }

    FullyImplicitFlowSolver::Result
    FullyImplicitFlowSolver::solveBackwardEulerStep(
        CBSStateSI& s,
        const Real physical_dt)
    {
        return solveImpl(s, true, physical_dt);
    }
}

#else

namespace cbs
{
    FullyImplicitFlowSolver::Result
    FullyImplicitFlowSolver::solve(CBSStateSI&)
    {
        throw std::runtime_error(
            "Fully implicit CBS requires MPI + PETSc");
    }

    FullyImplicitFlowSolver::Result
    FullyImplicitFlowSolver::solveBackwardEulerStep(
        CBSStateSI&,
        Real)
    {
        throw std::runtime_error(
            "Fully implicit CBS backward Euler requires MPI + PETSc");
    }
}

#endif
