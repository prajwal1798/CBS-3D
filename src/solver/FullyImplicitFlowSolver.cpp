//=============================================================================
// CBS3D++_SI
//
// Fully implicit steady CBS flow solve.
//
// The original semi-implicit CBS iteration is
//
//   M_dt (u* - u) = r_m(u)
//   A_p p          = B u*
//   M_dt (u - u*) + G p = 0
//
// where A_p = sum_e tau_e H_e and M_dt is the lumped momentum mass divided by
// the same local CBS pseudo-time. Instead of marching these equations millions
// of times, this solver treats u, u* and p as simultaneous unknowns. Picard
// linearisation freezes the nonlinear convection/stabilisation coefficients on
// the current corrected velocity. Every Picard iteration therefore solves the
// exact linearised CBS Step-1/2/3 fixed-point equations in one PETSc system.
//
// Per active fluid node the unknown ordering is
//
//   [u, v, w, u*, v*, w*, p].
//
// Temperature is intentionally not included: for the present constant-property
// forced-flow blanket problem there is no T -> (u,p) feedback. The converged
// flow is followed by the existing implicit steady CHT solve.
//=============================================================================

#include "cbs/solver/FullyImplicitFlowSolver.hpp"

#include <stdexcept>

#if defined(CBS3D_USE_MPI) && defined(CBS3D_USE_PETSC)

#include "cbs/boundary/Boundary.hpp"
#include "cbs/linalg/PetscHeaders.hpp"
#include "cbs/parallel/HaloExchange.hpp"

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
                throw std::runtime_error(
                    std::string("PETSc failure in ") + what +
                    " (code " + std::to_string(static_cast<int>(code)) + ")");
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

        Real element_cbs_tau(const CBSStateSI& s, const Int ie)
        {
            const Real h = s.alen_e(ie);
            const Real nu = fluid_kinematic_viscosity(s, ie);

            if (!(h > 0.0) || !std::isfinite(h) ||
                !(nu >= 0.0) || !std::isfinite(nu))
            {
                throw std::runtime_error(
                    "Fully implicit CBS: invalid h/nu for stabilisation");
            }

            Real speed = 0.0;
            for (Int a = 1; a <= 4; ++a)
            {
                const Int ip = s.intma(a, ie);
                const Real ux = s.unkno(1, ip);
                const Real uy = s.unkno(2, ip);
                const Real uz = s.unkno(3, ip);
                speed = std::max(
                    speed,
                    std::sqrt(ux * ux + uy * uy + uz * uz));
            }

            const Real dt_adv =
                h / std::max(speed, static_cast<Real>(1.0e-14));
            const Real dt_diff =
                nu > 0.0
                    ? h * h / (2.0 * nu)
                    : std::numeric_limits<Real>::max();

            const Real scale =
                env_real("CBS3D_IMPLICIT_FLOW_TAU_SCALE", 1.0);
            if (!(scale > 0.0))
            {
                throw std::runtime_error(
                    "CBS3D_IMPLICIT_FLOW_TAU_SCALE must be positive");
            }

            const Real tau =
                scale * s.cfg.csafm * std::min(dt_adv, dt_diff);

            if (!(tau > 0.0) || !std::isfinite(tau))
            {
                throw std::runtime_error(
                    "Fully implicit CBS: invalid local tau");
            }
            return tau;
        }

        Real tetra_mass(const Real volume, const Int a, const Int b)
        {
            return volume * (a == b ? 2.0 : 1.0) / 20.0;
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
            return m;
        }
    }

    FullyImplicitFlowSolver::Result
    FullyImplicitFlowSolver::solve(CBSStateSI& s)
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
                << "Fully implicit steady CBS flow\n"
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

            Real local_tau_min = std::numeric_limits<Real>::max();
            Real local_tau_max = 0.0;

            // Assemble M_dt = sum_e (V_e/4)/tau_e on the owner nodes.
            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                if (s.mat_elem(ie) != 0)
                {
                    continue;
                }

                const Real volume = s.detJ(ie) / 6.0;
                const Real tau = element_cbs_tau(s, ie);
                local_tau_min = std::min(local_tau_min, tau);
                local_tau_max = std::max(local_tau_max, tau);

                const Real contribution = (volume * 0.25) / tau;
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
                const Real tau = element_cbs_tau(s, ie);

                Real w[5][4] = {};
                Real div_w = 0.0;
                for (Int a = 1; a <= 4; ++a)
                {
                    const Int ip = s.intma(a, ie);
                    for (Int dim = 1; dim <= 3; ++dim)
                    {
                        w[a][dim] = s.unkno(dim, ip);
                        div_w +=
                            w[a][dim] * grad(s, ie, dim, a);
                    }
                }

                // g_b(node) = div(w N_b) evaluated at the four vertices.
                Real g[5][5] = {};
                for (Int basis = 1; basis <= 4; ++basis)
                {
                    for (Int node = 1; node <= 4; ++node)
                    {
                        Real value = 0.0;
                        for (Int dim = 1; dim <= 3; ++dim)
                        {
                            value +=
                                w[node][dim] *
                                grad(s, ie, dim, basis);
                        }
                        if (node == basis)
                        {
                            value += div_w;
                        }
                        g[basis][node] = value;
                    }
                }

                PetscInt rows[28] = {};
                PetscInt cols[28] = {};
                PetscScalar values[28 * 28] = {};

                for (Int a = 1; a <= 4; ++a)
                {
                    const PetscInt base =
                        static_cast<PetscInt>(
                            base_dof(s.intma(a, ie)));
                    for (Int q = 0; q < dofs_per_active_node; ++q)
                    {
                        rows[(a - 1) * 7 + q] = base + q;
                        cols[(a - 1) * 7 + q] = base + q;
                    }
                }

                for (Int a = 1; a <= 4; ++a)
                {
                    for (Int bnode = 1; bnode <= 4; ++bnode)
                    {
                        Real convection = 0.0;
                        for (Int c = 1; c <= 4; ++c)
                        {
                            convection +=
                                tetra_mass(volume, a, c) *
                                g[bnode][c];
                        }

                        Real characteristic = 0.0;
                        for (Int c = 1; c <= 4; ++c)
                        {
                            for (Int d = 1; d <= 4; ++d)
                            {
                                characteristic +=
                                    g[a][c] *
                                    tetra_mass(volume, c, d) *
                                    g[bnode][d];
                            }
                        }
                        characteristic *= 0.5 * tau;

                        // Step 1 equation:
                        // A_u(w) u + M_dt(u* - u) = 0.
                        for (Int i = 1; i <= 3; ++i)
                        {
                            const Int row =
                                (a - 1) * 7 +
                                corrected_offset + (i - 1);

                            for (Int k = 1; k <= 3; ++k)
                            {
                                Real diffusion = 0.0;
                                for (Int j = 1; j <= 3; ++j)
                                {
                                    if (i == k)
                                    {
                                        diffusion +=
                                            grad(s, ie, j, a) *
                                            grad(s, ie, j, bnode);
                                    }
                                }

                                diffusion +=
                                    grad(s, ie, k, a) *
                                    grad(s, ie, i, bnode);
                                diffusion -=
                                    (2.0 / 3.0) *
                                    grad(s, ie, i, a) *
                                    grad(s, ie, k, bnode);
                                diffusion *= nu * volume;

                                const Int col =
                                    (bnode - 1) * 7 +
                                    corrected_offset + (k - 1);

                                values[row * 28 + col] +=
                                    static_cast<PetscScalar>(
                                        diffusion +
                                        (i == k
                                            ? convection +
                                              characteristic
                                            : 0.0));
                            }
                        }

                        // Step 2 equation:
                        // tau H p - B u* = 0.
                        const Int pressure_row =
                            (a - 1) * 7 + pressure_offset;

                        for (Int i = 1; i <= 3; ++i)
                        {
                            const Int predictor_col =
                                (bnode - 1) * 7 +
                                predictor_offset + (i - 1);

                            values[
                                pressure_row * 28 +
                                predictor_col] -=
                                static_cast<PetscScalar>(
                                    volume * 0.25 *
                                    grad(s, ie, i, a));
                        }

                        Real h_ab = 0.0;
                        for (Int j = 1; j <= 3; ++j)
                        {
                            h_ab +=
                                grad(s, ie, j, a) *
                                grad(s, ie, j, bnode);
                        }

                        const Int pressure_col =
                            (bnode - 1) * 7 + pressure_offset;
                        values[
                            pressure_row * 28 +
                            pressure_col] +=
                            static_cast<PetscScalar>(
                                tau * volume * h_ab);

                        // Step 3 equation pressure gradient:
                        // M_dt(u-u*) + G p = 0.
                        for (Int i = 1; i <= 3; ++i)
                        {
                            const Int correction_row =
                                (a - 1) * 7 +
                                predictor_offset + (i - 1);

                            values[
                                correction_row * 28 +
                                pressure_col] +=
                                static_cast<PetscScalar>(
                                    volume * 0.25 *
                                    grad(s, ie, i, bnode));
                        }
                    }
                }

                // Production Step-1 natural deviatoric traction.
                for (Int face = 1; face <= 4; ++face)
                {
                    if (s.fedge(face, ie) == 0)
                    {
                        continue;
                    }

                    for (Int fa = 1; fa <= 3; ++fa)
                    {
                        const Int a = s.ippn1(face, fa);

                        for (Int bnode = 1; bnode <= 4; ++bnode)
                        {
                            for (Int i = 1; i <= 3; ++i)
                            {
                                const Int row =
                                    (a - 1) * 7 +
                                    corrected_offset + (i - 1);

                                for (Int k = 1; k <= 3; ++k)
                                {
                                    Real traction = 0.0;

                                    for (Int j = 1; j <= 3; ++j)
                                    {
                                        if (i == k)
                                        {
                                            traction +=
                                                grad(s, ie, j, bnode) *
                                                s.annxf(j, face, ie);
                                        }
                                    }

                                    traction +=
                                        grad(s, ie, i, bnode) *
                                        s.annxf(k, face, ie);

                                    traction -=
                                        (2.0 / 3.0) *
                                        grad(s, ie, k, bnode) *
                                        s.annxf(i, face, ie);

                                    traction *= nu / 3.0;

                                    const Int col =
                                        (bnode - 1) * 7 +
                                        corrected_offset + (k - 1);

                                    // r_m contains + traction; A_u=-r_m.
                                    values[row * 28 + col] -=
                                        static_cast<PetscScalar>(traction);
                                }
                            }
                        }
                    }

                    // Production Step-2 weak divergence boundary term:
                    // B u* = volume term - boundary flux.
                    if (s.fedge(face, ie) == 1)
                    {
                        for (Int fa = 1; fa <= 3; ++fa)
                        {
                            const Int a = s.ippn1(face, fa);
                            const Int pressure_row =
                                (a - 1) * 7 + pressure_offset;

                            for (Int fb = 1; fb <= 3; ++fb)
                            {
                                const Int bnode =
                                    s.ippn1(face, fb);

                                const Real face_mass_factor =
                                    s.cfg.fcon[2] *
                                    (a == bnode ? 2.0 : 1.0);

                                // Production Step 2 is deliberately
                                // mixed at the boundary: its volume divergence
                                // uses the predictor u*, while the face flux
                                // uses unkn1 (the previous corrected velocity).
                                // At the steady fixed point unkn1 == u, so the
                                // monolithic equation must couple this boundary
                                // term to corrected u, not to u*.
                                for (Int i = 1; i <= 3; ++i)
                                {
                                    const Int corrected_col =
                                        (bnode - 1) * 7 +
                                        corrected_offset + (i - 1);

                                    values[
                                        pressure_row * 28 +
                                        corrected_col] +=
                                        static_cast<PetscScalar>(
                                            face_mass_factor *
                                            s.annxf(i, face, ie));
                                }
                            }
                        }
                    }
                }

                check_petsc(MatSetValues(
                    A, 28, rows, 28, cols, values, ADD_VALUES),
                    "MatSetValues(fully implicit CBS element)");
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

                    check_petsc(MatSetValue(
                        A, step1_row, u_col,
                        -mdt, ADD_VALUES),
                        "MatSetValue(Step1 -Mdt)");
                    check_petsc(MatSetValue(
                        A, step1_row, predictor_col,
                        mdt, ADD_VALUES),
                        "MatSetValue(Step1 +Mdt)");
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
            check_petsc(VecSet(b, 0.0), "VecSet(flow b)");

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

            Real tau_min = 0.0;
            Real tau_max = 0.0;
            check_mpi(MPI_Allreduce(
                &local_tau_min, &tau_min, 1, MPI_DOUBLE, MPI_MIN,
                MPI_COMM_WORLD), "flow tau min");
            check_mpi(MPI_Allreduce(
                &local_tau_max, &tau_max, 1, MPI_DOUBLE, MPI_MAX,
                MPI_COMM_WORLD), "flow tau max");

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
}

#endif
