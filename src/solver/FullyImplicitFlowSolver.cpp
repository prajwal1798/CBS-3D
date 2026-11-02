//=============================================================================
// CBS3D++_SI
//
// Fully implicit steady laminar flow solve for the distributed P1/TET4 system.
//
// Nonlinearity is treated by Picard/Oseen iteration. At nonlinear iteration k,
// convection and CBS characteristic stabilisation are frozen on u^k while
// u^(k+1) and p^(k+1) are solved simultaneously:
//
//   [ A_Oseen(u^k)   G ] [u] = 0
//   [      -D       tauH] [p] = 0.
//
// tauH is the existing CBS equal-order pressure stabilisation. The velocity
// block uses the same conservative convection and full deviatoric viscous
// operator as MomentumAssemblyLiu.cpp. Strong constraints are imposed with
// MatZeroRowsColumnsIS, so prescribed inlet/wall/interface values also enter
// neighbouring free equations consistently.
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
        constexpr Int velocity_components = 3;
        constexpr Int dofs_per_active_node = 4;

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
                        "Fully implicit flow: invalid fluid density");
                }

                const Real mu =
                    s.cfg.turbulence_on > 0 ? s.mu_eff_e(ie) : s.mu_e(ie);

                if (!(mu >= 0.0) || !std::isfinite(mu))
                {
                    throw std::runtime_error(
                        "Fully implicit flow: invalid fluid viscosity");
                }
                return mu / s.rho_e(ie);
            }

            if (!(s.cfg.ani >= 0.0) || !std::isfinite(s.cfg.ani))
            {
                throw std::runtime_error(
                    "Fully implicit flow: invalid nondimensional viscosity");
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
                    "Fully implicit flow: invalid h/nu for CBS stabilisation");
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

            const Real advective =
                h / std::max(speed, static_cast<Real>(1.0e-14));
            const Real diffusive =
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

            const Real tau = scale * s.cfg.csafm *
                std::min(advective, diffusive);

            if (!(tau > 0.0) || !std::isfinite(tau))
            {
                throw std::runtime_error(
                    "Fully implicit flow: invalid CBS stabilisation tau");
            }
            return tau;
        }

        Real tetra_mass(const Real volume, const Int a, const Int b)
        {
            return volume * (a == b ? 2.0 : 1.0) / 20.0;
        }

        void set_default_option(
            const char* name,
            const char* value)
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
            set_default_option("-flow_ksp_max_it", "400");
            set_default_option("-flow_pc_type", "fieldsplit");
            set_default_option("-flow_pc_fieldsplit_type", "schur");
            set_default_option(
                "-flow_pc_fieldsplit_schur_fact_type", "full");
            set_default_option(
                "-flow_pc_fieldsplit_schur_precondition", "a11");

            // Oseen velocity block: local overlapping ILU is robust for the
            // nonsymmetric convection-diffusion operator at the present scale.
            set_default_option(
                "-flow_fieldsplit_velocity_ksp_type", "gmres");
            set_default_option(
                "-flow_fieldsplit_velocity_ksp_rtol", "1e-2");
            set_default_option(
                "-flow_fieldsplit_velocity_ksp_max_it", "30");
            set_default_option(
                "-flow_fieldsplit_velocity_pc_type", "asm");
            set_default_option(
                "-flow_fieldsplit_velocity_pc_asm_overlap", "1");
            set_default_option(
                "-flow_fieldsplit_velocity_sub_pc_type", "ilu");
            set_default_option(
                "-flow_fieldsplit_velocity_sub_pc_factor_levels", "1");

            // Pressure Schur approximation: the CBS pressure block is an
            // elliptic tau*Laplacian, for which BoomerAMG is appropriate.
            set_default_option(
                "-flow_fieldsplit_pressure_ksp_type", "gmres");
            set_default_option(
                "-flow_fieldsplit_pressure_ksp_rtol", "1e-2");
            set_default_option(
                "-flow_fieldsplit_pressure_ksp_max_it", "30");
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

                const PetscInt expected =
                    local_start + static_cast<PetscInt>(cursor);
                (void)expected;

                for (Int dim = 1; dim <= 3; ++dim)
                {
                    const Real old_value = s.unkno(dim, ip);
                    const Real linear_value =
                        static_cast<Real>(values[cursor + (dim - 1)]);
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
                    static_cast<Real>(values[cursor + 3]);
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
                    "Fully implicit flow: local vector/DOF ownership mismatch");
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
                        div_u += s.unkno(dim, ip) * grad(s, ie, dim, a);
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

            check_mpi(
                MPI_Allreduce(
                    &local_du2, &global_du2, 1, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD),
                "implicit flow du norm");
            check_mpi(
                MPI_Allreduce(
                    &local_u2, &global_u2, 1, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD),
                "implicit flow u norm");
            check_mpi(
                MPI_Allreduce(
                    &local_dp2, &global_dp2, 1, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD),
                "implicit flow dp norm");
            check_mpi(
                MPI_Allreduce(
                    &local_p2, &global_p2, 1, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD),
                "implicit flow p norm");
            check_mpi(
                MPI_Allreduce(
                    &local_div2, &global_div2, 1, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD),
                "implicit flow divergence norm");
            check_mpi(
                MPI_Allreduce(
                    &local_volume, &global_volume, 1, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD),
                "implicit flow fluid volume");
            check_mpi(
                MPI_Allreduce(
                    &local_umax, &global_umax, 1, MPI_DOUBLE, MPI_MAX,
                    MPI_COMM_WORLD),
                "implicit flow Umax");

            Metrics m;
            m.rel_u = std::sqrt(global_du2) /
                std::max(std::sqrt(global_u2), static_cast<Real>(1.0e-30));
            m.rel_p = std::sqrt(global_dp2) /
                std::max(std::sqrt(global_p2), static_cast<Real>(1.0));
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
                "Fully implicit flow requires distributed MPI execution");
        }
        if (s.cfg.ndim != 3 || s.cfg.nep != 4)
        {
            throw std::runtime_error(
                "Fully implicit flow currently requires 3D P1/TET4");
        }
        if (s.cfg.turbulence_on > 0)
        {
            throw std::runtime_error(
                "Fully implicit flow baseline is laminar only; implicit SA "
                "coupling has not yet been added");
        }

        for (const Int ip : s.owned_nodes)
        {
            if (s.node_symmetry(ip) != 0)
            {
                throw std::runtime_error(
                    "Fully implicit flow currently supports strong velocity "
                    "conditions but not distributed symmetry constraints");
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
                "Invalid fully implicit flow nonlinear controls");
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
        check_mpi(
            MPI_Allreduce(
                &local_active, &global_active, 1, MPI_LONG_LONG, MPI_SUM,
                MPI_COMM_WORLD),
            "implicit flow active node count");
        check_mpi(
            MPI_Exscan(
                &local_active, &active_offset, 1, MPI_LONG_LONG, MPI_SUM,
                MPI_COMM_WORLD),
            "implicit flow active node offset");
        if (s.mpi_rank == 0)
        {
            active_offset = 0;
        }

        const long long global_dofs =
            dofs_per_active_node * global_active;
        if (global_dofs > std::numeric_limits<Int>::max())
        {
            throw std::runtime_error(
                "Fully implicit flow global DOF count exceeds CBS Int range");
        }

        Array1D<Int> base_dof(s.cfg.npoin);
        base_dof.fill(-1);
        long long local_index = 0;
        for (const Int ip : s.owned_nodes)
        {
            if (!touches_fluid(s, ip))
            {
                continue;
            }
            base_dof(ip) = static_cast<Int>(
                dofs_per_active_node * (active_offset + local_index));
            ++local_index;
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
                        "Fully implicit flow: fluid node has no global DOF");
                }
            }
        }

        const PetscInt local_dofs =
            static_cast<PetscInt>(dofs_per_active_node * local_active);
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

        check_petsc(
            MatCreateAIJ(
                MPI_COMM_WORLD,
                local_dofs, local_dofs,
                n_global, n_global,
                64, nullptr, 64, nullptr,
                &A),
            "MatCreateAIJ(implicit flow)");
        check_petsc(
            MatSetOption(
                A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE),
            "MatSetOption(implicit flow)");

        check_petsc(MatCreateVecs(A, &x, &b), "MatCreateVecs");
        check_petsc(VecDuplicate(b, &residual), "VecDuplicate");

        std::vector<PetscInt> velocity_indices;
        std::vector<PetscInt> pressure_indices;
        std::vector<PetscInt> fixed_indices;

        velocity_indices.reserve(
            static_cast<std::size_t>(3 * local_active));
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
            for (Int dim = 0; dim < 3; ++dim)
            {
                velocity_indices.push_back(base + dim);
                if (velocity_fixed(s, ip))
                {
                    fixed_indices.push_back(base + dim);
                }
            }
            pressure_indices.push_back(base + 3);
            if (s.node_pressure_fixed(ip) != 0)
            {
                fixed_indices.push_back(base + 3);
            }
        }

        check_petsc(
            ISCreateGeneral(
                MPI_COMM_WORLD,
                static_cast<PetscInt>(velocity_indices.size()),
                velocity_indices.data(),
                PETSC_COPY_VALUES,
                &velocity_is),
            "ISCreateGeneral(velocity)");
        check_petsc(
            ISCreateGeneral(
                MPI_COMM_WORLD,
                static_cast<PetscInt>(pressure_indices.size()),
                pressure_indices.data(),
                PETSC_COPY_VALUES,
                &pressure_is),
            "ISCreateGeneral(pressure)");
        check_petsc(
            ISCreateGeneral(
                MPI_COMM_WORLD,
                static_cast<PetscInt>(fixed_indices.size()),
                fixed_indices.data(),
                PETSC_COPY_VALUES,
                &fixed_is),
            "ISCreateGeneral(fixed)");

        check_petsc(KSPCreate(MPI_COMM_WORLD, &ksp), "KSPCreate");
        check_petsc(
            KSPSetOptionsPrefix(ksp, "flow_"),
            "KSPSetOptionsPrefix(flow)");
        check_petsc(KSPSetType(ksp, KSPFGMRES), "KSPSetType(FGMRES)");
        check_petsc(
            KSPSetTolerances(
                ksp, 1.0e-8, 0.0, PETSC_DEFAULT, 400),
            "KSPSetTolerances(flow)");
        check_petsc(
            KSPSetInitialGuessNonzero(ksp, PETSC_TRUE),
            "KSPSetInitialGuessNonzero(flow)");

        PC pc = nullptr;
        check_petsc(KSPGetPC(ksp, &pc), "KSPGetPC(flow)");
        check_petsc(PCSetType(pc, PCFIELDSPLIT), "PCSetType(fieldsplit)");
        check_petsc(
            PCFieldSplitSetIS(pc, "velocity", velocity_is),
            "PCFieldSplitSetIS(velocity)");
        check_petsc(
            PCFieldSplitSetIS(pc, "pressure", pressure_is),
            "PCFieldSplitSetIS(pressure)");
        check_petsc(
            PCFieldSplitSetType(pc, PC_COMPOSITE_SCHUR),
            "PCFieldSplitSetType(Schur)");
        check_petsc(
            PCFieldSplitSetSchurFactType(
                pc, PC_FIELDSPLIT_SCHUR_FACT_FULL),
            "PCFieldSplitSetSchurFactType");
        check_petsc(
            PCFieldSplitSetSchurPre(
                pc, PC_FIELDSPLIT_SCHUR_PRE_A11, nullptr),
            "PCFieldSplitSetSchurPre");

        install_default_petsc_options();
        check_petsc(KSPSetFromOptions(ksp), "KSPSetFromOptions(flow)");

        if (s.mpi_rank == 0)
        {
            std::cout
                << "Fully implicit steady flow\n"
                << "  active fluid nodes : " << global_active << "\n"
                << "  coupled DOFs       : " << global_dofs
                << " (u,v,w,p)\n"
                << "  nonlinear method   : Picard/Oseen\n"
                << "  linear method      : PETSc FGMRES + "
                   "PCFIELDSPLIT Schur\n"
                << "  nonlinear tol      : " << nonlinear_tol << "\n"
                << "  pressure tol       : " << pressure_tol << "\n"
                << "  divergence tol     : " << divergence_tol << " 1/s\n"
                << "  max Picard iters   : " << max_nonlinear << "\n"
                << "  initial relaxation : "
                << requested_relaxation << "\n";
        }

        Real relaxation = requested_relaxation;
        Real previous_rel_u = std::numeric_limits<Real>::max();

        for (Int nonlinear = 1;
             nonlinear <= max_nonlinear;
             ++nonlinear)
        {
            result.nonlinear_iterations = nonlinear;

            HaloExchange::broadcastOwnedToGhosts(
                s.unkno, s.partition_metadata, MPI_COMM_WORLD);
            HaloExchange::broadcastOwnedToGhosts(
                s.pres, s.partition_metadata, MPI_COMM_WORLD);

            check_petsc(MatZeroEntries(A), "MatZeroEntries(flow)");

            Real local_tau_min = std::numeric_limits<Real>::max();
            Real local_tau_max = 0.0;

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
                        "Fully implicit flow: invalid fluid element volume");
                }

                const Real nu = fluid_kinematic_viscosity(s, ie);
                const Real tau = element_cbs_tau(s, ie);
                local_tau_min = std::min(local_tau_min, tau);
                local_tau_max = std::max(local_tau_max, tau);

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

                PetscInt rows[16] = {};
                PetscInt cols[16] = {};
                PetscScalar values[16 * 16] = {};

                for (Int a = 1; a <= 4; ++a)
                {
                    const PetscInt base =
                        static_cast<PetscInt>(
                            base_dof(s.intma(a, ie)));
                    for (Int q = 0; q < 4; ++q)
                    {
                        rows[(a - 1) * 4 + q] = base + q;
                        cols[(a - 1) * 4 + q] = base + q;
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

                        for (Int i = 1; i <= 3; ++i)
                        {
                            const Int row =
                                (a - 1) * 4 + (i - 1);

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
                                    (bnode - 1) * 4 + (k - 1);
                                values[row * 16 + col] +=
                                    static_cast<PetscScalar>(
                                        diffusion +
                                        (i == k
                                            ? convection +
                                              characteristic
                                            : 0.0));
                            }

                            // Kinematic pressure gradient:
                            // integral N_a * dp/dx_i dV.
                            const Int pcol =
                                (bnode - 1) * 4 + 3;
                            values[row * 16 + pcol] +=
                                static_cast<PetscScalar>(
                                    volume * 0.25 *
                                    grad(s, ie, i, bnode));
                        }

                        // CBS stabilised continuity:
                        // tau H p - [volume weak divergence] = 0.
                        const Int prow = (a - 1) * 4 + 3;
                        for (Int i = 1; i <= 3; ++i)
                        {
                            const Int ucol =
                                (bnode - 1) * 4 + (i - 1);
                            values[prow * 16 + ucol] -=
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
                        const Int pcol =
                            (bnode - 1) * 4 + 3;
                        values[prow * 16 + pcol] +=
                            static_cast<PetscScalar>(
                                tau * volume * h_ab);
                    }
                }

                // Match the production weak boundary operators.
                for (Int face = 1; face <= 4; ++face)
                {
                    if (s.fedge(face, ie) == 0)
                    {
                        continue;
                    }

                    // Momentum: subtract the natural deviatoric traction
                    // operator because the production residual adds +tau.n.
                    for (Int fa = 1; fa <= 3; ++fa)
                    {
                        const Int a = s.ippn1(face, fa);
                        for (Int bnode = 1; bnode <= 4; ++bnode)
                        {
                            for (Int i = 1; i <= 3; ++i)
                            {
                                const Int row =
                                    (a - 1) * 4 + (i - 1);
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
                                        (bnode - 1) * 4 + (k - 1);
                                    values[row * 16 + col] -=
                                        static_cast<PetscScalar>(traction);
                                }
                            }
                        }
                    }

                    // Continuity boundary flux is present only on the
                    // production fedge==1 faces.
                    if (s.fedge(face, ie) == 1)
                    {
                        for (Int fa = 1; fa <= 3; ++fa)
                        {
                            const Int a = s.ippn1(face, fa);
                            const Int prow = (a - 1) * 4 + 3;

                            for (Int fb = 1; fb <= 3; ++fb)
                            {
                                const Int bnode =
                                    s.ippn1(face, fb);
                                const Real face_mass_factor =
                                    s.cfg.fcon[2] *
                                    (a == bnode ? 2.0 : 1.0);

                                for (Int i = 1; i <= 3; ++i)
                                {
                                    const Int ucol =
                                        (bnode - 1) * 4 + (i - 1);
                                    values[prow * 16 + ucol] +=
                                        static_cast<PetscScalar>(
                                            face_mass_factor *
                                            s.annxf(i, face, ie));
                                }
                            }
                        }
                    }
                }

                check_petsc(
                    MatSetValues(
                        A, 16, rows, 16, cols, values, ADD_VALUES),
                    "MatSetValues(implicit flow element)");
            }

            check_petsc(
                MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY),
                "MatAssemblyBegin(implicit flow)");
            check_petsc(
                MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY),
                "MatAssemblyEnd(implicit flow)");

            // Pack current iterate. Fixed entries carry their prescribed value
            // and are used by MatZeroRowsColumnsIS to shift free-equation RHS.
            check_petsc(VecSet(x, 0.0), "VecSet(x)");
            check_petsc(VecSet(b, 0.0), "VecSet(b)");

            PetscScalar* x_values = nullptr;
            check_petsc(VecGetArray(x, &x_values), "VecGetArray(x)");
            long long cursor = 0;
            for (const Int ip : s.owned_nodes)
            {
                if (!touches_fluid(s, ip))
                {
                    continue;
                }
                x_values[cursor + 0] = s.unkno(1, ip);
                x_values[cursor + 1] = s.unkno(2, ip);
                x_values[cursor + 2] = s.unkno(3, ip);
                x_values[cursor + 3] = s.pres(ip);
                cursor += 4;
            }
            check_petsc(VecRestoreArray(x, &x_values), "VecRestoreArray(x)");
            check_petsc(
                VecAssemblyBegin(x), "VecAssemblyBegin(x)");
            check_petsc(
                VecAssemblyEnd(x), "VecAssemblyEnd(x)");
            check_petsc(
                VecAssemblyBegin(b), "VecAssemblyBegin(b)");
            check_petsc(
                VecAssemblyEnd(b), "VecAssemblyEnd(b)");

            check_petsc(
                MatZeroRowsColumnsIS(A, fixed_is, 1.0, x, b),
                "MatZeroRowsColumnsIS(flow BC)");

            check_petsc(KSPSetOperators(ksp, A, A), "KSPSetOperators(flow)");
            check_petsc(KSPSolve(ksp, b, x), "KSPSolve(flow)");

            KSPConvergedReason reason;
            PetscInt linear_iterations = 0;
            check_petsc(
                KSPGetConvergedReason(ksp, &reason),
                "KSPGetConvergedReason(flow)");
            check_petsc(
                KSPGetIterationNumber(ksp, &linear_iterations),
                "KSPGetIterationNumber(flow)");

            if (reason <= 0)
            {
                std::ostringstream message;
                message
                    << "Fully implicit flow linear solve failed at Picard "
                    << nonlinear
                    << ": KSP reason=" << static_cast<int>(reason)
                    << ", iterations=" << linear_iterations;
                throw std::runtime_error(message.str());
            }

            result.total_krylov_iterations +=
                static_cast<Int>(linear_iterations);

            check_petsc(
                MatMult(A, x, residual), "MatMult(flow audit)");
            check_petsc(
                VecAYPX(residual, -1.0, b), "VecAYPX(flow audit)");

            PetscReal rn = 0.0;
            PetscReal bn = 0.0;
            check_petsc(
                VecNorm(residual, NORM_2, &rn),
                "VecNorm(flow residual)");
            check_petsc(
                VecNorm(b, NORM_2, &bn),
                "VecNorm(flow rhs)");
            result.linear_relative_residual =
                bn > 0.0 ? rn / bn : rn;

            if (!std::isfinite(result.linear_relative_residual))
            {
                throw std::runtime_error(
                    "Fully implicit flow produced nonfinite linear residual");
            }

            const Metrics metrics =
                update_state_from_solution(s, x, relaxation);

            result.relative_velocity_change = metrics.rel_u;
            result.relative_pressure_change = metrics.rel_p;
            result.divergence_rms = metrics.div_rms;
            result.maximum_velocity = metrics.umax;

            Real tau_min = 0.0;
            Real tau_max = 0.0;
            check_mpi(
                MPI_Allreduce(
                    &local_tau_min, &tau_min, 1, MPI_DOUBLE, MPI_MIN,
                    MPI_COMM_WORLD),
                "implicit flow tau min");
            check_mpi(
                MPI_Allreduce(
                    &local_tau_max, &tau_max, 1, MPI_DOUBLE, MPI_MAX,
                    MPI_COMM_WORLD),
                "implicit flow tau max");

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
                    << "  tau=[" << tau_min << "," << tau_max << "]\n"
                    << std::flush;
            }

            if (metrics.rel_u < nonlinear_tol &&
                metrics.rel_p < pressure_tol &&
                metrics.div_rms < divergence_tol)
            {
                result.converged = true;
                break;
            }

            // Conservative automatic damping if Picard starts to move away
            // from its fixed point. Never increase above the requested value.
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

        check_petsc(ISDestroy(&fixed_is), "ISDestroy(fixed)");
        check_petsc(ISDestroy(&pressure_is), "ISDestroy(pressure)");
        check_petsc(ISDestroy(&velocity_is), "ISDestroy(velocity)");
        check_petsc(KSPDestroy(&ksp), "KSPDestroy(flow)");
        check_petsc(VecDestroy(&residual), "VecDestroy(flow residual)");
        check_petsc(VecDestroy(&b), "VecDestroy(flow b)");
        check_petsc(VecDestroy(&x), "VecDestroy(flow x)");
        check_petsc(MatDestroy(&A), "MatDestroy(flow)");

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
            "Fully implicit flow requires MPI + PETSc");
    }
}

#endif
