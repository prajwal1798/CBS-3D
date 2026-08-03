#include "cbs/solver/Solver.hpp"

#include "cbs/assembly/PressureAssembly.hpp"
#include "cbs/linalg/PetscPressureSolver.hpp"
#include "cbs/parallel/HaloExchange.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#ifdef CBS3D_USE_MPI
#include <mpi.h>
#endif

namespace cbs
{
#ifdef CBS3D_USE_MPI
    namespace
    {
        void check_step2_mpi(int ierr, const char* operation)
        {
            if (ierr != MPI_SUCCESS)
            {
                throw std::runtime_error(
                    std::string("Solver::runDistributedStep2 - MPI failure in ")
                    + operation);
            }
        }
    }
#endif

    void Solver::runDistributedStep2()
    {
#if defined(CBS3D_USE_MPI) && defined(CBS3D_USE_PETSC)
        if (!s_.mpi_enabled)
        {
            throw std::runtime_error(
                "Solver::runDistributedStep2 requires more than one MPI rank");
        }

        // Step 1 leaves u* owner/ghost consistent.
        runDistributedStep1();

        // Geometry-only P1 pressure stiffness coefficients.
        PressureAssembly::buildElementPressureTerms(s_);

        // Rank-owned tetrahedra assemble the weak divergence locally.
        PressureAssembly::assembleStep2Rhs(s_);

        // Complete the distributed FE vector assembly on shared-node owners.
        HaloExchange::sumGhostContributionsToOwners(
            s_.rhs1,
            s_.partition_metadata,
            MPI_COMM_WORLD);

        const ConjugateGradient::Result result =
            PetscPressureSolver::solveDistributedPressure(s_);

        s_.last_cg_iterations = result.iterations;
        s_.last_cg_initial_l2 = result.initial_l2;
        s_.last_cg_final_l2 = result.final_l2;
        s_.last_cg_relative_l2 = result.final_relative_l2;
        s_.last_cg_max_abs = result.final_max_abs;

        if (!result.converged)
        {
            throw std::runtime_error(
                "Solver::runDistributedStep2 - distributed PETSc pressure "
                "solve failed to converge");
        }

        // PETSc writes only owner pressure values; update every ghost copy.
        HaloExchange::broadcastOwnedToGhosts(
            s_.pres,
            s_.partition_metadata,
            MPI_COMM_WORLD);

        Real local_pressure_l2_sq = 0.0;
        Real local_pressure_max = 0.0;

        for (const Int ip : s_.owned_nodes)
        {
            if ((s_.node_material_mask(ip) &
                 CBSStateSI::node_touches_fluid) == 0)
            {
                continue;
            }

            const Real pressure = s_.pres(ip);

            if (!std::isfinite(pressure))
            {
                throw std::runtime_error(
                    "Solver::runDistributedStep2 produced non-finite pressure");
            }

            local_pressure_l2_sq += pressure * pressure;
            local_pressure_max =
                std::max(local_pressure_max, std::abs(pressure));
        }

        Real global_pressure_l2_sq = 0.0;
        Real global_pressure_max = 0.0;

        check_step2_mpi(
            MPI_Allreduce(
                &local_pressure_l2_sq,
                &global_pressure_l2_sq,
                1,
                MPI_DOUBLE,
                MPI_SUM,
                MPI_COMM_WORLD),
            "MPI_Allreduce pressure L2 norm");

        check_step2_mpi(
            MPI_Allreduce(
                &local_pressure_max,
                &global_pressure_max,
                1,
                MPI_DOUBLE,
                MPI_MAX,
                MPI_COMM_WORLD),
            "MPI_Allreduce pressure maximum");

        if (s_.mpi_rank == 0)
        {
            std::cout
                << "============================================================\n"
                << "CBS3D DISTRIBUTED STEP 2 COMPLETE\n"
                << "============================================================\n"
                << "MPI ranks              : " << s_.mpi_size << "\n"
                << "PETSc CG iterations    : " << result.iterations << "\n"
                << "initial residual L2    : " << result.initial_l2 << "\n"
                << "final residual L2      : " << result.final_l2 << "\n"
                << "relative residual L2   : " << result.final_relative_l2 << "\n"
                << "maximum residual       : " << result.final_max_abs << "\n"
                << "pressure L2 norm       : "
                << std::sqrt(global_pressure_l2_sq) << "\n"
                << "maximum |pressure|     : " << global_pressure_max << "\n"
                << "next solver stage      : distributed velocity correction\n"
                << "============================================================\n";
        }
#else
        throw std::runtime_error(
            "Solver::runDistributedStep2 requires MPI and PETSc support");
#endif
    }
}
