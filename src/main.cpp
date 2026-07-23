//=============================================================================
// CBS3D++_SI
//
// Program entry point for serial, OpenMP and MPI execution.
//
// Serial/OpenMP execution:
//
//     cbs3dpp_si <case_name>
//
// MPI persistent distributed-loop execution:
//
//     cbs3dpp_si <case_name> [partition_root]
//
// For an MPI run with N ranks, the default partition root is:
//
//     cbs_partitions_N
//
// Rank r reads the local case:
//
//     <partition_root>/rank_rrrr/<case_name>_rank_rrrr
//
// where rrrr is the four-digit MPI rank number.
//=============================================================================

#include "cbs/solver/Solver.hpp"

#ifdef CBS3D_USE_PETSC
#include "cbs/linalg/PetscPressureSolver.hpp"
#endif

#ifdef CBS3D_USE_MPI
#include <mpi.h>
#endif

#include <exception>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace
{
#ifdef CBS3D_USE_MPI
    // Broadcasts one text string from rank zero to every MPI rank.
    void broadcastString(std::string& value, const int rank)
    {
        int length = 0;

        if (rank == 0)
        {
            length = static_cast<int>(value.size());
        }

        if (MPI_Bcast(&length, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
        {
            throw std::runtime_error(
                "MPI_Bcast failed while broadcasting a string length");
        }

        if (length < 0)
        {
            throw std::runtime_error(
                "A negative string length was received through MPI");
        }

        if (rank != 0)
        {
            value.resize(static_cast<std::size_t>(length));
        }

        if (length > 0)
        {
            if (MPI_Bcast(
                    value.data(),
                    length,
                    MPI_CHAR,
                    0,
                    MPI_COMM_WORLD) != MPI_SUCCESS)
            {
                throw std::runtime_error(
                    "MPI_Bcast failed while broadcasting a string");
            }
        }
    }


    // Returns the MPI rank as a four-digit text field.
    std::string rankText(const int rank)
    {
        std::ostringstream text;
        text << std::setw(4) << std::setfill('0') << rank;
        return text.str();
    }


    // Builds the base name of the rank-local CBS input files.
    std::string localCaseName(
        const std::string& case_name,
        const std::string& partition_root,
        const int rank)
    {
        const std::string rank_text = rankText(rank);

        return partition_root
             + "/rank_" + rank_text
             + "/" + case_name
             + "_rank_" + rank_text;
    }
#endif
}


int main(int argc, char** argv)
{
    int mpi_rank = 0;
    int mpi_size = 1;

#ifdef CBS3D_USE_MPI
    int mpi_thread_level = MPI_THREAD_SINGLE;

    // Only the main thread performs MPI communication. OpenMP worker threads
    // remain inside local finite-element loops.
    if (MPI_Init_thread(
            &argc,
            &argv,
            MPI_THREAD_FUNNELED,
            &mpi_thread_level) != MPI_SUCCESS)
    {
        std::cerr << "ERROR: MPI_Init_thread failed\n";
        return 1;
    }

    if (mpi_thread_level < MPI_THREAD_FUNNELED)
    {
        std::cerr
            << "ERROR: MPI_THREAD_FUNNELED is not available\n";

        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }

    if (MPI_Comm_set_errhandler(
            MPI_COMM_WORLD,
            MPI_ERRORS_RETURN) != MPI_SUCCESS)
    {
        std::cerr << "ERROR: unable to set the MPI error handler\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }

    if (MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank) != MPI_SUCCESS ||
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size) != MPI_SUCCESS)
    {
        std::cerr
            << "ERROR: unable to obtain the MPI rank or communicator size\n";

        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
#endif

    try
    {
        std::string case_name;
        std::string partition_root;

        // Rank zero reads the command line. The selected values are then
        // broadcast to all MPI ranks.
        if (mpi_rank == 0)
        {
            if (argc >= 2)
            {
                case_name = argv[1];
            }
            else
            {
                std::cout << "Enter Problem Name?\n";
                std::cin >> case_name;
            }

            if (case_name.empty())
            {
                throw std::runtime_error(
                    "The problem name must not be empty");
            }

            if (argc >= 3)
            {
                partition_root = argv[2];
            }

#ifdef CBS3D_USE_MPI
            if (mpi_size > 1 && partition_root.empty())
            {
                partition_root =
                    "cbs_partitions_" + std::to_string(mpi_size);
            }
#endif
        }

#ifdef CBS3D_USE_MPI
        broadcastString(case_name, mpi_rank);
        broadcastString(partition_root, mpi_rank);
#endif

        std::string solver_case_name = case_name;

#ifdef CBS3D_USE_MPI
        if (mpi_size > 1)
        {
            solver_case_name = localCaseName(
                case_name,
                partition_root,
                mpi_rank);
        }
#endif

        if (mpi_rank == 0)
        {
            std::cout << "Problem Name: " << case_name << "\n";

#ifdef CBS3D_USE_MPI
            if (mpi_size > 1)
            {
                std::cout
                    << "MPI ranks: " << mpi_size << "\n"
                    << "Partition root: " << partition_root << "\n";
            }
#endif
        }

        cbs::Solver solver(solver_case_name);

        solver.setMpiContext(
            static_cast<cbs::Int>(mpi_rank),
            static_cast<cbs::Int>(mpi_size));

#ifdef CBS3D_USE_MPI
        if (mpi_size > 1)
        {
            // Preprocess once, build the distributed PETSc matrix and AMG
            // hierarchy once, then advance the complete distributed CBS loop.
            solver.runDistributedLoop();
        }
        else
#endif
        {
            solver.run();
        }

#ifdef CBS3D_USE_PETSC
        // Releases the legacy serial PETSc cache when that path was used.
        cbs::PetscPressureSolver::shutdown();
#endif

#ifdef CBS3D_USE_MPI
        if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS)
        {
            throw std::runtime_error(
                "MPI_Barrier failed before finalisation");
        }

        if (MPI_Finalize() != MPI_SUCCESS)
        {
            std::cerr << "ERROR: MPI_Finalize failed\n";
            return 1;
        }
#endif
    }
    catch (const std::exception& e)
    {
        std::cerr
            << "\nERROR on MPI rank "
            << mpi_rank
            << ": "
            << e.what()
            << "\n";

#ifdef CBS3D_USE_MPI
        MPI_Abort(MPI_COMM_WORLD, 1);
#endif
        return 1;
    }

    return 0;
}
