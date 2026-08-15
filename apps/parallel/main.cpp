//=============================================================================
// CBS3D++_SI distributed application
//
// Supports the validated MPI domain-decomposition path with PETSc pressure
// solution. OpenMP may also be enabled at build time for future hybrid runs.
//
// Usage:
//
//     cbs3d_parallel <case_name> [partition_root]
//
// Rank r reads:
//
//     <partition_root>/rank_rrrr/<case_name>_rank_rrrr.*
//=============================================================================

#include "cbs/solver/Solver.hpp"

#ifdef CBS3D_USE_PETSC
#include "cbs/linalg/PetscPressureSolver.hpp"
#endif

#include <mpi.h>

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <streambuf>
#include <string>

namespace
{
    class NullStreamBuffer final : public std::streambuf
    {
    protected:
        int overflow(const int character) override
        {
            return traits_type::not_eof(character);
        }
    };


    class ScopedStdoutSilence final
    {
    public:
        explicit ScopedStdoutSilence(const bool suppress)
        {
            if (suppress)
            {
                previous_ = std::cout.rdbuf(&null_buffer_);
            }
        }

        ~ScopedStdoutSilence()
        {
            if (previous_ != nullptr)
            {
                std::cout.rdbuf(previous_);
            }
        }

        ScopedStdoutSilence(const ScopedStdoutSilence&) = delete;
        ScopedStdoutSilence& operator=(const ScopedStdoutSilence&) = delete;

    private:
        NullStreamBuffer null_buffer_;
        std::streambuf* previous_ = nullptr;
    };


    bool environment_flag_enabled(const char* name)
    {
        const char* value = std::getenv(name);

        return value != nullptr &&
               value[0] != '\0' &&
               std::string(value) != "0";
    }


    void check_mpi(const int error_code, const char* operation)
    {
        if (error_code != MPI_SUCCESS)
        {
            throw std::runtime_error(
                std::string("MPI failure in ") + operation);
        }
    }


    void broadcast_string(std::string& value, const int rank)
    {
        int length = 0;

        if (rank == 0)
        {
            length = static_cast<int>(value.size());
        }

        check_mpi(
            MPI_Bcast(
                &length,
                1,
                MPI_INT,
                0,
                MPI_COMM_WORLD),
            "broadcasting a string length");

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
            check_mpi(
                MPI_Bcast(
                    value.data(),
                    length,
                    MPI_CHAR,
                    0,
                    MPI_COMM_WORLD),
                "broadcasting a string");
        }
    }


    std::string rank_text(const int rank)
    {
        std::ostringstream text;
        text << std::setw(4) << std::setfill('0') << rank;
        return text.str();
    }


    std::string local_case_name(
        const std::string& case_name,
        const std::string& partition_root,
        const int rank)
    {
        const std::string rank_id = rank_text(rank);

        return partition_root
             + "/rank_" + rank_id
             + "/" + case_name
             + "_rank_" + rank_id;
    }


    std::string lower_copy(std::string text)
    {
        std::transform(
            text.begin(),
            text.end(),
            text.begin(),
            [](const unsigned char character)
            {
                return static_cast<char>(std::tolower(character));
            });

        return text;
    }


    bool blank_or_comment(const std::string& line)
    {
        const auto first = line.find_first_not_of(" \t\r\n");

        return first == std::string::npos ||
               line[first] == '#' ||
               line[first] == '!';
    }


    // Temporary production safety check.
    //
    // The distributed production loop is currently laminar. The repository
    // already contains serial Spalart-Allmaras kernels, so a turbulent MPI case
    // could otherwise be mistaken for a supported distributed calculation.
    // Older laminar parameter files do not contain the optional SA block and
    // therefore return turbulence_on=0.
    int read_turbulence_flag(const std::string& parameter_path)
    {
        std::ifstream input(parameter_path);

        if (!input)
        {
            throw std::runtime_error(
                "Unable to open rank-local parameter file: "
                + parameter_path);
        }

        std::string line;

        while (std::getline(input, line))
        {
            const std::string lower = lower_copy(line);

            if (lower.find("turbulence_on") == std::string::npos)
            {
                continue;
            }

            while (std::getline(input, line))
            {
                if (blank_or_comment(line))
                {
                    continue;
                }

                std::istringstream values(line);
                int turbulence_on = 0;

                if (!(values >> turbulence_on))
                {
                    throw std::runtime_error(
                        "Invalid Spalart-Allmaras control data after the "
                        "turbulence_on label in " + parameter_path);
                }

                return turbulence_on;
            }

            throw std::runtime_error(
                "Missing Spalart-Allmaras control data after the "
                "turbulence_on label in " + parameter_path);
        }

        return 0;
    }


    void check_consistent_turbulence_flag(
        const std::string& solver_case_name,
        const int mpi_size)
    {
        if (mpi_size <= 1)
        {
            return;
        }

        const int local_flag =
            read_turbulence_flag(solver_case_name + ".par");

        int minimum_flag = 0;
        int maximum_flag = 0;

        check_mpi(
            MPI_Allreduce(
                &local_flag,
                &minimum_flag,
                1,
                MPI_INT,
                MPI_MIN,
                MPI_COMM_WORLD),
            "checking the minimum turbulence flag");

        check_mpi(
            MPI_Allreduce(
                &local_flag,
                &maximum_flag,
                1,
                MPI_INT,
                MPI_MAX,
                MPI_COMM_WORLD),
            "checking the maximum turbulence flag");

        if (minimum_flag != maximum_flag)
        {
            throw std::runtime_error(
                "Rank-local parameter files contain inconsistent "
                "turbulence_on values");
        }

        // Distributed Spalart-Allmaras is implemented in the production MPI
        // loop: the wall-distance search gathers the global wall surface, the
        // SA node classification is reduced over shared nodes, and the
        // transport step sums its element-assembled accumulators onto their
        // owners before the nodal update.  The blanket rejection that used to
        // stand here is therefore removed.
        //
        // The remaining check is a consistency check only.  It verifies that
        // every rank read the same turbulence_on value from its parameter file,
        // because a mismatch would make some ranks execute the SA step and its
        // collectives while others skipped them, which deadlocks rather than
        // failing.  That check is performed above and applies whatever the
        // value is.
        //
        // Model selection is validated in the solver, not here, so that the
        // serial and distributed paths cannot disagree about which models
        // exist.
    }


    void shutdown_optional_petsc_pressure_path() noexcept
    {
#ifdef CBS3D_USE_PETSC
        try
        {
            cbs::PetscPressureSolver::shutdown();
        }
        catch (...)
        {
            // Do not hide the original solver result during process teardown.
        }
#endif
    }
}


int main(int argc, char** argv)
{
    int mpi_rank = 0;
    int mpi_size = 1;
    int mpi_thread_level = MPI_THREAD_SINGLE;

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

    const ScopedStdoutSilence rank_stdout(
        mpi_rank != 0 &&
        !environment_flag_enabled("CBS3D_ALL_RANK_OUTPUT"));

    try
    {
        std::string case_name;
        std::string partition_root;

        if (mpi_rank == 0)
        {
            if (argc > 3)
            {
                throw std::runtime_error(
                    "Usage: cbs3d_parallel <case_name> [partition_root]");
            }

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

            if (mpi_size > 1 && partition_root.empty())
            {
                partition_root =
                    "cbs_partitions_" + std::to_string(mpi_size);
            }
        }

        broadcast_string(case_name, mpi_rank);
        broadcast_string(partition_root, mpi_rank);

        std::string solver_case_name = case_name;

        if (mpi_size > 1)
        {
            solver_case_name = local_case_name(
                case_name,
                partition_root,
                mpi_rank);

            check_consistent_turbulence_flag(
                solver_case_name,
                mpi_size);
        }

        if (mpi_rank == 0)
        {
            std::cout
                << "CBS3D++_SI parallel | case="
                << case_name
                << " | MPI ranks="
                << mpi_size;

            if (mpi_size > 1)
            {
                std::cout
                    << " | partitions="
                    << partition_root;
            }

            std::cout << "\n";
        }

        cbs::Solver solver(solver_case_name);
        solver.setMpiContext(
            static_cast<cbs::Int>(mpi_rank),
            static_cast<cbs::Int>(mpi_size));

        if (mpi_size > 1)
        {
            solver.runDistributedProductionLoop();
        }
        else
        {
            solver.run();
        }

        shutdown_optional_petsc_pressure_path();

        check_mpi(
            MPI_Barrier(MPI_COMM_WORLD),
            "the final barrier");

        if (MPI_Finalize() != MPI_SUCCESS)
        {
            std::cerr << "ERROR: MPI_Finalize failed\n";
            return 1;
        }
    }
    catch (const std::exception& error)
    {
        std::cerr
            << "\nERROR on MPI rank "
            << mpi_rank
            << ": "
            << error.what()
            << "\n";

        shutdown_optional_petsc_pressure_path();
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }

    return 0;
}
