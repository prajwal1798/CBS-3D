//=============================================================================
// CBS3D++_SI serial application
//
// Supports:
//   - deterministic single-thread execution;
//   - OpenMP shared-memory execution when enabled at build time;
//   - the internal pressure CG solver;
//   - the serial PETSc pressure path when enabled at build time.
//
// Usage:
//
//     cbs3d_serial <case_name>
//=============================================================================

#include "cbs/solver/Solver.hpp"

#ifdef CBS3D_USE_PETSC
#include "cbs/linalg/PetscPressureSolver.hpp"
#endif

#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{
    std::string read_case_name(const int argc, char** argv)
    {
        if (argc > 2)
        {
            throw std::runtime_error(
                "Usage: cbs3d_serial <case_name>");
        }

        if (argc == 2)
        {
            const std::string case_name = argv[1];

            if (case_name.empty())
            {
                throw std::runtime_error(
                    "The problem name must not be empty");
            }

            return case_name;
        }

        std::cout << "Enter Problem Name?\n";

        std::string case_name;
        std::cin >> case_name;

        if (case_name.empty())
        {
            throw std::runtime_error(
                "The problem name must not be empty");
        }

        return case_name;
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
    int return_code = 0;

    try
    {
        const std::string case_name = read_case_name(argc, argv);

        std::cout
            << "CBS3D++_SI serial | case="
            << case_name
            << "\n";

        cbs::Solver solver(case_name);
        solver.setMpiContext(0, 1);
        solver.run();
    }
    catch (const std::exception& error)
    {
        std::cerr
            << "\nERROR: "
            << error.what()
            << "\n";

        return_code = 1;
    }

    shutdown_optional_petsc_pressure_path();
    return return_code;
}
