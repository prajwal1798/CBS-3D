#!/usr/bin/env python3
"""Apply concise production-console output to the distributed CBS3D solver.

This transformation deliberately does not remove numerical validation, runtime
checks, exceptions, residual CSV output, or classic VTU/PVTU/PVD output.

Default behaviour after the patch:

* only MPI rank zero writes ordinary stdout;
* rank-local diagnostic chatter remains suppressed;
* distributed setup internals remain silent;
* the user sees one startup summary, interval-controlled iteration lines and one
  completion summary;
* stderr and fatal exceptions remain visible on every rank.

Debug opt-ins:

* CBS3D_VERBOSE=1         expose rank-zero setup diagnostics;
* CBS3D_ALL_RANK_OUTPUT=1 expose stdout from every MPI rank.

Compatible with Python 3.6 on Swansea Sunbird.
"""

from pathlib import Path
import sys


MAIN_PATH = Path("src/main.cpp")
PRODUCTION_PATH = Path("src/solver/SolverDistributedProductionLoop.cpp")


def replace_once(text, old, new, label):
    count = text.count(old)
    if count != 1:
        raise RuntimeError(
            "{}: expected one source match, found {}".format(label, count)
        )
    return text.replace(old, new, 1)


def patch_main(text):
    if "CBS3D_ALL_RANK_OUTPUT" in text:
        return text, False

    text = replace_once(
        text,
        "#include <exception>\n#include <iomanip>\n",
        "#include <cstdlib>\n#include <exception>\n#include <iomanip>\n",
        "main include insertion",
    )

    text = replace_once(
        text,
        "#include <stdexcept>\n#include <string>\n",
        "#include <stdexcept>\n#include <streambuf>\n#include <string>\n",
        "main streambuf include insertion",
    )

    helper = r'''namespace
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


    bool environmentFlagEnabled(const char* name)
    {
        const char* value = std::getenv(name);

        return value != nullptr &&
               value[0] != '\0' &&
               std::string(value) != "0";
    }

#ifdef CBS3D_USE_MPI
'''

    text = replace_once(
        text,
        "namespace\n{\n#ifdef CBS3D_USE_MPI\n",
        helper,
        "main console helper insertion",
    )

    text = replace_once(
        text,
        "#endif\n\n    try\n",
        "#endif\n\n    const ScopedStdoutSilence rank_stdout(\n        mpi_rank != 0 &&\n        !environmentFlagEnabled(\"CBS3D_ALL_RANK_OUTPUT\"));\n\n    try\n",
        "main rank-output suppression",
    )

    old_startup = r'''        if (mpi_rank == 0)
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
'''

    new_startup = r'''        if (mpi_rank == 0)
        {
            std::cout
                << "CBS3D++_SI | case=" << case_name;

#ifdef CBS3D_USE_MPI
            if (mpi_size > 1)
            {
                std::cout
                    << " | MPI ranks=" << mpi_size
                    << " | partitions=" << partition_root;
            }
#endif

            std::cout << "\n";
        }
'''

    text = replace_once(
        text,
        old_startup,
        new_startup,
        "main concise startup",
    )

    return text, True


def patch_production_loop(text):
    if "CBS3D_VERBOSE" in text:
        return text, False

    text = replace_once(
        text,
        "#include <cmath>\n#include <iostream>\n",
        "#include <cmath>\n#include <cstdlib>\n#include <iostream>\n",
        "production include insertion",
    )

    text = replace_once(
        text,
        "#include <stdexcept>\n#include <string>\n",
        "#include <stdexcept>\n#include <streambuf>\n#include <string>\n",
        "production streambuf include insertion",
    )

    helper = r'''        constexpr Real residual_epsilon = 1.0e-30;

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
'''

    text = replace_once(
        text,
        "        constexpr Real residual_epsilon = 1.0e-30;\n",
        helper,
        "production console helper insertion",
    )

    old_setup = r'''        const double setup_start = MPI_Wtime();

        runDistributedPreprocessing();

        Preprocess::classifyFaceEdges(s_);
        Preprocess::elementSize(s_);
        Coloring::build(s_);

        reject_unsupported_options(s_);

        PressureAssembly::buildElementPressureTerms(s_);

        ThermalBoundaryState thermal_boundary(s_.cfg.npoin);

        if (s_.cfg.temp_calc > 0)
        {
            thermal_boundary = build_thermal_boundary_state(s_);
        }

        Boundary::applyOwnedVelocityConstraints(s_);
        apply_owned_pressure_constraints(s_);

        if (s_.cfg.temp_calc > 0)
        {
            apply_owned_temperature_constraints(s_, thermal_boundary);
        }

        HaloExchange::broadcastOwnedToGhosts(
            s_.unkno,
            s_.partition_metadata,
            MPI_COMM_WORLD);

        HaloExchange::broadcastOwnedToGhosts(
            s_.pres,
            s_.partition_metadata,
            MPI_COMM_WORLD);

        if (s_.cfg.temp_calc > 0)
        {
            HaloExchange::broadcastOwnedToGhosts(
                s_.temperature,
                s_.partition_metadata,
                MPI_COMM_WORLD);
        }

        updateVelocityMagnitude();

        PetscPersistentDistributedPressureSystem pressure_system;
        pressure_system.initialise(s_);
'''

    new_setup = r'''        const double setup_start = MPI_Wtime();

        ThermalBoundaryState thermal_boundary(s_.cfg.npoin);
        PetscPersistentDistributedPressureSystem pressure_system;

        {
            const ScopedStdoutSilence setup_stdout(
                s_.mpi_rank != 0 ||
                !environment_flag_enabled("CBS3D_VERBOSE"));

            runDistributedPreprocessing();

            Preprocess::classifyFaceEdges(s_);
            Preprocess::elementSize(s_);
            Coloring::build(s_);

            reject_unsupported_options(s_);

            PressureAssembly::buildElementPressureTerms(s_);

            if (s_.cfg.temp_calc > 0)
            {
                thermal_boundary = build_thermal_boundary_state(s_);
            }

            Boundary::applyOwnedVelocityConstraints(s_);
            apply_owned_pressure_constraints(s_);

            if (s_.cfg.temp_calc > 0)
            {
                apply_owned_temperature_constraints(s_, thermal_boundary);
            }

            HaloExchange::broadcastOwnedToGhosts(
                s_.unkno,
                s_.partition_metadata,
                MPI_COMM_WORLD);

            HaloExchange::broadcastOwnedToGhosts(
                s_.pres,
                s_.partition_metadata,
                MPI_COMM_WORLD);

            if (s_.cfg.temp_calc > 0)
            {
                HaloExchange::broadcastOwnedToGhosts(
                    s_.temperature,
                    s_.partition_metadata,
                    MPI_COMM_WORLD);
            }

            updateVelocityMagnitude();
            pressure_system.initialise(s_);
        }
'''

    text = replace_once(
        text,
        old_setup,
        new_setup,
        "production setup suppression",
    )

    old_ready = r'''        if (s_.mpi_rank == 0)
        {
            std::cout
                << "============================================================\n"
                << "CBS3D DD-4B/DD-4C PRODUCTION LOOP READY\n"
                << "============================================================\n"
                << "case                       : "
                << distributed_case_name << "\n"
                << "MPI ranks                  : " << s_.mpi_size << "\n"
                << "iterations requested       : " << s_.cfg.ntime << "\n"
                << "energy equation            : "
                << (s_.cfg.temp_calc > 0 ? "ON" : "OFF") << "\n"
                << "pressure matrix builds     : 1\n"
                << "AMG hierarchy builds       : 1\n"
                << "distributed residual CSV   : "
                << (s_.cfg.residual_log_enabled > 0 ? "ON" : "OFF") << "\n"
                << "PVTU/PVD output            : "
                << (s_.cfg.vtu_output_enabled > 0 ? "ON" : "OFF") << "\n"
                << "maximum setup time [s]     : "
                << global_setup_seconds << "\n"
                << "============================================================\n";
        }
'''

    new_ready = r'''        if (s_.mpi_rank == 0)
        {
            std::cout
                << "Run configuration\n"
                << "  mesh       : "
                << s_.partition_metadata.global_nelem << " tetrahedra, "
                << s_.partition_metadata.global_npoin << " nodes\n"
                << "  parallel   : " << s_.mpi_size << " MPI ranks\n"
                << "  iterations : " << s_.cfg.ntime << "\n"
                << "  energy     : "
                << (s_.cfg.temp_calc > 0 ? "ON" : "OFF") << "\n"
                << "  residuals  : "
                << (s_.cfg.residual_log_enabled > 0 ? "ON" : "OFF") << "\n"
                << "  VTU/PVD    : "
                << (s_.cfg.vtu_output_enabled > 0 ? "ON" : "OFF") << "\n"
                << "  setup time : " << global_setup_seconds << " s\n";
        }
'''

    text = replace_once(
        text,
        old_ready,
        new_ready,
        "production concise ready summary",
    )

    text = replace_once(
        text,
        '                    << "DD iteration " << iteration\n',
        '                    << "Iteration " << iteration\n',
        "production iteration label",
    )

    old_complete = r'''        if (s_.mpi_rank == 0)
        {
            std::cout
                << "============================================================\n"
                << "CBS3D DD-4B/DD-4C PRODUCTION LOOP COMPLETE\n"
                << "============================================================\n"
                << "MPI ranks                  : " << s_.mpi_size << "\n"
                << "completed iterations       : " << last_iteration << "\n"
                << "stop reason                : " << stop_reason << "\n"
                << "final physical time        : " << s_.cfg.rtime << "\n"
                << "pressure matrix builds     : 1\n"
                << "AMG hierarchy builds       : 1\n"
                << "distributed output case    : "
                << distributed_case_name << "\n"
                << "maximum loop time [s]      : "
                << global_loop_seconds << "\n"
                << "next development stage     : rank-count validation\n"
                << "============================================================\n";
        }
'''

    new_complete = r'''        if (s_.mpi_rank == 0)
        {
            std::cout
                << "Run complete\n"
                << "  iterations    : " << last_iteration << "\n"
                << "  stop reason   : " << stop_reason << "\n"
                << "  physical time : " << s_.cfg.rtime << "\n"
                << "  loop time     : " << global_loop_seconds << " s\n"
                << "  output case   : " << distributed_case_name << "\n";
        }
'''

    text = replace_once(
        text,
        old_complete,
        new_complete,
        "production concise completion summary",
    )

    return text, True


def main():
    if not MAIN_PATH.is_file() or not PRODUCTION_PATH.is_file():
        raise RuntimeError(
            "Run this tool from the CBS3D repository root"
        )

    main_text = MAIN_PATH.read_text()
    production_text = PRODUCTION_PATH.read_text()

    main_text, main_changed = patch_main(main_text)
    production_text, production_changed = patch_production_loop(production_text)

    if main_changed:
        MAIN_PATH.write_text(main_text)

    if production_changed:
        PRODUCTION_PATH.write_text(production_text)

    if not main_changed and not production_changed:
        print("Production-console cleanup: ALREADY APPLIED")
        return 0

    print("Production-console cleanup: APPLIED")
    print("modified: {}".format(MAIN_PATH))
    print("modified: {}".format(PRODUCTION_PATH))
    print("default stdout: rank-zero concise production output")
    print("debug opt-in: CBS3D_VERBOSE=1")
    print("all-rank opt-in: CBS3D_ALL_RANK_OUTPUT=1")
    print("numerical checks and classic VTU/PVTU/PVD output: unchanged")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print("ERROR: {}".format(error), file=sys.stderr)
        raise SystemExit(2)
