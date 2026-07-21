#include "cbs/solver/Solver.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef CBS3D_USE_MPI
#include <mpi.h>
#endif

namespace cbs
{
#ifdef CBS3D_USE_MPI
    namespace
    {
        void check_pressure_audit_mpi(
            const int error_code,
            const char* operation)
        {
            if (error_code != MPI_SUCCESS)
            {
                throw std::runtime_error(
                    std::string(
                        "Distributed pressure audit MPI failure in ")
                    + operation);
            }
        }


        bool pressure_audit_nearly_equal(
            const Real a,
            const Real b)
        {
            const Real scale =
                std::fmax(
                    1.0,
                    std::fmax(std::fabs(a), std::fabs(b)));

            return std::fabs(a - b) <= 1.0e-12 * scale;
        }
    }
#endif


    void Solver::auditDistributedPressureBoundary() const
    {
#ifdef CBS3D_USE_MPI
        const Int global_npoin =
            static_cast<Int>(
                s_.partition_metadata.global_npoin);

        if (global_npoin < 1)
        {
            throw std::runtime_error(
                "Distributed pressure audit found invalid "
                "global node count");
        }

        const Size global_storage =
            static_cast<Size>(global_npoin) + 1U;

        // Independent global references reconstructed without using
        // node_pressure_fixed.
        std::vector<Int> local_outlet_reference(
            global_storage,
            0);

        std::vector<Int> global_outlet_reference(
            global_storage,
            0);

        std::vector<Int> local_fluid_reference(
            global_storage,
            0);

        std::vector<Int> global_fluid_reference(
            global_storage,
            0);


        const auto global_node_id =
            [this, global_npoin](const Int local_node) -> Int
            {
                if (local_node < 1 ||
                    local_node > s_.cfg.npoin)
                {
                    throw std::runtime_error(
                        "Distributed pressure audit found "
                        "local node out of range");
                }

                const Size local_index =
                    static_cast<Size>(local_node);

                if (local_index >=
                    s_.local_to_global_node.size())
                {
                    throw std::runtime_error(
                        "Distributed pressure audit found "
                        "incomplete local-to-global node map");
                }

                const Int global_node =
                    s_.local_to_global_node[local_index];

                if (global_node < 1 ||
                    global_node > global_npoin)
                {
                    throw std::runtime_error(
                        "Distributed pressure audit found "
                        "global node out of range");
                }

                return global_node;
            };


        // Independently reconstruct fluid connectivity from owned fluid
        // tetrahedra rather than consuming node_material_mask.
        for (Int ie = 1; ie <= s_.cfg.nelem; ++ie)
        {
            if (s_.mat_elem(ie) != 0)
            {
                continue;
            }

            for (Int in = 1; in <= s_.cfg.nep; ++in)
            {
                const Int local_node =
                    s_.intma(in, ie);

                const Int global_node =
                    global_node_id(local_node);

                local_fluid_reference[
                    static_cast<Size>(global_node)] = 1;
            }
        }


        // Independently reconstruct explicit prescribed-pressure state from
        // this rank's physical BC 520 boundary faces.
        for (Int ib = 1; ib <= s_.cfg.nboun; ++ib)
        {
            const Int bc =
                s_.iside(s_.cfg.bsid, ib);

            if (bc != s_.cfg.bc_pressure_outlet)
            {
                continue;
            }

            for (Int in = 1;
                 in <= s_.cfg.nsidp;
                 ++in)
            {
                const Int local_node =
                    s_.iside(in, ib);

                const Int global_node =
                    global_node_id(local_node);

                local_outlet_reference[
                    static_cast<Size>(global_node)] = 1;
            }
        }


        check_pressure_audit_mpi(
            MPI_Allreduce(
                local_fluid_reference.data(),
                global_fluid_reference.data(),
                global_npoin + 1,
                MPI_INT,
                MPI_BOR,
                MPI_COMM_WORLD),
            "MPI_Allreduce fluid-node reference");


        check_pressure_audit_mpi(
            MPI_Allreduce(
                local_outlet_reference.data(),
                global_outlet_reference.data(),
                global_npoin + 1,
                MPI_INT,
                MPI_BOR,
                MPI_COMM_WORLD),
            "MPI_Allreduce pressure-outlet reference");


        Int explicit_outlet_nodes = 0;
        Int nonfluid_outlet_nodes = 0;

        for (Int global_node = 1;
             global_node <= global_npoin;
             ++global_node)
        {
            const Size global_index =
                static_cast<Size>(global_node);

            if (global_outlet_reference[global_index] == 0)
            {
                continue;
            }

            ++explicit_outlet_nodes;

            if (global_fluid_reference[global_index] == 0)
            {
                ++nonfluid_outlet_nodes;
            }
        }

        if (nonfluid_outlet_nodes != 0)
        {
            throw std::runtime_error(
                "Distributed pressure audit found a BC 520 node "
                "without fluid connectivity");
        }


        std::vector<Int> global_pressure_reference =
            global_outlet_reference;

        bool fallback_used = false;
        Int fallback_global_node = 0;


        // Independently reproduce the deterministic fallback rule only when
        // the global mesh has no explicit BC 520 pressure outlet.
        if (explicit_outlet_nodes < 1)
        {
            fallback_used = true;

            const Int local_requested_node =
                s_.cfg.pnode;

            Int requested_min = 0;
            Int requested_max = 0;

            check_pressure_audit_mpi(
                MPI_Allreduce(
                    &local_requested_node,
                    &requested_min,
                    1,
                    MPI_INT,
                    MPI_MIN,
                    MPI_COMM_WORLD),
                "MPI_Allreduce requested pnode minimum");

            check_pressure_audit_mpi(
                MPI_Allreduce(
                    &local_requested_node,
                    &requested_max,
                    1,
                    MPI_INT,
                    MPI_MAX,
                    MPI_COMM_WORLD),
                "MPI_Allreduce requested pnode maximum");

            if (requested_min != requested_max)
            {
                throw std::runtime_error(
                    "Distributed pressure audit found inconsistent "
                    "pnode values across MPI ranks");
            }

            if (requested_min >= 1 &&
                requested_min <= global_npoin &&
                global_fluid_reference[
                    static_cast<Size>(requested_min)] != 0)
            {
                fallback_global_node =
                    requested_min;
            }
            else
            {
                for (Int global_node = 1;
                     global_node <= global_npoin;
                     ++global_node)
                {
                    if (global_fluid_reference[
                            static_cast<Size>(global_node)] != 0)
                    {
                        fallback_global_node =
                            global_node;
                        break;
                    }
                }
            }

            if (fallback_global_node < 1 ||
                fallback_global_node > global_npoin)
            {
                throw std::runtime_error(
                    "Distributed pressure audit found no "
                    "fluid-connected fallback node");
            }

            global_pressure_reference[
                static_cast<Size>(fallback_global_node)] = 1;
        }


        Int reference_fixed_nodes = 0;

        for (Int global_node = 1;
             global_node <= global_npoin;
             ++global_node)
        {
            if (global_pressure_reference[
                    static_cast<Size>(global_node)] != 0)
            {
                ++reference_fixed_nodes;
            }
        }


        // Compare every owner and ghost copy against the independently
        // reconstructed global pressure state.
        Int local_flag_failures = 0;
        Int local_invalid_copies = 0;

        for (Int ip = 1; ip <= s_.cfg.npoin; ++ip)
        {
            const Int flag =
                s_.node_pressure_fixed(ip);

            if (flag != 0 && flag != 1)
            {
                ++local_invalid_copies;
                continue;
            }

            const Int global_node =
                global_node_id(ip);

            const Int expected =
                global_pressure_reference[
                    static_cast<Size>(global_node)];

            if (flag != expected)
            {
                ++local_flag_failures;
            }
        }


        Int global_flag_failures = 0;
        Int global_invalid_copies = 0;

        check_pressure_audit_mpi(
            MPI_Allreduce(
                &local_flag_failures,
                &global_flag_failures,
                1,
                MPI_INT,
                MPI_SUM,
                MPI_COMM_WORLD),
            "MPI_Allreduce pressure-flag failures");

        check_pressure_audit_mpi(
            MPI_Allreduce(
                &local_invalid_copies,
                &global_invalid_copies,
                1,
                MPI_INT,
                MPI_SUM,
                MPI_COMM_WORLD),
            "MPI_Allreduce invalid pressure copies");


        // Count only owner copies for the unique global pressure-node count.
        Int local_owned_fixed_nodes = 0;
        Int local_owned_invalid_flags = 0;

        for (const Int ip : s_.owned_nodes)
        {
            const Int flag =
                s_.node_pressure_fixed(ip);

            if (flag == 1)
            {
                ++local_owned_fixed_nodes;
            }
            else if (flag != 0)
            {
                ++local_owned_invalid_flags;
            }
        }


        Int global_owned_fixed_nodes = 0;
        Int global_owned_invalid_flags = 0;

        check_pressure_audit_mpi(
            MPI_Allreduce(
                &local_owned_fixed_nodes,
                &global_owned_fixed_nodes,
                1,
                MPI_INT,
                MPI_SUM,
                MPI_COMM_WORLD),
            "MPI_Allreduce owned pressure-fixed count");

        check_pressure_audit_mpi(
            MPI_Allreduce(
                &local_owned_invalid_flags,
                &global_owned_invalid_flags,
                1,
                MPI_INT,
                MPI_SUM,
                MPI_COMM_WORLD),
            "MPI_Allreduce owned invalid pressure flags");


        // Validate the rank-local bc_list and bc_values rebuilt from the
        // reconciled node_pressure_fixed array.
        Int local_list_failures = 0;
        Int expected_local_fixed_nodes = 0;

        for (Int ip = 1; ip <= s_.cfg.npoin; ++ip)
        {
            if (s_.node_pressure_fixed(ip) == 1)
            {
                ++expected_local_fixed_nodes;
            }
        }

        Int list_limit =
            s_.cfg.bc_fixed;

        if (list_limit < 0 ||
            list_limit > s_.cfg.npoin)
        {
            ++local_list_failures;
            list_limit = 0;
        }

        if (s_.cfg.bc_fixed !=
            expected_local_fixed_nodes)
        {
            ++local_list_failures;
        }

        std::vector<Int> listed(
            static_cast<Size>(s_.cfg.npoin) + 1U,
            0);

        for (Int ibc = 1;
             ibc <= list_limit;
             ++ibc)
        {
            const Int ip =
                s_.bc_list(ibc);

            if (ip < 1 ||
                ip > s_.cfg.npoin)
            {
                ++local_list_failures;
                continue;
            }

            const Size local_index =
                static_cast<Size>(ip);

            if (listed[local_index] != 0 ||
                s_.node_pressure_fixed(ip) != 1)
            {
                ++local_list_failures;
                continue;
            }

            listed[local_index] = 1;

            const Real prescribed_pressure =
                s_.bc_values(ibc);

            if (!std::isfinite(prescribed_pressure) ||
                !pressure_audit_nearly_equal(
                    prescribed_pressure,
                    s_.cfg.outlet_pressure_gauge))
            {
                ++local_list_failures;
            }
        }

        for (Int ip = 1; ip <= s_.cfg.npoin; ++ip)
        {
            if (s_.node_pressure_fixed(ip) == 1 &&
                listed[static_cast<Size>(ip)] == 0)
            {
                ++local_list_failures;
            }
        }


        Int global_list_failures = 0;

        check_pressure_audit_mpi(
            MPI_Allreduce(
                &local_list_failures,
                &global_list_failures,
                1,
                MPI_INT,
                MPI_SUM,
                MPI_COMM_WORLD),
            "MPI_Allreduce pressure-list failures");


        const bool reference_count_matches =
            global_owned_fixed_nodes ==
            reference_fixed_nodes;

        if (global_flag_failures != 0 ||
            global_invalid_copies != 0 ||
            global_owned_invalid_flags != 0 ||
            global_list_failures != 0 ||
            !reference_count_matches)
        {
            throw std::runtime_error(
                "Distributed pressure-boundary reconciliation failed");
        }


        if (s_.mpi_rank == 0)
        {
            std::cout
                << "============================================================\n"
                << "CBS3D DISTRIBUTED PRESSURE BOUNDARY\n"
                << "============================================================\n"
                << "MPI ranks                         : "
                << s_.mpi_size << "\n"
                << "explicit pressure-outlet nodes    : "
                << explicit_outlet_nodes << "\n"
                << "pressure-fixed owned nodes         : "
                << global_owned_fixed_nodes << "\n"
                << "invalid pressure flags             : "
                << global_owned_invalid_flags << "\n"
                << "owner/ghost pressure-flag agreement: PASS\n"
                << "local pressure-list consistency    : PASS\n"
                << "independent global reference       : PASS\n";

            if (fallback_used)
            {
                std::cout
                    << "pressure reference fallback        : "
                    << fallback_global_node << "\n";
            }
            else
            {
                std::cout
                    << "pressure reference fallback        : "
                    << "NOT USED\n";
            }

            std::cout
                << "CBS Steps 1 to 4                    : NOT STARTED\n"
                << "RESULT                               : PASS\n"
                << "============================================================\n";
        }
#else
        throw std::runtime_error(
            "Solver::auditDistributedPressureBoundary "
            "requires an MPI build");
#endif
    }
}
