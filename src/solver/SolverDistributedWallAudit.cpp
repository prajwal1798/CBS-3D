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
        void check_wall_audit_mpi(
            const int error_code,
            const char* operation)
        {
            if (error_code != MPI_SUCCESS)
            {
                throw std::runtime_error(
                    std::string("Distributed wall audit MPI failure in ")
                    + operation);
            }
        }

        bool wall_audit_nearly_equal(
            const Real a,
            const Real b)
        {
            const Real scale =
                std::fmax(
                    1.0,
                    std::fmax(std::fabs(a), std::fabs(b)));

            return std::fabs(a - b) <= 1.0e-10 * scale;
        }
    }
#endif

    void Solver::auditDistributedWallClassification() const
    {
#ifdef CBS3D_USE_MPI
        const Int global_npoin =
            static_cast<Int>(s_.partition_metadata.global_npoin);

        if (global_npoin < 1)
        {
            throw std::runtime_error(
                "Distributed wall audit found invalid global node count");
        }

        const Int physical_bit =
            CBSStateSI::node_on_physical_wall;

        const Int interface_bit =
            CBSStateSI::node_on_material_interface;

        const Int maximum_valid_mask =
            physical_bit | interface_bit;

        const Size global_storage =
            static_cast<Size>(global_npoin) + 1U;

        std::vector<Int> local_reference_mask(
            global_storage,
            0);

        std::vector<Int> global_reference_mask(
            global_storage,
            0);

        std::vector<Real> local_reference_normal(
            3U * global_storage,
            0.0);

        std::vector<Real> global_reference_normal(
            3U * global_storage,
            0.0);

        const auto global_node_id =
            [this, global_npoin](const Int local_node) -> Int
            {
                if (local_node < 1 ||
                    local_node > s_.cfg.npoin)
                {
                    throw std::runtime_error(
                        "Distributed wall audit found local node out of range");
                }

                const Size local_index =
                    static_cast<Size>(local_node);

                if (local_index >=
                    s_.local_to_global_node.size())
                {
                    throw std::runtime_error(
                        "Distributed wall audit found incomplete node map");
                }

                const Int global_node =
                    s_.local_to_global_node[local_index];

                if (global_node < 1 ||
                    global_node > global_npoin)
                {
                    throw std::runtime_error(
                        "Distributed wall audit found global node out of range");
                }

                return global_node;
            };

        const auto is_physical_no_slip_boundary =
            [this](const Int bc) -> bool
            {
                return
                    bc == s_.cfg.bc_temperature_one_noslip ||
                    bc == s_.cfg.bc_temperature_zero_noslip ||
                    bc == s_.cfg.bc_noslip_adiabatic_wall ||
                    bc == s_.cfg.bc_noslip_heatflux_wall ||
                    bc == s_.cfg.bc_cht_interface;
            };

        // Independent physical-wall reference from rank-local physical
        // boundary faces. Artificial partition interfaces are absent.
        for (Int ib = 1; ib <= s_.cfg.nboun; ++ib)
        {
            const Int bc =
                s_.iside(s_.cfg.bsid, ib);

            if (!is_physical_no_slip_boundary(bc))
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

                local_reference_mask[
                    static_cast<Size>(global_node)] |= physical_bit;

                const Size base =
                    3U * static_cast<Size>(global_node);

                local_reference_normal[base] +=
                    s_.face_norm(1, ib);

                local_reference_normal[base + 1U] +=
                    s_.face_norm(2, ib);

                local_reference_normal[base + 2U] +=
                    s_.face_norm(3, ib);
            }
        }

        const Int material_interface_mask =
            CBSStateSI::node_touches_fluid |
            CBSStateSI::node_touches_solid;

        // Independent interface reference. Duplicate local copies are harmless
        // because the global reduction uses bitwise OR.
        for (Int ip = 1; ip <= s_.cfg.npoin; ++ip)
        {
            if (s_.node_material_mask(ip) !=
                material_interface_mask)
            {
                continue;
            }

            const Int global_node =
                global_node_id(ip);

            local_reference_mask[
                static_cast<Size>(global_node)] |= interface_bit;
        }

        check_wall_audit_mpi(
            MPI_Allreduce(
                local_reference_mask.data(),
                global_reference_mask.data(),
                global_npoin + 1,
                MPI_INT,
                MPI_BOR,
                MPI_COMM_WORLD),
            "MPI_Allreduce wall-mask reference");

        check_wall_audit_mpi(
            MPI_Allreduce(
                local_reference_normal.data(),
                global_reference_normal.data(),
                static_cast<int>(
                    global_reference_normal.size()),
                MPI_DOUBLE,
                MPI_SUM,
                MPI_COMM_WORLD),
            "MPI_Allreduce wall-normal reference");

        Int local_mask_failures = 0;
        Int local_normal_failures = 0;

        for (Int ip = 1; ip <= s_.cfg.npoin; ++ip)
        {
            const Int global_node =
                global_node_id(ip);

            const Size global_index =
                static_cast<Size>(global_node);

            if (s_.node_wall_mask(ip) !=
                global_reference_mask[global_index])
            {
                ++local_mask_failures;
            }

            const Size base =
                3U * global_index;

            for (Int dim = 1; dim <= 3; ++dim)
            {
                const Real expected =
                    global_reference_normal[
                        base + static_cast<Size>(dim - 1)];

                if (!wall_audit_nearly_equal(
                        s_.node_wall_normal_sum(dim, ip),
                        expected))
                {
                    ++local_normal_failures;
                }
            }
        }

        Int global_mask_failures = 0;
        Int global_normal_failures = 0;

        check_wall_audit_mpi(
            MPI_Allreduce(
                &local_mask_failures,
                &global_mask_failures,
                1,
                MPI_INT,
                MPI_SUM,
                MPI_COMM_WORLD),
            "MPI_Allreduce wall-mask failures");

        check_wall_audit_mpi(
            MPI_Allreduce(
                &local_normal_failures,
                &global_normal_failures,
                1,
                MPI_INT,
                MPI_SUM,
                MPI_COMM_WORLD),
            "MPI_Allreduce wall-normal failures");

        Int local_owned_counts[4] = {0, 0, 0, 0};
        Int local_invalid_masks = 0;

        for (const Int ip : s_.owned_nodes)
        {
            const Int mask =
                s_.node_wall_mask(ip);

            if (mask >= 0 &&
                mask <= maximum_valid_mask)
            {
                ++local_owned_counts[mask];
            }
            else
            {
                ++local_invalid_masks;
            }
        }

        Int global_owned_counts[4] = {0, 0, 0, 0};
        Int global_invalid_masks = 0;

        check_wall_audit_mpi(
            MPI_Allreduce(
                local_owned_counts,
                global_owned_counts,
                4,
                MPI_INT,
                MPI_SUM,
                MPI_COMM_WORLD),
            "MPI_Allreduce owned wall-mask counts");

        check_wall_audit_mpi(
            MPI_Allreduce(
                &local_invalid_masks,
                &global_invalid_masks,
                1,
                MPI_INT,
                MPI_SUM,
                MPI_COMM_WORLD),
            "MPI_Allreduce invalid wall masks");

        Int reference_counts[4] = {0, 0, 0, 0};

        for (Int global_node = 1;
             global_node <= global_npoin;
             ++global_node)
        {
            const Int mask =
                global_reference_mask[
                    static_cast<Size>(global_node)];

            if (mask < 0 ||
                mask > maximum_valid_mask)
            {
                throw std::runtime_error(
                    "Distributed wall audit produced invalid reference mask");
            }

            ++reference_counts[mask];
        }

        Int local_list_failures = 0;

        if (s_.cfg.npoin_wall < 0 ||
            s_.cfg.npoin_wall > s_.cfg.npoin)
        {
            ++local_list_failures;
        }

        Int expected_local_wall_nodes = 0;

        for (Int ip = 1; ip <= s_.cfg.npoin; ++ip)
        {
            if (s_.node_wall_mask(ip) != 0)
            {
                ++expected_local_wall_nodes;
            }
        }

        if (s_.cfg.npoin_wall !=
            expected_local_wall_nodes)
        {
            ++local_list_failures;
        }

        std::vector<Int> listed(
            static_cast<Size>(s_.cfg.npoin) + 1U,
            0);

        for (Int iw = 1;
             iw <= s_.cfg.npoin_wall;
             ++iw)
        {
            const Int ip =
                s_.wall_node_list(iw);

            if (ip < 1 || ip > s_.cfg.npoin)
            {
                ++local_list_failures;
                continue;
            }

            if (listed[static_cast<Size>(ip)] != 0 ||
                s_.node_wall_mask(ip) == 0)
            {
                ++local_list_failures;
                continue;
            }

            listed[static_cast<Size>(ip)] = 1;

            const Real nx =
                s_.node_wall_normal_sum(1, ip);

            const Real ny =
                s_.node_wall_normal_sum(2, ip);

            const Real nz =
                s_.node_wall_normal_sum(3, ip);

            const Real length =
                std::sqrt(nx * nx + ny * ny + nz * nz);

            const bool physical_wall =
                (s_.node_wall_mask(ip) & physical_bit) != 0;

            const Real expected_x =
                physical_wall && length > 0.0
                    ? nx / length
                    : 0.0;

            const Real expected_y =
                physical_wall && length > 0.0
                    ? ny / length
                    : 0.0;

            const Real expected_z =
                physical_wall && length > 0.0
                    ? nz / length
                    : 0.0;

            if (!wall_audit_nearly_equal(
                    s_.wall_node_norm(1, iw),
                    expected_x) ||
                !wall_audit_nearly_equal(
                    s_.wall_node_norm(2, iw),
                    expected_y) ||
                !wall_audit_nearly_equal(
                    s_.wall_node_norm(3, iw),
                    expected_z))
            {
                ++local_list_failures;
            }
        }

        for (Int ip = 1; ip <= s_.cfg.npoin; ++ip)
        {
            if (s_.node_wall_mask(ip) != 0 &&
                listed[static_cast<Size>(ip)] == 0)
            {
                ++local_list_failures;
            }
        }

        Int global_list_failures = 0;

        check_wall_audit_mpi(
            MPI_Allreduce(
                &local_list_failures,
                &global_list_failures,
                1,
                MPI_INT,
                MPI_SUM,
                MPI_COMM_WORLD),
            "MPI_Allreduce wall-list failures");

        const bool counts_match =
            global_owned_counts[0] == reference_counts[0] &&
            global_owned_counts[1] == reference_counts[1] &&
            global_owned_counts[2] == reference_counts[2] &&
            global_owned_counts[3] == reference_counts[3];

        if (global_mask_failures != 0 ||
            global_normal_failures != 0 ||
            global_list_failures != 0 ||
            global_invalid_masks != 0 ||
            !counts_match)
        {
            throw std::runtime_error(
                "Distributed wall/interface reconciliation failed");
        }

        const Int physical_only =
            global_owned_counts[physical_bit];

        const Int interface_only =
            global_owned_counts[interface_bit];

        const Int overlap =
            global_owned_counts[maximum_valid_mask];

        const Int physical_total =
            physical_only + overlap;

        const Int interface_total =
            interface_only + overlap;

        const Int wall_union =
            physical_only + interface_only + overlap;

        if (s_.mpi_rank == 0)
        {
            std::cout
                << "============================================================\n"
                << "CBS3D DISTRIBUTED WALL CLASSIFICATION\n"
                << "============================================================\n"
                << "MPI ranks                    : "
                << s_.mpi_size << "\n"
                << "physical-only owned nodes    : "
                << physical_only << "\n"
                << "interface-only owned nodes   : "
                << interface_only << "\n"
                << "physical/interface overlap   : "
                << overlap << "\n"
                << "physical wall owned nodes    : "
                << physical_total << "\n"
                << "material-interface nodes     : "
                << interface_total << "\n"
                << "wall/interface union nodes   : "
                << wall_union << "\n"
                << "invalid wall masks           : "
                << global_invalid_masks << "\n"
                << "owner/ghost wall-mask agreement  : PASS\n"
                << "owner/ghost wall-normal agreement: PASS\n"
                << "local wall-list consistency      : PASS\n"
                << "independent global reference     : PASS\n"
                << "CBS Steps 1 to 4                 : NOT STARTED\n"
                << "RESULT                            : PASS\n"
                << "============================================================\n";
        }
#else
        throw std::runtime_error(
            "Solver::auditDistributedWallClassification "
            "requires an MPI build");
#endif
    }
}
