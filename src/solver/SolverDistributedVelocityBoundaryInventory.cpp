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
        enum VelocityBoundaryInventoryBit : Int
        {
            inventory_bc_500 = 1 << 0,
            inventory_bc_501_502 = 1 << 1,
            inventory_bc_503 = 1 << 2,
            inventory_bc_506 = 1 << 3,
            inventory_bc_507 = 1 << 4,
            inventory_bc_508 = 1 << 5,
            inventory_bc_510 = 1 << 6,
            inventory_bc_511 = 1 << 7,
            inventory_bc_520 = 1 << 8,
            inventory_bc_530_532_901 = 1 << 9,
            inventory_bc_504 = 1 << 10,
            inventory_bc_902 = 1 << 11,
            inventory_material_solid = 1 << 12
        };


        void check_velocity_inventory_mpi(
            const int error_code,
            const char* operation)
        {
            if (error_code != MPI_SUCCESS)
            {
                throw std::runtime_error(
                    std::string(
                        "Distributed velocity-boundary inventory MPI failure in ")
                    + operation);
            }
        }


        Int count_set_bits(Int value)
        {
            Int count = 0;

            while (value != 0)
            {
                count += value & 1;
                value >>= 1;
            }

            return count;
        }


        bool finite_face_normal(
            const CBSStateSI& s,
            const Int boundary_face)
        {
            return
                std::isfinite(s.face_norm(1, boundary_face)) &&
                std::isfinite(s.face_norm(2, boundary_face)) &&
                std::isfinite(s.face_norm(3, boundary_face)) &&
                std::isfinite(s.face_norm(4, boundary_face)) &&
                s.face_norm(4, boundary_face) > 0.0;
        }


        Int inventory_bit_for_boundary(
            const CBSStateSI& s,
            const Int bc)
        {
            if (bc == s.cfg.bc_adiabatic_prescribed_velocity)
            {
                return inventory_bc_500;
            }

            if (bc == s.cfg.bc_temperature_one_noslip ||
                bc == s.cfg.bc_temperature_zero_noslip)
            {
                return inventory_bc_501_502;
            }

            if (bc == s.cfg.bc_temperature_zero_prescribed_velocity)
            {
                return inventory_bc_503;
            }

            if (bc == s.cfg.bc_pressure)
            {
                return inventory_bc_504;
            }

            if (bc == s.cfg.bc_symmetry_no_flux)
            {
                return inventory_bc_506;
            }

            if (bc == s.cfg.bc_bfs_parabolic_inlet)
            {
                return inventory_bc_507;
            }

            if (bc == s.cfg.bc_parabolic_inlet)
            {
                return inventory_bc_508;
            }

            if (bc == s.cfg.bc_velocity_temperature_inlet)
            {
                return inventory_bc_510;
            }

            if (bc == s.cfg.bc_massflow_temperature_inlet)
            {
                return inventory_bc_511;
            }

            if (bc == s.cfg.bc_pressure_outlet)
            {
                return inventory_bc_520;
            }

            if (bc == s.cfg.bc_noslip_adiabatic_wall ||
                bc == s.cfg.bc_noslip_heatflux_wall ||
                bc == s.cfg.bc_cht_interface)
            {
                return inventory_bc_530_532_901;
            }

            if (bc == s.cfg.bc_heatflux_marker)
            {
                return inventory_bc_902;
            }

            throw std::runtime_error(
                "Distributed velocity-boundary inventory found "
                "an unsupported boundary identifier");
        }


        bool mask_has(
            const Int mask,
            const Int bit)
        {
            return (mask & bit) != 0;
        }
    }
#endif


    void Solver::auditDistributedVelocityBoundaryInventory() const
    {
#ifdef CBS3D_USE_MPI
        const Int global_npoin =
            static_cast<Int>(
                s_.partition_metadata.global_npoin);

        if (global_npoin < 1)
        {
            throw std::runtime_error(
                "Distributed velocity-boundary inventory found "
                "an invalid global node count");
        }

        const Size global_storage =
            static_cast<Size>(global_npoin) + 1U;

        std::vector<Int> local_boundary_mask(
            global_storage,
            0);

        std::vector<Int> global_boundary_mask(
            global_storage,
            0);

        std::vector<Int> local_inlet_face_count(
            global_storage,
            0);

        std::vector<Int> global_inlet_face_count(
            global_storage,
            0);

        // Three consecutive entries are stored for each global node:
        //
        //     [3*i + 0]  sum(A_f n_x)
        //     [3*i + 1]  sum(A_f n_y)
        //     [3*i + 2]  sum(A_f n_z)
        std::vector<Real> local_inlet_normal_sum(
            3U * global_storage,
            0.0);

        std::vector<Real> global_inlet_normal_sum(
            3U * global_storage,
            0.0);


        const auto global_node_id =
            [this, global_npoin](const Int local_node) -> Int
            {
                if (local_node < 1 ||
                    local_node > s_.cfg.npoin)
                {
                    throw std::runtime_error(
                        "Distributed velocity-boundary inventory found "
                        "a local node outside the rank-local range");
                }

                const Size local_index =
                    static_cast<Size>(local_node);

                if (local_index >=
                    s_.local_to_global_node.size())
                {
                    throw std::runtime_error(
                        "Distributed velocity-boundary inventory found "
                        "an incomplete local-to-global node map");
                }

                const Int global_node =
                    s_.local_to_global_node[local_index];

                if (global_node < 1 ||
                    global_node > global_npoin)
                {
                    throw std::runtime_error(
                        "Distributed velocity-boundary inventory found "
                        "a global node outside the mesh range");
                }

                return global_node;
            };


        // -------------------------------------------------------------
        // Build the rank-local contribution directly from physical
        // boundary faces. Artificial MPI partition faces are absent from
        // the boundary-face files and therefore cannot enter this audit.
        // -------------------------------------------------------------
        for (Int ib = 1; ib <= s_.cfg.nboun; ++ib)
        {
            const Int bc =
                s_.iside(s_.cfg.bsid, ib);

            const Int boundary_bit =
                inventory_bit_for_boundary(s_, bc);

            if (!finite_face_normal(s_, ib))
            {
                throw std::runtime_error(
                    "Distributed velocity-boundary inventory found "
                    "an invalid physical boundary-face normal");
            }

            for (Int in = 1;
                 in <= s_.cfg.nsidp;
                 ++in)
            {
                const Int local_node =
                    s_.iside(in, ib);

                const Int global_node =
                    global_node_id(local_node);

                const Size global_index =
                    static_cast<Size>(global_node);

                local_boundary_mask[global_index] |=
                    boundary_bit;

                if (bc == s_.cfg.bc_massflow_temperature_inlet)
                {
                    ++local_inlet_face_count[global_index];

                    const Size vector_index =
                        3U * global_index;

                    local_inlet_normal_sum[vector_index + 0U] +=
                        s_.face_norm(1, ib);

                    local_inlet_normal_sum[vector_index + 1U] +=
                        s_.face_norm(2, ib);

                    local_inlet_normal_sum[vector_index + 2U] +=
                        s_.face_norm(3, ib);
                }
            }
        }


        // -------------------------------------------------------------
        // Add material-solid adjacency as a separate constraint. The
        // material mask has already been independently validated by DD-2.
        // -------------------------------------------------------------
        for (Int ip = 1; ip <= s_.cfg.npoin; ++ip)
        {
            const Int material_mask =
                s_.node_material_mask(ip);

            const Int valid_material_mask =
                CBSStateSI::node_touches_fluid |
                CBSStateSI::node_touches_solid;

            if (material_mask < CBSStateSI::node_touches_fluid ||
                material_mask > valid_material_mask)
            {
                throw std::runtime_error(
                    "Distributed velocity-boundary inventory found "
                    "an invalid material node mask");
            }

            if ((material_mask &
                 CBSStateSI::node_touches_solid) == 0)
            {
                continue;
            }

            const Int global_node =
                global_node_id(ip);

            local_boundary_mask[
                static_cast<Size>(global_node)] |=
                inventory_material_solid;
        }


        check_velocity_inventory_mpi(
            MPI_Allreduce(
                local_boundary_mask.data(),
                global_boundary_mask.data(),
                global_npoin + 1,
                MPI_INT,
                MPI_BOR,
                MPI_COMM_WORLD),
            "MPI_Allreduce global boundary masks");


        check_velocity_inventory_mpi(
            MPI_Allreduce(
                local_inlet_face_count.data(),
                global_inlet_face_count.data(),
                global_npoin + 1,
                MPI_INT,
                MPI_SUM,
                MPI_COMM_WORLD),
            "MPI_Allreduce inlet-face multiplicity");


        check_velocity_inventory_mpi(
            MPI_Allreduce(
                local_inlet_normal_sum.data(),
                global_inlet_normal_sum.data(),
                3 * (global_npoin + 1),
                MPI_DOUBLE,
                MPI_SUM,
                MPI_COMM_WORLD),
            "MPI_Allreduce inlet normal sums");


        const Int prescribed_nonzero_bits =
            inventory_bc_500 |
            inventory_bc_503 |
            inventory_bc_507 |
            inventory_bc_508 |
            inventory_bc_510 |
            inventory_bc_511;

        const Int physical_noslip_bits =
            inventory_bc_501_502 |
            inventory_bc_530_532_901;

        const Int strong_velocity_bits =
            prescribed_nonzero_bits |
            physical_noslip_bits |
            inventory_material_solid;


        Int moving_wall_nodes = 0;
        Int temperature_noslip_nodes = 0;
        Int legacy_fixed_x_nodes = 0;
        Int symmetry_nodes = 0;
        Int bfs_inlet_nodes = 0;
        Int parabolic_inlet_nodes = 0;
        Int prescribed_velocity_inlet_nodes = 0;
        Int massflow_inlet_nodes = 0;
        Int pressure_boundary_nodes = 0;
        Int pressure_outlet_nodes = 0;
        Int physical_noslip_nodes = 0;
        Int thermal_marker_nodes = 0;
        Int material_solid_nodes = 0;

        Int inlet_wall_overlap = 0;
        Int inlet_material_overlap = 0;
        Int inlet_outlet_overlap = 0;
        Int wall_outlet_overlap = 0;
        Int moving_wall_wall_overlap = 0;
        Int symmetry_strong_overlap = 0;

        Int multiple_prescribed_family_nodes = 0;
        Int multiple_strong_family_nodes = 0;
        Int multiple_inlet_face_nodes = 0;
        Int zero_resultant_inlet_normals = 0;
        Int invalid_global_masks = 0;


        for (Int global_node = 1;
             global_node <= global_npoin;
             ++global_node)
        {
            const Size global_index =
                static_cast<Size>(global_node);

            const Int mask =
                global_boundary_mask[global_index];

            if (mask < 0)
            {
                ++invalid_global_masks;
                continue;
            }

            const bool has_500 =
                mask_has(mask, inventory_bc_500);

            const bool has_501_502 =
                mask_has(mask, inventory_bc_501_502);

            const bool has_503 =
                mask_has(mask, inventory_bc_503);

            const bool has_504 =
                mask_has(mask, inventory_bc_504);

            const bool has_506 =
                mask_has(mask, inventory_bc_506);

            const bool has_507 =
                mask_has(mask, inventory_bc_507);

            const bool has_508 =
                mask_has(mask, inventory_bc_508);

            const bool has_510 =
                mask_has(mask, inventory_bc_510);

            const bool has_511 =
                mask_has(mask, inventory_bc_511);

            const bool has_520 =
                mask_has(mask, inventory_bc_520);

            const bool has_physical_noslip =
                (mask & physical_noslip_bits) != 0;

            const bool has_902 =
                mask_has(mask, inventory_bc_902);

            const bool has_material_solid =
                mask_has(mask, inventory_material_solid);

            moving_wall_nodes += has_500 ? 1 : 0;
            temperature_noslip_nodes += has_501_502 ? 1 : 0;
            legacy_fixed_x_nodes += has_503 ? 1 : 0;
            pressure_boundary_nodes += has_504 ? 1 : 0;
            symmetry_nodes += has_506 ? 1 : 0;
            bfs_inlet_nodes += has_507 ? 1 : 0;
            parabolic_inlet_nodes += has_508 ? 1 : 0;
            prescribed_velocity_inlet_nodes += has_510 ? 1 : 0;
            massflow_inlet_nodes += has_511 ? 1 : 0;
            pressure_outlet_nodes += has_520 ? 1 : 0;
            physical_noslip_nodes += has_physical_noslip ? 1 : 0;
            thermal_marker_nodes += has_902 ? 1 : 0;
            material_solid_nodes += has_material_solid ? 1 : 0;


            if (has_511 && has_physical_noslip)
            {
                ++inlet_wall_overlap;
            }

            if (has_511 && has_material_solid)
            {
                ++inlet_material_overlap;
            }

            if (has_511 && has_520)
            {
                ++inlet_outlet_overlap;
            }

            if (has_physical_noslip && has_520)
            {
                ++wall_outlet_overlap;
            }

            if (has_500 && has_physical_noslip)
            {
                ++moving_wall_wall_overlap;
            }

            if (has_506 &&
                (mask & strong_velocity_bits) != 0)
            {
                ++symmetry_strong_overlap;
            }


            const Int prescribed_family_count =
                count_set_bits(
                    mask & prescribed_nonzero_bits);

            if (prescribed_family_count > 1)
            {
                ++multiple_prescribed_family_nodes;
            }


            const Int strong_family_count =
                count_set_bits(
                    mask & strong_velocity_bits);

            if (strong_family_count > 1)
            {
                ++multiple_strong_family_nodes;
            }


            if (global_inlet_face_count[global_index] > 1)
            {
                ++multiple_inlet_face_nodes;
            }


            if (has_511)
            {
                const Size vector_index =
                    3U * global_index;

                const Real nx =
                    global_inlet_normal_sum[
                        vector_index + 0U];

                const Real ny =
                    global_inlet_normal_sum[
                        vector_index + 1U];

                const Real nz =
                    global_inlet_normal_sum[
                        vector_index + 2U];

                const Real magnitude =
                    std::sqrt(
                        nx * nx +
                        ny * ny +
                        nz * nz);

                if (!std::isfinite(magnitude) ||
                    magnitude <= 1.0e-14)
                {
                    ++zero_resultant_inlet_normals;
                }
            }
        }


        if (invalid_global_masks != 0 ||
            zero_resultant_inlet_normals != 0)
        {
            throw std::runtime_error(
                "Distributed velocity-boundary inventory failed");
        }


        if (s_.mpi_rank == 0)
        {
            std::cout
                << "============================================================\n"
                << "CBS3D DISTRIBUTED VELOCITY-BOUNDARY INVENTORY\n"
                << "============================================================\n"
                << "MPI ranks                           : "
                << s_.mpi_size << "\n"
                << "BC 500 moving/prescribed nodes      : "
                << moving_wall_nodes << "\n"
                << "BC 501/502 no-slip nodes            : "
                << temperature_noslip_nodes << "\n"
                << "BC 503 fixed-x velocity nodes       : "
                << legacy_fixed_x_nodes << "\n"
                << "BC 504 pressure-only nodes          : "
                << pressure_boundary_nodes << "\n"
                << "BC 506 symmetry nodes               : "
                << symmetry_nodes << "\n"
                << "BC 507 BFS inlet nodes              : "
                << bfs_inlet_nodes << "\n"
                << "BC 508 parabolic inlet nodes        : "
                << parabolic_inlet_nodes << "\n"
                << "BC 510 velocity-inlet nodes         : "
                << multiple_prescribed_family_nodes
                << "BC 511 mass-flow inlet nodes        : "
                << massflow_inlet_nodes << "\n"
                << "BC 520 pressure-outlet nodes        : "
                << pressure_outlet_nodes << "\n"
                << "physical no-slip nodes              : "
                << physical_noslip_nodes << "\n"
                << "BC 902 thermal-marker nodes         : "
                << thermal_marker_nodes << "\n"
                << "material-solid constrained nodes    : "
                << material_solid_nodes << "\n"
                << "BC 511 / physical-wall overlap      : "
                << inlet_wall_overlap << "\n"
                << "BC 511 / material-solid overlap     : "
                << inlet_material_overlap << "\n"
                << "BC 511 / outlet overlap             : "
                << inlet_outlet_overlap << "\n"
                << "physical-wall / outlet overlap      : "
                << wall_outlet_overlap << "\n"
                << "moving-wall / wall overlap          : "
                << moving_wall_wall_overlap << "\n"
                << "symmetry / strong-BC overlap        : "
                << symmetry_strong_overlap << "\n"
                << "multiple prescribed-velocity families: "
                << multiple_prescribed_family_nodes << "\n"
                << "multiple strong-velocity families   : "
                << multiple_strong_family_nodes << "\n"
                << "nodes touched by multiple BC 511 faces: "
                << multiple_inlet_face_nodes << "\n"
                << "zero-resultant BC 511 normals        : "
                << zero_resultant_inlet_normals << "\n"
                << "invalid global boundary masks       : "
                << invalid_global_masks << "\n"
                << "CBS Steps 1 to 4                    : NOT STARTED\n"
                << "RESULT                               : PASS\n"
                << "============================================================\n";
        }
#else
        throw std::runtime_error(
            "Solver::auditDistributedVelocityBoundaryInventory "
            "requires an MPI build");
#endif
    }
}
