#pragma once

#include "cbs/core/Types.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace cbs
{
    struct PartitionNeighbour
    {
        Int rank = -1;
        std::vector<Int> send_local_nodes;
        std::vector<Int> recv_local_nodes;
    };

    struct PartitionMetadata
    {
        Int version = 0;
        Int partition_id = 0;
        Int mpi_rank = 0;
        Int mpi_size = 1;

        std::int64_t global_nelem = 0;
        std::int64_t global_npoin = 0;
        std::int64_t global_nboun = 0;

        Int local_nelem = 0;
        Int local_npoin = 0;
        Int owned_nodes = 0;
        Int ghost_nodes = 0;
        Int local_nboun = 0;

        // One-based local node/element IDs map to original global Gmsh IDs.
        std::vector<std::int64_t> local_to_global_node;
        std::vector<std::int64_t> local_to_global_element;

        // One-based local node ID -> owning MPI rank.
        std::vector<Int> node_owner_rank;

        std::vector<PartitionNeighbour> neighbours;
    };

    class PartitionMetadataIO
    {
    public:
        static PartitionMetadata read(const std::string& path);
    };
}
