#include "cbs/parallel/PartitionMetadata.hpp"

#include <fstream>
#include <stdexcept>
#include <string>

namespace cbs
{
    namespace
    {
        void requireToken(std::istream& in, const char* expected)
        {
            std::string token;
            if (!(in >> token) || token != expected)
            {
                throw std::runtime_error(
                    std::string("PartitionMetadataIO - expected token '")
                    + expected + "'");
            }
        }
    }

    PartitionMetadata PartitionMetadataIO::read(const std::string& path)
    {
        std::ifstream in(path);
        if (!in)
        {
            throw std::runtime_error(
                "PartitionMetadataIO - cannot open file: " + path);
        }

        PartitionMetadata metadata;

        requireToken(in, "CBS3D_MPI_PARTITION");
        in >> metadata.version;
        if (metadata.version != 1)
        {
            throw std::runtime_error(
                "PartitionMetadataIO - unsupported metadata version");
        }

        requireToken(in, "PARTITION");
        in >> metadata.partition_id
           >> metadata.mpi_rank
           >> metadata.mpi_size;

        requireToken(in, "GLOBAL");
        in >> metadata.global_nelem
           >> metadata.global_npoin
           >> metadata.global_nboun;

        requireToken(in, "LOCAL");
        in >> metadata.local_nelem
           >> metadata.local_npoin
           >> metadata.owned_nodes
           >> metadata.ghost_nodes
           >> metadata.local_nboun;

        requireToken(in, "NODES");
        metadata.local_to_global_node.assign(
            static_cast<Size>(metadata.local_npoin) + 1U,
            0);
        metadata.node_owner_rank.assign(
            static_cast<Size>(metadata.local_npoin) + 1U,
            -1);

        for (Int i = 1; i <= metadata.local_npoin; ++i)
        {
            Int local_id = 0;
            std::int64_t global_id = 0;
            Int owner_rank = -1;
            in >> local_id >> global_id >> owner_rank;

            if (local_id < 1 || local_id > metadata.local_npoin)
            {
                throw std::runtime_error(
                    "PartitionMetadataIO - local node ID out of range");
            }

            metadata.local_to_global_node[static_cast<Size>(local_id)] = global_id;
            metadata.node_owner_rank[static_cast<Size>(local_id)] = owner_rank;
        }

        requireToken(in, "ELEMENTS");
        metadata.local_to_global_element.assign(
            static_cast<Size>(metadata.local_nelem) + 1U,
            0);

        for (Int i = 1; i <= metadata.local_nelem; ++i)
        {
            Int local_id = 0;
            std::int64_t global_id = 0;
            in >> local_id >> global_id;

            if (local_id < 1 || local_id > metadata.local_nelem)
            {
                throw std::runtime_error(
                    "PartitionMetadataIO - local element ID out of range");
            }

            metadata.local_to_global_element[static_cast<Size>(local_id)] = global_id;
        }

        requireToken(in, "NEIGHBOURS");
        Int neighbour_count = 0;
        in >> neighbour_count;
        if (neighbour_count < 0)
        {
            throw std::runtime_error(
                "PartitionMetadataIO - negative neighbour count");
        }

        metadata.neighbours.reserve(static_cast<Size>(neighbour_count));

        for (Int i = 0; i < neighbour_count; ++i)
        {
            PartitionNeighbour neighbour;
            Int send_count = 0;
            Int recv_count = 0;

            requireToken(in, "NEIGHBOUR");
            in >> neighbour.rank >> send_count >> recv_count;

            if (send_count < 0 || recv_count < 0)
            {
                throw std::runtime_error(
                    "PartitionMetadataIO - negative halo count");
            }

            requireToken(in, "SEND");
            neighbour.send_local_nodes.resize(static_cast<Size>(send_count));
            for (Int& local_id : neighbour.send_local_nodes)
            {
                in >> local_id;
                if (local_id < 1 || local_id > metadata.local_npoin)
                {
                    throw std::runtime_error(
                        "PartitionMetadataIO - SEND node ID out of range");
                }
            }

            requireToken(in, "RECV");
            neighbour.recv_local_nodes.resize(static_cast<Size>(recv_count));
            for (Int& local_id : neighbour.recv_local_nodes)
            {
                in >> local_id;
                if (local_id < 1 || local_id > metadata.local_npoin)
                {
                    throw std::runtime_error(
                        "PartitionMetadataIO - RECV node ID out of range");
                }
            }

            metadata.neighbours.push_back(std::move(neighbour));
        }

        requireToken(in, "END");

        if (!in)
        {
            throw std::runtime_error(
                "PartitionMetadataIO - malformed metadata file: " + path);
        }

        if (metadata.owned_nodes + metadata.ghost_nodes != metadata.local_npoin)
        {
            throw std::runtime_error(
                "PartitionMetadataIO - owned+ghost node count mismatch");
        }

        return metadata;
    }
}
