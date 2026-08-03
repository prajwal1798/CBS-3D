#include "cbs/parallel/PartitionMetadata.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#ifdef CBS3D_USE_MPI
#include <mpi.h>
#endif

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

#ifdef CBS3D_USE_MPI
        void checkMpi(const int error_code, const char* operation)
        {
            if (error_code != MPI_SUCCESS)
            {
                throw std::runtime_error(
                    std::string("PartitionMetadataIO MPI failure in ")
                    + operation);
            }
        }

        long double squareSumTo(const long double value)
        {
            if (value <= 0.0L)
            {
                return 0.0L;
            }

            return value * (value + 1.0L) *
                   (2.0L * value + 1.0L) / 6.0L;
        }

        bool nearlyEqualLongDouble(
            const long double lhs,
            const long double rhs)
        {
            const long double scale =
                std::max(
                    1.0L,
                    std::max(std::fabs(lhs), std::fabs(rhs)));

            return std::fabs(lhs - rhs) <= 1.0e-12L * scale;
        }

        // Older partition exports preserved the original Gmsh volume-element
        // tags. When Gmsh numbered boundary triangles before tetrahedra, the
        // tetrahedral tags form one contiguous globally shifted interval rather
        // than 1..global_nelem. The solver only requires a unique dense global
        // tetrahedron numbering, so detect that interval collectively and
        // remove its constant offset on every MPI rank.
        void normaliseGlobalElementIds(PartitionMetadata& metadata)
        {
            int mpi_initialised = 0;
            checkMpi(
                MPI_Initialized(&mpi_initialised),
                "MPI_Initialized");

            if (mpi_initialised == 0)
            {
                return;
            }

            long long local_count = metadata.local_nelem;
            long long global_count = 0;

            std::int64_t local_min =
                std::numeric_limits<std::int64_t>::max();
            std::int64_t local_max =
                std::numeric_limits<std::int64_t>::min();

            long double local_sum = 0.0L;
            long double local_square_sum = 0.0L;

            for (Int local_element = 1;
                 local_element <= metadata.local_nelem;
                 ++local_element)
            {
                const std::int64_t global_id =
                    metadata.local_to_global_element[
                        static_cast<Size>(local_element)];

                if (global_id < 1)
                {
                    throw std::runtime_error(
                        "PartitionMetadataIO - global element ID must be positive");
                }

                local_min = std::min(local_min, global_id);
                local_max = std::max(local_max, global_id);

                const long double value =
                    static_cast<long double>(global_id);

                local_sum += value;
                local_square_sum += value * value;
            }

            std::int64_t global_min = 0;
            std::int64_t global_max = 0;
            long double global_sum = 0.0L;
            long double global_square_sum = 0.0L;

            checkMpi(
                MPI_Allreduce(
                    &local_count,
                    &global_count,
                    1,
                    MPI_LONG_LONG,
                    MPI_SUM,
                    MPI_COMM_WORLD),
                "MPI_Allreduce global element count");

            checkMpi(
                MPI_Allreduce(
                    &local_min,
                    &global_min,
                    1,
                    MPI_LONG_LONG,
                    MPI_MIN,
                    MPI_COMM_WORLD),
                "MPI_Allreduce minimum global element ID");

            checkMpi(
                MPI_Allreduce(
                    &local_max,
                    &global_max,
                    1,
                    MPI_LONG_LONG,
                    MPI_MAX,
                    MPI_COMM_WORLD),
                "MPI_Allreduce maximum global element ID");

            checkMpi(
                MPI_Allreduce(
                    &local_sum,
                    &global_sum,
                    1,
                    MPI_LONG_DOUBLE,
                    MPI_SUM,
                    MPI_COMM_WORLD),
                "MPI_Allreduce global element ID sum");

            checkMpi(
                MPI_Allreduce(
                    &local_square_sum,
                    &global_square_sum,
                    1,
                    MPI_LONG_DOUBLE,
                    MPI_SUM,
                    MPI_COMM_WORLD),
                "MPI_Allreduce global element ID square sum");

            if (global_count != metadata.global_nelem)
            {
                throw std::runtime_error(
                    "PartitionMetadataIO - global element record count does not "
                    "match GLOBAL nelem");
            }

            const std::int64_t interval_count =
                global_max - global_min + 1;

            const long double expected_sum =
                static_cast<long double>(global_count) *
                static_cast<long double>(global_min + global_max) / 2.0L;

            const long double expected_square_sum =
                squareSumTo(static_cast<long double>(global_max)) -
                squareSumTo(static_cast<long double>(global_min - 1));

            const bool contiguous_global_interval =
                interval_count == global_count &&
                nearlyEqualLongDouble(global_sum, expected_sum) &&
                nearlyEqualLongDouble(
                    global_square_sum,
                    expected_square_sum);

            if (!contiguous_global_interval)
            {
                throw std::runtime_error(
                    "PartitionMetadataIO - global element IDs are neither a "
                    "dense nor a constant-offset contiguous numbering");
            }

            const std::int64_t global_offset = global_min - 1;

            for (Int local_element = 1;
                 local_element <= metadata.local_nelem;
                 ++local_element)
            {
                std::int64_t& global_id =
                    metadata.local_to_global_element[
                        static_cast<Size>(local_element)];

                global_id -= global_offset;

                if (global_id < 1 ||
                    global_id > metadata.global_nelem)
                {
                    throw std::runtime_error(
                        "PartitionMetadataIO - normalised global element ID is "
                        "outside 1..global_nelem");
                }
            }

            if (metadata.mpi_rank == 0 && global_offset != 0)
            {
                std::cout
                    << "Normalised legacy global element IDs by offset "
                    << global_offset
                    << " to dense range 1.."
                    << metadata.global_nelem
                    << "\n";
            }
        }
#endif
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

#ifdef CBS3D_USE_MPI
        normaliseGlobalElementIds(metadata);
#endif

        return metadata;
    }
}
