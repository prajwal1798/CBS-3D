#include "cbs/parallel/HaloExchange.hpp"

#ifdef CBS3D_USE_MPI

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

namespace cbs
{
    namespace
    {
        void check_mpi(const int error_code, const char* context)
        {
            if (error_code != MPI_SUCCESS)
            {
                throw std::runtime_error(
                    std::string("HaloExchange extrema - MPI failure in ")
                    + context);
            }
        }

        template <typename Combine>
        void reduce_ghost_values_to_owners(
            Array1D<Real>& values,
            const PartitionMetadata& metadata,
            const int tag,
            Combine combine,
            MPI_Comm communicator)
        {
            std::vector<std::vector<Real>> send_buffers(
                metadata.neighbours.size());
            std::vector<std::vector<Real>> recv_buffers(
                metadata.neighbours.size());
            std::vector<MPI_Request> requests;
            requests.reserve(2U * metadata.neighbours.size());

            // Reverse halo direction: ghost copies send to the node owner.
            for (Size i = 0; i < metadata.neighbours.size(); ++i)
            {
                const PartitionNeighbour& neighbour = metadata.neighbours[i];
                const std::vector<Int>& receive_nodes =
                    neighbour.send_local_nodes;

                recv_buffers[i].resize(receive_nodes.size());

                MPI_Request request{};
                check_mpi(
                    MPI_Irecv(
                        recv_buffers[i].data(),
                        static_cast<int>(recv_buffers[i].size()),
                        MPI_DOUBLE,
                        neighbour.rank,
                        tag,
                        communicator,
                        &request),
                    "MPI_Irecv");
                requests.push_back(request);
            }

            for (Size i = 0; i < metadata.neighbours.size(); ++i)
            {
                const PartitionNeighbour& neighbour = metadata.neighbours[i];
                const std::vector<Int>& send_nodes =
                    neighbour.recv_local_nodes;

                send_buffers[i].resize(send_nodes.size());
                for (Size j = 0; j < send_nodes.size(); ++j)
                {
                    send_buffers[i][j] = values(send_nodes[j]);
                }

                MPI_Request request{};
                check_mpi(
                    MPI_Isend(
                        send_buffers[i].data(),
                        static_cast<int>(send_buffers[i].size()),
                        MPI_DOUBLE,
                        neighbour.rank,
                        tag,
                        communicator,
                        &request),
                    "MPI_Isend");
                requests.push_back(request);
            }

            if (!requests.empty())
            {
                check_mpi(
                    MPI_Waitall(
                        static_cast<int>(requests.size()),
                        requests.data(),
                        MPI_STATUSES_IGNORE),
                    "MPI_Waitall");
            }

            for (Size i = 0; i < metadata.neighbours.size(); ++i)
            {
                const PartitionNeighbour& neighbour = metadata.neighbours[i];
                const std::vector<Int>& receive_nodes =
                    neighbour.send_local_nodes;

                for (Size j = 0; j < receive_nodes.size(); ++j)
                {
                    const Int node = receive_nodes[j];
                    values(node) = combine(values(node), recv_buffers[i][j]);
                }
            }
        }
    }

    void HaloExchange::minGhostContributionsToOwners(
        Array1D<Real>& values,
        const PartitionMetadata& metadata,
        MPI_Comm communicator)
    {
        reduce_ghost_values_to_owners(
            values,
            metadata,
            816,
            [](const Real a, const Real b)
            {
                return std::min(a, b);
            },
            communicator);
    }

    void HaloExchange::maxGhostContributionsToOwners(
        Array1D<Real>& values,
        const PartitionMetadata& metadata,
        MPI_Comm communicator)
    {
        reduce_ghost_values_to_owners(
            values,
            metadata,
            817,
            [](const Real a, const Real b)
            {
                return std::max(a, b);
            },
            communicator);
    }
}

#endif
