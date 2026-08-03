#include "cbs/parallel/HaloExchange.hpp"

#ifdef CBS3D_USE_MPI

#include <stdexcept>
#include <vector>

namespace cbs
{
    namespace
    {
        void checkMpi(const int error_code, const char* context)
        {
            if (error_code != MPI_SUCCESS)
            {
                throw std::runtime_error(
                    std::string("HaloExchange - MPI failure in ") + context);
            }
        }

        template <typename PackSend, typename ApplyRecv>
        void exchange(
            const PartitionMetadata& metadata,
            const Int values_per_node,
            const int tag,
            const bool reverse,
            PackSend pack_send,
            ApplyRecv apply_recv,
            MPI_Comm communicator)
        {
            std::vector<std::vector<Real>> send_buffers(metadata.neighbours.size());
            std::vector<std::vector<Real>> recv_buffers(metadata.neighbours.size());
            std::vector<MPI_Request> requests;
            requests.reserve(2U * metadata.neighbours.size());

            for (Size i = 0; i < metadata.neighbours.size(); ++i)
            {
                const PartitionNeighbour& neighbour = metadata.neighbours[i];
                const std::vector<Int>& recv_nodes = reverse
                    ? neighbour.send_local_nodes
                    : neighbour.recv_local_nodes;

                recv_buffers[i].resize(
                    static_cast<Size>(values_per_node)
                    * recv_nodes.size());

                MPI_Request request{};
                checkMpi(
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
                const std::vector<Int>& send_nodes = reverse
                    ? neighbour.recv_local_nodes
                    : neighbour.send_local_nodes;

                pack_send(send_nodes, send_buffers[i]);

                MPI_Request request{};
                checkMpi(
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
                checkMpi(
                    MPI_Waitall(
                        static_cast<int>(requests.size()),
                        requests.data(),
                        MPI_STATUSES_IGNORE),
                    "MPI_Waitall");
            }

            for (Size i = 0; i < metadata.neighbours.size(); ++i)
            {
                const PartitionNeighbour& neighbour = metadata.neighbours[i];
                const std::vector<Int>& recv_nodes = reverse
                    ? neighbour.send_local_nodes
                    : neighbour.recv_local_nodes;
                apply_recv(recv_nodes, recv_buffers[i]);
            }
        }

        // Integer counterpart of exchange(). It is kept separate from the
        // established Real/MPI_DOUBLE path so that DD-1 behaviour is unchanged.
        template <typename PackSend, typename ApplyRecv>
        void exchangeInt(
            const PartitionMetadata& metadata,
            const int tag,
            const bool reverse,
            PackSend pack_send,
            ApplyRecv apply_recv,
            MPI_Comm communicator)
        {
            std::vector<std::vector<Int>> send_buffers(metadata.neighbours.size());
            std::vector<std::vector<Int>> recv_buffers(metadata.neighbours.size());
            std::vector<MPI_Request> requests;
            requests.reserve(2U * metadata.neighbours.size());

            for (Size i = 0; i < metadata.neighbours.size(); ++i)
            {
                const PartitionNeighbour& neighbour = metadata.neighbours[i];

                const std::vector<Int>& recv_nodes = reverse
                    ? neighbour.send_local_nodes
                    : neighbour.recv_local_nodes;

                recv_buffers[i].resize(recv_nodes.size());

                MPI_Request request{};
                checkMpi(
                    MPI_Irecv(
                        recv_buffers[i].data(),
                        static_cast<int>(recv_buffers[i].size()),
                        MPI_INT,
                        neighbour.rank,
                        tag,
                        communicator,
                        &request),
                    "MPI_Irecv integer halo");

                requests.push_back(request);
            }

            for (Size i = 0; i < metadata.neighbours.size(); ++i)
            {
                const PartitionNeighbour& neighbour = metadata.neighbours[i];

                const std::vector<Int>& send_nodes = reverse
                    ? neighbour.recv_local_nodes
                    : neighbour.send_local_nodes;

                pack_send(send_nodes, send_buffers[i]);

                MPI_Request request{};
                checkMpi(
                    MPI_Isend(
                        send_buffers[i].data(),
                        static_cast<int>(send_buffers[i].size()),
                        MPI_INT,
                        neighbour.rank,
                        tag,
                        communicator,
                        &request),
                    "MPI_Isend integer halo");

                requests.push_back(request);
            }

            if (!requests.empty())
            {
                checkMpi(
                    MPI_Waitall(
                        static_cast<int>(requests.size()),
                        requests.data(),
                        MPI_STATUSES_IGNORE),
                    "MPI_Waitall integer halo");
            }

            for (Size i = 0; i < metadata.neighbours.size(); ++i)
            {
                const PartitionNeighbour& neighbour = metadata.neighbours[i];

                const std::vector<Int>& recv_nodes = reverse
                    ? neighbour.send_local_nodes
                    : neighbour.recv_local_nodes;

                apply_recv(recv_nodes, recv_buffers[i]);
            }
        }
    }

    void HaloExchange::broadcastOwnedToGhosts(
        Array1D<Real>& values,
        const PartitionMetadata& metadata,
        MPI_Comm communicator)
    {
        exchange(
            metadata,
            1,
            810,
            false,
            [&values](const std::vector<Int>& nodes, std::vector<Real>& buffer)
            {
                buffer.resize(nodes.size());
                for (Size i = 0; i < nodes.size(); ++i)
                {
                    buffer[i] = values(nodes[i]);
                }
            },
            [&values](const std::vector<Int>& nodes, const std::vector<Real>& buffer)
            {
                for (Size i = 0; i < nodes.size(); ++i)
                {
                    values(nodes[i]) = buffer[i];
                }
            },
            communicator);
    }

    void HaloExchange::broadcastOwnedToGhosts(
        Array1D<Int>& values,
        const PartitionMetadata& metadata,
        MPI_Comm communicator)
    {
        exchangeInt(
            metadata,
            814,
            false,
            [&values](
                const std::vector<Int>& nodes,
                std::vector<Int>& buffer)
            {
                buffer.resize(nodes.size());

                for (Size i = 0; i < nodes.size(); ++i)
                {
                    buffer[i] = values(nodes[i]);
                }
            },
            [&values](
                const std::vector<Int>& nodes,
                const std::vector<Int>& buffer)
            {
                for (Size i = 0; i < nodes.size(); ++i)
                {
                    values(nodes[i]) = buffer[i];
                }
            },
            communicator);
    }


    void HaloExchange::broadcastOwnedToGhosts(
        Array2D<Real>& values,
        const PartitionMetadata& metadata,
        MPI_Comm communicator)
    {
        const Int components = values.dim1();
        exchange(
            metadata,
            components,
            811,
            false,
            [&values, components](const std::vector<Int>& nodes, std::vector<Real>& buffer)
            {
                buffer.resize(static_cast<Size>(components) * nodes.size());
                Size k = 0;
                for (const Int node : nodes)
                {
                    for (Int component = 1; component <= components; ++component)
                    {
                        buffer[k++] = values(component, node);
                    }
                }
            },
            [&values, components](const std::vector<Int>& nodes, const std::vector<Real>& buffer)
            {
                Size k = 0;
                for (const Int node : nodes)
                {
                    for (Int component = 1; component <= components; ++component)
                    {
                        values(component, node) = buffer[k++];
                    }
                }
            },
            communicator);
    }

    void HaloExchange::sumGhostContributionsToOwners(
        Array1D<Real>& values,
        const PartitionMetadata& metadata,
        MPI_Comm communicator)
    {
        exchange(
            metadata,
            1,
            812,
            true,
            [&values](const std::vector<Int>& nodes, std::vector<Real>& buffer)
            {
                buffer.resize(nodes.size());
                for (Size i = 0; i < nodes.size(); ++i)
                {
                    buffer[i] = values(nodes[i]);
                }
            },
            [&values](const std::vector<Int>& nodes, const std::vector<Real>& buffer)
            {
                for (Size i = 0; i < nodes.size(); ++i)
                {
                    values(nodes[i]) += buffer[i];
                }
            },
            communicator);
    }

    void HaloExchange::orGhostMasksToOwners(
        Array1D<Int>& values,
        const PartitionMetadata& metadata,
        MPI_Comm communicator)
    {
        exchangeInt(
            metadata,
            815,
            true,
            [&values](
                const std::vector<Int>& nodes,
                std::vector<Int>& buffer)
            {
                buffer.resize(nodes.size());

                for (Size i = 0; i < nodes.size(); ++i)
                {
                    buffer[i] = values(nodes[i]);
                }
            },
            [&values](
                const std::vector<Int>& nodes,
                const std::vector<Int>& buffer)
            {
                for (Size i = 0; i < nodes.size(); ++i)
                {
                    values(nodes[i]) |= buffer[i];
                }
            },
            communicator);
    }


    void HaloExchange::sumGhostContributionsToOwners(
        Array2D<Real>& values,
        const PartitionMetadata& metadata,
        MPI_Comm communicator)
    {
        const Int components = values.dim1();
        exchange(
            metadata,
            components,
            813,
            true,
            [&values, components](const std::vector<Int>& nodes, std::vector<Real>& buffer)
            {
                buffer.resize(static_cast<Size>(components) * nodes.size());
                Size k = 0;
                for (const Int node : nodes)
                {
                    for (Int component = 1; component <= components; ++component)
                    {
                        buffer[k++] = values(component, node);
                    }
                }
            },
            [&values, components](const std::vector<Int>& nodes, const std::vector<Real>& buffer)
            {
                Size k = 0;
                for (const Int node : nodes)
                {
                    for (Int component = 1; component <= components; ++component)
                    {
                        values(component, node) += buffer[k++];
                    }
                }
            },
            communicator);
    }
}

#endif
