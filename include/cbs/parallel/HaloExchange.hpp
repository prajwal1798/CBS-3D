#pragma once

#include "cbs/core/Array.hpp"
#include "cbs/parallel/PartitionMetadata.hpp"

#ifdef CBS3D_USE_MPI
#include <mpi.h>
#endif

namespace cbs
{
    class HaloExchange
    {
    public:
#ifdef CBS3D_USE_MPI
        static void broadcastOwnedToGhosts(
            Array1D<Real>& values,
            const PartitionMetadata& metadata,
            MPI_Comm communicator = MPI_COMM_WORLD);

        static void broadcastOwnedToGhosts(
            Array2D<Real>& values,
            const PartitionMetadata& metadata,
            MPI_Comm communicator = MPI_COMM_WORLD);

        static void sumGhostContributionsToOwners(
            Array1D<Real>& values,
            const PartitionMetadata& metadata,
            MPI_Comm communicator = MPI_COMM_WORLD);

        static void sumGhostContributionsToOwners(
            Array2D<Real>& values,
            const PartitionMetadata& metadata,
            MPI_Comm communicator = MPI_COMM_WORLD);
#endif
    };
}
