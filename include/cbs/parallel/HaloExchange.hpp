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

        // Copies integer owner values to all local ghost copies.
        static void broadcastOwnedToGhosts(
            Array1D<Int>& values,
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

        // Combines ghost material-mask bits on the owning node.
        // Reduces a nodal quantity to its minimum over every rank that holds a
        // copy of the node, leaving the result on the owner.
        //
        // Needed for local time stepping.  A node on a partition interface is
        // surrounded by elements on both sides of the cut, and each rank can
        // only see its own, so each rank's candidate timestep for that node is
        // the minimum over a subset and is therefore too large.  Without this
        // reduction two ranks would advance the same node with different dt,
        // which makes the assembled pressure operator rank-dependent.
        //
        // Follow with broadcastOwnedToGhosts to publish the reduced value back
        // to the ghost copies.
        static void minGhostContributionsToOwners(
            Array1D<Real>& values,
            const PartitionMetadata& metadata,
            MPI_Comm communicator);

        static void orGhostMasksToOwners(
            Array1D<Int>& values,
            const PartitionMetadata& metadata,
            MPI_Comm communicator = MPI_COMM_WORLD);

        static void sumGhostContributionsToOwners(
            Array2D<Real>& values,
            const PartitionMetadata& metadata,
            MPI_Comm communicator = MPI_COMM_WORLD);
#endif
    };
}
