#pragma once

//=============================================================================
// CBS3D++_SI
//
// Persistent distributed symmetry/slip projector for BC 506.
//
// A shared node can belong to symmetry faces carried by different MPI ranks and
// can also lie at the intersection of two or three symmetry planes.  Projecting
// independently with one rank-local face normal is therefore partition
// dependent and incorrect at edges/corners.
//
// For every symmetry node we assemble the symmetric positive-semidefinite
// metric
//
//     G_i = sum_{f incident on i} A_f n_f n_f^T
//
// over all BC 506 faces.  Contributions are reverse-summed onto the node owner
// and broadcast to ghosts.  The range of G_i is the span of all independent
// symmetry normals at the node.  If Q contains an orthonormal basis of that
// range, the admissible tangential-velocity projector is
//
//     P_i = I - Q Q^T.
//
// Applying u <- P_i u therefore enforces all incident impermeability
// constraints simultaneously while preserving every admissible tangential
// component.  The projector is built once during distributed setup and reused
// after each velocity update.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"
#include "cbs/parallel/HaloExchange.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#ifdef CBS3D_USE_MPI
#include <mpi.h>
#endif

namespace cbs
{
    class DistributedSymmetryProjector
    {
    public:
        void initialise(CBSStateSI& s)
        {
            projector_.resize(9, s.cfg.npoin);
            projector_.fill(0.0);

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                set_identity(ip);
            }

            bool have_symmetry = false;

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                if (s.node_symmetry(ip) != 0)
                {
                    have_symmetry = true;
                    break;
                }
            }

#ifdef CBS3D_USE_MPI
            if (s.mpi_enabled && s.mpi_size > 1)
            {
                int local = have_symmetry ? 1 : 0;
                int global = 0;

                if (MPI_Allreduce(
                        &local,
                        &global,
                        1,
                        MPI_INT,
                        MPI_MAX,
                        MPI_COMM_WORLD) != MPI_SUCCESS)
                {
                    throw std::runtime_error(
                        "DistributedSymmetryProjector - MPI_Allreduce failed");
                }

                have_symmetry = global != 0;
            }
#endif

            if (!have_symmetry)
            {
                ready_ = true;
                return;
            }

#ifndef CBS3D_USE_MPI
            if (s.mpi_enabled)
            {
                throw std::runtime_error(
                    "DistributedSymmetryProjector requires an MPI-enabled build");
            }
#endif

            Array2D<Real> metric;
            metric.resize(9, s.cfg.npoin);
            metric.fill(0.0);

            for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
            {
                if (s.iside(s.cfg.bsid, ib) != s.cfg.bc_symmetry_no_flux)
                {
                    continue;
                }

                const Real area = s.face_norm(4, ib);

                if (!(area > 0.0) || !std::isfinite(area))
                {
                    throw std::runtime_error(
                        "DistributedSymmetryProjector - invalid BC 506 face area");
                }

                const std::array<Real, 3> area_normal =
                {
                    s.face_norm(1, ib),
                    s.face_norm(2, ib),
                    s.face_norm(3, ib)
                };

                for (const Real value : area_normal)
                {
                    if (!std::isfinite(value))
                    {
                        throw std::runtime_error(
                            "DistributedSymmetryProjector - non-finite BC 506 normal");
                    }
                }

                // face_norm stores A*n, hence
                //
                //     (A n)(A n)^T / A = A n n^T.
                for (Int in = 1; in <= s.cfg.nsidp; ++in)
                {
                    const Int ip = s.iside(in, ib);

                    if (ip < 1 || ip > s.cfg.npoin)
                    {
                        throw std::runtime_error(
                            "DistributedSymmetryProjector - BC 506 node out of range");
                    }

                    for (Int row = 0; row < 3; ++row)
                    {
                        for (Int col = 0; col < 3; ++col)
                        {
                            metric(3 * row + col + 1, ip) +=
                                area_normal[static_cast<std::size_t>(row)] *
                                area_normal[static_cast<std::size_t>(col)] /
                                area;
                        }
                    }
                }
            }

#ifdef CBS3D_USE_MPI
            if (s.mpi_enabled && s.mpi_size > 1)
            {
                HaloExchange::sumGhostContributionsToOwners(
                    metric,
                    s.partition_metadata,
                    MPI_COMM_WORLD);

                HaloExchange::broadcastOwnedToGhosts(
                    metric,
                    s.partition_metadata,
                    MPI_COMM_WORLD);
            }
#endif

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                if (s.node_symmetry(ip) == 0)
                {
                    continue;
                }

                build_node_projector(metric, ip);
            }

#ifdef CBS3D_USE_MPI
            if (s.mpi_enabled && s.mpi_size > 1)
            {
                // Owner-authoritative final copy.  Every rank should derive the
                // same projector from the reconciled metric; this broadcast also
                // removes any last-bit roundoff difference at shared nodes.
                HaloExchange::broadcastOwnedToGhosts(
                    projector_,
                    s.partition_metadata,
                    MPI_COMM_WORLD);
            }
#endif

            ready_ = true;
        }

        void applyOwned(CBSStateSI& s) const
        {
            if (!ready_)
            {
                throw std::runtime_error(
                    "DistributedSymmetryProjector used before initialise");
            }

            if (!s.mpi_enabled)
            {
                throw std::runtime_error(
                    "DistributedSymmetryProjector::applyOwned requires MPI state");
            }

            for (const Int ip : s.owned_nodes)
            {
                if (s.node_symmetry(ip) == 0)
                {
                    continue;
                }

                const std::array<Real, 3> old_velocity =
                {
                    s.unkno(1, ip),
                    s.unkno(2, ip),
                    s.unkno(3, ip)
                };

                std::array<Real, 3> projected = {0.0, 0.0, 0.0};

                for (Int row = 0; row < 3; ++row)
                {
                    for (Int col = 0; col < 3; ++col)
                    {
                        projected[static_cast<std::size_t>(row)] +=
                            projector_(3 * row + col + 1, ip) *
                            old_velocity[static_cast<std::size_t>(col)];
                    }
                }

                if (!std::isfinite(projected[0]) ||
                    !std::isfinite(projected[1]) ||
                    !std::isfinite(projected[2]))
                {
                    throw std::runtime_error(
                        "DistributedSymmetryProjector produced non-finite velocity");
                }

                s.unkno(1, ip) = projected[0];
                s.unkno(2, ip) = projected[1];
                s.unkno(3, ip) = projected[2];
            }
        }

    private:
        Array2D<Real> projector_;
        bool ready_ = false;

        void set_identity(const Int ip)
        {
            for (Int row = 0; row < 3; ++row)
            {
                for (Int col = 0; col < 3; ++col)
                {
                    projector_(3 * row + col + 1, ip) =
                        row == col ? 1.0 : 0.0;
                }
            }
        }

        void build_node_projector(
            const Array2D<Real>& metric,
            const Int ip)
        {
            Real scale = 0.0;

            for (Int entry = 1; entry <= 9; ++entry)
            {
                const Real value = metric(entry, ip);

                if (!std::isfinite(value))
                {
                    throw std::runtime_error(
                        "DistributedSymmetryProjector found non-finite symmetry metric");
                }

                scale = std::max(scale, std::abs(value));
            }

            if (!(scale > 0.0) || !std::isfinite(scale))
            {
                throw std::runtime_error(
                    "DistributedSymmetryProjector - symmetry node has zero global normal metric");
            }

            std::array<std::array<Real, 3>, 3> basis{};
            Int rank = 0;
            const Real tolerance = 1.0e-11 * scale;

            // Modified Gram-Schmidt on the columns of the symmetric metric.
            // The accepted columns span range(G), i.e. all independent symmetry
            // normal directions incident on this node.
            for (Int col = 0; col < 3; ++col)
            {
                std::array<Real, 3> v =
                {
                    metric(col + 1, ip),
                    metric(3 + col + 1, ip),
                    metric(6 + col + 1, ip)
                };

                for (Int q = 0; q < rank; ++q)
                {
                    Real projection = 0.0;

                    for (Int k = 0; k < 3; ++k)
                    {
                        projection +=
                            v[static_cast<std::size_t>(k)] *
                            basis[static_cast<std::size_t>(q)]
                                 [static_cast<std::size_t>(k)];
                    }

                    for (Int k = 0; k < 3; ++k)
                    {
                        v[static_cast<std::size_t>(k)] -=
                            projection *
                            basis[static_cast<std::size_t>(q)]
                                 [static_cast<std::size_t>(k)];
                    }
                }

                Real norm_squared = 0.0;

                for (const Real value : v)
                {
                    norm_squared += value * value;
                }

                const Real norm = std::sqrt(norm_squared);

                if (norm <= tolerance)
                {
                    continue;
                }

                if (rank >= 3)
                {
                    throw std::runtime_error(
                        "DistributedSymmetryProjector - invalid symmetry metric rank");
                }

                for (Int k = 0; k < 3; ++k)
                {
                    basis[static_cast<std::size_t>(rank)]
                         [static_cast<std::size_t>(k)] =
                        v[static_cast<std::size_t>(k)] / norm;
                }

                ++rank;
            }

            if (rank < 1 || rank > 3)
            {
                throw std::runtime_error(
                    "DistributedSymmetryProjector could not recover symmetry-normal span");
            }

            set_identity(ip);

            for (Int q = 0; q < rank; ++q)
            {
                for (Int row = 0; row < 3; ++row)
                {
                    for (Int col = 0; col < 3; ++col)
                    {
                        projector_(3 * row + col + 1, ip) -=
                            basis[static_cast<std::size_t>(q)]
                                 [static_cast<std::size_t>(row)] *
                            basis[static_cast<std::size_t>(q)]
                                 [static_cast<std::size_t>(col)];
                    }
                }
            }

            // Defensive algebraic audit: P must be finite, symmetric and
            // idempotent to roundoff.  A failed projector is much safer to stop
            // during setup than to contaminate the velocity field silently.
            Real symmetry_error = 0.0;
            Real idempotence_error = 0.0;

            for (Int row = 0; row < 3; ++row)
            {
                for (Int col = 0; col < 3; ++col)
                {
                    const Real p_rc =
                        projector_(3 * row + col + 1, ip);

                    const Real p_cr =
                        projector_(3 * col + row + 1, ip);

                    if (!std::isfinite(p_rc))
                    {
                        throw std::runtime_error(
                            "DistributedSymmetryProjector produced non-finite projector");
                    }

                    symmetry_error =
                        std::max(symmetry_error, std::abs(p_rc - p_cr));

                    Real p2_rc = 0.0;

                    for (Int k = 0; k < 3; ++k)
                    {
                        p2_rc +=
                            projector_(3 * row + k + 1, ip) *
                            projector_(3 * k + col + 1, ip);
                    }

                    idempotence_error =
                        std::max(
                            idempotence_error,
                            std::abs(p2_rc - p_rc));
                }
            }

            if (symmetry_error > 1.0e-9 ||
                idempotence_error > 1.0e-9)
            {
                throw std::runtime_error(
                    "DistributedSymmetryProjector failed projector algebra audit");
            }
        }
    };
}
