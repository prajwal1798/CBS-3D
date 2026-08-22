#pragma once

//=============================================================================
// CBS3D++_SI
//
// Authoritative turbulence-wall topology for Spalart-Allmaras preprocessing.
//
// A thermal boundary face is not automatically a turbulence wall. In CHT,
// BC 532 may belong to the exterior of a solid and therefore must not enter the
// fluid SA wall-distance geometry. Conversely, the conformal fluid-solid
// interface is an internal tetrahedral face and is not guaranteed to be present
// in the .plt boundary-face list.
//
// This module therefore defines the turbulence wall from two sources:
//
//   1. explicit wall-type boundary faces whose parent tetrahedron is fluid;
//   2. every tetrahedral face shared by one fluid and one solid element.
//
// In MPI calculations, an interface face can lie on a partition cut. Each rank
// first identifies locally unmatched faces whose three nodes are known material-
// interface nodes, then all ranks match those faces by sorted global node IDs.
// Only a fluid/solid pair is accepted as a reconstructed CHT wall. Fluid/fluid
// and solid/solid partition faces are communication surfaces, not walls.
//
// The same material/topology rules are used to reconcile SA wall-node flags.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"
#include "cbs/turbulence/WallDistance.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef CBS3D_USE_MPI
#include <mpi.h>
#endif

namespace cbs
{
    namespace turbulence
    {
        class TurbulenceWallTopology
        {
        private:
            struct KeyedTriangle
            {
                std::array<long long, 3> key{};
                WallTriangle triangle{};
            };

            struct FaceRecord
            {
                std::array<Int, 3> local_key{};
                KeyedTriangle keyed_triangle{};
                bool fluid = false;
            };

            static bool isExplicitWallBc(const CBSStateSI& s, const Int bc)
            {
                return bc == s.cfg.bc_noslip_adiabatic_wall
                    || bc == s.cfg.bc_noslip_heatflux_wall
                    || bc == s.cfg.bc_cht_interface;
            }

            static bool isMaterialInterfaceNode(
                const CBSStateSI& s,
                const Int ip)
            {
                if (ip < 1 || ip > s.cfg.npoin)
                {
                    throw std::runtime_error(
                        "TurbulenceWallTopology - node index outside mesh");
                }

                const Int interface_mask =
                    CBSStateSI::node_touches_fluid |
                    CBSStateSI::node_touches_solid;

                return s.node_material_mask(ip) == interface_mask;
            }

            static long long globalNodeId(
                const CBSStateSI& s,
                const Int ip)
            {
#ifdef CBS3D_USE_MPI
                if (s.mpi_enabled && s.mpi_size > 1)
                {
                    if (s.local_to_global_node.size()
                        < static_cast<Size>(s.cfg.npoin) + 1U)
                    {
                        throw std::runtime_error(
                            "TurbulenceWallTopology - MPI global-node map is missing");
                    }

                    const Int global_id =
                        s.local_to_global_node[static_cast<Size>(ip)];

                    if (global_id < 1)
                    {
                        throw std::runtime_error(
                            "TurbulenceWallTopology - invalid MPI global-node ID");
                    }

                    return static_cast<long long>(global_id);
                }
#else
                (void)s;
#endif

                return static_cast<long long>(ip);
            }

            static std::array<Real, 3> point(
                const CBSStateSI& s,
                const Int ip)
            {
                const std::array<Real, 3> p =
                {
                    s.coord(1, ip),
                    s.coord(2, ip),
                    s.coord(3, ip)
                };

                if (!std::isfinite(p[0]) ||
                    !std::isfinite(p[1]) ||
                    !std::isfinite(p[2]))
                {
                    throw std::runtime_error(
                        "TurbulenceWallTopology - non-finite wall coordinate");
                }

                return p;
            }

            static KeyedTriangle keyedTriangle(
                const CBSStateSI& s,
                const std::array<Int, 3>& nodes)
            {
                struct NodeKey
                {
                    long long global_id = 0;
                    Int local_id = 0;
                };

                std::array<NodeKey, 3> order{};

                for (Size i = 0; i < 3; ++i)
                {
                    const Int ip = nodes[i];

                    if (ip < 1 || ip > s.cfg.npoin)
                    {
                        throw std::runtime_error(
                            "TurbulenceWallTopology - face node outside mesh");
                    }

                    order[i].global_id = globalNodeId(s, ip);
                    order[i].local_id = ip;
                }

                std::sort(
                    order.begin(),
                    order.end(),
                    [](const NodeKey& lhs, const NodeKey& rhs)
                    {
                        return lhs.global_id < rhs.global_id;
                    });

                if (order[0].global_id == order[1].global_id ||
                    order[1].global_id == order[2].global_id)
                {
                    throw std::runtime_error(
                        "TurbulenceWallTopology - degenerate face node IDs");
                }

                KeyedTriangle result;

                for (Size i = 0; i < 3; ++i)
                {
                    result.key[i] = order[i].global_id;
                }

                result.triangle.a = point(s, order[0].local_id);
                result.triangle.b = point(s, order[1].local_id);
                result.triangle.c = point(s, order[2].local_id);

                return result;
            }

            static bool sameKey(
                const std::array<long long, 3>& a,
                const std::array<long long, 3>& b)
            {
                return a[0] == b[0] &&
                       a[1] == b[1] &&
                       a[2] == b[2];
            }

            static bool sameLocalKey(
                const std::array<Int, 3>& a,
                const std::array<Int, 3>& b)
            {
                return a[0] == b[0] &&
                       a[1] == b[1] &&
                       a[2] == b[2];
            }

            static void validateDuplicateGeometry(
                const KeyedTriangle& a,
                const KeyedTriangle& b)
            {
                const std::array<std::array<Real, 3>, 3> pa =
                {
                    a.triangle.a,
                    a.triangle.b,
                    a.triangle.c
                };

                const std::array<std::array<Real, 3>, 3> pb =
                {
                    b.triangle.a,
                    b.triangle.b,
                    b.triangle.c
                };

                Real scale = 1.0;
                Real maximum_difference = 0.0;

                for (Size i = 0; i < 3; ++i)
                {
                    for (Size k = 0; k < 3; ++k)
                    {
                        scale = std::max(scale, std::fabs(pa[i][k]));
                        scale = std::max(scale, std::fabs(pb[i][k]));
                        maximum_difference = std::max(
                            maximum_difference,
                            std::fabs(pa[i][k] - pb[i][k]));
                    }
                }

                if (maximum_difference > 1.0e-11 * scale)
                {
                    throw std::runtime_error(
                        "TurbulenceWallTopology - duplicate global face has inconsistent coordinates");
                }
            }

            static void deduplicate(std::vector<KeyedTriangle>& triangles)
            {
                std::sort(
                    triangles.begin(),
                    triangles.end(),
                    [](const KeyedTriangle& lhs, const KeyedTriangle& rhs)
                    {
                        return lhs.key < rhs.key;
                    });

                std::vector<KeyedTriangle> unique;
                unique.reserve(triangles.size());

                for (const KeyedTriangle& triangle : triangles)
                {
                    if (!unique.empty() && sameKey(unique.back().key, triangle.key))
                    {
                        validateDuplicateGeometry(unique.back(), triangle);
                        continue;
                    }

                    unique.push_back(triangle);
                }

                triangles.swap(unique);
            }

            static FaceRecord elementFace(
                const CBSStateSI& s,
                const Int ie,
                const Int face)
            {
                if (ie < 1 || ie > s.cfg.nelem ||
                    face < 1 || face > s.cfg.nsid)
                {
                    throw std::runtime_error(
                        "TurbulenceWallTopology - invalid element face");
                }

                std::array<Int, 3> nodes =
                {
                    s.intma(s.ippn1(face, 1), ie),
                    s.intma(s.ippn1(face, 2), ie),
                    s.intma(s.ippn1(face, 3), ie)
                };

                FaceRecord record;
                record.local_key = nodes;
                std::sort(record.local_key.begin(), record.local_key.end());
                record.keyed_triangle = keyedTriangle(s, nodes);
                record.fluid = s.mat_elem(ie) == 0;

                return record;
            }

            static bool faceCanBeMaterialInterface(
                const CBSStateSI& s,
                const FaceRecord& face)
            {
                for (const Int ip : face.local_key)
                {
                    if (!isMaterialInterfaceNode(s, ip))
                    {
                        return false;
                    }
                }

                return true;
            }

            static std::vector<KeyedTriangle> collectLocal(
                const CBSStateSI& s,
                std::vector<FaceRecord>& partition_candidates)
            {
                std::vector<KeyedTriangle> triangles;
                partition_candidates.clear();

                // Explicit boundary faces are turbulence walls only when their
                // parent tetrahedron belongs to the fluid. This excludes, for
                // example, a plasma-facing BC 532 on the exterior of EUROFER.
                for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
                {
                    const Int bc = s.iside(s.cfg.bsid, ib);

                    if (!isExplicitWallBc(s, bc))
                    {
                        continue;
                    }

                    const Int parent = s.iside(s.cfg.nsidpe, ib);

                    if (parent < 1 || parent > s.cfg.nelem)
                    {
                        throw std::runtime_error(
                            "TurbulenceWallTopology - wall boundary face has invalid parent tetrahedron");
                    }

                    if (s.mat_elem(parent) != 0)
                    {
                        continue;
                    }

                    const std::array<Int, 3> nodes =
                    {
                        s.iside(1, ib),
                        s.iside(2, ib),
                        s.iside(3, ib)
                    };

                    triangles.push_back(keyedTriangle(s, nodes));
                }

                // The material interface is internal and may be absent from
                // iside. Only faces whose three nodes are globally known to
                // touch both fluid and solid can possibly be CHT interfaces;
                // this filter keeps the face inventory small on production
                // meshes while preserving exactness.
                std::vector<FaceRecord> candidate_faces;

                for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
                {
                    for (Int face = 1; face <= s.cfg.nsid; ++face)
                    {
                        FaceRecord record = elementFace(s, ie, face);

                        if (faceCanBeMaterialInterface(s, record))
                        {
                            candidate_faces.push_back(record);
                        }
                    }
                }

                std::sort(
                    candidate_faces.begin(),
                    candidate_faces.end(),
                    [](const FaceRecord& lhs, const FaceRecord& rhs)
                    {
                        return lhs.local_key < rhs.local_key;
                    });

                Size begin = 0;

                while (begin < candidate_faces.size())
                {
                    Size end = begin + 1U;

                    while (end < candidate_faces.size() &&
                           sameLocalKey(
                               candidate_faces[begin].local_key,
                               candidate_faces[end].local_key))
                    {
                        ++end;
                    }

                    const Size multiplicity = end - begin;

                    if (multiplicity > 2U)
                    {
                        throw std::runtime_error(
                            "TurbulenceWallTopology - non-manifold tetrahedral face detected");
                    }

                    if (multiplicity == 2U)
                    {
                        const FaceRecord& first = candidate_faces[begin];
                        const FaceRecord& second = candidate_faces[begin + 1U];

                        if (first.fluid != second.fluid)
                        {
                            triangles.push_back(
                                first.fluid
                                    ? first.keyed_triangle
                                    : second.keyed_triangle);
                        }
                    }
                    else
                    {
#ifdef CBS3D_USE_MPI
                        if (s.mpi_enabled && s.mpi_size > 1)
                        {
                            partition_candidates.push_back(candidate_faces[begin]);
                        }
#endif
                    }

                    begin = end;
                }

                deduplicate(triangles);
                return triangles;
            }

#ifdef CBS3D_USE_MPI
            static void checkMpi(const int code, const char* operation)
            {
                if (code != MPI_SUCCESS)
                {
                    throw std::runtime_error(
                        std::string("TurbulenceWallTopology MPI failure in ")
                        + operation);
                }
            }

            static std::vector<KeyedTriangle> collectCrossRankInterfaces(
                const CBSStateSI& s,
                const std::vector<FaceRecord>& local_candidates)
            {
                if (!s.mpi_enabled || s.mpi_size <= 1)
                {
                    return {};
                }

                if (local_candidates.size()
                    > static_cast<Size>(std::numeric_limits<int>::max()))
                {
                    throw std::runtime_error(
                        "TurbulenceWallTopology - too many local partition-face candidates for MPI");
                }

                const int local_count =
                    static_cast<int>(local_candidates.size());

                std::vector<int> counts(static_cast<Size>(s.mpi_size), 0);

                checkMpi(
                    MPI_Allgather(
                        &local_count,
                        1,
                        MPI_INT,
                        counts.data(),
                        1,
                        MPI_INT,
                        MPI_COMM_WORLD),
                    "MPI_Allgather partition-face counts");

                std::vector<int> offsets(static_cast<Size>(s.mpi_size), 0);
                long long total_faces_ll = 0;

                for (Int rank = 0; rank < s.mpi_size; ++rank)
                {
                    offsets[static_cast<Size>(rank)] =
                        static_cast<int>(total_faces_ll);
                    total_faces_ll += counts[static_cast<Size>(rank)];
                }

                if (total_faces_ll > std::numeric_limits<int>::max())
                {
                    throw std::runtime_error(
                        "TurbulenceWallTopology - global partition-face candidate count exceeds MPI limit");
                }

                const int total_faces = static_cast<int>(total_faces_ll);

                std::vector<long long> local_keys(
                    static_cast<Size>(local_count) * 3U);
                std::vector<int> local_material(
                    static_cast<Size>(local_count));
                std::vector<double> local_coordinates(
                    static_cast<Size>(local_count) * 9U);

                for (int i = 0; i < local_count; ++i)
                {
                    const FaceRecord& face =
                        local_candidates[static_cast<Size>(i)];

                    for (Size k = 0; k < 3; ++k)
                    {
                        local_keys[static_cast<Size>(i) * 3U + k] =
                            face.keyed_triangle.key[k];
                    }

                    local_material[static_cast<Size>(i)] =
                        face.fluid ? 1 : 0;

                    const std::array<std::array<Real, 3>, 3> vertices =
                    {
                        face.keyed_triangle.triangle.a,
                        face.keyed_triangle.triangle.b,
                        face.keyed_triangle.triangle.c
                    };

                    for (Size corner = 0; corner < 3; ++corner)
                    {
                        for (Size k = 0; k < 3; ++k)
                        {
                            local_coordinates[
                                static_cast<Size>(i) * 9U +
                                corner * 3U + k] = vertices[corner][k];
                        }
                    }
                }

                std::vector<int> key_counts(static_cast<Size>(s.mpi_size));
                std::vector<int> key_offsets(static_cast<Size>(s.mpi_size));
                std::vector<int> coordinate_counts(static_cast<Size>(s.mpi_size));
                std::vector<int> coordinate_offsets(static_cast<Size>(s.mpi_size));

                for (Int rank = 0; rank < s.mpi_size; ++rank)
                {
                    const long long count = counts[static_cast<Size>(rank)];
                    const long long offset = offsets[static_cast<Size>(rank)];

                    if (count * 9LL > std::numeric_limits<int>::max() ||
                        offset * 9LL > std::numeric_limits<int>::max())
                    {
                        throw std::runtime_error(
                            "TurbulenceWallTopology - MPI partition-face packing exceeds integer range");
                    }

                    key_counts[static_cast<Size>(rank)] =
                        static_cast<int>(count * 3LL);
                    key_offsets[static_cast<Size>(rank)] =
                        static_cast<int>(offset * 3LL);
                    coordinate_counts[static_cast<Size>(rank)] =
                        static_cast<int>(count * 9LL);
                    coordinate_offsets[static_cast<Size>(rank)] =
                        static_cast<int>(offset * 9LL);
                }

                std::vector<long long> all_keys(
                    static_cast<Size>(total_faces) * 3U);
                std::vector<int> all_material(
                    static_cast<Size>(total_faces));
                std::vector<double> all_coordinates(
                    static_cast<Size>(total_faces) * 9U);

                checkMpi(
                    MPI_Allgatherv(
                        local_keys.data(),
                        local_count * 3,
                        MPI_LONG_LONG,
                        all_keys.data(),
                        key_counts.data(),
                        key_offsets.data(),
                        MPI_LONG_LONG,
                        MPI_COMM_WORLD),
                    "MPI_Allgatherv partition-face keys");

                checkMpi(
                    MPI_Allgatherv(
                        local_material.data(),
                        local_count,
                        MPI_INT,
                        all_material.data(),
                        counts.data(),
                        offsets.data(),
                        MPI_INT,
                        MPI_COMM_WORLD),
                    "MPI_Allgatherv partition-face material classes");

                checkMpi(
                    MPI_Allgatherv(
                        local_coordinates.data(),
                        local_count * 9,
                        MPI_DOUBLE,
                        all_coordinates.data(),
                        coordinate_counts.data(),
                        coordinate_offsets.data(),
                        MPI_DOUBLE,
                        MPI_COMM_WORLD),
                    "MPI_Allgatherv partition-face coordinates");

                struct GatheredFace
                {
                    KeyedTriangle keyed_triangle;
                    bool fluid = false;
                };

                std::vector<GatheredFace> gathered;
                gathered.resize(static_cast<Size>(total_faces));

                for (int i = 0; i < total_faces; ++i)
                {
                    GatheredFace& face = gathered[static_cast<Size>(i)];

                    for (Size k = 0; k < 3; ++k)
                    {
                        face.keyed_triangle.key[k] =
                            all_keys[static_cast<Size>(i) * 3U + k];
                    }

                    std::array<std::array<Real, 3>, 3> vertices{};

                    for (Size corner = 0; corner < 3; ++corner)
                    {
                        for (Size k = 0; k < 3; ++k)
                        {
                            vertices[corner][k] =
                                all_coordinates[
                                    static_cast<Size>(i) * 9U +
                                    corner * 3U + k];
                        }
                    }

                    face.keyed_triangle.triangle.a = vertices[0];
                    face.keyed_triangle.triangle.b = vertices[1];
                    face.keyed_triangle.triangle.c = vertices[2];
                    face.fluid = all_material[static_cast<Size>(i)] != 0;
                }

                std::sort(
                    gathered.begin(),
                    gathered.end(),
                    [](const GatheredFace& lhs, const GatheredFace& rhs)
                    {
                        return lhs.keyed_triangle.key < rhs.keyed_triangle.key;
                    });

                std::vector<KeyedTriangle> interfaces;
                Size begin = 0;

                while (begin < gathered.size())
                {
                    Size end = begin + 1U;

                    while (end < gathered.size() &&
                           sameKey(
                               gathered[begin].keyed_triangle.key,
                               gathered[end].keyed_triangle.key))
                    {
                        ++end;
                    }

                    const Size multiplicity = end - begin;

                    if (multiplicity > 2U)
                    {
                        throw std::runtime_error(
                            "TurbulenceWallTopology - non-manifold global partition face detected");
                    }

                    if (multiplicity == 2U)
                    {
                        const GatheredFace& first = gathered[begin];
                        const GatheredFace& second = gathered[begin + 1U];

                        validateDuplicateGeometry(
                            first.keyed_triangle,
                            second.keyed_triangle);

                        if (first.fluid != second.fluid)
                        {
                            interfaces.push_back(
                                first.fluid
                                    ? first.keyed_triangle
                                    : second.keyed_triangle);
                        }
                    }

                    begin = end;
                }

                deduplicate(interfaces);
                return interfaces;
            }

            static std::vector<KeyedTriangle> gatherTriangles(
                const CBSStateSI& s,
                const std::vector<KeyedTriangle>& local)
            {
                if (!s.mpi_enabled || s.mpi_size <= 1)
                {
                    return local;
                }

                if (local.size()
                    > static_cast<Size>(std::numeric_limits<int>::max()))
                {
                    throw std::runtime_error(
                        "TurbulenceWallTopology - too many local wall triangles for MPI");
                }

                const int local_count = static_cast<int>(local.size());
                std::vector<int> counts(static_cast<Size>(s.mpi_size), 0);

                checkMpi(
                    MPI_Allgather(
                        &local_count,
                        1,
                        MPI_INT,
                        counts.data(),
                        1,
                        MPI_INT,
                        MPI_COMM_WORLD),
                    "MPI_Allgather wall-triangle counts");

                std::vector<int> offsets(static_cast<Size>(s.mpi_size), 0);
                long long total_ll = 0;

                for (Int rank = 0; rank < s.mpi_size; ++rank)
                {
                    offsets[static_cast<Size>(rank)] = static_cast<int>(total_ll);
                    total_ll += counts[static_cast<Size>(rank)];
                }

                if (total_ll > std::numeric_limits<int>::max())
                {
                    throw std::runtime_error(
                        "TurbulenceWallTopology - global wall-triangle count exceeds MPI limit");
                }

                const int total = static_cast<int>(total_ll);

                std::vector<long long> local_keys(
                    static_cast<Size>(local_count) * 3U);
                std::vector<double> local_coordinates(
                    static_cast<Size>(local_count) * 9U);

                for (int i = 0; i < local_count; ++i)
                {
                    const KeyedTriangle& item = local[static_cast<Size>(i)];

                    for (Size k = 0; k < 3; ++k)
                    {
                        local_keys[static_cast<Size>(i) * 3U + k] = item.key[k];
                    }

                    const std::array<std::array<Real, 3>, 3> vertices =
                    {
                        item.triangle.a,
                        item.triangle.b,
                        item.triangle.c
                    };

                    for (Size corner = 0; corner < 3; ++corner)
                    {
                        for (Size k = 0; k < 3; ++k)
                        {
                            local_coordinates[
                                static_cast<Size>(i) * 9U +
                                corner * 3U + k] = vertices[corner][k];
                        }
                    }
                }

                std::vector<int> key_counts(static_cast<Size>(s.mpi_size));
                std::vector<int> key_offsets(static_cast<Size>(s.mpi_size));
                std::vector<int> coordinate_counts(static_cast<Size>(s.mpi_size));
                std::vector<int> coordinate_offsets(static_cast<Size>(s.mpi_size));

                for (Int rank = 0; rank < s.mpi_size; ++rank)
                {
                    const long long count = counts[static_cast<Size>(rank)];
                    const long long offset = offsets[static_cast<Size>(rank)];

                    if (count * 9LL > std::numeric_limits<int>::max() ||
                        offset * 9LL > std::numeric_limits<int>::max())
                    {
                        throw std::runtime_error(
                            "TurbulenceWallTopology - MPI wall-triangle packing exceeds integer range");
                    }

                    key_counts[static_cast<Size>(rank)] = static_cast<int>(count * 3LL);
                    key_offsets[static_cast<Size>(rank)] = static_cast<int>(offset * 3LL);
                    coordinate_counts[static_cast<Size>(rank)] = static_cast<int>(count * 9LL);
                    coordinate_offsets[static_cast<Size>(rank)] = static_cast<int>(offset * 9LL);
                }

                std::vector<long long> all_keys(
                    static_cast<Size>(total) * 3U);
                std::vector<double> all_coordinates(
                    static_cast<Size>(total) * 9U);

                checkMpi(
                    MPI_Allgatherv(
                        local_keys.data(),
                        local_count * 3,
                        MPI_LONG_LONG,
                        all_keys.data(),
                        key_counts.data(),
                        key_offsets.data(),
                        MPI_LONG_LONG,
                        MPI_COMM_WORLD),
                    "MPI_Allgatherv wall-triangle keys");

                checkMpi(
                    MPI_Allgatherv(
                        local_coordinates.data(),
                        local_count * 9,
                        MPI_DOUBLE,
                        all_coordinates.data(),
                        coordinate_counts.data(),
                        coordinate_offsets.data(),
                        MPI_DOUBLE,
                        MPI_COMM_WORLD),
                    "MPI_Allgatherv wall-triangle coordinates");

                std::vector<KeyedTriangle> gathered;
                gathered.resize(static_cast<Size>(total));

                for (int i = 0; i < total; ++i)
                {
                    KeyedTriangle& item = gathered[static_cast<Size>(i)];

                    for (Size k = 0; k < 3; ++k)
                    {
                        item.key[k] = all_keys[static_cast<Size>(i) * 3U + k];
                    }

                    for (Size k = 0; k < 3; ++k)
                    {
                        item.triangle.a[k] =
                            all_coordinates[static_cast<Size>(i) * 9U + k];
                        item.triangle.b[k] =
                            all_coordinates[static_cast<Size>(i) * 9U + 3U + k];
                        item.triangle.c[k] =
                            all_coordinates[static_cast<Size>(i) * 9U + 6U + k];
                    }
                }

                deduplicate(gathered);
                return gathered;
            }
#endif

            static std::vector<KeyedTriangle> collectGlobalKeyed(
                const CBSStateSI& s,
                Size& local_count)
            {
                std::vector<FaceRecord> partition_candidates;
                std::vector<KeyedTriangle> local =
                    collectLocal(s, partition_candidates);

                local_count = local.size();

#ifdef CBS3D_USE_MPI
                if (s.mpi_enabled && s.mpi_size > 1)
                {
                    std::vector<KeyedTriangle> cross_rank =
                        collectCrossRankInterfaces(s, partition_candidates);

                    local.insert(
                        local.end(),
                        cross_rank.begin(),
                        cross_rank.end());

                    deduplicate(local);
                    return gatherTriangles(s, local);
                }
#endif

                return local;
            }

        public:
            // Rebuilds only the SA wall-node flag. Active/inlet flags remain
            // owned by TurbulenceBoundary. Interface nodes are determined from
            // the reconciled material mask, so this remains correct even when
            // no explicit BC 901 face exists in .plt.
            static void reconcileWallNodeClassification(CBSStateSI& s)
            {
                s.sa_wall_node.fill(0);

                for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
                {
                    if (s.sa_active_node(ip) != 0 &&
                        isMaterialInterfaceNode(s, ip))
                    {
                        s.sa_wall_node(ip) = 1;
                    }
                }

                for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
                {
                    const Int bc = s.iside(s.cfg.bsid, ib);

                    if (!isExplicitWallBc(s, bc))
                    {
                        continue;
                    }

                    const Int parent = s.iside(s.cfg.nsidpe, ib);

                    if (parent < 1 || parent > s.cfg.nelem)
                    {
                        throw std::runtime_error(
                            "TurbulenceWallTopology - explicit wall has invalid parent tetrahedron");
                    }

                    if (s.mat_elem(parent) != 0)
                    {
                        continue;
                    }

                    for (Int corner = 1; corner <= s.cfg.nsidp; ++corner)
                    {
                        const Int ip = s.iside(corner, ib);

                        if (s.sa_active_node(ip) != 0)
                        {
                            s.sa_wall_node(ip) = 1;
                        }
                    }
                }
            }

            // Returns the exact fluid turbulence-wall surface visible to this
            // rank. Under MPI this is the globally gathered and deduplicated
            // surface, including material interfaces reconstructed across cuts.
            static std::vector<WallTriangle> collectFluidWallTriangles(
                const CBSStateSI& s)
            {
                Size local_count = 0;
                std::vector<KeyedTriangle> keyed =
                    collectGlobalKeyed(s, local_count);

                std::vector<WallTriangle> triangles;
                triangles.reserve(keyed.size());

                for (const KeyedTriangle& item : keyed)
                {
                    triangles.push_back(item.triangle);
                }

                return triangles;
            }

            // Computes the SA distance using the corrected CHT topology. The
            // BVH remains exact; only the wall-surface inventory is changed.
            static void computeWallDistance(
                CBSStateSI& s,
                WallDistanceStats& stats)
            {
                const auto gather_start = std::chrono::steady_clock::now();

                Size local_count = 0;
                std::vector<KeyedTriangle> keyed =
                    collectGlobalKeyed(s, local_count);

                stats.local_wall_triangles = local_count;
                stats.global_wall_triangles = keyed.size();
                stats.searched_wall_triangles = keyed.size();
                stats.gather_seconds =
                    std::chrono::duration<Real>(
                        std::chrono::steady_clock::now() - gather_start).count();

                if (keyed.empty())
                {
                    throw std::runtime_error(
                        "TurbulenceWallTopology - no fluid turbulence-wall triangles found");
                }

                std::vector<WallTriangle> triangles;
                triangles.reserve(keyed.size());

                for (const KeyedTriangle& item : keyed)
                {
                    triangles.push_back(item.triangle);
                }

                const auto build_start = std::chrono::steady_clock::now();
                WallTriangleBVH bvh;
                bvh.build(triangles);

                stats.bvh_nodes = bvh.nodeCount();
                stats.bvh_depth = bvh.depth();
                stats.build_seconds =
                    std::chrono::duration<Real>(
                        std::chrono::steady_clock::now() - build_start).count();

                if (s.wall_distance.size() != s.cfg.npoin)
                {
                    throw std::runtime_error(
                        "TurbulenceWallTopology - wall_distance array has wrong size");
                }

                s.wall_distance.fill(0.0);
                const auto query_start = std::chrono::steady_clock::now();

#ifdef CBS3D_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
                for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
                {
                    if (s.sa_active_node(ip) == 0 ||
                        s.sa_wall_node(ip) != 0)
                    {
                        s.wall_distance(ip) = 0.0;
                        continue;
                    }

                    const std::array<Real, 3> p = point(s, ip);
                    const Real distance_squared =
                        bvh.nearestDistanceSquared(
                            p,
                            std::numeric_limits<Real>::max());

                    if (!(distance_squared >= 0.0) ||
                        !std::isfinite(distance_squared))
                    {
                        throw std::runtime_error(
                            "TurbulenceWallTopology - invalid nearest-wall distance");
                    }

                    s.wall_distance(ip) = std::sqrt(distance_squared);
                }

                stats.query_seconds =
                    std::chrono::duration<Real>(
                        std::chrono::steady_clock::now() - query_start).count();

                Real local_min = std::numeric_limits<Real>::max();
                Real local_max = 0.0;
                long long local_queries = 0;

                for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
                {
                    if (s.sa_active_node(ip) == 0)
                    {
                        continue;
                    }

                    const Real distance = s.wall_distance(ip);

                    if (!(distance >= 0.0) || !std::isfinite(distance))
                    {
                        throw std::runtime_error(
                            "TurbulenceWallTopology - non-finite stored wall distance");
                    }

                    local_min = std::min(local_min, distance);
                    local_max = std::max(local_max, distance);
                    ++local_queries;
                }

                if (local_queries == 0)
                {
                    local_min = 0.0;
                }

#ifdef CBS3D_USE_MPI
                if (s.mpi_enabled && s.mpi_size > 1)
                {
                    Real global_min = 0.0;
                    Real global_max = 0.0;
                    long long global_queries = 0;
                    Real global_query_seconds = 0.0;

                    checkMpi(
                        MPI_Allreduce(
                            &local_min,
                            &global_min,
                            1,
                            MPI_DOUBLE,
                            MPI_MIN,
                            MPI_COMM_WORLD),
                        "MPI_Allreduce wall-distance minimum");

                    checkMpi(
                        MPI_Allreduce(
                            &local_max,
                            &global_max,
                            1,
                            MPI_DOUBLE,
                            MPI_MAX,
                            MPI_COMM_WORLD),
                        "MPI_Allreduce wall-distance maximum");

                    checkMpi(
                        MPI_Allreduce(
                            &local_queries,
                            &global_queries,
                            1,
                            MPI_LONG_LONG,
                            MPI_SUM,
                            MPI_COMM_WORLD),
                        "MPI_Allreduce wall-distance query count");

                    checkMpi(
                        MPI_Allreduce(
                            &stats.query_seconds,
                            &global_query_seconds,
                            1,
                            MPI_DOUBLE,
                            MPI_MAX,
                            MPI_COMM_WORLD),
                        "MPI_Allreduce wall-distance query time");

                    stats.min_distance = global_min;
                    stats.max_distance = global_max;
                    stats.queried_nodes =
                        global_queries > std::numeric_limits<Int>::max()
                            ? std::numeric_limits<Int>::max()
                            : static_cast<Int>(global_queries);
                    stats.query_seconds = global_query_seconds;

                    return;
                }
#endif

                stats.min_distance = local_min;
                stats.max_distance = local_max;
                stats.queried_nodes =
                    local_queries > std::numeric_limits<Int>::max()
                        ? std::numeric_limits<Int>::max()
                        : static_cast<Int>(local_queries);
            }
        };
    }
}
