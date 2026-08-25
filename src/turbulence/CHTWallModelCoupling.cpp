//=============================================================================
// CBS3D++_SI
//
// Conformal CHT wall-model coupling.
//
// This translation unit deliberately owns all mutable cache state.  Keeping the
// implementation out of the public header avoids ODR/compiler problems and
// makes the MPI collectives explicit in one place.
//=============================================================================

#include "cbs/turbulence/CHTWallModelCoupling.hpp"

#include "cbs/parallel/HaloExchange.hpp"
#include "cbs/turbulence/ThermalWallTreatment.hpp"
#include "cbs/turbulence/WallTreatment.hpp"

#ifdef CBS3D_USE_MPI
#include <mpi.h>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace cbs
{
    namespace turbulence
    {
        namespace
        {
            struct ModelFace
            {
                Int parent = 0;
                Int local_face = 0;
                Int sample_node = 0;
                std::array<Int, 3> nodes = {0, 0, 0};
                std::array<long long, 3> global_key = {0, 0, 0};
                std::array<Real, 3> normal = {0.0, 0.0, 0.0};
                Real area = 0.0;
                Real sample_height = 0.0;
            };

            struct CandidateFace
            {
                Int parent = 0;
                Int local_face = 0;
                bool fluid = false;
                std::array<Int, 3> local_key = {0, 0, 0};
                std::array<long long, 3> global_key = {0, 0, 0};
            };

            struct CouplingCache
            {
                const CBSStateSI* state = nullptr;
                Int npoin = 0;
                Int nelem = 0;
                bool ready = false;
                long long global_faces = 0;
                std::vector<ModelFace> faces;
                std::vector<Int> wall_nodes;
                std::vector<std::array<Real, 9>> normal_projector;
                std::vector<Real> bulk_conductivity;
            };

            CouplingCache cache;

            bool environment_flag(const char* name)
            {
                const char* value = std::getenv(name);
                return value != nullptr && value[0] != '\0' &&
                    std::string(value) != "0";
            }

            Real environment_real(
                const char* name,
                const Real default_value)
            {
                const char* value = std::getenv(name);

                if (value == nullptr || value[0] == '\0')
                {
                    return default_value;
                }

                char* end = nullptr;
                const Real parsed = std::strtod(value, &end);

                if (end == value || *end != '\0' || !std::isfinite(parsed))
                {
                    throw std::runtime_error(
                        std::string(
                            "CHTWallModelCoupling - invalid real environment value for ")
                        + name);
                }

                return parsed;
            }

            bool requested()
            {
                return environment_flag("CBS3D_CHT_WALL_TREATMENT");
            }

            void check_mpi(const int code, const char* operation)
            {
#ifdef CBS3D_USE_MPI
                if (code != MPI_SUCCESS)
                {
                    throw std::runtime_error(
                        std::string("CHTWallModelCoupling MPI failure in ")
                        + operation);
                }
#else
                (void)code;
                (void)operation;
#endif
            }

            bool is_material_interface_node(
                const CBSStateSI& s,
                const Int ip)
            {
                if (ip < 1 || ip > s.cfg.npoin)
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - node outside local mesh");
                }

                const Int required =
                    CBSStateSI::node_touches_fluid |
                    CBSStateSI::node_touches_solid;

                return s.node_material_mask(ip) == required;
            }

            long long global_node_id(
                const CBSStateSI& s,
                const Int ip)
            {
                if (ip < 1 || ip > s.cfg.npoin)
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - node outside local mesh");
                }

#ifdef CBS3D_USE_MPI
                if (s.mpi_enabled && s.mpi_size > 1)
                {
                    if (static_cast<Size>(ip) >= s.local_to_global_node.size())
                    {
                        throw std::runtime_error(
                            "CHTWallModelCoupling - missing local-to-global node map");
                    }

                    const Int gid =
                        s.local_to_global_node[static_cast<Size>(ip)];

                    if (gid < 1)
                    {
                        throw std::runtime_error(
                            "CHTWallModelCoupling - invalid global node id");
                    }

                    return static_cast<long long>(gid);
                }
#endif

                return static_cast<long long>(ip);
            }

            Int dnkdx_index(
                const CBSStateSI& s,
                const Int ie,
                const Int dim,
                const Int local_node)
            {
                return (ie - 1) * s.cfg.ndim * s.cfg.nep
                    + (dim - 1) * s.cfg.nep
                    + local_node;
            }

            Real grad(
                const CBSStateSI& s,
                const Int ie,
                const Int dim,
                const Int local_node)
            {
                return s.dNkdx(dnkdx_index(s, ie, dim, local_node));
            }

            void validate_configuration(const CBSStateSI& s)
            {
                if (!requested())
                {
                    return;
                }

                if (environment_flag("CBS3D_SA_WALL_TREATMENT"))
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - do not enable CBS3D_SA_WALL_TREATMENT and CBS3D_CHT_WALL_TREATMENT together");
                }

                if (s.cfg.turbulence_on < 1 || s.cfg.turbulence_model != 0)
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - CHT wall treatment requires standard Spalart-Allmaras");
                }

                if (s.cfg.temp_calc < 1)
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - CHT wall treatment requires the energy equation");
                }

                if (s.cfg.dimensional_mode < 1 ||
                    s.cfg.material_properties_enabled < 1)
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - dimensional material properties are required");
                }

                if (s.cfg.turbulent_thermal_diffusivity_on < 1)
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - turbulent thermal diffusivity must be enabled for the bulk fluid");
                }

                if (!(s.cfg.sa_prandtl_t > 0.0) ||
                    !std::isfinite(s.cfg.sa_prandtl_t))
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - sa_prandtl_t must be positive and finite");
                }

                if (s.cfg.ndim != 3 || s.cfg.nep != 4 ||
                    s.cfg.nsid != 4 || s.cfg.nsidp != 3)
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - requires 3-D P1 TET4/TRI3 topology");
                }

                int local_fluid = 0;
                int local_solid = 0;

                for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
                {
                    if (s.mat_elem(ie) == 0)
                    {
                        local_fluid = 1;
                    }
                    else
                    {
                        local_solid = 1;
                    }
                }

#ifdef CBS3D_USE_MPI
                if (s.mpi_enabled && s.mpi_size > 1)
                {
                    int global_fluid = 0;
                    int global_solid = 0;

                    check_mpi(
                        MPI_Allreduce(
                            &local_fluid,
                            &global_fluid,
                            1,
                            MPI_INT,
                            MPI_MAX,
                            MPI_COMM_WORLD),
                        "MPI_Allreduce fluid-domain audit");

                    check_mpi(
                        MPI_Allreduce(
                            &local_solid,
                            &global_solid,
                            1,
                            MPI_INT,
                            MPI_MAX,
                            MPI_COMM_WORLD),
                        "MPI_Allreduce solid-domain audit");

                    local_fluid = global_fluid;
                    local_solid = global_solid;
                }
#endif

                if (local_fluid == 0 || local_solid == 0)
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - both fluid and solid materials are required");
                }
            }

            CandidateFace candidate_face(
                const CBSStateSI& s,
                const Int ie,
                const Int face)
            {
                CandidateFace result;
                result.parent = ie;
                result.local_face = face;
                result.fluid = s.mat_elem(ie) == 0;

                for (Int corner = 1; corner <= 3; ++corner)
                {
                    const Int ip =
                        s.intma(s.ippn1(face, corner), ie);

                    result.local_key[static_cast<Size>(corner - 1)] = ip;
                    result.global_key[static_cast<Size>(corner - 1)] =
                        global_node_id(s, ip);
                }

                std::sort(
                    result.local_key.begin(),
                    result.local_key.end());

                std::sort(
                    result.global_key.begin(),
                    result.global_key.end());

                return result;
            }

            bool can_be_interface(
                const CBSStateSI& s,
                const CandidateFace& face)
            {
                for (const Int ip : face.local_key)
                {
                    if (!is_material_interface_node(s, ip))
                    {
                        return false;
                    }
                }

                return true;
            }

            ModelFace make_model_face(
                const CBSStateSI& s,
                const CandidateFace& candidate)
            {
                if (!candidate.fluid)
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - attempted to model a solid-side face");
                }

                const Int ie = candidate.parent;
                const Int face = candidate.local_face;
                const Real area = s.annxf(4, face, ie);

                if (!(area > 0.0) || !std::isfinite(area) ||
                    !(s.detJ(ie) > 0.0) || !std::isfinite(s.detJ(ie)))
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - invalid interface face geometry");
                }

                // A conformal material interface is internal to the complete CHT
                // mesh.  If it also appears in fedge, the same physical face has
                // been classified simultaneously as an exterior boundary and a
                // material interface, which is a mesh/BC error rather than a wall
                // model choice.
                if (s.fedge(face, ie) != 0)
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - reconstructed material interface is also marked as an exterior fedge");
                }

                ModelFace result;
                result.parent = ie;
                result.local_face = face;
                result.sample_node = s.intma(face, ie);
                result.global_key = candidate.global_key;
                result.area = area;

                // TET4 local face f is opposite local node f.  The exact
                // perpendicular altitude is h = 3V/A = detJ/(2A).
                result.sample_height = 0.5 * s.detJ(ie) / area;

                if (!(result.sample_height > 0.0) ||
                    !std::isfinite(result.sample_height))
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - invalid first-fluid-tet wall height");
                }

                Real normal_norm2 = 0.0;

                for (Int dim = 1; dim <= 3; ++dim)
                {
                    const Real component =
                        s.annxf(dim, face, ie) / area;

                    result.normal[static_cast<Size>(dim - 1)] = component;
                    normal_norm2 += component * component;
                }

                const Real normal_norm = std::sqrt(normal_norm2);

                if (!(normal_norm > 0.0) || !std::isfinite(normal_norm))
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - invalid fluid-side interface normal");
                }

                for (Real& component : result.normal)
                {
                    component /= normal_norm;
                }

                for (Int corner = 1; corner <= 3; ++corner)
                {
                    const Int ip =
                        s.intma(s.ippn1(face, corner), ie);

                    result.nodes[static_cast<Size>(corner - 1)] = ip;

                    if (ip == result.sample_node)
                    {
                        throw std::runtime_error(
                            "CHTWallModelCoupling - opposite sample node lies on interface face");
                    }
                }

                return result;
            }

            void collect_local_faces(
                const CBSStateSI& s,
                std::vector<ModelFace>& accepted,
                std::vector<CandidateFace>& partition_candidates)
            {
#ifndef CBS3D_USE_MPI
                (void)partition_candidates;
#endif
                std::vector<CandidateFace> candidates;

                for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
                {
                    for (Int face = 1; face <= s.cfg.nsid; ++face)
                    {
                        CandidateFace candidate =
                            candidate_face(s, ie, face);

                        if (can_be_interface(s, candidate))
                        {
                            candidates.push_back(candidate);
                        }
                    }
                }

                std::sort(
                    candidates.begin(),
                    candidates.end(),
                    [](const CandidateFace& lhs, const CandidateFace& rhs)
                    {
                        return lhs.local_key < rhs.local_key;
                    });

                Size begin = 0;

                while (begin < candidates.size())
                {
                    Size end = begin + 1U;

                    while (end < candidates.size() &&
                           candidates[end].local_key ==
                               candidates[begin].local_key)
                    {
                        ++end;
                    }

                    const Size multiplicity = end - begin;

                    if (multiplicity > 2U)
                    {
                        throw std::runtime_error(
                            "CHTWallModelCoupling - non-manifold local tetrahedral face");
                    }

                    if (multiplicity == 2U)
                    {
                        const CandidateFace& a = candidates[begin];
                        const CandidateFace& b = candidates[begin + 1U];

                        if (a.fluid != b.fluid)
                        {
                            accepted.push_back(
                                make_model_face(s, a.fluid ? a : b));
                        }
                    }
                    else
                    {
#ifdef CBS3D_USE_MPI
                        if (s.mpi_enabled && s.mpi_size > 1)
                        {
                            partition_candidates.push_back(candidates[begin]);
                        }
#endif
                    }

                    begin = end;
                }
            }

            void collect_cross_rank_faces(
                const CBSStateSI& s,
                const std::vector<CandidateFace>& local_candidates,
                std::vector<ModelFace>& accepted)
            {
#ifdef CBS3D_USE_MPI
                if (!s.mpi_enabled || s.mpi_size <= 1)
                {
                    return;
                }

                if (local_candidates.size() >
                    static_cast<Size>(std::numeric_limits<int>::max()))
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - too many partition interface candidates");
                }

                const int local_count =
                    static_cast<int>(local_candidates.size());

                std::vector<int> counts(
                    static_cast<Size>(s.mpi_size),
                    0);

                check_mpi(
                    MPI_Allgather(
                        &local_count,
                        1,
                        MPI_INT,
                        counts.data(),
                        1,
                        MPI_INT,
                        MPI_COMM_WORLD),
                    "MPI_Allgather interface candidate counts");

                std::vector<int> offsets(
                    static_cast<Size>(s.mpi_size),
                    0);

                long long total_ll = 0;

                for (Int rank = 0; rank < s.mpi_size; ++rank)
                {
                    if (total_ll > std::numeric_limits<int>::max())
                    {
                        throw std::runtime_error(
                            "CHTWallModelCoupling - global candidate count exceeds MPI integer range");
                    }

                    offsets[static_cast<Size>(rank)] =
                        static_cast<int>(total_ll);
                    total_ll += counts[static_cast<Size>(rank)];
                }

                if (total_ll > std::numeric_limits<int>::max())
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - global candidate count exceeds MPI integer range");
                }

                const int total = static_cast<int>(total_ll);

                std::vector<long long> send_keys(
                    static_cast<Size>(local_count) * 3U);

                std::vector<int> send_material(
                    static_cast<Size>(local_count));

                for (int i = 0; i < local_count; ++i)
                {
                    const CandidateFace& face =
                        local_candidates[static_cast<Size>(i)];

                    send_material[static_cast<Size>(i)] =
                        face.fluid ? 1 : 0;

                    for (Size k = 0; k < 3; ++k)
                    {
                        send_keys[static_cast<Size>(i) * 3U + k] =
                            face.global_key[k];
                    }
                }

                std::vector<int> counts3(
                    static_cast<Size>(s.mpi_size));

                std::vector<int> offsets3(
                    static_cast<Size>(s.mpi_size));

                for (Int rank = 0; rank < s.mpi_size; ++rank)
                {
                    if (counts[static_cast<Size>(rank)] >
                        std::numeric_limits<int>::max() / 3)
                    {
                        throw std::runtime_error(
                            "CHTWallModelCoupling - interface candidate key count exceeds MPI integer range");
                    }

                    counts3[static_cast<Size>(rank)] =
                        3 * counts[static_cast<Size>(rank)];
                    offsets3[static_cast<Size>(rank)] =
                        3 * offsets[static_cast<Size>(rank)];
                }

                std::vector<long long> all_keys(
                    static_cast<Size>(total) * 3U);

                std::vector<int> all_material(
                    static_cast<Size>(total));

                check_mpi(
                    MPI_Allgatherv(
                        send_keys.data(),
                        local_count * 3,
                        MPI_LONG_LONG,
                        all_keys.data(),
                        counts3.data(),
                        offsets3.data(),
                        MPI_LONG_LONG,
                        MPI_COMM_WORLD),
                    "MPI_Allgatherv interface candidate keys");

                check_mpi(
                    MPI_Allgatherv(
                        send_material.data(),
                        local_count,
                        MPI_INT,
                        all_material.data(),
                        counts.data(),
                        offsets.data(),
                        MPI_INT,
                        MPI_COMM_WORLD),
                    "MPI_Allgatherv interface material classes");

                struct Gathered
                {
                    std::array<long long, 3> key = {0, 0, 0};
                    bool fluid = false;
                };

                std::vector<Gathered> gathered(
                    static_cast<Size>(total));

                for (int i = 0; i < total; ++i)
                {
                    Gathered& item =
                        gathered[static_cast<Size>(i)];

                    item.fluid =
                        all_material[static_cast<Size>(i)] != 0;

                    for (Size k = 0; k < 3; ++k)
                    {
                        item.key[k] =
                            all_keys[static_cast<Size>(i) * 3U + k];
                    }
                }

                std::sort(
                    gathered.begin(),
                    gathered.end(),
                    [](const Gathered& lhs, const Gathered& rhs)
                    {
                        return lhs.key < rhs.key;
                    });

                std::set<std::array<long long, 3>> interface_keys;
                Size begin = 0;

                while (begin < gathered.size())
                {
                    Size end = begin + 1U;

                    while (end < gathered.size() &&
                           gathered[end].key == gathered[begin].key)
                    {
                        ++end;
                    }

                    const Size multiplicity = end - begin;

                    if (multiplicity > 2U)
                    {
                        throw std::runtime_error(
                            "CHTWallModelCoupling - non-manifold cross-rank tetrahedral face");
                    }

                    if (multiplicity == 2U &&
                        gathered[begin].fluid !=
                            gathered[begin + 1U].fluid)
                    {
                        interface_keys.insert(gathered[begin].key);
                    }

                    begin = end;
                }

                for (const CandidateFace& candidate : local_candidates)
                {
                    if (candidate.fluid &&
                        interface_keys.find(candidate.global_key) !=
                            interface_keys.end())
                    {
                        accepted.push_back(
                            make_model_face(s, candidate));
                    }
                }
#else
                (void)s;
                (void)local_candidates;
                (void)accepted;
#endif
            }

            std::array<Real, 9> projector_from_gram(
                const std::array<Real, 9>& gram)
            {
                std::array<Real, 9> projector{};
                Real frobenius2 = 0.0;

                for (const Real value : gram)
                {
                    frobenius2 += value * value;
                }

                const Real frobenius = std::sqrt(frobenius2);

                if (!(frobenius > 0.0))
                {
                    return projector;
                }

                const Real tolerance = 1.0e-12 * frobenius;
                std::array<std::array<Real, 3>, 3> basis{};
                Int rank = 0;

                for (Int column = 0; column < 3; ++column)
                {
                    std::array<Real, 3> vector =
                    {
                        gram[static_cast<Size>(column)],
                        gram[static_cast<Size>(3 + column)],
                        gram[static_cast<Size>(6 + column)]
                    };

                    for (Int q = 0; q < rank; ++q)
                    {
                        const Real dot =
                            vector[0] * basis[static_cast<Size>(q)][0] +
                            vector[1] * basis[static_cast<Size>(q)][1] +
                            vector[2] * basis[static_cast<Size>(q)][2];

                        for (Int k = 0; k < 3; ++k)
                        {
                            vector[static_cast<Size>(k)] -=
                                dot * basis[static_cast<Size>(q)][static_cast<Size>(k)];
                        }
                    }

                    const Real norm = std::sqrt(
                        vector[0] * vector[0] +
                        vector[1] * vector[1] +
                        vector[2] * vector[2]);

                    if (norm <= tolerance)
                    {
                        continue;
                    }

                    for (Int k = 0; k < 3; ++k)
                    {
                        basis[static_cast<Size>(rank)][static_cast<Size>(k)] =
                            vector[static_cast<Size>(k)] / norm;
                    }

                    ++rank;
                }

                for (Int q = 0; q < rank; ++q)
                {
                    for (Int i = 0; i < 3; ++i)
                    {
                        for (Int j = 0; j < 3; ++j)
                        {
                            projector[static_cast<Size>(3 * i + j)] +=
                                basis[static_cast<Size>(q)][static_cast<Size>(i)] *
                                basis[static_cast<Size>(q)][static_cast<Size>(j)];
                        }
                    }
                }

                return projector;
            }

            bool nonzero_projector(
                const std::array<Real, 9>& projector)
            {
                return projector[0] + projector[4] + projector[8] > 0.5;
            }

            void build_cache(const CBSStateSI& s)
            {
                validate_configuration(s);

                cache = CouplingCache{};
                cache.state = &s;
                cache.npoin = s.cfg.npoin;
                cache.nelem = s.cfg.nelem;

                if (!requested())
                {
                    cache.ready = true;
                    return;
                }

                std::vector<CandidateFace> partition_candidates;
                collect_local_faces(
                    s,
                    cache.faces,
                    partition_candidates);

                collect_cross_rank_faces(
                    s,
                    partition_candidates,
                    cache.faces);

                std::sort(
                    cache.faces.begin(),
                    cache.faces.end(),
                    [](const ModelFace& lhs, const ModelFace& rhs)
                    {
                        return lhs.global_key < rhs.global_key;
                    });

                if (std::adjacent_find(
                        cache.faces.begin(),
                        cache.faces.end(),
                        [](const ModelFace& lhs, const ModelFace& rhs)
                        {
                            return lhs.global_key == rhs.global_key;
                        }) != cache.faces.end())
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - duplicate local CHT interface face");
                }

                const long long local_faces =
                    static_cast<long long>(cache.faces.size());

                cache.global_faces = local_faces;

#ifdef CBS3D_USE_MPI
                if (s.mpi_enabled && s.mpi_size > 1)
                {
                    check_mpi(
                        MPI_Allreduce(
                            &local_faces,
                            &cache.global_faces,
                            1,
                            MPI_LONG_LONG,
                            MPI_SUM,
                            MPI_COMM_WORLD),
                        "MPI_Allreduce CHT wall-face count");
                }
#endif

                if (cache.global_faces <= 0)
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - no conformal fluid-solid interface faces found globally");
                }

                Array2D<Real> gram;
                gram.resize(9, s.cfg.npoin);
                gram.fill(0.0);

                for (const ModelFace& face : cache.faces)
                {
                    for (const Int ip : face.nodes)
                    {
                        for (Int i = 0; i < 3; ++i)
                        {
                            for (Int j = 0; j < 3; ++j)
                            {
                                gram(3 * i + j + 1, ip) +=
                                    face.area *
                                    face.normal[static_cast<Size>(i)] *
                                    face.normal[static_cast<Size>(j)];
                            }
                        }
                    }
                }

#ifdef CBS3D_USE_MPI
                if (s.mpi_enabled && s.mpi_size > 1)
                {
                    HaloExchange::sumGhostContributionsToOwners(
                        gram,
                        s.partition_metadata,
                        MPI_COMM_WORLD);

                    HaloExchange::broadcastOwnedToGhosts(
                        gram,
                        s.partition_metadata,
                        MPI_COMM_WORLD);
                }
#endif

                cache.normal_projector.resize(
                    static_cast<Size>(s.cfg.npoin) + 1U);

                for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
                {
                    std::array<Real, 9> local_gram{};

                    for (Int entry = 0; entry < 9; ++entry)
                    {
                        local_gram[static_cast<Size>(entry)] =
                            gram(entry + 1, ip);
                    }

                    cache.normal_projector[static_cast<Size>(ip)] =
                        projector_from_gram(local_gram);

                    if (nonzero_projector(
                            cache.normal_projector[
                                static_cast<Size>(ip)]))
                    {
                        cache.wall_nodes.push_back(ip);
                    }
                }

                const long long local_wall_node_occurrences =
                    static_cast<long long>(cache.wall_nodes.size());

                long long global_wall_node_occurrences =
                    local_wall_node_occurrences;

#ifdef CBS3D_USE_MPI
                if (s.mpi_enabled && s.mpi_size > 1)
                {
                    check_mpi(
                        MPI_Allreduce(
                            &local_wall_node_occurrences,
                            &global_wall_node_occurrences,
                            1,
                            MPI_LONG_LONG,
                            MPI_SUM,
                            MPI_COMM_WORLD),
                        "MPI_Allreduce CHT wall-node audit");
                }
#endif

                if (global_wall_node_occurrences <= 0)
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - no interface wall-node projector survived globally");
                }

                cache.bulk_conductivity.assign(
                    static_cast<Size>(s.cfg.nelem) + 1U,
                    0.0);

                cache.ready = true;

                if (!s.mpi_enabled || s.mpi_rank == 0)
                {
                    WallTreatmentOptions options;
                    options.kappa = environment_real(
                        "CBS3D_WALL_KAPPA",
                        options.kappa);
                    options.log_intercept = environment_real(
                        "CBS3D_WALL_B",
                        options.log_intercept);

                    std::cout
                        << "  SA CHT wall treatment: ON\n"
                        << "    momentum model      : Spalding smooth wall\n"
                        << "    thermal model       : Kader continuous T+\n"
                        << "    kappa / B           : "
                        << options.kappa << " / "
                        << options.log_intercept << "\n"
                        << "    global CHT faces    : "
                        << cache.global_faces << "\n"
                        << "    coupling            : conformal fluid-solid topology\n";
                }
            }

            void ensure_cache(const CBSStateSI& s)
            {
                if (!cache.ready ||
                    cache.state != &s ||
                    cache.npoin != s.cfg.npoin ||
                    cache.nelem != s.cfg.nelem)
                {
                    build_cache(s);
                }
            }

            bool owned_node(const CBSStateSI& s, const Int ip)
            {
                if (!s.mpi_enabled)
                {
                    return true;
                }

                return std::find(
                    s.owned_nodes.begin(),
                    s.owned_nodes.end(),
                    ip) != s.owned_nodes.end();
            }

            std::array<Real, 3> project_tangent(
                const std::array<Real, 3>& velocity,
                const std::array<Real, 9>& projector)
            {
                std::array<Real, 3> normal_component{};

                for (Int i = 0; i < 3; ++i)
                {
                    for (Int j = 0; j < 3; ++j)
                    {
                        normal_component[static_cast<Size>(i)] +=
                            projector[static_cast<Size>(3 * i + j)] *
                            velocity[static_cast<Size>(j)];
                    }
                }

                return
                {
                    velocity[0] - normal_component[0],
                    velocity[1] - normal_component[1],
                    velocity[2] - normal_component[2]
                };
            }

            WallTreatmentOptions wall_options()
            {
                WallTreatmentOptions options;
                options.kappa = environment_real(
                    "CBS3D_WALL_KAPPA",
                    options.kappa);
                options.log_intercept = environment_real(
                    "CBS3D_WALL_B",
                    options.log_intercept);
                return options;
            }

            void molecular_properties(
                const CBSStateSI& s,
                const Int ie,
                Real& density,
                Real& molecular_nu)
            {
                density = s.rho_e(ie);
                const Real mu = s.mu_e(ie);

                if (!(density > 0.0) || !std::isfinite(density) ||
                    !(mu > 0.0) || !std::isfinite(mu))
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - invalid fluid molecular properties");
                }

                molecular_nu = mu / density;
            }

            Real bulk_thermal_conductivity(
                const CBSStateSI& s,
                const Int ie)
            {
                if (s.mat_elem(ie) != 0)
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - bulk fluid conductivity requested on solid");
                }

                Real conductivity = s.k_e(ie);

                if (s.cfg.turbulence_on > 0 &&
                    s.cfg.turbulent_thermal_diffusivity_on > 0)
                {
                    if (s.nu_t_e(ie) < 0.0 ||
                        !std::isfinite(s.nu_t_e(ie)))
                    {
                        throw std::runtime_error(
                            "CHTWallModelCoupling - invalid element eddy viscosity");
                    }

                    conductivity +=
                        s.rho_cp_e(ie) *
                        s.nu_t_e(ie) /
                        s.cfg.sa_prandtl_t;
                }

                if (!(conductivity > 0.0) ||
                    !std::isfinite(conductivity))
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - invalid bulk fluid conductivity");
                }

                return conductivity;
            }

            WallTreatmentResult evaluate_momentum_wall(
                const CBSStateSI& s,
                const ModelFace& face,
                const bool current_velocity)
            {
                Real density = 0.0;
                Real molecular_nu = 0.0;
                molecular_properties(
                    s,
                    face.parent,
                    density,
                    molecular_nu);

                const std::array<Real, 3> velocity =
                {
                    current_velocity
                        ? s.unkno(1, face.sample_node)
                        : s.unkn1(1, face.sample_node),
                    current_velocity
                        ? s.unkno(2, face.sample_node)
                        : s.unkn1(2, face.sample_node),
                    current_velocity
                        ? s.unkno(3, face.sample_node)
                        : s.unkn1(3, face.sample_node)
                };

                return WallTreatment::evaluateSpalding(
                    velocity,
                    face.normal,
                    face.sample_height,
                    density,
                    molecular_nu,
                    wall_options());
            }

            ThermalWallTreatmentResult evaluate_thermal_wall(
                const CBSStateSI& s,
                const ModelFace& face,
                const bool current_velocity)
            {
                const WallTreatmentResult momentum =
                    evaluate_momentum_wall(
                        s,
                        face,
                        current_velocity);

                const Int ie = face.parent;
                const Real density = s.rho_e(ie);

                if (!(density > 0.0) || !std::isfinite(density) ||
                    !(s.rho_cp_e(ie) > 0.0) ||
                    !std::isfinite(s.rho_cp_e(ie)))
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - invalid fluid rho/rhoCp for thermal wall law");
                }

                const Real cp = s.rho_cp_e(ie) / density;

                return ThermalWallTreatment::evaluateKader(
                    momentum.friction_velocity,
                    face.sample_height,
                    density,
                    cp,
                    s.mu_e(ie),
                    s.k_e(ie));
            }

            void add_characteristic_boundary_term(
                CBSStateSI& s,
                const ModelFace& face)
            {
                const Int ie = face.parent;
                Real flux_gradient[4][4] = {};
                Real mean_velocity[4] = {};

                // This reproduces only the boundary half of
                // MomentumAssembly::step1_characteristic_correction for a face
                // that is a boundary of the FLUID subdomain but is internal to
                // the complete CHT mesh.  The volume term already exists.
                for (Int i = 1; i <= 3; ++i)
                {
                    for (Int j = 1; j <= 3; ++j)
                    {
                        for (Int a = 1; a <= s.cfg.nep; ++a)
                        {
                            const Int ip = s.intma(a, ie);
                            flux_gradient[i][j] +=
                                grad(s, ie, j, a) *
                                s.unkn1(j, ip) *
                                s.unkn1(i, ip);
                        }
                    }
                }

                for (Int k = 1; k <= 3; ++k)
                {
                    for (Int a = 1; a <= s.cfg.nep; ++a)
                    {
                        mean_velocity[k] +=
                            s.unkn1(k, s.intma(a, ie));
                    }

                    mean_velocity[k] /=
                        static_cast<Real>(s.cfg.nep);
                }

                const Real coefficient =
                    0.5 * s.delte(ie) * s.cfg.fcon[4];

                if (!(coefficient >= 0.0) ||
                    !std::isfinite(coefficient))
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - invalid characteristic wall coefficient");
                }

                for (Int i = 1; i <= 3; ++i)
                {
                    Real divergence = 0.0;

                    for (Int j = 1; j <= 3; ++j)
                    {
                        divergence += flux_gradient[i][j];
                    }

                    Real face_value = 0.0;

                    for (Int k = 1; k <= 3; ++k)
                    {
                        face_value +=
                            mean_velocity[k] *
                            divergence *
                            s.annxf(k, face.local_face, ie) *
                            coefficient;
                    }

                    for (const Int ip : face.nodes)
                    {
                        s.rhs(i, ip) += face_value;
                    }
                }
            }
        }

        bool CHTWallModelCoupling::enabled(const CBSStateSI& s)
        {
            if (!requested())
            {
                return false;
            }

            ensure_cache(s);
            return true;
        }

        long long CHTWallModelCoupling::globalWallFaceCount(
            const CBSStateSI& s)
        {
            return enabled(s) ? cache.global_faces : 0;
        }

        bool CHTWallModelCoupling::isModelWallNode(
            const CBSStateSI& s,
            const Int ip)
        {
            if (!enabled(s) || ip < 1 || ip > s.cfg.npoin)
            {
                return false;
            }

            return nonzero_projector(
                cache.normal_projector[static_cast<Size>(ip)]);
        }

        CHTCapturedWallVelocity CHTWallModelCoupling::captureVelocity(
            const CBSStateSI& s,
            const bool owned_only)
        {
            CHTCapturedWallVelocity captured;

            if (!enabled(s))
            {
                return captured;
            }

            captured.nodes.reserve(cache.wall_nodes.size());
            captured.values.reserve(cache.wall_nodes.size());

            for (const Int ip : cache.wall_nodes)
            {
                if (owned_only && !owned_node(s, ip))
                {
                    continue;
                }

                captured.nodes.push_back(ip);
                captured.values.push_back(
                    {
                        s.unkno(1, ip),
                        s.unkno(2, ip),
                        s.unkno(3, ip)
                    });
            }

            return captured;
        }

        void CHTWallModelCoupling::restoreTangentialAndEnforceImpermeability(
            CBSStateSI& s,
            const CHTCapturedWallVelocity& captured)
        {
            if (!enabled(s))
            {
                return;
            }

            if (captured.nodes.size() != captured.values.size())
            {
                throw std::runtime_error(
                    "CHTWallModelCoupling - inconsistent captured velocity arrays");
            }

            for (Size entry = 0; entry < captured.nodes.size(); ++entry)
            {
                const Int ip = captured.nodes[entry];

                if (!isModelWallNode(s, ip))
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - invalid captured interface node");
                }

                if (s.node_velocity_bc_type(ip) ==
                    CBSStateSI::velocity_bc_moving_wall)
                {
                    continue;
                }

                const std::array<Real, 3> tangent =
                    project_tangent(
                        captured.values[entry],
                        cache.normal_projector[static_cast<Size>(ip)]);

                s.unkno(1, ip) = tangent[0];
                s.unkno(2, ip) = tangent[1];
                s.unkno(3, ip) = tangent[2];
            }
        }

        void CHTWallModelCoupling::enforceImpermeability(
            CBSStateSI& s,
            const bool owned_only)
        {
            if (!enabled(s))
            {
                return;
            }

            for (const Int ip : cache.wall_nodes)
            {
                if (owned_only && !owned_node(s, ip))
                {
                    continue;
                }

                if (s.node_velocity_bc_type(ip) ==
                    CBSStateSI::velocity_bc_moving_wall)
                {
                    continue;
                }

                const std::array<Real, 3> current =
                {
                    s.unkno(1, ip),
                    s.unkno(2, ip),
                    s.unkno(3, ip)
                };

                const std::array<Real, 3> tangent =
                    project_tangent(
                        current,
                        cache.normal_projector[static_cast<Size>(ip)]);

                s.unkno(1, ip) = tangent[0];
                s.unkno(2, ip) = tangent[1];
                s.unkno(3, ip) = tangent[2];
            }
        }

        CHTWallMomentumDiagnostics CHTWallModelCoupling::addMomentumWallFlux(
            CBSStateSI& s)
        {
            CHTWallMomentumDiagnostics diagnostics;

            if (!enabled(s))
            {
                return diagnostics;
            }

            diagnostics.minimum_y_plus =
                std::numeric_limits<Real>::max();

            for (const ModelFace& face : cache.faces)
            {
                // The complete CHT mesh treats this as an internal face, so the
                // base momentum assembly does not add its characteristic
                // boundary term.  It is nevertheless a boundary of the FLUID
                // momentum subdomain and therefore receives that CBS term here.
                add_characteristic_boundary_term(s, face);

                Real density = 0.0;
                Real molecular_nu = 0.0;
                molecular_properties(
                    s,
                    face.parent,
                    density,
                    molecular_nu);

                const WallTreatmentResult wall =
                    evaluate_momentum_wall(s, face, false);

                const Real face_weight = face.area / 3.0;

                for (const Int ip : face.nodes)
                {
                    for (Int component = 0; component < 3; ++component)
                    {
                        const Real load =
                            face_weight *
                            wall.kinematic_wall_traction[
                                static_cast<Size>(component)];

                        s.rhs(component + 1, ip) += load;

                        diagnostics.assembled_nodal_load[
                            static_cast<Size>(component)] += load;
                    }
                }

                for (Int component = 0; component < 3; ++component)
                {
                    diagnostics.modeled_surface_load[
                        static_cast<Size>(component)] +=
                        face.area *
                        wall.kinematic_wall_traction[
                            static_cast<Size>(component)];
                }

                const Real work =
                    face.area *
                    (wall.kinematic_wall_traction[0] *
                         wall.tangential_velocity[0] +
                     wall.kinematic_wall_traction[1] *
                         wall.tangential_velocity[1] +
                     wall.kinematic_wall_traction[2] *
                         wall.tangential_velocity[2]);

                const Real work_scale =
                    face.area *
                    std::max(
                        1.0,
                        wall.wall_shear_magnitude / density *
                            wall.tangential_speed);

                if (work > 1.0e-13 * work_scale)
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - wall traction performs positive work on the fluid");
                }

                diagnostics.modeled_wall_work += work;
                diagnostics.local_area += face.area;
                ++diagnostics.local_faces;

                if (wall.y_plus > 0.0)
                {
                    diagnostics.minimum_y_plus =
                        std::min(
                            diagnostics.minimum_y_plus,
                            wall.y_plus);

                    diagnostics.maximum_y_plus =
                        std::max(
                            diagnostics.maximum_y_plus,
                            wall.y_plus);
                }
            }

            if (diagnostics.minimum_y_plus ==
                std::numeric_limits<Real>::max())
            {
                diagnostics.minimum_y_plus = 0.0;
            }

            for (Int component = 0; component < 3; ++component)
            {
                const Real surface =
                    diagnostics.modeled_surface_load[
                        static_cast<Size>(component)];

                const Real nodal =
                    diagnostics.assembled_nodal_load[
                        static_cast<Size>(component)];

                const Real tolerance =
                    1.0e-12 *
                    std::max(
                        {1.0, std::abs(surface), std::abs(nodal)});

                if (std::abs(surface - nodal) > tolerance)
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - TRI3 Spalding load conservation failed");
                }
            }

            return diagnostics;
        }

        void CHTWallModelCoupling::prepareThermalStabilityConductivity(
            CBSStateSI& s)
        {
            if (!enabled(s))
            {
                return;
            }

            std::vector<Real> maximum_kn(
                static_cast<Size>(s.cfg.nelem) + 1U,
                0.0);

            // Reconstruct the ordinary bulk turbulent conductivity from the
            // molecular value and the CURRENT SA eddy viscosity.  Do not use
            // the previous value of k_eff_e here because it may still contain
            // the wall-model stability inflation from the preceding iteration.
            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                if (s.mat_elem(ie) == 0)
                {
                    cache.bulk_conductivity[static_cast<Size>(ie)] =
                        bulk_thermal_conductivity(s, ie);

                    s.k_eff_e(ie) =
                        cache.bulk_conductivity[static_cast<Size>(ie)];
                }
            }

            for (const ModelFace& face : cache.faces)
            {
                // Current corrected velocity is used only to establish a
                // conservative NEXT-step stability coefficient.  Step 4 itself
                // uses the lagged beginning-of-iteration velocity, matching the
                // explicit CBS time level used by the wall shear.
                const ThermalWallTreatmentResult thermal =
                    evaluate_thermal_wall(s, face, true);

                maximum_kn[static_cast<Size>(face.parent)] =
                    std::max(
                        maximum_kn[static_cast<Size>(face.parent)],
                        thermal.wall_normal_conductivity);
            }

            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                if (s.mat_elem(ie) != 0)
                {
                    continue;
                }

                const Real bulk =
                    cache.bulk_conductivity[static_cast<Size>(ie)];

                const Real wall =
                    maximum_kn[static_cast<Size>(ie)];

                if (wall > 0.0)
                {
                    s.k_eff_e(ie) = std::max(bulk, wall);
                }
            }
        }

        CHTWallThermalDiagnostics
        CHTWallModelCoupling::correctThermalWallDiffusion(CBSStateSI& s)
        {
            CHTWallThermalDiagnostics diagnostics;

            if (!enabled(s))
            {
                return diagnostics;
            }

            diagnostics.minimum_y_plus =
                std::numeric_limits<Real>::max();
            diagnostics.minimum_prandtl =
                std::numeric_limits<Real>::max();
            diagnostics.minimum_kn_over_k =
                std::numeric_limits<Real>::max();

            std::set<Int> corrected_parents;

            // EnergyAssembly has just used the deliberately conservative
            // isotropic k_eff_e.  First remove that temporary isotropic
            // inflation and recover the ordinary SA bulk operator exactly once
            // per wall-adjacent fluid tetrahedron.
            for (const ModelFace& face : cache.faces)
            {
                const Int ie = face.parent;

                if (!corrected_parents.insert(ie).second)
                {
                    continue;
                }

                const Real k_bulk =
                    cache.bulk_conductivity[static_cast<Size>(ie)];

                const Real k_assembled = s.k_eff_e(ie);

                if (!(k_bulk > 0.0) || !std::isfinite(k_bulk) ||
                    !(k_assembled > 0.0) ||
                    !std::isfinite(k_assembled))
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - invalid thermal conductivity during correction");
                }

                Real grad_t[3] = {0.0, 0.0, 0.0};

                for (Int a = 1; a <= s.cfg.nep; ++a)
                {
                    const Int ip = s.intma(a, ie);

                    for (Int dim = 1; dim <= 3; ++dim)
                    {
                        grad_t[dim - 1] +=
                            grad(s, ie, dim, a) *
                            s.temperature1(ip);
                    }
                }

                const Real volume = s.detJ(ie) / 6.0;
                const Real delta_iso = k_bulk - k_assembled;

                if (delta_iso != 0.0)
                {
                    for (Int a = 1; a <= s.cfg.nep; ++a)
                    {
                        const Int ip = s.intma(a, ie);
                        Real grad_dot = 0.0;

                        for (Int dim = 1; dim <= 3; ++dim)
                        {
                            grad_dot +=
                                grad(s, ie, dim, a) *
                                grad_t[dim - 1];
                        }

                        // Base residual: -k_assembled V gradNa.gradT.
                        // Add -(k_bulk-k_assembled)V gradNa.gradT so the net
                        // coefficient becomes exactly k_bulk.
                        s.rhs1(ip) -=
                            delta_iso * volume * grad_dot;
                    }
                }
            }

            // Add the intended wall-normal rank-one conductivity correction:
            //
            //   K = k_bulk I + (k_n-k_bulk) n tensor n.
            //
            // The sum of its four P1 nodal residuals is zero because
            // sum_a grad(N_a)=0, so no artificial interface heat source is
            // created.  The solid operator is untouched.
            for (const ModelFace& face : cache.faces)
            {
                const Int ie = face.parent;

                const Real k_bulk =
                    cache.bulk_conductivity[static_cast<Size>(ie)];

                // Lag the wall-law velocity to the beginning of the CBS
                // iteration.  The timestep entering this explicit Step 4 was
                // established from the preceding completed state, so using
                // unkn1 here keeps the coefficient at the same time level.
                const ThermalWallTreatmentResult thermal =
                    evaluate_thermal_wall(s, face, false);

                const WallTreatmentResult momentum =
                    evaluate_momentum_wall(s, face, false);

                Real grad_t[3] = {0.0, 0.0, 0.0};

                for (Int a = 1; a <= s.cfg.nep; ++a)
                {
                    const Int ip = s.intma(a, ie);

                    for (Int dim = 1; dim <= 3; ++dim)
                    {
                        grad_t[dim - 1] +=
                            grad(s, ie, dim, a) *
                            s.temperature1(ip);
                    }
                }

                const Real normal_grad_t =
                    grad_t[0] * face.normal[0] +
                    grad_t[1] * face.normal[1] +
                    grad_t[2] * face.normal[2];

                const Real delta_k =
                    thermal.wall_normal_conductivity - k_bulk;

                const Real volume = s.detJ(ie) / 6.0;
                Real face_correction_sum = 0.0;

                for (Int a = 1; a <= s.cfg.nep; ++a)
                {
                    const Int ip = s.intma(a, ie);

                    const Real grad_na_n =
                        grad(s, ie, 1, a) * face.normal[0] +
                        grad(s, ie, 2, a) * face.normal[1] +
                        grad(s, ie, 3, a) * face.normal[2];

                    const Real correction =
                        -delta_k *
                        volume *
                        grad_na_n *
                        normal_grad_t;

                    s.rhs1(ip) += correction;
                    face_correction_sum += correction;
                }

                const Real correction_scale =
                    std::max(
                        1.0,
                        std::abs(delta_k * volume * normal_grad_t) /
                            std::max(face.sample_height, 1.0e-30));

                if (std::abs(face_correction_sum) >
                    5.0e-12 * correction_scale)
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - thermal wall correction violates element energy conservation");
                }

                diagnostics.residual_correction_sum +=
                    face_correction_sum;
                diagnostics.local_area += face.area;
                ++diagnostics.local_faces;

                if (momentum.y_plus > 0.0)
                {
                    diagnostics.minimum_y_plus =
                        std::min(
                            diagnostics.minimum_y_plus,
                            momentum.y_plus);

                    diagnostics.maximum_y_plus =
                        std::max(
                            diagnostics.maximum_y_plus,
                            momentum.y_plus);
                }

                diagnostics.minimum_prandtl =
                    std::min(
                        diagnostics.minimum_prandtl,
                        thermal.molecular_prandtl);

                diagnostics.maximum_prandtl =
                    std::max(
                        diagnostics.maximum_prandtl,
                        thermal.molecular_prandtl);

                const Real ratio =
                    thermal.wall_normal_conductivity /
                    s.k_e(ie);

                diagnostics.minimum_kn_over_k =
                    std::min(
                        diagnostics.minimum_kn_over_k,
                        ratio);

                diagnostics.maximum_kn_over_k =
                    std::max(
                        diagnostics.maximum_kn_over_k,
                        ratio);
            }

            if (diagnostics.minimum_y_plus ==
                std::numeric_limits<Real>::max())
            {
                diagnostics.minimum_y_plus = 0.0;
            }

            if (diagnostics.minimum_prandtl ==
                std::numeric_limits<Real>::max())
            {
                diagnostics.minimum_prandtl = 0.0;
            }

            if (diagnostics.minimum_kn_over_k ==
                std::numeric_limits<Real>::max())
            {
                diagnostics.minimum_kn_over_k = 0.0;
            }

            return diagnostics;
        }
    }
}
