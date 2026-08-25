#pragma once

//=============================================================================
// CBS3D++_SI
//
// Model-consistent coarse-wall coupling for conformal fluid/solid CHT.
//
// Activated only by:
//
//     CBS3D_CHT_WALL_TREATMENT=1
//
// The module is intentionally separate from WallModelCoupling, which remains
// the verified pure-fluid/explicit-wall Spalding path.  CHT interfaces are
// internal tetrahedral faces and therefore need different topology and weak
// forms.
//
// Momentum
// --------
// The fluid side of every fluid/solid face receives the Spalding wall traction.
// Internal CHT faces have fedge=0, so MomentumAssembly has not assembled a
// natural viscous wall term there; consequently no viscous boundary term is
// subtracted.  The CBS characteristic face term that an equivalent explicit
// wall would receive is restored explicitly, then the Spalding traction is
// added.  Strong material no-slip is intercepted afterwards: tangential wall
// velocity is restored while the span of incident wall normals is projected out
// exactly, preserving impermeability at planes, edges and corners.
//
// Energy
// ------
// The conformal CHT discretisation has one shared temperature unknown on the
// interface.  Adding equal/opposite Robin loads at that same DOF would cancel,
// while adding only a fluid load would create energy.  Instead the unresolved
// thermal wall layer is represented by replacing only the FIRST FLUID TET's
// wall-normal conductivity with a Kader-law conductivity:
//
//     K = k_bulk I + sum_f (k_n,f - k_bulk) n_f (x) n_f
//
// where k_bulk is the ordinary SA turbulent conductivity and k_n,f follows from
// the Kader temperature law using the same u_tau obtained from Spalding.  Solid
// conduction and tangential fluid diffusion are untouched.  The correction is
// conservative because sum_a grad(N_a)=0.
//
// Explicit thermal stability is protected without changing TimeStep.cpp: before
// timestep evaluation, wall-adjacent k_eff_e is inflated to max(k_bulk,k_n).
// EnergyAssembly temporarily sees that conservative isotropic upper bound; the
// residual correction below then removes the artificial tangential inflation
// and installs the intended anisotropic wall-normal tensor.  After every SA
// update the same preparation is repeated, so the NEXT timestep always sees a
// conductivity at least as restrictive as the wall model used in Step 4.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"
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
        struct CHTCapturedWallVelocity
        {
            std::vector<Int> nodes;
            std::vector<std::array<Real, 3>> values;
        };

        struct CHTWallMomentumDiagnostics
        {
            long long local_faces = 0;
            Real local_area = 0.0;
            Real minimum_y_plus = 0.0;
            Real maximum_y_plus = 0.0;
            Real modeled_wall_work = 0.0;
            std::array<Real, 3> modeled_surface_load = {0.0, 0.0, 0.0};
            std::array<Real, 3> assembled_nodal_load = {0.0, 0.0, 0.0};
        };

        struct CHTWallThermalDiagnostics
        {
            long long local_faces = 0;
            Real local_area = 0.0;
            Real minimum_y_plus = 0.0;
            Real maximum_y_plus = 0.0;
            Real minimum_prandtl = 0.0;
            Real maximum_prandtl = 0.0;
            Real minimum_kn_over_k = 0.0;
            Real maximum_kn_over_k = 0.0;
            Real residual_correction_sum = 0.0;
        };

        class CHTWallModelCoupling
        {
        private:
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

            struct Cache
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

            inline static Cache cache_;

            static bool environmentFlag(const char* name)
            {
                const char* value = std::getenv(name);
                return value != nullptr && value[0] != '\0' &&
                    std::string(value) != "0";
            }

            static Real environmentReal(
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
                        std::string("CHTWallModelCoupling - invalid real environment value for ") + name);
                }

                return parsed;
            }

            static bool requested()
            {
                return environmentFlag("CBS3D_CHT_WALL_TREATMENT");
            }

            static void checkMpi(const int code, const char* operation)
            {
#ifdef CBS3D_USE_MPI
                if (code != MPI_SUCCESS)
                {
                    throw std::runtime_error(
                        std::string("CHTWallModelCoupling MPI failure in ") + operation);
                }
#else
                (void)code;
                (void)operation;
#endif
            }

            static bool isMaterialInterfaceNode(
                const CBSStateSI& s,
                const Int ip)
            {
                if (ip < 1 || ip > s.cfg.npoin)
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - node outside local mesh");
                }

                const Int mask =
                    CBSStateSI::node_touches_fluid |
                    CBSStateSI::node_touches_solid;

                return s.node_material_mask(ip) == mask;
            }

            static long long globalNodeId(
                const CBSStateSI& s,
                const Int ip)
            {
#ifdef CBS3D_USE_MPI
                if (s.mpi_enabled && s.mpi_size > 1)
                {
                    if (static_cast<Size>(ip) >= s.local_to_global_node.size())
                    {
                        throw std::runtime_error(
                            "CHTWallModelCoupling - missing local-to-global node map");
                    }

                    const Int gid = s.local_to_global_node[static_cast<Size>(ip)];

                    if (gid < 1)
                    {
                        throw std::runtime_error(
                            "CHTWallModelCoupling - invalid global node id");
                    }

                    return static_cast<long long>(gid);
                }
#else
                (void)s;
#endif

                return static_cast<long long>(ip);
            }

            static Int dNkdxIndex(
                const CBSStateSI& s,
                const Int ie,
                const Int dim,
                const Int local_node)
            {
                return (ie - 1) * s.cfg.ndim * s.cfg.nep
                    + (dim - 1) * s.cfg.nep
                    + local_node;
            }

            static Real grad(
                const CBSStateSI& s,
                const Int ie,
                const Int dim,
                const Int local_node)
            {
                return s.dNkdx(dNkdxIndex(s, ie, dim, local_node));
            }

            static void validateConfiguration(const CBSStateSI& s)
            {
                if (!requested())
                {
                    return;
                }

                if (environmentFlag("CBS3D_SA_WALL_TREATMENT"))
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
                    checkMpi(
                        MPI_Allreduce(&local_fluid, &global_fluid, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD),
                        "MPI_Allreduce fluid-domain audit");
                    checkMpi(
                        MPI_Allreduce(&local_solid, &global_solid, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD),
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

            static CandidateFace candidateFace(
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
                    const Int ip = s.intma(s.ippn1(face, corner), ie);
                    result.local_key[static_cast<Size>(corner - 1)] = ip;
                    result.global_key[static_cast<Size>(corner - 1)] =
                        globalNodeId(s, ip);
                }

                std::sort(result.local_key.begin(), result.local_key.end());
                std::sort(result.global_key.begin(), result.global_key.end());
                return result;
            }

            static bool canBeInterface(
                const CBSStateSI& s,
                const CandidateFace& face)
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

            static ModelFace makeModelFace(
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
                    const Real component = s.annxf(dim, face, ie) / area;
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
                    const Int ip = s.intma(s.ippn1(face, corner), ie);
                    result.nodes[static_cast<Size>(corner - 1)] = ip;

                    if (ip == result.sample_node)
                    {
                        throw std::runtime_error(
                            "CHTWallModelCoupling - opposite sample node lies on interface face");
                    }
                }

                return result;
            }

            static void collectLocalFaces(
                const CBSStateSI& s,
                std::vector<ModelFace>& accepted,
                std::vector<CandidateFace>& partition_candidates)
            {
                std::vector<CandidateFace> candidates;

                for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
                {
                    for (Int face = 1; face <= s.cfg.nsid; ++face)
                    {
                        CandidateFace candidate = candidateFace(s, ie, face);

                        if (canBeInterface(s, candidate))
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
                           candidates[end].local_key == candidates[begin].local_key)
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
                            accepted.push_back(makeModelFace(s, a.fluid ? a : b));
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

            static void collectCrossRankFaces(
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

                const int local_count = static_cast<int>(local_candidates.size());
                std::vector<int> counts(static_cast<Size>(s.mpi_size), 0);

                checkMpi(
                    MPI_Allgather(
                        &local_count, 1, MPI_INT,
                        counts.data(), 1, MPI_INT,
                        MPI_COMM_WORLD),
                    "MPI_Allgather interface candidate counts");

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
                        "CHTWallModelCoupling - global candidate count exceeds MPI integer range");
                }

                const int total = static_cast<int>(total_ll);
                std::vector<long long> send_keys(static_cast<Size>(local_count) * 3U);
                std::vector<int> send_material(static_cast<Size>(local_count));

                for (int i = 0; i < local_count; ++i)
                {
                    const CandidateFace& face = local_candidates[static_cast<Size>(i)];
                    send_material[static_cast<Size>(i)] = face.fluid ? 1 : 0;

                    for (Size k = 0; k < 3; ++k)
                    {
                        send_keys[static_cast<Size>(i) * 3U + k] = face.global_key[k];
                    }
                }

                std::vector<int> counts3(static_cast<Size>(s.mpi_size));
                std::vector<int> offsets3(static_cast<Size>(s.mpi_size));

                for (Int rank = 0; rank < s.mpi_size; ++rank)
                {
                    counts3[static_cast<Size>(rank)] =
                        3 * counts[static_cast<Size>(rank)];
                    offsets3[static_cast<Size>(rank)] =
                        3 * offsets[static_cast<Size>(rank)];
                }

                std::vector<long long> all_keys(static_cast<Size>(total) * 3U);
                std::vector<int> all_material(static_cast<Size>(total));

                checkMpi(
                    MPI_Allgatherv(
                        send_keys.data(), local_count * 3, MPI_LONG_LONG,
                        all_keys.data(), counts3.data(), offsets3.data(), MPI_LONG_LONG,
                        MPI_COMM_WORLD),
                    "MPI_Allgatherv interface candidate keys");

                checkMpi(
                    MPI_Allgatherv(
                        send_material.data(), local_count, MPI_INT,
                        all_material.data(), counts.data(), offsets.data(), MPI_INT,
                        MPI_COMM_WORLD),
                    "MPI_Allgatherv interface material classes");

                struct Gathered
                {
                    std::array<long long, 3> key = {0, 0, 0};
                    bool fluid = false;
                };

                std::vector<Gathered> gathered(static_cast<Size>(total));

                for (int i = 0; i < total; ++i)
                {
                    Gathered& item = gathered[static_cast<Size>(i)];
                    item.fluid = all_material[static_cast<Size>(i)] != 0;

                    for (Size k = 0; k < 3; ++k)
                    {
                        item.key[k] = all_keys[static_cast<Size>(i) * 3U + k];
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
                        gathered[begin].fluid != gathered[begin + 1U].fluid)
                    {
                        interface_keys.insert(gathered[begin].key);
                    }

                    begin = end;
                }

                for (const CandidateFace& candidate : local_candidates)
                {
                    if (candidate.fluid &&
                        interface_keys.find(candidate.global_key) != interface_keys.end())
                    {
                        accepted.push_back(makeModelFace(s, candidate));
                    }
                }
#else
                (void)s;
                (void)local_candidates;
                (void)accepted;
#endif
            }

            static std::array<Real, 9> projectorFromGram(
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

            static bool nonzeroProjector(const std::array<Real, 9>& projector)
            {
                return projector[0] + projector[4] + projector[8] > 0.5;
            }

            static void buildCache(const CBSStateSI& s)
            {
                validateConfiguration(s);

                cache_ = Cache{};
                cache_.state = &s;
                cache_.npoin = s.cfg.npoin;
                cache_.nelem = s.cfg.nelem;

                if (!requested())
                {
                    cache_.ready = true;
                    return;
                }

                std::vector<CandidateFace> partition_candidates;
                collectLocalFaces(s, cache_.faces, partition_candidates);
                collectCrossRankFaces(s, partition_candidates, cache_.faces);

                std::sort(
                    cache_.faces.begin(),
                    cache_.faces.end(),
                    [](const ModelFace& lhs, const ModelFace& rhs)
                    {
                        return lhs.global_key < rhs.global_key;
                    });

                if (std::adjacent_find(
                        cache_.faces.begin(),
                        cache_.faces.end(),
                        [](const ModelFace& lhs, const ModelFace& rhs)
                        {
                            return lhs.global_key == rhs.global_key;
                        }) != cache_.faces.end())
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - duplicate local CHT interface face");
                }

                long long local_faces =
                    static_cast<long long>(cache_.faces.size());
                cache_.global_faces = local_faces;

#ifdef CBS3D_USE_MPI
                if (s.mpi_enabled && s.mpi_size > 1)
                {
                    checkMpi(
                        MPI_Allreduce(
                            &local_faces,
                            &cache_.global_faces,
                            1,
                            MPI_LONG_LONG,
                            MPI_SUM,
                            MPI_COMM_WORLD),
                        "MPI_Allreduce CHT wall-face count");
                }
#endif

                if (cache_.global_faces <= 0)
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - no conformal fluid-solid interface faces found globally");
                }

                Array2D<Real> gram;
                gram.resize(9, s.cfg.npoin);
                gram.fill(0.0);

                for (const ModelFace& face : cache_.faces)
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

                cache_.normal_projector.resize(
                    static_cast<Size>(s.cfg.npoin) + 1U);

                for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
                {
                    std::array<Real, 9> local_gram{};

                    for (Int entry = 0; entry < 9; ++entry)
                    {
                        local_gram[static_cast<Size>(entry)] =
                            gram(entry + 1, ip);
                    }

                    cache_.normal_projector[static_cast<Size>(ip)] =
                        projectorFromGram(local_gram);

                    if (nonzeroProjector(
                            cache_.normal_projector[static_cast<Size>(ip)]))
                    {
                        cache_.wall_nodes.push_back(ip);
                    }
                }

                long long wall_node_occurrences =
                    static_cast<long long>(cache_.wall_nodes.size());
                long long global_wall_node_occurrences = wall_node_occurrences;

#ifdef CBS3D_USE_MPI
                if (s.mpi_enabled && s.mpi_size > 1)
                {
                    checkMpi(
                        MPI_Allreduce(
                            &wall_node_occurrences,
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

                cache_.bulk_conductivity.assign(
                    static_cast<Size>(s.cfg.nelem) + 1U,
                    0.0);
                cache_.ready = true;

                if (!s.mpi_enabled || s.mpi_rank == 0)
                {
                    WallTreatmentOptions options;
                    options.kappa =
                        environmentReal("CBS3D_WALL_KAPPA", options.kappa);
                    options.log_intercept =
                        environmentReal("CBS3D_WALL_B", options.log_intercept);

                    std::cout
                        << "  SA CHT wall treatment: ON\n"
                        << "    momentum model      : Spalding smooth wall\n"
                        << "    thermal model       : Kader continuous T+\n"
                        << "    kappa / B           : "
                        << options.kappa << " / " << options.log_intercept << "\n"
                        << "    global CHT faces    : " << cache_.global_faces << "\n"
                        << "    coupling            : internal fluid-solid topology\n";
                }
            }

            static void ensureCache(const CBSStateSI& s)
            {
                if (!cache_.ready ||
                    cache_.state != &s ||
                    cache_.npoin != s.cfg.npoin ||
                    cache_.nelem != s.cfg.nelem)
                {
                    buildCache(s);
                }
            }

            static bool ownedNode(const CBSStateSI& s, const Int ip)
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

            static std::array<Real, 3> projectTangent(
                const std::array<Real, 3>& velocity,
                const std::array<Real, 9>& projector)
            {
                std::array<Real, 3> normal{};

                for (Int i = 0; i < 3; ++i)
                {
                    for (Int j = 0; j < 3; ++j)
                    {
                        normal[static_cast<Size>(i)] +=
                            projector[static_cast<Size>(3 * i + j)] *
                            velocity[static_cast<Size>(j)];
                    }
                }

                return
                {
                    velocity[0] - normal[0],
                    velocity[1] - normal[1],
                    velocity[2] - normal[2]
                };
            }

            static WallTreatmentOptions wallOptions()
            {
                WallTreatmentOptions options;
                options.kappa = environmentReal("CBS3D_WALL_KAPPA", options.kappa);
                options.log_intercept = environmentReal("CBS3D_WALL_B", options.log_intercept);
                return options;
            }

            static void molecularProperties(
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

            static Real bulkThermalConductivity(
                const CBSStateSI& s,
                const Int ie)
            {
                if (s.mat_elem(ie) != 0)
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - bulk fluid conductivity requested on solid");
                }

                Real k = s.k_e(ie);

                if (s.cfg.turbulence_on > 0 &&
                    s.cfg.turbulent_thermal_diffusivity_on > 0)
                {
                    if (s.nu_t_e(ie) < 0.0 || !std::isfinite(s.nu_t_e(ie)))
                    {
                        throw std::runtime_error(
                            "CHTWallModelCoupling - invalid element eddy viscosity");
                    }

                    k += s.rho_cp_e(ie) *
                        s.nu_t_e(ie) / s.cfg.sa_prandtl_t;
                }

                if (!(k > 0.0) || !std::isfinite(k))
                {
                    throw std::runtime_error(
                        "CHTWallModelCoupling - invalid bulk fluid conductivity");
                }

                return k;
            }

            static WallTreatmentResult evaluateMomentumWall(
                const CBSStateSI& s,
                const ModelFace& face,
                const bool current_velocity)
            {
                Real density = 0.0;
                Real molecular_nu = 0.0;
                molecularProperties(s, face.parent, density, molecular_nu);

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
                    wallOptions());
            }

            static ThermalWallTreatmentResult evaluateThermalWall(
                const CBSStateSI& s,
                const ModelFace& face,
                const bool current_velocity)
            {
                const WallTreatmentResult momentum =
                    evaluateMomentumWall(s, face, current_velocity);

                const Int ie = face.parent;
                const Real density = s.rho_e(ie);
                const Real cp = s.rho_cp_e(ie) / density;

                return ThermalWallTreatment::evaluateKader(
                    momentum.friction_velocity,
                    face.sample_height,
                    density,
                    cp,
                    s.mu_e(ie),
                    s.k_e(ie));
            }

            static void addCharacteristicBoundaryTerm(
                CBSStateSI& s,
                const ModelFace& face,
                CHTWallMomentumDiagnostics& diagnostics)
            {
                const Int ie = face.parent;
                Real flux_gradient[4][4] = {};
                Real mean_velocity[4] = {};

                // d(u_j u_i)/dx_j, using the same nodal old velocity and P1
                // gradients as MomentumAssembly::step1_characteristic_correction.
                for (Int i = 1; i <= 3; ++i)
                {
                    for (Int j = 1; j <= 3; ++j)
                    {
                        for (Int a = 1; a <= s.cfg.nep; ++a)
                        {
                            const Int ip = s.intma(a, ie);
                            flux_gradient[i][j] +=
                                grad(s, ie, j, a) *
                                s.unkn1(j, ip) * s.unkn1(i, ip);
                        }
                    }
                }

                for (Int k = 1; k <= 3; ++k)
                {
                    for (Int a = 1; a <= s.cfg.nep; ++a)
                    {
                        mean_velocity[k] += s.unkn1(k, s.intma(a, ie));
                    }
                    mean_velocity[k] /= static_cast<Real>(s.cfg.nep);
                }

                const Real coefficient =
                    0.5 * s.delte(ie) * s.cfg.fcon[4];

                if (!(coefficient >= 0.0) || !std::isfinite(coefficient))
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
                            mean_velocity[k] * divergence *
                            s.annxf(k, face.local_face, ie) * coefficient;
                    }

                    for (const Int ip : face.nodes)
                    {
                        s.rhs(i, ip) += face_value;
                    }
                }

                (void)diagnostics;
            }

        public:
            static bool enabled(const CBSStateSI& s)
            {
                if (!requested())
                {
                    return false;
                }

                ensureCache(s);
                return true;
            }

            static long long globalWallFaceCount(const CBSStateSI& s)
            {
                return enabled(s) ? cache_.global_faces : 0;
            }

            static bool isModelWallNode(const CBSStateSI& s, const Int ip)
            {
                if (!enabled(s) || ip < 1 || ip > s.cfg.npoin)
                {
                    return false;
                }

                return nonzeroProjector(
                    cache_.normal_projector[static_cast<Size>(ip)]);
            }

            static CHTCapturedWallVelocity captureVelocity(
                const CBSStateSI& s,
                const bool owned_only)
            {
                CHTCapturedWallVelocity captured;

                if (!enabled(s))
                {
                    return captured;
                }

                captured.nodes.reserve(cache_.wall_nodes.size());
                captured.values.reserve(cache_.wall_nodes.size());

                for (const Int ip : cache_.wall_nodes)
                {
                    if (owned_only && !ownedNode(s, ip))
                    {
                        continue;
                    }

                    captured.nodes.push_back(ip);
                    captured.values.push_back(
                        {s.unkno(1, ip), s.unkno(2, ip), s.unkno(3, ip)});
                }

                return captured;
            }

            static void restoreTangentialAndEnforceImpermeability(
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
                        projectTangent(
                            captured.values[entry],
                            cache_.normal_projector[static_cast<Size>(ip)]);

                    s.unkno(1, ip) = tangent[0];
                    s.unkno(2, ip) = tangent[1];
                    s.unkno(3, ip) = tangent[2];
                }
            }

            static void enforceImpermeability(
                CBSStateSI& s,
                const bool owned_only)
            {
                if (!enabled(s))
                {
                    return;
                }

                for (const Int ip : cache_.wall_nodes)
                {
                    if (owned_only && !ownedNode(s, ip))
                    {
                        continue;
                    }

                    const std::array<Real, 3> current =
                    {s.unkno(1, ip), s.unkno(2, ip), s.unkno(3, ip)};
                    const std::array<Real, 3> tangent =
                        projectTangent(
                            current,
                            cache_.normal_projector[static_cast<Size>(ip)]);

                    s.unkno(1, ip) = tangent[0];
                    s.unkno(2, ip) = tangent[1];
                    s.unkno(3, ip) = tangent[2];
                }
            }

            static CHTWallMomentumDiagnostics addMomentumWallFlux(CBSStateSI& s)
            {
                CHTWallMomentumDiagnostics diagnostics;

                if (!enabled(s))
                {
                    return diagnostics;
                }

                diagnostics.minimum_y_plus =
                    std::numeric_limits<Real>::max();

                for (const ModelFace& face : cache_.faces)
                {
                    // Restore the characteristic boundary term that the base
                    // CBS assembly omits because an internal CHT face has
                    // fedge=0.  Convective normal flux needs no analogous term:
                    // impermeability makes u.n identically zero on the face.
                    addCharacteristicBoundaryTerm(s, face, diagnostics);

                    Real density = 0.0;
                    Real molecular_nu = 0.0;
                    molecularProperties(s, face.parent, density, molecular_nu);

                    const WallTreatmentResult wall =
                        evaluateMomentumWall(s, face, false);
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

                    const Real work = face.area *
                        (wall.kinematic_wall_traction[0] * wall.tangential_velocity[0] +
                         wall.kinematic_wall_traction[1] * wall.tangential_velocity[1] +
                         wall.kinematic_wall_traction[2] * wall.tangential_velocity[2]);

                    if (work > 1.0e-13 * face.area *
                        std::max(1.0, wall.wall_shear_magnitude / density *
                                      wall.tangential_speed))
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
                            std::min(diagnostics.minimum_y_plus, wall.y_plus);
                        diagnostics.maximum_y_plus =
                            std::max(diagnostics.maximum_y_plus, wall.y_plus);
                    }
                }

                if (diagnostics.minimum_y_plus ==
                    std::numeric_limits<Real>::max())
                {
                    diagnostics.minimum_y_plus = 0.0;
                }

                // The characteristic correction is conservative in the CBS
                // formulation but is deliberately excluded from this check:
                // this audit concerns only the newly imposed wall shear load.
                for (Int component = 0; component < 3; ++component)
                {
                    const Real surface =
                        diagnostics.modeled_surface_load[static_cast<Size>(component)];
                    const Real nodal =
                        diagnostics.assembled_nodal_load[static_cast<Size>(component)];
                    const Real tolerance =
                        1.0e-12 * std::max({1.0, std::abs(surface), std::abs(nodal)});

                    if (std::abs(surface - nodal) > tolerance)
                    {
                        throw std::runtime_error(
                            "CHTWallModelCoupling - TRI3 Spalding load conservation failed");
                    }
                }

                return diagnostics;
            }

            // Rebuild the ordinary SA turbulent conductivity exactly from the
            // current nu_t_e, then inflate only wall-adjacent fluid elements to
            // a conservative isotropic stability bound.  TimeStep therefore
            // sees alpha >= every Kader wall-normal alpha used by Step 4.
            static void prepareThermalStabilityConductivity(CBSStateSI& s)
            {
                if (!enabled(s))
                {
                    return;
                }

                std::vector<Real> maximum_kn(
                    static_cast<Size>(s.cfg.nelem) + 1U,
                    0.0);

                for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
                {
                    if (s.mat_elem(ie) == 0)
                    {
                        cache_.bulk_conductivity[static_cast<Size>(ie)] =
                            bulkThermalConductivity(s, ie);
                        s.k_eff_e(ie) =
                            cache_.bulk_conductivity[static_cast<Size>(ie)];
                    }
                }

                for (const ModelFace& face : cache_.faces)
                {
                    const ThermalWallTreatmentResult thermal =
                        evaluateThermalWall(s, face, true);
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
                        cache_.bulk_conductivity[static_cast<Size>(ie)];
                    const Real wall = maximum_kn[static_cast<Size>(ie)];

                    if (wall > 0.0)
                    {
                        s.k_eff_e(ie) = std::max(bulk, wall);
                    }
                }
            }

            static CHTWallThermalDiagnostics correctThermalWallDiffusion(
                CBSStateSI& s)
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

                // First remove the deliberately inflated ISOTROPIC stability
                // conductivity.  EnergyAssembly used k_eff_e; the intended bulk
                // operator uses the ordinary SA k_bulk.
                for (const ModelFace& face : cache_.faces)
                {
                    const Int ie = face.parent;

                    if (!corrected_parents.insert(ie).second)
                    {
                        continue;
                    }

                    const Real k_bulk =
                        cache_.bulk_conductivity[static_cast<Size>(ie)];
                    const Real k_assembled = s.k_eff_e(ie);

                    if (!(k_bulk > 0.0) || !std::isfinite(k_bulk) ||
                        !(k_assembled > 0.0) || !std::isfinite(k_assembled))
                    {
                        throw std::runtime_error(
                            "CHTWallModelCoupling - invalid thermal conductivity during correction");
                    }

                    Real gradT[3] = {0.0, 0.0, 0.0};

                    for (Int a = 1; a <= s.cfg.nep; ++a)
                    {
                        const Int ip = s.intma(a, ie);
                        for (Int dim = 1; dim <= 3; ++dim)
                        {
                            gradT[dim - 1] +=
                                grad(s, ie, dim, a) * s.temperature1(ip);
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
                                    grad(s, ie, dim, a) * gradT[dim - 1];
                            }
                            s.rhs1(ip) -= delta_iso * volume * grad_dot;
                        }
                    }
                }

                // Then install one Kader wall-normal rank-one correction for
                // every incident CHT face.
                for (const ModelFace& face : cache_.faces)
                {
                    const Int ie = face.parent;
                    const Real k_bulk =
                        cache_.bulk_conductivity[static_cast<Size>(ie)];
                    const ThermalWallTreatmentResult thermal =
                        evaluateThermalWall(s, face, false);
                    const WallTreatmentResult momentum =
                        evaluateMomentumWall(s, face, false);

                    Real gradT[3] = {0.0, 0.0, 0.0};

                    for (Int a = 1; a <= s.cfg.nep; ++a)
                    {
                        const Int ip = s.intma(a, ie);
                        for (Int dim = 1; dim <= 3; ++dim)
                        {
                            gradT[dim - 1] +=
                                grad(s, ie, dim, a) * s.temperature1(ip);
                        }
                    }

                    const Real normal_gradT =
                        gradT[0] * face.normal[0] +
                        gradT[1] * face.normal[1] +
                        gradT[2] * face.normal[2];

                    const Real delta_k =
                        thermal.wall_normal_conductivity - k_bulk;
                    const Real volume = s.detJ(ie) / 6.0;
                    Real face_correction_sum = 0.0;

                    for (Int a = 1; a <= s.cfg.nep; ++a)
                    {
                        const Int ip = s.intma(a, ie);
                        const Real gradNa_n =
                            grad(s, ie, 1, a) * face.normal[0] +
                            grad(s, ie, 2, a) * face.normal[1] +
                            grad(s, ie, 3, a) * face.normal[2];

                        const Real correction =
                            -delta_k * volume * gradNa_n * normal_gradT;

                        s.rhs1(ip) += correction;
                        face_correction_sum += correction;
                    }

                    const Real correction_scale =
                        std::max(
                            1.0,
                            std::abs(delta_k * volume * normal_gradT) /
                                std::max(face.sample_height, 1.0e-30));

                    if (std::abs(face_correction_sum) >
                        5.0e-12 * correction_scale)
                    {
                        throw std::runtime_error(
                            "CHTWallModelCoupling - thermal wall correction violates element energy conservation");
                    }

                    diagnostics.residual_correction_sum += face_correction_sum;
                    diagnostics.local_area += face.area;
                    ++diagnostics.local_faces;

                    if (momentum.y_plus > 0.0)
                    {
                        diagnostics.minimum_y_plus =
                            std::min(diagnostics.minimum_y_plus, momentum.y_plus);
                        diagnostics.maximum_y_plus =
                            std::max(diagnostics.maximum_y_plus, momentum.y_plus);
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
                        thermal.wall_normal_conductivity / s.k_e(ie);
                    diagnostics.minimum_kn_over_k =
                        std::min(diagnostics.minimum_kn_over_k, ratio);
                    diagnostics.maximum_kn_over_k =
                        std::max(diagnostics.maximum_kn_over_k, ratio);
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
        };
    }
}
