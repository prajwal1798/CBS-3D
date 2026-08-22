//=============================================================================
// CBS3D++_SI
//
// Production momentum coupling for the smooth-wall Spalding treatment.
//=============================================================================

#include "cbs/turbulence/WallModelCoupling.hpp"

#include "cbs/parallel/HaloExchange.hpp"
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
                std::array<Real, 3> normal = {0.0, 0.0, 0.0};
                Real area = 0.0;
                Real sample_height = 0.0;
            };

            struct CouplingCache
            {
                const CBSStateSI* state = nullptr;
                Int npoin = 0;
                Int nelem = 0;
                Int nboun = 0;
                bool ready = false;
                long long global_faces = 0;
                std::vector<ModelFace> faces;
                std::vector<Int> wall_nodes;
                std::vector<std::array<Real, 9>> normal_projector;
            };

            CouplingCache cache;

            bool environment_flag(const char* name)
            {
                const char* value = std::getenv(name);
                return value != nullptr && value[0] != '\0' && std::string(value) != "0";
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
                        std::string("WallModelCoupling - invalid real environment value for ")
                        + name);
                }

                return parsed;
            }

            bool requested()
            {
                return environment_flag("CBS3D_SA_WALL_TREATMENT");
            }

            bool modelled_bc(const CBSStateSI& s, const Int bc)
            {
                return
                    bc == s.cfg.bc_noslip_adiabatic_wall ||
                    bc == s.cfg.bc_noslip_heatflux_wall;
            }

            void check_mpi(const int code, const char* operation)
            {
#ifdef CBS3D_USE_MPI
                if (code != MPI_SUCCESS)
                {
                    throw std::runtime_error(
                        std::string("WallModelCoupling MPI failure in ") + operation);
                }
#else
                (void)code;
                (void)operation;
#endif
            }

            bool global_contains_solid(const CBSStateSI& s)
            {
                int local = 0;

                for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
                {
                    if (s.mat_elem(ie) != 0)
                    {
                        local = 1;
                        break;
                    }
                }

#ifdef CBS3D_USE_MPI
                if (s.mpi_enabled && s.mpi_size > 1)
                {
                    int global = 0;
                    check_mpi(
                        MPI_Allreduce(
                            &local,
                            &global,
                            1,
                            MPI_INT,
                            MPI_MAX,
                            MPI_COMM_WORLD),
                        "MPI_Allreduce material-domain audit");
                    return global != 0;
                }
#endif

                return local != 0;
            }

            void validate_configuration(const CBSStateSI& s)
            {
                if (!requested())
                {
                    return;
                }

                if (s.cfg.turbulence_on < 1)
                {
                    throw std::runtime_error(
                        "WallModelCoupling - CBS3D_SA_WALL_TREATMENT requires turbulence_on=1");
                }

                if (s.cfg.turbulence_model != 0)
                {
                    throw std::runtime_error(
                        "WallModelCoupling - production wall treatment currently supports standard SA only");
                }

                if (global_contains_solid(s))
                {
                    throw std::runtime_error(
                        "WallModelCoupling - the current production coupling is intentionally restricted to pure-fluid SA. "
                        "CHT/interface wall modelling is deferred until the momentum model passes flat-plate verification "
                        "and the thermal wall treatment is derived.");
                }

                if (s.cfg.ndim != 3 || s.cfg.nep != 4 ||
                    s.cfg.nsid != 4 || s.cfg.nsidp != 3)
                {
                    throw std::runtime_error(
                        "WallModelCoupling - production coupling requires 3-D P1 TET4/TRI3 topology");
                }
            }

            long long global_node_id(const CBSStateSI& s, const Int ip)
            {
                if (ip < 1 || ip > s.cfg.npoin)
                {
                    throw std::runtime_error(
                        "WallModelCoupling - node outside local mesh");
                }

#ifdef CBS3D_USE_MPI
                if (s.mpi_enabled && s.mpi_size > 1)
                {
                    if (static_cast<Size>(ip) >= s.local_to_global_node.size())
                    {
                        throw std::runtime_error(
                            "WallModelCoupling - missing local-to-global node map");
                    }

                    const Int gid = s.local_to_global_node[static_cast<Size>(ip)];

                    if (gid < 1)
                    {
                        throw std::runtime_error(
                            "WallModelCoupling - invalid global node id");
                    }

                    return static_cast<long long>(gid);
                }
#endif

                return static_cast<long long>(ip);
            }

            long long audit_global_face_uniqueness(
                const CBSStateSI& s,
                const std::vector<ModelFace>& faces)
            {
                std::vector<std::array<long long, 3>> local_keys;
                local_keys.reserve(faces.size());

                for (const ModelFace& face : faces)
                {
                    std::array<long long, 3> key =
                    {
                        global_node_id(s, face.nodes[0]),
                        global_node_id(s, face.nodes[1]),
                        global_node_id(s, face.nodes[2])
                    };
                    std::sort(key.begin(), key.end());
                    local_keys.push_back(key);
                }

                std::sort(local_keys.begin(), local_keys.end());

                if (std::adjacent_find(local_keys.begin(), local_keys.end()) != local_keys.end())
                {
                    throw std::runtime_error(
                        "WallModelCoupling - duplicate model wall face on one rank");
                }

#ifdef CBS3D_USE_MPI
                if (s.mpi_enabled && s.mpi_size > 1)
                {
                    if (local_keys.size() >
                        static_cast<Size>(std::numeric_limits<int>::max()))
                    {
                        throw std::runtime_error(
                            "WallModelCoupling - local model-wall count exceeds MPI integer range");
                    }

                    const int local_count = static_cast<int>(local_keys.size());
                    std::vector<int> counts(static_cast<Size>(s.mpi_size), 0);

                    check_mpi(
                        MPI_Allgather(
                            &local_count,
                            1,
                            MPI_INT,
                            counts.data(),
                            1,
                            MPI_INT,
                            MPI_COMM_WORLD),
                        "MPI_Allgather wall-face counts");

                    std::vector<int> offsets(static_cast<Size>(s.mpi_size), 0);
                    long long total = 0;

                    for (Int rank = 0; rank < s.mpi_size; ++rank)
                    {
                        offsets[static_cast<Size>(rank)] = static_cast<int>(total);
                        total += counts[static_cast<Size>(rank)];
                    }

                    if (total > std::numeric_limits<int>::max())
                    {
                        throw std::runtime_error(
                            "WallModelCoupling - global model-wall count exceeds MPI integer range");
                    }

                    std::vector<long long> send(
                        static_cast<Size>(local_count) * 3U);

                    for (int i = 0; i < local_count; ++i)
                    {
                        for (Size k = 0; k < 3; ++k)
                        {
                            send[static_cast<Size>(i) * 3U + k] =
                                local_keys[static_cast<Size>(i)][k];
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

                    std::vector<long long> receive(
                        static_cast<Size>(total) * 3U);

                    check_mpi(
                        MPI_Allgatherv(
                            send.data(),
                            local_count * 3,
                            MPI_LONG_LONG,
                            receive.data(),
                            counts3.data(),
                            offsets3.data(),
                            MPI_LONG_LONG,
                            MPI_COMM_WORLD),
                        "MPI_Allgatherv wall-face keys");

                    std::vector<std::array<long long, 3>> global_keys(
                        static_cast<Size>(total));

                    for (int i = 0; i < static_cast<int>(total); ++i)
                    {
                        for (Size k = 0; k < 3; ++k)
                        {
                            global_keys[static_cast<Size>(i)][k] =
                                receive[static_cast<Size>(i) * 3U + k];
                        }
                    }

                    std::sort(global_keys.begin(), global_keys.end());

                    if (std::adjacent_find(global_keys.begin(), global_keys.end()) !=
                        global_keys.end())
                    {
                        throw std::runtime_error(
                            "WallModelCoupling - one physical model wall face is represented on more than one MPI rank");
                    }

                    return total;
                }
#endif

                return static_cast<long long>(faces.size());
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
                        gram[static_cast<Size>(0 * 3 + column)],
                        gram[static_cast<Size>(1 * 3 + column)],
                        gram[static_cast<Size>(2 * 3 + column)]
                    };

                    for (Int q = 0; q < rank; ++q)
                    {
                        const Real projection =
                            vector[0] * basis[static_cast<Size>(q)][0] +
                            vector[1] * basis[static_cast<Size>(q)][1] +
                            vector[2] * basis[static_cast<Size>(q)][2];

                        for (Int k = 0; k < 3; ++k)
                        {
                            vector[static_cast<Size>(k)] -=
                                projection * basis[static_cast<Size>(q)][static_cast<Size>(k)];
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
                            projector[static_cast<Size>(i * 3 + j)] +=
                                basis[static_cast<Size>(q)][static_cast<Size>(i)] *
                                basis[static_cast<Size>(q)][static_cast<Size>(j)];
                        }
                    }
                }

                return projector;
            }

            bool nonzero_projector(const std::array<Real, 9>& projector)
            {
                const Real trace = projector[0] + projector[4] + projector[8];
                return trace > 0.5;
            }

            void build_cache(const CBSStateSI& s)
            {
                validate_configuration(s);

                cache = CouplingCache{};
                cache.state = &s;
                cache.npoin = s.cfg.npoin;
                cache.nelem = s.cfg.nelem;
                cache.nboun = s.cfg.nboun;

                if (!requested())
                {
                    cache.ready = true;
                    return;
                }

                Array2D<Real> gram;
                gram.resize(9, s.cfg.npoin);
                gram.fill(0.0);

                for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
                {
                    const Int bc = s.iside(s.cfg.bsid, ib);

                    if (!modelled_bc(s, bc))
                    {
                        continue;
                    }

                    const Int parent = s.iside(s.cfg.nsidpe, ib);
                    const Int local_face = s.iside(s.cfg.nsidpl, ib);

                    if (parent < 1 || parent > s.cfg.nelem ||
                        local_face < 1 || local_face > s.cfg.nsid)
                    {
                        throw std::runtime_error(
                            "WallModelCoupling - invalid model-wall parent/local-face record");
                    }

                    if (s.mat_elem(parent) != 0)
                    {
                        continue;
                    }

                    const Real area = s.face_norm(4, ib);

                    if (!(area > 0.0) || !std::isfinite(area))
                    {
                        throw std::runtime_error(
                            "WallModelCoupling - invalid model-wall TRI3 area");
                    }

                    std::array<Real, 3> normal =
                    {
                        s.face_norm(1, ib) / area,
                        s.face_norm(2, ib) / area,
                        s.face_norm(3, ib) / area
                    };

                    const Real normal_norm = std::sqrt(
                        normal[0] * normal[0] +
                        normal[1] * normal[1] +
                        normal[2] * normal[2]);

                    if (!(normal_norm > 0.0) || !std::isfinite(normal_norm))
                    {
                        throw std::runtime_error(
                            "WallModelCoupling - invalid model-wall normal");
                    }

                    for (Real& value : normal)
                    {
                        value /= normal_norm;
                    }

                    ModelFace face;
                    face.parent = parent;
                    face.local_face = local_face;
                    face.area = area;
                    face.normal = normal;
                    face.sample_node = s.intma(local_face, parent);

                    for (Int in = 1; in <= 3; ++in)
                    {
                        const Int ip = s.iside(in, ib);

                        if (ip < 1 || ip > s.cfg.npoin)
                        {
                            throw std::runtime_error(
                                "WallModelCoupling - model-wall node outside local mesh");
                        }

                        face.nodes[static_cast<Size>(in - 1)] = ip;

                        if (ip == face.sample_node)
                        {
                            throw std::runtime_error(
                                "WallModelCoupling - opposite/sample node lies on model wall face");
                        }

                        for (Int i = 0; i < 3; ++i)
                        {
                            for (Int j = 0; j < 3; ++j)
                            {
                                gram(i * 3 + j + 1, ip) +=
                                    area * normal[static_cast<Size>(i)] *
                                    normal[static_cast<Size>(j)];
                            }
                        }
                    }

                    if (!(s.detJ(parent) > 0.0) || !std::isfinite(s.detJ(parent)))
                    {
                        throw std::runtime_error(
                            "WallModelCoupling - invalid parent tetrahedron volume");
                    }

                    face.sample_height = 0.5 * s.detJ(parent) / area;

                    if (!(face.sample_height > 0.0) ||
                        !std::isfinite(face.sample_height))
                    {
                        throw std::runtime_error(
                            "WallModelCoupling - invalid off-wall sample height");
                    }

                    if (s.fedge(local_face, parent) == 0)
                    {
                        throw std::runtime_error(
                            "WallModelCoupling - model wall face is absent from the Step-1 boundary flux inventory");
                    }

                    cache.faces.push_back(face);
                }

                if (cache.faces.empty())
                {
                    throw std::runtime_error(
                        "WallModelCoupling - wall treatment requested but no fluid BC 530/532 faces were found");
                }

                cache.global_faces = audit_global_face_uniqueness(s, cache.faces);

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
                            cache.normal_projector[static_cast<Size>(ip)]))
                    {
                        cache.wall_nodes.push_back(ip);
                    }
                }

                if (cache.wall_nodes.empty())
                {
                    throw std::runtime_error(
                        "WallModelCoupling - no model-wall nodes survived projector construction");
                }

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
                        << "  SA production wall treatment: ON\n"
                        << "    model              : Spalding smooth wall\n"
                        << "    kappa / B          : "
                        << options.kappa << " / "
                        << options.log_intercept << "\n"
                        << "    global wall faces  : "
                        << cache.global_faces << "\n"
                        << "    local wall nodes   : "
                        << cache.wall_nodes.size() << "\n"
                        << "    CHT coupling       : DISABLED pending flat-plate + thermal validation\n";
                }
            }

            void ensure_cache(const CBSStateSI& s)
            {
                if (!cache.ready ||
                    cache.state != &s ||
                    cache.npoin != s.cfg.npoin ||
                    cache.nelem != s.cfg.nelem ||
                    cache.nboun != s.cfg.nboun)
                {
                    build_cache(s);
                }
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
                            projector[static_cast<Size>(i * 3 + j)] *
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

            Int dNkdx_index(
                const CBSStateSI& s,
                const Int ie,
                const Int dim,
                const Int local_node)
            {
                return (ie - 1) * s.cfg.ndim * s.cfg.nep
                    + (dim - 1) * s.cfg.nep
                    + local_node;
            }

            Real momentum_diffusivity(
                const CBSStateSI& s,
                const Int ie)
            {
                if (s.cfg.dimensional_mode > 0 &&
                    s.cfg.material_properties_enabled > 0)
                {
                    const Real rho = s.rho_e(ie);
                    const Real mu =
                        s.cfg.turbulence_on > 0
                            ? s.mu_eff_e(ie)
                            : s.mu_e(ie);

                    if (!(rho > 0.0) || !std::isfinite(rho) ||
                        mu < 0.0 || !std::isfinite(mu))
                    {
                        throw std::runtime_error(
                            "WallModelCoupling - invalid effective momentum diffusivity");
                    }

                    return mu / rho;
                }

                if (!(s.cfg.ani > 0.0) || !std::isfinite(s.cfg.ani))
                {
                    throw std::runtime_error(
                        "WallModelCoupling - invalid non-dimensional molecular diffusivity");
                }

                return s.cfg.ani;
            }

            void molecular_properties(
                const CBSStateSI& s,
                const Int ie,
                Real& density,
                Real& molecular_nu)
            {
                if (s.cfg.dimensional_mode > 0 &&
                    s.cfg.material_properties_enabled > 0)
                {
                    density = s.rho_e(ie);
                    const Real mu = s.mu_e(ie);

                    if (!(density > 0.0) || !std::isfinite(density) ||
                        !(mu > 0.0) || !std::isfinite(mu))
                    {
                        throw std::runtime_error(
                            "WallModelCoupling - invalid molecular material properties");
                    }

                    molecular_nu = mu / density;
                    return;
                }

                density = 1.0;
                molecular_nu = s.cfg.ani;

                if (!(molecular_nu > 0.0) || !std::isfinite(molecular_nu))
                {
                    throw std::runtime_error(
                        "WallModelCoupling - invalid non-dimensional molecular viscosity");
                }
            }
        }

        bool WallModelCoupling::enabled(const CBSStateSI& s)
        {
            if (!requested())
            {
                return false;
            }

            ensure_cache(s);
            return true;
        }

        CapturedWallVelocity WallModelCoupling::captureVelocity(
            const CBSStateSI& s,
            const bool owned_only)
        {
            CapturedWallVelocity captured;

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

        void WallModelCoupling::restoreTangentialAndEnforceImpermeability(
            CBSStateSI& s,
            const CapturedWallVelocity& captured)
        {
            if (!enabled(s))
            {
                return;
            }

            if (captured.nodes.size() != captured.values.size())
            {
                throw std::runtime_error(
                    "WallModelCoupling - captured wall velocity arrays are inconsistent");
            }

            for (Size entry = 0; entry < captured.nodes.size(); ++entry)
            {
                const Int ip = captured.nodes[entry];

                if (ip < 1 || ip > s.cfg.npoin ||
                    !isModelWallNode(s, ip))
                {
                    throw std::runtime_error(
                        "WallModelCoupling - invalid captured model-wall node");
                }

                // A moving-wall condition has deliberately higher priority than
                // a stationary model wall at a geometric edge.
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

        void WallModelCoupling::enforceImpermeability(
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

        WallModelMomentumDiagnostics WallModelCoupling::replaceMomentumWallFlux(
            CBSStateSI& s)
        {
            WallModelMomentumDiagnostics diagnostics;

            if (!enabled(s))
            {
                return diagnostics;
            }

            WallTreatmentOptions options;
            options.kappa = environment_real(
                "CBS3D_WALL_KAPPA",
                options.kappa);
            options.log_intercept = environment_real(
                "CBS3D_WALL_B",
                options.log_intercept);

            diagnostics.minimum_y_plus =
                std::numeric_limits<Real>::max();

            for (const ModelFace& face : cache.faces)
            {
                const Int ie = face.parent;
                const Real old_diffusivity = momentum_diffusivity(s, ie);
                Real gradient[3][3] = {};

                for (Int component = 0; component < 3; ++component)
                {
                    for (Int a = 1; a <= s.cfg.nep; ++a)
                    {
                        const Int ip = s.intma(a, ie);

                        for (Int dim = 0; dim < 3; ++dim)
                        {
                            gradient[component][dim] +=
                                s.unkn1(component + 1, ip) *
                                s.dNkdx(
                                    dNkdx_index(
                                        s,
                                        ie,
                                        dim + 1,
                                        a));
                        }
                    }
                }

                // Remove exactly the natural boundary term assembled by the
                // current componentwise diffusion operator.  Use annxf and
                // fdif[2] exactly as MomentumAssembly does, rather than a
                // reconstructed equivalent, so there is no double counting.
                std::array<Real, 3> old_natural{};

                for (Int component = 0; component < 3; ++component)
                {
                    for (Int dim = 0; dim < 3; ++dim)
                    {
                        old_natural[static_cast<Size>(component)] +=
                            gradient[component][dim] *
                            s.annxf(dim + 1, face.local_face, ie) *
                            old_diffusivity * s.cfg.fdif[2];
                    }
                }

                Real density = 1.0;
                Real molecular_nu = 0.0;
                molecular_properties(
                    s,
                    ie,
                    density,
                    molecular_nu);

                const std::array<Real, 3> sample_velocity =
                {
                    s.unkn1(1, face.sample_node),
                    s.unkn1(2, face.sample_node),
                    s.unkn1(3, face.sample_node)
                };

                const WallTreatmentResult wall =
                    WallTreatment::evaluateSpalding(
                        sample_velocity,
                        face.normal,
                        face.sample_height,
                        density,
                        molecular_nu,
                        options);

                const Real face_weight = face.area / 3.0;

                for (const Int ip : face.nodes)
                {
                    for (Int component = 0; component < 3; ++component)
                    {
                        s.rhs(component + 1, ip) -=
                            old_natural[static_cast<Size>(component)];

                        const Real model_contribution =
                            face_weight *
                            wall.kinematic_wall_traction[
                                static_cast<Size>(component)];

                        s.rhs(component + 1, ip) += model_contribution;

                        diagnostics.assembled_nodal_load[
                            static_cast<Size>(component)] +=
                            model_contribution;
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
                        "WallModelCoupling - modeled wall traction performs positive work on the fluid");
                }

                diagnostics.modeled_wall_work += work;
                diagnostics.local_area += face.area;
                ++diagnostics.local_faces;

                if (wall.wall_shear_magnitude > 0.0)
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
                const Real a = diagnostics.modeled_surface_load[
                    static_cast<Size>(component)];
                const Real b = diagnostics.assembled_nodal_load[
                    static_cast<Size>(component)];
                const Real tolerance =
                    1.0e-12 * std::max({1.0, std::abs(a), std::abs(b)});

                if (std::abs(a - b) > tolerance)
                {
                    throw std::runtime_error(
                        "WallModelCoupling - TRI3 wall-load conservation check failed");
                }
            }

            return diagnostics;
        }

        long long WallModelCoupling::globalWallFaceCount(const CBSStateSI& s)
        {
            if (!enabled(s))
            {
                return 0;
            }

            return cache.global_faces;
        }

        bool WallModelCoupling::isModelWallNode(
            const CBSStateSI& s,
            const Int ip)
        {
            if (!enabled(s))
            {
                return false;
            }

            if (ip < 1 || ip > s.cfg.npoin)
            {
                return false;
            }

            return nonzero_projector(
                cache.normal_projector[static_cast<Size>(ip)]);
        }
    }
}
