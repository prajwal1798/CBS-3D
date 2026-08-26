//=============================================================================
// CBS3D++_SI
//
// Verified P1/TET4 Spalart-Allmaras finite-element assembly.
//
// This implementation replaces the historical centroid treatment of the
// nonlinear SA production/destruction source by a positive four-point
// tetrahedral quadrature rule.  Linear P1 operators remain integrated in their
// exact closed forms.
//
// Important discrete choices:
//   * q = nu_tilde is P1 on each tetrahedron.
//   * grad(q) and grad(u) are element-constant.
//   * Galerkin advection is integrated exactly with the P1 consistent mass
//     coefficients V/10 and V/20.
//   * div[(nu+q) grad(q)] is integrated exactly for constant element molecular
//     viscosity by using the exact element mean q_bar.
//   * P(q,d) and D(q,d) are nonlinear and are evaluated at four positive
//     degree-2 tetrahedron quadrature points, including the test function N_a.
//   * positive destruction is linearised as C_D(q_old,d) q_old q_new and then
//     row-sum lumped: L_a = integral N_a C_D q_old dV.
//   * the CBS scalar characteristic term includes both its element-volume and
//     physical-boundary contributions, matching the transferred momentum form.
//   * nodal nu_t/mu_t output is a volume-weighted projection of the element
//     quantities instead of an unweighted incident-element count average.
//=============================================================================

#include "cbs/assembly/SpalartAllmarasAssembly.hpp"

#include "cbs/parallel/HaloExchange.hpp"
#include "cbs/turbulence/SpalartAllmaras.hpp"

#ifdef CBS3D_USE_MPI
#include <mpi.h>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <ios>
#include <sstream>
#include <stdexcept>
#include <string>

namespace cbs
{
    namespace
    {
        constexpr Real tetra_q_a = 0.58541019662496845446;
        constexpr Real tetra_q_b = 0.13819660112501051518;
        constexpr Real tetra_q_w = 0.25;

        bool is_fluid_element(const CBSStateSI& s, const Int ie)
        {
            return s.mat_elem(ie) == 0;
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

        Real grad(
            const CBSStateSI& s,
            const Int ie,
            const Int dim,
            const Int local_node)
        {
            return s.dNkdx(dNkdx_index(s, ie, dim, local_node));
        }

        void validate_transport_dimensions(const CBSStateSI& s)
        {
            if (s.cfg.ndim != 3 || s.cfg.nep != 4 ||
                s.cfg.nsid != 4 || s.cfg.nsidp != 3)
            {
                throw std::runtime_error(
                    "SpalartAllmarasAssembly - SA transport requires ndim=3, nep=4, nsid=4, nsidp=3");
            }
        }

        Real molecular_kinematic_viscosity(
            const CBSStateSI& s,
            const Int ie,
            bool& ok)
        {
            ok = true;

            if (s.cfg.dimensional_mode > 0 &&
                s.cfg.material_properties_enabled > 0)
            {
                const Real rho = s.rho_e(ie);
                const Real mu = s.mu_e(ie);

                if (!(rho > 0.0) || !std::isfinite(rho) ||
                    mu < 0.0 || !std::isfinite(mu))
                {
                    ok = false;
                    return 0.0;
                }

                return mu / rho;
            }

            if (!(s.cfg.ani > 0.0) || !std::isfinite(s.cfg.ani))
            {
                ok = false;
                return 0.0;
            }

            return s.cfg.ani;
        }

        Real element_vorticity_magnitude(
            const CBSStateSI& s,
            const Int ie)
        {
            Real duidxj[4][4] = {};

            for (Int i = 1; i <= 3; ++i)
            {
                for (Int a = 1; a <= 4; ++a)
                {
                    const Int ip = s.intma(a, ie);

                    for (Int j = 1; j <= 3; ++j)
                    {
                        duidxj[i][j] +=
                            s.unkno(i, ip) * grad(s, ie, j, a);
                    }
                }
            }

            const Real wx = duidxj[3][2] - duidxj[2][3];
            const Real wy = duidxj[1][3] - duidxj[3][1];
            const Real wz = duidxj[2][1] - duidxj[1][2];

            return std::sqrt(wx * wx + wy * wy + wz * wz);
        }

        void compute_sa_gradient(
            const CBSStateSI& s,
            const Int ie,
            const Real q_node[5],
            Real grad_q[4])
        {
            grad_q[1] = 0.0;
            grad_q[2] = 0.0;
            grad_q[3] = 0.0;

            for (Int a = 1; a <= 4; ++a)
            {
                for (Int j = 1; j <= 3; ++j)
                {
                    grad_q[j] += q_node[a] * grad(s, ie, j, a);
                }
            }
        }

        void compute_velocity_sum(
            const CBSStateSI& s,
            const Int ie,
            Real velocity_sum[4])
        {
            velocity_sum[1] = 0.0;
            velocity_sum[2] = 0.0;
            velocity_sum[3] = 0.0;

            for (Int a = 1; a <= 4; ++a)
            {
                const Int ip = s.intma(a, ie);
                velocity_sum[1] += s.unkno(1, ip);
                velocity_sum[2] += s.unkno(2, ip);
                velocity_sum[3] += s.unkno(3, ip);
            }
        }

        std::array<Real, 4> quadrature_shape_values(const Int qp)
        {
            std::array<Real, 4> n =
            {
                tetra_q_b,
                tetra_q_b,
                tetra_q_b,
                tetra_q_b
            };

            n[static_cast<Size>(qp)] = tetra_q_a;
            return n;
        }

        std::string describe_sa_node_failure(
            const CBSStateSI& s,
            const Int ip,
            const char* reason)
        {
            std::ostringstream text;
            text.setf(std::ios::scientific);
            text.precision(8);
            text << "SpalartAllmarasAssembly - " << reason
                 << "\n    rank                 " << s.mpi_rank
                 << "\n    local node           " << ip
                 << "\n    coordinates          ("
                 << s.coord(1, ip) << ", "
                 << s.coord(2, ip) << ", "
                 << s.coord(3, ip) << ")"
                 << "\n    wall distance        " << s.wall_distance(ip)
                 << "\n    nu_tilde old         " << s.nu_tilde1(ip)
                 << "\n    nu_tilde new         " << s.nu_tilde(ip)
                 << "\n    sa_rhs               " << s.sa_rhs(ip)
                 << "\n    sa_production        " << s.sa_production(ip)
                 << "\n    sa_destruction       " << s.sa_destruction(ip)
                 << "\n    sa_diffusion         " << s.sa_diffusion(ip)
                 << "\n    sa_source            " << s.sa_source(ip)
                 << "\n    sa_destruction_lhs   " << s.sa_destruction_lhs(ip)
                 << "\n    elcoe2               " << s.elcoe2(ip);
            return text.str();
        }
    }

    void SpalartAllmarasAssembly::resetEffectiveProperties(CBSStateSI& s)
    {
#pragma omp parallel for schedule(static)
        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            s.nu_tilde_e(ie) = 0.0;
            s.nu_t_e(ie) = 0.0;
            s.mu_t_e(ie) = 0.0;
            s.mu_eff_e(ie) = s.mu_e(ie);
            s.k_eff_e(ie) = s.k_e(ie);
        }

        s.nu_t.fill(0.0);
        s.mu_t.fill(0.0);
    }

    void SpalartAllmarasAssembly::updateEddyViscosity(CBSStateSI& s)
    {
        resetEffectiveProperties(s);

        if (s.cfg.turbulence_on < 1)
        {
            return;
        }

        turbulence::SpalartAllmarasConstants constants;
        s.nu_t.fill(0.0);
        s.mu_t.fill(0.0);
        s.sa_nodal_weight.fill(0.0);

        int bad_material = 0;
        Int diverged_element = 0;
        Real diverged_ratio = 0.0;

        for (Int colour = 0; colour < s.ncolor; ++colour)
        {
            const Int cbeg = s.color_ptr[static_cast<Size>(colour)];
            const Int cend = s.color_ptr[static_cast<Size>(colour) + 1];

#pragma omp parallel for schedule(static)
            for (Int k = cbeg; k < cend; ++k)
            {
                const Int ie = s.color_elem[static_cast<Size>(k)];

                if (!is_fluid_element(s, ie))
                {
                    continue;
                }

                if (!(s.detJ(ie) > 0.0) || !std::isfinite(s.detJ(ie)))
                {
#pragma omp atomic write
                    bad_material = 1;
                    continue;
                }

                Real q_average = 0.0;
                for (Int a = 1; a <= 4; ++a)
                {
                    const Int ip = s.intma(a, ie);
                    q_average += s.sa_wall_node(ip) != 0
                        ? 0.0
                        : std::max(s.cfg.sa_nu_tilde_floor, s.nu_tilde(ip));
                }
                q_average *= 0.25;

                bool material_ok = true;
                const Real molecular_nu =
                    molecular_kinematic_viscosity(s, ie, material_ok);

                if (!material_ok)
                {
#pragma omp atomic write
                    bad_material = 1;
                    continue;
                }

                if (s.cfg.sa_nu_tilde_ceiling_ratio > 0.0 &&
                    molecular_nu > 0.0 &&
                    q_average > s.cfg.sa_nu_tilde_ceiling_ratio * molecular_nu)
                {
#pragma omp critical(sa_ceiling_failure)
                    {
                        if (diverged_element == 0)
                        {
                            diverged_element = ie;
                            diverged_ratio = q_average / molecular_nu;
                        }
                    }
                }

                const Real nu_t = turbulence::eddyKinematicViscosity(
                    q_average,
                    molecular_nu,
                    constants);
                const Real mu_t = s.rho_e(ie) * nu_t;

                s.nu_tilde_e(ie) = q_average;
                s.nu_t_e(ie) = nu_t;
                s.mu_t_e(ie) = mu_t;
                s.mu_eff_e(ie) = s.mu_e(ie) + mu_t;

                if (s.cfg.turbulent_thermal_diffusivity_on > 0)
                {
                    s.k_eff_e(ie) = s.k_e(ie)
                        + s.rho_cp_e(ie) * nu_t / s.cfg.sa_prandtl_t;
                }

                // Volume-weighted projection is partition independent and does
                // not give a tiny near-wall tetrahedron the same nodal weight as
                // a much larger neighbour.
                const Real element_volume = s.detJ(ie) / 6.0;
                for (Int a = 1; a <= 4; ++a)
                {
                    const Int ip = s.intma(a, ie);
                    s.nu_t(ip) += element_volume * nu_t;
                    s.mu_t(ip) += element_volume * mu_t;
                    s.sa_nodal_weight(ip) += element_volume;
                }
            }
        }

        if (bad_material != 0)
        {
            throw std::runtime_error(
                "SpalartAllmarasAssembly::updateEddyViscosity - invalid fluid material or element geometry");
        }

        if (diverged_element != 0)
        {
            std::ostringstream text;
            text << "SpalartAllmarasAssembly::updateEddyViscosity - nu_tilde/nu = "
                 << diverged_ratio
                 << " exceeded sa_nu_tilde_ceiling_ratio at local element "
                 << diverged_element;
            throw std::runtime_error(text.str());
        }

#ifdef CBS3D_USE_MPI
        if (s.mpi_enabled && s.mpi_size > 1)
        {
            HaloExchange::sumGhostContributionsToOwners(
                s.nu_t, s.partition_metadata, MPI_COMM_WORLD);
            HaloExchange::sumGhostContributionsToOwners(
                s.mu_t, s.partition_metadata, MPI_COMM_WORLD);
            HaloExchange::sumGhostContributionsToOwners(
                s.sa_nodal_weight, s.partition_metadata, MPI_COMM_WORLD);
        }
#endif

#pragma omp parallel for schedule(static)
        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            const Real weight = s.sa_nodal_weight(ip);
            if (weight > 0.0)
            {
                s.nu_t(ip) /= weight;
                s.mu_t(ip) /= weight;
            }
            else
            {
                s.nu_t(ip) = 0.0;
                s.mu_t(ip) = 0.0;
            }
        }

#ifdef CBS3D_USE_MPI
        if (s.mpi_enabled && s.mpi_size > 1)
        {
            HaloExchange::broadcastOwnedToGhosts(
                s.nu_t, s.partition_metadata, MPI_COMM_WORLD);
            HaloExchange::broadcastOwnedToGhosts(
                s.mu_t, s.partition_metadata, MPI_COMM_WORLD);
        }
#endif
    }

    void SpalartAllmarasAssembly::assembleTransportRhs(CBSStateSI& s)
    {
        validate_transport_dimensions(s);

        s.sa_rhs.fill(0.0);
        s.sa_source.fill(0.0);
        s.sa_production.fill(0.0);
        s.sa_destruction.fill(0.0);
        s.sa_diffusion.fill(0.0);
        s.sa_residual.fill(0.0);
        s.sa_destruction_lhs.fill(0.0);

        if (s.cfg.turbulence_on < 1)
        {
            return;
        }

        turbulence::SpalartAllmarasConstants constants;

        int bad_detj = 0;
        int bad_material = 0;
        int bad_wall_distance = 0;
        int bad_dt = 0;

        for (Int colour = 0; colour < s.ncolor; ++colour)
        {
            const Int cbeg = s.color_ptr[static_cast<Size>(colour)];
            const Int cend = s.color_ptr[static_cast<Size>(colour) + 1];

#pragma omp parallel for schedule(static)
            for (Int k = cbeg; k < cend; ++k)
            {
                const Int ie = s.color_elem[static_cast<Size>(k)];

                if (!is_fluid_element(s, ie))
                {
                    continue;
                }

                if (!(s.detJ(ie) > 0.0) || !std::isfinite(s.detJ(ie)))
                {
#pragma omp atomic write
                    bad_detj = 1;
                    continue;
                }

                const Real volume = s.detJ(ie) / 6.0;
                const Real nodal_volume = volume / 4.0;
                const Real advective_mass_factor = volume / 20.0;

                Real q_node[5] = {};
                Real d_node[5] = {};
                Real q_sum = 0.0;

                for (Int a = 1; a <= 4; ++a)
                {
                    const Int ip = s.intma(a, ie);
                    q_node[a] = s.sa_wall_node(ip) != 0
                        ? 0.0
                        : std::max(s.cfg.sa_nu_tilde_floor, s.nu_tilde1(ip));
                    q_sum += q_node[a];

                    d_node[a] = s.wall_distance(ip);
                    if (!(d_node[a] > 0.0) || !std::isfinite(d_node[a]))
                    {
#pragma omp atomic write
                        bad_wall_distance = 1;
                    }
                }

                if (bad_wall_distance != 0)
                {
                    continue;
                }

                const Real q_average = 0.25 * q_sum;

                bool material_ok = true;
                const Real molecular_nu =
                    molecular_kinematic_viscosity(s, ie, material_ok);
                if (!material_ok)
                {
#pragma omp atomic write
                    bad_material = 1;
                    continue;
                }

                Real grad_q[4] = {};
                compute_sa_gradient(s, ie, q_node, grad_q);
                const Real grad_q_squared =
                    grad_q[1] * grad_q[1] +
                    grad_q[2] * grad_q[2] +
                    grad_q[3] * grad_q[3];

                Real velocity_sum[4] = {};
                compute_velocity_sum(s, ie, velocity_sum);
                Real mean_velocity[4] = {};
                for (Int j = 1; j <= 3; ++j)
                {
                    mean_velocity[j] = 0.25 * velocity_sum[j];
                }

                const Real omega = element_vorticity_magnitude(s, ie);
                const Real diffusion_coefficient =
                    (molecular_nu + q_average) / constants.sigma;
                const Real nonlinear_gradient_source =
                    constants.cb2 * grad_q_squared / constants.sigma;

                const Real element_dt = s.delte(ie);
                if (!(element_dt > 0.0) || !std::isfinite(element_dt))
                {
#pragma omp atomic write
                    bad_dt = 1;
                    continue;
                }

                const Real advective_derivative =
                    mean_velocity[1] * grad_q[1] +
                    mean_velocity[2] * grad_q[2] +
                    mean_velocity[3] * grad_q[3];
                const Real characteristic_factor =
                    0.5 * element_dt * volume;

                Real local_rhs[5] = {};
                Real local_production[5] = {};
                Real local_destruction[5] = {};
                Real local_diffusion[5] = {};
                Real local_source[5] = {};
                Real local_destruction_lhs[5] = {};

                // Exact P1 Galerkin advection, exact linear diffusion, and the
                // element-volume half of the CBS characteristic correction.
                for (Int a = 1; a <= 4; ++a)
                {
                    const Int ip = s.intma(a, ie);
                    const Real weighted_u = velocity_sum[1] + s.unkno(1, ip);
                    const Real weighted_v = velocity_sum[2] + s.unkno(2, ip);
                    const Real weighted_w = velocity_sum[3] + s.unkno(3, ip);

                    const Real advection = -advective_mass_factor *
                        (weighted_u * grad_q[1] +
                         weighted_v * grad_q[2] +
                         weighted_w * grad_q[3]);

                    const Real diffusion = -volume * diffusion_coefficient *
                        (grad(s, ie, 1, a) * grad_q[1] +
                         grad(s, ie, 2, a) * grad_q[2] +
                         grad(s, ie, 3, a) * grad_q[3]);

                    const Real nonlinear_rhs =
                        nodal_volume * nonlinear_gradient_source;

                    const Real characteristic =
                        -characteristic_factor * advective_derivative *
                        (mean_velocity[1] * grad(s, ie, 1, a) +
                         mean_velocity[2] * grad(s, ie, 2, a) +
                         mean_velocity[3] * grad(s, ie, 3, a));

                    local_rhs[a] +=
                        advection + diffusion + nonlinear_rhs + characteristic;
                    local_diffusion[a] += diffusion + nonlinear_rhs;
                }

                // Complete the transferred CBS characteristic weak form on
                // physical faces. annxf is the area-weighted outward normal and
                // fcon[4] is the same TRI3 integration factor used by Step 1.
                for (Int is = 1; is <= s.cfg.nsid; ++is)
                {
                    if (s.fedge(is, ie) == 0)
                    {
                        continue;
                    }

                    Real u_dot_area_normal = 0.0;
                    for (Int j = 1; j <= 3; ++j)
                    {
                        u_dot_area_normal +=
                            mean_velocity[j] * s.annxf(j, is, ie);
                    }

                    const Real face_characteristic =
                        0.5 * element_dt * s.cfg.fcon[4] *
                        u_dot_area_normal * advective_derivative;

                    for (Int face_node = 1;
                         face_node <= s.cfg.nsidp;
                         ++face_node)
                    {
                        const Int a = s.ippn1(is, face_node);
                        local_rhs[a] += face_characteristic;
                    }
                }

                // Positive four-point quadrature for the nonlinear SA source.
                // The rule is exact through polynomial degree two and, unlike a
                // centroid evaluation, samples the near-wall q/d variation.
                for (Int qp = 0; qp < 4; ++qp)
                {
                    const std::array<Real, 4> n =
                        quadrature_shape_values(qp);

                    Real q_qp = 0.0;
                    Real d_qp = 0.0;
                    for (Int a = 1; a <= 4; ++a)
                    {
                        q_qp += n[static_cast<Size>(a - 1)] * q_node[a];
                        d_qp += n[static_cast<Size>(a - 1)] * d_node[a];
                    }
                    d_qp = std::max(d_qp, s.cfg.sa_min_wall_distance);

                    Real s_bar = 0.0;
                    if (q_qp > 0.0)
                    {
                        s_bar = turbulence::sBar(
                            q_qp, molecular_nu, d_qp, constants);
                    }

                    Real s_tilde = s.cfg.sa_use_stilde_limiter > 0
                        ? turbulence::limitedSTilde(omega, s_bar, constants)
                        : omega + s_bar;
                    s_tilde = std::max(s_tilde, s.cfg.sa_min_stilde);

                    const Real chi_value =
                        turbulence::chi(q_qp, molecular_nu);
                    const Real ft2_value =
                        turbulence::ft2(chi_value, constants);
                    const Real r_value = turbulence::rFunction(
                        q_qp, s_tilde, d_qp, constants);
                    const Real fw_value =
                        turbulence::fw(r_value, constants);
                    const Real production = turbulence::productionTerm(
                        q_qp, s_tilde, ft2_value, constants);
                    const Real destruction_coefficient =
                        turbulence::destructionCoefficient(
                            d_qp, fw_value, ft2_value, constants);
                    const Real destruction =
                        destruction_coefficient * q_qp * q_qp;

                    const Real qp_volume = volume * tetra_q_w;

                    for (Int a = 1; a <= 4; ++a)
                    {
                        const Real test_weight =
                            qp_volume * n[static_cast<Size>(a - 1)];
                        const Real production_rhs =
                            test_weight * production;
                        const Real destruction_rhs =
                            test_weight * destruction;

                        local_rhs[a] += production_rhs;

                        if (s.cfg.sa_implicit_destruction > 0 &&
                            destruction_coefficient > 0.0 &&
                            q_qp > 0.0)
                        {
                            // Row-sum lumping of
                            // integral N_a C_D q_old N_b dV.
                            local_destruction_lhs[a] +=
                                test_weight * destruction_coefficient * q_qp;
                        }
                        else
                        {
                            local_rhs[a] -= destruction_rhs;
                        }

                        local_production[a] += production_rhs;
                        local_destruction[a] += destruction_rhs;
                        local_source[a] += production_rhs - destruction_rhs;
                    }
                }

                for (Int a = 1; a <= 4; ++a)
                {
                    const Int ip = s.intma(a, ie);
                    s.sa_rhs(ip) += local_rhs[a];
                    s.sa_production(ip) += local_production[a];
                    s.sa_destruction(ip) += local_destruction[a];
                    s.sa_diffusion(ip) += local_diffusion[a];
                    s.sa_source(ip) += local_source[a];
                    s.sa_destruction_lhs(ip) += local_destruction_lhs[a];
                }
            }
        }

        if (bad_detj != 0)
        {
            throw std::runtime_error(
                "SpalartAllmarasAssembly::assembleTransportRhs - invalid detJ");
        }
        if (bad_material != 0)
        {
            throw std::runtime_error(
                "SpalartAllmarasAssembly::assembleTransportRhs - invalid fluid material properties");
        }
        if (bad_wall_distance != 0)
        {
            throw std::runtime_error(
                "SpalartAllmarasAssembly::assembleTransportRhs - invalid wall distance");
        }
        if (bad_dt != 0)
        {
            throw std::runtime_error(
                "SpalartAllmarasAssembly::assembleTransportRhs - invalid element timestep");
        }
    }

    void SpalartAllmarasAssembly::updateNuTilde(CBSStateSI& s)
    {
        if (s.cfg.turbulence_on < 1)
        {
            return;
        }

        Int failed_node = 0;
        const char* failure_reason = "";

        const bool have_owned = !s.owned_nodes.empty();
        const Int update_count = have_owned
            ? static_cast<Int>(s.owned_nodes.size())
            : s.cfg.npoin;

#pragma omp parallel for schedule(static)
        for (Int index = 0; index < update_count; ++index)
        {
            const Int ip = have_owned
                ? s.owned_nodes[static_cast<Size>(index)]
                : index + 1;

            if (s.sa_active_node(ip) == 0)
            {
                s.nu_tilde(ip) = 0.0;
                s.sa_residual(ip) = 0.0;
                continue;
            }

            if (s.sa_wall_node(ip) != 0)
            {
                s.nu_tilde(ip) = 0.0;
                s.sa_residual(ip) = -s.nu_tilde1(ip);
                continue;
            }

            if (!(s.elcoe2(ip) > 0.0) || !std::isfinite(s.elcoe2(ip)))
            {
                s.nu_tilde(ip) = s.nu_tilde1(ip);
                s.sa_residual(ip) = 0.0;
                continue;
            }

            const Real old_value = s.nu_tilde1(ip);
            Real denominator = 1.0;
            if (s.cfg.sa_implicit_destruction > 0)
            {
                denominator +=
                    s.elcoe2(ip) * s.sa_destruction_lhs(ip);
            }

            if (!(denominator > 0.0) || !std::isfinite(denominator))
            {
#pragma omp critical(sa_update_failure)
                {
                    if (failed_node == 0)
                    {
                        failed_node = ip;
                        failure_reason =
                            "invalid semi-implicit destruction denominator";
                    }
                }
                continue;
            }

            Real new_value =
                (old_value + s.elcoe2(ip) * s.sa_rhs(ip)) / denominator;

            if (!std::isfinite(new_value))
            {
#pragma omp critical(sa_update_failure)
                {
                    if (failed_node == 0)
                    {
                        failed_node = ip;
                        failure_reason = "non-finite nu_tilde update";
                    }
                }
                continue;
            }

            new_value = std::max(s.cfg.sa_nu_tilde_floor, new_value);
            s.nu_tilde(ip) = new_value;
            s.sa_residual(ip) = new_value - old_value;
        }

        if (failed_node != 0)
        {
            throw std::runtime_error(
                describe_sa_node_failure(s, failed_node, failure_reason));
        }
    }
}
