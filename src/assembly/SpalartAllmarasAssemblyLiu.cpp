//=============================================================================
// CBS3D++_SI
//
// Spalart-Allmaras assembly aligned with Chun-Bin Liu (2005), Chapters 3-4 and
// Appendix B, for the Swansea/Nithiarasu CBS formulation.
//
// The production discretisation implemented here follows the thesis matrices:
//
//   M_nu      Appendix B Eq. B.48
//   C_nu      Appendix B Eqs. B.50-B.51
//   K_nu      Appendix B Eqs. B.52-B.53
//   f_nuOmega Appendix B Eqs. B.54-B.55
//   source    Appendix B Eqs. B.56-B.57
//   K_u_nu    Appendix B Eqs. B.60-B.61
//
// Consequences that differ from the historical C++ SA port:
//
//   * the transported convection is conservative: div(u nu_tilde), not merely
//     u.grad(nu_tilde);
//   * the CBS characteristic term is the exact P1/TET4 K_u_nu product using
//     div(u N_a), not an element-mean-velocity surrogate;
//   * the production/destruction vector uses the consistent P1 mass weighting
//     of Appendix B, not equal V/4 centroid forcing and not an unrelated
//     four-point nonlinear source quadrature;
//   * the Liu model is the fully turbulent SA-noft2 form shown in Eq. 3.23;
//   * destruction remains in the explicit Step-4 residual, as in Eq. 4.25;
//   * wall nu_tilde=0 is imposed strongly.  The present outlet uses homogeneous
//     natural SA flux, so the generic f_nuGamma term of Eq. B.58 is zero.
//
// Modern owner/ghost exchange, finite-value checks and the final non-negative
// projection are retained as implementation safeguards; they do not redefine
// the finite-element residual above.
//=============================================================================

#include "cbs/assembly/SpalartAllmarasAssembly.hpp"
#include "cbs/assembly/LiuNithiarasuKernels.hpp"
#include "cbs/parallel/HaloExchange.hpp"
#include "cbs/turbulence/SpalartAllmaras.hpp"

#ifdef CBS3D_USE_MPI
#include <mpi.h>
#endif

#include <algorithm>
#include <cmath>
#include <ios>
#include <sstream>
#include <stdexcept>
#include <string>

namespace cbs
{
    namespace
    {
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

        void validate_dimensions(const CBSStateSI& s)
        {
            if (s.cfg.ndim != 3 || s.cfg.nep != 4 ||
                s.cfg.nsid != 4 || s.cfg.nsidp != 3)
            {
                throw std::runtime_error(
                    "SpalartAllmarasAssembly - Liu/Nithiarasu SA requires "
                    "3-D P1 TET4/TRI3 topology");
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

        Real vorticity_magnitude(const CBSStateSI& s, const Int ie)
        {
            Real du[4][4] = {};
            for (Int i = 1; i <= 3; ++i)
            {
                for (Int a = 1; a <= 4; ++a)
                {
                    const Int ip = s.intma(a, ie);
                    for (Int j = 1; j <= 3; ++j)
                    {
                        du[i][j] += s.unkno(i, ip) * grad(s, ie, j, a);
                    }
                }
            }

            const Real wx = du[3][2] - du[2][3];
            const Real wy = du[1][3] - du[3][1];
            const Real wz = du[2][1] - du[1][2];
            return std::sqrt(wx * wx + wy * wy + wz * wz);
        }

        void compute_sa_gradient(
            const CBSStateSI& s,
            const Int ie,
            const Real q[5],
            Real grad_q[4])
        {
            for (Int j = 1; j <= 3; ++j)
            {
                grad_q[j] = 0.0;
                for (Int a = 1; a <= 4; ++a)
                {
                    grad_q[j] += q[a] * grad(s, ie, j, a);
                }
            }
        }

        std::string describe_node_failure(
            const CBSStateSI& s,
            const Int ip,
            const char* reason)
        {
            std::ostringstream text;
            text.setf(std::ios::scientific);
            text.precision(8);
            text << "SpalartAllmarasAssembly - " << reason
                 << "\n  rank           " << s.mpi_rank
                 << "\n  node           " << ip
                 << "\n  xyz            "
                 << s.coord(1, ip) << " "
                 << s.coord(2, ip) << " "
                 << s.coord(3, ip)
                 << "\n  wall distance  " << s.wall_distance(ip)
                 << "\n  old nu_tilde   " << s.nu_tilde1(ip)
                 << "\n  rhs            " << s.sa_rhs(ip)
                 << "\n  production     " << s.sa_production(ip)
                 << "\n  destruction    " << s.sa_destruction(ip)
                 << "\n  diffusion      " << s.sa_diffusion(ip)
                 << "\n  dt/M           " << s.elcoe2(ip);
            return text.str();
        }
    }

    void SpalartAllmarasAssembly::resetEffectiveProperties(CBSStateSI& s)
    {
#ifdef CBS3D_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
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
        s.sa_nodal_weight.fill(0.0);

        int bad_material = 0;
        Int diverged_element = 0;
        Real diverged_ratio = 0.0;

        for (Int colour = 0; colour < s.ncolor; ++colour)
        {
            const Int cbeg = s.color_ptr[static_cast<Size>(colour)];
            const Int cend = s.color_ptr[static_cast<Size>(colour) + 1];

#ifdef CBS3D_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (Int k = cbeg; k < cend; ++k)
            {
                const Int ie = s.color_elem[static_cast<Size>(k)];
                if (!is_fluid_element(s, ie))
                {
                    continue;
                }

                if (!(s.detJ(ie) > 0.0) || !std::isfinite(s.detJ(ie)))
                {
#ifdef CBS3D_USE_OPENMP
#pragma omp atomic write
#endif
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
#ifdef CBS3D_USE_OPENMP
#pragma omp atomic write
#endif
                    bad_material = 1;
                    continue;
                }

                if (s.cfg.sa_nu_tilde_ceiling_ratio > 0.0 &&
                    molecular_nu > 0.0 &&
                    q_average > s.cfg.sa_nu_tilde_ceiling_ratio * molecular_nu)
                {
#ifdef CBS3D_USE_OPENMP
#pragma omp critical(sa_liu_ceiling)
#endif
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

                const Real volume = s.detJ(ie) / 6.0;
                for (Int a = 1; a <= 4; ++a)
                {
                    const Int ip = s.intma(a, ie);
                    s.nu_t(ip) += volume * nu_t;
                    s.mu_t(ip) += volume * mu_t;
                    s.sa_nodal_weight(ip) += volume;
                }
            }
        }

        if (bad_material != 0)
        {
            throw std::runtime_error(
                "SpalartAllmarasAssembly::updateEddyViscosity - invalid material or TET4 geometry");
        }
        if (diverged_element != 0)
        {
            throw std::runtime_error(
                "SpalartAllmarasAssembly::updateEddyViscosity - nu_tilde/nu ceiling exceeded at element "
                + std::to_string(diverged_element)
                + ", ratio=" + std::to_string(diverged_ratio));
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

#ifdef CBS3D_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
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
        validate_dimensions(s);

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

#ifdef CBS3D_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (Int k = cbeg; k < cend; ++k)
            {
                const Int ie = s.color_elem[static_cast<Size>(k)];
                if (!is_fluid_element(s, ie))
                {
                    continue;
                }

                if (!(s.detJ(ie) > 0.0) || !std::isfinite(s.detJ(ie)))
                {
#ifdef CBS3D_USE_OPENMP
#pragma omp atomic write
#endif
                    bad_detj = 1;
                    continue;
                }

                const Real volume = s.detJ(ie) / 6.0;
                const Real nodal_volume = volume / 4.0;

                Real q[5] = {};
                Real velocity[5][4] = {};
                Real q_sum = 0.0;
                Real d_sum = 0.0;

                bool local_distance_ok = true;
                for (Int a = 1; a <= 4; ++a)
                {
                    const Int ip = s.intma(a, ie);
                    q[a] = s.sa_wall_node(ip) != 0
                        ? 0.0
                        : std::max(s.cfg.sa_nu_tilde_floor, s.nu_tilde1(ip));
                    q_sum += q[a];

                    const Real da = s.wall_distance(ip);
                    if (!(da > 0.0) || !std::isfinite(da))
                    {
                        local_distance_ok = false;
                    }
                    d_sum += da;

                    for (Int i = 1; i <= 3; ++i)
                    {
                        velocity[a][i] = s.unkno(i, ip);
                    }
                }

                if (!local_distance_ok)
                {
#ifdef CBS3D_USE_OPENMP
#pragma omp atomic write
#endif
                    bad_wall_distance = 1;
                    continue;
                }

                const Real q_average = 0.25 * q_sum;
                const Real wall_distance = std::max(
                    0.25 * d_sum,
                    s.cfg.sa_min_wall_distance);

                bool material_ok = true;
                const Real molecular_nu =
                    molecular_kinematic_viscosity(s, ie, material_ok);
                if (!material_ok)
                {
#ifdef CBS3D_USE_OPENMP
#pragma omp atomic write
#endif
                    bad_material = 1;
                    continue;
                }

                const Real dt = s.delte(ie);
                if (!(dt > 0.0) || !std::isfinite(dt))
                {
#ifdef CBS3D_USE_OPENMP
#pragma omp atomic write
#endif
                    bad_dt = 1;
                    continue;
                }

                Real grad_q[4] = {};
                compute_sa_gradient(s, ie, q, grad_q);
                const Real grad_q_squared =
                    grad_q[1] * grad_q[1]
                    + grad_q[2] * grad_q[2]
                    + grad_q[3] * grad_q[3];

                Real div_u = 0.0;
                for (Int a = 1; a <= 4; ++a)
                {
                    for (Int i = 1; i <= 3; ++i)
                    {
                        div_u += velocity[a][i] * grad(s, ie, i, a);
                    }
                }

                Real div_uq[5] = {};
                liu_nithiarasu::conservative_scalar_divergence_nodes(
                    velocity,
                    q,
                    grad_q,
                    div_u,
                    div_uq);

                const Real omega = vorticity_magnitude(s, ie);
                const Real s_bar = q_average > 0.0
                    ? turbulence::sBar(
                        q_average,
                        molecular_nu,
                        wall_distance,
                        constants)
                    : 0.0;

                // Liu Eq. 3.26.  No ft2 term and no AJS S-tilde branch appears
                // in the thesis formulation.  The tiny floor is only a division
                // guard for r and does not act in ordinary attached flow.
                const Real s_tilde = std::max(
                    omega + s_bar,
                    s.cfg.sa_min_stilde);

                const Real r_value = turbulence::rFunction(
                    q_average,
                    s_tilde,
                    wall_distance,
                    constants);
                const Real fw_value = turbulence::fw(r_value, constants);

                // Passing ft2=0 selects the fully turbulent Liu/SA-noft2 model
                // while retaining the well-tested scalar helper functions.
                const Real production = turbulence::productionTerm(
                    q_average,
                    s_tilde,
                    0.0,
                    constants);
                const Real destruction_coefficient =
                    turbulence::destructionCoefficient(
                        wall_distance,
                        fw_value,
                        0.0,
                        constants);

                const Real diffusion_coefficient =
                    (molecular_nu + q_average) / constants.sigma;
                const Real nonlinear_source =
                    constants.cb2 * grad_q_squared / constants.sigma;

                Real local_rhs[5] = {};
                Real local_production[5] = {};
                Real local_destruction[5] = {};
                Real local_diffusion[5] = {};
                Real local_source[5] = {};

                for (Int a = 1; a <= 4; ++a)
                {
                    // C_nu q: -int N_a div(u q) dV.
                    const Real advection =
                        -liu_nithiarasu::consistent_linear_load(
                            volume,
                            div_uq,
                            a);

                    // K_nu q: -int grad(N_a) . [(nu+qbar)/sigma grad(q)] dV.
                    const Real diffusion = -volume * diffusion_coefficient *
                        (grad(s, ie, 1, a) * grad_q[1]
                         + grad(s, ie, 2, a) * grad_q[2]
                         + grad(s, ie, 3, a) * grad_q[3]);

                    // f_nuOmega, second diffusion term.  grad(q) is constant on
                    // a P1 tetrahedron, so V/4 is exact.
                    const Real nonlinear = nodal_volume * nonlinear_source;

                    // K_u_nu: -dt/2 int div(u N_a) div(u q) dV.
                    Real grad_na[4] =
                    {
                        0.0,
                        grad(s, ie, 1, a),
                        grad(s, ie, 2, a),
                        grad(s, ie, 3, a)
                    };
                    Real div_u_na[5] = {};
                    liu_nithiarasu::test_function_divergence_nodes(
                        velocity,
                        grad_na,
                        div_u,
                        a,
                        div_u_na);
                    const Real characteristic = -0.5 * dt *
                        liu_nithiarasu::p1_product_integral(
                            volume,
                            div_u_na,
                            div_uq);

                    // Appendix B Eqs. B.56-B.57.  The closure coefficients are
                    // element values; the remaining q_h is integrated with the
                    // consistent P1 mass coefficients, not V/4 equal forcing.
                    const Real q_load =
                        liu_nithiarasu::consistent_scalar_load(
                            volume,
                            q,
                            a);
                    const Real production_rhs =
                        (q_average > 0.0 ? production / q_average : 0.0)
                        * q_load;
                    const Real destruction_rhs =
                        destruction_coefficient * q_average * q_load;

                    local_rhs[a] =
                        advection
                        + characteristic
                        + diffusion
                        + nonlinear
                        + production_rhs
                        - destruction_rhs;

                    local_production[a] = production_rhs;
                    local_destruction[a] = destruction_rhs;
                    local_diffusion[a] = diffusion + nonlinear;
                    local_source[a] = production_rhs - destruction_rhs;
                }

                // The generic thesis vector f_nuGamma permits a prescribed SA
                // diffusive flux.  CBS3D's current SA boundary set is instead:
                //   wall/inlet: strong Dirichlet;
                //   outlet:     homogeneous natural flux.
                // Therefore f_nuGamma is exactly zero for the supported cases.

                for (Int a = 1; a <= 4; ++a)
                {
                    const Int ip = s.intma(a, ie);
                    s.sa_rhs(ip) += local_rhs[a];
                    s.sa_production(ip) += local_production[a];
                    s.sa_destruction(ip) += local_destruction[a];
                    s.sa_diffusion(ip) += local_diffusion[a];
                    s.sa_source(ip) += local_source[a];
                }
            }
        }

        if (bad_detj != 0)
        {
            throw std::runtime_error(
                "SpalartAllmarasAssembly - invalid TET4 Jacobian");
        }
        if (bad_material != 0)
        {
            throw std::runtime_error(
                "SpalartAllmarasAssembly - invalid molecular viscosity/density");
        }
        if (bad_wall_distance != 0)
        {
            throw std::runtime_error(
                "SpalartAllmarasAssembly - invalid SA wall distance");
        }
        if (bad_dt != 0)
        {
            throw std::runtime_error(
                "SpalartAllmarasAssembly - invalid CBS element timestep");
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
        const Int count = have_owned
            ? static_cast<Int>(s.owned_nodes.size())
            : s.cfg.npoin;

#ifdef CBS3D_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (Int index = 0; index < count; ++index)
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
            Real new_value = old_value + s.elcoe2(ip) * s.sa_rhs(ip);

            if (!std::isfinite(new_value))
            {
#ifdef CBS3D_USE_OPENMP
#pragma omp critical(sa_liu_update_failure)
#endif
                {
                    if (failed_node == 0)
                    {
                        failed_node = ip;
                        failure_reason = "non-finite Liu SA explicit update";
                    }
                }
                continue;
            }

            // Liu's Step 4 is explicit.  sa_implicit_destruction is deliberately
            // not used by this production assembly.  The non-negative projection
            // remains a guard for the present non-SA-neg branch and must be
            // monitored during verification rather than mistaken for physics.
            new_value = std::max(s.cfg.sa_nu_tilde_floor, new_value);
            s.nu_tilde(ip) = new_value;
            s.sa_residual(ip) = new_value - old_value;
        }

        if (failed_node != 0)
        {
            throw std::runtime_error(
                describe_node_failure(s, failed_node, failure_reason));
        }
    }
}
