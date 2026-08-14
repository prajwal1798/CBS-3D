#include "cbs/assembly/SpalartAllmarasAssembly.hpp"

#include "cbs/turbulence/SpalartAllmaras.hpp"

#include "cbs/parallel/HaloExchange.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace cbs
{
    namespace
    {
        bool is_fluid_element(const CBSStateSI& s, Int ie)
        {
            return s.mat_elem(ie) == 0;
        }

        Int dNkdx_index(
            const CBSStateSI& s,
            Int ie,
            Int dim,
            Int local_node)
        {
            return (ie - 1) * s.cfg.ndim * s.cfg.nep
                + (dim - 1) * s.cfg.nep
                + local_node;
        }

        Real grad(
            const CBSStateSI& s,
            Int ie,
            Int dim,
            Int local_node)
        {
            return s.dNkdx(dNkdx_index(s, ie, dim, local_node));
        }

        Real positive_part(Real value)
        {
            return value > 0.0 ? value : 0.0;
        }

        // Returns the molecular kinematic viscosity of an element.
        //
        // This routine is called from inside an OpenMP parallel region, so it
        // must not throw: an exception that leaves a structured parallel block
        // terminates the process instead of unwinding, and under MPI it would
        // abort one rank while the others waited in the next collective.  Invalid
        // input is reported through ok and converted into a single exception on
        // the master thread after the parallel region has closed.
        Real molecular_kinematic_viscosity(
            const CBSStateSI& s,
            Int ie,
            bool& ok)
        {
            ok = true;

            if (s.cfg.dimensional_mode > 0 && s.cfg.material_properties_enabled > 0)
            {
                if (s.rho_e(ie) <= 0.0 || !std::isfinite(s.rho_e(ie)))
                {
                    ok = false;
                    return 0.0;
                }

                if (s.mu_e(ie) < 0.0 || !std::isfinite(s.mu_e(ie)))
                {
                    ok = false;
                    return 0.0;
                }

                return s.mu_e(ie) / s.rho_e(ie);
            }

            if (s.cfg.ani <= 0.0 || !std::isfinite(s.cfg.ani))
            {
                ok = false;
                return 0.0;
            }

            return s.cfg.ani;
        }

        // Element-averaged wall distance.  Non-throwing for the same reason as
        // molecular_kinematic_viscosity.
        Real element_wall_distance(
            const CBSStateSI& s,
            Int ie,
            bool& ok)
        {
            ok = true;

            Real d = 0.0;

            for (Int a = 1; a <= s.cfg.nep; ++a)
            {
                const Int ip = s.intma(a, ie);
                const Real node_distance = s.wall_distance(ip);

                if (node_distance <= 0.0 || !std::isfinite(node_distance))
                {
                    ok = false;
                    return s.cfg.sa_min_wall_distance;
                }

                d += node_distance;
            }

            d /= static_cast<Real>(s.cfg.nep);
            d = std::max(d, s.cfg.sa_min_wall_distance);

            return d;
        }

        void validate_transport_dimensions(const CBSStateSI& s)
        {
            if (s.cfg.ndim != 3 ||
                s.cfg.nep != 4 ||
                s.cfg.nsid != 4 ||
                s.cfg.nsidp != 3)
            {
                throw std::runtime_error(
                    "SpalartAllmarasAssembly - SA transport requires ndim=3, nep=4, nsid=4, nsidp=3");
            }
        }

        Real element_vorticity_magnitude(
            const CBSStateSI& s,
            Int ie)
        {
            Real duidxj[4][4] = {};

            for (Int i = 1; i <= s.cfg.ndim; ++i)
            {
                for (Int a = 1; a <= s.cfg.nep; ++a)
                {
                    const Int ip = s.intma(a, ie);

                    for (Int j = 1; j <= s.cfg.ndim; ++j)
                    {
                        duidxj[i][j] += s.unkno(i, ip) * grad(s, ie, j, a);
                    }
                }
            }

            const Real omega_x = duidxj[3][2] - duidxj[2][3];
            const Real omega_y = duidxj[1][3] - duidxj[3][1];
            const Real omega_z = duidxj[2][1] - duidxj[1][2];

            return std::sqrt(
                omega_x * omega_x +
                omega_y * omega_y +
                omega_z * omega_z);
        }

        void compute_sa_gradient(
            const CBSStateSI& s,
            Int ie,
            const Real q_node[5],
            Real gradient_q[4])
        {
            gradient_q[1] = 0.0;
            gradient_q[2] = 0.0;
            gradient_q[3] = 0.0;

            for (Int a = 1; a <= s.cfg.nep; ++a)
            {
                for (Int dim = 1; dim <= s.cfg.ndim; ++dim)
                {
                    gradient_q[dim] += q_node[a] * grad(s, ie, dim, a);
                }
            }
        }

        void compute_velocity_sum(
            const CBSStateSI& s,
            Int ie,
            Real velocity_sum[4])
        {
            velocity_sum[1] = 0.0;
            velocity_sum[2] = 0.0;
            velocity_sum[3] = 0.0;

            for (Int a = 1; a <= s.cfg.nep; ++a)
            {
                const Int ip = s.intma(a, ie);

                velocity_sum[1] += s.unkno(1, ip);
                velocity_sum[2] += s.unkno(2, ip);
                velocity_sum[3] += s.unkno(3, ip);
            }
        }
    }

    //=========================================================================
    // Resets all turbulence-derived properties to their laminar molecular values.
    //=========================================================================
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

    //=========================================================================
    // Converts the transported SA variable into eddy viscosity and effective
    // material properties.
    //=========================================================================
    void SpalartAllmarasAssembly::updateEddyViscosity(CBSStateSI& s)
    {
        resetEffectiveProperties(s);

        if (s.cfg.turbulence_on < 1)
        {
            return;
        }

        turbulence::SpalartAllmarasConstants constants;

        // The nodal averages are accumulated directly into state arrays rather
        // than local vectors, because under domain decomposition the sums and
        // the counts both have to cross partition interfaces before the division
        // is performed.  Averaging local contributions first and exchanging
        // afterwards would make the interface value depend on how the mesh was
        // partitioned.
        s.nu_t.fill(0.0);
        s.mu_t.fill(0.0);
        s.sa_nodal_weight.fill(0.0);

        bool bad_material = false;

        // The element loop scatters into nodal accumulators, so it is driven by
        // the same colour arrays used by every other assembly in the solver.
        // Within one colour no two elements share a node, so the scatter needs
        // neither atomics nor per-thread buffers.
        for (Int colour = 0; colour < s.ncolor; ++colour)
        {
            const Int cbeg = s.color_ptr[static_cast<std::size_t>(colour)];
            const Int cend = s.color_ptr[static_cast<std::size_t>(colour) + 1];

#pragma omp parallel for schedule(static)
            for (Int k = cbeg; k < cend; ++k)
            {
                const Int ie = s.color_elem[static_cast<std::size_t>(k)];

                if (!is_fluid_element(s, ie))
                {
                    continue;
                }

                Real nu_tilde_avg = 0.0;

                for (Int in = 1; in <= s.cfg.nep; ++in)
                {
                    const Int ip = s.intma(in, ie);

                    // A node on a no-slip wall carries nu_tilde = 0.  Using the
                    // floor there instead would leave a non-zero working
                    // variable on the wall and produce a spurious eddy viscosity
                    // in the first cell, which is exactly where the log law is
                    // most sensitive.  This matches the treatment in
                    // assembleTransportRhs.
                    const Real nodal_value = s.sa_wall_node(ip) != 0
                        ? 0.0
                        : std::max(s.cfg.sa_nu_tilde_floor, s.nu_tilde(ip));

                    nu_tilde_avg += nodal_value;
                }

                nu_tilde_avg /= static_cast<Real>(s.cfg.nep);

                bool material_ok = true;
                const Real molecular_nu =
                    molecular_kinematic_viscosity(s, ie, material_ok);

                if (!material_ok)
                {
#pragma omp atomic write
                    bad_material = true;
                    continue;
                }

                Real nu_t = 0.0;
                if (molecular_nu > 0.0)
                {
                    nu_t = turbulence::eddyKinematicViscosity(
                        nu_tilde_avg,
                        molecular_nu,
                        constants);
                }

                const Real mu_t = s.rho_e(ie) * nu_t;

                s.nu_tilde_e(ie) = nu_tilde_avg;
                s.nu_t_e(ie) = nu_t;
                s.mu_t_e(ie) = mu_t;
                s.mu_eff_e(ie) = s.mu_e(ie) + mu_t;

                if (s.cfg.turbulent_thermal_diffusivity_on > 0)
                {
                    s.k_eff_e(ie) = s.k_e(ie)
                        + s.rho_cp_e(ie) * nu_t / s.cfg.sa_prandtl_t;
                }

                for (Int in = 1; in <= s.cfg.nep; ++in)
                {
                    const Int ip = s.intma(in, ie);

                    s.nu_t(ip) += nu_t;
                    s.mu_t(ip) += mu_t;
                    s.sa_nodal_weight(ip) += 1.0;
                }
            }
        }

        if (bad_material)
        {
            throw std::runtime_error(
                "SpalartAllmarasAssembly::updateEddyViscosity - invalid density,"
                " molecular viscosity or non-dimensional ani at one or more elements");
        }

        // Complete the sums and the counts on the owning rank before dividing.
#ifdef CBS3D_USE_MPI
        if (s.mpi_enabled && s.mpi_size > 1)
        {
            HaloExchange::sumGhostContributionsToOwners(
                s.nu_t,
                s.partition_metadata,
                MPI_COMM_WORLD);

            HaloExchange::sumGhostContributionsToOwners(
                s.mu_t,
                s.partition_metadata,
                MPI_COMM_WORLD);

            HaloExchange::sumGhostContributionsToOwners(
                s.sa_nodal_weight,
                s.partition_metadata,
                MPI_COMM_WORLD);
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

        // Push the completed averages back onto the ghost copies so that every
        // rank sees the same nodal eddy-viscosity field.
#ifdef CBS3D_USE_MPI
        if (s.mpi_enabled && s.mpi_size > 1)
        {
            HaloExchange::broadcastOwnedToGhosts(
                s.nu_t,
                s.partition_metadata,
                MPI_COMM_WORLD);

            HaloExchange::broadcastOwnedToGhosts(
                s.mu_t,
                s.partition_metadata,
                MPI_COMM_WORLD);
        }
#endif
    }

    //=========================================================================
    // Assembles the finite-element residual of the SA transport equation.
    //
    // The transported scalar is
    //
    //     q = nu_tilde
    //
    // and the implemented strong form is
    //
    //     dq/dt + u . grad(q)
    //       = P - D
    //       + 1/sigma div[(nu + q) grad(q)]
    //       + cb2/sigma grad(q) . grad(q)
    //
    // The weak residual assembled here is used in
    //
    //     M_L/dt (q_new - q_old) = R_SA
    //
    // with q_old stored in nu_tilde1.  The corrected velocity from CBS Step 3,
    // stored in unkno, is used for advection and vorticity production.
    //
    // Boundary conditions are not imposed by adding boundary flux terms here.
    // Wall and inlet values are imposed strongly after the nodal update.  The
    // zero-gradient outlet condition is the natural boundary condition of the
    // diffusion term and therefore needs no explicit RHS contribution.
    //=========================================================================
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

        // Error flags raised inside the OpenMP region and acted on afterwards.
        // Every write stores the same value, true, so the only requirement is
        // that the store is not torn; the atomic write makes that explicit
        // rather than relying on it.
        bool bad_detj = false;
        bool bad_material = false;
        bool bad_wall_distance = false;

        for (Int colour = 0; colour < s.ncolor; ++colour)
        {
            const Int cbeg = s.color_ptr[static_cast<std::size_t>(colour)];
            const Int cend = s.color_ptr[static_cast<std::size_t>(colour) + 1];

#pragma omp parallel for schedule(static)
            for (Int k = cbeg; k < cend; ++k)
            {
                const Int ie = s.color_elem[static_cast<std::size_t>(k)];

                if (!is_fluid_element(s, ie))
                {
                    continue;
                }

                if (s.detJ(ie) <= 0.0 || !std::isfinite(s.detJ(ie)))
                {
#pragma omp atomic write
                    bad_detj = true;
                    continue;
                }

                const Real volume = s.detJ(ie) / 6.0;
                const Real nodal_volume = volume / static_cast<Real>(s.cfg.nep);
                const Real advective_mass_factor = volume / 20.0;

                Real q_node[5] = {};
                Real q_sum = 0.0;

                for (Int a = 1; a <= s.cfg.nep; ++a)
                {
                    const Int ip = s.intma(a, ie);

                    q_node[a] = std::max(
                        s.cfg.sa_nu_tilde_floor,
                        s.nu_tilde1(ip));

                    if (s.sa_wall_node(ip) != 0)
                    {
                        q_node[a] = 0.0;
                    }

                    q_sum += q_node[a];
                }

                const Real q_average = q_sum / static_cast<Real>(s.cfg.nep);

                bool material_ok = true;
                const Real molecular_nu =
                    molecular_kinematic_viscosity(s, ie, material_ok);

                if (!material_ok)
                {
#pragma omp atomic write
                    bad_material = true;
                    continue;
                }

                bool distance_ok = true;
                const Real wall_distance =
                    element_wall_distance(s, ie, distance_ok);

                if (!distance_ok)
                {
#pragma omp atomic write
                    bad_wall_distance = true;
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

                const Real omega = element_vorticity_magnitude(s, ie);

                Real s_bar = 0.0;
                Real s_tilde = 0.0;

                if (q_average > 0.0)
                {
                    s_bar = turbulence::sBar(
                        q_average,
                        molecular_nu,
                        wall_distance,
                        constants);
                }

                if (s.cfg.sa_use_stilde_limiter > 0)
                {
                    s_tilde = turbulence::limitedSTilde(omega, s_bar, constants);
                }
                else
                {
                    s_tilde = omega + s_bar;
                }

                s_tilde = std::max(s_tilde, s.cfg.sa_min_stilde);

                const Real chi_value = turbulence::chi(q_average, molecular_nu);
                const Real ft2_value = turbulence::ft2(chi_value, constants);
                const Real r_value = turbulence::rFunction(
                    q_average,
                    s_tilde,
                    wall_distance,
                    constants);
                const Real fw_value = turbulence::fw(r_value, constants);

                const Real production = turbulence::productionTerm(
                    q_average,
                    s_tilde,
                    ft2_value,
                    constants);

                const Real destruction_coefficient =
                    turbulence::destructionCoefficient(
                        wall_distance,
                        fw_value,
                        ft2_value,
                        constants);

                const Real destruction = destruction_coefficient * q_average * q_average;

                const Real diffusion_coefficient =
                    (molecular_nu + q_average) / constants.sigma;

                const Real nonlinear_gradient_source =
                    constants.cb2 * grad_q_squared / constants.sigma;

                Real local_rhs[5] = {};
                Real local_production[5] = {};
                Real local_destruction[5] = {};
                Real local_diffusion[5] = {};
                Real local_source[5] = {};

                // local_destruction_lhs stores the nodal coefficient added to
                // the left-hand side when positive destruction is treated
                // semi-implicitly:
                //
                //     D = C_D q^2  ->  C_D q_old q_new
                //
                // It is accumulated into sa_destruction_lhs, which has this
                // single meaning.  It must not be stored in sa_residual: the
                // distributed solver sums sa_rhs and sa_destruction_lhs across
                // partition interfaces before the nodal update, and sa_residual
                // holds the post-update increment used by the convergence
                // monitor.  Overloading one array for both breaks the exchange.
                Real local_destruction_lhs[5] = {};

                for (Int a = 1; a <= s.cfg.nep; ++a)
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

                    const Real production_rhs = nodal_volume * production;
                    const Real destruction_rhs = nodal_volume * destruction;
                    const Real nonlinear_rhs = nodal_volume * nonlinear_gradient_source;

                    local_rhs[a] += advection;
                    local_rhs[a] += diffusion;
                    local_rhs[a] += production_rhs;
                    local_rhs[a] += nonlinear_rhs;

                    // Destruction is stiff near the wall because it contains
                    // (nu_tilde/d)^2.  If sa_implicit_destruction is enabled and
                    // the destruction coefficient is positive, use the standard
                    // first-order linearisation:
                    //
                    //     C_D q^2  ->  C_D q_old q_new
                    //
                    // The term C_D q_old is therefore added to the nodal
                    // mass/time diagonal.  Negative C_D is not a destruction
                    // sink; it is kept explicit in the RHS with the correct sign.
                    if (s.cfg.sa_implicit_destruction > 0 &&
                        destruction_coefficient > 0.0 &&
                        q_average > 0.0)
                    {
                        local_destruction_lhs[a] +=
                            nodal_volume * destruction_coefficient * q_average;
                    }
                    else
                    {
                        local_rhs[a] -= destruction_rhs;
                    }

                    local_production[a] += production_rhs;
                    local_destruction[a] += destruction_rhs;
                    local_diffusion[a] += diffusion + nonlinear_rhs;
                    local_source[a] += production_rhs - destruction_rhs;
                }

                for (Int a = 1; a <= s.cfg.nep; ++a)
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

        if (bad_detj)
        {
            throw std::runtime_error(
                "SpalartAllmarasAssembly::assembleTransportRhs - invalid detJ at one or more elements");
        }

        if (bad_material)
        {
            throw std::runtime_error(
                "SpalartAllmarasAssembly::assembleTransportRhs - invalid density,"
                " molecular viscosity or non-dimensional ani at one or more elements");
        }

        if (bad_wall_distance)
        {
            throw std::runtime_error(
                "SpalartAllmarasAssembly::assembleTransportRhs - invalid wall"
                " distance at one or more element nodes;"
                " check that WallDistance::compute ran during preprocessing");
        }
    }

    //=========================================================================
    // Applies the nodal SA update.
    //
    // The current milestone uses the explicit lumped update:
    //
    //     q_new = q_old + [M_L/dt]^(-1) R_SA
    //
    // where:
    //
    //     q       = nu_tilde
    //     q_old   = nu_tilde1
    //     R_SA    = sa_rhs
    //     elcoe2  = [M_L/dt]^(-1) on fluid-connected nodes
    //
    // Negative values are clipped to sa_nu_tilde_floor in the first standard
    // non-negative SA branch.  The SA-neg branch will later replace this clipping
    // with the proper negative-model equations.
    //=========================================================================
    void SpalartAllmarasAssembly::updateNuTilde(CBSStateSI& s)
    {
        if (s.cfg.turbulence_on < 1)
        {
            return;
        }

#pragma omp parallel for schedule(static)
        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
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

            if (s.elcoe2(ip) <= 0.0 || !std::isfinite(s.elcoe2(ip)))
            {
                s.nu_tilde(ip) = s.nu_tilde1(ip);
                s.sa_residual(ip) = 0.0;
                continue;
            }

            const Real old_value = s.nu_tilde1(ip);

            Real denominator = 1.0;
            if (s.cfg.sa_implicit_destruction > 0)
            {
                denominator += s.elcoe2(ip) * s.sa_destruction_lhs(ip);
            }

            if (denominator <= 0.0 || !std::isfinite(denominator))
            {
                denominator = 1.0;
            }

            Real new_value =
                (old_value + s.elcoe2(ip) * s.sa_rhs(ip)) / denominator;

            if (!std::isfinite(new_value))
            {
                new_value = old_value;
            }

            new_value = std::max(s.cfg.sa_nu_tilde_floor, new_value);

            s.nu_tilde(ip) = new_value;
            s.sa_residual(ip) = new_value - old_value;
        }
    }
}
