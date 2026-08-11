//=============================================================================
// CBS3D++_SI
//
// Energy-equation residual assembly for three-dimensional conjugate heat
// transfer.
//
// Governing equation:
//
//     rho cp [ dT/dt + u . grad(T) ]
//         = div(k grad(T)) + Q
//
// The assembled thermal residual is:
//
//     r_T = r_conv + r_stab + r_diff + r_source + r_flux
//
// where:
//
//     r_conv    fluid convection
//     r_stab    CBS/SUPG-style thermal stabilisation
//     r_diff    thermal diffusion in fluid and solid
//     r_source  volumetric heat generation
//     r_flux    prescribed surface heat flux
//
// The nodal update is performed later in Steps::step4Energy():
//
//     T^(n+1) = T^n + elcoe2p * rhs1
//
// The temperature stored in temperature1 is used during assembly, while the
// corrected velocity available after CBS Step 3 is used for fluid convection.
//=============================================================================

#include "cbs/assembly/EnergyAssembly.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace cbs
{
    namespace
    {
        // Returns the one-dimensional storage position of:
        //
        //     dN_local_node / dx_dim
        //
        // for tetrahedral element ie.
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


        // Fluid elements are identified by:
        //
        //     mat_elem(e) = 0
        bool is_fluid_element(
            const CBSStateSI& s,
            Int ie)
        {
            return s.mat_elem(ie) == 0;
        }


        // Checks the fixed dimensions required by the present three-dimensional
        // P1 tetrahedral energy formulation.
        void validate_energy_dimensions(const CBSStateSI& s)
        {
            if (s.cfg.ndim != 3 ||
                s.cfg.nep != 4 ||
                s.cfg.nsid != 4 ||
                s.cfg.nsidp != 3 ||
                s.cfg.ndim1 != 4)
            {
                throw std::runtime_error(
                    "EnergyAssembly - CBS3D energy assembly requires ndim=3, nep=4, nsid=4, nsidp=3, ndim1=4");
            }
        }


        // Returns one Cartesian derivative of a tetrahedral shape function:
        //
        //     grad(N_a)_dim = dN_a / dx_dim
        Real grad(
            const CBSStateSI& s,
            Int ie,
            Int dim,
            Int local_node)
        {
            return s.dNkdx(dNkdx_index(s, ie, dim, local_node));
        }


        //=====================================================================
        // Calculates the constant temperature gradient inside one P1
        // tetrahedral element.
        //
        // The finite-element interpolation is:
        //
        //     T(x,y,z) = sum_a N_a(x,y,z) T_a
        //
        // Therefore:
        //
        //     grad(T) = sum_a T_a grad(N_a)
        //
        // Since grad(N_a) is constant in a linear tetrahedron, grad(T) is also
        // constant inside the element.
        //
        // The temperature from the beginning of the CBS iteration,
        // temperature1, is used so that the thermal residual remains explicit
        // in temperature.
        //=====================================================================
        void compute_temperature_gradient(
            const CBSStateSI& s,
            Int ie,
            Real& dTdx,
            Real& dTdy,
            Real& dTdz)
        {
            dTdx = 0.0;
            dTdy = 0.0;
            dTdz = 0.0;

            // Step 4 uses the temperature saved at the start of the CBS
            // iteration.  This keeps the scalar update explicit in temperature,
            // as in CBS2D++_SI.
            for (Int a = 1; a <= s.cfg.nep; ++a)
            {
                const Int ip = s.intma(a, ie);
                const Real T = s.temperature1(ip);

                dTdx += grad(s, ie, 1, a) * T;
                dTdy += grad(s, ie, 2, a) * T;
                dTdz += grad(s, ie, 3, a) * T;
            }
        }


        // Checks the element thermal properties required by the energy
        // equation:
        //
        //     rho_cp_e = rho_e cp_e > 0
        //
        //     k_e > 0
        //
        //     Qvol_e is finite
        void validate_material_properties(
            const CBSStateSI& s,
            Int ie)
        {
            if (s.rho_cp_e(ie) <= 0.0 || !std::isfinite(s.rho_cp_e(ie)))
            {
                throw std::runtime_error(
                    "EnergyAssembly::assembleStep4Rhs - invalid rho*Cp at element "
                    + std::to_string(ie));
            }

            if (s.k_e(ie) <= 0.0 || !std::isfinite(s.k_e(ie)))
            {
                throw std::runtime_error(
                    "EnergyAssembly::assembleStep4Rhs - invalid thermal conductivity at element "
                    + std::to_string(ie));
            }

            if (!std::isfinite(s.Qvol_e(ie)))
            {
                throw std::runtime_error(
                    "EnergyAssembly::assembleStep4Rhs - invalid volumetric source at element "
                    + std::to_string(ie));
            }
        }


        //=====================================================================
        // Adds the Galerkin convection contribution for one fluid element.
        //
        // Strong thermal-advection term:
        //
        //     rho cp u . grad(T)
        //
        // Weak residual contribution:
        //
        //     r_conv,a^(e)
        //       = -integral(V_e) N_a rho cp
        //          [u . grad(T)] dV
        //
        // For P1 tetrahedra:
        //
        //     integral(V_e) N_a N_b dV
        //       = V_e/10    when a = b
        //       = V_e/20    when a != b
        //
        // Hence:
        //
        //     r_conv,a^(e)
        //       = -rho cp V_e/20
        //          [(sum_b u_b + u_a) . grad(T)]
        //
        // Since:
        //
        //     det(J_e) = 6 V_e
        //
        // the implemented factor is:
        //
        //     rho cp det(J_e) / 120
        //
        // The corrected velocity from CBS Step 3 is used.
        //=====================================================================
        void add_fluid_convection(
            const CBSStateSI& s,
            Int ie,
            Real dTdx,
            Real dTdy,
            Real dTdz,
            Real lrhs[5])
        {
            // Weak advection term:
            //
            //   - int_Omega N_a rhoCp (u.gradT) dOmega
            //
            // For P1 tetrahedra, u is linearly interpolated and gradT is
            // constant.  The exact mass-integrated nodal factor is:
            //
            //   int N_a N_b dOmega =
            //       V/10 if a=b
            //       V/20 if a!=b
            //
            // Therefore:
            //
            //   R_a_conv = -rhoCp * V/20 * (sum_b u_b + u_a).gradT
            //
            // With detJ = 6V, V/20 = detJ/120.
            const Real adv_factor = s.rho_cp_e(ie) * s.detJ(ie) / 120.0;

            Real u_sum = 0.0;
            Real v_sum = 0.0;
            Real w_sum = 0.0;

            for (Int b = 1; b <= s.cfg.nep; ++b)
            {
                const Int ip = s.intma(b, ie);

                // Use the corrected velocity available after Step 3.
                u_sum += s.unkno(1, ip);
                v_sum += s.unkno(2, ip);
                w_sum += s.unkno(3, ip);
            }

            for (Int a = 1; a <= s.cfg.nep; ++a)
            {
                const Int ip = s.intma(a, ie);

                const Real u_weight = u_sum + s.unkno(1, ip);
                const Real v_weight = v_sum + s.unkno(2, ip);
                const Real w_weight = w_sum + s.unkno(3, ip);

                const Real advT =
                    u_weight * dTdx +
                    v_weight * dTdy +
                    w_weight * dTdz;

                lrhs[a] -= adv_factor * advT;
            }
        }


        //=====================================================================
        // Adds the CBS convection-stabilisation contribution for one
        // fluid element.
        //
        // The implemented scalar characteristic term is:
        //
        //     r_stab,a^(e)
        //       = dt_e/2 integral(V_e)
        //          rho cp [u_bar . grad(N_a)]
        //                 [u_bar . grad(T)] dV
        //
        // where the element-average velocity is:
        //
        //     u_bar = (1/4) sum_a u_a
        //
        // Because the velocity averages and the P1 gradients are constant
        // within the element:
        //
        //     r_stab,a^(e)
        //       = dt_e/2 rho cp V_e
        //          [u_bar . grad(N_a)]
        //          [u_bar . grad(T)]
        //
        // This term is applied only to fluid elements.
        //=====================================================================
        void add_fluid_convection_stabilisation(
            const CBSStateSI& s,
            Int ie,
            Real dTdx,
            Real dTdy,
            Real dTdz,
            Real lrhs[5])
        {
            // Scalar CBS/SUPG-style characteristic correction for thermal
            // advection:
            //
            //   + dt/2 int_Omega rhoCp (u.gradN_a) (u.gradT) dOmega
            //
            // This is the scalar analogue of the Step-1 characteristic
            // correction used for momentum.  It is applied only in the fluid,
            // where thermal advection exists.  The element-average velocity is
            // used to keep the first CHT port robust and deterministic.
            const Real dt = s.delte(ie);

            if (dt <= 0.0 || !std::isfinite(dt))
            {
                throw std::runtime_error(
                    "EnergyAssembly::assembleStep4Rhs - element timestep must be positive for thermal stabilisation at element "
                    + std::to_string(ie));
            }

            Real ubar = 0.0;
            Real vbar = 0.0;
            Real wbar = 0.0;

            for (Int a = 1; a <= s.cfg.nep; ++a)
            {
                const Int ip = s.intma(a, ie);

                ubar += s.unkno(1, ip);
                vbar += s.unkno(2, ip);
                wbar += s.unkno(3, ip);
            }

            ubar /= static_cast<Real>(s.cfg.nep);
            vbar /= static_cast<Real>(s.cfg.nep);
            wbar /= static_cast<Real>(s.cfg.nep);

            const Real advT =
                ubar * dTdx +
                vbar * dTdy +
                wbar * dTdz;

            const Real volume = s.detJ(ie) / 6.0;
            const Real tau_factor = 0.5 * dt * s.rho_cp_e(ie) * volume;

            for (Int a = 1; a <= s.cfg.nep; ++a)
            {
                const Real u_gradNa =
                    ubar * grad(s, ie, 1, a) +
                    vbar * grad(s, ie, 2, a) +
                    wbar * grad(s, ie, 3, a);

                lrhs[a] += tau_factor * u_gradNa * advT;
            }
        }


        //=====================================================================
        // Adds the thermal-diffusion contribution for one element.
        //
        // Conductive term:
        //
        //     div(k grad(T))
        //
        // After integration by parts:
        //
        //     r_diff,a^(e)
        //       = -integral(V_e)
        //          k grad(N_a) . grad(T) dV
        //
        //       + integral(Gamma_e)
        //          N_a k grad(T) . n dGamma
        //
        // The volume contribution is assembled here:
        //
        //     r_diff,a^(e)
        //       = -k_e V_e grad(N_a) . grad(T)
        //
        // Prescribed external heat flux is added separately by
        // add_prescribed_heat_flux().
        //
        // This term is assembled in both fluid and solid elements.
        //=====================================================================
        void add_thermal_diffusion(
            const CBSStateSI& s,
            Int ie,
            Real dTdx,
            Real dTdy,
            Real dTdz,
            Real lrhs[5])
        {
            // Diffusion weak form after integration by parts:
            //
            //   - int_Omega k grad(N_a).grad(T) dOmega
            //
            // External Neumann heat flux is not included here.  It is added
            // explicitly only for BC 532 in add_prescribed_heat_flux().
            const Real volume = s.detJ(ie) / 6.0;

            // Laminar and solid calculation:
            //
            //     k_used = k_e
            //
            // Turbulent fluid heat-transfer calculation:
            //
            //     k_used = k_eff_e
            //
            // where:
            //
            //     k_eff_e = k_e + rho cp nu_t / Pr_t
            //
            // The turbulent addition is never applied in solid elements.
            Real k = s.k_e(ie);

            if (s.cfg.turbulence_on > 0 &&
                s.cfg.turbulent_thermal_diffusivity_on > 0 &&
                is_fluid_element(s, ie))
            {
                k = s.k_eff_e(ie);
            }

            if (k <= 0.0 || !std::isfinite(k))
            {
                throw std::runtime_error(
                    "EnergyAssembly::assembleStep4Rhs - invalid effective thermal conductivity at element "
                    + std::to_string(ie));
            }

            for (Int a = 1; a <= s.cfg.nep; ++a)
            {
                const Real gradNi_dot_gradT =
                    grad(s, ie, 1, a) * dTdx +
                    grad(s, ie, 2, a) * dTdy +
                    grad(s, ie, 3, a) * dTdz;

                lrhs[a] -= k * volume * gradNi_dot_gradT;
            }
        }


        //=====================================================================
        // Resolves the element volumetric heat source.
        //
        // Input contract:
        //
        //     fluid element (mat_elem == 0)
        //         Q_e = Qvol_e from .matprop
        //
        //     solid element (mat_elem != 0)
        //         source_solid != 0  -> Q_e = source_solid from .par
        //         source_solid == 0  -> Q_e = Qvol_e from .matprop
        //
        // A non-zero .par source and non-zero solid .matprop source are not
        // added together.  That ambiguous double specification is rejected.
        //=====================================================================
        Real element_volumetric_source(
            const CBSStateSI& s,
            Int ie)
        {
            const Real material_source = s.Qvol_e(ie);

            if (!std::isfinite(material_source))
            {
                throw std::runtime_error(
                    "EnergyAssembly::assembleStep4Rhs - non-finite .matprop volumetric source at element "
                    + std::to_string(ie));
            }

            if (is_fluid_element(s, ie))
            {
                return material_source;
            }

            const Real parameter_source = s.cfg.source_solid;

            if (!std::isfinite(parameter_source))
            {
                throw std::runtime_error(
                    "EnergyAssembly::assembleStep4Rhs - non-finite source_solid in .par");
            }

            if (parameter_source != 0.0 && material_source != 0.0)
            {
                throw std::runtime_error(
                    "EnergyAssembly::assembleStep4Rhs - solid volumetric heat source is specified in both .par source_solid and .matprop Qvol; use exactly one source definition");
            }

            return parameter_source != 0.0
                ? parameter_source
                : material_source;
        }


        //=====================================================================
        // Adds a uniform element volumetric heat source.
        //
        // Weak source contribution:
        //
        //     r_source,a^(e)
        //       = integral(V_e) N_a Q_e dV
        //
        // For a P1 tetrahedron:
        //
        //     integral(V_e) N_a dV = V_e/4
        //
        // therefore:
        //
        //     r_source,a^(e) = Q_e V_e/4
        //
        // Since det(J_e) = 6V_e:
        //
        //     V_e/4 = det(J_e)/24
        //
        // which is represented by mass_factor.
        //=====================================================================
        void add_volumetric_source(
            const CBSStateSI& s,
            Int ie,
            Real lrhs[5])
        {
            // Dimensional source term:
            //
            //   int_Omega N_a Q dOmega = Q * V/4 = Q * detJ/24.
            //
            // Solid source_solid is read from .par. Per-material Qvol_e from
            // .matprop remains available when source_solid is zero and is the
            // only volumetric-source mechanism for fluid elements.
            const Real qvol = element_volumetric_source(s, ie);

            if (qvol == 0.0)
            {
                return;
            }

            const Real source = qvol * s.detJ(ie) * s.cfg.mass_factor;

            for (Int a = 1; a <= s.cfg.nep; ++a)
            {
                lrhs[a] += source;
            }
        }


        //=====================================================================
        // Adds the prescribed external heat-flux contribution on BC 532.
        //
        // The adopted sign convention is:
        //
        //     positive heat_flux_bc = heat entering the thermal domain
        //
        // For a triangular P1 boundary face:
        //
        //     r_flux,a^(f)
        //       = integral(Gamma_f) N_a q'' dGamma
        //
        // Since:
        //
        //     integral(Gamma_f) N_a dGamma = A_f/3
        //
        // each face node receives:
        //
        //     r_flux,a^(f) = q'' A_f/3
        //
        // No heat-flux term is applied on BC 901. That boundary is the
        // conformal fluid-solid interface.
        //=====================================================================
        void add_prescribed_heat_flux(CBSStateSI& s)
        {
            // Prescribed heat flux is applied only on BC 532.
            //
            // For a triangular P1 face:
            //
            //   int_Gamma N_a q'' dGamma = q'' * A/3.
            //
            // Positive heat_flux_bc means heat entering the computational
            // thermal domain.
            if (s.cfg.heat_flux_bc == 0.0)
            {
                return;
            }

            for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
            {
                const Int bc = s.iside(s.cfg.bsid, ib);

                if (bc != s.cfg.bc_noslip_heatflux_wall)
                {
                    continue;
                }

                const Real area = s.face_norm(4, ib);

                if (area <= 0.0 || !std::isfinite(area))
                {
                    throw std::runtime_error(
                        "EnergyAssembly::assembleStep4Rhs - invalid heat-flux boundary face area");
                }

                const Real contribution =
                    s.cfg.heat_flux_bc * area / static_cast<Real>(s.cfg.nsidp);

                for (Int in = 1; in <= s.cfg.nsidp; ++in)
                {
                    const Int ip = s.iside(in, ib);

                    if (ip < 1 || ip > s.cfg.npoin)
                    {
                        throw std::runtime_error(
                            "EnergyAssembly::assembleStep4Rhs - heat-flux boundary node out of range");
                    }

                    s.rhs1(ip) += contribution;
                }
            }
        }
    }


    //=========================================================================
    // Assembles the complete CBS Step 4 thermal residual.
    //
    // For every tetrahedral element:
    //
    //     1. Validate geometry and material properties.
    //     2. Calculate grad(T) from temperature1.
    //     3. Add convection and stabilisation for fluid elements.
    //     4. Add diffusion for fluid and solid elements.
    //     5. Add volumetric heat generation.
    //     6. Scatter the element residual into rhs1.
    //
    // After all volume terms are assembled, prescribed heat flux is added on
    // BC 532 boundary faces.
    //
    // Fluid-solid interface continuity requires no separate boundary term
    // because the conformal mesh shares the same nodal temperature unknowns
    // across BC 901.
    //
    // Inputs:
    //     temperature1  temperature at the beginning of the CBS iteration
    //     unkno         corrected velocity from CBS Step 3
    //     rho_cp_e      element volumetric heat capacity
    //     k_e           element thermal conductivity
    //     Qvol_e        optional per-material volumetric heat source
    //     source_solid  optional uniform solid source from .par
    //     delte         element time step used by stabilisation
    //
    // Output:
    //     rhs1          assembled global thermal residual
    //=========================================================================
    void EnergyAssembly::assembleStep4Rhs(CBSStateSI& s)
    {
        validate_energy_dimensions(s);

        s.rhs1.fill(0.0);

        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            if (s.detJ(ie) <= 0.0 || !std::isfinite(s.detJ(ie)))
            {
                throw std::runtime_error(
                    "EnergyAssembly::assembleStep4Rhs - invalid detJ at element "
                    + std::to_string(ie));
            }

            validate_material_properties(s, ie);

            Real dTdx = 0.0;
            Real dTdy = 0.0;
            Real dTdz = 0.0;

            compute_temperature_gradient(s, ie, dTdx, dTdy, dTdz);

            Real lrhs[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };

            if (is_fluid_element(s, ie))
            {
                add_fluid_convection(s, ie, dTdx, dTdy, dTdz, lrhs);
                add_fluid_convection_stabilisation(s, ie, dTdx, dTdy, dTdz, lrhs);
            }

            // Diffusion and volumetric source are assembled in both fluid and
            // solid materials. Interface continuity is handled naturally by
            // the conformal shared-node FEM mesh, so BC 901 has no boundary RHS.
            // In MPI production, rhs1 is reverse-added from ghosts to owners
            // immediately after this shared assembly, which also parallelises
            // the volumetric-source contribution without a separate collective.
            add_thermal_diffusion(s, ie, dTdx, dTdy, dTdz, lrhs);
            add_volumetric_source(s, ie, lrhs);

            for (Int a = 1; a <= s.cfg.nep; ++a)
            {
                const Int ip = s.intma(a, ie);
                s.rhs1(ip) += lrhs[a];
            }
        }

        add_prescribed_heat_flux(s);
    }


    //=========================================================================
    // Reserved hook for real-time or BDF thermal-history contributions.
    //
    // The present implementation intentionally performs no operation. The
    // current Step 4 update is:
    //
    //     T^(n+1) = T^n + elcoe2p * rhs1
    //
    // using the lumped thermal capacitance assembled in
    // Preprocess::massMatrix().
    //=========================================================================
    void EnergyAssembly::applyRealTimeEnergyTerm(CBSStateSI& s)
    {
        // Real-time/BDF scalar history terms are intentionally not included in
        // the first 3D CHT assembly.  The current semi-implicit driver updates:
        //
        //   temperature = temperature1 + rhs1 * elcoe2p
        //
        // using the thermal capacitance assembled in Preprocess::massMatrix().
        (void)s;
    }
}
