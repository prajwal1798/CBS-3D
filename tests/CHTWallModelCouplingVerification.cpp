#include "cbs/turbulence/CHTWallModelCoupling.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>

namespace
{
    using cbs::CBSStateSI;
    using cbs::Int;
    using cbs::Real;
    using cbs::turbulence::CHTWallModelCoupling;
    using cbs::turbulence::ThermalWallTreatment;
    using cbs::turbulence::WallTreatment;

    bool close(
        const Real a,
        const Real b,
        const Real relative = 5.0e-11,
        const Real absolute = 5.0e-12)
    {
        return std::abs(a - b) <=
            absolute + relative * std::max(std::abs(a), std::abs(b));
    }

    void set_face_map(CBSStateSI& s)
    {
        s.ippn1.resize(4, 3);
        s.ippn1(1, 1) = 2; s.ippn1(1, 2) = 3; s.ippn1(1, 3) = 4;
        s.ippn1(2, 1) = 1; s.ippn1(2, 2) = 4; s.ippn1(2, 3) = 3;
        s.ippn1(3, 1) = 1; s.ippn1(3, 2) = 2; s.ippn1(3, 3) = 4;
        s.ippn1(4, 1) = 1; s.ippn1(4, 2) = 3; s.ippn1(4, 3) = 2;
    }

    Int grad_index(
        const CBSStateSI& s,
        const Int ie,
        const Int dim,
        const Int a)
    {
        return (ie - 1) * s.cfg.ndim * s.cfg.nep
            + (dim - 1) * s.cfg.nep + a;
    }

    CBSStateSI make_state()
    {
        CBSStateSI s;
        s.cfg.npoin = 5;
        s.cfg.nelem = 2;
        s.cfg.nboun = 0;
        s.cfg.temp_calc = 1;
        s.cfg.turbulence_on = 1;
        s.cfg.turbulence_model = 0;
        s.cfg.turbulent_thermal_diffusivity_on = 1;
        s.cfg.dimensional_mode = 1;
        s.cfg.material_properties_enabled = 1;
        s.cfg.sa_prandtl_t = 0.90;

        set_face_map(s);

        // Fluid tet above z=0; solid tet below z=0.  They share nodes 1,2,3.
        s.intma.resize(4, 2);
        s.intma(1, 1) = 1; s.intma(2, 1) = 2;
        s.intma(3, 1) = 3; s.intma(4, 1) = 4;
        s.intma(1, 2) = 1; s.intma(2, 2) = 3;
        s.intma(3, 2) = 2; s.intma(4, 2) = 5;

        s.mat_elem.resize(2);
        s.mat_elem(1) = 0;
        s.mat_elem(2) = 1;

        s.node_material_mask.resize(5);
        s.node_material_mask(1) = 3;
        s.node_material_mask(2) = 3;
        s.node_material_mask(3) = 3;
        s.node_material_mask(4) = CBSStateSI::node_touches_fluid;
        s.node_material_mask(5) = CBSStateSI::node_touches_solid;

        s.detJ.resize(2);
        s.detJ(1) = 1.0e-3; // A=0.5 -> interface altitude = 1 mm
        s.detJ(2) = 1.0e-3;

        s.annxf.resize(4, 4, 2);
        s.annxf.fill(0.0);
        // Shared face is local face 4 on both tets.
        s.annxf(3, 4, 1) = -0.5;
        s.annxf(4, 4, 1) =  0.5;
        s.annxf(3, 4, 2) =  0.5;
        s.annxf(4, 4, 2) =  0.5;

        s.fedge.resize(4, 2);
        s.fedge.fill(0);

        s.dNkdx.resize(24);
        s.dNkdx.fill(0.0);

        // Fluid tet: x=N2, y=N3, z=0.001*N4.
        const Real gfluid[3][4] =
        {
            {-1.0, 1.0, 0.0, 0.0},
            {-1.0, 0.0, 1.0, 0.0},
            {-1000.0, 0.0, 0.0, 1000.0}
        };

        for (Int dim = 1; dim <= 3; ++dim)
        {
            for (Int a = 1; a <= 4; ++a)
            {
                s.dNkdx(grad_index(s, 1, dim, a)) =
                    gfluid[dim - 1][a - 1];
            }
        }

        s.rho_e.resize(2);
        s.cp_e.resize(2);
        s.rho_cp_e.resize(2);
        s.mu_e.resize(2);
        s.k_e.resize(2);
        s.nu_t_e.resize(2);
        s.k_eff_e.resize(2);

        s.rho_e(1) = 6.7;
        s.cp_e(1) = 5200.0;
        s.rho_cp_e(1) = s.rho_e(1) * s.cp_e(1);
        s.mu_e(1) = 3.1e-5;
        s.k_e(1) = 0.24;
        s.nu_t_e(1) = 1.0e-5;
        s.k_eff_e(1) = s.k_e(1)
            + s.rho_cp_e(1) * s.nu_t_e(1) / s.cfg.sa_prandtl_t;

        s.rho_e(2) = 7650.0;
        s.cp_e(2) = 560.0;
        s.rho_cp_e(2) = s.rho_e(2) * s.cp_e(2);
        s.mu_e(2) = 0.0;
        s.k_e(2) = 28.3;
        s.nu_t_e(2) = 0.0;
        s.k_eff_e(2) = s.k_e(2);

        s.unkn1.resize(3, 5);
        s.unkno.resize(3, 5);
        s.unkn1.fill(0.0);
        s.unkno.fill(0.0);
        // Opposite fluid sample moves tangentially in +x.
        s.unkn1(1, 4) = 10.0;
        s.unkno(1, 4) = 10.0;

        s.rhs.resize(3, 5);
        s.rhs.fill(0.0);

        s.temperature1.resize(5);
        s.temperature.resize(5);
        s.rhs1.resize(5);
        s.rhs1.fill(0.0);

        s.temperature1(1) = 650.0;
        s.temperature1(2) = 650.0;
        s.temperature1(3) = 650.0;
        s.temperature1(4) = 600.0;
        s.temperature1(5) = 700.0;
        for (Int ip = 1; ip <= 5; ++ip)
        {
            s.temperature(ip) = s.temperature1(ip);
        }

        s.delte.resize(2);
        s.delte(1) = 1.0e-7;
        s.delte(2) = 1.0e-7;

        s.node_velocity_bc_type.resize(5);
        s.node_velocity_bc_type.fill(CBSStateSI::velocity_bc_free);
        for (const Int ip : {1, 2, 3})
        {
            s.node_velocity_bc_type(ip) = CBSStateSI::velocity_bc_noslip;
        }

        return s;
    }
}

int main()
{
#if defined(_WIN32)
    _putenv_s("CBS3D_SA_WALL_TREATMENT", "");
    _putenv_s("CBS3D_CHT_WALL_TREATMENT", "1");
#else
    unsetenv("CBS3D_SA_WALL_TREATMENT");
    setenv("CBS3D_CHT_WALL_TREATMENT", "1", 1);
#endif

    try
    {
        CBSStateSI s = make_state();

        if (CHTWallModelCoupling::globalWallFaceCount(s) != 1)
        {
            std::printf("FAIL: expected exactly one conformal CHT face\n");
            return 1;
        }

        for (const Int ip : {1, 2, 3})
        {
            if (!CHTWallModelCoupling::isModelWallNode(s, ip))
            {
                std::printf("FAIL: interface node %d not classified\n", ip);
                return 1;
            }
        }

        if (CHTWallModelCoupling::isModelWallNode(s, 4) ||
            CHTWallModelCoupling::isModelWallNode(s, 5))
        {
            std::printf("FAIL: non-interface node classified as CHT wall\n");
            return 1;
        }

        const auto momentum = CHTWallModelCoupling::addMomentumWallFlux(s);

        if (momentum.local_faces != 1 ||
            !(momentum.modeled_surface_load[0] < 0.0) ||
            !(momentum.modeled_wall_work <= 0.0) ||
            !(momentum.minimum_y_plus > 0.0))
        {
            std::printf("FAIL: CHT Spalding momentum diagnostics\n");
            return 1;
        }

        // Verify tangent-space restoration: the interface normal is z, so x/y
        // survive while z is removed.
        for (const Int ip : {1, 2, 3})
        {
            s.unkno(1, ip) = 3.0;
            s.unkno(2, ip) = 4.0;
            s.unkno(3, ip) = 5.0;
        }

        const auto captured = CHTWallModelCoupling::captureVelocity(s, false);
        for (const Int ip : {1, 2, 3})
        {
            s.unkno(1, ip) = 0.0;
            s.unkno(2, ip) = 0.0;
            s.unkno(3, ip) = 0.0;
        }
        CHTWallModelCoupling::restoreTangentialAndEnforceImpermeability(
            s, captured);

        for (const Int ip : {1, 2, 3})
        {
            if (!close(s.unkno(1, ip), 3.0) ||
                !close(s.unkno(2, ip), 4.0) ||
                !close(s.unkno(3, ip), 0.0))
            {
                std::printf("FAIL: CHT interface tangent projection\n");
                return 1;
            }
        }

        // Re-establish the representative old velocity used by the wall laws.
        s.unkno.fill(0.0);
        s.unkno(1, 4) = 10.0;

        CHTWallModelCoupling::prepareThermalStabilityConductivity(s);

        const Real k_bulk = s.k_e(1)
            + s.rho_cp_e(1) * s.nu_t_e(1) / s.cfg.sa_prandtl_t;
        const Real k_assembled = s.k_eff_e(1);

        if (!(k_assembled >= k_bulk) || !std::isfinite(k_assembled))
        {
            std::printf("FAIL: thermal stability conductivity not conservative\n");
            return 1;
        }

        // Assemble exactly the fluid isotropic diffusion term that
        // EnergyAssembly sees after the stability preparation.
        const Real volume = s.detJ(1) / 6.0;
        const Real gradTz = -50000.0; // (600-650)/0.001

        for (Int a = 1; a <= 4; ++a)
        {
            const Int ip = s.intma(a, 1);
            const Real gradNa_z = s.dNkdx(grad_index(s, 1, 3, a));
            s.rhs1(ip) = -k_assembled * volume * gradNa_z * gradTz;
        }

        // Sentinel on the solid-only node. The CHT wall correction may alter
        // shared interface residuals (that is the fluid contribution), but must
        // never touch the solid-only node or solid-element operator.
        s.rhs1(5) = 123.456;

        const auto thermal = CHTWallModelCoupling::correctThermalWallDiffusion(s);

        const auto wall = WallTreatment::evaluateSpalding(
            {10.0, 0.0, 0.0},
            {0.0, 0.0, -1.0},
            1.0e-3,
            s.rho_e(1),
            s.mu_e(1) / s.rho_e(1));

        const auto kader = ThermalWallTreatment::evaluateKader(
            wall.friction_velocity,
            1.0e-3,
            s.rho_e(1),
            s.cp_e(1),
            s.mu_e(1),
            s.k_e(1));

        for (Int a = 1; a <= 4; ++a)
        {
            const Int ip = s.intma(a, 1);
            const Real gradNa_z = s.dNkdx(grad_index(s, 1, 3, a));
            const Real expected =
                -kader.wall_normal_conductivity *
                volume * gradNa_z * gradTz;

            if (!close(s.rhs1(ip), expected, 2.0e-10, 2.0e-9))
            {
                std::printf(
                    "FAIL: thermal tensor node=%d got=% .17e expected=% .17e\n",
                    ip, s.rhs1(ip), expected);
                return 1;
            }
        }

        if (!close(s.rhs1(5), 123.456) ||
            !close(thermal.residual_correction_sum, 0.0, 0.0, 2.0e-9) ||
            thermal.local_faces != 1 ||
            !(thermal.maximum_kn_over_k >= thermal.minimum_kn_over_k) ||
            !(thermal.minimum_prandtl > 0.0))
        {
            std::printf("FAIL: thermal conservation/material isolation diagnostics\n");
            return 1;
        }
    }
    catch (const std::exception& error)
    {
        std::printf("FAIL: unexpected exception: %s\n", error.what());
        return 1;
    }

    std::printf("PASS: conformal CHT Spalding/Kader wall coupling\n");
    return 0;
}
