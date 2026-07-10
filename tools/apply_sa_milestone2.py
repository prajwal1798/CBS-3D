#!/usr/bin/env python3
"""
Apply Spalart-Allmaras Milestone 2 scaffolding to CBS3D++_SI.

This script is intentionally conservative:
  - it creates the turbulence/boundary/assembly scaffold files;
  - it adds turbulence controls to RunConfig;
  - it adds turbulence fields to CBSStateSI;
  - it allocates and initialises the new fields;
  - it adds the new source files to CMake;
  - it keeps turbulence disabled by default, so laminar behaviour remains unchanged.

Run from the repository root:

    python tools/apply_sa_milestone2.py

Then build normally and commit the resulting source changes.
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read(path: str) -> str:
    return (ROOT / path).read_text(encoding="utf-8")


def write(path: str, text: str) -> None:
    target = ROOT / path
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(text, encoding="utf-8", newline="\n")


def replace_once(text: str, old: str, new: str, path: str) -> str:
    if old not in text:
        raise RuntimeError(f"Required patch marker not found in {path}: {old!r}")
    return text.replace(old, new, 1)


def patch_run_config() -> None:
    path = "include/cbs/core/RunConfig.hpp"
    text = read(path)

    if "turbulence_on" in text:
        print(f"[skip] {path}: turbulence controls already present")
        return

    marker = "        Real heat_flux_bc = 0.0;\n"
    insertion = marker + r'''

        // --------------------------------------------------------------------
        // Spalart-Allmaras turbulence controls
        // --------------------------------------------------------------------
        // turbulence_on:
        //     0  laminar solver path; all turbulence arrays remain passive
        //     1  solve and apply the selected one-equation turbulence model
        //
        // turbulence_model:
        //     0  standard Spalart-Allmaras model
        //     1  SA-neg compatible branch, reserved for the robust variant
        //
        // turbulent_thermal_diffusivity_on:
        //     0  turbulent viscosity affects momentum only
        //     1  fluid thermal conductivity receives rho*cp*nu_t/Pr_t
        Int turbulence_on = 0;
        Int turbulence_model = 0;
        Int turbulent_thermal_diffusivity_on = 0;

        // Fully turbulent inlet/farfield working-variable ratio:
        //     nu_tilde_inlet = sa_inlet_ratio * nu_molecular
        Real sa_inlet_ratio = 3.0;

        // Turbulent Prandtl number used only for fluid heat-transfer coupling.
        Real sa_prandtl_t = 0.90;

        // Numerical protections for the first non-negative SA implementation.
        Real sa_min_wall_distance = 1.0e-14;
        Real sa_min_stilde = 1.0e-14;
        Real sa_nu_tilde_floor = 0.0;

        Int sa_use_stilde_limiter = 1;
        Int sa_implicit_destruction = 1;
'''
    text = replace_once(text, marker, insertion, path)
    write(path, text)
    print(f"[patch] {path}: added turbulence controls")


def patch_cbs_state() -> None:
    path = "include/cbs/core/CBSStateSI.hpp"
    text = read(path)

    if "Array1D<Real> nu_tilde;" not in text:
        marker = "        Array1D<Real> Qvol_e;\n"
        insertion = marker + r'''

        // --------------------------------------------------------------------
        // Spalart-Allmaras turbulence state
        // --------------------------------------------------------------------
        // nu_tilde is the transported SA working variable.  nu_t and mu_t are
        // derived quantities and must never be prescribed independently.
        Array1D<Real> nu_tilde;
        Array1D<Real> nu_tilde1;
        Array1D<Real> nu_t;
        Array1D<Real> mu_t;

        // True minimum distance from a fluid node to the nearest physical
        // no-slip wall triangle.  This is geometric data, not a mesh-line
        // index or nearest-wall-node approximation.
        Array1D<Real> wall_distance;

        // SA nodal assembly and diagnostics.
        Array1D<Real> sa_rhs;
        Array1D<Real> sa_source;
        Array1D<Real> sa_production;
        Array1D<Real> sa_destruction;
        Array1D<Real> sa_diffusion;
        Array1D<Real> sa_residual;

        // Nodal turbulence classification flags.
        Array1D<Int> sa_active_node;
        Array1D<Int> sa_wall_node;
        Array1D<Int> sa_inlet_node;

        // Element-averaged turbulence quantities and effective properties used
        // by momentum and energy assemblies.  Molecular values remain stored in
        // mu_e and k_e and are not overwritten.
        Array1D<Real> nu_tilde_e;
        Array1D<Real> nu_t_e;
        Array1D<Real> mu_t_e;
        Array1D<Real> mu_eff_e;
        Array1D<Real> k_eff_e;
'''
        text = replace_once(text, marker, insertion, path)
        print(f"[patch] {path}: added turbulence fields")
    else:
        print(f"[skip] {path}: turbulence fields already present")

    if "nu_tilde.resize(cfg.npoin);" not in text:
        marker = "            Qvol_e.resize(cfg.nelem);\n"
        insertion = marker + r'''

            // Spalart-Allmaras nodal fields.
            nu_tilde.resize(cfg.npoin);
            nu_tilde1.resize(cfg.npoin);
            nu_t.resize(cfg.npoin);
            mu_t.resize(cfg.npoin);
            wall_distance.resize(cfg.npoin);

            sa_rhs.resize(cfg.npoin);
            sa_source.resize(cfg.npoin);
            sa_production.resize(cfg.npoin);
            sa_destruction.resize(cfg.npoin);
            sa_diffusion.resize(cfg.npoin);
            sa_residual.resize(cfg.npoin);

            sa_active_node.resize(cfg.npoin);
            sa_wall_node.resize(cfg.npoin);
            sa_inlet_node.resize(cfg.npoin);

            // Spalart-Allmaras element fields.
            nu_tilde_e.resize(cfg.nelem);
            nu_t_e.resize(cfg.nelem);
            mu_t_e.resize(cfg.nelem);
            mu_eff_e.resize(cfg.nelem);
            k_eff_e.resize(cfg.nelem);
'''
        text = replace_once(text, marker, insertion, path)
        print(f"[patch] {path}: added turbulence allocation")
    else:
        print(f"[skip] {path}: turbulence allocation already present")

    if "nu_tilde.fill(0.0);" not in text:
        marker = "            Qvol_e.fill(0.0);\n"
        insertion = marker + r'''

            nu_tilde.fill(0.0);
            nu_tilde1.fill(0.0);
            nu_t.fill(0.0);
            mu_t.fill(0.0);
            wall_distance.fill(1.0e300);

            sa_rhs.fill(0.0);
            sa_source.fill(0.0);
            sa_production.fill(0.0);
            sa_destruction.fill(0.0);
            sa_diffusion.fill(0.0);
            sa_residual.fill(0.0);

            sa_active_node.fill(0);
            sa_wall_node.fill(0);
            sa_inlet_node.fill(0);

            nu_tilde_e.fill(0.0);
            nu_t_e.fill(0.0);
            mu_t_e.fill(0.0);
            mu_eff_e.fill(0.0);
            k_eff_e.fill(1.0);
'''
        text = replace_once(text, marker, insertion, path)
        print(f"[patch] {path}: added turbulence default initialisation")
    else:
        print(f"[skip] {path}: turbulence default initialisation already present")

    write(path, text)


def patch_mesh_io() -> None:
    path = "src/io/MeshIO.cpp"
    text = read(path)

    if "s.mu_eff_e(ie) = s.mu_e(ie);" not in text:
        marker = "                s.Qvol_e(ie) = 0.0;\n"
        insertion = marker + r'''
                s.mu_eff_e(ie) = s.mu_e(ie);
                s.k_eff_e(ie) = s.k_e(ie);
                s.nu_tilde_e(ie) = 0.0;
                s.nu_t_e(ie) = 0.0;
                s.mu_t_e(ie) = 0.0;
'''
        text = text.replace(marker, insertion)
        print(f"[patch] {path}: effective default properties tied to molecular values")
    else:
        print(f"[skip] {path}: effective default properties already patched")

    if "s.mu_eff_e(ie) = s.mu_e(ie);\n                s.k_eff_e(ie) = s.k_e(ie);" not in text:
        marker = "            s.alpha_e(ie) = s.k_e(ie) / s.rho_cp_e(ie);\n"
        insertion = marker + r'''
            s.mu_eff_e(ie) = s.mu_e(ie);
            s.k_eff_e(ie) = s.k_e(ie);
            s.nu_tilde_e(ie) = 0.0;
            s.nu_t_e(ie) = 0.0;
            s.mu_t_e(ie) = 0.0;
'''
        text = replace_once(text, marker, insertion, path)
        print(f"[patch] {path}: effective material properties after .matprop load")

    write(path, text)


def patch_cmake() -> None:
    path = "CMakeLists.txt"
    text = read(path)

    if "src/turbulence/SpalartAllmaras.cpp" in text:
        print(f"[skip] {path}: SA sources already present")
        return

    marker = "    src/assembly/EnergyAssembly.cpp\n"
    insertion = marker + r'''
    src/assembly/SpalartAllmarasAssembly.cpp

    src/turbulence/SpalartAllmaras.cpp
    src/turbulence/WallDistance.cpp

    src/boundary/TurbulenceBoundary.cpp
'''
    text = replace_once(text, marker, insertion, path)
    write(path, text)
    print(f"[patch] {path}: added SA source files")


def write_sa_core_files() -> None:
    write("include/cbs/turbulence/SpalartAllmaras.hpp", r'''#pragma once

//=============================================================================
// CBS3D++_SI
//
// Pure algebraic functions and constants for the Spalart-Allmaras turbulence
// model.  These routines do not access mesh storage and are therefore suitable
// for deterministic unit testing before the transport equation is coupled into
// the CBS solver.
//=============================================================================

#include "cbs/core/Types.hpp"

namespace cbs::turbulence
{
    struct SpalartAllmarasConstants
    {
        Real cb1 = 0.1355;
        Real sigma = 2.0 / 3.0;
        Real cb2 = 0.622;
        Real kappa = 0.41;
        Real cw2 = 0.3;
        Real cw3 = 2.0;
        Real cv1 = 7.1;
        Real ct3 = 1.2;
        Real ct4 = 0.5;

        // Allmaras-Johnson-Spalart S_tilde limiter constants.
        Real c2 = 0.7;
        Real c3 = 0.9;

        // SA-neg reserved constant.  The negative branch is not activated in
        // Milestone 2 but the value is kept here so the public API is stable.
        Real cn1 = 16.0;

        [[nodiscard]] Real cw1() const noexcept;
    };

    [[nodiscard]] Real chi(Real nu_tilde, Real molecular_nu);
    [[nodiscard]] Real fv1(Real chi_value, const SpalartAllmarasConstants& c = {});
    [[nodiscard]] Real fv2(Real chi_value, Real fv1_value);
    [[nodiscard]] Real ft2(Real chi_value, const SpalartAllmarasConstants& c = {});

    [[nodiscard]] Real sBar(
        Real nu_tilde,
        Real molecular_nu,
        Real wall_distance,
        const SpalartAllmarasConstants& c = {});

    [[nodiscard]] Real limitedSTilde(
        Real omega,
        Real s_bar,
        const SpalartAllmarasConstants& c = {});

    [[nodiscard]] Real rFunction(
        Real nu_tilde,
        Real s_tilde,
        Real wall_distance,
        const SpalartAllmarasConstants& c = {});

    [[nodiscard]] Real fw(Real r_value, const SpalartAllmarasConstants& c = {});

    [[nodiscard]] Real eddyKinematicViscosity(
        Real nu_tilde,
        Real molecular_nu,
        const SpalartAllmarasConstants& c = {});

    [[nodiscard]] Real negativeBranchFn(
        Real chi_value,
        const SpalartAllmarasConstants& c = {});

    [[nodiscard]] Real productionTerm(
        Real nu_tilde,
        Real s_tilde,
        Real ft2_value,
        const SpalartAllmarasConstants& c = {});

    [[nodiscard]] Real destructionCoefficient(
        Real wall_distance,
        Real fw_value,
        Real ft2_value,
        const SpalartAllmarasConstants& c = {});

    [[nodiscard]] Real destructionTerm(
        Real nu_tilde,
        Real wall_distance,
        Real fw_value,
        Real ft2_value,
        const SpalartAllmarasConstants& c = {});
}
''')

    write("src/turbulence/SpalartAllmaras.cpp", r'''#include "cbs/turbulence/SpalartAllmaras.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cbs::turbulence
{
    namespace
    {
        constexpr Real tiny = 1.0e-30;

        [[nodiscard]] Real safe_positive(const Real value, const char* name)
        {
            if (!std::isfinite(value) || value <= 0.0)
            {
                throw std::runtime_error(
                    std::string("SpalartAllmaras - non-positive or non-finite ") + name);
            }
            return value;
        }
    }

    Real SpalartAllmarasConstants::cw1() const noexcept
    {
        return cb1 / (kappa * kappa) + (1.0 + cb2) / sigma;
    }

    Real chi(const Real nu_tilde, const Real molecular_nu)
    {
        return nu_tilde / safe_positive(molecular_nu, "molecular_nu");
    }

    Real fv1(const Real chi_value, const SpalartAllmarasConstants& c)
    {
        const Real chi3 = chi_value * chi_value * chi_value;
        const Real cv13 = c.cv1 * c.cv1 * c.cv1;
        const Real denominator = chi3 + cv13;

        if (std::abs(denominator) <= tiny)
        {
            return 0.0;
        }

        return chi3 / denominator;
    }

    Real fv2(const Real chi_value, const Real fv1_value)
    {
        const Real denominator = 1.0 + chi_value * fv1_value;

        if (std::abs(denominator) <= tiny)
        {
            return 0.0;
        }

        return 1.0 - chi_value / denominator;
    }

    Real ft2(const Real chi_value, const SpalartAllmarasConstants& c)
    {
        return c.ct3 * std::exp(-c.ct4 * chi_value * chi_value);
    }

    Real sBar(
        const Real nu_tilde,
        const Real molecular_nu,
        const Real wall_distance,
        const SpalartAllmarasConstants& c)
    {
        const Real d = safe_positive(wall_distance, "wall_distance");
        const Real chi_value = chi(nu_tilde, molecular_nu);
        const Real fv1_value = fv1(chi_value, c);
        const Real fv2_value = fv2(chi_value, fv1_value);

        return nu_tilde * fv2_value / (c.kappa * c.kappa * d * d);
    }

    Real limitedSTilde(
        const Real omega,
        const Real s_bar,
        const SpalartAllmarasConstants& c)
    {
        if (!std::isfinite(omega) || omega < 0.0)
        {
            throw std::runtime_error("SpalartAllmaras - omega must be non-negative and finite");
        }

        if (!std::isfinite(s_bar))
        {
            throw std::runtime_error("SpalartAllmaras - s_bar must be finite");
        }

        if (s_bar >= -c.c2 * omega)
        {
            return std::max(omega + s_bar, 0.0);
        }

        const Real denominator = (c.c3 - 2.0 * c.c2) * omega - s_bar;

        if (std::abs(denominator) <= tiny)
        {
            return std::max(omega, 0.0);
        }

        const Real limited = omega
            + omega * (c.c2 * c.c2 * omega + c.c3 * s_bar) / denominator;

        return std::max(limited, 0.0);
    }

    Real rFunction(
        const Real nu_tilde,
        const Real s_tilde,
        const Real wall_distance,
        const SpalartAllmarasConstants& c)
    {
        const Real d = safe_positive(wall_distance, "wall_distance");

        if (!std::isfinite(s_tilde) || s_tilde <= tiny)
        {
            return 10.0;
        }

        const Real denominator = s_tilde * c.kappa * c.kappa * d * d;

        if (denominator <= tiny)
        {
            return 10.0;
        }

        return std::clamp(nu_tilde / denominator, 0.0, 10.0);
    }

    Real fw(const Real r_value, const SpalartAllmarasConstants& c)
    {
        const Real r = std::clamp(r_value, 0.0, 10.0);
        const Real r6 = std::pow(r, 6.0);
        const Real g = r + c.cw2 * (r6 - r);
        const Real g6 = std::pow(g, 6.0);
        const Real cw36 = std::pow(c.cw3, 6.0);

        return g * std::pow((1.0 + cw36) / (g6 + cw36), 1.0 / 6.0);
    }

    Real eddyKinematicViscosity(
        const Real nu_tilde,
        const Real molecular_nu,
        const SpalartAllmarasConstants& c)
    {
        if (nu_tilde <= 0.0)
        {
            return 0.0;
        }

        const Real chi_value = chi(nu_tilde, molecular_nu);
        return nu_tilde * fv1(chi_value, c);
    }

    Real negativeBranchFn(
        const Real chi_value,
        const SpalartAllmarasConstants& c)
    {
        const Real chi3 = chi_value * chi_value * chi_value;
        const Real denominator = c.cn1 - chi3;

        if (std::abs(denominator) <= tiny)
        {
            throw std::runtime_error("SpalartAllmaras - singular SA-neg fn denominator");
        }

        return (c.cn1 + chi3) / denominator;
    }

    Real productionTerm(
        const Real nu_tilde,
        const Real s_tilde,
        const Real ft2_value,
        const SpalartAllmarasConstants& c)
    {
        return c.cb1 * (1.0 - ft2_value) * s_tilde * nu_tilde;
    }

    Real destructionCoefficient(
        const Real wall_distance,
        const Real fw_value,
        const Real ft2_value,
        const SpalartAllmarasConstants& c)
    {
        const Real d = safe_positive(wall_distance, "wall_distance");
        return (c.cw1() * fw_value - c.cb1 * ft2_value / (c.kappa * c.kappa)) / (d * d);
    }

    Real destructionTerm(
        const Real nu_tilde,
        const Real wall_distance,
        const Real fw_value,
        const Real ft2_value,
        const SpalartAllmarasConstants& c)
    {
        return destructionCoefficient(wall_distance, fw_value, ft2_value, c)
            * nu_tilde * nu_tilde;
    }
}
''')

    print("[write] SA algebra files")


def write_boundary_files() -> None:
    write("include/cbs/boundary/TurbulenceBoundary.hpp", r'''#pragma once

//=============================================================================
// CBS3D++_SI
//
// Boundary and initial-value handling for the Spalart-Allmaras working variable.
// Partition boundaries are never treated as turbulence walls.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"

namespace cbs
{
    class TurbulenceBoundary
    {
    public:
        static void classifyNodes(CBSStateSI& s);
        static void initialiseNuTilde(CBSStateSI& s);
        static void applyWallValues(CBSStateSI& s);
        static void applyInletValues(CBSStateSI& s);
    };
}
''')

    write("src/boundary/TurbulenceBoundary.cpp", r'''#include "cbs/boundary/TurbulenceBoundary.hpp"

#include <algorithm>
#include <vector>

namespace cbs
{
    namespace
    {
        [[nodiscard]] bool is_fluid_element(const CBSStateSI& s, const Int ie)
        {
            return s.mat_elem(ie) == 0;
        }

        [[nodiscard]] bool is_wall_bc(const CBSStateSI& s, const Int bc)
        {
            return bc == s.cfg.bc_noslip_adiabatic_wall
                || bc == s.cfg.bc_noslip_heatflux_wall
                || bc == s.cfg.bc_cht_interface;
        }

        [[nodiscard]] bool is_inlet_bc(const CBSStateSI& s, const Int bc)
        {
            return bc == s.cfg.bc_velocity_temperature_inlet
                || bc == s.cfg.bc_massflow_temperature_inlet;
        }

        [[nodiscard]] Real nodal_molecular_nu(
            const CBSStateSI& s,
            const Int ip,
            const std::vector<Int>& count)
        {
            Real sum = 0.0;

            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                if (!is_fluid_element(s, ie))
                {
                    continue;
                }

                bool touches_node = false;
                for (Int in = 1; in <= s.cfg.nep; ++in)
                {
                    touches_node = touches_node || (s.intma(in, ie) == ip);
                }

                if (touches_node && s.rho_e(ie) > 0.0)
                {
                    sum += s.mu_e(ie) / s.rho_e(ie);
                }
            }

            const Int n = count[static_cast<std::size_t>(ip)];
            return n > 0 ? sum / static_cast<Real>(n) : 0.0;
        }
    }

    void TurbulenceBoundary::classifyNodes(CBSStateSI& s)
    {
        s.sa_active_node.fill(0);
        s.sa_wall_node.fill(0);
        s.sa_inlet_node.fill(0);

        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            if (!is_fluid_element(s, ie))
            {
                continue;
            }

            for (Int in = 1; in <= s.cfg.nep; ++in)
            {
                s.sa_active_node(s.intma(in, ie)) = 1;
            }
        }

        for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
        {
            const Int bc = s.iside(s.cfg.bsid, ib);

            if (!is_wall_bc(s, bc) && !is_inlet_bc(s, bc))
            {
                continue;
            }

            for (Int i = 1; i <= s.cfg.nsidp; ++i)
            {
                const Int ip = s.iside(i, ib);

                if (is_wall_bc(s, bc))
                {
                    s.sa_wall_node(ip) = 1;
                }

                if (is_inlet_bc(s, bc))
                {
                    s.sa_inlet_node(ip) = 1;
                }
            }
        }
    }

    void TurbulenceBoundary::initialiseNuTilde(CBSStateSI& s)
    {
        classifyNodes(s);

        s.nu_tilde.fill(0.0);
        s.nu_tilde1.fill(0.0);
        s.nu_t.fill(0.0);
        s.mu_t.fill(0.0);

        if (s.cfg.turbulence_on < 1)
        {
            return;
        }

        std::vector<Int> fluid_touch_count(static_cast<std::size_t>(s.cfg.npoin + 1), 0);

        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            if (!is_fluid_element(s, ie))
            {
                continue;
            }

            for (Int in = 1; in <= s.cfg.nep; ++in)
            {
                ++fluid_touch_count[static_cast<std::size_t>(s.intma(in, ie))];
            }
        }

        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            if (s.sa_active_node(ip) == 0 || s.sa_wall_node(ip) != 0)
            {
                continue;
            }

            const Real nu = nodal_molecular_nu(s, ip, fluid_touch_count);
            s.nu_tilde(ip) = std::max(0.0, s.cfg.sa_inlet_ratio * nu);
            s.nu_tilde1(ip) = s.nu_tilde(ip);
        }

        applyWallValues(s);
        applyInletValues(s);
    }

    void TurbulenceBoundary::applyWallValues(CBSStateSI& s)
    {
        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            if (s.sa_wall_node(ip) != 0)
            {
                s.nu_tilde(ip) = 0.0;
                s.nu_tilde1(ip) = 0.0;
                s.nu_t(ip) = 0.0;
                s.mu_t(ip) = 0.0;
            }
        }
    }

    void TurbulenceBoundary::applyInletValues(CBSStateSI& s)
    {
        if (s.cfg.turbulence_on < 1)
        {
            return;
        }

        std::vector<Int> fluid_touch_count(static_cast<std::size_t>(s.cfg.npoin + 1), 0);

        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            if (!is_fluid_element(s, ie))
            {
                continue;
            }

            for (Int in = 1; in <= s.cfg.nep; ++in)
            {
                ++fluid_touch_count[static_cast<std::size_t>(s.intma(in, ie))];
            }
        }

        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            if (s.sa_inlet_node(ip) != 0 && s.sa_wall_node(ip) == 0)
            {
                const Real nu = nodal_molecular_nu(s, ip, fluid_touch_count);
                s.nu_tilde(ip) = std::max(0.0, s.cfg.sa_inlet_ratio * nu);
                s.nu_tilde1(ip) = s.nu_tilde(ip);
            }
        }
    }
}
''')

    print("[write] turbulence boundary files")


def write_wall_distance_files() -> None:
    write("include/cbs/turbulence/WallDistance.hpp", r'''#pragma once

//=============================================================================
// CBS3D++_SI
//
// Wall-distance computation for the Spalart-Allmaras model.  The required
// quantity is the true minimum Euclidean distance from each active fluid node to
// the nearest physical no-slip wall triangle.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"

#include <array>
#include <vector>

namespace cbs::turbulence
{
    struct WallTriangle
    {
        std::array<Real, 3> a{};
        std::array<Real, 3> b{};
        std::array<Real, 3> c{};
    };

    class WallDistance
    {
    public:
        static std::vector<WallTriangle> collectPhysicalWallTriangles(const CBSStateSI& s);
        static Real pointTriangleDistance(const std::array<Real, 3>& p, const WallTriangle& tri);
        static void computeSerial(CBSStateSI& s);
    };
}
''')

    write("src/turbulence/WallDistance.cpp", r'''#include "cbs/turbulence/WallDistance.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cbs::turbulence
{
    namespace
    {
        [[nodiscard]] std::array<Real, 3> node_point(const CBSStateSI& s, const Int ip)
        {
            return {s.coord(1, ip), s.coord(2, ip), s.coord(3, ip)};
        }

        [[nodiscard]] bool is_wall_bc(const CBSStateSI& s, const Int bc)
        {
            return bc == s.cfg.bc_noslip_adiabatic_wall
                || bc == s.cfg.bc_noslip_heatflux_wall
                || bc == s.cfg.bc_cht_interface;
        }

        [[nodiscard]] std::array<Real, 3> sub(
            const std::array<Real, 3>& a,
            const std::array<Real, 3>& b)
        {
            return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
        }

        [[nodiscard]] std::array<Real, 3> add_scaled(
            const std::array<Real, 3>& a,
            const std::array<Real, 3>& b,
            const Real scale)
        {
            return {a[0] + scale * b[0], a[1] + scale * b[1], a[2] + scale * b[2]};
        }

        [[nodiscard]] Real dot(
            const std::array<Real, 3>& a,
            const std::array<Real, 3>& b)
        {
            return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        }

        [[nodiscard]] Real distance_squared(
            const std::array<Real, 3>& a,
            const std::array<Real, 3>& b)
        {
            const auto d = sub(a, b);
            return dot(d, d);
        }
    }

    std::vector<WallTriangle> WallDistance::collectPhysicalWallTriangles(const CBSStateSI& s)
    {
        std::vector<WallTriangle> triangles;
        triangles.reserve(static_cast<std::size_t>(s.cfg.nboun));

        for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
        {
            const Int bc = s.iside(s.cfg.bsid, ib);

            if (!is_wall_bc(s, bc))
            {
                continue;
            }

            WallTriangle tri;
            tri.a = node_point(s, s.iside(1, ib));
            tri.b = node_point(s, s.iside(2, ib));
            tri.c = node_point(s, s.iside(3, ib));
            triangles.push_back(tri);
        }

        return triangles;
    }

    Real WallDistance::pointTriangleDistance(
        const std::array<Real, 3>& p,
        const WallTriangle& tri)
    {
        // Real-Time Collision Detection, Christer Ericson, point-triangle
        // closest-point test.  The formula is purely geometric and does not
        // assume structured wall-normal mesh lines.
        const auto ab = sub(tri.b, tri.a);
        const auto ac = sub(tri.c, tri.a);
        const auto ap = sub(p, tri.a);

        const Real d1 = dot(ab, ap);
        const Real d2 = dot(ac, ap);
        if (d1 <= 0.0 && d2 <= 0.0)
        {
            return std::sqrt(distance_squared(p, tri.a));
        }

        const auto bp = sub(p, tri.b);
        const Real d3 = dot(ab, bp);
        const Real d4 = dot(ac, bp);
        if (d3 >= 0.0 && d4 <= d3)
        {
            return std::sqrt(distance_squared(p, tri.b));
        }

        const Real vc = d1 * d4 - d3 * d2;
        if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
        {
            const Real v = d1 / (d1 - d3);
            const auto closest = add_scaled(tri.a, ab, v);
            return std::sqrt(distance_squared(p, closest));
        }

        const auto cp = sub(p, tri.c);
        const Real d5 = dot(ab, cp);
        const Real d6 = dot(ac, cp);
        if (d6 >= 0.0 && d5 <= d6)
        {
            return std::sqrt(distance_squared(p, tri.c));
        }

        const Real vb = d5 * d2 - d1 * d6;
        if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
        {
            const Real w = d2 / (d2 - d6);
            const auto closest = add_scaled(tri.a, ac, w);
            return std::sqrt(distance_squared(p, closest));
        }

        const Real va = d3 * d6 - d5 * d4;
        if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
        {
            const Real w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
            const auto bc = sub(tri.c, tri.b);
            const auto closest = add_scaled(tri.b, bc, w);
            return std::sqrt(distance_squared(p, closest));
        }

        const Real denom = 1.0 / (va + vb + vc);
        const Real v = vb * denom;
        const Real w = vc * denom;
        const auto closest = add_scaled(add_scaled(tri.a, ab, v), ac, w);
        return std::sqrt(distance_squared(p, closest));
    }

    void WallDistance::computeSerial(CBSStateSI& s)
    {
        const auto wall_triangles = collectPhysicalWallTriangles(s);

        if (wall_triangles.empty())
        {
            throw std::runtime_error(
                "WallDistance::computeSerial - no physical no-slip wall triangles found");
        }

        s.wall_distance.fill(std::numeric_limits<Real>::max());

        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            if (s.sa_active_node(ip) == 0)
            {
                continue;
            }

            if (s.sa_wall_node(ip) != 0)
            {
                s.wall_distance(ip) = s.cfg.sa_min_wall_distance;
                continue;
            }

            const auto p = node_point(s, ip);
            Real d_min = std::numeric_limits<Real>::max();

            for (const WallTriangle& tri : wall_triangles)
            {
                d_min = std::min(d_min, pointTriangleDistance(p, tri));
            }

            s.wall_distance(ip) = std::max(d_min, s.cfg.sa_min_wall_distance);
        }
    }
}
''')

    print("[write] wall-distance files")


def write_assembly_files() -> None:
    write("include/cbs/assembly/SpalartAllmarasAssembly.hpp", r'''#pragma once

//=============================================================================
// CBS3D++_SI
//
// Spalart-Allmaras turbulence assembly and effective-property update entry
// points.  Milestone 2 intentionally provides the state-safe scaffolding and
// algebraic eddy-viscosity update, while the full transport equation is added in
// a later milestone.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"

namespace cbs
{
    class SpalartAllmarasAssembly
    {
    public:
        static void resetEffectiveProperties(CBSStateSI& s);
        static void updateEddyViscosity(CBSStateSI& s);
    };
}
''')

    write("src/assembly/SpalartAllmarasAssembly.cpp", r'''#include "cbs/assembly/SpalartAllmarasAssembly.hpp"

#include "cbs/turbulence/SpalartAllmaras.hpp"

#include <algorithm>
#include <vector>

namespace cbs
{
    namespace
    {
        [[nodiscard]] bool is_fluid_element(const CBSStateSI& s, const Int ie)
        {
            return s.mat_elem(ie) == 0;
        }
    }

    void SpalartAllmarasAssembly::resetEffectiveProperties(CBSStateSI& s)
    {
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

        std::vector<Real> nodal_nu_t_sum(static_cast<std::size_t>(s.cfg.npoin + 1), 0.0);
        std::vector<Real> nodal_mu_t_sum(static_cast<std::size_t>(s.cfg.npoin + 1), 0.0);
        std::vector<Int> nodal_count(static_cast<std::size_t>(s.cfg.npoin + 1), 0);

        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            if (!is_fluid_element(s, ie))
            {
                continue;
            }

            Real nu_tilde_avg = 0.0;
            for (Int in = 1; in <= s.cfg.nep; ++in)
            {
                nu_tilde_avg += std::max(s.cfg.sa_nu_tilde_floor, s.nu_tilde(s.intma(in, ie)));
            }
            nu_tilde_avg /= static_cast<Real>(s.cfg.nep);

            const Real molecular_nu = s.rho_e(ie) > 0.0
                ? s.mu_e(ie) / s.rho_e(ie)
                : 0.0;

            const Real nu_t = molecular_nu > 0.0
                ? turbulence::eddyKinematicViscosity(nu_tilde_avg, molecular_nu)
                : 0.0;

            const Real mu_t = s.rho_e(ie) * nu_t;

            s.nu_tilde_e(ie) = nu_tilde_avg;
            s.nu_t_e(ie) = nu_t;
            s.mu_t_e(ie) = mu_t;
            s.mu_eff_e(ie) = s.mu_e(ie) + mu_t;

            if (s.cfg.turbulent_thermal_diffusivity_on > 0)
            {
                s.k_eff_e(ie) = s.k_e(ie) + s.rho_cp_e(ie) * nu_t / s.cfg.sa_prandtl_t;
            }

            for (Int in = 1; in <= s.cfg.nep; ++in)
            {
                const Int ip = s.intma(in, ie);
                nodal_nu_t_sum[static_cast<std::size_t>(ip)] += nu_t;
                nodal_mu_t_sum[static_cast<std::size_t>(ip)] += mu_t;
                ++nodal_count[static_cast<std::size_t>(ip)];
            }
        }

        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            const Int n = nodal_count[static_cast<std::size_t>(ip)];
            if (n > 0)
            {
                s.nu_t(ip) = nodal_nu_t_sum[static_cast<std::size_t>(ip)] / static_cast<Real>(n);
                s.mu_t(ip) = nodal_mu_t_sum[static_cast<std::size_t>(ip)] / static_cast<Real>(n);
            }
        }
    }
}
''')

    print("[write] SA assembly files")


def main() -> None:
    patch_run_config()
    patch_cbs_state()
    patch_mesh_io()
    patch_cmake()
    write_sa_core_files()
    write_boundary_files()
    write_wall_distance_files()
    write_assembly_files()
    print("\nSA Milestone 2 scaffolding applied.")
    print("Next commands:")
    print("  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCBS3D_ENABLE_PETSC=OFF")
    print("  cmake --build build -j")


if __name__ == "__main__":
    main()
