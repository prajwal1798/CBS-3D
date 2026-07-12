//=============================================================================
// CBS3D++_SI
//
// High-level execution of the semi-implicit CBS steps and the optional
// Spalart-Allmaras turbulence transport step.
//
// The sequence is:
//
//     Step 1   Build the momentum residual and calculate the predicted velocity.
//
//     Step 2   Build and solve the pressure equation.
//
//     Step 3   Correct the velocity using the new pressure gradient.
//
//     Step SA  Advance the Spalart-Allmaras working variable when turbulence is
//              enabled.  This step uses the corrected velocity from Step 3.
//
//     Step 4   Advance the temperature equation when thermal calculation is
//              enabled.
//
// Step SA is not called Step 5 because it is not a fifth CBS splitting step.  It
// is a transported turbulence scalar inserted between velocity correction and
// energy update.
//=============================================================================

#include "cbs/solver/Steps.hpp"

#include "cbs/assembly/EnergyAssembly.hpp"
#include "cbs/assembly/MomentumAssembly.hpp"
#include "cbs/assembly/PressureAssembly.hpp"
#include "cbs/assembly/SpalartAllmarasAssembly.hpp"
#include "cbs/boundary/Boundary.hpp"
#include "cbs/boundary/TurbulenceBoundary.hpp"
#include "cbs/linalg/ConjugateGradient.hpp"

#ifdef CBS3D_USE_PETSC
#include "cbs/linalg/PetscPressureSolver.hpp"
#endif

#include <chrono>
#include <cmath>
#include <iostream>
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


        // Pressure and momentum are assembled only on fluid elements.
        bool is_fluid_element(
            const CBSStateSI& s,
            Int ie)
        {
            return s.mat_elem(ie) == 0;
        }


        // Checks the fixed dimensions required by the current CBS3D
        // P1-tetrahedral formulation.
        void validate_step_dimensions(const CBSStateSI& s)
        {
            if (s.cfg.ndim != 3 ||
                s.cfg.nep != 4 ||
                s.cfg.gdim != 13)
            {
                throw std::runtime_error(
                    "Steps - CBS3D solver steps require ndim=3, nep=4, gdim=13");
            }
        }


        // Returns one Cartesian derivative of a tetrahedral shape function:
        //
        //     grad(N_a)_dim
        Real grad(
            const CBSStateSI& s,
            Int ie,
            Int dim,
            Int local_node)
        {
            return s.dNkdx(dNkdx_index(s, ie, dim, local_node));
        }


        // Enforces the zero-velocity condition in the solid domain:
        //
        //     u = v = w = 0
        //
        // Nodes shared by fluid and solid elements are also set to zero,
        // which gives the no-slip condition at the conformal CHT interface.
        void enforce_zero_velocity_in_solid(CBSStateSI& s)
        {
            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                if (is_fluid_element(s, ie))
                {
                    continue;
                }

                for (Int in = 1; in <= s.cfg.nep; ++in)
                {
                    const Int ip = s.intma(in, ie);

                    for (Int idim = 1; idim <= s.cfg.ndim; ++idim)
                    {
                        s.unkno(idim, ip) = 0.0;
                    }
                }
            }
        }


        // Applies the Step 1 nodal predictor:
        //
        //     u* = u^n + D_u^(-1) r_m
        //
        // where:
        //
        //     D_u^(-1) = elcoe2
        //     r_m       = rhs
        void update_velocity_from_rhs_using_predictor_mass(CBSStateSI& s)
        {
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                for (Int idim = 1; idim <= s.cfg.ndim; ++idim)
                {
                    s.unkno(idim, ip) =
                        s.unkn1(idim, ip) + s.rhs(idim, ip) * s.elcoe2(ip);
                }
            }
        }


        // Applies the Step 3 nodal velocity correction:
        //
        //     u^(n+1) = u* + D_u^(-1) r_p
        //
        // where r_p contains the assembled pressure-gradient contribution.
        void add_velocity_correction_from_rhs(CBSStateSI& s)
        {
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                for (Int idim = 1; idim <= s.cfg.ndim; ++idim)
                {
                    s.unkno(idim, ip) += s.rhs(idim, ip) * s.elcoe2(ip);
                }
            }
        }


        // Reapplies all velocity constraints after Steps 1 and 3.
        //
        // The order is important:
        //
        //     1. symmetry projection
        //     2. prescribed velocity and wall conditions
        //     3. outlet backflow control
        //     4. zero velocity in the material solid
        void apply_velocity_boundary_package(CBSStateSI& s)
        {
            Boundary::applySymmetry(s);
            Boundary::applyVelocity(s);
            Boundary::applyOutletBackflowControl(s);
            enforce_zero_velocity_in_solid(s);
        }
    }


    //=========================================================================
    // Stops execution when the unavailable explicit CBS path is requested.
    //=========================================================================
    void Steps::rejectExplicitMode(const char* step_name)
    {
        throw std::runtime_error(
            std::string(step_name)
            + " - fully explicit CBS3D mode was requested, but explicit step kernels are not yet ported.");
    }


    //=========================================================================
    // Selects and executes CBS Step 1.
    //=========================================================================
    void Steps::step1(CBSStateSI& s)
    {
        if (s.cfg.cbs_scheme == 1)
        {
            step1SemiImplicit(s);
            return;
        }

        if (s.cfg.cbs_scheme == 0)
        {
            rejectExplicitMode("Steps::step1");
        }

        throw std::runtime_error("Steps::step1 - cbs_scheme must be 1 or 0");
    }


    //=========================================================================
    // Selects and executes CBS Step 2.
    //=========================================================================
    void Steps::step2(CBSStateSI& s)
    {
        if (s.cfg.cbs_scheme == 1)
        {
            step2SemiImplicit(s);
            return;
        }

        if (s.cfg.cbs_scheme == 0)
        {
            rejectExplicitMode("Steps::step2");
        }

        throw std::runtime_error("Steps::step2 - cbs_scheme must be 1 or 0");
    }


    //=========================================================================
    // Selects and executes CBS Step 3.
    //=========================================================================
    void Steps::step3(CBSStateSI& s)
    {
        if (s.cfg.cbs_scheme == 1)
        {
            step3SemiImplicit(s);
            return;
        }

        if (s.cfg.cbs_scheme == 0)
        {
            rejectExplicitMode("Steps::step3");
        }

        throw std::runtime_error("Steps::step3 - cbs_scheme must be 1 or 0");
    }


    //=========================================================================
    // Executes the optional Spalart-Allmaras turbulence transport step.
    //
    // This step is deliberately placed after Step 3.  The SA equation requires
    // the corrected velocity field for both advection and vorticity production.
    // After nu_tilde is updated, the eddy viscosity and effective material
    // properties are refreshed before Step 4 uses the thermal conductivity.
    //=========================================================================
    void Steps::stepSpalartAllmaras(CBSStateSI& s)
    {
        if (s.cfg.turbulence_on < 1)
        {
            return;
        }

        validate_step_dimensions(s);

        SpalartAllmarasAssembly::assembleTransportRhs(s);
        SpalartAllmarasAssembly::updateNuTilde(s);

        TurbulenceBoundary::applyWallValues(s);
        TurbulenceBoundary::applyInletValues(s);

        SpalartAllmarasAssembly::updateEddyViscosity(s);
    }


    //=========================================================================
    // Executes CBS Step 4 when temperature calculation is enabled.
    //=========================================================================
    void Steps::step4(CBSStateSI& s)
    {
        step4Energy(s);
    }


    //=========================================================================
    // CBS Step 1: momentum predictor.
    //=========================================================================
    void Steps::step1SemiImplicit(CBSStateSI& s)
    {
        validate_step_dimensions(s);

        MomentumAssembly::assembleStep1Rhs(s);

        update_velocity_from_rhs_using_predictor_mass(s);

        apply_velocity_boundary_package(s);
    }


    //=========================================================================
    // CBS Step 2: pressure-system assembly and solution.
    //=========================================================================
    void Steps::step2SemiImplicit(CBSStateSI& s)
    {
        validate_step_dimensions(s);

        PressureAssembly::assembleStep2Rhs(s);

        ConjugateGradient::Result result;

        if (s.cfg.solver_opt == 1)
        {
            result = ConjugateGradient::solvePressure(s);
        }
        else if (s.cfg.solver_opt == 2)
        {
            throw std::runtime_error(
                "Steps::step2SemiImplicit - banded pressure solver path is not ported for CBS3D yet");
        }
        else if (s.cfg.solver_opt == 3)
        {
#ifdef CBS3D_USE_PETSC
            result = PetscPressureSolver::solvePressure(s);
#else
            throw std::runtime_error(
                "Steps::step2SemiImplicit - solver_opt=3 requests PETSc, but this executable was built without CBS3D_USE_PETSC");
#endif
        }
        else
        {
            throw std::runtime_error(
                "Steps::step2SemiImplicit - unknown pressure solver option");
        }

        s.last_cg_iterations  = result.iterations;
        s.last_cg_initial_l2  = result.initial_l2;
        s.last_cg_final_l2    = result.final_l2;
        s.last_cg_relative_l2 = result.final_relative_l2;
        s.last_cg_max_abs     = result.final_max_abs;

        if (!result.converged)
        {
            throw std::runtime_error(
                "Steps::step2SemiImplicit - pressure solver failed to converge");
        }

        Boundary::applyPressure(s);
    }


    //=========================================================================
    // CBS Step 3: pressure-gradient velocity correction.
    //=========================================================================
    void Steps::step3SemiImplicit(CBSStateSI& s)
    {
        validate_step_dimensions(s);

        s.rhs.fill(0.0);

        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            if (!is_fluid_element(s, ie))
            {
                continue;
            }

            if (s.detJ(ie) <= 0.0 || !std::isfinite(s.detJ(ie)))
            {
                throw std::runtime_error(
                    "Steps::step3SemiImplicit - invalid detJ at element "
                    + std::to_string(ie));
            }

            Real grad_pres[4] = { 0.0, 0.0, 0.0, 0.0 };

            for (Int idim = 1; idim <= s.cfg.ndim; ++idim)
            {
                for (Int a = 1; a <= s.cfg.nep; ++a)
                {
                    const Int ip = s.intma(a, ie);
                    grad_pres[idim] += grad(s, ie, idim, a) * s.pres(ip);
                }
            }

            const Real vol4 = s.detJ(ie) * s.cfg.fcon[1];

            for (Int a = 1; a <= s.cfg.nep; ++a)
            {
                const Int ip = s.intma(a, ie);

                for (Int idim = 1; idim <= s.cfg.ndim; ++idim)
                {
                    s.rhs(idim, ip) -= vol4 * grad_pres[idim];
                }
            }
        }

        add_velocity_correction_from_rhs(s);

        apply_velocity_boundary_package(s);
    }


    //=========================================================================
    // CBS Step 4: temperature update.
    //=========================================================================
    void Steps::step4Energy(CBSStateSI& s)
    {
        if (s.cfg.temp_calc < 1)
        {
            return;
        }

        EnergyAssembly::assembleStep4Rhs(s);

#pragma omp parallel for schedule(static)
        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            s.temperature(ip) =
                s.temperature1(ip) + s.rhs1(ip) * s.elcoe2p(ip);
        }

        Boundary::applyTemperature(s);
    }
}
