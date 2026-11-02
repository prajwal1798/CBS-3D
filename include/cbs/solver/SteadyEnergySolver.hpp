#pragma once

//=============================================================================
// CBS3D++_SI
//
// FILE      : SteadyEnergySolver.hpp
// PURPOSE   : IMPLICIT STEADY CONJUGATE ENERGY SOLVE FOR A FROZEN VELOCITY.
//
//             THE EXPLICIT STEP 4 MARCHES  T <- T + C^-1 R(T)  UNTIL R = 0,
//             WHERE R(T) = M T + c + q IS AFFINE IN T FOR A FIXED VELOCITY:
//               M  ELEMENT OPERATOR (CONVECTION + div(u) CORRECTION +
//                  CHARACTERISTIC STABILISATION IN THE FLUID, CONDUCTION
//                  EVERYWHERE)          -> EnergyAssembly::elementOperator
//               c  T_ref SHIFT AND VOLUMETRIC SOURCES (R AT T = 0)
//               q  PRESCRIBED WALL HEAT FLUX (BC 532)
//             THIS SOLVER FINDS THE SAME STEADY STATE DIRECTLY:
//
//                 -M T = c + q   (FREE NODES),   T = T_D   (DIRICHLET NODES)
//
//             ONE DISTRIBUTED PETSc SOLVE (GMRES + ASM/ILU(1) BY DEFAULT;
//             CLASSICAL BOOMERAMG BREAKS DOWN ON THE ADVECTION OPERATOR),
//             OPTIONS PREFIX energy_, E.G.
//             PETSC_OPTIONS="-energy_ksp_monitor -energy_pc_type hypre".
//             CONSTANT PROPERTIES AND NO BUOYANCY: THE FLOW DOES NOT DEPEND
//             ON T, SO THIS IS EXACT FOR THE GIVEN VELOCITY FIELD.
//
//             WITH CBS3D_TEMPERATURE_AFC=1 THE NONLINEAR DEFECT-CORRECTION
//             PATH USES THE TRANSIENT ZALESAK LIMITER; THIS METHOD FAILED
//             TO CONVERGE ON THE ORIGINAL U-CHANNEL (G12B/G12C) AND IS NOT
//             A VALIDATED IMPLICIT HIGH-ORDER CHT SOLVER.
//
//             CBS3D_STEADY_ENERGY_LOW_ORDER_ONLY=1 REQUIRES AFC=1 ONLY
//             TO ENABLE FLUID GRAPH VISCOSITY, AND SOLVES A_LOW*T=b
//             WITHOUT NONLINEAR ANTIDIFFUSIVE CORRECTION. G12D CASE PASS,
//             BUT IT IS A MORE DIFFUSIVE EXPERIMENTAL BASELINE.
//
//             CBS3D_IMPLICIT_ENERGY_BE_DT > 0 [s] AND
//             CBS3D_IMPLICIT_ENERGY_BE_STEPS=N ENABLE *PHYSICAL-TIME*
//             BACKWARD EULER FOR THE FROZEN VELOCITY USING LOW_ORDER_ONLY.
//             THE LUMPED CAPACITY IS ASSEMBLED FROM RHO_CP*V/4, NOT FROM
//             LOCAL CBS PSEUDO-TIME. N PHYSICAL STEPS REUSE ONE PETSc KSP
//             MATRIX. STEPWISE TIME REPORTS ARE RELATIVE TO THE IMPORTED
//             INITIAL TEMPERATURE STATE. NO TEMPERATURE-FLOW FEEDBACK OR
//             HIGH-ORDER BOUNDED BE IS IMPLEMENTED.
//
//             HIGH-ORDER AND UNLIMITED STEADY MODES CHECK THE ASSEMBLED
//             ENERGY RESIDUAL. LOW-ORDER STEADY / PHYSICAL BE MODES CHECK
//             THEIR OWN PETSC MATRIX RESIDUAL. THESE ARE DISTINCT CHECKS.
//
// USAGE     : CBS3D_STEADY_ENERGY=1 IN THE MPI PRODUCTION RUN (SEE
//             SolverDistributedProductionLoop.cpp); THE VELOCITY COMES FROM A
//             RESTART (CBS3D_RESTART_VELOCITY_ARRAY=mean_velocity FOR A
//             TIME-AVERAGED FIELD).
//=============================================================================

#include "cbs/core/Array.hpp"
#include "cbs/core/CBSStateSI.hpp"

namespace cbs
{
    class SteadyEnergySolver
    {
    public:
        struct Result
        {
            Int iterations = 0;
            Real ksp_relative_residual = 0.0;
            Real residual_check = 0.0;     // max |R_AFC(T)| / max |c + q| OVER FREE NODES (AFC OR UNLIMITED)
            Real seconds = 0.0;
            bool converged = false;
            Real t_min = 0.0;
            Real t_max = 0.0;
            bool low_order = false;        // experimental graph-viscosity FEM
            bool transient_be = false;      // true physical-time backward Euler
            Int be_steps = 0;
            Real be_dt = 0.0;               // physical global timestep in seconds
            bool afc = false;              // true only for nonlinear flux-corrected solution
            Int afc_iterations = 0;        // DEFECT-CORRECTION ITERATIONS
            Real afc_change = 0.0;         // LAST max |T^{k+1} - T^k| [K]
        };

        // fixed / value / heat_flux_load: THE THERMAL BOUNDARY STATE OF THE
        // PRODUCTION LOOP (DIRICHLET MASK AND VALUES, NODAL BC 532 LOAD).
        static Result solve(CBSStateSI& s,
                            const Array1D<Int>& fixed,
                            const Array1D<Real>& value,
                            const Array1D<Real>& heat_flux_load);
    };
}
