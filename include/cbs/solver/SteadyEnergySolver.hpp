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
//             WITH CBS3D_TEMPERATURE_AFC=1 THE SAME BOUNDED (NO UNDER/
//             OVERSHOOT) SCHEME IS SOLVED: LOW-ORDER MATRIX (GRAPH VISCOSITY
//             d_ab IN FLUID ELEMENTS) PLUS ZALESAK-LIMITED ANTIDIFFUSIVE
//             FLUXES, BY DEFECT CORRECTION  A_L T^{k+1} = b + g(T^k)  WITH ONE
//             PRECONDITIONER FOR ALL ITERATIONS (CBS3D_STEADY_ENERGY_TOL [K],
//             _MAXIT, _RELAX).
//
//             AFTER THE SOLVE THE FULL EXPLICIT RESIDUAL R(T) IS ASSEMBLED
//             AGAIN WITH THE ORIGINAL CODE AND REPORTED (MUST BE ~ SOLVER
//             TOLERANCE): AN INDEPENDENT CHECK OF THE MATRIX.
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
            bool afc = false;              // CBS3D_TEMPERATURE_AFC=1: LIMITED SOLUTION
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
