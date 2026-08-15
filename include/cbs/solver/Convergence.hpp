#pragma once

//=============================================================================
// CBS3D++_SI
//
// Convergence measures for the CBS iteration loop.
//
// For a nodal variable phi, the relative change between two CBS iterations is
//
//     R_phi =
//         sqrt[ sum_i (phi_i^(n+1) - phi_i^n)^2
//               /
//               (sum_i (phi_i^(n+1))^2 + epsilon) ]
//
// The current-field L2 norm is
//
//     ||phi^(n+1)||_2 = sqrt[ sum_i (phi_i^(n+1))^2 ]
//
// For velocity and temperature, an additional residual is reconstructed on
// the assembled RHS scale by dividing the nodal update by the inverse diagonal:
//
//     r_u,i = (u_i^(n+1) - u_i^n) / elcoe2(i)
//
//     r_T,i = (T_i^(n+1) - T_i^n) / elcoe2p(i)
//
// Material-domain treatment:
//
//     velocity     fluid-only nodes that do not touch the solid
//     pressure     all fluid-connected nodes, including interface nodes
//     temperature  complete fluid-solid thermal domain
//
// Residual storage in CBSStateSI::hb:
//
//     hb[0]   relative residual of u
//     hb[1]   L2 norm of u
//     hb[2]   RHS-scale residual of u
//
//     hb[3]   relative residual of v
//     hb[4]   L2 norm of v
//     hb[5]   RHS-scale residual of v
//
//     hb[6]   relative residual of w
//     hb[7]   L2 norm of w
//     hb[8]   RHS-scale residual of w
//
//     hb[9]   relative pressure residual
//     hb[10]  L2 norm of pressure
//     hb[11]  absolute pressure-update norm
//
//     hb[12]  relative temperature residual
//     hb[13]  L2 norm of temperature
//     hb[14]  RHS-scale temperature residual
//
// CBSStateSI::hb must therefore contain 15 entries.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"

namespace cbs
{
    class Convergence
    {
    public:
        // Calculates all velocity, pressure and temperature convergence
        // measures for the current CBS iteration.
        static void evaluate(CBSStateSI& s);

        // Returns true when every enabled steady-state criterion is satisfied.
        static bool steadyStateReached(const CBSStateSI& s);

        // Returns the largest relative residual among u, v and w.
        static Real velocityResidual(const CBSStateSI& s);

        // Returns the relative pressure residual.
        static Real pressureResidual(const CBSStateSI& s);

        // Returns the relative temperature residual.
        static Real temperatureResidual(const CBSStateSI& s);

        // Returns the relative L2 residual of the Spalart-Allmaras working
        // variable:
        //
        //     R_nu = ||nu_tilde^(n+1) - nu_tilde^n|| / max(||nu_tilde^(n+1)||, eps)
        //
        // Under MPI the two sums are accumulated over owned nodes only and then
        // reduced with MPI_Allreduce, so the reported value is the global
        // residual and does not change with the number of partitions.  Summing
        // over all local nodes instead would count every interface node once per
        // rank that holds a copy.
        static Real turbulenceResidual(const CBSStateSI& s);

        // Global extrema of the Spalart-Allmaras fields.
        //
        // All values are reduced with MPI_Allreduce over owned nodes and local
        // elements, so they are independent of the partition count and can be
        // compared directly between a serial and a distributed run.  They exist
        // to make a diverging SA field visible in the iteration monitor while it
        // is still growing, rather than after it has produced a non-finite
        // effective viscosity.
        struct TurbulenceDiagnostics
        {
            Real nu_tilde_min = 0.0;
            Real nu_tilde_max = 0.0;
            Real nu_t_min = 0.0;
            Real nu_t_max = 0.0;
            Real mu_t_max = 0.0;
            Real mu_eff_max = 0.0;
            Real chi_max = 0.0;
            Real residual = 0.0;
        };

        static TurbulenceDiagnostics turbulenceDiagnostics(const CBSStateSI& s);
    };
}
