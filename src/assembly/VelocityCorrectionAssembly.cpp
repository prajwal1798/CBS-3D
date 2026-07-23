//=============================================================================
// CBS3D++_SI
//
// Rank-local finite-element assembly of the semi-implicit CBS Step 3
// pressure-gradient velocity correction.
//=============================================================================

#include "cbs/assembly/VelocityCorrectionAssembly.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace cbs
{
    namespace
    {
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
            return s.dNkdx(
                dNkdx_index(s, ie, dim, local_node));
        }


        void validate_inputs(const CBSStateSI& s)
        {
            if (s.cfg.ndim != 3 ||
                s.cfg.nep != 4 ||
                s.cfg.gdim != 13)
            {
                throw std::runtime_error(
                    "VelocityCorrectionAssembly requires the CBS3D "
                    "P1-tetrahedral dimensions ndim=3, nep=4, gdim=13");
            }
        }
    }


    //=========================================================================
    // Assembles the Step 3 pressure-gradient contribution from rank-owned
    // fluid tetrahedra:
    //
    //     grad(p_h)|_e = sum_b p_b grad(N_b)
    //
    //     r_{p,a}^{(e)} = - (V_e / 4) grad(p_h)|_e
    //
    // Since every rank-local tetrahedron is uniquely owned, each element is
    // integrated exactly once globally. Contributions accumulated at ghost
    // nodes are completed later by the reverse halo exchange.
    //=========================================================================
    void VelocityCorrectionAssembly::assembleStep3Rhs(CBSStateSI& s)
    {
        validate_inputs(s);

        s.rhs.fill(0.0);

        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            if (s.mat_elem(ie) != 0)
            {
                continue;
            }

            if (s.detJ(ie) <= 0.0 ||
                !std::isfinite(s.detJ(ie)))
            {
                throw std::runtime_error(
                    "VelocityCorrectionAssembly found invalid detJ at element "
                    + std::to_string(ie));
            }

            Real pressure_gradient[4] =
                {0.0, 0.0, 0.0, 0.0};

            for (Int dim = 1;
                 dim <= s.cfg.ndim;
                 ++dim)
            {
                for (Int a = 1;
                     a <= s.cfg.nep;
                     ++a)
                {
                    const Int ip = s.intma(a, ie);
                    const Real pressure = s.pres(ip);

                    if (!std::isfinite(pressure))
                    {
                        throw std::runtime_error(
                            "VelocityCorrectionAssembly found non-finite "
                            "pressure at local node "
                            + std::to_string(ip));
                    }

                    pressure_gradient[dim] +=
                        grad(s, ie, dim, a) * pressure;
                }
            }

            const Real volume_quarter =
                s.detJ(ie) * s.cfg.fcon[1];

            if (volume_quarter <= 0.0 ||
                !std::isfinite(volume_quarter))
            {
                throw std::runtime_error(
                    "VelocityCorrectionAssembly found invalid V/4 at element "
                    + std::to_string(ie));
            }

            for (Int a = 1;
                 a <= s.cfg.nep;
                 ++a)
            {
                const Int ip = s.intma(a, ie);

                for (Int dim = 1;
                     dim <= s.cfg.ndim;
                     ++dim)
                {
                    s.rhs(dim, ip) -=
                        volume_quarter * pressure_gradient[dim];
                }
            }
        }
    }
}
