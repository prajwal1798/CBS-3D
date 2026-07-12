#include "cbs/turbulence/TurbulencePreprocess.hpp"

#include "cbs/assembly/SpalartAllmarasAssembly.hpp"
#include "cbs/boundary/TurbulenceBoundary.hpp"
#include "cbs/turbulence/WallDistance.hpp"

#include <stdexcept>

namespace cbs
{
    //=========================================================================
    // Prepares the Spalart-Allmaras turbulence model before time advancement.
    //
    // This routine performs only preprocessing and initialisation.  It does not
    // assemble the SA transport equation and it does not modify the CBS Step 1,
    // Step 2, Step 3 or Step 4 algorithms.
    //
    // Operations performed when turbulence_on is enabled:
    //
    //     1. classify SA-active, SA-wall and SA-inlet nodes;
    //     2. compute wall_distance(node) using the threaded wall-distance search;
    //     3. initialise the transported SA working variable nu_tilde;
    //     4. compute the first eddy-viscosity and effective-property fields.
    //
    // The wall-distance calculation in this milestone is deliberately not an MPI
    // algorithm.  It assumes that the current process has the complete wall-face
    // list required to compute the true nearest-wall distance.  The distributed
    // version should later read precomputed distances from the partition files or
    // use a proper distributed nearest-wall search.
    //=========================================================================
    void TurbulencePreprocess::prepareSpalartAllmaras(CBSStateSI& s)
    {
        SpalartAllmarasAssembly::resetEffectiveProperties(s);

        if (s.cfg.turbulence_on < 1)
        {
            return;
        }

        if (s.mpi_enabled)
        {
            throw std::runtime_error(
                "Spalart-Allmaras wall-distance preprocessing is currently "
                "implemented for the complete serial/OpenMP mesh only");
        }

        TurbulenceBoundary::classifyNodes(s);

        turbulence::WallDistance::compute(s);

        TurbulenceBoundary::initialiseNuTilde(s);

        SpalartAllmarasAssembly::updateEddyViscosity(s);
    }
}
