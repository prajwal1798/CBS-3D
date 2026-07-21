#pragma once

//=============================================================================
// CBS3D++_SI
//
// Finite-element preprocessing routines for the three-dimensional
// semi-implicit Characteristic-Based Split solver.
//
// The routines in this module prepare the tetrahedral mesh and nodal data
// required before the CBS iterations begin. They:
//
//     1. Validate the mesh dimensions and boundary identifiers.
//     2. Calculate P1 tetrahedral shape-function gradients and element volumes.
//     3. Match each boundary triangle to its parent tetrahedral face.
//     4. Calculate area-weighted outward face normals.
//     5. Assemble lumped momentum mass and thermal capacitance terms.
//     6. Classify boundary faces and determine an element length scale.
//     7. Build wall, interface and pressure-boundary node lists.
//     8. Convert a prescribed inlet mass-flow rate into velocity magnitude.
//
// Current CHT boundary identifiers:
//
//     511  mass-flow inlet
//     520  pressure outlet
//     530  no-slip adiabatic wall
//     532  no-slip prescribed heat-flux wall
//     901  fluid-solid interface
//
// Boundary identifier 902 is retained only as a legacy heat-flux marker.
// It must not be interpreted as a fluid-solid interface.
//=============================================================================

#include "cbs/core/CBSStateSI.hpp"

namespace cbs
{
    class Preprocess
    {
    public:
        // Check mesh dimensions and report unsupported boundary identifiers.
        static void validateBoundaryFlags(CBSStateSI& s);

        // Calculate constant P1 shape-function gradients and det(J) for every
        // tetrahedral element.
        static void shapeFunctionDerivatives(CBSStateSI& s);

        // Match every boundary triangle to the corresponding local face of its
        // parent tetrahedral element.
        static void assignBoundaryFaceNumbers(CBSStateSI& s);

        // Calculate area-weighted outward normals and areas for all
        // tetrahedral and boundary faces.
        static void getNormals(CBSStateSI& s);

        // Assemble nodal lumped mass, thermal capacitance and the
        // lumped-minus-consistent mass correction.
        static void massMatrix(CBSStateSI& s);

        // Build and reconcile each node's fluid/solid connectivity mask.
        static void buildMaterialNodeMasks(CBSStateSI& s);

        // Mark exterior faces as prescribed-velocity or other boundary faces.
        static void classifyFaceEdges(CBSStateSI& s);

        // Calculate the minimum tetrahedral altitude used as the element
        // characteristic length.
        static void elementSize(CBSStateSI& s);

        // Build the no-slip wall and conformal fluid-solid interface node list.
        static void wallDetermination(CBSStateSI& s);

        // Convert the prescribed inlet mass-flow rate into velocity magnitude.
        static void computeMassFlowInletVelocity(CBSStateSI& s);

        // Initialise the nodal velocity magnitude and its previous value.
        static void initialiseVelocityMagnitude(CBSStateSI& s);

        // Build the pressure-outlet node list or select one reference node.
        static void detectPressureBoundaryNodes(CBSStateSI& s);
    };
}
