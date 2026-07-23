#pragma once

//=============================================================================
// CBS3D++_SI
//
// Main driver for the three-dimensional semi-implicit
// Characteristic-Based Split finite-element solver.
//
// The Solver class controls the complete calculation. It reads the problem
// data, prepares the finite-element quantities, advances the CBS equations,
// checks convergence and writes the numerical results.
//
// The detailed numerical calculations are performed by the specialised
// preprocessing, assembly, boundary-condition, time-step and output modules.
//=============================================================================

#include "cbs/core/Types.hpp"
#include "cbs/core/CBSStateSI.hpp"
#include "cbs/utils/SolverProfiler.hpp"

#include <string>

namespace cbs
{
    class Solver
    {
    public:
        // Creates a solver for the selected case.
        //
        // The case name is used as the common base name of the input and
        // output files.
        explicit Solver(std::string case_name);

        // Executes the complete CBS solution procedure.
        void run();

        // Reads and validates one rank-local partition without advancing the
        // CBS equations. This is the first distributed-MPI integration stage.
        void runPartitionInitialisation();

        // Executes rank-local geometric preprocessing and reconciles shared
        // nodal mass quantities across MPI partition interfaces. CBS Steps
        // 1 to 4 and the pressure solve are not advanced by this milestone.
        void runDistributedPreprocessing();

        // Preprocesses the distributed mesh and executes one complete
        // owner/ghost-consistent CBS momentum-predictor step.
        void runDistributedStep1();

        // Executes distributed Step 1 followed by the MPI/PETSc pressure step.
        void runDistributedStep2();

        // Executes distributed Steps 1 to 3 and measures the post-correction
        // weak-divergence residual on unconstrained pressure rows.
        void runDistributedStep3();

        // Returns read-only access to the complete solver state.
        [[nodiscard]] const CBSStateSI& state() const noexcept;

        // Returns modifiable access to the complete solver state.
        [[nodiscard]] CBSStateSI& state() noexcept;

        // Stores the MPI rank and number of MPI processes in the solver state.
        //
        // This method already exists in Solver.cpp and is declared here so
        // that the header and source file remain consistent.
        void setMpiContext(Int rank, Int size) noexcept;

    private:
        // Base name of the problem input and output files.
        std::string case_name_;

        // Complete mesh, physical parameters, solution variables and
        // numerical work arrays used by the CBS solver.
        CBSStateSI s_;

        // Measures the computational time spent in the main solver stages.
        SolverProfiler profiler_;

        // Reads the input files and prepares the mesh, geometry, boundary
        // conditions and initial solution fields.
        void initialise();

        // Reads the rank-local .mpi file and builds ownership/global maps.
        void readPartitionMetadata();

        // Verifies the integrated owner-to-ghost halo communication.
        void auditPartitionHalo() const;

        // Verifies that distributed nodal mass and thermal capacitance agree
        // with independently integrated rank-local element contributions.
        void auditDistributedPreprocessing() const;

        // Verifies the globally reconciled fluid/solid nodal masks.
        void auditDistributedMaterialMasks() const;

        // Independently verifies distributed wall/interface classification and
        // area-weighted wall-normal reconciliation.
        void auditDistributedWallClassification() const;

        // Independently reconstructs distributed prescribed-pressure state and
        // validates owner/ghost flags and each rank-local pressure list.
        void auditDistributedPressureBoundary() const;

        // Inventories global velocity-relevant boundary classifications,
        // overlaps and inlet-normal multiplicity without modifying solver state.
        void auditDistributedVelocityBoundaryInventory() const;

        // Independently validates the persistent distributed nodal velocity-
        // boundary types, priorities, values, topology flags and inlet normals.
        void auditDistributedVelocityBoundaryState() const;

        // Builds the element and global pressure operators used in CBS Step 2.
        void preparePressureSystem();

        // Advances the solution through one complete CBS iteration.
        void advanceOneStep(Int iitime);

        // Calculates the velocity magnitude at every mesh node.
        void updateVelocityMagnitude();

        // Checks the steady-state convergence criterion after the prescribed
        // minimum number of CBS iterations.
        [[nodiscard]] bool steadyStateReached(Int iitime) const;

        // Checks whether the requested transient end time has been reached.
        [[nodiscard]] bool transientEndTimeReached() const;
    };
}
