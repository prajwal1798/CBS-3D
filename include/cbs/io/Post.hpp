#pragma once

//=============================================================================
// CBS3D++_SI
//
// Post-processing, terminal monitoring and solution-output control.
//
// The module provides:
//
//     1. Program banner and run summary.
//     2. Stage-by-stage terminal messages.
//     3. Per-iteration residual CSV output.
//     4. VTU solution files for ParaView.
//     5. A PVD collection file for time-dependent visualisation.
//     6. Console progress, elapsed time and estimated time remaining.
//     7. Optional launch of the live residual plotting script.
//
// The solver stores mesh indices using one-based numbering. VTU connectivity is
// written with zero-based numbering, as required by the VTK XML format.
//
// The output routines preserve the physical domains of the CHT problem:
//
//     pressure     written only on fluid-connected nodes
//     velocity     written only on fluid-only flow nodes
//     temperature  written on the complete fluid-solid thermal domain
//
// Boundary triangles are included as separate VTK cells so that material and
// boundary-condition identifiers can be inspected directly in ParaView.
//=============================================================================

#include "cbs/core/Types.hpp"
#include "cbs/core/CBSStateSI.hpp"

#include <string>

namespace cbs
{
    class Post
    {
    public:
        // Prints the program title and solver capabilities.
        static void printBanner();

        // Prints the mesh, numerical controls and output configuration.
        static void printRunSummary(
            const CBSStateSI& s,
            const std::string& case_name);

        // Prints the start of one preprocessing or solution stage.
        static void printStage(
            const std::string& stage,
            const std::string& detail = std::string{});

        // Prints successful completion of one stage.
        static void printStageDone(
            const std::string& stage,
            const std::string& detail = std::string{});

        // Creates output directories, residual files and the initial solution
        // output before the CBS iteration loop begins.
        static void initialiseRunOutputs(
            CBSStateSI& s,
            const std::string& case_name);

        // Appends one CBS-iteration row to the residual CSV file.
        static void writeResidualRow(
            const CBSStateSI& s,
            const std::string& case_name,
            Int iitime);

        // Returns true when a VTU file is required at the current iteration.
        static bool shouldWriteVTU(
            CBSStateSI& s,
            Int iitime);

        // Writes one VTU solution and updates the PVD time collection.
        static void writeSolution(
            CBSStateSI& s,
            const std::string& case_name,
            Int iitime);

        // Writes one VTK XML unstructured-grid solution file.
        static void writeVTU(
            const CBSStateSI& s,
            const std::string& case_name,
            Int iitime);

        // Writes the ParaView time-series collection file.
        static void writePVD(
            const CBSStateSI& s,
            const std::string& case_name);

        // Prints the current CBS iteration, residuals and timing information.
        static void printProgressLine(
            const CBSStateSI& s,
            Int iitime,
            bool force_newline);

        // Prints the final stopping reason and convergence summary.
        static void printFinalSummary(
            const CBSStateSI& s,
            Int last_iteration,
            const std::string& stop_reason);

    private:
        // Launches the optional external Python residual monitor.
        static void launchLiveResidualPlotterIfRequested(
            const CBSStateSI& s,
            const std::string& case_name);
    };
}
