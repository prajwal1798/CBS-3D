//=============================================================================
// CBS3D++_SI
//
// Main driver for the three-dimensional semi-implicit
// Characteristic-Based Split finite-element solver.
//
// This file controls the order of the complete numerical calculation:
//
//     1. Read the case files.
//     2. Preprocess the mesh and boundary data.
//     3. Assemble the pressure operator.
//     4. Advance the CBS equations through Steps 1 to 4.
//     5. Evaluate convergence.
//     6. Write residuals and solution files.
//
// The detailed finite-element calculations are performed by the specialised
// modules called from this driver.
//=============================================================================

#include "cbs/solver/Solver.hpp"

#include "cbs/assembly/PressureAssembly.hpp"
#include "cbs/boundary/Boundary.hpp"
#include "cbs/io/MeshIO.hpp"
#include "cbs/io/Post.hpp"
#include "cbs/parallel/Coloring.hpp"
#include "cbs/parallel/HaloExchange.hpp"
#include "cbs/parallel/PartitionMetadata.hpp"
#include "cbs/preprocess/Preprocess.hpp"
#include "cbs/solver/Convergence.hpp"
#include "cbs/solver/Steps.hpp"
#include "cbs/timestep/TimeStep.hpp"
#include "cbs/turbulence/TurbulencePreprocess.hpp"

#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

#ifdef CBS3D_USE_MPI
#include <mpi.h>
#endif

namespace cbs
{
#ifdef CBS3D_USE_MPI
    namespace
    {
        void checkMpi(const int error_code, const char* operation)
        {
            if (error_code != MPI_SUCCESS)
            {
                throw std::runtime_error(
                    std::string("Solver MPI failure in ") + operation);
            }
        }
    }
#endif
    //=========================================================================
    // Stores the case name used to locate the solver input and output files.
    //=========================================================================
    Solver::Solver(std::string case_name)
        : case_name_(std::move(case_name))
    {
    }


    //=========================================================================
    // Returns read-only access to the complete solver state.
    //=========================================================================
    const CBSStateSI& Solver::state() const noexcept
    {
        return s_;
    }


    //=========================================================================
    // Returns modifiable access to the complete solver state.
    //=========================================================================
    CBSStateSI& Solver::state() noexcept
    {
        return s_;
    }


    //=========================================================================
    // Stores the MPI process information in the solver state.
    //
    // The numerical solver remains in serial mode when size is equal to one.
    //=========================================================================
    void Solver::setMpiContext(const Int rank, const Int size) noexcept
    {
        s_.mpi_rank = rank;
        s_.mpi_size = size;
        s_.mpi_enabled = size > 1;
    }




    //=========================================================================
    // Reads one rank-local CBS case and verifies the associated MPI metadata.
    //
    // This is the first distributed-solver integration stage. The routine does
    // not preprocess the finite-element operators and does not advance CBS
    // Steps 1 to 4. It verifies only:
    //
    //     rank-local input files
    //     global/local ownership maps
    //     global partition totals
    //     owner-to-ghost halo communication
    //=========================================================================
    void Solver::runPartitionInitialisation()
    {
        if (!s_.mpi_enabled)
        {
            throw std::runtime_error(
                "Solver::runPartitionInitialisation requires more than one MPI rank");
        }

        MeshIO::readAll(case_name_, s_);
        readPartitionMetadata();
        auditPartitionHalo();
    }



    //=========================================================================
    // Executes the fir
    //
    // Geometric quantities are calculated independently from each rank's owned
    // tetrahedra. Shared nodal mass quantities are reconciled by massMatrix().
    //
    // This milestone intentionally stops before pressure assembly and CBS
    // Steps 1 to 4.
    //=========================================================================
    void Solver::runDistributedPreprocessing()
    {
        if (!s_.mpi_enabled)
        {
            throw std::runtime_error(
                "Solver::runDistributedPreprocessing requires more than one MPI rank");
        }

        MeshIO::readAll(case_name_, s_);
        readPartitionMetadata();

        // Validate partition ownership and the basic forward halo operation.
        auditPartitionHalo();

        // These operations use only owned tetrahedra or rank-local physical
        // boundary faces.
        Preprocess::validateBoundaryFlags(s_);
        Preprocess::shapeFunctionDerivatives(s_);
        Preprocess::assignBoundaryFaceNumbers(s_);
        Preprocess::getNormals(s_);

        // The first numerically meaningful reverse-add and forward halo stage.
        Preprocess::massMatrix(s_);

        // Reconcile the fluid/solid classification of every shared node.
        Preprocess::buildMaterialNodeMasks(s_);

        auditDistributedPreprocessing();
        auditDistributedMaterialMasks();
    }


    //=========================================================================
    // Executes the complete CBS solution procedure.
    //
    // The solver is first initialised and the pressure system is assembled.
    // The main loop then advances the solution by one CBS iteration at a time,
    // writes the requested output and checks the stopping conditions.
    //=========================================================================
    void Solver::run()
    {
        // Print the program title and solver identification.
        Post::printBanner();

        // Read and preprocess the complete problem definition.
        initialise();

        // Assemble the pressure operator before entering the CBS loop.
        preparePressureSystem();

        // Print the complete run configuration and initialise output files.
        Post::printRunSummary(s_, case_name_);
        Post::initialiseRunOutputs(s_, case_name_);

        // Store the final completed iteration and the reason for termination.
        Int last_iteration = 0;
        std::string stop_reason = "maximum iteration count reached";

        Post::printStageDone("Solver loop", "starting CBS iterations");

        // Main CBS iteration loop.
        for (Int iitime = 1; iitime <= s_.cfg.ntime; ++iitime)
        {
            last_iteration = iitime;

            // Reserved location for rereading selected runtime parameters.
            if (s_.cfg.runtime_mod > 0 && (iitime % s_.cfg.iwrite) == 0)
            {
                // Runtime .par reread hook.  Kept visible until the runtime
                // modification grammar is ported.
            }

            // Complete CBS Steps 1 to 4 and evaluate the new residuals.
            advanceOneStep(iitime);

            // Write residuals, console information and VTU solution files.
            {
                auto timer = profiler_.time(SolverProfiler::Section::PostOutput);

                Post::writeResidualRow(s_, case_name_, iitime);
                Post::printProgressLine(s_, iitime, false);

                if (Post::shouldWriteVTU(s_, iitime))
                {
                    Post::writeSolution(s_, case_name_, iitime);
                }
            }

            // Print the timing information for the completed CBS iteration.
            profiler_.printIteration(std::cout, iitime, s_.cfg.console_log_every);

            // Stop a transient calculation when the requested physical time
            // has been reached.
            if (s_.cfg.transient_on > 0 && transientEndTimeReached())
            {
                stop_reason = "transient end time reached";
                break;
            }

            // Stop a steady calculation when all convergence criteria have
            // been satisfied.
            if (s_.cfg.transient_on < 1 && steadyStateReached(iitime))
            {
                stop_reason = "steady-state convergence reached";
                break;
            }
        }

        // Print the final residual information even when the usual console
        // output interval does not coincide with the last iteration.
        Post::printProgressLine(s_, last_iteration, true);

        // Ensure that the final solution is written exactly once.
        if (s_.cfg.vtu_output_enabled > 0)
        {
            const bool already_written =
                !s_.output_iterations.empty() &&
                s_.output_iterations.back() == last_iteration;

            if (!already_written)
            {
                Post::writeSolution(s_, case_name_, last_iteration);
            }
        }

        // Print the final convergence and termination summary.
        Post::printFinalSummary(s_, last_iteration, stop_reason);
    }


    //=========================================================================
    // Reads and preprocesses the complete numerical problem.
    //
    // The routine performs the following operations:
    //
    //     1. Read mesh, boundary-condition, parameter and material files.
    //     2. Validate the boundary-condition identifiers.
    //     3. Calculate tetrahedral shape-function derivatives.
    //     4. Identify boundary faces and calculate outward normals.
    //     5. Assemble the lumped momentum and thermal mass terms.
    //     6. Classify walls, pressure boundaries and special face edges.
    //     7. Build OpenMP element-colouring groups.
    //     8. Apply the initial velocity, pressure and temperature conditions.
//     9. If requested, preprocess the Spalart-Allmaras turbulence model.
    //=========================================================================
    void Solver::initialise()
    {
        // Read the complete problem definition.
        Post::printStage("Input", "reading .plt/.bco/.par/material controls");
        MeshIO::readAll(case_name_, s_);

        if (s_.mpi_enabled)
        {
            readPartitionMetadata();
        }

        Post::printStageDone("Input", "case files loaded");

        // Check that every boundary identifier is supported by the solver.
        Post::printStage("Preprocess", "validating BC flags");
        Preprocess::validateBoundaryFlags(s_);
        Post::printStageDone("Preprocess", "BC flags accepted");

        // Calculate element volumes and Cartesian derivatives of the
        // tetrahedral shape functions.
        Post::printStage("Geometry", "shape derivatives");
        Preprocess::shapeFunctionDerivatives(s_);
        Post::printStageDone("Geometry", "shape derivatives ready");

        // Match every boundary triangle to its tetrahedral face and calculate
        // the corresponding unit normal.
        Post::printStage("Boundary geometry", "face numbers and normals");
        Preprocess::assignBoundaryFaceNumbers(s_);
        Preprocess::getNormals(s_);
        Post::printStageDone("Boundary geometry", "normals ready");

        // Assemble the lumped nodal mass used by the momentum equation and the
        // thermal capacitance used by the energy equation.
        Post::printStage("Mass matrix", "momentum and thermal capacitance");
        Preprocess::massMatrix(s_);
        Post::printStageDone("Mass matrix", "lumped diagonals ready");

        // Identify all nodal and facial boundary groups required by the
        // momentum, pressure and energy equations.
        Post::printStage("Boundary classification", "fedge, walls, pressure nodes");
        Preprocess::classifyFaceEdges(s_);
        Preprocess::elementSize(s_);
        Preprocess::wallDetermination(s_);
        Preprocess::computeMassFlowInletVelocity(s_);
        Preprocess::initialiseVelocityMagnitude(s_);
        Preprocess::detectPressureBoundaryNodes(s_);

        // Divide the elements into independent colour groups so that OpenMP
        // threads can assemble nodal contributions without write conflicts.
        Coloring::build(s_);   // element coloring for race-free OpenMP scatter
        std::cout << "  colors     : " << s_.ncolor << " (parallel scatter groups)\n";
        Post::printStageDone("Boundary classification", "solver boundary state ready");

        // Apply the prescribed initial and boundary values to the solution.
        Post::printStage("Initial boundary values", "velocity, pressure, temperature");
        Boundary::applyTemperature(s_);
        Boundary::applyVelocity(s_);
        Boundary::applyPressure(s_);
        updateVelocityMagnitude();
        Post::printStageDone("Initial boundary values", "initial field constrained");


        // Precompute all geometry-dependent Spalart-Allmaras quantities before
        // the CBS time loop.  In the current milestone this includes the
        // OpenMP wall-distance search and the initial eddy-viscosity field.
        if (s_.cfg.turbulence_on > 0)
        {
            Post::printStage(
                "Turbulence preprocessing",
                "Spalart-Allmaras wall distance");

            TurbulencePreprocess::prepareSpalartAllmaras(s_);

            Post::printStageDone(
                "Turbulence preprocessing",
                "SA wall distance ready");
        }
    }




    //=========================================================================
    // Reads and validates the rank-local partition metadata.
    //=========================================================================
    void Solver::readPartitionMetadata()
    {
        const std::string metadata_path = case_name_ + ".mpi";
        s_.partition_metadata = PartitionMetadataIO::read(metadata_path);

        const PartitionMetadata& metadata = s_.partition_metadata;

        if (metadata.mpi_rank != s_.mpi_rank)
        {
            throw std::runtime_error(
                "Partition MPI rank does not match the executing MPI rank");
        }

        if (metadata.mpi_size != s_.mpi_size)
        {
            throw std::runtime_error(
                "Partition MPI size does not match the executing communicator size");
        }

        if (metadata.partition_id != s_.mpi_rank + 1)
        {
            throw std::runtime_error(
                "Partition ID does not satisfy partition_id = MPI rank + 1");
        }

        if (metadata.local_nelem != s_.cfg.nelem ||
            metadata.local_npoin != s_.cfg.npoin ||
            metadata.local_nboun != s_.cfg.nboun)
        {
            throw std::runtime_error(
                "Partition metadata does not match the rank-local CBS mesh sizes");
        }

        if (metadata.global_nelem > std::numeric_limits<Int>::max() ||
            metadata.global_npoin > std::numeric_limits<Int>::max() ||
            metadata.global_nboun > std::numeric_limits<Int>::max())
        {
            throw std::runtime_error(
                "Global partition size exceeds the CBS Int index range");
        }

        // All tetrahedra stored in a rank-local .plt file are owned by that
        // MPI rank. The exporter does not duplicate ghost tetrahedra.
        s_.owned_elements.resize(static_cast<Size>(s_.cfg.nelem));
        std::iota(s_.owned_elements.begin(), s_.owned_elements.end(), 1);
        s_.ghost_elements.clear();

        // The exporter writes owned nodes first and ghost nodes afterwards.
        s_.owned_nodes.resize(static_cast<Size>(metadata.owned_nodes));
        std::iota(s_.owned_nodes.begin(), s_.owned_nodes.end(), 1);

        s_.ghost_nodes.resize(static_cast<Size>(metadata.ghost_nodes));
        std::iota(
            s_.ghost_nodes.begin(),
            s_.ghost_nodes.end(),
            metadata.owned_nodes + 1);

        s_.node_owner_rank = metadata.node_owner_rank;

        s_.local_to_global_node.assign(
            static_cast<Size>(metadata.local_npoin) + 1U,
            0);

        s_.local_to_global_element.assign(
            static_cast<Size>(metadata.local_nelem) + 1U,
            0);

        s_.global_to_local_node.assign(
            static_cast<Size>(metadata.global_npoin) + 1U,
            -1);

        for (Int local_node = 1;
             local_node <= metadata.local_npoin;
             ++local_node)
        {
            const std::int64_t global_node_64 =
                metadata.local_to_global_node[static_cast<Size>(local_node)];

            if (global_node_64 < 1 ||
                global_node_64 > metadata.global_npoin)
            {
                throw std::runtime_error(
                    "Global node ID is outside the partition range");
            }

            const Int global_node = static_cast<Int>(global_node_64);
            s_.local_to_global_node[static_cast<Size>(local_node)] = global_node;
            s_.global_to_local_node[static_cast<Size>(global_node)] = local_node;

            const Int owner_rank =
                s_.node_owner_rank[static_cast<Size>(local_node)];

            if (owner_rank < 0 || owner_rank >= s_.mpi_size)
            {
                throw std::runtime_error(
                    "Node owner rank is outside the MPI communicator");
            }
        }

        for (Int local_element = 1;
             local_element <= metadata.local_nelem;
             ++local_element)
        {
            const std::int64_t global_element_64 =
                metadata.local_to_global_element[static_cast<Size>(local_element)];

            if (global_element_64 < 1 ||
                global_element_64 > metadata.global_nelem)
            {
                throw std::runtime_error(
                    "Global element ID is outside the partition range");
            }

            s_.local_to_global_element[static_cast<Size>(local_element)] =
                static_cast<Int>(global_element_64);
        }

        s_.element_owner_rank.assign(
            static_cast<Size>(s_.cfg.nelem) + 1U,
            s_.mpi_rank);

        s_.pressure_global_to_local.assign(
            static_cast<Size>(metadata.global_npoin) + 1U,
            -1);
    }


    //=========================================================================
    // Verifies global partition totals and owner-to-ghost communication.
    //=========================================================================
    void Solver::auditPartitionHalo() const
    {
#ifdef CBS3D_USE_MPI
        const PartitionMetadata& metadata = s_.partition_metadata;

        Array1D<Real> global_id_value(s_.cfg.npoin);
        global_id_value.fill(0.0);

        for (const Int local_node : s_.owned_nodes)
        {
            global_id_value(local_node) = static_cast<Real>(
                s_.local_to_global_node[static_cast<Size>(local_node)]);
        }

        HaloExchange::broadcastOwnedToGhosts(
            global_id_value,
            metadata,
            MPI_COMM_WORLD);

        long long local_halo_failures = 0;

        for (Int local_node = 1;
             local_node <= s_.cfg.npoin;
             ++local_node)
        {
            const Real expected = static_cast<Real>(
                s_.local_to_global_node[static_cast<Size>(local_node)]);

            if (global_id_value(local_node) != expected)
            {
                ++local_halo_failures;
            }
        }

        const long long local_owned_elements = s_.cfg.nelem;
        const long long local_owned_nodes = metadata.owned_nodes;
        const long long local_boundary_faces = s_.cfg.nboun;

        long long global_owned_elements = 0;
        long long global_owned_nodes = 0;
        long long global_boundary_faces = 0;
        long long global_halo_failures = 0;

        checkMpi(
            MPI_Allreduce(
                &local_owned_elements,
                &global_owned_elements,
                1,
                MPI_LONG_LONG,
                MPI_SUM,
                MPI_COMM_WORLD),
            "MPI_Allreduce owned elements");

        checkMpi(
            MPI_Allreduce(
                &local_owned_nodes,
                &global_owned_nodes,
                1,
                MPI_LONG_LONG,
                MPI_SUM,
                MPI_COMM_WORLD),
            "MPI_Allreduce owned nodes");

        checkMpi(
            MPI_Allreduce(
                &local_boundary_faces,
                &global_boundary_faces,
                1,
                MPI_LONG_LONG,
                MPI_SUM,
                MPI_COMM_WORLD),
            "MPI_Allreduce boundary faces");

        checkMpi(
            MPI_Allreduce(
                &local_halo_failures,
                &global_halo_failures,
                1,
                MPI_LONG_LONG,
                MPI_SUM,
                MPI_COMM_WORLD),
            "MPI_Allreduce halo failures");

        int local_cells = s_.cfg.nelem;
        int min_cells = 0;
        int max_cells = 0;

        checkMpi(
            MPI_Allreduce(
                &local_cells,
                &min_cells,
                1,
                MPI_INT,
                MPI_MIN,
                MPI_COMM_WORLD),
            "MPI_Allreduce minimum cells");

        checkMpi(
            MPI_Allreduce(
                &local_cells,
                &max_cells,
                1,
                MPI_INT,
                MPI_MAX,
                MPI_COMM_WORLD),
            "MPI_Allreduce maximum cells");

        const bool totals_match =
            global_owned_elements == metadata.global_nelem &&
            global_owned_nodes == metadata.global_npoin &&
            global_boundary_faces == metadata.global_nboun;

        if (!totals_match || global_halo_failures != 0)
        {
            throw std::runtime_error(
                "Integrated MPI partition initialisation failed");
        }

        if (s_.mpi_rank == 0)
        {
            std::cout
                << "============================================================\n"
                << "CBS3D MPI PARTITION INITIALISATION\n"
                << "============================================================\n"
                << "MPI ranks                    : " << s_.mpi_size << "\n"
                << "owned tetrahedra             : " << global_owned_elements << "\n"
                << "unique owned nodes           : " << global_owned_nodes << "\n"
                << "physical boundary triangles  : " << global_boundary_faces << "\n"
                << "cells per rank min/max        : "
                << min_cells << " / " << max_cells << "\n"
                << "halo global-ID exchange      : PASS\n"
                << "RESULT                        : PASS\n"
                << "============================================================\n";
        }
#else
        throw std::runtime_error(
            "Solver::auditPartitionHalo requires an MPI build");
#endif
    }



    //=========================================================================
    // Verifies distributed lumped-mass and thermal-capacitance reconciliation.
    //
    // Each element contributes:
    //
    //     nep * detJ * mass_factor
    //
    // to the total scalar lumped mass. Summing only owned nodal entries avoids
    // counting replicated ghost values.
    //=========================================================================
    void Solver::auditDistributedPreprocessing() const
    {
#ifdef CBS3D_USE_MPI
        Real local_element_mass = 0.0;
        Real local_element_capacity = 0.0;

        for (Int ie = 1; ie <= s_.cfg.nelem; ++ie)
        {
            if (s_.detJ(ie) <= 0.0 ||
                !std::isfinite(s_.detJ(ie)))
            {
                throw std::runtime_error(
                    "Distributed preprocessing audit found invalid detJ");
            }

            const Real element_lumped_mass =
                static_cast<Real>(s_.cfg.nep)
                * s_.detJ(ie)
                * s_.cfg.mass_factor;

            local_element_mass += element_lumped_mass;

            local_element_capacity +=
                s_.rho_cp_e(ie) * element_lumped_mass;
        }

        Real local_owned_mass = 0.0;
        Real local_owned_capacity = 0.0;

        for (const Int ip : s_.owned_nodes)
        {
            if (s_.Mdiag_real(ip) <= 0.0 ||
                !std::isfinite(s_.Mdiag_real(ip)) ||
                s_.elcoe2p(ip) <= 0.0 ||
                !std::isfinite(s_.elcoe2p(ip)))
            {
                throw std::runtime_error(
                    "Distributed preprocessing audit found an invalid owned-node coefficient");
            }

            local_owned_mass +=
                s_.Mdiag_real(ip);

            local_owned_capacity +=
                1.0 / s_.elcoe2p(ip);
        }

        Real global_element_mass = 0.0;
        Real global_element_capacity = 0.0;
        Real global_owned_mass = 0.0;
        Real global_owned_capacity = 0.0;

        checkMpi(
            MPI_Allreduce(
                &local_element_mass,
                &global_element_mass,
                1,
                MPI_DOUBLE,
                MPI_SUM,
                MPI_COMM_WORLD),
            "MPI_Allreduce element lumped mass");

        checkMpi(
            MPI_Allreduce(
                &local_element_capacity,
                &global_element_capacity,
                1,
                MPI_DOUBLE,
                MPI_SUM,
                MPI_COMM_WORLD),
            "MPI_Allreduce element thermal capacity");

        checkMpi(
            MPI_Allreduce(
                &local_owned_mass,
                &global_owned_mass,
                1,
                MPI_DOUBLE,
                MPI_SUM,
                MPI_COMM_WORLD),
            "MPI_Allreduce owned nodal mass");

        checkMpi(
            MPI_Allreduce(
                &local_owned_capacity,
                &global_owned_capacity,
                1,
                MPI_DOUBLE,
                MPI_SUM,
                MPI_COMM_WORLD),
            "MPI_Allreduce owned nodal thermal capacity");

        const auto nearly_equal = [](const Real a, const Real b)
        {
            const Real scale =
                std::fmax(
                    1.0,
                    std::fmax(std::fabs(a), std::fabs(b)));

            return std::fabs(a - b) <= 1.0e-9 * scale;
        };

        if (!nearly_equal(global_owned_mass, global_element_mass) ||
            !nearly_equal(
                global_owned_capacity,
                global_element_capacity))
        {
            throw std::runtime_error(
                "Distributed preprocessing mass reconciliation failed");
        }

        if (s_.mpi_rank == 0)
        {
            std::cout
                << "============================================================\n"
                << "CBS3D DISTRIBUTED PREPROCESSING\n"
                << "============================================================\n"
                << "MPI ranks                    : " << s_.mpi_size << "\n"
                << "element lumped mass          : " << global_element_mass << "\n"
                << "unique owned nodal mass      : " << global_owned_mass << "\n"
                << "element thermal capacity     : " << global_element_capacity << "\n"
                << "owned nodal thermal capacity : " << global_owned_capacity << "\n"
                << "mass reconciliation          : PASS\n"
                << "thermal reconciliation       : PASS\n"
                << "CBS Steps 1 to 4             : NOT STARTED\n"
                << "RESULT                        : PASS\n"
                << "============================================================\n";
        }
#else
        throw std::runtime_error(
            "Solver::auditDistributedPreprocessing requires an MPI build");
#endif
    }


    //=========================================================================
    // Independently verifies the distributed nodal material masks.
    //
    // A global reference mask is constructed directly from the owned
    // tetrahedra using MPI_BOR over global node IDs. This reference does not
    // depend on the neighbour halo maps, so it provides an independent check
    // of the reverse-OR and forward-broadcast implementation.
    //=========================================================================
    void Solver::auditDistributedMaterialMasks() const
    {
#ifdef CBS3D_USE_MPI
        const Int global_npoin =
            static_cast<Int>(s_.partition_metadata.global_npoin);

        if (global_npoin < 1)
        {
            throw std::runtime_error(
                "Distributed material-mask audit found invalid global node count");
        }

        std::vector<Int> local_reference(
            static_cast<Size>(global_npoin) + 1U,
            0);

        std::vector<Int> global_reference(
            static_cast<Size>(global_npoin) + 1U,
            0);

        // Construct the independent rank-local contribution in global-node
        // numbering directly from the owned tetrahedra.
        for (Int ie = 1; ie <= s_.cfg.nelem; ++ie)
        {
            const Int material_bit =
                s_.mat_elem(ie) == 0
                    ? CBSStateSI::node_touches_fluid
                    : CBSStateSI::node_touches_solid;

            for (Int in = 1; in <= s_.cfg.nep; ++in)
            {
                const Int ip = s_.intma(in, ie);

                if (ip < 1 || ip > s_.cfg.npoin)
                {
                    throw std::runtime_error(
                        "Distributed material-mask audit found invalid local node");
                }

                const Size local_index = static_cast<Size>(ip);

                if (local_index >= s_.local_to_global_node.size())
                {
                    throw std::runtime_error(
                        "Distributed material-mask audit found incomplete "
                        "local-to-global node map");
                }

                const Int global_id =
                    s_.local_to_global_node[local_index];

                if (global_id < 1 || global_id > global_npoin)
                {
                    throw std::runtime_error(
                        "Distributed material-mask audit found invalid global node ID");
                }

                local_reference[static_cast<Size>(global_id)] |=
                    material_bit;
            }
        }

        checkMpi(
            MPI_Allreduce(
                local_reference.data(),
                global_reference.data(),
                global_npoin + 1,
                MPI_INT,
                MPI_BOR,
                MPI_COMM_WORLD),
            "MPI_Allreduce global material-mask reference");

        Int local_copy_mismatches = 0;

        for (Int ip = 1; ip <= s_.cfg.npoin; ++ip)
        {
            const Size local_index = static_cast<Size>(ip);

            if (local_index >= s_.local_to_global_node.size())
            {
                ++local_copy_mismatches;
                continue;
            }

            const Int global_id =
                s_.local_to_global_node[local_index];

            if (global_id < 1 || global_id > global_npoin)
            {
                ++local_copy_mismatches;
                continue;
            }

            if (s_.node_material_mask(ip) !=
                global_reference[static_cast<Size>(global_id)])
            {
                ++local_copy_mismatches;
            }
        }

        Int global_copy_mismatches = 0;

        checkMpi(
            MPI_Allreduce(
                &local_copy_mismatches,
                &global_copy_mismatches,
                1,
                MPI_INT,
                MPI_SUM,
                MPI_COMM_WORLD),
            "MPI_Allreduce material-mask mismatches");

        Int local_owned_counts[4] = {0, 0, 0, 0};

        for (const Int ip : s_.owned_nodes)
        {
            const Int mask = s_.node_material_mask(ip);

            if (mask >= 0 && mask <= 3)
            {
                ++local_owned_counts[mask];
            }
            else
            {
                ++local_owned_counts[0];
            }
        }

        Int global_owned_counts[4] = {0, 0, 0, 0};

        checkMpi(
            MPI_Allreduce(
                local_owned_counts,
                global_owned_counts,
                4,
                MPI_INT,
                MPI_SUM,
                MPI_COMM_WORLD),
            "MPI_Allreduce owned material-mask counts");

        Int reference_counts[4] = {0, 0, 0, 0};

        for (Int global_id = 1;
             global_id <= global_npoin;
             ++global_id)
        {
            const Int mask =
                global_reference[static_cast<Size>(global_id)];

            if (mask >= 0 && mask <= 3)
            {
                ++reference_counts[mask];
            }
            else
            {
                ++reference_counts[0];
            }
        }

        const Int global_owned_total =
            global_owned_counts[1]
            + global_owned_counts[2]
            + global_owned_counts[3];

        const bool counts_match =
            global_owned_counts[0] == 0
            && reference_counts[0] == 0
            && global_owned_counts[1] == reference_counts[1]
            && global_owned_counts[2] == reference_counts[2]
            && global_owned_counts[3] == reference_counts[3]
            && global_owned_total == global_npoin;

        if (global_copy_mismatches != 0 || !counts_match)
        {
            throw std::runtime_error(
                "Distributed material-mask reconciliation failed");
        }

        if (s_.mpi_rank == 0)
        {
            std::cout
                << "============================================================\n"
                << "CBS3D DISTRIBUTED MATERIAL MASKS\n"
                << "============================================================\n"
                << "MPI ranks                    : " << s_.mpi_size << "\n"
                << "fluid-only owned nodes       : "
                << global_owned_counts[1] << "\n"
                << "solid-only owned nodes       : "
                << global_owned_counts[2] << "\n"
                << "fluid-solid interface nodes  : "
                << global_owned_counts[3] << "\n"
                << "unique owned nodes           : "
                << global_owned_total << "\n"
                << "invalid material masks       : "
                << global_owned_counts[0] << "\n"
                << "owner/ghost mask agreement   : PASS\n"
                << "independent MPI_BOR reference: PASS\n"
                << "CBS Steps 1 to 4             : NOT STARTED\n"
                << "RESULT                        : PASS\n"
                << "============================================================\n";
        }
#else
        throw std::runtime_error(
            "Solver::auditDistributedMaterialMasks requires an MPI build");
#endif
    }


    //=========================================================================
    // Builds the pressure operator required by CBS Step 2.
    //
    // Element pressure terms are calculated first and then assembled into the
    // global pressure system used by the pressure-correction solver.
    //=========================================================================
    void Solver::preparePressureSystem()
    {
        Post::printStage("Pressure operator", "building element/global pressure terms");
        PressureAssembly::buildElementPressureTerms(s_);
        PressureAssembly::buildGlobalPressureTerms(s_);
        Post::printStageDone("Pressure operator", "pressure system ready");
    }


    //=========================================================================
    // Advances the solution through one complete CBS iteration.
    //
    // The calculation follows the standard semi-implicit CBS sequence:
    //
    //     Step 1  Calculate the intermediate velocity increment.
    //     Step 2  Solve the pressure-correction equation.
    //     Step 3  Correct the velocity field.
    //     Step 4  Advance the temperature field when energy is enabled.
    //
    // The previous solution is stored before the new iteration begins so that
    // residuals can be evaluated after all CBS steps are complete.
    //=========================================================================
    void Solver::advanceOneStep(const Int iitime)
    {
        profiler_.resetIteration();

        // Store the previous iteration values for residual evaluation.
        s_.unkn1 = s_.unkno;
        s_.pres1 = s_.pres;
        s_.temperature1 = s_.temperature;

        // Store the SA working variable at the beginning of the current
        // iteration.  The SA residual assembly is explicit in nu_tilde and
        // therefore uses nu_tilde1 as q^n.
        if (s_.cfg.turbulence_on > 0)
        {
            s_.nu_tilde1 = s_.nu_tilde;
        }

        // Calculate the physical or pseudo-time step used by this iteration.
        {
            auto timer = profiler_.time(SolverProfiler::Section::TimeStepCompute);
            TimeStep::computeTimeStep(s_, iitime);
        }

        // Advance the accumulated physical time for transient calculations.
        if (s_.cfg.transient_on > 0)
        {
            s_.cfg.rtime += s_.cfg.dtreal;
        }

        // Update the diagonal left-hand-side coefficients that depend on the
        // current time-step size.
        {
            auto timer = profiler_.time(SolverProfiler::Section::LhsDiagonal);
            TimeStep::updateLhsDiagonal(s_);
        }

        // CBS Step 1: calculate the intermediate momentum contribution.
        {
            auto timer = profiler_.time(SolverProfiler::Section::Step1Momentum);
            Steps::step1(s_);
        }

        // Optionally correct the time step using the Step 2 pressure condition.
        if (s_.cfg.step2_check > 0)
        {
            auto timer = profiler_.time(SolverProfiler::Section::Step2TimeStepCorrection);
            TimeStep::applyStep2PressureTimeStepCorrection(s_);
            TimeStep::updateLhsDiagonal(s_);
        }

        // CBS Step 2: assemble and solve the pressure-correction equation.
        {
            auto timer = profiler_.time(SolverProfiler::Section::Step2PressureRhs);
            Steps::step2(s_);
        }

        // CBS Step 3: correct the nodal velocity using the pressure increment.
        {
            auto timer = profiler_.time(SolverProfiler::Section::Step3VelocityCorrection);
            Steps::step3(s_);
        }

        // Optional Spalart-Allmaras transport step.  It is placed after
        // Step 3 because the SA advection and production terms use the corrected
        // velocity field.  It is placed before Step 4 so that the energy equation
        // can use the updated turbulent thermal conductivity.
        if (s_.cfg.turbulence_on > 0)
        {
            auto timer = profiler_.time(SolverProfiler::Section::StepSpalartAllmaras);
            Steps::stepSpalartAllmaras(s_);
        }

        // CBS Step 4: advance the energy equation when temperature is enabled.
        {
            auto timer = profiler_.time(SolverProfiler::Section::Step4Energy);
            Steps::step4(s_);
        }

        // Update the scalar velocity magnitude used by output and monitoring.
        {
            auto timer = profiler_.time(SolverProfiler::Section::VelocityMagnitude);
            updateVelocityMagnitude();
        }

        // Calculate the relative and absolute residual measures.
        {
            auto timer = profiler_.time(SolverProfiler::Section::Convergence);
            Convergence::evaluate(s_);
        }
    }


    //=========================================================================
    // Calculates the magnitude of the three-dimensional nodal velocity.
    //
    // A very small positive value is included inside the square root to avoid
    // an exactly zero argument in later calculations that may use velocity.
    //=========================================================================
    void Solver::updateVelocityMagnitude()
    {
        for (Int ip = 1; ip <= s_.cfg.npoin; ++ip)
        {
            const Real u = s_.unkno(1, ip);
            const Real v = s_.unkno(2, ip);
            const Real w = s_.unkno(3, ip);

            s_.velocity(ip) = std::sqrt(u * u + v * v + w * w + 1.0e-16);
        }
    }


    //=========================================================================
    // Checks whether the steady-state solution has converged.
    //
    // Convergence is not accepted until the specified minimum number of CBS
    // iterations has been completed.
    //=========================================================================
    bool Solver::steadyStateReached(const Int iitime) const
    {
        if (iitime < s_.cfg.steady_min_iterations)
        {
            return false;
        }

        return Convergence::steadyStateReached(s_);
    }


    //=========================================================================
    // Checks whether the requested physical end time has been reached.
    //=========================================================================
    bool Solver::transientEndTimeReached() const
    {
        return s_.cfg.end_rtime > 0.0 &&
               s_.cfg.rtime >= s_.cfg.end_rtime;
    }
}
