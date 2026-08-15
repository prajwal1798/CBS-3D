//=============================================================================
// CBS3D++_SI
//
// Distributed-memory residual and ParaView output implementation.
//=============================================================================

#include "cbs/io/DistributedPost.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

#ifdef CBS3D_USE_MPI
#include <mpi.h>
#endif

namespace cbs
{
#if defined(CBS3D_USE_MPI)
    namespace
    {
        namespace fs = std::filesystem;

        void check_mpi(const int error_code, const char* operation)
        {
            if (error_code != MPI_SUCCESS)
            {
                throw std::runtime_error(
                    std::string("DistributedPost MPI failure in ") + operation);
            }
        }

        Real safe_value(const Real value)
        {
            return std::isfinite(value) ? value : 0.0;
        }

        std::string step_tag(const Int iteration)
        {
            std::ostringstream text;
            text << "step_" << std::setw(8) << std::setfill('0') << iteration;
            return text.str();
        }

        std::string rank_tag(const Int rank)
        {
            std::ostringstream text;
            text << "rank_" << std::setw(4) << std::setfill('0') << rank;
            return text.str();
        }

        std::string case_name_from_local_path(
            const std::string& rank_local_case_name)
        {
            const fs::path local_path(rank_local_case_name);
            std::string name = local_path.filename().string();

            const std::string marker = "_rank_";
            const std::size_t marker_position = name.rfind(marker);

            if (marker_position != std::string::npos)
            {
                const std::string suffix =
                    name.substr(marker_position + marker.size());

                const bool numeric_suffix =
                    !suffix.empty() &&
                    std::all_of(
                        suffix.begin(),
                        suffix.end(),
                        [](const unsigned char character)
                        {
                            return character >= '0' && character <= '9';
                        });

                if (numeric_suffix)
                {
                    name.resize(marker_position);
                }
            }

            if (name.empty())
            {
                throw std::runtime_error(
                    "DistributedPost obtained an empty distributed case name");
            }

            return name;
        }

        fs::path case_output_directory(const std::string& case_name)
        {
            return fs::path("output") / case_name;
        }

        std::string piece_file_name(
            const std::string& case_name,
            const Int iteration,
            const Int rank)
        {
            return case_name + "_" + step_tag(iteration)
                 + "_" + rank_tag(rank) + ".vtu";
        }

        std::string parallel_file_name(
            const std::string& case_name,
            const Int iteration)
        {
            return case_name + "_" + step_tag(iteration) + ".pvtu";
        }

        fs::path residual_file_path(const std::string& case_name)
        {
            return case_output_directory(case_name) /
                (case_name + "_distributed_residuals.csv");
        }

        fs::path pvd_file_path(const std::string& case_name)
        {
            return case_output_directory(case_name) /
                (case_name + "_distributed.pvd");
        }

        // Molecular dynamic viscosity used to normalise mu_t at a node.
        //
        // There is no node-to-element adjacency in the state, so a single sweep
        // builds the nodal value from the fluid elements that contain the node.
        // The result is cached because the VTU writer needs it for every node
        // and the sweep is over the whole element list.
        Real nodal_reference_mu(const CBSStateSI& s, const Int ip)
        {
            static thread_local const CBSStateSI* cached_state = nullptr;
            static thread_local std::vector<Real> cached_mu;

            if (cached_state != &s ||
                cached_mu.size() != static_cast<std::size_t>(s.cfg.npoin) + 1U)
            {
                cached_mu.assign(
                    static_cast<std::size_t>(s.cfg.npoin) + 1U,
                    0.0);

                for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
                {
                    if (s.mat_elem(ie) != 0)
                    {
                        continue;
                    }

                    for (Int a = 1; a <= s.cfg.nep; ++a)
                    {
                        const Int node = s.intma(a, ie);

                        cached_mu[static_cast<std::size_t>(node)] = s.mu_e(ie);
                    }
                }

                cached_state = &s;
            }

            return cached_mu[static_cast<std::size_t>(ip)];
        }


        bool touches_fluid(const CBSStateSI& s, const Int ip)
        {
            return
                (s.node_material_mask(ip) &
                 CBSStateSI::node_touches_fluid) != 0;
        }

        bool touches_solid(const CBSStateSI& s, const Int ip)
        {
            return
                (s.node_material_mask(ip) &
                 CBSStateSI::node_touches_solid) != 0;
        }

        bool velocity_active(const CBSStateSI& s, const Int ip)
        {
            return touches_fluid(s, ip) && !touches_solid(s, ip);
        }

        long long global_node_id(const CBSStateSI& s, const Int ip)
        {
            const std::size_t index = static_cast<std::size_t>(ip);

            if (index >= s.local_to_global_node.size())
            {
                throw std::runtime_error(
                    "DistributedPost found an incomplete node map");
            }

            return static_cast<long long>(s.local_to_global_node[index]);
        }

        long long global_element_id(const CBSStateSI& s, const Int ie)
        {
            const std::size_t index = static_cast<std::size_t>(ie);

            if (index >= s.local_to_global_element.size())
            {
                throw std::runtime_error(
                    "DistributedPost found an incomplete element map");
            }

            return static_cast<long long>(s.local_to_global_element[index]);
        }

        void write_residual_header(const std::string& case_name)
        {
            std::ofstream output(
                residual_file_path(case_name),
                std::ios::trunc);

            if (!output)
            {
                throw std::runtime_error(
                    "DistributedPost cannot create the distributed residual CSV");
            }

            output
                << "iteration,time,dt,"
                << "u_rel,u_norm,u_abs,"
                << "v_rel,v_norm,v_abs,"
                << "w_rel,w_norm,w_abs,"
                << "p_rel,p_norm,p_abs,"
                << "T_rel,T_norm,T_abs,"
                << "velocity_rel_max,"
                << "continuity_rms,continuity_max,"
                << "maximum_velocity,maximum_velocity_correction,"
                << "cg_iterations,cg_initial_l2,cg_final_l2,"
                << "cg_relative_l2,cg_max_abs,iteration_wall_seconds,"
                << "sa_rel,sa_nu_tilde_min,sa_nu_tilde_max,sa_chi_max,"
                << "sa_nu_t_min,sa_nu_t_max,sa_mu_t_max,sa_mu_eff_max\n";
        }

        void write_piece(
            const CBSStateSI& s,
            const std::string& case_name,
            const Int iteration)
        {
            const fs::path output_path =
                case_output_directory(case_name) /
                piece_file_name(case_name, iteration, s.mpi_rank);

            std::ofstream output(output_path, std::ios::trunc);

            if (!output)
            {
                throw std::runtime_error(
                    "DistributedPost cannot create rank-local VTU piece");
            }

            const Int total_cells = s.cfg.nelem + s.cfg.nboun;

            output << std::setprecision(16);
            output << "<?xml version=\"1.0\"?>\n";
            output << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
            output << "  <UnstructuredGrid>\n";
            output << "    <Piece NumberOfPoints=\"" << s.cfg.npoin
                   << "\" NumberOfCells=\"" << total_cells << "\">\n";

            output << "      <Points>\n";
            output << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                output << "          "
                       << safe_value(s.coord(1, ip)) << ' '
                       << safe_value(s.coord(2, ip)) << ' '
                       << safe_value(s.coord(3, ip)) << '\n';
            }
            output << "        </DataArray>\n";
            output << "      </Points>\n";

            output << "      <Cells>\n";
            output << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                output << "          "
                       << s.intma(1, ie) - 1 << ' '
                       << s.intma(2, ie) - 1 << ' '
                       << s.intma(3, ie) - 1 << ' '
                       << s.intma(4, ie) - 1 << '\n';
            }
            for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
            {
                output << "          "
                       << s.iside(1, ib) - 1 << ' '
                       << s.iside(2, ib) - 1 << ' '
                       << s.iside(3, ib) - 1 << '\n';
            }
            output << "        </DataArray>\n";

            output << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
            long long offset = 0;
            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                offset += 4;
                output << "          " << offset << '\n';
            }
            for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
            {
                offset += 3;
                output << "          " << offset << '\n';
            }
            output << "        </DataArray>\n";

            output << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                output << "          10\n";
            }
            for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
            {
                output << "          5\n";
            }
            output << "        </DataArray>\n";
            output << "      </Cells>\n";

            output << "      <PointData Scalars=\"temperature\" Vectors=\"velocity\">\n";

            output << "        <DataArray type=\"Float64\" Name=\"pressure\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                output << "          "
                       << (touches_fluid(s, ip)
                           ? safe_value(s.pres(ip))
                           : 0.0)
                       << '\n';
            }
            output << "        </DataArray>\n";

            output << "        <DataArray type=\"Float64\" Name=\"temperature\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                output << "          " << safe_value(s.temperature(ip)) << '\n';
            }
            output << "        </DataArray>\n";

            output << "        <DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                if (velocity_active(s, ip))
                {
                    output << "          "
                           << safe_value(s.unkno(1, ip)) << ' '
                           << safe_value(s.unkno(2, ip)) << ' '
                           << safe_value(s.unkno(3, ip)) << '\n';
                }
                else
                {
                    output << "          0 0 0\n";
                }
            }
            output << "        </DataArray>\n";

            output << "        <DataArray type=\"Float64\" Name=\"velocity_magnitude\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                output << "          "
                       << (velocity_active(s, ip)
                           ? safe_value(s.velocity(ip))
                           : 0.0)
                       << '\n';
            }
            output << "        </DataArray>\n";

            // Spalart-Allmaras nodal fields.
            //
            // These are written whenever turbulence is enabled so that the
            // distributed output carries the same SA state as the serial path.
            // Without them the post-processing has to infer nu_t from nu_tilde
            // and a viscosity quoted elsewhere, which silently breaks for a
            // variable-property case, and the SA budget terms cannot be
            // examined at all.
            //
            // Every array is written on all local nodes including the ghost
            // layer.  Ghost entries carry the owner's value because the SA step
            // broadcasts after the update, and vtkGhostType is already written
            // above so a reader can drop duplicates.
            if (s.cfg.turbulence_on > 0)
            {
                struct SaPointField
                {
                    const char* name;
                    const Array1D<Real>* values;
                };

                const SaPointField sa_fields[] =
                {
                    { "nu_tilde",       &s.nu_tilde },
                    { "nu_t",           &s.nu_t },
                    { "mu_t",           &s.mu_t },
                    { "wall_distance",  &s.wall_distance },
                    { "sa_rhs",         &s.sa_rhs },
                    { "sa_residual",    &s.sa_residual },
                    { "sa_production",  &s.sa_production },
                    { "sa_destruction", &s.sa_destruction },
                    { "sa_diffusion",   &s.sa_diffusion },
                    { "sa_source",      &s.sa_source }
                };

                for (const SaPointField& field : sa_fields)
                {
                    output << "        <DataArray type=\"Float64\" Name=\""
                           << field.name
                           << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";

                    for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
                    {
                        output << "          "
                               << safe_value((*field.values)(ip)) << '\n';
                    }

                    output << "        </DataArray>\n";
                }

                // mu_t/mu is the quantity the NASA TMR reference plots use, so
                // it is exported directly rather than left to be reconstructed
                // by a post-processor that would have to guess the molecular
                // viscosity at each node.
                output << "        <DataArray type=\"Float64\" Name=\"mu_t_over_mu\""
                          " NumberOfComponents=\"1\" format=\"ascii\">\n";

                for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
                {
                    const Real reference_mu = nodal_reference_mu(s, ip);

                    output << "          "
                           << (reference_mu > 0.0
                               ? safe_value(s.mu_t(ip) / reference_mu)
                               : 0.0)
                           << '\n';
                }

                output << "        </DataArray>\n";

                struct SaFlagField
                {
                    const char* name;
                    const Array1D<Int>* values;
                };

                const SaFlagField sa_flags[] =
                {
                    { "sa_active_node", &s.sa_active_node },
                    { "sa_wall_node",   &s.sa_wall_node },
                    { "sa_inlet_node",  &s.sa_inlet_node }
                };

                for (const SaFlagField& flag : sa_flags)
                {
                    output << "        <DataArray type=\"Int32\" Name=\""
                           << flag.name
                           << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";

                    for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
                    {
                        output << "          " << (*flag.values)(ip) << '\n';
                    }

                    output << "        </DataArray>\n";
                }
            }

            output << "        <DataArray type=\"Int64\" Name=\"global_node_id\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                output << "          " << global_node_id(s, ip) << '\n';
            }
            output << "        </DataArray>\n";

            output << "        <DataArray type=\"Int32\" Name=\"owner_rank\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                output << "          "
                       << s.node_owner_rank[static_cast<std::size_t>(ip)]
                       << '\n';
            }
            output << "        </DataArray>\n";

            output << "        <DataArray type=\"UInt8\" Name=\"is_owned\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                output << "          "
                       << (s.node_owner_rank[static_cast<std::size_t>(ip)] == s.mpi_rank
                           ? 1
                           : 0)
                       << '\n';
            }
            output << "        </DataArray>\n";

            // VTK duplicate-point flag. Ghost copies are required by local
            // tetrahedral connectivity but are not independent global unknowns.
            output << "        <DataArray type=\"UInt8\" Name=\"vtkGhostType\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                output << "          "
                       << (s.node_owner_rank[static_cast<std::size_t>(ip)] == s.mpi_rank
                           ? 0
                           : 1)
                       << '\n';
            }
            output << "        </DataArray>\n";

            output << "        <DataArray type=\"Int32\" Name=\"node_domain_kind\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                output << "          " << s.node_material_mask(ip) << '\n';
            }
            output << "        </DataArray>\n";

            output << "        <DataArray type=\"Int32\" Name=\"velocity_bc_type\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                output << "          " << s.node_velocity_bc_type(ip) << '\n';
            }
            output << "        </DataArray>\n";

            output << "        <DataArray type=\"UInt8\" Name=\"pressure_fixed\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                output << "          " << (s.node_pressure_fixed(ip) != 0 ? 1 : 0) << '\n';
            }
            output << "        </DataArray>\n";
            output << "      </PointData>\n";

            output << "      <CellData Scalars=\"material_id\">\n";

            output << "        <DataArray type=\"UInt8\" Name=\"cell_kind\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                output << "          1\n";
            }
            for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
            {
                output << "          2\n";
            }
            output << "        </DataArray>\n";

            output << "        <DataArray type=\"Int64\" Name=\"global_element_id\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                output << "          " << global_element_id(s, ie) << '\n';
            }
            for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
            {
                output << "          -1\n";
            }
            output << "        </DataArray>\n";

            output << "        <DataArray type=\"Int32\" Name=\"material_id\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                output << "          " << s.mat_elem(ie) << '\n';
            }
            for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
            {
                output << "          -1\n";
            }
            output << "        </DataArray>\n";

            // Spalart-Allmaras element fields.
            //
            // The eddy viscosity actually used by the momentum and energy
            // assemblies is the element quantity mu_eff_e, not the nodal
            // average written above, so it is exported directly.  Boundary
            // faces are appended to the cell list and carry -1.
            if (s.cfg.turbulence_on > 0)
            {
                struct SaCellField
                {
                    const char* name;
                    const Array1D<Real>* values;
                };

                const SaCellField sa_cells[] =
                {
                    { "nu_tilde_e", &s.nu_tilde_e },
                    { "nu_t_e",     &s.nu_t_e },
                    { "mu_t_e",     &s.mu_t_e },
                    { "mu_eff_e",   &s.mu_eff_e },
                    { "rho_e",      &s.rho_e },
                    { "mu_e",       &s.mu_e }
                };

                for (const SaCellField& field : sa_cells)
                {
                    output << "        <DataArray type=\"Float64\" Name=\""
                           << field.name
                           << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";

                    for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
                    {
                        output << "          "
                               << safe_value((*field.values)(ie)) << '\n';
                    }

                    for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
                    {
                        output << "          -1\n";
                    }

                    output << "        </DataArray>\n";
                }
            }

            output << "        <DataArray type=\"Int32\" Name=\"bc_id\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                output << "          0\n";
            }
            for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
            {
                output << "          " << s.iside(s.cfg.bsid, ib) << '\n';
            }
            output << "        </DataArray>\n";

            output << "        <DataArray type=\"Int64\" Name=\"parent_global_element\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                output << "          " << global_element_id(s, ie) << '\n';
            }
            for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
            {
                const Int parent = s.iside(s.cfg.nsidpe, ib);
                output << "          " << global_element_id(s, parent) << '\n';
            }
            output << "        </DataArray>\n";

            output << "      </CellData>\n";
            output << "    </Piece>\n";
            output << "  </UnstructuredGrid>\n";
            output << "</VTKFile>\n";
        }

        void write_parallel_descriptor(
            const CBSStateSI& s,
            const std::string& case_name,
            const Int iteration)
        {
            const fs::path output_path =
                case_output_directory(case_name) /
                parallel_file_name(case_name, iteration);

            std::ofstream output(output_path, std::ios::trunc);

            if (!output)
            {
                throw std::runtime_error(
                    "DistributedPost cannot create PVTU descriptor");
            }

            output << "<?xml version=\"1.0\"?>\n";
            output << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
            output << "  <PUnstructuredGrid GhostLevel=\"0\">\n";
            output << "    <PPointData Scalars=\"temperature\" Vectors=\"velocity\">\n";
            output << "      <PDataArray type=\"Float64\" Name=\"pressure\" NumberOfComponents=\"1\"/>\n";
            output << "      <PDataArray type=\"Float64\" Name=\"temperature\" NumberOfComponents=\"1\"/>\n";
            output << "      <PDataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\"/>\n";
            output << "      <PDataArray type=\"Float64\" Name=\"velocity_magnitude\" NumberOfComponents=\"1\"/>\n";
            // The SA declarations must match the piece files exactly, or a
            // reader silently drops the arrays it cannot find in every piece.
            if (s.cfg.turbulence_on > 0)
            {
                const char* sa_real_names[] =
                {
                    "nu_tilde", "nu_t", "mu_t", "wall_distance",
                    "sa_rhs", "sa_residual", "sa_production",
                    "sa_destruction", "sa_diffusion", "sa_source",
                    "mu_t_over_mu"
                };

                for (const char* name : sa_real_names)
                {
                    output << "      <PDataArray type=\"Float64\" Name=\""
                           << name << "\" NumberOfComponents=\"1\"/>\n";
                }

                const char* sa_flag_names[] =
                {
                    "sa_active_node", "sa_wall_node", "sa_inlet_node"
                };

                for (const char* name : sa_flag_names)
                {
                    output << "      <PDataArray type=\"Int32\" Name=\""
                           << name << "\" NumberOfComponents=\"1\"/>\n";
                }
            }

            output << "      <PDataArray type=\"Int64\" Name=\"global_node_id\" NumberOfComponents=\"1\"/>\n";
            output << "      <PDataArray type=\"Int32\" Name=\"owner_rank\" NumberOfComponents=\"1\"/>\n";
            output << "      <PDataArray type=\"UInt8\" Name=\"is_owned\" NumberOfComponents=\"1\"/>\n";
            output << "      <PDataArray type=\"UInt8\" Name=\"vtkGhostType\" NumberOfComponents=\"1\"/>\n";
            output << "      <PDataArray type=\"Int32\" Name=\"node_domain_kind\" NumberOfComponents=\"1\"/>\n";
            output << "      <PDataArray type=\"Int32\" Name=\"velocity_bc_type\" NumberOfComponents=\"1\"/>\n";
            output << "      <PDataArray type=\"UInt8\" Name=\"pressure_fixed\" NumberOfComponents=\"1\"/>\n";
            output << "    </PPointData>\n";
            output << "    <PCellData Scalars=\"material_id\">\n";
            output << "      <PDataArray type=\"UInt8\" Name=\"cell_kind\" NumberOfComponents=\"1\"/>\n";
            output << "      <PDataArray type=\"Int64\" Name=\"global_element_id\" NumberOfComponents=\"1\"/>\n";
            output << "      <PDataArray type=\"Int32\" Name=\"material_id\" NumberOfComponents=\"1\"/>\n";
            if (s.cfg.turbulence_on > 0)
            {
                const char* sa_cell_names[] =
                {
                    "nu_tilde_e", "nu_t_e", "mu_t_e",
                    "mu_eff_e", "rho_e", "mu_e"
                };

                for (const char* name : sa_cell_names)
                {
                    output << "      <PDataArray type=\"Float64\" Name=\""
                           << name << "\" NumberOfComponents=\"1\"/>\n";
                }
            }

            output << "      <PDataArray type=\"Int32\" Name=\"bc_id\" NumberOfComponents=\"1\"/>\n";
            output << "      <PDataArray type=\"Int64\" Name=\"parent_global_element\" NumberOfComponents=\"1\"/>\n";
            output << "    </PCellData>\n";
            output << "    <PPoints>\n";
            output << "      <PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>\n";
            output << "    </PPoints>\n";

            for (Int rank = 0; rank < s.mpi_size; ++rank)
            {
                output << "    <Piece Source=\""
                       << piece_file_name(case_name, iteration, rank)
                       << "\"/>\n";
            }

            output << "  </PUnstructuredGrid>\n";
            output << "</VTKFile>\n";
        }

        void write_pvd(
            const CBSStateSI& s,
            const std::string& case_name)
        {
            if (s.output_iterations.size() != s.output_times.size())
            {
                throw std::runtime_error(
                    "DistributedPost output history is inconsistent");
            }

            std::ofstream output(pvd_file_path(case_name), std::ios::trunc);

            if (!output)
            {
                throw std::runtime_error(
                    "DistributedPost cannot create distributed PVD file");
            }

            output << std::setprecision(16);
            output << "<?xml version=\"1.0\"?>\n";
            output << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
            output << "  <Collection>\n";

            for (std::size_t index = 0;
                 index < s.output_iterations.size();
                 ++index)
            {
                output << "    <DataSet timestep=\""
                       << safe_value(s.output_times[index])
                       << "\" group=\"\" part=\"0\" file=\""
                       << parallel_file_name(
                              case_name,
                              s.output_iterations[index])
                       << "\"/>\n";
            }

            output << "  </Collection>\n";
            output << "</VTKFile>\n";
        }
    }
#endif

    std::string DistributedPost::distributedCaseName(
        const std::string& rank_local_case_name)
    {
#if defined(CBS3D_USE_MPI)
        return case_name_from_local_path(rank_local_case_name);
#else
        (void)rank_local_case_name;
        throw std::runtime_error(
            "DistributedPost requires an MPI-enabled build");
#endif
    }

    void DistributedPost::initialise(
        CBSStateSI& s,
        const std::string& rank_local_case_name)
    {
#if defined(CBS3D_USE_MPI)
        if (!s.mpi_enabled)
        {
            throw std::runtime_error(
                "DistributedPost::initialise requires more than one MPI rank");
        }

        const std::string case_name =
            case_name_from_local_path(rank_local_case_name);

        if (s.mpi_rank == 0)
        {
            fs::create_directories(case_output_directory(case_name));
        }

        check_mpi(
            MPI_Barrier(MPI_COMM_WORLD),
            "MPI_Barrier output directory");

        s.output_iterations.clear();
        s.output_times.clear();
        s.next_vtu_output_time = s.cfg.vtu_output_every_sim_time;

        if (s.mpi_rank == 0 && s.cfg.residual_log_enabled > 0)
        {
            write_residual_header(case_name);
        }

        check_mpi(
            MPI_Barrier(MPI_COMM_WORLD),
            "MPI_Barrier residual initialisation");

        if (s.cfg.vtu_output_enabled > 0)
        {
            writeSolution(s, rank_local_case_name, 0);
        }
#else
        (void)s;
        (void)rank_local_case_name;
        throw std::runtime_error(
            "DistributedPost::initialise requires an MPI-enabled build");
#endif
    }

    bool DistributedPost::shouldWriteSolution(
        CBSStateSI& s,
        const Int iteration)
    {
        if (s.cfg.vtu_output_enabled < 1)
        {
            return false;
        }

        if (iteration == 0)
        {
            return true;
        }

        const bool by_iteration =
            s.cfg.vtu_output_every_iterations > 0 &&
            (iteration % s.cfg.vtu_output_every_iterations) == 0;

        bool by_time = false;

        if (s.cfg.transient_on > 0 &&
            s.cfg.vtu_output_every_sim_time > 0.0 &&
            s.cfg.rtime + 1.0e-14 >= s.next_vtu_output_time)
        {
            by_time = true;

            while (s.next_vtu_output_time <= s.cfg.rtime + 1.0e-14)
            {
                s.next_vtu_output_time +=
                    s.cfg.vtu_output_every_sim_time;
            }
        }

        return by_iteration || by_time;
    }

    void DistributedPost::writeSolution(
        CBSStateSI& s,
        const std::string& rank_local_case_name,
        const Int iteration)
    {
#if defined(CBS3D_USE_MPI)
        if (!s.mpi_enabled)
        {
            throw std::runtime_error(
                "DistributedPost::writeSolution requires more than one MPI rank");
        }

        const std::string case_name =
            case_name_from_local_path(rank_local_case_name);

        write_piece(s, case_name, iteration);

        check_mpi(
            MPI_Barrier(MPI_COMM_WORLD),
            "MPI_Barrier VTU pieces");

        s.output_iterations.push_back(iteration);
        s.output_times.push_back(s.cfg.rtime);

        if (s.mpi_rank == 0)
        {
            write_parallel_descriptor(s, case_name, iteration);
            write_pvd(s, case_name);
        }

        check_mpi(
            MPI_Barrier(MPI_COMM_WORLD),
            "MPI_Barrier PVTU/PVD metadata");
#else
        (void)s;
        (void)rank_local_case_name;
        (void)iteration;
        throw std::runtime_error(
            "DistributedPost::writeSolution requires an MPI-enabled build");
#endif
    }

    void DistributedPost::writeResidualRow(
        const CBSStateSI& s,
        const std::string& rank_local_case_name,
        const Int iteration,
        const Real continuity_rms,
        const Real continuity_max,
        const Real maximum_velocity,
        const Real maximum_velocity_correction,
        const Real iteration_wall_seconds,
        const Convergence::TurbulenceDiagnostics& turbulence)
    {
#if defined(CBS3D_USE_MPI)
        if (s.mpi_rank != 0 || s.cfg.residual_log_enabled < 1)
        {
            return;
        }

        if (s.cfg.residual_log_every > 1 &&
            (iteration % s.cfg.residual_log_every) != 0)
        {
            return;
        }

        const std::string case_name =
            case_name_from_local_path(rank_local_case_name);

        std::ofstream output(
            residual_file_path(case_name),
            std::ios::app);

        if (!output)
        {
            throw std::runtime_error(
                "DistributedPost cannot append the distributed residual CSV");
        }

        output << std::setprecision(16)
               << iteration << ','
               << safe_value(s.cfg.rtime) << ','
               << safe_value(s.cfg.dtreal) << ',';

        for (const Real value : s.hb)
        {
            output << safe_value(value) << ',';
        }

        output
            << safe_value(std::max({s.hb[0], s.hb[3], s.hb[6]})) << ','
            << safe_value(continuity_rms) << ','
            << safe_value(continuity_max) << ','
            << safe_value(maximum_velocity) << ','
            << safe_value(maximum_velocity_correction) << ','
            << s.last_cg_iterations << ','
            << safe_value(s.last_cg_initial_l2) << ','
            << safe_value(s.last_cg_final_l2) << ','
            << safe_value(s.last_cg_relative_l2) << ','
            << safe_value(s.last_cg_max_abs) << ','
            << safe_value(iteration_wall_seconds) << ','
            << safe_value(turbulence.residual) << ','
            << safe_value(turbulence.nu_tilde_min) << ','
            << safe_value(turbulence.nu_tilde_max) << ','
            << safe_value(turbulence.chi_max) << ','
            << safe_value(turbulence.nu_t_min) << ','
            << safe_value(turbulence.nu_t_max) << ','
            << safe_value(turbulence.mu_t_max) << ','
            << safe_value(turbulence.mu_eff_max)
            << '\n';
#else
        (void)s;
        (void)rank_local_case_name;
        (void)iteration;
        (void)continuity_rms;
        (void)continuity_max;
        (void)maximum_velocity;
        (void)maximum_velocity_correction;
        (void)iteration_wall_seconds;
        (void)turbulence;
        throw std::runtime_error(
            "DistributedPost::writeResidualRow requires an MPI-enabled build");
#endif
    }

    bool DistributedPost::solutionAlreadyWritten(
        const CBSStateSI& s,
        const Int iteration)
    {
        return
            std::find(
                s.output_iterations.begin(),
                s.output_iterations.end(),
                iteration) != s.output_iterations.end();
    }
}
