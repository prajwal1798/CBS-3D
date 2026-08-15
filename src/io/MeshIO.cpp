//=============================================================================
// CBS3D++_SI
// Mesh, boundary-condition, parameter and material-file reader
//=============================================================================

#include "cbs/io/MeshIO.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace cbs
{
    namespace
    {
        //-------------------------------------------------------------------------
        // Removes leading and trailing spaces from an input line.
        //-------------------------------------------------------------------------
        std::string trim_line(const std::string& line)
        {
            const std::size_t first = line.find_first_not_of(" \t\r\n");

            if (first == std::string::npos)
            {
                return "";
            }

            const std::size_t last = line.find_last_not_of(" \t\r\n");

            return line.substr(first, last - first + 1);
        }

        //-------------------------------------------------------------------------
        // Reads the next non-empty, non-comment line from an input file.
        // Lines beginning with '#' or '!' are treated as comments.
        //-------------------------------------------------------------------------
        bool read_data_line(std::istream& input, std::string& line)
        {
            std::string raw_line;

            while (std::getline(input, raw_line))
            {
                raw_line = trim_line(raw_line);

                if (raw_line.empty())
                {
                    continue;
                }

                if (raw_line[0] == '#' || raw_line[0] == '!')
                {
                    continue;
                }

                line = raw_line;
                return true;
            }

            return false;
        }

        //-------------------------------------------------------------------------
        // Reads the next required data line from an input file.
        //-------------------------------------------------------------------------
        std::string next_data_line(
            std::istream& input,
            const std::string& file_name)
        {
            std::string line;

            if (!read_data_line(input, line))
            {
                throw std::runtime_error(
                    "Unexpected end of file: " + file_name);
            }

            return line;
        }

        //-------------------------------------------------------------------------
        // Reads all non-empty, non-comment lines from a file.
        //-------------------------------------------------------------------------
        std::vector<std::string> read_data_lines(
            const std::string& file_name)
        {
            std::ifstream input(file_name);

            if (!input)
            {
                throw std::runtime_error(
                    "Cannot open file: " + file_name);
            }

            std::vector<std::string> lines;
            std::string line;

            while (read_data_line(input, line))
            {
                lines.push_back(line);
            }

            return lines;
        }

        void set_freestream_initial_condition(CBSStateSI& s)
        {
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                s.unkno(1, ip) = s.cfg.cinf[0];
                s.unkno(2, ip) = s.cfg.cinf[1];
                s.unkno(3, ip) = s.cfg.cinf[2];

                s.pres(ip) = s.cfg.cinf[3];
                s.temperature(ip) = s.cfg.cinf[4];
            }

            s.cfg.istart = 1;
            s.cfg.rtime = 0.0;
        }

        void read_optional_alpha_source_flux(std::istream& in, CBSStateSI& s)
        {
            std::string label;
            std::string line;

            if (!read_data_line(in, label))
            {
                return;
            }

            if (!read_data_line(in, line))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - missing alpha/source/flux data line");
            }

            std::istringstream iss(line);

            if (!(iss >> s.cfg.alpha_sf
                >> s.cfg.k_ratio
                >> s.cfg.source_solid
                >> s.cfg.heat_flux_bc))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - invalid alpha/source/flux line");
            }
        }

        void read_cbs3d_timestep_controls(std::istream& in, CBSStateSI& s)
        {
            std::string label;
            std::string line;

            if (!read_data_line(in, label))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - missing CBS3D timestep-control label line");
            }

            if (!read_data_line(in, line))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - missing CBS3D timestep-control data line");
            }

            std::istringstream iss(line);

            if (!(iss >> s.cfg.cbs_scheme
                >> s.cfg.ilots
                >> s.cfg.htype
                >> s.cfg.step2_check
                >> s.cfg.rem_deltp
                >> s.cfg.beta_opt
                >> s.cfg.dtfix_end
                >> s.cfg.csafm2
                >> s.cfg.epsilon1
                >> s.cfg.deltr))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - invalid CBS3D timestep-control line; "
                    "expected cbs_scheme ilots htype step2_check rem_deltp beta_opt "
                    "dtfix_end csafm2 epsilon1 deltr");
            }
        }

        void read_optional_dimensional_controls(std::istream& in, CBSStateSI& s)
        {
            std::string label;
            std::string line;

            if (!read_data_line(in, label))
            {
                return;
            }

            if (!read_data_line(in, line))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - missing dimensional-control data line");
            }

            {
                std::istringstream iss(line);

                if (!(iss >> s.cfg.dimensional_mode
                    >> s.cfg.material_properties_enabled
                    >> s.cfg.mass_flow_inlet_enabled))
                {
                    throw std::runtime_error(
                        "MeshIO::readParameterFile - invalid dimensional-control line");
                }
            }

            if (!read_data_line(in, label))
            {
                return;
            }

            if (!read_data_line(in, line))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - missing reference/inlet/outlet data line");
            }

            {
                std::istringstream iss(line);

                if (!(iss >> s.cfg.p_ref
                    >> s.cfg.model_depth
                    >> s.cfg.inlet_mass_flow_rate
                    >> s.cfg.inlet_density
                    >> s.cfg.inlet_temperature
                    >> s.cfg.outlet_pressure_gauge))
                {
                    throw std::runtime_error(
                        "MeshIO::readParameterFile - invalid reference/inlet/outlet line");
                }
            }

            if (!read_data_line(in, label))
            {
                return;
            }

            if (!read_data_line(in, line))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - missing prescribed inlet velocity line");
            }

            {
                std::istringstream iss(line);

                if (!(iss >> s.cfg.inlet_u >> s.cfg.inlet_v >> s.cfg.inlet_w))
                {
                    throw std::runtime_error(
                        "MeshIO::readParameterFile - invalid prescribed inlet velocity line; "
                        "CBS3D expects inlet_u inlet_v inlet_w");
                }
            }

            if (!read_data_line(in, label))
            {
                return;
            }

            if (!read_data_line(in, line))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - missing Nusselt reference line");
            }

            {
                std::istringstream iss(line);

                if (!(iss >> s.cfg.nusselt_Tinf
                    >> s.cfg.nusselt_Tref
                    >> s.cfg.nusselt_diameter))
                {
                    throw std::runtime_error(
                        "MeshIO::readParameterFile - invalid Nusselt reference line");
                }
            }
        }


        void read_pressure_cg_controls(std::istream& in, CBSStateSI& s)
        {
            std::string label;
            std::string line;

            if (!read_data_line(in, label))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - missing pressure-CG-control label line");
            }

            if (!read_data_line(in, line))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - missing pressure-CG-control data line");
            }

            std::istringstream iss(line);

            if (!(iss >> s.cfg.cg_preconditioner
                >> s.cfg.cg_conv_test
                >> s.cfg.cg_max_iter
                >> s.cfg.l2norm_pres_tolerance))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - invalid pressure-CG-control line; "
                    "expected cg_preconditioner cg_conv_test cg_max_iter l2norm_pres_tolerance");
            }
        }

        void read_output_monitor_controls(std::istream& in, CBSStateSI& s)
        {
            std::string label;
            std::string line;

            if (!read_data_line(in, label))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - missing CBS3D output/monitor-control label line");
            }

            if (!read_data_line(in, line))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - missing CBS3D output/monitor-control data line");
            }

            std::istringstream iss(line);

            if (!(iss >> s.cfg.residual_log_enabled
                >> s.cfg.residual_log_every
                >> s.cfg.console_log_every
                >> s.cfg.live_residual_plot
                >> s.cfg.vtu_output_enabled
                >> s.cfg.vtu_output_every_iterations
                >> s.cfg.vtu_output_every_sim_time
                >> s.cfg.write_boundary_debug_arrays
                >> s.cfg.steady_min_iterations))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - invalid CBS3D output/monitor-control line; "
                    "expected residual_log_enabled residual_log_every console_log_every live_residual_plot "
                    "vtu_output_enabled vtu_output_every_iterations vtu_output_every_sim_time "
                    "write_boundary_debug_arrays steady_min_iterations");
            }
        }

        void read_artificial_diffusion_control(std::istream& in, CBSStateSI& s)
        {
            std::string label;
            std::string line;

            if (!read_data_line(in, label))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - missing artificial-diffusion-control label line");
            }

            if (!read_data_line(in, line))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - missing artificial-diffusion-control data line");
            }

            std::istringstream iss(line);

            if (!(iss >> s.cfg.art_diff))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - invalid artificial-diffusion-control line; expected art_diff");
            }
        }


        //-------------------------------------------------------------------------
        // Reads the optional Spalart-Allmaras turbulence-control block.
        //
        // The block is deliberately optional so that all older laminar and CHT
        // parameter files remain valid.  If the file ends after the artificial
        // diffusion line, the turbulence controls keep their RunConfig defaults:
        //
        //     turbulence_on = 0
        //
        // New SA cases append four labelled data blocks:
        //
        //     turbulence_on turbulence_model turbulent_thermal_diffusivity_on
        //     sa_inlet_ratio sa_prandtl_t
        //     sa_min_wall_distance sa_min_stilde sa_nu_tilde_floor
        //     sa_use_stilde_limiter sa_implicit_destruction [sa_nu_tilde_ceiling_ratio]
        //-------------------------------------------------------------------------
        void read_optional_spalart_allmaras_controls(std::istream& in, CBSStateSI& s)
        {
            std::string label;
            std::string line;

            if (!read_data_line(in, label))
            {
                return;
            }

            if (!read_data_line(in, line))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - missing Spalart-Allmaras control data line");
            }

            {
                std::istringstream iss(line);

                if (!(iss >> s.cfg.turbulence_on
                    >> s.cfg.turbulence_model
                    >> s.cfg.turbulent_thermal_diffusivity_on))
                {
                    throw std::runtime_error(
                        "MeshIO::readParameterFile - invalid Spalart-Allmaras control line; "
                        "expected turbulence_on turbulence_model turbulent_thermal_diffusivity_on");
                }
            }

            if (!read_data_line(in, label))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - missing Spalart-Allmaras inlet/thermal label line");
            }

            if (!read_data_line(in, line))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - missing Spalart-Allmaras inlet/thermal data line");
            }

            {
                std::istringstream iss(line);

                if (!(iss >> s.cfg.sa_inlet_ratio >> s.cfg.sa_prandtl_t))
                {
                    throw std::runtime_error(
                        "MeshIO::readParameterFile - invalid Spalart-Allmaras inlet/thermal line; "
                        "expected sa_inlet_ratio sa_prandtl_t");
                }
            }

            if (!read_data_line(in, label))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - missing Spalart-Allmaras numerical-floor label line");
            }

            if (!read_data_line(in, line))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - missing Spalart-Allmaras numerical-floor data line");
            }

            {
                std::istringstream iss(line);

                if (!(iss >> s.cfg.sa_min_wall_distance
                    >> s.cfg.sa_min_stilde
                    >> s.cfg.sa_nu_tilde_floor))
                {
                    throw std::runtime_error(
                        "MeshIO::readParameterFile - invalid Spalart-Allmaras numerical-floor line; "
                        "expected sa_min_wall_distance sa_min_stilde sa_nu_tilde_floor");
                }
            }

            if (!read_data_line(in, label))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - missing Spalart-Allmaras switch label line");
            }

            if (!read_data_line(in, line))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - missing Spalart-Allmaras switch data line");
            }

            {
                std::istringstream iss(line);

                if (!(iss >> s.cfg.sa_use_stilde_limiter
                    >> s.cfg.sa_implicit_destruction))
                {
                    throw std::runtime_error(
                        "MeshIO::readParameterFile - invalid Spalart-Allmaras switch line; "
                        "expected sa_use_stilde_limiter sa_implicit_destruction");
                }

                // Optional third field: sa_nu_tilde_ceiling_ratio.
                //
                // It is read only if present so that existing parameter files
                // remain valid and keep the built-in default.  Zero disables the
                // divergence bound.
                Real ceiling_ratio = 0.0;

                if (iss >> ceiling_ratio)
                {
                    if (ceiling_ratio < 0.0 || !std::isfinite(ceiling_ratio))
                    {
                        throw std::runtime_error(
                            "MeshIO::readParameterFile - sa_nu_tilde_ceiling_ratio "
                            "must be zero (disabled) or positive");
                    }

                    s.cfg.sa_nu_tilde_ceiling_ratio = ceiling_ratio;
                }
            }
        }

        bool material_files_required(const CBSStateSI& s)
        {
            // Strict rule:
            //
            //   material_properties_enabled = 1
            //       .material and .matprop are required.
            //
            //   material_properties_enabled = 0
            //       no material files are required.  Every element is treated
            //       as fluid material 0 with safe default properties.
            //
            // This allows pure flow cases such as 3D lid-driven cavity to use
            // only .plt, .bco and .par when temp_calc=0 and no material model is
            // requested.
            return s.cfg.material_properties_enabled > 0;
        }

        void initialise_default_flow_only_materials(CBSStateSI& s)
        {
            // Safe default material state for non-CHT / no-material-file runs.
            //
            // mat_elem = 0 identifies fluid elements throughout the CBS3D C++
            // solver.  The default thermal values are intentionally positive so
            // that optional diagnostics or accidental temperature allocation do
            // not encounter zero capacitance/conductivity.  They are not used
            // by the momentum equation unless material_properties_enabled > 0.
            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                s.mat_elem(ie) = 0;

                s.rho_e(ie) = 1.0;
                s.cp_e(ie) = 1.0;
                s.rho_cp_e(ie) = 1.0;

                s.mu_e(ie) = 0.0;
                s.k_e(ie) = 1.0;
                s.alpha_e(ie) = 1.0;
                s.Qvol_e(ie) = 0.0;

                s.mu_eff_e(ie) = s.mu_e(ie);
                s.k_eff_e(ie) = s.k_e(ie);
                s.nu_tilde_e(ie) = 0.0;
                s.nu_t_e(ie) = 0.0;
                s.mu_t_e(ie) = 0.0;
            }
        }

        void check_material_property_values(
            Int id,
            Real rho,
            Real cp,
            Real k,
            Real mu)
        {
            if (rho <= 0.0 || cp <= 0.0 || k <= 0.0)
            {
                throw std::runtime_error(
                    "MeshIO::readMaterialPropertyFile - material " + std::to_string(id) +
                    " has non-positive rho/cp/k");
            }

            if (mu < 0.0)
            {
                throw std::runtime_error(
                    "MeshIO::readMaterialPropertyFile - material " + std::to_string(id) +
                    " has negative viscosity");
            }
        }

        bool supported_solver_bc_id(Int bc)
        {
            switch (bc)
            {
            case 500:
            case 501:
            case 502:
            case 503:
            case 506:
            case 507:
            case 508:
            case 510:
            case 511:
            case 520:
            case 530:
            case 532:
            case 901:
            case 902:
                return true;
            default:
                return false;
            }
        }

        void enforce_initial_material_domain_fields(CBSStateSI& s)
        {
            if (s.cfg.material_properties_enabled < 1)
            {
                return;
            }

            std::vector<char> touches_fluid(
                static_cast<std::size_t>(s.cfg.npoin + 1), 0);
            std::vector<char> touches_solid(
                static_cast<std::size_t>(s.cfg.npoin + 1), 0);

            Int fluid_elements = 0;
            Int solid_elements = 0;

            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                const bool fluid_element = (s.mat_elem(ie) == 0);

                if (fluid_element)
                {
                    ++fluid_elements;
                }
                else
                {
                    ++solid_elements;
                }

                for (Int inod = 1; inod <= s.cfg.nep; ++inod)
                {
                    const Int ip = s.intma(inod, ie);

                    if (fluid_element)
                    {
                        touches_fluid[static_cast<std::size_t>(ip)] = 1;
                    }
                    else
                    {
                        touches_solid[static_cast<std::size_t>(ip)] = 1;
                    }
                }
            }

            Int fluid_only_nodes = 0;
            Int solid_only_nodes = 0;
            Int interface_nodes = 0;

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                const bool fluid =
                    touches_fluid[static_cast<std::size_t>(ip)] != 0;
                const bool solid =
                    touches_solid[static_cast<std::size_t>(ip)] != 0;

                if (fluid && solid)
                {
                    ++interface_nodes;
                }
                else if (fluid)
                {
                    ++fluid_only_nodes;
                }
                else if (solid)
                {
                    ++solid_only_nodes;
                }
                else
                {
                    throw std::runtime_error(
                        "MeshIO::initialiseFields - orphan node not touched by any material element: "
                        + std::to_string(ip));
                }

                if (solid)
                {
                    for (Int idim = 1; idim <= s.cfg.ndim; ++idim)
                    {
                        s.unkno(idim, ip) = 0.0;
                        s.unkn1(idim, ip) = 0.0;
                    }
                }

                if (solid && !fluid)
                {
                    s.pres(ip) = 0.0;
                    s.pres1(ip) = 0.0;
                }
            }

            std::cout
                << "Material domain audit:\n"
                << "  Fluid elements              : " << fluid_elements << "\n"
                << "  Solid elements              : " << solid_elements << "\n"
                << "  Fluid-only nodes            : " << fluid_only_nodes << "\n"
                << "  Solid-only nodes            : " << solid_only_nodes << "\n"
                << "  Fluid-solid interface nodes : " << interface_nodes << "\n";
        }
    }

    // =========================================================================

    void MeshIO::readAll(const std::string& case_name, CBSStateSI& s)
    {
        s.cfg.case_name = case_name;
        s.initialise_local_topology();

        const CaseFiles files(case_name);

        readSizes(files, s);

        s.set_problem_sizes(
            s.cfg.nelem,
            s.cfg.npoin,
            s.cfg.nboun,
            s.cfg.nflag);

        readMeshFile(files, s);

        readBoundaryFile(files, s);

        readParameterFile(files, s);

        if (material_files_required(s))
        {
            readMaterialFile(files, s);
            readMaterialPropertyFile(files, s);
        }
        else
        {
            initialise_default_flow_only_materials(s);

            std::cout
                << "Material files disabled by .par controls: "
                << "using default all-fluid material state "
                << "(mat_elem=0) for all elements.\n";
        }

        initialiseFields(s);
    }

    // =========================================================================

    void MeshIO::readSizes(const CaseFiles& files, CBSStateSI& s)
    {
        {
            std::ifstream in(files.par);
            if (!in)
            {
                throw std::runtime_error(
                    "MeshIO::readSizes - cannot open parameter file: " + files.par);
            }

            next_data_line(in, files.par);
            {
                const std::string line = next_data_line(in, files.par);
                std::istringstream iss(line);

                if (!(iss >> s.cfg.solver_opt))
                {
                    throw std::runtime_error("MeshIO::readSizes - invalid SOLVER_OPT");
                }
            }

            if (s.cfg.solver_opt < 1 || s.cfg.solver_opt > 3)
            {
                throw std::runtime_error(
                    "MeshIO::readSizes - SOLVER_OPT must be 1, 2, or 3");
            }
        }

        {
            std::ifstream in(files.plt);
            if (!in)
            {
                throw std::runtime_error(
                    "MeshIO::readSizes - cannot open mesh file: " + files.plt);
            }

            const std::string line = next_data_line(in, files.plt);
            std::istringstream iss(line);

            if (s.cfg.solver_opt == 1 || s.cfg.solver_opt == 3)
            {
                if (!(iss >> s.cfg.nelem >> s.cfg.npoin >> s.cfg.nboun))
                {
                    throw std::runtime_error("MeshIO::readSizes - invalid mesh header");
                }

                s.cfg.nbw = 0;
            }
            else
            {
                if (!(iss >> s.cfg.nelem >> s.cfg.npoin >> s.cfg.nboun >> s.cfg.nbw))
                {
                    throw std::runtime_error(
                        "MeshIO::readSizes - invalid banded mesh header");
                }
            }

            if (s.cfg.nelem < 1 || s.cfg.npoin < 1 || s.cfg.nboun < 0)
            {
                throw std::runtime_error(
                    "MeshIO::readSizes - non-positive mesh size in .plt header");
            }
        }

        {
            std::ifstream in(files.bco);
            if (!in)
            {
                throw std::runtime_error(
                    "MeshIO::readSizes - cannot open boundary file: " + files.bco);
            }

            const std::string line = next_data_line(in, files.bco);
            std::istringstream iss(line);

            Int mass_check = 0;

            if (!(iss >> s.cfg.nflag))
            {
                throw std::runtime_error("MeshIO::readSizes - invalid .bco header");
            }

            // 3D CBS .bco generally uses: nflag mass_check.
            // The CHT-first C++ path reads mass_check for compatibility but does
            // not use mass-group data yet.
            if (!(iss >> mass_check))
            {
                mass_check = 0;
            }

            if (s.cfg.nflag < 1)
            {
                throw std::runtime_error(
                    "MeshIO::readSizes - .bco must contain at least one BC mapping");
            }
        }
    }

    // =========================================================================

    void MeshIO::readMeshFile(const CaseFiles& files, CBSStateSI& s)
    {
        std::ifstream in(files.plt);
        if (!in)
        {
            throw std::runtime_error(
                "MeshIO::readMeshFile - cannot open mesh file: " + files.plt);
        }

        {
            const std::string line = next_data_line(in, files.plt);
            std::istringstream iss(line);

            Int nelem = 0;
            Int npoin = 0;
            Int nboun = 0;
            Int nbw = 0;

            if (s.cfg.solver_opt == 1 || s.cfg.solver_opt == 3)
            {
                if (!(iss >> nelem >> npoin >> nboun))
                {
                    throw std::runtime_error(
                        "MeshIO::readMeshFile - invalid mesh header");
                }
            }
            else
            {
                if (!(iss >> nelem >> npoin >> nboun >> nbw))
                {
                    throw std::runtime_error(
                        "MeshIO::readMeshFile - invalid banded mesh header");
                }
            }

            if (nelem != s.cfg.nelem ||
                npoin != s.cfg.npoin ||
                nboun != s.cfg.nboun)
            {
                throw std::runtime_error(
                    "MeshIO::readMeshFile - mesh header mismatch with allocated sizes");
            }
        }

        // ---------------------------------------------------------------------
        // Tetrahedral connectivity:
        //     element_id n1 n2 n3 n4
        // ---------------------------------------------------------------------
        for (Int i = 1; i <= s.cfg.nelem; ++i)
        {
            Int ie = 0;

            if (!(in >> ie))
            {
                throw std::runtime_error(
                    "MeshIO::readMeshFile - cannot read element index");
            }

            if (ie < 1 || ie > s.cfg.nelem)
            {
                throw std::runtime_error(
                    "MeshIO::readMeshFile - element index out of range");
            }

            for (Int j = 1; j <= s.cfg.nep; ++j)
            {
                if (!(in >> s.intma(j, ie)))
                {
                    throw std::runtime_error(
                        "MeshIO::readMeshFile - cannot read tetrahedral connectivity");
                }

                if (s.intma(j, ie) < 1 || s.intma(j, ie) > s.cfg.npoin)
                {
                    throw std::runtime_error(
                        "MeshIO::readMeshFile - connectivity node out of range");
                }
            }
        }

        // ---------------------------------------------------------------------
        // Nodal coordinates:
        //     node_id x y z
        // ---------------------------------------------------------------------
        for (Int i = 1; i <= s.cfg.npoin; ++i)
        {
            Int ip = 0;

            if (!(in >> ip))
            {
                throw std::runtime_error(
                    "MeshIO::readMeshFile - cannot read node index");
            }

            if (ip < 1 || ip > s.cfg.npoin)
            {
                throw std::runtime_error(
                    "MeshIO::readMeshFile - node index out of range");
            }

            for (Int j = 1; j <= s.cfg.ndim; ++j)
            {
                if (!(in >> s.coord(j, ip)))
                {
                    throw std::runtime_error(
                        "MeshIO::readMeshFile - cannot read 3D coordinate");
                }
            }
        }

        // ---------------------------------------------------------------------
        // Boundary faces:
        //     f1 f2 f3 parent_element boundary_side_id
        //
        // Internal legacy-compatible storage:
        //     iside(1:3,ib) = triangular face nodes
        //     iside(4,ib)   = local face number; filled later by preprocess
        //     iside(5,ib)   = parent tetrahedron
        //     iside(6,ib)   = raw boundary-side id, later mapped through .bco
        // ---------------------------------------------------------------------
        for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
        {
            Int n1 = 0;
            Int n2 = 0;
            Int n3 = 0;
            Int parent = 0;
            Int raw_bc = 0;

            if (!(in >> n1 >> n2 >> n3 >> parent >> raw_bc))
            {
                throw std::runtime_error(
                    "MeshIO::readMeshFile - cannot read boundary-face record");
            }

            if (n1 < 1 || n1 > s.cfg.npoin ||
                n2 < 1 || n2 > s.cfg.npoin ||
                n3 < 1 || n3 > s.cfg.npoin)
            {
                throw std::runtime_error(
                    "MeshIO::readMeshFile - boundary-face node out of range");
            }

            if (parent < 1 || parent > s.cfg.nelem)
            {
                throw std::runtime_error(
                    "MeshIO::readMeshFile - boundary-face parent element out of range");
            }

            s.iside(1, ib) = n1;
            s.iside(2, ib) = n2;
            s.iside(3, ib) = n3;
            s.iside(s.cfg.nsidpl, ib) = 0;
            s.iside(s.cfg.nsidpe, ib) = parent;
            s.iside(s.cfg.bsid, ib) = raw_bc;
        }
    }

    // =========================================================================

    void MeshIO::readBoundaryFile(const CaseFiles& files, CBSStateSI& s)
    {
        std::ifstream in(files.bco);
        if (!in)
        {
            throw std::runtime_error(
                "MeshIO::readBoundaryFile - cannot open boundary file: " + files.bco);
        }

        Int nflag = 0;
        Int mass_check = 0;

        {
            const std::string line = next_data_line(in, files.bco);
            std::istringstream iss(line);

            if (!(iss >> nflag))
            {
                throw std::runtime_error(
                    "MeshIO::readBoundaryFile - invalid .bco header");
            }

            if (!(iss >> mass_check))
            {
                mass_check = 0;
            }

            if (nflag != s.cfg.nflag)
            {
                throw std::runtime_error(
                    "MeshIO::readBoundaryFile - boundary flag count mismatch");
            }
        }

        std::unordered_map<Int, Int> raw_to_solver_bc;
        raw_to_solver_bc.reserve(static_cast<std::size_t>(s.cfg.nflag));

        for (Int i = 1; i <= s.cfg.nflag; ++i)
        {
            const std::string line = next_data_line(in, files.bco);
            std::istringstream iss(line);

            Int raw_side_id = 0;
            Int solver_bc = 0;
            Int ignored_group = 0;

            if (!(iss >> raw_side_id >> solver_bc))
            {
                throw std::runtime_error(
                    "MeshIO::readBoundaryFile - cannot read boundary flag record");
            }

            if (mass_check > 0)
            {
                // Compatibility with legacy .bco files:
                //     side_id solver_bc group
                // The current CHT-first CBS3D path does not use group IDs yet.
                iss >> ignored_group;
            }

            if (raw_side_id <= 0)
            {
                throw std::runtime_error(
                    "MeshIO::readBoundaryFile - raw boundary-side ID must be positive");
            }

            if (!supported_solver_bc_id(solver_bc))
            {
                throw std::runtime_error(
                    "MeshIO::readBoundaryFile - unsupported solver BC ID "
                    + std::to_string(solver_bc)
                    + " for raw boundary-side ID "
                    + std::to_string(raw_side_id));
            }

            if (solver_bc == 902)
            {
                std::cout
                    << "WARNING: .bco uses BC 902 on raw boundary-side ID "
                    << raw_side_id
                    << ". In the conformal CHT workflow, 902 must not be used "
                    << "as a fluid-solid interface marker. Use 532 only for "
                    << "external heat-flux walls and material adjacency for "
                    << "the internal interface.\n";
            }

            s.flag_list(1, i) = raw_side_id;
            s.flag_list(2, i) = solver_bc;

            if (!raw_to_solver_bc.emplace(raw_side_id, solver_bc).second)
            {
                throw std::runtime_error(
                    "MeshIO::readBoundaryFile - duplicate raw boundary-side ID in .bco");
            }
        }

        for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
        {
            const Int raw_id = s.iside(s.cfg.bsid, ib);
            const auto it = raw_to_solver_bc.find(raw_id);

            if (it == raw_to_solver_bc.end())
            {
                throw std::runtime_error(
                    "MeshIO::readBoundaryFile - boundary-face raw ID " +
                    std::to_string(raw_id) + " not found in .bco file");
            }

            s.iside(s.cfg.bsid, ib) = it->second;
        }
    }

    // =========================================================================

    void MeshIO::readParameterFile(const CaseFiles& files, CBSStateSI& s)
    {
        std::ifstream in(files.par);
        if (!in)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - cannot open parameter file: " + files.par);
        }

        next_data_line(in, files.par);
        {
            const std::string line = next_data_line(in, files.par);
            std::istringstream iss(line);

            if (!(iss >> s.cfg.solver_opt))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - invalid SOLVER_OPT");
            }
        }

        next_data_line(in, files.par);
        {
            const std::string line = next_data_line(in, files.par);
            std::istringstream iss(line);

            if (!(iss >> s.cfg.restart_opt))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - invalid RESTART_OPT");
            }
        }

        next_data_line(in, files.par);
        {
            const std::string line = next_data_line(in, files.par);
            std::istringstream iss(line);

            if (!(iss >> s.cfg.temp_calc))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - invalid TEMP_CALC");
            }
        }

        next_data_line(in, files.par);
        {
            const std::string line = next_data_line(in, files.par);
            std::istringstream iss(line);

            // 3D free-stream values are:
            //     Ux Uy Uz P T
            if (!(iss >> s.cfg.cinf[0]
                >> s.cfg.cinf[1]
                >> s.cfg.cinf[2]
                >> s.cfg.cinf[3]
                >> s.cfg.cinf[4]))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - invalid 3D freestream line; "
                    "expected Ux Uy Uz P T");
            }
        }

        next_data_line(in, files.par);
        {
            const std::string line = next_data_line(in, files.par);
            std::istringstream iss(line);

            if (!(iss >> s.cfg.ntime
                >> s.cfg.transient_on
                >> s.cfg.dtfixed
                >> s.cfg.dtfix
                >> s.cfg.iwrite))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - invalid timestep-control line");
            }
        }

        next_data_line(in, files.par);
        {
            const std::string line = next_data_line(in, files.par);
            std::istringstream iss(line);

            if (!(iss >> s.cfg.csafm >> s.cfg.theta[0]))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - invalid CSAFM/THETA line");
            }

            s.cfg.theta[1] = 1.0;
        }

        read_cbs3d_timestep_controls(in, s);

        next_data_line(in, files.par);
        {
            const std::string line = next_data_line(in, files.par);
            std::istringstream iss(line);

            if (!(iss >> s.cfg.re >> s.cfg.pr >> s.cfg.ra >> s.cfg.ri))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - invalid Re/Pr/Ra/Ri line");
            }
        }

        next_data_line(in, files.par);
        {
            const std::string line = next_data_line(in, files.par);
            std::istringstream iss(line);

            if (!(iss >> s.cfg.convection_type))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - invalid CONVECTION_TYPE");
            }
        }

        next_data_line(in, files.par);
        {
            const std::string line = next_data_line(in, files.par);
            std::istringstream iss(line);

            if (!(iss >> s.cfg.pnode))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - invalid PNODE");
            }
        }

        next_data_line(in, files.par);
        {
            const std::string line = next_data_line(in, files.par);
            std::istringstream iss(line);

            if (!(iss >> s.cfg.write_output
                >> s.cfg.write_time_output
                >> s.cfg.time_output_interval
                >> s.cfg.end_rtime))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - invalid output-control line");
            }
        }

        next_data_line(in, files.par);
        {
            const std::string line = next_data_line(in, files.par);
            std::istringstream iss(line);

            if (!(iss >> s.cfg.relToler >> s.cfg.absToler))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - invalid CG tolerance line");
            }
        }

        next_data_line(in, files.par);
        {
            const std::string line = next_data_line(in, files.par);
            std::istringstream iss(line);

            if (!(iss >> s.cfg.vel_check
                >> s.cfg.l2norm_vel_tolerance
                >> s.cfg.temp_check
                >> s.cfg.l2norm_temp_tolerance))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - invalid steady-state tolerance line");
            }

            // Optional trailing pair: sa_check and l2norm_sa_tolerance.
            //
            // Read only when present so that existing parameter files stay
            // valid and keep the previous stopping behaviour, in which the SA
            // residual is reported but does not gate the stop.
            Int sa_check = 0;
            Real sa_tolerance = 0.0;

            if (iss >> sa_check)
            {
                if (sa_check < 0 || sa_check > 1)
                {
                    throw std::runtime_error(
                        "MeshIO::readParameterFile - sa_check must be 0 or 1");
                }

                s.cfg.sa_check = sa_check;

                if (iss >> sa_tolerance)
                {
                    if (!(sa_tolerance > 0.0) || !std::isfinite(sa_tolerance))
                    {
                        throw std::runtime_error(
                            "MeshIO::readParameterFile - l2norm_sa_tolerance "
                            "must be positive");
                    }

                    s.cfg.l2norm_sa_tolerance = sa_tolerance;
                }
                else if (sa_check > 0)
                {
                    throw std::runtime_error(
                        "MeshIO::readParameterFile - sa_check=1 requires "
                        "l2norm_sa_tolerance on the same line");
                }
            }
        }

        next_data_line(in, files.par);
        {
            const std::string line = next_data_line(in, files.par);
            std::istringstream iss(line);

            if (!(iss >> s.cfg.paraview_output
                >> s.cfg.tecplot_output
                >> s.cfg.nusselt_calc
                >> s.cfg.nusselt_flag))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - invalid output-format line");
            }
        }

        next_data_line(in, files.par);
        {
            const std::string line = next_data_line(in, files.par);
            std::istringstream iss(line);

            if (!(iss >> s.cfg.runtime_mod))
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - invalid RUNTIME_MOD");
            }
        }

        read_optional_alpha_source_flux(in, s);
        read_optional_dimensional_controls(in, s);
        read_pressure_cg_controls(in, s);
        read_output_monitor_controls(in, s);
        read_artificial_diffusion_control(in, s);
        read_optional_spalart_allmaras_controls(in, s);

        if (s.cfg.solver_opt < 1 || s.cfg.solver_opt > 3)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - SOLVER_OPT must be 1, 2, or 3");
        }

        if (s.cfg.restart_opt < 0 || s.cfg.restart_opt > 1)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - RESTART_OPT must be 0 or 1");
        }

        if (s.cfg.temp_calc < 0 || s.cfg.temp_calc > 1)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - TEMP_CALC must be 0 or 1");
        }

        if (s.cfg.ntime < 1)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - NTIME must be >= 1");
        }

        if (s.cfg.transient_on < 0 || s.cfg.transient_on > 1)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - TRANSIENT_ON must be 0 or 1");
        }

        if (s.cfg.dtfixed < -1 || s.cfg.dtfixed > 1)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - DTFIXED must be -1, 0 or 1");
        }

        if (s.cfg.dtfix <= 0.0)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - DTFIX must be > 0");
        }

        if (s.cfg.iwrite < 1)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - IWRITE must be >= 1");
        }

        if (s.cfg.csafm <= 0.0)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - CSAFM must be > 0");
        }

        if (s.cfg.theta[0] < 0.5 || s.cfg.theta[0] > 1.0)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - THETA(1) must satisfy 0.5 <= theta <= 1.0");
        }

        if (s.cfg.cbs_scheme < 0 || s.cfg.cbs_scheme > 1)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - CBS_SCHEME must be 0 explicit or 1 semi-implicit");
        }

        if (s.cfg.htype < 1 || s.cfg.htype > 3)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - HTYPE must be 1, 2, or 3");
        }

        if (s.cfg.step2_check < 0 || s.cfg.step2_check > 1)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - STEP2_CHECK must be 0 or 1");
        }

        if (s.cfg.rem_deltp < 0 || s.cfg.rem_deltp > 1)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - REM_DELTP must be 0 or 1");
        }

        if (s.cfg.beta_opt < 0 || s.cfg.beta_opt > 1)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - BETA_OPT must be 0 or 1");
        }

        if (s.cfg.dtfix_end < 0)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - DTFIX_END must be >= 0");
        }

        if (s.cfg.csafm2 <= 0.0)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - CSAFM2 must be > 0");
        }

        if (s.cfg.epsilon1 <= 0.0)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - EPSILON1 must be > 0");
        }

        if (s.cfg.deltr <= 0.0)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - DELTR must be > 0");
        }

        if (s.cfg.cg_preconditioner < 0 || s.cfg.cg_preconditioner > 1)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - CG_PRECONDITIONER must be 0 or 1");
        }

        if (s.cfg.cg_conv_test < 1 || s.cfg.cg_conv_test > 3)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - CG_CONV_TEST must be 1, 2, or 3");
        }

        if (s.cfg.cg_max_iter < 1)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - CG_MAX_ITER must be >= 1");
        }

        if (s.cfg.l2norm_pres_tolerance < 0.0)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - L2NORM_PRES_TOLERANCE must be >= 0");
        }

        if (s.cfg.residual_log_enabled < 0 || s.cfg.residual_log_enabled > 1)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - RESIDUAL_LOG_ENABLED must be 0 or 1");
        }

        if (s.cfg.residual_log_every < 1)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - RESIDUAL_LOG_EVERY must be >= 1");
        }

        if (s.cfg.console_log_every < 0)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - CONSOLE_LOG_EVERY must be >= 0");
        }

        if (s.cfg.live_residual_plot < 0 || s.cfg.live_residual_plot > 1)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - LIVE_RESIDUAL_PLOT must be 0 or 1");
        }

        if (s.cfg.vtu_output_enabled < 0 || s.cfg.vtu_output_enabled > 1)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - VTU_OUTPUT_ENABLED must be 0 or 1");
        }

        if (s.cfg.vtu_output_every_iterations < 0)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - VTU_OUTPUT_EVERY_ITERATIONS must be >= 0");
        }

        if (s.cfg.vtu_output_every_sim_time < 0.0)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - VTU_OUTPUT_EVERY_SIM_TIME must be >= 0");
        }

        if (s.cfg.write_boundary_debug_arrays < 0 || s.cfg.write_boundary_debug_arrays > 1)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - WRITE_BOUNDARY_DEBUG_ARRAYS must be 0 or 1");
        }

        if (s.cfg.steady_min_iterations < 1)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - STEADY_MIN_ITERATIONS must be >= 1");
        }

        if (s.cfg.art_diff < 0 || s.cfg.art_diff > 1)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - ART_DIFF must be 0 or 1");
        }

        if (s.cfg.convection_type < 0 || s.cfg.convection_type > 1)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - CONVECTION_TYPE must be 0 or 1");
        }

        if (s.cfg.pnode < 1 || s.cfg.pnode > s.cfg.npoin)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - PNODE out of range");
        }

        if (s.cfg.model_depth <= 0.0)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - MODEL_DEPTH must be > 0");
        }

        if (s.cfg.material_properties_enabled < 0 ||
            s.cfg.material_properties_enabled > 1)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - MATERIAL_PROPERTIES_ENABLED must be 0 or 1");
        }

        if (s.cfg.mass_flow_inlet_enabled < 0 ||
            s.cfg.mass_flow_inlet_enabled > 1)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - MASS_FLOW_INLET_ENABLED must be 0 or 1");
        }

        if (s.cfg.inlet_density <= 0.0)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - INLET_DENSITY must be > 0");
        }

        if (s.cfg.nusselt_diameter <= 0.0)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - NUSSELT_DIAMETER must be > 0");
        }


        if (s.cfg.turbulence_on < 0 || s.cfg.turbulence_on > 1)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - TURBULENCE_ON must be 0 or 1");
        }

        // turbulence_model selects the SA variant:
        //
        //     0  standard Spalart-Allmaras
        //     1  SA-neg, reserved
        //
        // SA-neg is not implemented.  The transported variable is clipped at
        // sa_nu_tilde_floor rather than following the negative branch of the
        // model, so accepting turbulence_model = 1 would silently run standard
        // SA under an SA-neg label and invalidate any result reported as SA-neg.
        if (s.cfg.turbulence_model != 0)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - turbulence_model must be 0."
                " Value 1 is reserved for SA-neg, which is not implemented:"
                " the solver clips nu_tilde at sa_nu_tilde_floor instead of"
                " integrating the negative branch of the model.");
        }

        if (s.cfg.turbulent_thermal_diffusivity_on < 0 ||
            s.cfg.turbulent_thermal_diffusivity_on > 1)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - TURBULENT_THERMAL_DIFFUSIVITY_ON must be 0 or 1");
        }

        if (s.cfg.sa_inlet_ratio < 0.0)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - SA_INLET_RATIO must be non-negative");
        }

        if (s.cfg.sa_prandtl_t <= 0.0)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - SA_PRANDTL_T must be positive");
        }

        if (s.cfg.sa_min_wall_distance <= 0.0)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - SA_MIN_WALL_DISTANCE must be positive");
        }

        if (s.cfg.sa_min_stilde <= 0.0)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - SA_MIN_STILDE must be positive");
        }

        if (s.cfg.sa_nu_tilde_floor < 0.0)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - SA_NU_TILDE_FLOOR must be non-negative");
        }

        if (s.cfg.sa_use_stilde_limiter < 0 || s.cfg.sa_use_stilde_limiter > 1)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - SA_USE_STILDE_LIMITER must be 0 or 1");
        }

        if (s.cfg.sa_implicit_destruction < 0 || s.cfg.sa_implicit_destruction > 1)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - SA_IMPLICIT_DESTRUCTION must be 0 or 1");
        }

        if (s.cfg.turbulence_on > 0)
        {
            if (s.cfg.dimensional_mode < 1 || s.cfg.material_properties_enabled < 1)
            {
                throw std::runtime_error(
                    "MeshIO::readParameterFile - Spalart-Allmaras requires dimensional_mode=1 "
                    "and material_properties_enabled=1 so that nu = mu/rho is available");
            }
        }

        s.cfg.istart = 1;
        s.cfg.rtime = 0.0;
    }

    // =========================================================================

    void MeshIO::readMaterialFile(const CaseFiles& files, CBSStateSI& s)
    {
        const std::vector<std::string> lines = read_data_lines(files.material);

        if (lines.empty())
        {
            throw std::runtime_error(
                "MeshIO::readMaterialFile - empty material file: " + files.material);
        }

        Size first_record = 0;

        {
            std::istringstream iss(lines[0]);
            std::vector<Int> vals;
            Int v = 0;

            while (iss >> v)
            {
                vals.push_back(v);
            }

            if (vals.size() == 1)
            {
                if (vals[0] != s.cfg.nelem)
                {
                    throw std::runtime_error(
                        "MeshIO::readMaterialFile - material header NELEM mismatch");
                }

                first_record = 1;
            }
        }

        const Size expected_records = static_cast<Size>(s.cfg.nelem);

        if (lines.size() - first_record != expected_records)
        {
            throw std::runtime_error(
                "MeshIO::readMaterialFile - number of material records does not match NELEM");
        }

        for (Int i = 1; i <= s.cfg.nelem; ++i)
        {
            const std::string& line = lines[first_record + static_cast<Size>(i - 1)];
            std::istringstream iss(line);

            Int ie = 0;
            Int n1 = 0;
            Int n2 = 0;
            Int n3 = 0;
            Int n4 = 0;
            Int mat = 0;

            if (!(iss >> ie >> n1 >> n2 >> n3 >> n4 >> mat))
            {
                throw std::runtime_error(
                    "MeshIO::readMaterialFile - invalid material record; "
                    "expected eid n1 n2 n3 n4 material_id");
            }

            if (ie < 1 || ie > s.cfg.nelem)
            {
                throw std::runtime_error(
                    "MeshIO::readMaterialFile - element index out of range");
            }

            if (s.intma(1, ie) != n1 ||
                s.intma(2, ie) != n2 ||
                s.intma(3, ie) != n3 ||
                s.intma(4, ie) != n4)
            {
                throw std::runtime_error(
                    "MeshIO::readMaterialFile - connectivity mismatch with mesh file");
            }

            s.mat_elem(ie) = mat;
        }
    }

    // =========================================================================

    void MeshIO::readMaterialPropertyFile(const CaseFiles& files, CBSStateSI& s)
    {
        std::ifstream in(files.matprop);

        if (!in)
        {
            throw std::runtime_error(
                "MeshIO::readMaterialPropertyFile - cannot open .matprop file: " +
                files.matprop);
        }

        Int nmat = 0;

        if (!(in >> nmat))
        {
            throw std::runtime_error(
                "MeshIO::readMaterialPropertyFile - cannot read material count");
        }

        if (nmat < 1)
        {
            throw std::runtime_error(
                "MeshIO::readMaterialPropertyFile - invalid number of materials");
        }

        struct MatData
        {
            Real rho = 0.0;
            Real cp = 0.0;
            Real k = 0.0;
            Real mu = 0.0;
            Real Qvol = 0.0;
        };

        std::unordered_map<Int, MatData> mat;
        mat.reserve(static_cast<std::size_t>(nmat));

        for (Int i = 1; i <= nmat; ++i)
        {
            Int id = -1;
            std::string name;
            std::string phase;

            MatData data;

            if (!(in >> id
                >> name
                >> phase
                >> data.rho
                >> data.cp
                >> data.k
                >> data.mu
                >> data.Qvol))
            {
                throw std::runtime_error(
                    "MeshIO::readMaterialPropertyFile - invalid .matprop record");
            }

            check_material_property_values(id, data.rho, data.cp, data.k, data.mu);

            if (!mat.emplace(id, data).second)
            {
                throw std::runtime_error(
                    "MeshIO::readMaterialPropertyFile - duplicate material ID");
            }
        }

        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            const Int mid = s.mat_elem(ie);
            const auto it = mat.find(mid);

            if (it == mat.end())
            {
                throw std::runtime_error(
                    "MeshIO::readMaterialPropertyFile - element " +
                    std::to_string(ie) + " uses material ID " +
                    std::to_string(mid) + " missing from .matprop");
            }

            const MatData& m = it->second;

            s.rho_e(ie) = m.rho;
            s.cp_e(ie) = m.cp;
            s.k_e(ie) = m.k;
            s.mu_e(ie) = m.mu;
            s.Qvol_e(ie) = m.Qvol;

            s.rho_cp_e(ie) = s.rho_e(ie) * s.cp_e(ie);

            if (s.rho_cp_e(ie) <= 0.0)
            {
                throw std::runtime_error(
                    "MeshIO::readMaterialPropertyFile - invalid rho*cp in element");
            }

            s.alpha_e(ie) = s.k_e(ie) / s.rho_cp_e(ie);
        }

        std::cout << "Material properties loaded successfully\n";
    }

    // =========================================================================

    void MeshIO::initialiseFields(CBSStateSI& s)
    {
        s.unkno.fill(0.0);
        s.unkn1.fill(0.0);

        s.pres.fill(0.0);
        s.pres1.fill(0.0);

        s.temperature.fill(0.0);
        s.temperature1.fill(0.0);

        s.unkno_temp.fill(0.0);
        s.pres_temp.fill(0.0);
        s.temperature_temp.fill(0.0);

        s.unknn1.fill(0.0);
        s.unknn2.fill(0.0);
        s.tempert1.fill(0.0);
        s.tempert2.fill(0.0);

        set_freestream_initial_condition(s);

        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            s.unkn1(1, ip) = s.unkno(1, ip);
            s.unkn1(2, ip) = s.unkno(2, ip);
            s.unkn1(3, ip) = s.unkno(3, ip);

            s.pres1(ip) = s.pres(ip);
            s.temperature1(ip) = s.temperature(ip);
        }

        enforce_initial_material_domain_fields(s);

        s.cfg.rtime_output = s.cfg.rtime + s.cfg.time_output_interval;

        if (s.cfg.convection_type == 1)
        {
            s.cfg.ani = s.cfg.pr;
            s.cfg.temp_calc = 1;
        }
        else
        {
            s.cfg.ani = 1.0 / s.cfg.re;

            if (s.cfg.ri > 0.0)
            {
                s.cfg.temp_calc = 1;
            }
        }

        if (s.cfg.material_properties_enabled < 1)
        {
            // Now that ani has been derived from Re/Pr/convection mode, refresh
            // the default material fields with a consistent nondimensional
            // fluid viscosity.  These values are primarily safety/default data;
            // dimensional CHT uses explicit .material/.matprop files.
            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                s.mat_elem(ie) = 0;
                s.rho_e(ie) = 1.0;
                s.cp_e(ie) = 1.0;
                s.rho_cp_e(ie) = 1.0;
                s.mu_e(ie) = s.cfg.ani;
                s.k_e(ie) = 1.0;
                s.alpha_e(ie) = 1.0;
                s.Qvol_e(ie) = 0.0;

                s.mu_eff_e(ie) = s.mu_e(ie);
                s.k_eff_e(ie) = s.k_e(ie);
                s.nu_tilde_e(ie) = 0.0;
                s.nu_t_e(ie) = 0.0;
                s.mu_t_e(ie) = 0.0;
            }
        }

        if (s.cfg.temp_calc < 1)
        {
            s.cfg.nusselt_calc = 0;
        }
    }

}
