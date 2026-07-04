//=============================================================================
// CBS3D++_SI
//
// Post-processing, terminal monitoring and solution-output control.
//
// The terminal output is intentionally ASCII-only for reliable use in Windows
// consoles, MSYS/UCRT64 terminals, SSH sessions and HPC batch logs.
//
// The module writes:
//
//     residual CSV files
//     VTU unstructured-grid solution files
//     PVD time-series collection files
//
// It also provides the console progress display and the optional live residual
// plotting hook.
//=============================================================================

#include "cbs/io/Post.hpp"

#include "cbs/solver/Convergence.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace cbs
{
    namespace
    {
        namespace fs = std::filesystem;
        using Clock = std::chrono::steady_clock;

        Clock::time_point run_wall_clock_start = Clock::now();

        constexpr Int panel_width = 96;
        constexpr Int label_width = 28;

        constexpr const char* reset   = "\033[0m";
        constexpr const char* dim     = "\033[2m";
        constexpr const char* bold    = "\033[1m";
        constexpr const char* cyan    = "\033[36m";
        constexpr const char* green   = "\033[32m";
        constexpr const char* yellow  = "\033[33m";
        constexpr const char* magenta = "\033[35m";
        constexpr const char* red     = "\033[31m";
        constexpr const char* blue    = "\033[34m";

        // Returns the root directory used by all solver output.
        fs::path output_root()
        {
            return fs::path("output");
        }

        // Returns the case-specific output directory:
        //
        //     output/<case_name>
        fs::path case_output_dir(const std::string& case_name)
        {
            return output_root() / case_name;
        }

        // Returns the residual CSV file path.
        fs::path residual_file_path(const std::string& case_name)
        {
            return case_output_dir(case_name) / (case_name + "_residuals.csv");
        }

        // Returns the ParaView collection-file path.
        fs::path pvd_file_path(const std::string& case_name)
        {
            return case_output_dir(case_name) / (case_name + ".pvd");
        }

        // Converts the iteration number to a fixed-width solution tag:
        //
        //     step_00000042
        std::string step_tag(const Int iitime)
        {
            std::ostringstream oss;
            oss << "step_" << std::setw(8) << std::setfill('0') << iitime;
            return oss.str();
        }

        // Returns the VTU file name for one CBS iteration.
        std::string vtu_file_name(
            const std::string& case_name,
            const Int iitime)
        {
            return case_name + "_" + step_tag(iitime) + ".vtu";
        }

        // Returns the full path of one VTU solution file.
        fs::path vtu_file_path(
            const std::string& case_name,
            const Int iitime)
        {
            return case_output_dir(case_name) / vtu_file_name(case_name, iitime);
        }

        // Creates the case output directory when it does not already exist.
        void ensure_output_directory(const std::string& case_name)
        {
            fs::create_directories(case_output_dir(case_name));
        }

        // Replaces NaN or infinite output values by zero so that CSV and VTK
        // files remain readable.
        Real safe_value(const Real value)
        {
            return std::isfinite(value) ? value : 0.0;
        }

        // Formats one scalar in scientific notation for terminal output.
        std::string sci(const Real value, const Int precision = 2)
        {
            std::ostringstream oss;
            oss << std::scientific << std::setprecision(precision) << safe_value(value);
            return oss.str();
        }

        // Converts an integer control flag to ON or OFF.
        std::string on_off(const Int flag)
        {
            return flag > 0 ? "ON" : "OFF";
        }

        // Returns a readable name for the selected pressure solver.
        std::string pressure_solver_name(const Int solver_opt)
        {
            if (solver_opt == 1)
            {
                return "Conjugate Gradient";
            }

            if (solver_opt == 2)
            {
                return "Banded Gaussian";
            }

            if (solver_opt == 3)
            {
                return "PETSc KSP + AMG";
            }

            return "Unknown";
        }

        // Returns a readable name for the selected CBS formulation.
        std::string scheme_name(const Int cbs_scheme)
        {
            if (cbs_scheme == 1)
            {
                return "Semi-Implicit CBS";
            }

            if (cbs_scheme == 0)
            {
                return "Fully Explicit CBS";
            }

            return "Unknown";
        }

        // Returns the elapsed wall-clock time since output initialisation.
        Real elapsed_wall_seconds()
        {
            return std::chrono::duration<Real>(Clock::now() - run_wall_clock_start).count();
        }

        // Formats elapsed time as seconds, minutes or hours.
        std::string format_seconds(const Real seconds)
        {
            std::ostringstream oss;

            if (seconds < 60.0)
            {
                oss << std::fixed << std::setprecision(1) << seconds << "s";
                return oss.str();
            }

            const auto total = static_cast<long long>(seconds);
            const long long h = total / 3600;
            const long long m = (total % 3600) / 60;
            const long long s = total % 60;

            if (h > 0)
            {
                oss << h << "h " << m << "m " << s << "s";
            }
            else
            {
                oss << m << "m " << s << "s";
            }

            return oss.str();
        }

        // Converts a completion fraction in [0,1] to a percentage string.
        std::string percent_string(const Real fraction)
        {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(1)
                << 100.0 * std::clamp(fraction, 0.0, 1.0)
                << "%";
            return oss.str();
        }

        // Builds the fixed-width ASCII progress bar used by the terminal
        // monitor.
        std::string progress_bar(
            const Real fraction,
            const Int width)
        {
            const Real f = std::clamp(fraction, 0.0, 1.0);
            const Int filled =
                static_cast<Int>(std::round(f * static_cast<Real>(width)));

            std::string bar;
            bar.reserve(static_cast<std::size_t>(width));

            for (Int i = 0; i < width; ++i)
            {
                if (i < filled)
                {
                    bar.push_back('#');
                }
                else if (i == filled && filled < width)
                {
                    bar.push_back('>');
                }
                else
                {
                    bar.push_back('-');
                }
            }

            return bar;
        }

        // Selects one frame of the ASCII progress spinner.
        const char* spinner(const Int iitime)
        {
            static constexpr std::array<const char*, 4> frames =
            {
                "|", "/", "-", "\\"
            };

            return frames[static_cast<std::size_t>(iitime) % frames.size()];
        }

        // Shortens text to the available panel width.
        std::string clip(
            const std::string& input,
            const std::size_t max_len)
        {
            if (input.size() <= max_len)
            {
                return input;
            }

            if (max_len <= 3)
            {
                return input.substr(0, max_len);
            }

            return input.substr(0, max_len - 3) + "...";
        }

        // Pads one text field to a fixed terminal width.
        std::string pad_right(
            const std::string& input,
            const std::size_t width)
        {
            const std::string clipped = clip(input, width);

            if (clipped.size() >= width)
            {
                return clipped;
            }

            return clipped + std::string(width - clipped.size(), ' ');
        }

        // Prints one horizontal panel border.
        void hr(const char* colour = cyan)
        {
            std::cout
                << colour
                << "+"
                << std::string(panel_width - 2, '-')
                << "+"
                << reset
                << "\n";
        }

        // Prints the lower panel border.
        void footer(const char* colour = cyan)
        {
            hr(colour);
        }

        // Prints an internal panel separator.
        void separator_row(const char* colour = cyan)
        {
            std::cout
                << colour
                << "+"
                << std::string(panel_width - 2, '-')
                << "+"
                << reset
                << "\n";
        }

        // Prints one centred row inside the terminal panel.
        void centre_row(
            const std::string& text,
            const char* colour = cyan)
        {
            const std::size_t inner = static_cast<std::size_t>(panel_width - 2);
            const std::string clipped = clip(text, inner);
            const std::size_t left =
                (inner > clipped.size()) ? (inner - clipped.size()) / 2 : 0;
            const std::size_t right = inner - left - clipped.size();

            std::cout
                << cyan << "|" << reset
                << colour
                << std::string(left, ' ')
                << clipped
                << std::string(right, ' ')
                << reset
                << cyan << "|" << reset
                << "\n";
        }

        // Prints one labelled value inside the terminal panel.
        void key_value_row(
            const std::string& key,
            const std::string& value,
            const char* key_colour = dim)
        {
            const std::size_t inner = static_cast<std::size_t>(panel_width - 2);
            const std::string prefix =
                "  " + pad_right(key, label_width) + " ";

            const std::size_t available_value_width =
                inner > prefix.size() ? inner - prefix.size() : 0;

            const std::string clipped_value =
                clip(value, available_value_width);

            std::cout
                << cyan << "|" << reset
                << key_colour << prefix << reset
                << clipped_value
                << std::string(
                    available_value_width > clipped_value.size()
                        ? available_value_width - clipped_value.size()
                        : 0,
                    ' ')
                << cyan << "|" << reset
                << "\n";
        }

        // Quotes one path before it is passed to an external shell command.
        std::string quote_path_for_command(const fs::path& p)
        {
            std::ostringstream oss;
            oss << '"' << p.string() << '"';
            return oss.str();
        }

        // Builds a diagnostic boundary-condition identifier for every node.
        //
        // When a node belongs to several boundary faces, the identifier from
        // the last visited face is retained. This array is intended only for
        // ParaView debugging and not for numerical boundary enforcement.
        std::vector<Int> build_node_bc_debug(const CBSStateSI& s)
        {
            std::vector<Int> node_bc(static_cast<std::size_t>(s.cfg.npoin + 1), 0);

            for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
            {
                const Int bc = s.iside(s.cfg.bsid, ib);

                for (Int in = 1; in <= s.cfg.nsidp; ++in)
                {
                    const Int ip = s.iside(in, ib);

                    if (ip >= 1 && ip <= s.cfg.npoin)
                    {
                        node_bc[static_cast<std::size_t>(ip)] = bc;
                    }
                }
            }

            return node_bc;
        }

        // Nodal material-domain classification used when writing physical
        // fields to VTU.
        struct NodePhysicsMasks
        {
            std::vector<char> touches_fluid;
            std::vector<char> touches_solid;
            std::vector<char> pressure_active;
            std::vector<char> velocity_active;
            std::vector<char> thermal_active;
            std::vector<Int> node_domain_kind;
        };

        // Builds nodal fluid/solid connectivity masks.
        //
        // Domain codes written to ParaView:
        //
        //     0  orphan/unclassified node
        //     1  fluid-only node
        //     2  solid-only node
        //     3  conformal fluid-solid interface node
        //
        // Field activity:
        //
        //     pressure_active = touches fluid
        //
        //     velocity_active = touches fluid and does not touch solid
        //
        //     thermal_active  = touches fluid or solid
        NodePhysicsMasks build_node_physics_masks(const CBSStateSI& s)
        {
            NodePhysicsMasks masks;
            const std::size_t n = static_cast<std::size_t>(s.cfg.npoin + 1);

            masks.touches_fluid.assign(n, 0);
            masks.touches_solid.assign(n, 0);
            masks.pressure_active.assign(n, 0);
            masks.velocity_active.assign(n, 0);
            masks.thermal_active.assign(n, 0);
            masks.node_domain_kind.assign(n, 0);

            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                const bool fluid_element = (s.mat_elem(ie) == 0);

                for (Int in = 1; in <= s.cfg.nep; ++in)
                {
                    const Int ip = s.intma(in, ie);

                    if (ip < 1 || ip > s.cfg.npoin)
                    {
                        throw std::runtime_error(
                            "Post::writeVTU - element connectivity node out of range");
                    }

                    if (fluid_element)
                    {
                        masks.touches_fluid[static_cast<std::size_t>(ip)] = 1;
                    }
                    else
                    {
                        masks.touches_solid[static_cast<std::size_t>(ip)] = 1;
                    }
                }
            }

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                const bool fluid =
                    masks.touches_fluid[static_cast<std::size_t>(ip)] != 0;
                const bool solid =
                    masks.touches_solid[static_cast<std::size_t>(ip)] != 0;

                masks.pressure_active[static_cast<std::size_t>(ip)] =
                    fluid ? 1 : 0;

                masks.velocity_active[static_cast<std::size_t>(ip)] =
                    (fluid && !solid) ? 1 : 0;

                masks.thermal_active[static_cast<std::size_t>(ip)] =
                    (fluid || solid) ? 1 : 0;

                if (fluid && solid)
                {
                    masks.node_domain_kind[static_cast<std::size_t>(ip)] = 3;
                }
                else if (fluid)
                {
                    masks.node_domain_kind[static_cast<std::size_t>(ip)] = 1;
                }
                else if (solid)
                {
                    masks.node_domain_kind[static_cast<std::size_t>(ip)] = 2;
                }
                else
                {
                    masks.node_domain_kind[static_cast<std::size_t>(ip)] = 0;
                }
            }

            return masks;
        }

        // Tests whether pressure is physically active at one node.
        bool pressure_active_node(
            const NodePhysicsMasks& masks,
            const Int ip)
        {
            return masks.pressure_active[static_cast<std::size_t>(ip)] != 0;
        }

        // Tests whether velocity is physically active at one node.
        bool velocity_active_node(
            const NodePhysicsMasks& masks,
            const Int ip)
        {
            return masks.velocity_active[static_cast<std::size_t>(ip)] != 0;
        }

        // Creates the residual CSV file and writes its column names.
        //
        // The file contains relative, field-norm and absolute/RHS-scale
        // residuals for velocity, pressure and temperature, followed by the
        // Pressure CG convergence data.
        void write_residual_header(const std::string& case_name)
        {
            std::ofstream out(residual_file_path(case_name), std::ios::trunc);

            if (!out)
            {
                throw std::runtime_error(
                    "Post::initialiseRunOutputs - cannot open residual CSV file");
            }

            out
                << "iteration,"
                << "time,"
                << "dt,"
                << "u_rel,u_norm,u_abs,"
                << "v_rel,v_norm,v_abs,"
                << "w_rel,w_norm,w_abs,"
                << "p_rel,p_norm,p_abs,"
                << "T_rel,T_norm,T_abs,"
                << "velocity_rel_max,"
                << "cg_iterations,"
                << "cg_initial_l2,"
                << "cg_final_l2,"
                << "cg_relative_l2,"
                << "cg_max_abs\n";
        }

        // Returns true when the iteration-based VTU output interval has been
        // reached.
        bool write_by_iteration(
            const CBSStateSI& s,
            const Int iitime)
        {
            return s.cfg.vtu_output_every_iterations > 0 &&
                   (iitime % s.cfg.vtu_output_every_iterations == 0);
        }

        // Returns true when the transient simulated-time output interval has
        // been reached.
        //
        // next_vtu_output_time is advanced beyond the current physical time so
        // that large time steps do not repeatedly trigger the same interval.
        bool write_by_simulated_time(
            CBSStateSI& s,
            const Int iitime)
        {
            if (s.cfg.transient_on < 1 ||
                s.cfg.vtu_output_every_sim_time <= 0.0)
            {
                return false;
            }

            if (iitime == 0)
            {
                s.next_vtu_output_time = s.cfg.vtu_output_every_sim_time;
                return true;
            }

            if (s.cfg.rtime + 1.0e-14 < s.next_vtu_output_time)
            {
                return false;
            }

            while (s.next_vtu_output_time <= s.cfg.rtime + 1.0e-14)
            {
                s.next_vtu_output_time += s.cfg.vtu_output_every_sim_time;
            }

            return true;
        }
    }

    //=========================================================================
    // Prints the solver title and the principal capabilities of the current
    // executable.
    //=========================================================================
    void Post::printBanner()
    {
        std::cout << "\n";
        std::cout << bold << cyan;
        std::cout << "   ____ ____ ____  _____ ____    _     ____ ___\n";
        std::cout << "  / ___| __ ) ___||___ /|  _ \\ _| |_  / ___|_ _|\n";
        std::cout << " | |   |  _ \\___ \\  |_ \\| | | |_   _| \\___ \\| |\n";
        std::cout << " | |___| |_) |__) |___) | |_| | |_|    ___) | |\n";
        std::cout << "  \\____|____/____/|____/|____/        |____/___|\n";
        std::cout << reset << "\n";

        hr();
        centre_row("CBS3D++_SI - Semi-Implicit Characteristic-Based Split Solver", bold);
        separator_row();
        centre_row("3D tetrahedral incompressible flow | optional conjugate heat transfer", reset);
        centre_row("VTU/PVD output | residual CSV | live residual plot hook", reset);
        footer();
        std::cout << "\n";
    }

    //=========================================================================
    // Prints the numerical and output configuration before the CBS loop starts.
    //
    // Reported information includes mesh size, CBS formulation, pressure
    // solver, transient and energy controls, output intervals and target number
    // of CBS iterations.
    //=========================================================================
    void Post::printRunSummary(
        const CBSStateSI& s,
        const std::string& case_name)
    {
        hr(magenta);
        centre_row("RUN DASHBOARD", bold);
        separator_row(magenta);

        key_value_row("Case", case_name);
        key_value_row(
            "Mesh",
            std::to_string(s.cfg.nelem) + " tets | "
            + std::to_string(s.cfg.npoin) + " nodes | "
            + std::to_string(s.cfg.nboun) + " boundary faces");
        key_value_row("Pressure solver", pressure_solver_name(s.cfg.solver_opt));
        key_value_row("CBS scheme", scheme_name(s.cfg.cbs_scheme));
        key_value_row("Transient mode", on_off(s.cfg.transient_on));
        key_value_row("Energy equation", on_off(s.cfg.temp_calc));
        key_value_row(
            "Material mode",
            s.cfg.material_properties_enabled > 0
                ? "material files required"
                : "flow-only all-fluid default");
        key_value_row(
            "Residual CSV",
            on_off(s.cfg.residual_log_enabled) + std::string(" | every ")
            + std::to_string(s.cfg.residual_log_every) + " iter");
        key_value_row(
            "Console monitor",
            "every " + std::to_string(s.cfg.console_log_every) + " iter");
        key_value_row(
            "VTU/PVD output",
            on_off(s.cfg.vtu_output_enabled) + std::string(" | every ")
            + std::to_string(s.cfg.vtu_output_every_iterations) + " iter");
        key_value_row("Live plotter", on_off(s.cfg.live_residual_plot));
        key_value_row("Target", std::to_string(s.cfg.ntime) + " iterations");
        footer(magenta);

        std::cout
            << dim
            << "  Residual legend: RelU/RelV/RelW/RelP/RelT are convergence residuals, not field values.\n\n"
            << reset;
    }

    //=========================================================================
    // Prints the beginning of one input, preprocessing or solver stage.
    //=========================================================================
    void Post::printStage(
        const std::string& stage,
        const std::string& detail)
    {
        std::cout
            << "  " << blue << "[..]" << reset << " "
            << bold << pad_right(stage, 30) << reset;

        if (!detail.empty())
        {
            std::cout << dim << detail << reset;
        }

        std::cout << "\n";
    }

    //=========================================================================
    // Prints successful completion of one stage.
    //=========================================================================
    void Post::printStageDone(
        const std::string& stage,
        const std::string& detail)
    {
        std::cout
            << "  " << green << "[OK]" << reset << " "
            << bold << pad_right(stage, 30) << reset;

        if (!detail.empty())
        {
            std::cout << green << detail << reset;
        }

        std::cout << "\n";
    }

    //=========================================================================
    // Initialises all run-time output.
    //
    // The routine:
    //
    //     1. resets the wall-clock timer;
    //     2. creates the case output directory;
    //     3. clears the VTU/PVD output history;
    //     4. creates the residual CSV file when enabled;
    //     5. writes the initial solution at CBS iteration zero;
    //     6. launches the optional live residual monitor.
    //=========================================================================
    void Post::initialiseRunOutputs(
        CBSStateSI& s,
        const std::string& case_name)
    {
        run_wall_clock_start = Clock::now();

        ensure_output_directory(case_name);

        s.output_iterations.clear();
        s.output_times.clear();
        s.next_vtu_output_time = s.cfg.vtu_output_every_sim_time;

        if (s.cfg.residual_log_enabled > 0)
        {
            write_residual_header(case_name);
        }

        if (s.cfg.vtu_output_enabled > 0)
        {
            writeSolution(s, case_name, 0);
        }

        launchLiveResidualPlotterIfRequested(s, case_name);
    }

    //=========================================================================
    // Appends one residual record to the case CSV file.
    //
    // The row contains:
    //
    //     iteration, physical time and time step
    //     u, v and w residual triplets
    //     pressure residual triplet
    //     temperature residual triplet
    //     maximum relative velocity residual
    //     Pressure CG iterations and residual measures
    //=========================================================================
    void Post::writeResidualRow(
        const CBSStateSI& s,
        const std::string& case_name,
        const Int iitime)
    {
        if (s.cfg.residual_log_enabled < 1)
        {
            return;
        }

        if (s.cfg.residual_log_every > 1 &&
            (iitime % s.cfg.residual_log_every) != 0)
        {
            return;
        }

        std::ofstream out(residual_file_path(case_name), std::ios::app);

        if (!out)
        {
            throw std::runtime_error(
                "Post::writeResidualRow - cannot open residual CSV file");
        }

        out << std::setprecision(16)
            << iitime << ','
            << safe_value(s.cfg.rtime) << ','
            << safe_value(s.cfg.dtreal) << ','
            << safe_value(s.hb[0]) << ',' << safe_value(s.hb[1]) << ',' << safe_value(s.hb[2]) << ','
            << safe_value(s.hb[3]) << ',' << safe_value(s.hb[4]) << ',' << safe_value(s.hb[5]) << ','
            << safe_value(s.hb[6]) << ',' << safe_value(s.hb[7]) << ',' << safe_value(s.hb[8]) << ','
            << safe_value(s.hb[9]) << ',' << safe_value(s.hb[10]) << ',' << safe_value(s.hb[11]) << ','
            << safe_value(s.hb[12]) << ',' << safe_value(s.hb[13]) << ',' << safe_value(s.hb[14]) << ','
            << safe_value(Convergence::velocityResidual(s)) << ','
            << s.last_cg_iterations << ','
            << safe_value(s.last_cg_initial_l2) << ','
            << safe_value(s.last_cg_final_l2) << ','
            << safe_value(s.last_cg_relative_l2) << ','
            << safe_value(s.last_cg_max_abs)
            << '\n';
    }

    //=========================================================================
    // Tests whether a VTU solution is required.
    //
    // Output may be triggered by:
    //
    //     1. iteration interval;
    //     2. transient simulated-time interval.
    //
    // The initial state at iteration zero is always eligible when VTU output
    // is enabled.
    //=========================================================================
    bool Post::shouldWriteVTU(
        CBSStateSI& s,
        const Int iitime)
    {
        if (s.cfg.vtu_output_enabled < 1)
        {
            return false;
        }

        if (iitime == 0)
        {
            return true;
        }

        return write_by_iteration(s, iitime) ||
               write_by_simulated_time(s, iitime);
    }

    //=========================================================================
    // Writes one VTU solution and records it in the PVD time-series history.
    //=========================================================================
    void Post::writeSolution(
        CBSStateSI& s,
        const std::string& case_name,
        const Int iitime)
    {
        if (s.cfg.vtu_output_enabled < 1)
        {
            return;
        }

        writeVTU(s, case_name, iitime);

        s.output_iterations.push_back(iitime);
        s.output_times.push_back(s.cfg.rtime);

        writePVD(s, case_name);
    }

    //=========================================================================
    // Writes one VTK XML unstructured-grid file.
    //
    // Mesh cells:
    //
    //     tetrahedra          VTK cell type 10
    //     boundary triangles  VTK cell type 5
    //
    // Internal mesh numbering is one-based. VTK connectivity is written with
    // zero-based node indices.
    //
    // Point data:
    //
    //     pressure
    //     temperature
    //     velocity
    //     velocity magnitude
    //     material-domain masks
    //     optional boundary-condition debug identifier
    //
    // Cell data:
    //
    //     cell kind
    //     material identifier
    //     boundary-condition identifier
    //     parent tetrahedral element
    //=========================================================================
    void Post::writeVTU(
        const CBSStateSI& s,
        const std::string& case_name,
        const Int iitime)
    {
        ensure_output_directory(case_name);

        std::ofstream out(vtu_file_path(case_name, iitime), std::ios::trunc);

        if (!out)
        {
            throw std::runtime_error("Post::writeVTU - cannot open VTU file");
        }

        const Int total_cells = s.cfg.nelem + s.cfg.nboun;
        const std::vector<Int> node_bc =
            s.cfg.write_boundary_debug_arrays > 0
                ? build_node_bc_debug(s)
                : std::vector<Int>{};

        const NodePhysicsMasks masks = build_node_physics_masks(s);

        out << std::setprecision(16);
        out << "<?xml version=\"1.0\"?>\n";
        out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        out << "  <UnstructuredGrid>\n";
        out << "    <Piece NumberOfPoints=\"" << s.cfg.npoin
            << "\" NumberOfCells=\"" << total_cells << "\">\n";

        out << "      <Points>\n";
        out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            out << "          "
                << s.coord(1, ip) << ' '
                << s.coord(2, ip) << ' '
                << s.coord(3, ip) << '\n';
        }
        out << "        </DataArray>\n";
        out << "      </Points>\n";

        out << "      <Cells>\n";

        out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            out << "          "
                << s.intma(1, ie) - 1 << ' '
                << s.intma(2, ie) - 1 << ' '
                << s.intma(3, ie) - 1 << ' '
                << s.intma(4, ie) - 1 << '\n';
        }

        for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
        {
            out << "          "
                << s.iside(1, ib) - 1 << ' '
                << s.iside(2, ib) - 1 << ' '
                << s.iside(3, ib) - 1 << '\n';
        }
        out << "        </DataArray>\n";

        out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
        Int offset = 0;

        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            offset += 4;
            out << "          " << offset << '\n';
        }

        for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
        {
            offset += 3;
            out << "          " << offset << '\n';
        }
        out << "        </DataArray>\n";

        out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            out << "          10\n";
        }

        for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
        {
            out << "          5\n";
        }
        out << "        </DataArray>\n";

        out << "      </Cells>\n";

        out << "      <PointData Scalars=\"temperature\" Vectors=\"velocity\">\n";

        out << "        <DataArray type=\"Float64\" Name=\"pressure\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            const Real pressure =
                pressure_active_node(masks, ip)
                    ? safe_value(s.pres(ip))
                    : 0.0;

            out << "          " << pressure << '\n';
        }
        out << "        </DataArray>\n";

        out << "        <DataArray type=\"Float64\" Name=\"temperature\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            out << "          " << safe_value(s.temperature(ip)) << '\n';
        }
        out << "        </DataArray>\n";

        out << "        <DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            if (velocity_active_node(masks, ip))
            {
                out << "          "
                    << safe_value(s.unkno(1, ip)) << ' '
                    << safe_value(s.unkno(2, ip)) << ' '
                    << safe_value(s.unkno(3, ip)) << '\n';
            }
            else
            {
                out << "          0 0 0\n";
            }
        }
        out << "        </DataArray>\n";

        out << "        <DataArray type=\"Float64\" Name=\"velocity_magnitude\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            const Real velocity_magnitude =
                velocity_active_node(masks, ip)
                    ? safe_value(s.velocity(ip))
                    : 0.0;

            out << "          " << velocity_magnitude << '\n';
        }
        out << "        </DataArray>\n";

        out << "        <DataArray type=\"Int32\" Name=\"node_domain_kind\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            out << "          " << masks.node_domain_kind[static_cast<std::size_t>(ip)] << '\n';
        }
        out << "        </DataArray>\n";

        out << "        <DataArray type=\"Int32\" Name=\"pressure_active\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            out << "          " << static_cast<Int>(masks.pressure_active[static_cast<std::size_t>(ip)]) << '\n';
        }
        out << "        </DataArray>\n";

        out << "        <DataArray type=\"Int32\" Name=\"velocity_active\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            out << "          " << static_cast<Int>(masks.velocity_active[static_cast<std::size_t>(ip)]) << '\n';
        }
        out << "        </DataArray>\n";

        out << "        <DataArray type=\"Int32\" Name=\"thermal_active\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            out << "          " << static_cast<Int>(masks.thermal_active[static_cast<std::size_t>(ip)]) << '\n';
        }
        out << "        </DataArray>\n";

        if (s.cfg.write_boundary_debug_arrays > 0)
        {
            out << "        <DataArray type=\"Int32\" Name=\"node_bc_id_debug\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                out << "          " << node_bc[static_cast<std::size_t>(ip)] << '\n';
            }
            out << "        </DataArray>\n";
        }

        out << "      </PointData>\n";

        out << "      <CellData Scalars=\"mat_id\">\n";

        out << "        <DataArray type=\"Int32\" Name=\"cell_kind\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            out << "          1\n";
        }

        for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
        {
            out << "          2\n";
        }
        out << "        </DataArray>\n";

        out << "        <DataArray type=\"Int32\" Name=\"mat_id\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            out << "          " << s.mat_elem(ie) << '\n';
        }

        for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
        {
            out << "          -1\n";
        }
        out << "        </DataArray>\n";

        out << "        <DataArray type=\"Int32\" Name=\"bc_id\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            out << "          0\n";
        }

        for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
        {
            out << "          " << s.iside(s.cfg.bsid, ib) << '\n';
        }
        out << "        </DataArray>\n";

        out << "        <DataArray type=\"Int32\" Name=\"parent_element\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
        {
            out << "          " << ie << '\n';
        }

        for (Int ib = 1; ib <= s.cfg.nboun; ++ib)
        {
            out << "          " << s.iside(s.cfg.nsidpe, ib) << '\n';
        }
        out << "        </DataArray>\n";

        out << "      </CellData>\n";

        out << "    </Piece>\n";
        out << "  </UnstructuredGrid>\n";
        out << "</VTKFile>\n";
    }

    //=========================================================================
    // Rewrites the ParaView collection file from the stored output history.
    //
    // Each DataSet entry associates:
    //
    //     physical time <-> VTU file
    //=========================================================================
    void Post::writePVD(
        const CBSStateSI& s,
        const std::string& case_name)
    {
        ensure_output_directory(case_name);

        if (s.output_iterations.size() != s.output_times.size())
        {
            throw std::runtime_error(
                "Post::writePVD - output history arrays are inconsistent");
        }

        std::ofstream out(pvd_file_path(case_name), std::ios::trunc);

        if (!out)
        {
            throw std::runtime_error("Post::writePVD - cannot open PVD file");
        }

        out << std::setprecision(16);
        out << "<?xml version=\"1.0\"?>\n";
        out << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        out << "  <Collection>\n";

        for (std::size_t i = 0; i < s.output_iterations.size(); ++i)
        {
            out << "    <DataSet timestep=\""
                << s.output_times[i]
                << "\" group=\"\" part=\"0\" file=\""
                << vtu_file_name(case_name, s.output_iterations[i])
                << "\"/>\n";
        }

        out << "  </Collection>\n";
        out << "</VTKFile>\n";
    }

    //=========================================================================
    // Prints one in-place CBS progress line.
    //
    // Completion fraction:
    //
    //     f = iteration / requested_iterations
    //
    // Estimated remaining wall time:
    //
    //     ETA = elapsed * (1 - f) / f
    //
    // The line also reports physical time, dt, field residuals, Pressure CG
    // iterations and the Pressure CG relative residual.
    //=========================================================================
    void Post::printProgressLine(
        const CBSStateSI& s,
        const Int iitime,
        const bool force_newline)
    {
        if (s.cfg.console_log_every < 1 && !force_newline)
        {
            return;
        }

        if (!force_newline &&
            iitime != 1 &&
            (iitime % s.cfg.console_log_every) != 0)
        {
            return;
        }

        constexpr Int bar_width = 34;

        const Real fraction =
            std::clamp(
                static_cast<Real>(iitime) /
                static_cast<Real>(std::max<Int>(s.cfg.ntime, 1)),
                0.0,
                1.0);

        const Real elapsed = elapsed_wall_seconds();
        const Real eta =
            (fraction > 1.0e-12)
                ? elapsed * (1.0 - fraction) / fraction
                : 0.0;

        std::cout
            << '\r'
            << "  " << magenta << spinner(iitime) << reset << " "
            << cyan << "[" << green << progress_bar(fraction, bar_width) << cyan << "]" << reset
            << " " << bold << std::setw(6) << percent_string(fraction) << reset
            << "  " << dim << "it" << reset << "=" << iitime << "/" << s.cfg.ntime
            << "  " << dim << "t" << reset << "=" << sci(s.cfg.rtime, 2)
            << "  " << dim << "dt" << reset << "=" << sci(s.cfg.dtreal, 2)
            << "  " << yellow << "RelU" << reset << "=" << sci(s.hb[0], 2)
            << "  " << yellow << "RelV" << reset << "=" << sci(s.hb[3], 2)
            << "  " << yellow << "RelW" << reset << "=" << sci(s.hb[6], 2)
            << "  " << yellow << "RelP" << reset << "=" << sci(s.hb[9], 2);

        if (s.cfg.temp_calc > 0)
        {
            std::cout << "  " << yellow << "RelT" << reset << "=" << sci(s.hb[12], 2);
        }

        std::cout
            << "  " << blue << "CG" << reset << "=" << s.last_cg_iterations
            << "  " << blue << "CGrel" << reset << "=" << sci(s.last_cg_relative_l2, 2)
            << "  " << dim << "elapsed" << reset << "=" << format_seconds(elapsed)
            << "  " << dim << "eta" << reset << "=" << format_seconds(eta)
            << "     "
            << std::flush;

        if (force_newline)
        {
            std::cout << "\n";
        }
    }

    //=========================================================================
    // Prints the final stopping reason, completed CBS iteration, physical time,
    // wall-clock time, field residuals and last Pressure CG information.
    //=========================================================================
    void Post::printFinalSummary(
        const CBSStateSI& s,
        const Int last_iteration,
        const std::string& stop_reason)
    {
        std::cout << "\n";
        hr(green);
        centre_row("CBS3D++_SI RUN COMPLETE", bold);
        separator_row(green);

        key_value_row("Stop reason", stop_reason);
        key_value_row("Last iteration", std::to_string(last_iteration));
        key_value_row("Final physical time", sci(s.cfg.rtime, 6));
        key_value_row("Final dt", sci(s.cfg.dtreal, 6));
        key_value_row("Wall-clock elapsed", format_seconds(elapsed_wall_seconds()));
        key_value_row("Final RelU", sci(s.hb[0], 6));
        key_value_row("Final RelV", sci(s.hb[3], 6));
        key_value_row("Final RelW", sci(s.hb[6], 6));
        key_value_row("Final RelP", sci(s.hb[9], 6));

        if (s.cfg.temp_calc > 0)
        {
            key_value_row("Final RelT", sci(s.hb[12], 6));
        }

        key_value_row("Last CG iterations", std::to_string(s.last_cg_iterations));
        key_value_row("Last CG relative L2", sci(s.last_cg_relative_l2, 6));
        footer(green);
    }

    //=========================================================================
    // Launches tools/live_residual_plot.py when requested.
    //
    // The plotting process reads the residual CSV produced by the solver.
    // Failure to launch the external plotter is reported as a warning and does
    // not stop the numerical calculation.
    //=========================================================================
    void Post::launchLiveResidualPlotterIfRequested(
        const CBSStateSI& s,
        const std::string& case_name)
    {
        if (s.cfg.live_residual_plot < 1)
        {
            return;
        }

        if (s.cfg.residual_log_enabled < 1)
        {
            std::cout
                << "  " << red << "[!!]" << reset
                << " Live residual plotting requested, but residual logging is disabled.\n";
            return;
        }

        const fs::path script = fs::path("tools") / "live_residual_plot.py";
        const fs::path csv = residual_file_path(case_name);

        if (!fs::exists(script))
        {
            std::cout
                << "  " << yellow << "[!!]" << reset
                << " Live residual plotter requested, but "
                << script.string()
                << " was not found. Start it manually after copying the script.\n";
            return;
        }

#ifdef _WIN32
        const std::string command =
            "start \"CBS3D residuals\" python "
            + quote_path_for_command(script)
            + " "
            + quote_path_for_command(csv);
#else
        const std::string command =
            "python3 "
            + quote_path_for_command(script)
            + " "
            + quote_path_for_command(csv)
            + " >/dev/null 2>&1 &";
#endif

        const int rc = std::system(command.c_str());

        if (rc != 0)
        {
            std::cout
                << "  " << yellow << "[!!]" << reset
                << " Live residual plotter launch returned code "
                << rc
                << ". The solver will continue.\n";
        }
        else
        {
            std::cout
                << "  " << green << "[OK]" << reset
                << " Live residual monitor launched.\n";
        }
    }
}
