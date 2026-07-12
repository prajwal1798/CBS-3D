from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return path.read_text(encoding="utf-8")


def write(path, text):
    path.write_text(text, encoding="utf-8")


def replace_once(text, old, new, label):
    if new in text:
        print(f"[skip] {label}")
        return text

    if old not in text:
        raise SystemExit(f"Patch anchor not found: {label}")

    print(f"[patch] {label}")
    return text.replace(old, new, 1)


def insert_after(text, anchor, insertion, label):
    if insertion in text:
        print(f"[skip] {label}")
        return text

    if anchor not in text:
        raise SystemExit(f"Patch anchor not found: {label}")

    print(f"[patch] {label}")
    return text.replace(anchor, anchor + insertion, 1)


def patch_meshio():
    path = ROOT / "src" / "io" / "MeshIO.cpp"
    text = read(path)

    anchor_function = """        void read_artificial_diffusion_control(std::istream& in, CBSStateSI& s)
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
"""

    sa_function = """

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
        //     sa_use_stilde_limiter sa_implicit_destruction
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
            }
        }
"""

    text = insert_after(
        text,
        anchor_function,
        sa_function,
        "MeshIO.cpp: add optional SA .par reader")

    old_call = """        read_pressure_cg_controls(in, s);
        read_output_monitor_controls(in, s);
        read_artificial_diffusion_control(in, s);
"""

    new_call = """        read_pressure_cg_controls(in, s);
        read_output_monitor_controls(in, s);
        read_artificial_diffusion_control(in, s);
        read_optional_spalart_allmaras_controls(in, s);
"""

    text = replace_once(
        text,
        old_call,
        new_call,
        "MeshIO.cpp: read optional SA controls")

    old_validation_anchor = """        s.cfg.istart = 1;
        s.cfg.rtime = 0.0;
"""

    sa_validation = """
        if (s.cfg.turbulence_on < 0 || s.cfg.turbulence_on > 1)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - TURBULENCE_ON must be 0 or 1");
        }

        if (s.cfg.turbulence_model < 0 || s.cfg.turbulence_model > 1)
        {
            throw std::runtime_error(
                "MeshIO::readParameterFile - TURBULENCE_MODEL must be 0 or 1");
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

"""

    text = replace_once(
        text,
        old_validation_anchor,
        sa_validation + old_validation_anchor,
        "MeshIO.cpp: validate SA controls")

    write(path, text)


def patch_boundary():
    path = ROOT / "src" / "boundary" / "Boundary.cpp"
    text = read(path)

    old_normal_guard = """            if (bc == s.cfg.bc_massflow_temperature_inlet)
            {
                // Mass-flow inlet in 3D uses the inward normal direction.
                // face_norm is outward from the domain, so inflow is -n.
                unit_outward_normal(s, ib, nx, ny, nz);
            }
"""

    new_normal_guard = """            if (bc == s.cfg.bc_massflow_temperature_inlet &&
                s.cfg.mass_flow_inlet_enabled > 0)
            {
                // Mass-flow inlet in 3D uses the inward normal direction.
                // face_norm is outward from the domain, so inflow is -n.
                unit_outward_normal(s, ib, nx, ny, nz);
            }
"""

    text = replace_once(
        text,
        old_normal_guard,
        new_normal_guard,
        "Boundary.cpp: compute inlet normal only for mass-flow inlet")

    old_inlet = """                else if (bc == s.cfg.bc_massflow_temperature_inlet)
                {
                    set_velocity(
                        s,
                        ip,
                        -s.cfg.inlet_u_from_massflow * nx,
                        -s.cfg.inlet_u_from_massflow * ny,
                        -s.cfg.inlet_u_from_massflow * nz);
                }
"""

    new_inlet = """                else if (bc == s.cfg.bc_massflow_temperature_inlet)
                {
                    if (s.cfg.mass_flow_inlet_enabled > 0)
                    {
                        set_velocity(
                            s,
                            ip,
                            -s.cfg.inlet_u_from_massflow * nx,
                            -s.cfg.inlet_u_from_massflow * ny,
                            -s.cfg.inlet_u_from_massflow * nz);
                    }
                    else
                    {
                        // Flow-only benchmarks, such as the turbulent flat plate,
                        // often use BC 511 as a uniform velocity inlet without a
                        // prescribed mass-flow rate.  In that case the velocity
                        // components are read directly from the .par file.
                        set_velocity(
                            s,
                            ip,
                            s.cfg.inlet_u,
                            s.cfg.inlet_v,
                            s.cfg.inlet_w);
                    }
                }
"""

    text = replace_once(
        text,
        old_inlet,
        new_inlet,
        "Boundary.cpp: allow prescribed velocity on BC 511")

    write(path, text)


def patch_post():
    path = ROOT / "src" / "io" / "Post.cpp"
    text = read(path)

    old_summary = """        key_value_row("Energy equation", on_off(s.cfg.temp_calc));
        key_value_row(
            "Material mode",
"""

    new_summary = """        key_value_row("Energy equation", on_off(s.cfg.temp_calc));
        key_value_row(
            "Turbulence model",
            s.cfg.turbulence_on > 0
                ? "Spalart-Allmaras"
                : "OFF");
        key_value_row(
            "Material mode",
"""

    text = replace_once(
        text,
        old_summary,
        new_summary,
        "Post.cpp: print turbulence model in run summary")

    marker = """        out << "        <DataArray type=\"Float64\" Name=\"velocity_magnitude\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            const Real velocity_magnitude =
                velocity_active_node(masks, ip)
                    ? safe_value(s.velocity(ip))
                    : 0.0;

            out << "          " << velocity_magnitude << '\n';
        }
        out << "        </DataArray>\n";
"""

    insertion = """

        if (s.cfg.turbulence_on > 0)
        {
            out << "        <DataArray type=\"Float64\" Name=\"nu_tilde\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                out << "          " << safe_value(s.nu_tilde(ip)) << '\n';
            }
            out << "        </DataArray>\n";

            out << "        <DataArray type=\"Float64\" Name=\"nu_t\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                out << "          " << safe_value(s.nu_t(ip)) << '\n';
            }
            out << "        </DataArray>\n";

            out << "        <DataArray type=\"Float64\" Name=\"mu_t\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                out << "          " << safe_value(s.mu_t(ip)) << '\n';
            }
            out << "        </DataArray>\n";

            out << "        <DataArray type=\"Float64\" Name=\"wall_distance\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                out << "          " << safe_value(s.wall_distance(ip)) << '\n';
            }
            out << "        </DataArray>\n";

            out << "        <DataArray type=\"Float64\" Name=\"sa_residual\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                out << "          " << safe_value(s.sa_residual(ip)) << '\n';
            }
            out << "        </DataArray>\n";

            out << "        <DataArray type=\"Float64\" Name=\"sa_production\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                out << "          " << safe_value(s.sa_production(ip)) << '\n';
            }
            out << "        </DataArray>\n";

            out << "        <DataArray type=\"Float64\" Name=\"sa_destruction\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                out << "          " << safe_value(s.sa_destruction(ip)) << '\n';
            }
            out << "        </DataArray>\n";

            out << "        <DataArray type=\"Float64\" Name=\"sa_diffusion\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                out << "          " << safe_value(s.sa_diffusion(ip)) << '\n';
            }
            out << "        </DataArray>\n";

            out << "        <DataArray type=\"Float64\" Name=\"sa_source\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                out << "          " << safe_value(s.sa_source(ip)) << '\n';
            }
            out << "        </DataArray>\n";

            out << "        <DataArray type=\"Int32\" Name=\"sa_active_node\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                out << "          " << s.sa_active_node(ip) << '\n';
            }
            out << "        </DataArray>\n";

            out << "        <DataArray type=\"Int32\" Name=\"sa_wall_node\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                out << "          " << s.sa_wall_node(ip) << '\n';
            }
            out << "        </DataArray>\n";

            out << "        <DataArray type=\"Int32\" Name=\"sa_inlet_node\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                out << "          " << s.sa_inlet_node(ip) << '\n';
            }
            out << "        </DataArray>\n";
        }
"""

    if insertion in text:
        print("[skip] Post.cpp: write SA point diagnostics")
    elif marker in text:
        print("[patch] Post.cpp: write SA point diagnostics")
        text = text.replace(marker, marker + insertion, 1)
    else:
        raise SystemExit("Patch anchor not found: Post.cpp: write SA point diagnostics")

    write(path, text)


def main():
    patch_meshio()
    patch_boundary()
    patch_post()
    print("SA input/output support patch applied.")


if __name__ == "__main__":
    main()
