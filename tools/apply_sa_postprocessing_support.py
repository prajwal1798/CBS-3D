from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return path.read_text(encoding="utf-8")


def write(path, text):
    path.write_text(text, encoding="utf-8")


def insert_after(text, anchor, insertion, label):
    if insertion in text:
        print(f"[skip] {label}")
        return text

    if anchor not in text:
        raise SystemExit(f"Patch anchor not found: {label}")

    print(f"[patch] {label}")
    return text.replace(anchor, anchor + insertion, 1)


def patch_post_cpp():
    path = ROOT / "src" / "io" / "Post.cpp"
    text = read(path)

    helper_anchor = """        // Tests whether velocity is physically active at one node.
        bool velocity_active_node(
            const NodePhysicsMasks& masks,
            const Int ip)
        {
            return masks.velocity_active[static_cast<std::size_t>(ip)] != 0;
        }
"""

    helper_insertion = """

        // Builds a nodal molecular-viscosity field from fluid-element material
        // values.  The array is used only for post-processing diagnostics such
        // as mu_t/mu.  It does not modify the governing-equation coefficients.
        //
        // Only fluid elements are included:
        //
        //     mat_elem(e) = 0
        //
        // This keeps CHT solids out of the turbulence ratio and makes the same
        // output path valid for both all-fluid SA cases and later CHT+turbulence
        // cases.  Nodes not touched by any fluid element receive zero and their
        // mu_t/mu value is written as zero.
        std::vector<Real> build_nodal_fluid_molecular_viscosity(const CBSStateSI& s)
        {
            const std::size_t n = static_cast<std::size_t>(s.cfg.npoin + 1);

            std::vector<Real> mu_sum(n, 0.0);
            std::vector<Int> count(n, 0);

            for (Int ie = 1; ie <= s.cfg.nelem; ++ie)
            {
                if (s.mat_elem(ie) != 0)
                {
                    continue;
                }

                const Real mu = s.mu_e(ie);

                for (Int in = 1; in <= s.cfg.nep; ++in)
                {
                    const Int ip = s.intma(in, ie);

                    if (ip < 1 || ip > s.cfg.npoin)
                    {
                        throw std::runtime_error(
                            "Post::writeVTU - element connectivity node out of range while building nodal viscosity");
                    }

                    mu_sum[static_cast<std::size_t>(ip)] += mu;
                    ++count[static_cast<std::size_t>(ip)];
                }
            }

            std::vector<Real> mu_node(n, 0.0);

            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                const Int c = count[static_cast<std::size_t>(ip)];

                if (c > 0)
                {
                    mu_node[static_cast<std::size_t>(ip)] =
                        mu_sum[static_cast<std::size_t>(ip)] / static_cast<Real>(c);
                }
            }

            return mu_node;
        }
"""

    text = insert_after(
        text,
        helper_anchor,
        helper_insertion,
        "Post.cpp: add nodal molecular-viscosity diagnostic helper")

    output_anchor = """        out << "        <DataArray type=\"Float64\" Name=\"velocity_magnitude\" NumberOfComponents=\"1\" format=\"ascii\">\n";
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

    output_insertion = """

        if (s.cfg.turbulence_on > 0)
        {
            const std::vector<Real> nodal_mu = build_nodal_fluid_molecular_viscosity(s);

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

            out << "        <DataArray type=\"Float64\" Name=\"mu_t_over_mu\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                const Real mu = nodal_mu[static_cast<std::size_t>(ip)];
                const Real ratio =
                    mu > 0.0
                        ? s.mu_t(ip) / mu
                        : 0.0;

                out << "          " << safe_value(ratio) << '\n';
            }
            out << "        </DataArray>\n";

            out << "        <DataArray type=\"Float64\" Name=\"wall_distance\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                out << "          " << safe_value(s.wall_distance(ip)) << '\n';
            }
            out << "        </DataArray>\n";

            out << "        <DataArray type=\"Float64\" Name=\"sa_rhs\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                out << "          " << safe_value(s.sa_rhs(ip)) << '\n';
            }
            out << "        </DataArray>\n";

            out << "        <DataArray type=\"Float64\" Name=\"sa_residual\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                out << "          " << safe_value(s.sa_residual(ip)) << '\n';
            }
            out << "        </DataArray>\n";

            out << "        <DataArray type=\"Float64\" Name=\"sa_source\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
            {
                out << "          " << safe_value(s.sa_source(ip)) << '\n';
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

    text = insert_after(
        text,
        output_anchor,
        output_insertion,
        "Post.cpp: write SA point diagnostics")

    write(path, text)


def main():
    patch_post_cpp()
    print("SA post-processing support patch applied.")


if __name__ == "__main__":
    main()
