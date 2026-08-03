#!/usr/bin/env python3
"""
Patch CBS3D++_SI so nodal SA eddy viscosity is forced to zero on
SA wall nodes after element-to-node averaging.

Why this patch exists
---------------------
TurbulenceBoundary::applyWallValues correctly sets

    nu_tilde = 0, nu_t = 0, mu_t = 0

on SA wall nodes.  However, SpalartAllmarasAssembly::updateEddyViscosity()
then recomputes nodal nu_t/mu_t by averaging neighbouring element eddy
viscosities.  A wall node belongs to near-wall fluid elements whose element
averaged nu_tilde can be nonzero because of the adjacent interior nodes.  The
averaging step can therefore reintroduce nonzero nodal nu_t/mu_t on a wall.

This patch keeps element effective viscosities unchanged, but makes the nodal
VTU/diagnostic fields obey the strong SA wall condition exactly:

    wall/inactive node => nu_t = 0 and mu_t = 0
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT / "src" / "assembly" / "SpalartAllmarasAssembly.cpp"

OLD = '''        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            const Int count = nodal_count[static_cast<std::size_t>(ip)];

            if (count > 0)
            {
                s.nu_t(ip) = nodal_nu_t_sum[static_cast<std::size_t>(ip)]
                    / static_cast<Real>(count);

                s.mu_t(ip) = nodal_mu_t_sum[static_cast<std::size_t>(ip)]
                    / static_cast<Real>(count);
            }
        }
'''

NEW = '''        for (Int ip = 1; ip <= s.cfg.npoin; ++ip)
        {
            // Wall and inactive nodes must not inherit eddy viscosity from
            // neighbouring element averages.  The transported SA wall value is
            // imposed strongly as nu_tilde = 0, and the derived nodal eddy
            // viscosity must remain exactly zero for VTU diagnostics and for
            // any later nodal use of the turbulence field.
            if (s.sa_active_node(ip) == 0 || s.sa_wall_node(ip) != 0)
            {
                s.nu_t(ip) = 0.0;
                s.mu_t(ip) = 0.0;
                continue;
            }

            const Int count = nodal_count[static_cast<std::size_t>(ip)];

            if (count > 0)
            {
                s.nu_t(ip) = nodal_nu_t_sum[static_cast<std::size_t>(ip)]
                    / static_cast<Real>(count);

                s.mu_t(ip) = nodal_mu_t_sum[static_cast<std::size_t>(ip)]
                    / static_cast<Real>(count);
            }
            else
            {
                s.nu_t(ip) = 0.0;
                s.mu_t(ip) = 0.0;
            }
        }
'''


def main() -> int:
    text = TARGET.read_text(encoding="utf-8")

    if NEW in text:
        print("[skip] SA wall eddy-viscosity cleanup is already present.")
        return 0

    if OLD not in text:
        raise SystemExit(
            "Patch anchor not found. SpalartAllmarasAssembly.cpp has changed; "
            "inspect updateEddyViscosity() manually."
        )

    text = text.replace(OLD, NEW, 1)
    TARGET.write_text(text, encoding="utf-8", newline="\n")

    print("[patch] SpalartAllmarasAssembly.cpp: enforce zero nodal nu_t/mu_t on SA wall nodes")
    print("SA wall eddy-viscosity cleanup patch applied.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
