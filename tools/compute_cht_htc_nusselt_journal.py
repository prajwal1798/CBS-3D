#!/usr/bin/env python3
"""Restyle CBS3D CHT SVG figures for journal-ready presentation.

This module changes only SVG presentation. It does not read or modify the
numerical CSV/VTP data produced by compute_cht_htc_nusselt.py. The
implementation uses only the Python standard library for Sunbird.
"""
from __future__ import print_function

import re
from pathlib import Path


FILES = {
    "cht_wall_and_bulk_temperature_vs_s.svg": {
        "title": "Wall and bulk-fluid temperature",
        "x": "Streamwise coordinate, s/Dh",
        "y": "Temperature [K]",
        "labels": {
            "Wall": "Wall temperature, T_w",
            "Bulk coolant": "Bulk-fluid temperature, T_b",
        },
    },
    "cht_wall_superheat_vs_s.svg": {
        "title": "Wall-to-bulk temperature difference",
        "x": "Streamwise coordinate, s/Dh",
        "y": "T_w - T_b [K]",
        "labels": {"Wall superheat": "T_w - T_b"},
    },
    "cht_heat_flux_vs_s.svg": {
        "title": "Interfacial heat-flux reconstruction",
        "x": "Streamwise coordinate, s/Dh",
        "y": "Heat flux, q'' [W m^-2]",
        "labels": {
            "Fluid-side P1": "Fluid-side P1 flux",
            "Solid-side P1": "Solid-side P1 flux",
            "Energy-constrained recovered":
                "Energy-constrained recovered flux",
        },
    },
    "cht_htc_vs_s.svg": {
        "title": "Local heat-transfer coefficient",
        "x": "Streamwise coordinate, s/Dh",
        "y": "h [W m^-2 K^-1]",
        "labels": {"Sectional HTC": "Section-averaged h"},
    },
    "cht_nusselt_vs_s.svg": {
        "title": "Local Nusselt number",
        "x": "Streamwise coordinate, s/Dh",
        "y": "Nu",
        "labels": {"Sectional Nu": "Section-averaged Nu"},
    },
    "cht_mass_flow_vs_s.svg": {
        "title": "Cross-sectional mass-flow conservation",
        "x": "Streamwise coordinate, s/Dh",
        "y": "Mass-flow rate [kg s^-1]",
        "labels": {
            "Net": "Net",
            "Forward": "Forward",
            "Reverse": "Reverse",
        },
    },
    "cht_flux_jump_vs_s.svg": {
        "title": "One-sided interfacial flux discontinuity",
        "x": "Streamwise coordinate, s/Dh",
        "y": "Relative flux jump",
        "labels": {"Area-weighted jump": "Area-weighted jump"},
    },
}


FONT_STACK = (
    "STIX Two Text, STIXGeneral, Times New Roman, "
    "DejaVu Serif, serif"
)


def _replace_text_element(svg, old_text, new_text):
    pattern = re.compile(
        r"(<text\b[^>]*>)" + re.escape(old_text) + r"(</text>)"
    )
    return pattern.sub(
        lambda match: match.group(1) + new_text + match.group(2),
        svg,
    )


def _restyle_one(path, spec):
    svg = path.read_text(encoding="utf-8")

    svg = svg.replace(
        '<svg xmlns="http://www.w3.org/2000/svg"',
        '<svg xmlns="http://www.w3.org/2000/svg" '
        'shape-rendering="geometricPrecision" '
        'text-rendering="geometricPrecision"',
        1,
    )
    svg = svg.replace(
        'font-family="Arial"',
        'font-family="{}"'.format(FONT_STACK),
    )
    svg = svg.replace(
        'font-size="24" font-weight="bold"',
        'font-size="23" font-weight="600"',
    )
    svg = svg.replace('font-size="18"', 'font-size="19"')
    svg = svg.replace('font-size="14"', 'font-size="15"')

    svg = svg.replace(
        'stroke="#dddddd" stroke-width="1"',
        'stroke="#e6e6e6" stroke-width="0.8"',
    )
    svg = svg.replace(
        'stroke="black" stroke-width="1.5"',
        'stroke="#222222" stroke-width="1.2"',
    )
    svg = svg.replace('stroke-width="2.5"', 'stroke-width="2.2"')
    svg = svg.replace('stroke-width="3"', 'stroke-width="2.4"')

    svg = svg.replace('stroke="#1f77b4"', 'stroke="#0072B2"')
    svg = svg.replace(
        'stroke="#d62728"',
        'stroke="#D55E00" stroke-dasharray="9 5"',
    )
    svg = svg.replace(
        'stroke="#2ca02c"',
        'stroke="#009E73" stroke-dasharray="2.5 3.5"',
    )
    svg = svg.replace(
        'stroke="#9467bd"',
        'stroke="#CC79A7" stroke-dasharray="12 4 2 4"',
    )

    svg = re.sub(r"^<circle\b[^>]*/>\n?", "", svg, flags=re.MULTILINE)

    title_match = re.search(
        r'(<text x="600\.0" y="38"[^>]*>)(.*?)(</text>)',
        svg,
    )
    if title_match:
        svg = (
            svg[:title_match.start(2)]
            + spec["title"]
            + svg[title_match.end(2):]
        )

    for old_text, new_text in spec["labels"].items():
        svg = _replace_text_element(svg, old_text, new_text)

    svg = re.sub(
        r'(<text x="635\.0" y="690"[^>]*>).*?(</text>)',
        lambda match: match.group(1) + spec["x"] + match.group(2),
        svg,
    )
    svg = re.sub(
        r'(<text x="28" y="348\.5"[^>]*'
        r'transform="rotate\(-90 28 348\.5\)"[^>]*>).*?(</text>)',
        lambda match: match.group(1) + spec["y"] + match.group(2),
        svg,
    )

    svg = svg.replace(
        '<polyline ',
        '<polyline stroke-linecap="round" '
        'stroke-linejoin="round" ',
    )

    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(svg, encoding="utf-8")
    temporary.replace(path)


def regenerate_journal_plots(output_dir):
    output_dir = Path(output_dir).resolve()
    if not output_dir.is_dir():
        raise RuntimeError(
            "Plot directory does not exist: {}".format(output_dir)
        )

    missing = [
        name for name in FILES
        if not (output_dir / name).is_file()
    ]
    if missing:
        raise RuntimeError(
            "Missing SVG files in {}: {}".format(
                output_dir,
                ", ".join(missing),
            )
        )

    for name, specification in FILES.items():
        _restyle_one(output_dir / name, specification)
        print("restyled: {}".format(name))

    print("PASS: journal-style CHT figures regenerated")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=Path)
    arguments = parser.parse_args()
    regenerate_journal_plots(arguments.output_dir)
