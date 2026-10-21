#!/usr/bin/env python3
"""Summarise and integrity-check a native CBS3D mesh case.

The native case is the set of files sharing one base name:

    <base>.plt       tetrahedra, node coordinates and boundary faces
    <base>.bco       declared boundary-condition codes
    <base>.material  element material flags (optional)
    <base>.mpi       rank-local partition metadata (optional)

This tool does not require Gmsh.  It is intended both for the assembled
single-domain case and for each rank-local partition written by the Gmsh
workflow.  All indices in the native files are one-based.

Exit status:
    0  every check passed
    1  at least one check failed
    2  the case could not be read

Examples:
    python tools/mesh_summary.py cases/flat_plate_local/flat_plate_local
    python tools/mesh_summary.py cases/flat_plate_p8/ranks/rank_0000/... \
        --json report.json
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Sequence, Tuple


NODE_RECORD = re.compile(
    r"^\s*(\d+)\s*"
    r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+"
    r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+"
    r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*$"
)


def read_nonempty_lines(path: Path) -> List[str]:
    lines: List[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for raw in stream:
            text = raw.strip()
            if text and text[0] not in "#!":
                lines.append(text)
    return lines


def parse_ints(line: str, count: int, context: str) -> List[int]:
    fields = line.split()
    if len(fields) < count:
        raise ValueError(f"{context}: expected {count} integers, got {len(fields)}")
    return [int(value) for value in fields[:count]]


def tetra_volume(
    a: Sequence[float],
    b: Sequence[float],
    c: Sequence[float],
    d: Sequence[float],
) -> float:
    ab = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
    ac = (c[0] - a[0], c[1] - a[1], c[2] - a[2])
    ad = (d[0] - a[0], d[1] - a[1], d[2] - a[2])
    cross = (
        ac[1] * ad[2] - ac[2] * ad[1],
        ac[2] * ad[0] - ac[0] * ad[2],
        ac[0] * ad[1] - ac[1] * ad[0],
    )
    return (ab[0] * cross[0] + ab[1] * cross[1] + ab[2] * cross[2]) / 6.0


def sorted_face(nodes: Sequence[int]) -> Tuple[int, int, int]:
    return tuple(sorted(int(value) for value in nodes))  # type: ignore[return-value]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base", type=Path, help="Case base path without extension")
    parser.add_argument("--json", type=Path, help="Write the full report as JSON")
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Only write the JSON report and print PASS/FAIL",
    )
    args = parser.parse_args()

    base = args.base.expanduser().resolve()
    plt_path = base.with_suffix(".plt")
    bco_path = base.with_suffix(".bco")
    material_path = base.with_suffix(".material")
    mpi_path = base.with_suffix(".mpi")

    for required in (plt_path, bco_path):
        if not required.is_file():
            print(f"ERROR: missing required file: {required}", file=sys.stderr)
            return 2

    try:
        plt_lines = read_nonempty_lines(plt_path)
        header = parse_ints(plt_lines[0], 3, f"{plt_path.name} header")
        nelem, npoin, nboun = header
        if min(nelem, npoin, nboun) <= 0:
            raise ValueError("non-positive mesh size in .plt header")

        expected = 1 + nelem + npoin + nboun
        if len(plt_lines) != expected:
            raise ValueError(
                f"{plt_path.name}: header implies {expected} records, "
                f"found {len(plt_lines)}"
            )

        cursor = 1
        elements: List[Tuple[int, int, int, int, int]] = []
        for _ in range(nelem):
            eid, n0, n1, n2, n3 = parse_ints(
                plt_lines[cursor], 5, "element record"
            )
            elements.append((eid, n0, n1, n2, n3))
            cursor += 1

        coords: Dict[int, Tuple[float, float, float]] = {}
        for _ in range(npoin):
            raw_record = plt_lines[cursor]
            match = NODE_RECORD.match(raw_record)
            if match is None:
                raise ValueError(
                    f"node record: cannot parse id and coordinates: {raw_record!r}"
                )
            node_id = int(match.group(1))
            coords[node_id] = (
                float(match.group(2)),
                float(match.group(3)),
                float(match.group(4)),
            )
            cursor += 1

        boundaries: List[Tuple[int, int, int, int, int]] = []
        for _ in range(nboun):
            n0, n1, n2, parent, bc = parse_ints(
                plt_lines[cursor], 5, "boundary record"
            )
            boundaries.append((n0, n1, n2, parent, bc))
            cursor += 1

    except (ValueError, IndexError) as error:
        print(f"ERROR: cannot parse {plt_path}: {error}", file=sys.stderr)
        return 2

    checks: Dict[str, bool] = {}
    issues: List[str] = []

    node_ids = set(coords)
    checks["node_ids_are_one_based"] = min(node_ids) == 1
    checks["node_ids_are_contiguous"] = len(node_ids) == npoin and max(node_ids) == npoin

    bad_connectivity = 0
    volumes: List[float] = []
    face_to_element: Dict[Tuple[int, int, int], int] = {}
    face_counts: Counter[Tuple[int, int, int]] = Counter()

    for eid, n0, n1, n2, n3 in elements:
        nodes = (n0, n1, n2, n3)
        if any(node not in coords for node in nodes):
            bad_connectivity += 1
            continue
        a, b, c, d = (coords[node] for node in nodes)
        volumes.append(abs(tetra_volume(a, b, c, d)))
        # Record the first owner of each face key.  Interior faces are shared
        # by exactly two elements; a count above two is non-manifold.
        for face in (
            (n1, n2, n3),
            (n0, n3, n2),
            (n0, n1, n3),
            (n0, n2, n1),
        ):
            key = sorted_face(face)
            face_counts[key] += 1
            face_to_element.setdefault(key, eid)

    checks["element_connectivity_resolves"] = bad_connectivity == 0
    checks["all_tetrahedra_have_positive_volume"] = (
        len(volumes) == nelem and all(value > 0.0 for value in volumes)
    )
    checks["no_nonmanifold_element_faces"] = (
        max(face_counts.values(), default=0) <= 2
    )

    xyz = list(coords.values())
    bbox = {
        "xmin": min(v[0] for v in xyz),
        "xmax": max(v[0] for v in xyz),
        "ymin": min(v[1] for v in xyz),
        "ymax": max(v[1] for v in xyz),
        "zmin": min(v[2] for v in xyz),
        "zmax": max(v[2] for v in xyz),
    }

    bc_counts = Counter(record[4] for record in boundaries)
    bad_boundary_faces = 0
    bad_parents = 0
    for n0, n1, n2, parent, _bc in boundaries:
        if n0 not in coords or n1 not in coords or n2 not in coords:
            bad_boundary_faces += 1
            continue
        key = sorted_face((n0, n1, n2))
        if key not in face_to_element:
            bad_boundary_faces += 1
        elif parent != face_to_element[key]:
            bad_parents += 1

    checks["boundary_faces_match_element_faces"] = bad_boundary_faces == 0
    checks["boundary_parent_is_owner"] = bad_parents == 0

    try:
        bco_lines = read_nonempty_lines(bco_path)
        nflag, _mass_check = parse_ints(bco_lines[0], 2, f"{bco_path.name} header")
        declared = {
            parse_ints(line, 2, "bco mapping")[0]
            for line in bco_lines[1 : 1 + nflag]
        }
        checks["bco_mapping_count"] = nflag == len(bco_lines) - 1
        checks["bco_matches_plt_bcs"] = declared == set(bc_counts)
    except (ValueError, IndexError) as error:
        declared = set()
        checks["bco_mapping_count"] = False
        checks["bco_matches_plt_bcs"] = False
        issues.append(f"cannot parse {bco_path}: {error}")

    material_counts: Dict[int, int] = {}
    if material_path.is_file():
        material_lines = read_nonempty_lines(material_path)
        checks["material_record_count"] = len(material_lines) == nelem
        materials: Counter[int] = Counter()
        material_ids = []
        for line in material_lines:
            fields = line.split()
            if len(fields) < 6:
                checks["material_record_count"] = False
                break
            material_ids.append(int(fields[0]))
            materials[int(fields[5])] += 1
        checks["material_ids_are_one_based"] = material_ids == list(
            range(1, nelem + 1)
        )
        material_counts = dict(sorted(materials.items()))

    mpi_report: Dict[str, object] = {}
    if mpi_path.is_file():
        try:
            mpi_lines = read_nonempty_lines(mpi_path)
            partition = parse_ints(mpi_lines[1].replace("PARTITION", ""), 3, "PARTITION")
            global_counts = parse_ints(mpi_lines[2].replace("GLOBAL", ""), 3, "GLOBAL")
            local_counts = parse_ints(mpi_lines[3].replace("LOCAL", ""), 5, "LOCAL")
            mpi_report = {
                "partition_id": partition[0],
                "mpi_rank": partition[1],
                "parts": partition[2],
                "global_tetrahedra": global_counts[0],
                "global_nodes": global_counts[1],
                "global_boundary_faces": global_counts[2],
                "local_tetrahedra": local_counts[0],
                "local_nodes": local_counts[1],
                "owned_nodes": local_counts[2],
                "ghost_nodes": local_counts[3],
                "local_boundary_faces": local_counts[4],
            }
            checks["mpi_local_counts_match_plt"] = (
                local_counts[0] == nelem
                and local_counts[1] == npoin
                and local_counts[4] == nboun
            )
        except (ValueError, IndexError) as error:
            checks["mpi_local_counts_match_plt"] = False
            issues.append(f"cannot parse {mpi_path}: {error}")

    volume_total = sum(volumes)
    report = {
        "base": str(base),
        "result": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "issues": issues,
        "mesh": {
            "tetrahedra": nelem,
            "nodes": npoin,
            "boundary_faces": nboun,
            "total_volume": volume_total,
            "bbox": bbox,
        },
        "boundary_conditions": {str(key): value for key, value in sorted(bc_counts.items())},
        "materials": {str(key): value for key, value in material_counts.items()},
        "partition": mpi_report,
    }

    if args.json:
        args.json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    if not args.quiet:
        print(f"Case          : {base}")
        print(f"Tetrahedra    : {nelem}")
        print(f"Nodes         : {npoin}")
        print(f"Boundary faces: {nboun}")
        print(f"Total volume  : {volume_total:.12e}")
        print(
            "Bounding box  : "
            f"x[{bbox['xmin']:.6g}, {bbox['xmax']:.6g}] "
            f"y[{bbox['ymin']:.6g}, {bbox['ymax']:.6g}] "
            f"z[{bbox['zmin']:.6g}, {bbox['zmax']:.6g}]"
        )
        print("BC face counts: " + ", ".join(
            f"{key}={value}" for key, value in sorted(bc_counts.items())
        ))
        if material_counts:
            print("Material tets : " + ", ".join(
                f"{key}={value}" for key, value in material_counts.items()
            ))
        if mpi_report:
            print(
                "Partition     : "
                f"id={mpi_report['partition_id']} rank={mpi_report['mpi_rank']} "
                f"parts={mpi_report['parts']} owned={mpi_report['owned_nodes']} "
                f"ghosts={mpi_report['ghost_nodes']}"
            )
        failed = [name for name, ok in checks.items() if not ok]
        if failed:
            print("Failed checks : " + ", ".join(failed))

    print(f"RESULT: {report['result']}")
    return 0 if report["result"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
