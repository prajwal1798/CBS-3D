#!/usr/bin/env python3
"""Audit a distributed CBS3D VTU/PVTU result before computing h and Nu.

The audit is intentionally dependency-free. It reads the ASCII VTK XML pieces
sequentially and verifies that the topology and identifiers required for an
accurate conjugate heat-transfer post-processor are present.

This script is compatible with the older Python interpreter currently provided
on Sunbird.
"""

import argparse
import collections
import hashlib
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Optional, Tuple


REQUIRED_POINT_ARRAYS = {
    "temperature",
    "velocity",
    "global_node_id",
    "vtkGhostType",
    "node_domain_kind",
}

REQUIRED_CELL_ARRAYS = {
    "cell_kind",
    "global_element_id",
    "material_id",
    "bc_id",
    "parent_global_element",
}


def parse_numbers(text, cast):
    """Convert one ASCII VTK array payload into a Python list."""
    if not text:
        return []
    return [cast(token) for token in text.split()]


def named_arrays(parent):
    """Return DataArray children indexed by their Name attribute."""
    if parent is None:
        return {}

    result = {}  # type: Dict[str, ET.Element]
    for array in parent.findall("DataArray"):
        name = array.attrib.get("Name")
        if name:
            result[name] = array
    return result


def find_piece_sources(pvtu):
    """Resolve all rank-local VTU files referenced by a PVTU descriptor."""
    root = ET.parse(str(pvtu)).getroot()
    grid = root.find("PUnstructuredGrid")
    if grid is None:
        raise RuntimeError("{} is not a PUnstructuredGrid file".format(pvtu))

    pieces = []  # type: List[Path]
    for piece in grid.findall("Piece"):
        source = piece.attrib.get("Source")
        if not source:
            raise RuntimeError("PVTU piece without Source in {}".format(pvtu))
        pieces.append((pvtu.parent / source).resolve())

    if not pieces:
        raise RuntimeError("No Piece entries found in {}".format(pvtu))

    return pieces


def audit_piece(path):
    """Audit one rank-local VTU piece and return aggregateable counters."""
    if not path.is_file():
        raise RuntimeError("Missing VTU piece: {}".format(path))

    root = ET.parse(str(path)).getroot()
    grid = root.find("UnstructuredGrid")
    if grid is None:
        raise RuntimeError("{} is not an UnstructuredGrid file".format(path))

    piece = grid.find("Piece")
    if piece is None:
        raise RuntimeError("No Piece element in {}".format(path))

    number_of_points = int(piece.attrib["NumberOfPoints"])
    number_of_cells = int(piece.attrib["NumberOfCells"])

    point_arrays = named_arrays(piece.find("PointData"))
    cell_arrays = named_arrays(piece.find("CellData"))

    missing_point = REQUIRED_POINT_ARRAYS.difference(point_arrays)
    missing_cell = REQUIRED_CELL_ARRAYS.difference(cell_arrays)
    if missing_point or missing_cell:
        raise RuntimeError(
            "{}: missing point arrays={}, cell arrays={}".format(
                path.name,
                sorted(missing_point),
                sorted(missing_cell),
            )
        )

    cell_kind = parse_numbers(cell_arrays["cell_kind"].text, int)
    global_element = parse_numbers(
        cell_arrays["global_element_id"].text,
        int,
    )
    material_id = parse_numbers(cell_arrays["material_id"].text, int)
    bc_id = parse_numbers(cell_arrays["bc_id"].text, int)
    parent_global = parse_numbers(
        cell_arrays["parent_global_element"].text,
        int,
    )

    arrays = [cell_kind, global_element, material_id, bc_id, parent_global]
    if any(len(values) != number_of_cells for values in arrays):
        lengths = [len(values) for values in arrays]
        raise RuntimeError(
            "{}: cell-array lengths {} do not match NumberOfCells={}".format(
                path.name,
                lengths,
                number_of_cells,
            )
        )

    volume_material = {}  # type: Dict[int, int]
    material_counts = collections.Counter()
    bc_counts = collections.Counter()
    parent_material_counts = collections.Counter()

    for kind, gid, material in zip(
        cell_kind,
        global_element,
        material_id,
    ):
        if kind == 1:
            volume_material[gid] = material
            material_counts[material] += 1

    for kind, bc, parent in zip(cell_kind, bc_id, parent_global):
        if kind != 2:
            continue

        bc_counts[bc] += 1
        parent_material = volume_material.get(parent)
        parent_material_counts[(bc, parent_material)] += 1

    ghost_values = parse_numbers(
        point_arrays["vtkGhostType"].text,
        int,
    )
    if len(ghost_values) != number_of_points:
        raise RuntimeError(
            "{}: vtkGhostType length {} does not match NumberOfPoints={}".format(
                path.name,
                len(ghost_values),
                number_of_points,
            )
        )

    return (
        material_counts,
        bc_counts,
        parent_material_counts,
        set(point_arrays),
        set(cell_arrays),
        number_of_points,
        sum(1 for value in ghost_values if value != 0),
    )


def data_lines(path, limit=30):
    """Read the first non-comment data lines from a text input file."""
    lines = []  # type: List[str]
    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for raw in stream:
            line = raw.strip()
            if not line or line.startswith("#") or line.startswith("!"):
                continue

            lines.append(line)
            if len(lines) >= limit:
                break

    return lines


def audit_material_files(partition_root):
    """Report whether all rank-local material property files are identical."""
    files = sorted(partition_root.rglob("*.matprop"))

    print("\n===== MATERIAL PROPERTY FILES =====")
    print("matprop files found : {}".format(len(files)))

    if not files:
        print("ERROR: no .matprop files found")
        return

    hashes = collections.defaultdict(list)
    for path in files:
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        hashes[digest].append(path)

    print("unique contents     : {}".format(len(hashes)))

    for index, paths in enumerate(hashes.values(), start=1):
        representative = paths[0]
        print("\nmatprop variant {}: {}".format(index, representative))
        for line in data_lines(representative):
            print("  {}".format(line))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pvtu", required=True, type=Path)
    parser.add_argument("--partition-root", required=True, type=Path)
    args = parser.parse_args()

    pvtu = args.pvtu.resolve()
    partition_root = args.partition_root.resolve()

    pieces = find_piece_sources(pvtu)
    total_material = collections.Counter()
    total_bc = collections.Counter()
    total_parent_material = collections.Counter()
    point_array_intersection = None  # type: Optional[set]
    cell_array_intersection = None  # type: Optional[set]
    total_local_points = 0
    total_ghost_points = 0

    print("CBS3D heat-transfer post-processing audit")
    print("PVTU             : {}".format(pvtu))
    print("VTU pieces       : {}".format(len(pieces)))
    print("Partition root   : {}".format(partition_root))

    for index, piece_path in enumerate(pieces, start=1):
        (
            material_counts,
            bc_counts,
            parent_material_counts,
            point_arrays,
            cell_arrays,
            local_points,
            ghost_points,
        ) = audit_piece(piece_path)

        total_material.update(material_counts)
        total_bc.update(bc_counts)
        total_parent_material.update(parent_material_counts)
        total_local_points += local_points
        total_ghost_points += ghost_points

        point_array_intersection = (
            point_arrays
            if point_array_intersection is None
            else point_array_intersection.intersection(point_arrays)
        )
        cell_array_intersection = (
            cell_arrays
            if cell_array_intersection is None
            else cell_array_intersection.intersection(cell_arrays)
        )

        print(
            "  audited piece {:02d}/{:02d}: {}".format(
                index,
                len(pieces),
                piece_path.name,
            )
        )

    print("\n===== COMMON ARRAYS =====")
    print("PointData:")
    for name in sorted(point_array_intersection or []):
        print("  {}".format(name))

    print("CellData:")
    for name in sorted(cell_array_intersection or []):
        print("  {}".format(name))

    print("\n===== VOLUME MATERIAL COUNTS =====")
    for material, count in sorted(total_material.items()):
        print("material_id {:6d} : {}".format(material, count))

    print("\n===== BOUNDARY FACE COUNTS =====")
    for bc, count in sorted(total_bc.items()):
        print("bc_id {:6d} : {}".format(bc, count))

    print("\n===== BOUNDARY PARENT MATERIAL COUNTS =====")
    for (bc, material), count in sorted(
        total_parent_material.items(),
        key=lambda item: (
            item[0][0],
            -999999 if item[0][1] is None else item[0][1],
        ),
    ):
        print(
            "bc_id {:6d}, parent material {:>6s} : {}".format(
                bc,
                str(material),
                count,
            )
        )

    print("\n===== POINT OWNERSHIP SUMMARY =====")
    print("sum of rank-local points : {}".format(total_local_points))
    print("sum of ghost points      : {}".format(total_ghost_points))

    audit_material_files(partition_root)

    errors = []
    if total_bc[901] == 0:
        errors.append("no BC 901 interface faces were found")
    if total_bc[532] == 0:
        errors.append("no BC 532 imposed-heat-flux faces were found")
    if not list(partition_root.rglob("*.matprop")):
        errors.append("no material-property file was found")

    print("\n===== AUDIT RESULT =====")
    if errors:
        for error in errors:
            print("FAIL: {}".format(error))
        return 1

    print(
        "PASS: topology and material inputs needed for CHT post-processing exist"
    )
    print(
        "Next: reconstruct P1 tetrahedral temperature gradients and interface fluxes"
    )
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        print("ERROR: {}".format(error), file=sys.stderr)
        sys.exit(2)
