#!/usr/bin/env python3
"""Audit a distributed CBS3D VTU/PVTU result before computing h and Nu.

The audit is intentionally dependency-free. It reads the ASCII VTK XML pieces
sequentially and verifies that the topology and identifiers required for an
accurate conjugate heat-transfer post-processor are present.
"""

from __future__ import annotations

import argparse
import collections
import hashlib
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


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


def parse_numbers(text: str | None, cast):
    if not text:
        return []
    return [cast(token) for token in text.split()]


def named_arrays(parent: ET.Element | None) -> Dict[str, ET.Element]:
    if parent is None:
        return {}
    result: Dict[str, ET.Element] = {}
    for array in parent.findall("DataArray"):
        name = array.attrib.get("Name")
        if name:
            result[name] = array
    return result


def find_piece_sources(pvtu: Path) -> List[Path]:
    root = ET.parse(pvtu).getroot()
    grid = root.find("PUnstructuredGrid")
    if grid is None:
        raise RuntimeError(f"{pvtu} is not a PUnstructuredGrid file")

    pieces = []
    for piece in grid.findall("Piece"):
        source = piece.attrib.get("Source")
        if not source:
            raise RuntimeError(f"PVTU piece without Source in {pvtu}")
        pieces.append((pvtu.parent / source).resolve())

    if not pieces:
        raise RuntimeError(f"No Piece entries found in {pvtu}")
    return pieces


def audit_piece(
    path: Path,
) -> Tuple[
    collections.Counter,
    collections.Counter,
    collections.Counter,
    set,
    set,
    int,
    int,
]:
    if not path.is_file():
        raise RuntimeError(f"Missing VTU piece: {path}")

    root = ET.parse(path).getroot()
    grid = root.find("UnstructuredGrid")
    if grid is None:
        raise RuntimeError(f"{path} is not an UnstructuredGrid file")
    piece = grid.find("Piece")
    if piece is None:
        raise RuntimeError(f"No Piece element in {path}")

    number_of_points = int(piece.attrib["NumberOfPoints"])
    number_of_cells = int(piece.attrib["NumberOfCells"])

    point_arrays = named_arrays(piece.find("PointData"))
    cell_arrays = named_arrays(piece.find("CellData"))

    missing_point = REQUIRED_POINT_ARRAYS.difference(point_arrays)
    missing_cell = REQUIRED_CELL_ARRAYS.difference(cell_arrays)
    if missing_point or missing_cell:
        raise RuntimeError(
            f"{path.name}: missing point arrays={sorted(missing_point)}, "
            f"cell arrays={sorted(missing_cell)}"
        )

    cell_kind = parse_numbers(cell_arrays["cell_kind"].text, int)
    global_element = parse_numbers(cell_arrays["global_element_id"].text, int)
    material_id = parse_numbers(cell_arrays["material_id"].text, int)
    bc_id = parse_numbers(cell_arrays["bc_id"].text, int)
    parent_global = parse_numbers(
        cell_arrays["parent_global_element"].text, int
    )

    arrays = [cell_kind, global_element, material_id, bc_id, parent_global]
    if any(len(values) != number_of_cells for values in arrays):
        lengths = [len(values) for values in arrays]
        raise RuntimeError(
            f"{path.name}: cell-array lengths {lengths} do not match "
            f"NumberOfCells={number_of_cells}"
        )

    volume_material: Dict[int, int] = {}
    material_counts = collections.Counter()
    bc_counts = collections.Counter()
    parent_material_counts = collections.Counter()

    for kind, gid, material in zip(
        cell_kind, global_element, material_id
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
        point_arrays["vtkGhostType"].text, int
    )
    if len(ghost_values) != number_of_points:
        raise RuntimeError(
            f"{path.name}: vtkGhostType length {len(ghost_values)} does not "
            f"match NumberOfPoints={number_of_points}"
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


def data_lines(path: Path, limit: int = 30) -> List[str]:
    lines: List[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for raw in stream:
            line = raw.strip()
            if not line or line.startswith("#") or line.startswith("!"):
                continue
            lines.append(line)
            if len(lines) >= limit:
                break
    return lines


def audit_material_files(partition_root: Path) -> None:
    files = sorted(partition_root.rglob("*.matprop"))
    print("\n===== MATERIAL PROPERTY FILES =====")
    print(f"matprop files found : {len(files)}")

    if not files:
        print("ERROR: no .matprop files found")
        return

    hashes = collections.defaultdict(list)
    for path in files:
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        hashes[digest].append(path)

    print(f"unique contents     : {len(hashes)}")
    for index, paths in enumerate(hashes.values(), start=1):
        representative = paths[0]
        print(f"\nmatprop variant {index}: {representative}")
        for line in data_lines(representative):
            print(f"  {line}")


def main() -> int:
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
    point_array_intersection = None
    cell_array_intersection = None
    total_local_points = 0
    total_ghost_points = 0

    print("CBS3D heat-transfer post-processing audit")
    print(f"PVTU             : {pvtu}")
    print(f"VTU pieces       : {len(pieces)}")
    print(f"Partition root   : {partition_root}")

    for index, piece in enumerate(pieces, start=1):
        (
            material_counts,
            bc_counts,
            parent_material_counts,
            point_arrays,
            cell_arrays,
            local_points,
            ghost_points,
        ) = audit_piece(piece)

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
            f"  audited piece {index:02d}/{len(pieces):02d}: "
            f"{piece.name}"
        )

    print("\n===== COMMON ARRAYS =====")
    print("PointData:")
    for name in sorted(point_array_intersection or []):
        print(f"  {name}")
    print("CellData:")
    for name in sorted(cell_array_intersection or []):
        print(f"  {name}")

    print("\n===== VOLUME MATERIAL COUNTS =====")
    for material, count in sorted(total_material.items()):
        print(f"material_id {material:6d} : {count}")

    print("\n===== BOUNDARY FACE COUNTS =====")
    for bc, count in sorted(total_bc.items()):
        print(f"bc_id {bc:6d} : {count}")

    print("\n===== BOUNDARY PARENT MATERIAL COUNTS =====")
    for (bc, material), count in sorted(
        total_parent_material.items(),
        key=lambda item: (item[0][0], -999999 if item[0][1] is None else item[0][1]),
    ):
        print(f"bc_id {bc:6d}, parent material {material!s:>6s} : {count}")

    print("\n===== POINT OWNERSHIP SUMMARY =====")
    print(f"sum of rank-local points : {total_local_points}")
    print(f"sum of ghost points      : {total_ghost_points}")

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
            print(f"FAIL: {error}")
        return 1

    print("PASS: topology and material inputs needed for CHT post-processing exist")
    print("Next: reconstruct P1 tetrahedral temperature gradients and interface fluxes")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(2)
