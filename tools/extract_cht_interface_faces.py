#!/usr/bin/env python3
"""Extract the implicit conformal fluid-solid interface from distributed VTU.

The CBS3D distributed output stores only external boundary triangles explicitly.
The fluid-solid interface is represented by conformal tetrahedra with different
material identifiers sharing the same three global node identifiers. This tool
reconstructs that interface exactly from the P1 tetrahedral topology.

The implementation is compatible with the older Python available on Sunbird.
It processes one VTU piece at a time, resolves faces internal to each partition
in memory, and externally sorts only the unmatched partition-boundary faces.
"""

import argparse
import collections
import csv
import math
import os
import re
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


TETRA_FACES = (
    (0, 1, 2),
    (0, 1, 3),
    (0, 2, 3),
    (1, 2, 3),
)

CSV_HEADER = (
    "node_1",
    "node_2",
    "node_3",
    "fluid_element_id",
    "solid_element_id",
    "solid_material_id",
    "fluid_rank",
    "fluid_local_cell",
    "solid_rank",
    "solid_local_cell",
    "area_m2",
    "centroid_x_m",
    "centroid_y_m",
    "centroid_z_m",
    "node_mask_ok",
    "pair_location",
)


def parse_numbers(text, cast):
    if not text:
        return []
    return [cast(token) for token in text.split()]


def named_arrays(parent):
    result = {}
    if parent is None:
        return result
    for array in parent.findall("DataArray"):
        name = array.attrib.get("Name")
        if name:
            result[name] = array
    return result


def find_piece_sources(pvtu):
    root = ET.parse(str(pvtu)).getroot()
    grid = root.find("PUnstructuredGrid")
    if grid is None:
        raise RuntimeError("{} is not a PUnstructuredGrid file".format(pvtu))

    pieces = []
    for piece in grid.findall("Piece"):
        source = piece.attrib.get("Source")
        if not source:
            raise RuntimeError("PVTU piece without Source in {}".format(pvtu))
        path = (pvtu.parent / source).resolve()
        if not path.is_file():
            raise RuntimeError("Missing VTU piece: {}".format(path))
        pieces.append(path)

    if not pieces:
        raise RuntimeError("No Piece entries found in {}".format(pvtu))
    return pieces


def rank_from_piece_name(path, fallback):
    match = re.search(r"_rank_(\d+)\.vtu$", path.name)
    if match:
        return int(match.group(1))
    return fallback


def vector_subtract(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def cross_product(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def face_geometry(key, coordinate_by_gid, mask_by_gid):
    try:
        p0 = coordinate_by_gid[key[0]]
        p1 = coordinate_by_gid[key[1]]
        p2 = coordinate_by_gid[key[2]]
    except KeyError as error:
        raise RuntimeError(
            "Missing coordinate for global face node {}".format(error)
        )

    edge1 = vector_subtract(p1, p0)
    edge2 = vector_subtract(p2, p0)
    cross = cross_product(edge1, edge2)
    area = 0.5 * math.sqrt(
        cross[0] * cross[0]
        + cross[1] * cross[1]
        + cross[2] * cross[2]
    )

    if not math.isfinite(area) or area <= 0.0:
        raise RuntimeError(
            "Non-positive interface-candidate face area for nodes {}".format(key)
        )

    centroid = (
        (p0[0] + p1[0] + p2[0]) / 3.0,
        (p0[1] + p1[1] + p2[1]) / 3.0,
        (p0[2] + p1[2] + p2[2]) / 3.0,
    )
    mask_ok = int(
        mask_by_gid.get(key[0]) == 3
        and mask_by_gid.get(key[1]) == 3
        and mask_by_gid.get(key[2]) == 3
    )
    return (area, centroid[0], centroid[1], centroid[2], mask_ok)


def write_spill_record(stream, key, record, geometry):
    material, element_id, rank, local_cell = record
    area, cx, cy, cz, mask_ok = geometry
    fields = (
        key[0], key[1], key[2],
        material, element_id, rank, local_cell,
        "{:.17g}".format(area),
        "{:.17g}".format(cx),
        "{:.17g}".format(cy),
        "{:.17g}".format(cz),
        mask_ok,
    )
    stream.write("\t".join(str(value) for value in fields) + "\n")


def parse_spill_record(line):
    fields = line.rstrip("\n").split("\t")
    if len(fields) != 12:
        raise RuntimeError(
            "Malformed unmatched-face record with {} fields".format(len(fields))
        )
    key = (int(fields[0]), int(fields[1]), int(fields[2]))
    record = (
        int(fields[3]),
        int(fields[4]),
        int(fields[5]),
        int(fields[6]),
    )
    geometry = (
        float(fields[7]),
        float(fields[8]),
        float(fields[9]),
        float(fields[10]),
        int(fields[11]),
    )
    return key, record, geometry


def geometry_consistent(a, b):
    area_scale = max(abs(a[0]), abs(b[0]), 1.0e-300)
    if abs(a[0] - b[0]) > 1.0e-11 * area_scale:
        return False
    coordinate_scale = max(
        abs(a[1]), abs(a[2]), abs(a[3]),
        abs(b[1]), abs(b[2]), abs(b[3]),
        1.0,
    )
    tolerance = 1.0e-12 * coordinate_scale
    return (
        abs(a[1] - b[1]) <= tolerance
        and abs(a[2] - b[2]) <= tolerance
        and abs(a[3] - b[3]) <= tolerance
    )


def write_interface_pair(stream, key, first, second, geometry, location):
    material_a, element_a, rank_a, cell_a = first
    material_b, element_b, rank_b, cell_b = second

    if material_a == 0 and material_b > 0:
        fluid = first
        solid = second
    elif material_b == 0 and material_a > 0:
        fluid = second
        solid = first
    else:
        raise RuntimeError(
            "Face {} is not a fluid-solid material pair: {} and {}".format(
                key, material_a, material_b
            )
        )

    _, fluid_element, fluid_rank, fluid_cell = fluid
    solid_material, solid_element, solid_rank, solid_cell = solid
    area, cx, cy, cz, mask_ok = geometry

    fields = (
        key[0], key[1], key[2],
        fluid_element, solid_element, solid_material,
        fluid_rank, fluid_cell, solid_rank, solid_cell,
        "{:.17g}".format(area),
        "{:.17g}".format(cx),
        "{:.17g}".format(cy),
        "{:.17g}".format(cz),
        mask_ok,
        location,
    )
    stream.write("\t".join(str(value) for value in fields) + "\n")


def run_sort(source, destination, temp_directory, parallelism):
    command = [
        "sort",
        "-T", str(temp_directory),
        "-S", "50%",
        "--parallel={}".format(max(1, parallelism)),
        "-t", "\t",
        "-k1,1n",
        "-k2,2n",
        "-k3,3n",
        str(source),
        "-o", str(destination),
    ]
    environment = dict(os.environ)
    environment["LC_ALL"] = "C"
    subprocess.check_call(command, env=environment)


def read_piece(path, fallback_rank, owned_global_nodes, owned_interface_nodes):
    root = ET.parse(str(path)).getroot()
    piece = root.find("./UnstructuredGrid/Piece")
    if piece is None:
        raise RuntimeError("No UnstructuredGrid/Piece in {}".format(path))

    number_of_points = int(piece.attrib["NumberOfPoints"])
    number_of_cells = int(piece.attrib["NumberOfCells"])

    points_array = piece.find("./Points/DataArray")
    if points_array is None:
        raise RuntimeError("Missing Points/DataArray in {}".format(path))
    point_values = parse_numbers(points_array.text, float)
    if len(point_values) != 3 * number_of_points:
        raise RuntimeError(
            "{}: coordinate count {} does not equal 3*{}".format(
                path.name, len(point_values), number_of_points
            )
        )

    point_arrays = named_arrays(piece.find("PointData"))
    cell_arrays = named_arrays(piece.find("CellData"))
    cells_arrays = named_arrays(piece.find("Cells"))

    required_point = {"global_node_id", "node_domain_kind", "is_owned"}
    required_cell = {
        "cell_kind", "global_element_id", "material_id"
    }
    required_cells = {"connectivity", "offsets", "types"}

    missing_point = required_point.difference(point_arrays)
    missing_cell = required_cell.difference(cell_arrays)
    missing_cells = required_cells.difference(cells_arrays)
    if missing_point or missing_cell or missing_cells:
        raise RuntimeError(
            "{}: missing PointData={}, CellData={}, Cells={}".format(
                path.name,
                sorted(missing_point),
                sorted(missing_cell),
                sorted(missing_cells),
            )
        )

    global_ids = parse_numbers(point_arrays["global_node_id"].text, int)
    domain_kind = parse_numbers(point_arrays["node_domain_kind"].text, int)
    is_owned = parse_numbers(point_arrays["is_owned"].text, int)
    if not (
        len(global_ids) == number_of_points
        and len(domain_kind) == number_of_points
        and len(is_owned) == number_of_points
    ):
        raise RuntimeError("{}: inconsistent point-array lengths".format(path.name))

    coordinates = []
    coordinate_by_gid = {}
    mask_by_gid = {}
    for index, gid in enumerate(global_ids):
        coordinate = (
            point_values[3 * index],
            point_values[3 * index + 1],
            point_values[3 * index + 2],
        )
        coordinates.append(coordinate)
        coordinate_by_gid[gid] = coordinate
        mask_by_gid[gid] = domain_kind[index]

        if is_owned[index] == 1:
            if gid in owned_global_nodes:
                raise RuntimeError(
                    "Global node {} has more than one owner".format(gid)
                )
            owned_global_nodes.add(gid)
            if domain_kind[index] == 3:
                owned_interface_nodes.add(gid)

    connectivity = parse_numbers(cells_arrays["connectivity"].text, int)
    offsets = parse_numbers(cells_arrays["offsets"].text, int)
    cell_types = parse_numbers(cells_arrays["types"].text, int)
    cell_kind = parse_numbers(cell_arrays["cell_kind"].text, int)
    element_ids = parse_numbers(cell_arrays["global_element_id"].text, int)
    material_ids = parse_numbers(cell_arrays["material_id"].text, int)

    for values, label in (
        (offsets, "offsets"),
        (cell_types, "types"),
        (cell_kind, "cell_kind"),
        (element_ids, "global_element_id"),
        (material_ids, "material_id"),
    ):
        if len(values) != number_of_cells:
            raise RuntimeError(
                "{}: {} length {} does not equal NumberOfCells={}".format(
                    path.name, label, len(values), number_of_cells
                )
            )

    rank = rank_from_piece_name(path, fallback_rank)
    face_buckets = collections.defaultdict(list)
    material_counts = collections.Counter()
    volume_cells = 0
    boundary_cells = 0
    previous_offset = 0

    for cell_index in range(number_of_cells):
        end_offset = offsets[cell_index]
        local_nodes = connectivity[previous_offset:end_offset]
        previous_offset = end_offset

        if cell_kind[cell_index] == 1:
            volume_cells += 1
            if cell_types[cell_index] != 10 or len(local_nodes) != 4:
                raise RuntimeError(
                    "{}: volume cell {} is not a four-node VTK tetrahedron".format(
                        path.name, cell_index + 1
                    )
                )

            material = material_ids[cell_index]
            element_id = element_ids[cell_index]
            if element_id < 0:
                raise RuntimeError(
                    "{}: invalid global element id at cell {}".format(
                        path.name, cell_index + 1
                    )
                )
            material_counts[material] += 1

            gids = tuple(global_ids[node] for node in local_nodes)
            record = (material, element_id, rank, cell_index + 1)
            for face in TETRA_FACES:
                key = tuple(sorted((gids[face[0]], gids[face[1]], gids[face[2]])))
                face_buckets[key].append(record)
        elif cell_kind[cell_index] == 2:
            boundary_cells += 1
        else:
            raise RuntimeError(
                "{}: unsupported cell_kind {} at cell {}".format(
                    path.name, cell_kind[cell_index], cell_index + 1
                )
            )

    if previous_offset != len(connectivity):
        raise RuntimeError(
            "{}: final cell offset {} does not equal connectivity length {}".format(
                path.name, previous_offset, len(connectivity)
            )
        )

    return {
        "rank": rank,
        "face_buckets": face_buckets,
        "coordinate_by_gid": coordinate_by_gid,
        "mask_by_gid": mask_by_gid,
        "material_counts": material_counts,
        "volume_cells": volume_cells,
        "boundary_cells": boundary_cells,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pvtu", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--sort-parallel",
        type=int,
        default=1,
        help="GNU sort worker count; normally match Slurm cpus-per-task",
    )
    parser.add_argument("--keep-temp", action="store_true")
    args = parser.parse_args()

    pvtu = args.pvtu.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    stem = pvtu.stem
    final_csv = output_dir / (stem + "_interface_faces.csv")
    summary_path = output_dir / (stem + "_interface_summary.txt")
    temp_dir = output_dir / ("." + stem + "_interface_tmp")

    if temp_dir.exists():
        shutil.rmtree(str(temp_dir))
    temp_dir.mkdir(parents=True)

    spill_unsorted = temp_dir / "unmatched_faces.tsv"
    spill_sorted = temp_dir / "unmatched_faces_sorted.tsv"
    interface_unsorted = temp_dir / "interface_faces.tsv"
    interface_sorted = temp_dir / "interface_faces_sorted.tsv"

    pieces = find_piece_sources(pvtu)
    owned_global_nodes = set()
    owned_interface_nodes = set()
    total_material = collections.Counter()
    total_volume_cells = 0
    total_boundary_cells = 0
    local_same_material = 0
    local_interface = 0
    local_nonmanifold = 0

    print("CBS3D conformal CHT interface extraction")
    print("PVTU          : {}".format(pvtu))
    print("VTU pieces    : {}".format(len(pieces)))
    print("Output        : {}".format(final_csv))
    print("Temporary dir : {}".format(temp_dir))

    with spill_unsorted.open("w", encoding="utf-8") as spill_stream, \
            interface_unsorted.open("w", encoding="utf-8") as interface_stream:
        for index, piece_path in enumerate(pieces):
            data = read_piece(
                piece_path,
                index,
                owned_global_nodes,
                owned_interface_nodes,
            )
            total_material.update(data["material_counts"])
            total_volume_cells += data["volume_cells"]
            total_boundary_cells += data["boundary_cells"]

            for key, records in data["face_buckets"].items():
                if len(records) == 1:
                    geometry = face_geometry(
                        key,
                        data["coordinate_by_gid"],
                        data["mask_by_gid"],
                    )
                    write_spill_record(spill_stream, key, records[0], geometry)
                elif len(records) == 2:
                    if records[0][0] == records[1][0]:
                        local_same_material += 1
                    else:
                        geometry = face_geometry(
                            key,
                            data["coordinate_by_gid"],
                            data["mask_by_gid"],
                        )
                        write_interface_pair(
                            interface_stream,
                            key,
                            records[0],
                            records[1],
                            geometry,
                            "local",
                        )
                        local_interface += 1
                else:
                    local_nonmanifold += 1

            print(
                "  processed piece {:02d}/{:02d}: rank {:04d}, "
                "tetrahedra={}, boundary_faces={}, unmatched_faces={}".format(
                    index + 1,
                    len(pieces),
                    data["rank"],
                    data["volume_cells"],
                    data["boundary_cells"],
                    sum(1 for records in data["face_buckets"].values() if len(records) == 1),
                )
            )

            del data

    print("Sorting unmatched partition faces...")
    run_sort(spill_unsorted, spill_sorted, temp_dir, args.sort_parallel)

    cross_same_material = 0
    cross_interface = 0
    external_faces = 0
    global_nonmanifold = 0
    geometry_mismatches = 0

    with spill_sorted.open("r", encoding="utf-8") as source, \
            interface_unsorted.open("a", encoding="utf-8") as interface_stream:
        current_key = None
        group = []

        def process_group(key, records):
            result = {
                "same": 0,
                "interface": 0,
                "external": 0,
                "nonmanifold": 0,
                "geometry_mismatch": 0,
            }
            if key is None:
                return result
            if len(records) == 1:
                result["external"] = 1
                return result
            if len(records) != 2:
                result["nonmanifold"] = 1
                return result

            first_record, first_geometry = records[0]
            second_record, second_geometry = records[1]
            if not geometry_consistent(first_geometry, second_geometry):
                result["geometry_mismatch"] = 1

            if first_record[0] == second_record[0]:
                result["same"] = 1
            else:
                combined_geometry = (
                    first_geometry[0],
                    first_geometry[1],
                    first_geometry[2],
                    first_geometry[3],
                    int(first_geometry[4] == 1 and second_geometry[4] == 1),
                )
                write_interface_pair(
                    interface_stream,
                    key,
                    first_record,
                    second_record,
                    combined_geometry,
                    "cross_rank",
                )
                result["interface"] = 1
            return result

        for line in source:
            key, record, geometry = parse_spill_record(line)
            if current_key is None:
                current_key = key
            if key != current_key:
                result = process_group(current_key, group)
                cross_same_material += result["same"]
                cross_interface += result["interface"]
                external_faces += result["external"]
                global_nonmanifold += result["nonmanifold"]
                geometry_mismatches += result["geometry_mismatch"]
                current_key = key
                group = []
            group.append((record, geometry))

        result = process_group(current_key, group)
        cross_same_material += result["same"]
        cross_interface += result["interface"]
        external_faces += result["external"]
        global_nonmanifold += result["nonmanifold"]
        geometry_mismatches += result["geometry_mismatch"]

    print("Sorting extracted interface faces...")
    run_sort(interface_unsorted, interface_sorted, temp_dir, args.sort_parallel)

    interface_faces = 0
    interface_nodes = set()
    interface_area = 0.0
    minimum_area = None
    maximum_area = 0.0
    bad_mask_faces = 0
    duplicate_interface_faces = 0
    previous_key = None

    with interface_sorted.open("r", encoding="utf-8") as source, \
            final_csv.open("w", encoding="utf-8", newline="") as destination:
        writer = csv.writer(destination)
        writer.writerow(CSV_HEADER)
        for line in source:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(CSV_HEADER):
                raise RuntimeError(
                    "Malformed interface record with {} fields".format(len(fields))
                )
            key = (int(fields[0]), int(fields[1]), int(fields[2]))
            if key == previous_key:
                duplicate_interface_faces += 1
            previous_key = key

            area = float(fields[10])
            mask_ok = int(fields[14])
            interface_faces += 1
            interface_nodes.update(key)
            interface_area += area
            minimum_area = area if minimum_area is None else min(minimum_area, area)
            maximum_area = max(maximum_area, area)
            if mask_ok != 1:
                bad_mask_faces += 1
            writer.writerow(fields)

    total_internal_same = local_same_material + cross_same_material
    total_interface = local_interface + cross_interface
    incidence_accounted = (
        2 * (total_internal_same + total_interface) + external_faces
    )
    incidence_expected = 4 * total_volume_cells
    incidence_residual = incidence_expected - incidence_accounted
    missing_interface_nodes = owned_interface_nodes.difference(interface_nodes)
    extra_interface_nodes = interface_nodes.difference(owned_interface_nodes)

    lines = []
    lines.append("CBS3D conformal CHT interface extraction summary")
    lines.append("PVTU                         : {}".format(pvtu))
    lines.append("VTU pieces                   : {}".format(len(pieces)))
    lines.append("global owned nodes           : {}".format(len(owned_global_nodes)))
    lines.append("owned mask-3 interface nodes : {}".format(len(owned_interface_nodes)))
    lines.append("")
    lines.append("Volume material counts")
    for material, count in sorted(total_material.items()):
        lines.append("  material {:6d}              : {}".format(material, count))
    lines.append("  total tetrahedra            : {}".format(total_volume_cells))
    lines.append("  explicit boundary triangles : {}".format(total_boundary_cells))
    lines.append("")
    lines.append("Face reconstruction")
    lines.append("  same-material local faces   : {}".format(local_same_material))
    lines.append("  same-material cross-rank    : {}".format(cross_same_material))
    lines.append("  interface local faces       : {}".format(local_interface))
    lines.append("  interface cross-rank faces  : {}".format(cross_interface))
    lines.append("  total interface faces       : {}".format(interface_faces))
    lines.append("  external unmatched faces    : {}".format(external_faces))
    lines.append("  local non-manifold groups   : {}".format(local_nonmanifold))
    lines.append("  global non-manifold groups  : {}".format(global_nonmanifold))
    lines.append("  geometry mismatches         : {}".format(geometry_mismatches))
    lines.append("  duplicate interface faces   : {}".format(duplicate_interface_faces))
    lines.append("  tetra-face incidence error  : {}".format(incidence_residual))
    lines.append("")
    lines.append("Interface geometry")
    lines.append("  unique interface nodes      : {}".format(len(interface_nodes)))
    lines.append("  missing mask-3 nodes        : {}".format(len(missing_interface_nodes)))
    lines.append("  unexpected interface nodes  : {}".format(len(extra_interface_nodes)))
    lines.append("  bad node-mask faces         : {}".format(bad_mask_faces))
    lines.append("  total interface area [m^2]  : {:.16e}".format(interface_area))
    lines.append("  minimum face area [m^2]     : {:.16e}".format(minimum_area or 0.0))
    lines.append("  maximum face area [m^2]     : {:.16e}".format(maximum_area))
    lines.append("")
    lines.append("Output CSV                    : {}".format(final_csv))

    errors = []
    if local_nonmanifold != 0 or global_nonmanifold != 0:
        errors.append("non-manifold tetrahedral face groups were found")
    if geometry_mismatches != 0:
        errors.append("cross-rank copies disagree geometrically")
    if duplicate_interface_faces != 0:
        errors.append("duplicate interface faces were written")
    if incidence_residual != 0:
        errors.append("tetrahedral face-incidence identity does not close")
    if external_faces != total_boundary_cells:
        errors.append(
            "reconstructed external face count {} differs from explicit boundary count {}".format(
                external_faces, total_boundary_cells
            )
        )
    if interface_faces == 0:
        errors.append("no fluid-solid interface faces were reconstructed")
    if bad_mask_faces != 0:
        errors.append("some reconstructed material interfaces do not use three mask-3 nodes")
    if missing_interface_nodes or extra_interface_nodes:
        errors.append("interface-face nodes do not exactly match owned mask-3 nodes")

    lines.append("")
    lines.append("===== EXTRACTION RESULT =====")
    if errors:
        for error in errors:
            lines.append("FAIL: {}".format(error))
    else:
        lines.append("PASS: implicit fluid-solid interface reconstructed exactly")
        lines.append("Next: reconstruct fluid- and solid-side P1 temperature gradients")

    summary = "\n".join(lines) + "\n"
    summary_path.write_text(summary, encoding="utf-8")
    print("\n" + summary)

    if not args.keep_temp:
        shutil.rmtree(str(temp_dir))

    return 1 if errors else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        print("ERROR: {}".format(error), file=sys.stderr)
        sys.exit(2)
