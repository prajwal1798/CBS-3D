#!/usr/bin/env python3
"""Reconstruct one-sided P1 temperature gradients and CHT interface heat flux.

The CBS3D distributed VTU output contains the conformal fluid-solid interface
implicitly. ``extract_cht_interface_faces.py`` first identifies each interface
triangle and records its adjacent fluid and solid tetrahedra. This tool then
loads those tetrahedra, reconstructs the constant P1 temperature gradient on
both sides, orients the interface unit normal from solid to fluid, and evaluates

    q''_fluid = -k_fluid grad(T_fluid) dot n_solid_to_fluid
    q''_solid = -k_solid grad(T_solid) dot n_solid_to_fluid

The two one-sided fluxes are deliberately retained separately. Standard
continuous P1 finite elements enforce temperature continuity but do not enforce
pointwise normal-flux continuity; their difference is therefore a discretisation
and gradient-recovery diagnostic, not automatically a solver failure.

The implementation avoids NumPy and modern type-hint syntax so it runs with the
older Python interpreter available on Swansea Sunbird.
"""

import argparse
import collections
import csv
import hashlib
import math
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


REQUIRED_INTERFACE_COLUMNS = {
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
    "pair_location",
}

OUTPUT_HEADER = (
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
    "pair_location",
    "area_m2",
    "centroid_x_m",
    "centroid_y_m",
    "centroid_z_m",
    "normal_solid_to_fluid_x",
    "normal_solid_to_fluid_y",
    "normal_solid_to_fluid_z",
    "wall_temperature_K",
    "wall_temperature_fluid_side_K",
    "wall_temperature_solid_side_K",
    "wall_temperature_side_mismatch_K",
    "gradT_fluid_x_K_per_m",
    "gradT_fluid_y_K_per_m",
    "gradT_fluid_z_K_per_m",
    "gradT_solid_x_K_per_m",
    "gradT_solid_y_K_per_m",
    "gradT_solid_z_K_per_m",
    "dTdn_fluid_K_per_m",
    "dTdn_solid_K_per_m",
    "k_fluid_W_per_mK",
    "k_solid_W_per_mK",
    "heat_flux_fluid_W_per_m2",
    "heat_flux_solid_W_per_m2",
    "heat_flux_arithmetic_mean_W_per_m2",
    "heat_flux_jump_fluid_minus_solid_W_per_m2",
    "heat_flux_relative_jump",
    "fluid_centroid_distance_m",
    "solid_centroid_distance_m",
    "area_relative_mismatch",
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


def vector_subtract(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def vector_add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def vector_scale(a, scalar):
    return (a[0] * scalar, a[1] * scalar, a[2] * scalar)


def dot_product(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross_product(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def vector_norm(a):
    return math.sqrt(dot_product(a, a))


def determinant3(
    a11, a12, a13,
    a21, a22, a23,
    a31, a32, a33,
):
    return (
        a11 * (a22 * a33 - a23 * a32)
        - a12 * (a21 * a33 - a23 * a31)
        + a13 * (a21 * a32 - a22 * a31)
    )


def rank_from_piece_name(path, fallback):
    match = re.search(r"_rank_(\d+)\.vtu$", path.name)
    if match:
        return int(match.group(1))
    return fallback


def find_piece_sources(pvtu):
    root = ET.parse(str(pvtu)).getroot()
    grid = root.find("PUnstructuredGrid")
    if grid is None:
        raise RuntimeError("{} is not a PUnstructuredGrid file".format(pvtu))

    result = {}
    for index, piece in enumerate(grid.findall("Piece")):
        source = piece.attrib.get("Source")
        if not source:
            raise RuntimeError("PVTU Piece without Source in {}".format(pvtu))
        path = (pvtu.parent / source).resolve()
        if not path.is_file():
            raise RuntimeError("Missing VTU piece: {}".format(path))
        rank = rank_from_piece_name(path, index)
        if rank in result:
            raise RuntimeError("Duplicate VTU rank {}".format(rank))
        result[rank] = path

    if not result:
        raise RuntimeError("No Piece entries found in {}".format(pvtu))
    return result


def read_material_properties(partition_root):
    files = sorted(partition_root.rglob("*.matprop"))
    if not files:
        raise RuntimeError("No .matprop file found under {}".format(partition_root))

    hashes = collections.defaultdict(list)
    for path in files:
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        hashes[digest].append(path)
    if len(hashes) != 1:
        raise RuntimeError(
            "Partition material-property files have {} distinct contents".format(
                len(hashes)
            )
        )

    representative = files[0]
    lines = []
    with representative.open("r", encoding="utf-8", errors="replace") as stream:
        for raw in stream:
            line = raw.strip()
            if not line or line.startswith("#") or line.startswith("!"):
                continue
            lines.append(line)

    if not lines:
        raise RuntimeError("Empty material-property file: {}".format(representative))

    try:
        expected = int(lines[0].split()[0])
    except Exception:
        raise RuntimeError(
            "Cannot read material count from {}".format(representative)
        )

    properties = {}
    for line in lines[1:]:
        fields = line.split()
        if len(fields) < 8:
            raise RuntimeError(
                "Malformed material-property row in {}: {}".format(
                    representative, line
                )
            )
        material_id = int(fields[0])
        properties[material_id] = {
            "name": fields[1],
            "phase": fields[2],
            "rho": float(fields[3]),
            "cp": float(fields[4]),
            "k": float(fields[5]),
            "mu": float(fields[6]),
            "source": float(fields[7]),
        }

    if len(properties) != expected:
        raise RuntimeError(
            "Material file declares {} records but {} were parsed".format(
                expected, len(properties)
            )
        )
    if 0 not in properties:
        raise RuntimeError("Material 0 (fluid) is absent from {}".format(representative))
    if properties[0]["phase"].lower() != "fluid":
        raise RuntimeError("Material 0 is not labelled as fluid")
    if properties[0]["k"] <= 0.0:
        raise RuntimeError("Fluid thermal conductivity must be positive")

    return properties, representative, len(files)


def read_interface_rows(path):
    rows = []
    needed_by_rank = collections.defaultdict(set)
    solid_material_ids = set()

    with path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise RuntimeError("Interface CSV has no header: {}".format(path))
        missing = REQUIRED_INTERFACE_COLUMNS.difference(reader.fieldnames)
        if missing:
            raise RuntimeError(
                "Interface CSV is missing columns: {}".format(sorted(missing))
            )

        for line_number, source in enumerate(reader, start=2):
            try:
                row = {
                    "nodes": (
                        int(source["node_1"]),
                        int(source["node_2"]),
                        int(source["node_3"]),
                    ),
                    "fluid_element_id": int(source["fluid_element_id"]),
                    "solid_element_id": int(source["solid_element_id"]),
                    "solid_material_id": int(source["solid_material_id"]),
                    "fluid_rank": int(source["fluid_rank"]),
                    "fluid_local_cell": int(source["fluid_local_cell"]),
                    "solid_rank": int(source["solid_rank"]),
                    "solid_local_cell": int(source["solid_local_cell"]),
                    "area": float(source["area_m2"]),
                    "centroid": (
                        float(source["centroid_x_m"]),
                        float(source["centroid_y_m"]),
                        float(source["centroid_z_m"]),
                    ),
                    "pair_location": source["pair_location"],
                }
            except Exception as error:
                raise RuntimeError(
                    "Malformed interface CSV row {}: {}".format(line_number, error)
                )

            if row["area"] <= 0.0 or not math.isfinite(row["area"]):
                raise RuntimeError(
                    "Non-positive interface area at CSV row {}".format(line_number)
                )
            if tuple(sorted(row["nodes"])) != row["nodes"]:
                raise RuntimeError(
                    "Interface node key is not sorted at CSV row {}".format(line_number)
                )

            rows.append(row)
            needed_by_rank[row["fluid_rank"]].add(row["fluid_local_cell"])
            needed_by_rank[row["solid_rank"]].add(row["solid_local_cell"])
            solid_material_ids.add(row["solid_material_id"])

    if not rows:
        raise RuntimeError("Interface CSV contains no faces: {}".format(path))
    return rows, needed_by_rank, solid_material_ids


def tetra_gradient(coordinates, temperatures, context):
    p1, p2, p3, p4 = coordinates

    j11 = p2[0] - p1[0]
    j21 = p2[1] - p1[1]
    j31 = p2[2] - p1[2]
    j12 = p3[0] - p1[0]
    j22 = p3[1] - p1[1]
    j32 = p3[2] - p1[2]
    j13 = p4[0] - p1[0]
    j23 = p4[1] - p1[1]
    j33 = p4[2] - p1[2]

    det_j = determinant3(
        j11, j12, j13,
        j21, j22, j23,
        j31, j32, j33,
    )
    if not math.isfinite(det_j) or det_j <= 0.0:
        raise RuntimeError(
            "Non-positive tetrahedron determinant {} for {}".format(det_j, context)
        )

    inv_det = 1.0 / det_j
    inv11 = (j22 * j33 - j23 * j32) * inv_det
    inv12 = (j13 * j32 - j12 * j33) * inv_det
    inv13 = (j12 * j23 - j13 * j22) * inv_det
    inv21 = (j23 * j31 - j21 * j33) * inv_det
    inv22 = (j11 * j33 - j13 * j31) * inv_det
    inv23 = (j13 * j21 - j11 * j23) * inv_det
    inv31 = (j21 * j32 - j22 * j31) * inv_det
    inv32 = (j12 * j31 - j11 * j32) * inv_det
    inv33 = (j11 * j22 - j12 * j21) * inv_det

    grad_n2 = (inv11, inv12, inv13)
    grad_n3 = (inv21, inv22, inv23)
    grad_n4 = (inv31, inv32, inv33)
    grad_n1 = (
        -(grad_n2[0] + grad_n3[0] + grad_n4[0]),
        -(grad_n2[1] + grad_n3[1] + grad_n4[1]),
        -(grad_n2[2] + grad_n3[2] + grad_n4[2]),
    )

    gradients = (grad_n1, grad_n2, grad_n3, grad_n4)
    gradient_t = (0.0, 0.0, 0.0)
    for temperature, gradient_n in zip(temperatures, gradients):
        gradient_t = vector_add(
            gradient_t,
            vector_scale(gradient_n, temperature),
        )

    if not all(math.isfinite(value) for value in gradient_t):
        raise RuntimeError("Non-finite temperature gradient for {}".format(context))

    centroid = (
        sum(point[0] for point in coordinates) / 4.0,
        sum(point[1] for point in coordinates) / 4.0,
        sum(point[2] for point in coordinates) / 4.0,
    )
    return gradient_t, centroid, det_j


def read_requested_tetrahedra(path, rank, needed_cells):
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
        raise RuntimeError("Coordinate length mismatch in {}".format(path))

    point_arrays = named_arrays(piece.find("PointData"))
    cell_arrays = named_arrays(piece.find("CellData"))
    cells_arrays = named_arrays(piece.find("Cells"))

    required_point = {"temperature", "global_node_id"}
    required_cell = {"cell_kind", "global_element_id", "material_id"}
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
    temperatures = parse_numbers(point_arrays["temperature"].text, float)
    if len(global_ids) != number_of_points or len(temperatures) != number_of_points:
        raise RuntimeError("PointData length mismatch in {}".format(path))

    coordinates = [
        (
            point_values[3 * index],
            point_values[3 * index + 1],
            point_values[3 * index + 2],
        )
        for index in range(number_of_points)
    ]

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

    result = {}
    for local_cell in sorted(needed_cells):
        cell_index = local_cell - 1
        if cell_index < 0 or cell_index >= number_of_cells:
            raise RuntimeError(
                "Rank {} requested local cell {} outside 1..{}".format(
                    rank, local_cell, number_of_cells
                )
            )
        if cell_kind[cell_index] != 1 or cell_types[cell_index] != 10:
            raise RuntimeError(
                "Rank {} local cell {} is not a tetrahedral volume cell".format(
                    rank, local_cell
                )
            )

        start = 0 if cell_index == 0 else offsets[cell_index - 1]
        end = offsets[cell_index]
        local_nodes = connectivity[start:end]
        if len(local_nodes) != 4:
            raise RuntimeError(
                "Rank {} local cell {} has {} nodes".format(
                    rank, local_cell, len(local_nodes)
                )
            )

        gids = tuple(global_ids[index] for index in local_nodes)
        coords = tuple(coordinates[index] for index in local_nodes)
        temps = tuple(temperatures[index] for index in local_nodes)
        context = "rank {} local cell {} global element {}".format(
            rank, local_cell, element_ids[cell_index]
        )
        gradient_t, centroid, det_j = tetra_gradient(coords, temps, context)

        result[local_cell] = {
            "global_element_id": element_ids[cell_index],
            "material_id": material_ids[cell_index],
            "gids": gids,
            "coordinates": coords,
            "temperatures": temps,
            "coordinate_by_gid": dict(zip(gids, coords)),
            "temperature_by_gid": dict(zip(gids, temps)),
            "gradient_t": gradient_t,
            "centroid": centroid,
            "det_j": det_j,
        }

    return result


def face_geometry_and_normal(nodes, fluid_tet, solid_tet):
    try:
        p0 = fluid_tet["coordinate_by_gid"][nodes[0]]
        p1 = fluid_tet["coordinate_by_gid"][nodes[1]]
        p2 = fluid_tet["coordinate_by_gid"][nodes[2]]
    except KeyError as error:
        raise RuntimeError(
            "Fluid tetrahedron does not contain interface node {}".format(error)
        )

    for gid in nodes:
        if gid not in solid_tet["coordinate_by_gid"]:
            raise RuntimeError(
                "Solid tetrahedron does not contain interface node {}".format(gid)
            )

    edge1 = vector_subtract(p1, p0)
    edge2 = vector_subtract(p2, p0)
    raw_normal = cross_product(edge1, edge2)
    raw_norm = vector_norm(raw_normal)
    if not math.isfinite(raw_norm) or raw_norm <= 0.0:
        raise RuntimeError("Degenerate interface face {}".format(nodes))

    area = 0.5 * raw_norm
    normal = vector_scale(raw_normal, 1.0 / raw_norm)
    centroid = (
        (p0[0] + p1[0] + p2[0]) / 3.0,
        (p0[1] + p1[1] + p2[1]) / 3.0,
        (p0[2] + p1[2] + p2[2]) / 3.0,
    )

    fluid_vector = vector_subtract(fluid_tet["centroid"], centroid)
    if dot_product(normal, fluid_vector) < 0.0:
        normal = vector_scale(normal, -1.0)

    fluid_distance = dot_product(
        vector_subtract(fluid_tet["centroid"], centroid), normal
    )
    solid_distance = dot_product(
        vector_subtract(solid_tet["centroid"], centroid), normal
    )

    if fluid_distance <= 0.0 or solid_distance >= 0.0:
        raise RuntimeError(
            "Cannot orient solid-to-fluid normal for face {}: "
            "fluid_distance={}, solid_distance={}".format(
                nodes, fluid_distance, solid_distance
            )
        )

    return area, centroid, normal, fluid_distance, solid_distance


def finite_min(current, value):
    return value if current is None else min(current, value)


def finite_max(current, value):
    return value if current is None else max(current, value)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pvtu", required=True, type=Path)
    parser.add_argument("--interface-csv", required=True, type=Path)
    parser.add_argument("--partition-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    pvtu = args.pvtu.resolve()
    interface_csv = args.interface_csv.resolve()
    partition_root = args.partition_root.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if not pvtu.is_file():
        raise RuntimeError("PVTU file does not exist: {}".format(pvtu))
    if not interface_csv.is_file():
        raise RuntimeError("Interface CSV does not exist: {}".format(interface_csv))
    if not partition_root.is_dir():
        raise RuntimeError(
            "Partition root does not exist: {}".format(partition_root)
        )

    pieces = find_piece_sources(pvtu)
    material_properties, matprop_path, matprop_count = read_material_properties(
        partition_root
    )
    rows, needed_by_rank, solid_material_ids = read_interface_rows(interface_csv)

    missing_ranks = sorted(set(needed_by_rank).difference(pieces))
    if missing_ranks:
        raise RuntimeError("Interface CSV references absent ranks {}".format(missing_ranks))

    for material_id in sorted(solid_material_ids):
        if material_id not in material_properties:
            raise RuntimeError(
                "Interface references missing solid material {}".format(material_id)
            )
        if material_properties[material_id]["phase"].lower() != "solid":
            raise RuntimeError(
                "Interface material {} is not labelled as solid".format(material_id)
            )
        if material_properties[material_id]["k"] <= 0.0:
            raise RuntimeError(
                "Solid material {} has non-positive conductivity".format(material_id)
            )

    stem = pvtu.stem
    output_csv = output_dir / (stem + "_interface_heat_flux.csv")
    summary_path = output_dir / (stem + "_interface_heat_flux_summary.txt")

    print("CBS3D one-sided P1 CHT interface heat-flux reconstruction")
    print("PVTU                 : {}".format(pvtu))
    print("Interface CSV        : {}".format(interface_csv))
    print("Interface faces      : {}".format(len(rows)))
    print("VTU pieces           : {}".format(len(pieces)))
    print("Material properties  : {}".format(matprop_path))
    print("Output CSV           : {}".format(output_csv))

    tetrahedra = {}
    total_requested = sum(len(cells) for cells in needed_by_rank.values())
    loaded_requested = 0
    for rank in sorted(needed_by_rank):
        local_data = read_requested_tetrahedra(
            pieces[rank], rank, needed_by_rank[rank]
        )
        for local_cell, data in local_data.items():
            tetrahedra[(rank, local_cell)] = data
        loaded_requested += len(local_data)
        print(
            "  loaded rank {:04d}: requested tetrahedra={}, cumulative={}/{}".format(
                rank, len(local_data), loaded_requested, total_requested
            )
        )

    fluid_k = material_properties[0]["k"]
    interface_area = 0.0
    q_fluid_integral = 0.0
    q_solid_integral = 0.0
    q_mean_integral = 0.0
    tw_integral = 0.0
    max_temperature_mismatch = 0.0
    max_area_relative_mismatch = 0.0
    maximum_relative_flux_jump = 0.0
    relative_flux_jump_integral = 0.0
    negative_q_fluid_area = 0.0
    negative_q_solid_area = 0.0
    local_faces = 0
    cross_rank_faces = 0
    min_tw = None
    max_tw = None
    min_q_fluid = None
    max_q_fluid = None
    min_q_solid = None
    max_q_solid = None

    with output_csv.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(OUTPUT_HEADER)

        for index, row in enumerate(rows, start=1):
            fluid_key = (row["fluid_rank"], row["fluid_local_cell"])
            solid_key = (row["solid_rank"], row["solid_local_cell"])
            fluid_tet = tetrahedra.get(fluid_key)
            solid_tet = tetrahedra.get(solid_key)
            if fluid_tet is None or solid_tet is None:
                raise RuntimeError(
                    "Missing requested tetrahedron at interface face {}".format(index)
                )

            if fluid_tet["global_element_id"] != row["fluid_element_id"]:
                raise RuntimeError(
                    "Fluid element ID mismatch at interface face {}".format(index)
                )
            if solid_tet["global_element_id"] != row["solid_element_id"]:
                raise RuntimeError(
                    "Solid element ID mismatch at interface face {}".format(index)
                )
            if fluid_tet["material_id"] != 0:
                raise RuntimeError(
                    "Fluid tetrahedron has material {} at interface face {}".format(
                        fluid_tet["material_id"], index
                    )
                )
            if solid_tet["material_id"] != row["solid_material_id"]:
                raise RuntimeError(
                    "Solid material mismatch at interface face {}".format(index)
                )

            area, centroid, normal, fluid_distance, solid_distance = (
                face_geometry_and_normal(row["nodes"], fluid_tet, solid_tet)
            )
            area_relative_mismatch = abs(area - row["area"]) / max(
                area, row["area"], 1.0e-300
            )
            centroid_scale = max(
                abs(centroid[0]), abs(centroid[1]), abs(centroid[2]), 1.0
            )
            centroid_error = vector_norm(vector_subtract(centroid, row["centroid"]))
            if centroid_error > 1.0e-11 * centroid_scale:
                raise RuntimeError(
                    "Interface centroid mismatch at face {}: {}".format(
                        index, centroid_error
                    )
                )

            fluid_face_temperatures = [
                fluid_tet["temperature_by_gid"][gid] for gid in row["nodes"]
            ]
            solid_face_temperatures = [
                solid_tet["temperature_by_gid"][gid] for gid in row["nodes"]
            ]
            tw_fluid = sum(fluid_face_temperatures) / 3.0
            tw_solid = sum(solid_face_temperatures) / 3.0
            temperature_mismatch = max(
                abs(a - b)
                for a, b in zip(fluid_face_temperatures, solid_face_temperatures)
            )
            tw = 0.5 * (tw_fluid + tw_solid)

            grad_fluid = fluid_tet["gradient_t"]
            grad_solid = solid_tet["gradient_t"]
            dtdn_fluid = dot_product(grad_fluid, normal)
            dtdn_solid = dot_product(grad_solid, normal)
            solid_k = material_properties[row["solid_material_id"]]["k"]
            q_fluid = -fluid_k * dtdn_fluid
            q_solid = -solid_k * dtdn_solid
            q_mean = 0.5 * (q_fluid + q_solid)
            q_jump = q_fluid - q_solid
            relative_jump = abs(q_jump) / max(
                abs(q_fluid), abs(q_solid), 1.0e-300
            )

            values_to_check = (
                area, centroid[0], centroid[1], centroid[2],
                normal[0], normal[1], normal[2], tw,
                grad_fluid[0], grad_fluid[1], grad_fluid[2],
                grad_solid[0], grad_solid[1], grad_solid[2],
                dtdn_fluid, dtdn_solid, q_fluid, q_solid,
                q_mean, q_jump, relative_jump,
            )
            if not all(math.isfinite(value) for value in values_to_check):
                raise RuntimeError(
                    "Non-finite reconstructed quantity at interface face {}".format(index)
                )

            interface_area += area
            q_fluid_integral += q_fluid * area
            q_solid_integral += q_solid * area
            q_mean_integral += q_mean * area
            tw_integral += tw * area
            relative_flux_jump_integral += relative_jump * area
            max_temperature_mismatch = max(
                max_temperature_mismatch, temperature_mismatch
            )
            max_area_relative_mismatch = max(
                max_area_relative_mismatch, area_relative_mismatch
            )
            maximum_relative_flux_jump = max(
                maximum_relative_flux_jump, relative_jump
            )
            if q_fluid < 0.0:
                negative_q_fluid_area += area
            if q_solid < 0.0:
                negative_q_solid_area += area

            min_tw = finite_min(min_tw, tw)
            max_tw = finite_max(max_tw, tw)
            min_q_fluid = finite_min(min_q_fluid, q_fluid)
            max_q_fluid = finite_max(max_q_fluid, q_fluid)
            min_q_solid = finite_min(min_q_solid, q_solid)
            max_q_solid = finite_max(max_q_solid, q_solid)

            if row["pair_location"] == "local":
                local_faces += 1
            elif row["pair_location"] == "cross_rank":
                cross_rank_faces += 1
            else:
                raise RuntimeError(
                    "Unknown pair_location '{}' at face {}".format(
                        row["pair_location"], index
                    )
                )

            writer.writerow((
                row["nodes"][0], row["nodes"][1], row["nodes"][2],
                row["fluid_element_id"], row["solid_element_id"],
                row["solid_material_id"], row["fluid_rank"],
                row["fluid_local_cell"], row["solid_rank"],
                row["solid_local_cell"], row["pair_location"],
                "{:.17g}".format(area),
                "{:.17g}".format(centroid[0]),
                "{:.17g}".format(centroid[1]),
                "{:.17g}".format(centroid[2]),
                "{:.17g}".format(normal[0]),
                "{:.17g}".format(normal[1]),
                "{:.17g}".format(normal[2]),
                "{:.17g}".format(tw),
                "{:.17g}".format(tw_fluid),
                "{:.17g}".format(tw_solid),
                "{:.17g}".format(temperature_mismatch),
                "{:.17g}".format(grad_fluid[0]),
                "{:.17g}".format(grad_fluid[1]),
                "{:.17g}".format(grad_fluid[2]),
                "{:.17g}".format(grad_solid[0]),
                "{:.17g}".format(grad_solid[1]),
                "{:.17g}".format(grad_solid[2]),
                "{:.17g}".format(dtdn_fluid),
                "{:.17g}".format(dtdn_solid),
                "{:.17g}".format(fluid_k),
                "{:.17g}".format(solid_k),
                "{:.17g}".format(q_fluid),
                "{:.17g}".format(q_solid),
                "{:.17g}".format(q_mean),
                "{:.17g}".format(q_jump),
                "{:.17g}".format(relative_jump),
                "{:.17g}".format(fluid_distance),
                "{:.17g}".format(solid_distance),
                "{:.17g}".format(area_relative_mismatch),
            ))

    area_weighted_tw = tw_integral / interface_area
    mean_q_fluid = q_fluid_integral / interface_area
    mean_q_solid = q_solid_integral / interface_area
    mean_q_mean = q_mean_integral / interface_area
    mean_relative_jump = relative_flux_jump_integral / interface_area
    integrated_flux_difference = q_fluid_integral - q_solid_integral
    integrated_relative_difference = abs(integrated_flux_difference) / max(
        abs(q_fluid_integral), abs(q_solid_integral), 1.0e-300
    )

    lines = []
    lines.append("CBS3D one-sided P1 CHT interface heat-flux summary")
    lines.append("PVTU                                      : {}".format(pvtu))
    lines.append("Interface CSV                             : {}".format(interface_csv))
    lines.append("Output CSV                                : {}".format(output_csv))
    lines.append("Material property file                    : {}".format(matprop_path))
    lines.append("Material property copies                  : {}".format(matprop_count))
    lines.append("")
    lines.append("Topology and geometry")
    lines.append("  interface faces                         : {}".format(len(rows)))
    lines.append("  local interface faces                   : {}".format(local_faces))
    lines.append("  cross-rank interface faces              : {}".format(cross_rank_faces))
    lines.append("  unique requested tetrahedra             : {}".format(len(tetrahedra)))
    lines.append("  total interface area [m^2]              : {:.17e}".format(interface_area))
    lines.append("  maximum area relative mismatch          : {:.17e}".format(max_area_relative_mismatch))
    lines.append("")
    lines.append("Material properties")
    lines.append("  fluid material                          : 0 {}".format(material_properties[0]["name"]))
    lines.append("  fluid conductivity [W/(m K)]            : {:.17e}".format(fluid_k))
    for material_id in sorted(solid_material_ids):
        props = material_properties[material_id]
        lines.append(
            "  solid material {}                        : {}".format(
                material_id, props["name"]
            )
        )
        lines.append(
            "  solid conductivity {} [W/(m K)]         : {:.17e}".format(
                material_id, props["k"]
            )
        )
    lines.append("")
    lines.append("Interface temperature")
    lines.append("  minimum wall temperature [K]            : {:.17e}".format(min_tw))
    lines.append("  maximum wall temperature [K]            : {:.17e}".format(max_tw))
    lines.append("  area-weighted wall temperature [K]      : {:.17e}".format(area_weighted_tw))
    lines.append("  maximum fluid/solid nodal mismatch [K]  : {:.17e}".format(max_temperature_mismatch))
    lines.append("")
    lines.append("One-sided raw P1 normal heat flux")
    lines.append("  normal convention                        : solid to fluid")
    lines.append("  positive flux convention                 : heat flowing solid to fluid")
    lines.append("  minimum fluid-side flux [W/m^2]         : {:.17e}".format(min_q_fluid))
    lines.append("  maximum fluid-side flux [W/m^2]         : {:.17e}".format(max_q_fluid))
    lines.append("  area-mean fluid-side flux [W/m^2]       : {:.17e}".format(mean_q_fluid))
    lines.append("  integrated fluid-side heat rate [W]     : {:.17e}".format(q_fluid_integral))
    lines.append("  minimum solid-side flux [W/m^2]         : {:.17e}".format(min_q_solid))
    lines.append("  maximum solid-side flux [W/m^2]         : {:.17e}".format(max_q_solid))
    lines.append("  area-mean solid-side flux [W/m^2]       : {:.17e}".format(mean_q_solid))
    lines.append("  integrated solid-side heat rate [W]     : {:.17e}".format(q_solid_integral))
    lines.append("  area-mean arithmetic flux [W/m^2]       : {:.17e}".format(mean_q_mean))
    lines.append("  integrated arithmetic heat rate [W]     : {:.17e}".format(q_mean_integral))
    lines.append("  integrated fluid-solid difference [W]   : {:.17e}".format(integrated_flux_difference))
    lines.append("  integrated relative difference          : {:.17e}".format(integrated_relative_difference))
    lines.append("  area-mean facewise relative jump         : {:.17e}".format(mean_relative_jump))
    lines.append("  maximum facewise relative jump           : {:.17e}".format(maximum_relative_flux_jump))
    lines.append("  negative fluid-side flux area fraction  : {:.17e}".format(negative_q_fluid_area / interface_area))
    lines.append("  negative solid-side flux area fraction  : {:.17e}".format(negative_q_solid_area / interface_area))
    lines.append("")

    hard_errors = []
    if len(rows) != local_faces + cross_rank_faces:
        hard_errors.append("interface face location count is inconsistent")
    if max_temperature_mismatch > 1.0e-9:
        hard_errors.append(
            "fluid/solid shared-node temperatures differ by more than 1e-9 K"
        )
    if max_area_relative_mismatch > 1.0e-10:
        hard_errors.append(
            "reconstructed face area differs from extraction by more than 1e-10 relative"
        )

    lines.append("Validation result")
    if hard_errors:
        for error in hard_errors:
            lines.append("  FAIL: {}".format(error))
    else:
        lines.append("  PASS: interface topology, geometry, temperature continuity, and P1 gradients reconstructed")
        lines.append("  NOTE: one-sided flux jump is diagnostic; continuous P1 elements do not impose pointwise flux continuity")
        lines.append("  NEXT: construct a conservative recovered wall flux before calculating h and Nu")

    summary = "\n".join(lines) + "\n"
    summary_path.write_text(summary, encoding="utf-8")
    print("\n" + summary)
    print("Output CSV      : {}".format(output_csv))
    print("Summary         : {}".format(summary_path))

    return 1 if hard_errors else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        print("ERROR: {}".format(error), file=sys.stderr)
        sys.exit(2)
