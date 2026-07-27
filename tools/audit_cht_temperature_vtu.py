#!/usr/bin/env python3
"""Audit a CBS3D CHT temperature field directly from serial or rank-local VTU files.

The audit is intentionally read-only. It checks the three mechanisms most likely
responsible for a mesh-imprinted heat-flux temperature field:

1. BC 532 coverage and duplicate heat-flux boundary triangles.
2. Variation of the lumped nodal heat-load / thermal-capacitance ratio.
3. Agreement of duplicated MPI point values and correlation of surface
   temperature with the load/capacitance ratio.

Only Python's standard library is required. CBS3D VTU output is ASCII XML, so no
VTK Python installation is needed.
"""

from __future__ import annotations

import argparse
import math
import statistics
import sys
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import DefaultDict, Dict, Iterable, List, Mapping, Sequence, Tuple

Vec3 = Tuple[float, float, float]
CoordKey = Tuple[float, float, float]


@dataclass(frozen=True)
class Material:
    rho: float
    cp: float
    k: float
    mu: float
    qvol: float

    @property
    def rho_cp(self) -> float:
        return self.rho * self.cp

    @property
    def alpha(self) -> float:
        return self.k / self.rho_cp


@dataclass
class Piece:
    path: Path
    points: List[Vec3]
    connectivity: List[int]
    offsets: List[int]
    types: List[int]
    temperature: List[float]
    mat_id: List[int]
    bc_id: List[int]


def fail(message: str) -> "NoReturn":
    raise RuntimeError(message)


def parse_float_array(node: ET.Element, context: str) -> List[float]:
    text = node.text or ""
    try:
        return [float(value) for value in text.split()]
    except ValueError as exc:
        raise RuntimeError(f"{context}: invalid floating-point DataArray") from exc


def parse_int_array(node: ET.Element, context: str) -> List[int]:
    text = node.text or ""
    try:
        return [int(value) for value in text.split()]
    except ValueError as exc:
        raise RuntimeError(f"{context}: invalid integer DataArray") from exc


def named_data_array(parent: ET.Element, name: str, context: str) -> ET.Element:
    for node in parent.findall("DataArray"):
        if node.get("Name") == name:
            return node
    fail(f"{context}: missing DataArray Name={name!r}")


def read_piece(path: Path) -> Piece:
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError as exc:
        raise RuntimeError(f"{path}: malformed VTU XML") from exc

    piece_node = root.find("./UnstructuredGrid/Piece")
    if piece_node is None:
        fail(f"{path}: missing UnstructuredGrid/Piece")

    points_parent = piece_node.find("Points")
    cells_parent = piece_node.find("Cells")
    point_data = piece_node.find("PointData")
    cell_data = piece_node.find("CellData")
    if None in (points_parent, cells_parent, point_data, cell_data):
        fail(f"{path}: incomplete VTU Piece")

    assert points_parent is not None
    assert cells_parent is not None
    assert point_data is not None
    assert cell_data is not None

    point_array = points_parent.find("DataArray")
    if point_array is None:
        fail(f"{path}: missing point coordinates")
    flat_points = parse_float_array(point_array, str(path))
    if len(flat_points) % 3 != 0:
        fail(f"{path}: point coordinate count is not divisible by three")
    points = [
        (flat_points[i], flat_points[i + 1], flat_points[i + 2])
        for i in range(0, len(flat_points), 3)
    ]

    connectivity = parse_int_array(
        named_data_array(cells_parent, "connectivity", str(path)), str(path)
    )
    offsets = parse_int_array(
        named_data_array(cells_parent, "offsets", str(path)), str(path)
    )
    types = parse_int_array(
        named_data_array(cells_parent, "types", str(path)), str(path)
    )
    temperature = parse_float_array(
        named_data_array(point_data, "temperature", str(path)), str(path)
    )
    mat_id = parse_int_array(
        named_data_array(cell_data, "mat_id", str(path)), str(path)
    )
    bc_id = parse_int_array(
        named_data_array(cell_data, "bc_id", str(path)), str(path)
    )

    npoints_declared = int(piece_node.get("NumberOfPoints", "-1"))
    ncells_declared = int(piece_node.get("NumberOfCells", "-1"))
    if npoints_declared != len(points):
        fail(f"{path}: NumberOfPoints mismatch")
    if ncells_declared != len(offsets):
        fail(f"{path}: NumberOfCells mismatch")
    if len(types) != len(offsets) or len(mat_id) != len(offsets) or len(bc_id) != len(offsets):
        fail(f"{path}: cell-data array length mismatch")
    if len(temperature) != len(points):
        fail(f"{path}: temperature array length mismatch")
    if offsets and offsets[-1] != len(connectivity):
        fail(f"{path}: final connectivity offset mismatch")

    return Piece(
        path=path,
        points=points,
        connectivity=connectivity,
        offsets=offsets,
        types=types,
        temperature=temperature,
        mat_id=mat_id,
        bc_id=bc_id,
    )


def iter_cells(piece: Piece) -> Iterable[Tuple[int, int, List[int], int, int]]:
    start = 0
    for cell_index, end in enumerate(piece.offsets):
        if end < start or end > len(piece.connectivity):
            fail(f"{piece.path}: invalid cell offset at cell {cell_index}")
        nodes = piece.connectivity[start:end]
        yield cell_index, piece.types[cell_index], nodes, piece.mat_id[cell_index], piece.bc_id[cell_index]
        start = end


def sub(a: Vec3, b: Vec3) -> Vec3:
    return a[0] - b[0], a[1] - b[1], a[2] - b[2]


def cross(a: Vec3, b: Vec3) -> Vec3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def dot(a: Vec3, b: Vec3) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def norm(a: Vec3) -> float:
    return math.sqrt(dot(a, a))


def triangle_area(a: Vec3, b: Vec3, c: Vec3) -> float:
    return 0.5 * norm(cross(sub(b, a), sub(c, a)))


def tetra_volume(a: Vec3, b: Vec3, c: Vec3, d: Vec3) -> float:
    return abs(dot(sub(b, a), cross(sub(c, a), sub(d, a)))) / 6.0


def coord_key(point: Vec3, digits: int) -> CoordKey:
    return round(point[0], digits), round(point[1], digits), round(point[2], digits)


def read_matprop(path: Path) -> Dict[int, Material]:
    tokens = path.read_text(encoding="utf-8").split()
    if not tokens:
        fail(f"{path}: empty material-property file")
    try:
        nmat = int(tokens[0])
    except ValueError as exc:
        raise RuntimeError(f"{path}: invalid material count") from exc

    expected = 1 + 8 * nmat
    if len(tokens) < expected:
        fail(f"{path}: expected {nmat} material records")

    materials: Dict[int, Material] = {}
    pos = 1
    for _ in range(nmat):
        try:
            material_id = int(tokens[pos])
            # name and phase are retained only for human-readable input.
            rho = float(tokens[pos + 3])
            cp = float(tokens[pos + 4])
            k = float(tokens[pos + 5])
            mu = float(tokens[pos + 6])
            qvol = float(tokens[pos + 7])
        except (ValueError, IndexError) as exc:
            raise RuntimeError(f"{path}: invalid material-property record") from exc
        pos += 8
        if material_id in materials:
            fail(f"{path}: duplicate material ID {material_id}")
        if rho <= 0.0 or cp <= 0.0 or k <= 0.0 or mu < 0.0:
            fail(f"{path}: non-physical material properties for ID {material_id}")
        materials[material_id] = Material(rho, cp, k, mu, qvol)
    return materials


def descriptive(values: Sequence[float]) -> Mapping[str, float]:
    if not values:
        return {}
    ordered = sorted(values)
    n = len(ordered)

    def percentile(p: float) -> float:
        if n == 1:
            return ordered[0]
        x = p * (n - 1)
        lo = int(math.floor(x))
        hi = int(math.ceil(x))
        if lo == hi:
            return ordered[lo]
        w = x - lo
        return ordered[lo] * (1.0 - w) + ordered[hi] * w

    mean = statistics.fmean(ordered)
    std = statistics.pstdev(ordered) if n > 1 else 0.0
    return {
        "min": ordered[0],
        "p05": percentile(0.05),
        "median": percentile(0.50),
        "p95": percentile(0.95),
        "max": ordered[-1],
        "mean": mean,
        "std": std,
        "cv": std / abs(mean) if mean != 0.0 else math.inf,
    }


def pearson(x: Sequence[float], y: Sequence[float]) -> float:
    if len(x) != len(y) or len(x) < 2:
        return math.nan
    mx = statistics.fmean(x)
    my = statistics.fmean(y)
    dx = [value - mx for value in x]
    dy = [value - my for value in y]
    denom = math.sqrt(sum(value * value for value in dx) * sum(value * value for value in dy))
    if denom == 0.0:
        return math.nan
    return sum(a * b for a, b in zip(dx, dy)) / denom


def print_stats(label: str, stats: Mapping[str, float], unit: str = "") -> None:
    if not stats:
        print(f"{label}: no values")
        return
    suffix = f" {unit}" if unit else ""
    print(f"{label}:")
    for name in ("min", "p05", "median", "p95", "max", "mean", "std", "cv"):
        value = stats[name]
        if name == "cv":
            print(f"  {name:>8s} = {value:.9e}")
        else:
            print(f"  {name:>8s} = {value:.9e}{suffix}")


def discover_vtu_files(directory: Path, step: int) -> List[Path]:
    if not directory.is_dir():
        fail(f"VTU directory does not exist: {directory}")
    tag = f"{step:08d}"
    files = sorted(path for path in directory.glob("*.vtu") if tag in path.name)
    if not files:
        # Accept a serial writer name that may not use the zero-padded step tag.
        files = sorted(path for path in directory.glob("*.vtu") if str(step) in path.stem)
    if not files:
        fail(f"No VTU files for step {step} were found in {directory}")
    return files


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vtu-dir", type=Path, required=True)
    parser.add_argument("--step", type=int, required=True)
    parser.add_argument("--matprop", type=Path, required=True)
    parser.add_argument("--heat-flux", type=float, required=True)
    parser.add_argument("--heat-flux-bc", type=int, default=532)
    parser.add_argument("--reference-temperature", type=float, default=300.0)
    parser.add_argument("--dt", type=float, default=1.0e-5)
    parser.add_argument("--completed-iterations", type=int, default=None)
    parser.add_argument("--coord-digits", type=int, default=12)
    args = parser.parse_args()

    if args.step < 0 or args.heat_flux <= 0.0 or args.dt <= 0.0:
        fail("step, heat flux and dt must be positive")
    completed_iterations = args.completed_iterations or args.step

    materials = read_matprop(args.matprop)
    files = discover_vtu_files(args.vtu_dir, args.step)
    print(f"VTU pieces                     : {len(files)}")
    print(f"Material-property file         : {args.matprop}")

    capacity: DefaultDict[CoordKey, float] = defaultdict(float)
    heat_load: DefaultDict[CoordKey, float] = defaultdict(float)
    temperatures: DefaultDict[CoordKey, List[float]] = defaultdict(list)
    bc_face_count: Counter[int] = Counter()
    bc_face_area: DefaultDict[int, float] = defaultdict(float)
    boundary_faces: List[Tuple[int, float, Tuple[CoordKey, CoordKey, CoordKey], Vec3]] = []
    heat_face_points: List[Vec3] = []
    face_occurrences: Counter[Tuple[CoordKey, CoordKey, CoordKey]] = Counter()
    material_tet_count: Counter[int] = Counter()
    material_volume: DefaultDict[int, float] = defaultdict(float)
    total_thermal_energy = 0.0

    global_min = [math.inf, math.inf, math.inf]
    global_max = [-math.inf, -math.inf, -math.inf]

    for path in files:
        piece = read_piece(path)
        keys = [coord_key(point, args.coord_digits) for point in piece.points]
        for point, key, temperature in zip(piece.points, keys, piece.temperature):
            temperatures[key].append(temperature)
            for axis in range(3):
                global_min[axis] = min(global_min[axis], point[axis])
                global_max[axis] = max(global_max[axis], point[axis])

        for _, vtk_type, nodes, material_id, bc_id in iter_cells(piece):
            if vtk_type == 10:
                if len(nodes) != 4:
                    fail(f"{path}: tetrahedron does not contain four nodes")
                if material_id not in materials:
                    fail(f"{path}: material ID {material_id} missing from {args.matprop}")
                points = [piece.points[index] for index in nodes]
                volume = tetra_volume(*points)
                if volume <= 0.0 or not math.isfinite(volume):
                    fail(f"{path}: invalid tetrahedral volume")
                material = materials[material_id]
                nodal_capacity = material.rho_cp * volume / 4.0
                material_tet_count[material_id] += 1
                material_volume[material_id] += volume
                element_temperature_sum = 0.0
                for index in nodes:
                    key = keys[index]
                    capacity[key] += nodal_capacity
                    element_temperature_sum += piece.temperature[index]
                total_thermal_energy += (
                    material.rho_cp
                    * volume
                    * (element_temperature_sum / 4.0 - args.reference_temperature)
                )
            elif vtk_type == 5:
                if len(nodes) != 3:
                    fail(f"{path}: boundary triangle does not contain three nodes")
                points = [piece.points[index] for index in nodes]
                area = triangle_area(*points)
                if area <= 0.0 or not math.isfinite(area):
                    fail(f"{path}: invalid boundary-triangle area")
                face_keys = tuple(sorted(keys[index] for index in nodes))
                centroid = (
                    sum(point[0] for point in points) / 3.0,
                    sum(point[1] for point in points) / 3.0,
                    sum(point[2] for point in points) / 3.0,
                )
                face_occurrences[face_keys] += 1
                bc_face_count[bc_id] += 1
                bc_face_area[bc_id] += area
                boundary_faces.append((bc_id, area, face_keys, centroid))
                if bc_id == args.heat_flux_bc:
                    contribution = args.heat_flux * area / 3.0
                    for index in nodes:
                        key = keys[index]
                        heat_load[key] += contribution
                        heat_face_points.append(piece.points[index])

    if not heat_face_points:
        fail(f"No boundary triangles with BC {args.heat_flux_bc} were found")

    duplicate_point_mismatch = 0.0
    duplicate_point_count = 0
    representative_temperature: Dict[CoordKey, float] = {}
    for key, values in temperatures.items():
        representative_temperature[key] = statistics.fmean(values)
        if len(values) > 1:
            duplicate_point_count += 1
            duplicate_point_mismatch = max(duplicate_point_mismatch, max(values) - min(values))

    duplicate_faces = sum(count - 1 for count in face_occurrences.values() if count > 1)
    heat_area = bc_face_area[args.heat_flux_bc]
    total_power = args.heat_flux * heat_area

    heat_ranges = [
        max(point[axis] for point in heat_face_points) - min(point[axis] for point in heat_face_points)
        for axis in range(3)
    ]
    plane_axis = min(range(3), key=lambda axis: heat_ranges[axis])
    plane_coordinate = statistics.median(point[plane_axis] for point in heat_face_points)
    domain_span = max(global_max[axis] - global_min[axis] for axis in range(3))
    plane_tolerance = max(1.0e-12, domain_span * 1.0e-9)
    plane_counts: Counter[int] = Counter()
    plane_areas: DefaultDict[int, float] = defaultdict(float)
    for bc_id, area, _, centroid in boundary_faces:
        if abs(centroid[plane_axis] - plane_coordinate) <= plane_tolerance:
            plane_counts[bc_id] += 1
            plane_areas[bc_id] += area

    rates: List[float] = []
    heat_temperatures: List[float] = []
    missing_capacity = 0
    for key, load in heat_load.items():
        cap = capacity.get(key, 0.0)
        if cap <= 0.0:
            missing_capacity += 1
            continue
        rates.append(load / cap)
        heat_temperatures.append(representative_temperature[key])

    nominal_time = completed_iterations * args.dt
    nominal_input_energy = total_power * nominal_time

    print("\n=== Boundary-condition audit ===")
    for bc_id in sorted(bc_face_count):
        print(
            f"BC {bc_id:4d}: faces={bc_face_count[bc_id]:9d} "
            f"area={bc_face_area[bc_id]:.12e} m^2"
        )
    print(f"Duplicate physical boundary faces: {duplicate_faces}")
    print(
        f"Heat-flux plane               : axis={'xyz'[plane_axis]}, "
        f"coordinate={plane_coordinate:.12e}, tolerance={plane_tolerance:.3e}"
    )
    for bc_id in sorted(plane_counts):
        print(
            f"  plane BC {bc_id:4d}: faces={plane_counts[bc_id]:9d} "
            f"area={plane_areas[bc_id]:.12e} m^2"
        )

    print("\n=== Material audit ===")
    for material_id in sorted(material_tet_count):
        material = materials[material_id]
        print(
            f"Material {material_id:4d}: tets={material_tet_count[material_id]:9d} "
            f"volume={material_volume[material_id]:.12e} m^3 "
            f"rhoCp={material.rho_cp:.12e} J/(m^3 K) "
            f"k={material.k:.12e} W/(m K) alpha={material.alpha:.12e} m^2/s"
        )

    print("\n=== Heat-load / capacity audit ===")
    print(f"BC {args.heat_flux_bc} area              : {heat_area:.12e} m^2")
    print(f"Applied heat flux              : {args.heat_flux:.12e} W/m^2")
    print(f"Applied thermal power          : {total_power:.12e} W")
    print(f"Heated unique nodes            : {len(heat_load)}")
    print(f"Heated nodes missing capacity  : {missing_capacity}")
    print_stats("Nodal initial dT/dt = F_i/C_i", descriptive(rates), "K/s")
    print_stats("Temperature on BC 532 nodes", descriptive(heat_temperatures), "K")
    print(
        "Pearson correlation T versus F_i/C_i: "
        f"{pearson(rates, heat_temperatures):.9e}"
    )

    all_temperatures = list(representative_temperature.values())
    print("\n=== Field and restart-safety audit ===")
    print_stats("Global nodal temperature", descriptive(all_temperatures), "K")
    print(f"Unique coordinate points       : {len(representative_temperature)}")
    print(f"Coordinates duplicated in VTUs : {duplicate_point_count}")
    print(f"Maximum duplicate dT           : {duplicate_point_mismatch:.12e} K")
    print(f"Nominal accumulated pseudo-time: {nominal_time:.12e} s")
    print(f"Nominal heat input qAt         : {nominal_input_energy:.12e} J")
    print(
        f"Stored sensible-energy change  : {total_thermal_energy:.12e} J "
        f"relative to {args.reference_temperature:.6f} K"
    )

    print("\n=== Interpretation flags ===")
    mixed_plane_bcs = [bc for bc in plane_counts if bc != args.heat_flux_bc]
    if mixed_plane_bcs:
        print(
            "FAIL: the geometric heat-flux plane contains non-532 boundary "
            f"IDs: {sorted(mixed_plane_bcs)}"
        )
    else:
        print("PASS: the detected heat-flux plane contains only BC 532 triangles")

    if duplicate_faces:
        print("FAIL: duplicate physical boundary triangles were detected")
    else:
        print("PASS: no duplicate physical boundary triangles were detected")

    rate_stats = descriptive(rates)
    if rate_stats and rate_stats["cv"] > 0.10:
        print(
            "WARNING: F_i/C_i varies by more than 10%; an explicit lumped "
            "thermal step will imprint the surface mesh at short pseudo-time"
        )
    else:
        print("PASS: the BC-load/capacity ratio is comparatively uniform")

    if duplicate_point_mismatch > 1.0e-10:
        print("WARNING: duplicated rank-point temperatures are not identical")
    else:
        print("PASS: duplicated rank-point temperatures agree")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
