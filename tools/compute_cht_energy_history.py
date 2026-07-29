#!/usr/bin/env python3
"""Compute CHT energy and boundary-power histories from distributed CBS3D VTU.

The tool processes every saved PVTU state one piece at a time. For each state it
computes exact P1 tetrahedral thermal energy, inlet/outlet mass and enthalpy
rates, imposed BC-532 heat input, and interval energy-balance diagnostics.

No NumPy, matplotlib, or modern type-hint syntax is required; the script is
compatible with the older Python available on Swansea Sunbird. It writes CSV,
plain-text summaries, and dependency-free SVG plots.
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


def determinant3(a11, a12, a13, a21, a22, a23, a31, a32, a33):
    return (
        a11 * (a22 * a33 - a23 * a32)
        - a12 * (a21 * a33 - a23 * a31)
        + a13 * (a21 * a32 - a22 * a31)
    )


def tetra_volume(points):
    p0, p1, p2, p3 = points
    a = vector_subtract(p1, p0)
    b = vector_subtract(p2, p0)
    c = vector_subtract(p3, p0)
    det = determinant3(
        a[0], b[0], c[0],
        a[1], b[1], c[1],
        a[2], b[2], c[2],
    )
    volume = abs(det) / 6.0
    if not math.isfinite(volume) or volume <= 0.0:
        raise RuntimeError("Encountered non-positive tetrahedral volume")
    return volume


def centroid(points):
    count = float(len(points))
    return (
        sum(point[0] for point in points) / count,
        sum(point[1] for point in points) / count,
        sum(point[2] for point in points) / count,
    )


def find_piece_sources(pvtu):
    root = ET.parse(str(pvtu)).getroot()
    grid = root.find("PUnstructuredGrid")
    if grid is None:
        raise RuntimeError("{} is not a PUnstructuredGrid file".format(pvtu))
    paths = []
    for piece in grid.findall("Piece"):
        source = piece.attrib.get("Source")
        if not source:
            raise RuntimeError("PVTU Piece without Source in {}".format(pvtu))
        path = (pvtu.parent / source).resolve()
        if not path.is_file():
            raise RuntimeError("Missing VTU piece: {}".format(path))
        paths.append(path)
    if not paths:
        raise RuntimeError("No Piece entries found in {}".format(pvtu))
    return paths


def read_material_properties(partition_root):
    files = sorted(partition_root.rglob("*.matprop"))
    if not files:
        raise RuntimeError("No .matprop files found below {}".format(partition_root))

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
    expected = int(lines[0].split()[0])
    properties = {}
    for line in lines[1:]:
        fields = line.split()
        if len(fields) < 8:
            raise RuntimeError("Malformed material-property row: {}".format(line))
        material = int(fields[0])
        properties[material] = {
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
            "Material file declares {} entries but {} were parsed".format(
                expected, len(properties)
            )
        )
    if 0 not in properties:
        raise RuntimeError("Fluid material 0 is absent")
    return properties, representative, len(files)


def iteration_from_name(path):
    match = re.search(r"step_(\d+)\.pvtu$", path.name)
    if not match:
        raise RuntimeError("Cannot parse iteration from {}".format(path.name))
    return int(match.group(1))


def interpolate_scalar(values, weights):
    return sum(value * weight for value, weight in zip(values, weights))


def interpolate_vector(values, weights):
    return (
        sum(value[0] * weight for value, weight in zip(values, weights)),
        sum(value[1] * weight for value, weight in zip(values, weights)),
        sum(value[2] * weight for value, weight in zip(values, weights)),
    )


def integrate_boundary_face(face_points, face_temperatures, face_velocities,
                            parent_points, rho, cp):
    edge1 = vector_subtract(face_points[1], face_points[0])
    edge2 = vector_subtract(face_points[2], face_points[0])
    area_vector = vector_scale(cross_product(edge1, edge2), 0.5)
    face_center = centroid(face_points)
    parent_center = centroid(parent_points)

    if dot_product(area_vector, vector_subtract(parent_center, face_center)) > 0.0:
        area_vector = vector_scale(area_vector, -1.0)

    area = vector_norm(area_vector)
    if not math.isfinite(area) or area <= 0.0:
        raise RuntimeError("Encountered non-positive boundary-face area")
    normal = vector_scale(area_vector, 1.0 / area)

    quadrature = (
        (2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0),
        (1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0),
        (1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0),
    )

    mass_rate = 0.0
    enthalpy_rate = 0.0
    for weights in quadrature:
        velocity = interpolate_vector(face_velocities, weights)
        temperature = interpolate_scalar(face_temperatures, weights)
        normal_velocity = dot_product(velocity, normal)
        weight = area / 3.0
        mass_rate += rho * normal_velocity * weight
        enthalpy_rate += rho * cp * normal_velocity * temperature * weight

    return area, mass_rate, enthalpy_rate


def process_piece(path, properties, heat_bc, inlet_bc, outlet_bc):
    root = ET.parse(str(path)).getroot()
    piece = root.find("./UnstructuredGrid/Piece")
    if piece is None:
        raise RuntimeError("No UnstructuredGrid/Piece in {}".format(path))

    number_of_points = int(piece.attrib["NumberOfPoints"])
    number_of_cells = int(piece.attrib["NumberOfCells"])

    point_values = parse_numbers(piece.find("./Points/DataArray").text, float)
    if len(point_values) != 3 * number_of_points:
        raise RuntimeError("{}: invalid coordinate count".format(path.name))
    coordinates = [
        (point_values[3 * i], point_values[3 * i + 1], point_values[3 * i + 2])
        for i in range(number_of_points)
    ]

    point_arrays = named_arrays(piece.find("PointData"))
    cell_arrays = named_arrays(piece.find("CellData"))
    cells_arrays = named_arrays(piece.find("Cells"))

    required_point = {"temperature", "velocity"}
    required_cell = {
        "cell_kind", "global_element_id", "material_id", "bc_id",
        "parent_global_element",
    }
    required_cells = {"connectivity", "offsets", "types"}
    missing = (
        required_point.difference(point_arrays),
        required_cell.difference(cell_arrays),
        required_cells.difference(cells_arrays),
    )
    if any(missing):
        raise RuntimeError(
            "{}: missing arrays PointData={}, CellData={}, Cells={}".format(
                path.name, sorted(missing[0]), sorted(missing[1]), sorted(missing[2])
            )
        )

    temperatures = parse_numbers(point_arrays["temperature"].text, float)
    velocity_values = parse_numbers(point_arrays["velocity"].text, float)
    if len(temperatures) != number_of_points or len(velocity_values) != 3 * number_of_points:
        raise RuntimeError("{}: inconsistent point-data lengths".format(path.name))
    velocities = [
        (velocity_values[3 * i], velocity_values[3 * i + 1], velocity_values[3 * i + 2])
        for i in range(number_of_points)
    ]

    connectivity = parse_numbers(cells_arrays["connectivity"].text, int)
    offsets = parse_numbers(cells_arrays["offsets"].text, int)
    types = parse_numbers(cells_arrays["types"].text, int)
    cell_kind = parse_numbers(cell_arrays["cell_kind"].text, int)
    global_element = parse_numbers(cell_arrays["global_element_id"].text, int)
    material_id = parse_numbers(cell_arrays["material_id"].text, int)
    bc_id = parse_numbers(cell_arrays["bc_id"].text, int)
    parent_global = parse_numbers(cell_arrays["parent_global_element"].text, int)

    for values in (offsets, types, cell_kind, global_element, material_id, bc_id, parent_global):
        if len(values) != number_of_cells:
            raise RuntimeError("{}: inconsistent cell-data lengths".format(path.name))

    cells = []
    start = 0
    for end in offsets:
        cells.append(connectivity[start:end])
        start = end
    if start != len(connectivity):
        raise RuntimeError("{}: connectivity offset mismatch".format(path.name))

    result = {
        "energy_by_material_J": collections.defaultdict(float),
        "volume_by_material_m3": collections.defaultdict(float),
        "temperature_volume_by_material_K_m3": collections.defaultdict(float),
        "tetrahedra_by_material": collections.Counter(),
        "boundary_faces": collections.Counter(),
        "boundary_area_m2": collections.defaultdict(float),
        "mass_rate_kg_s": collections.defaultdict(float),
        "enthalpy_rate_W": collections.defaultdict(float),
    }

    parent_by_global = {}
    for index in range(number_of_cells):
        if cell_kind[index] != 1:
            continue
        local_nodes = cells[index]
        if types[index] != 10 or len(local_nodes) != 4:
            raise RuntimeError("{}: invalid tetrahedron cell {}".format(path.name, index + 1))
        material = material_id[index]
        if material not in properties:
            raise RuntimeError("Unknown material {} in {}".format(material, path.name))
        points = [coordinates[node] for node in local_nodes]
        volume = tetra_volume(points)
        mean_temperature = sum(temperatures[node] for node in local_nodes) / 4.0
        rho_cp = properties[material]["rho"] * properties[material]["cp"]
        result["energy_by_material_J"][material] += rho_cp * volume * mean_temperature
        result["volume_by_material_m3"][material] += volume
        result["temperature_volume_by_material_K_m3"][material] += volume * mean_temperature
        result["tetrahedra_by_material"][material] += 1
        parent_by_global[global_element[index]] = local_nodes

    fluid = properties[0]
    for index in range(number_of_cells):
        if cell_kind[index] != 2:
            continue
        local_nodes = cells[index]
        if types[index] != 5 or len(local_nodes) != 3:
            raise RuntimeError("{}: invalid boundary triangle {}".format(path.name, index + 1))
        bc = bc_id[index]
        result["boundary_faces"][bc] += 1
        parent_nodes = parent_by_global.get(parent_global[index])
        if parent_nodes is None:
            raise RuntimeError(
                "{}: parent element {} unavailable for boundary cell {}".format(
                    path.name, parent_global[index], index + 1
                )
            )
        face_points = [coordinates[node] for node in local_nodes]
        face_temperatures = [temperatures[node] for node in local_nodes]
        face_velocities = [velocities[node] for node in local_nodes]
        parent_points = [coordinates[node] for node in parent_nodes]
        area, mass_rate, enthalpy_rate = integrate_boundary_face(
            face_points,
            face_temperatures,
            face_velocities,
            parent_points,
            fluid["rho"],
            fluid["cp"],
        )
        result["boundary_area_m2"][bc] += area
        if bc == inlet_bc or bc == outlet_bc:
            result["mass_rate_kg_s"][bc] += mass_rate
            result["enthalpy_rate_W"][bc] += enthalpy_rate

    return result


def process_state(pvtu, properties, heat_bc, inlet_bc, outlet_bc, heat_flux):
    total = {
        "energy_by_material_J": collections.defaultdict(float),
        "volume_by_material_m3": collections.defaultdict(float),
        "temperature_volume_by_material_K_m3": collections.defaultdict(float),
        "tetrahedra_by_material": collections.Counter(),
        "boundary_faces": collections.Counter(),
        "boundary_area_m2": collections.defaultdict(float),
        "mass_rate_kg_s": collections.defaultdict(float),
        "enthalpy_rate_W": collections.defaultdict(float),
    }

    pieces = find_piece_sources(pvtu)
    for index, path in enumerate(pieces, start=1):
        local = process_piece(path, properties, heat_bc, inlet_bc, outlet_bc)
        for field in (
            "energy_by_material_J", "volume_by_material_m3",
            "temperature_volume_by_material_K_m3", "boundary_area_m2",
            "mass_rate_kg_s", "enthalpy_rate_W",
        ):
            for key, value in local[field].items():
                total[field][key] += value
        total["tetrahedra_by_material"].update(local["tetrahedra_by_material"])
        total["boundary_faces"].update(local["boundary_faces"])
        print("    state {} piece {:02d}/{:02d}".format(pvtu.name, index, len(pieces)))

    fluid_energy = total["energy_by_material_J"].get(0, 0.0)
    solid_energy = sum(
        value for material, value in total["energy_by_material_J"].items()
        if material != 0
    )
    fluid_volume = total["volume_by_material_m3"].get(0, 0.0)
    solid_volume = sum(
        value for material, value in total["volume_by_material_m3"].items()
        if material != 0
    )
    fluid_temperature_mean = (
        total["temperature_volume_by_material_K_m3"].get(0, 0.0) / fluid_volume
        if fluid_volume > 0.0 else float("nan")
    )
    solid_temperature_mean = (
        sum(
            value for material, value in total["temperature_volume_by_material_K_m3"].items()
            if material != 0
        ) / solid_volume
        if solid_volume > 0.0 else float("nan")
    )

    inlet_mass = total["mass_rate_kg_s"].get(inlet_bc, 0.0)
    outlet_mass = total["mass_rate_kg_s"].get(outlet_bc, 0.0)
    inlet_enthalpy = total["enthalpy_rate_W"].get(inlet_bc, 0.0)
    outlet_enthalpy = total["enthalpy_rate_W"].get(outlet_bc, 0.0)
    net_outward_enthalpy = inlet_enthalpy + outlet_enthalpy
    heat_area = total["boundary_area_m2"].get(heat_bc, 0.0)
    applied_power = heat_flux * heat_area

    fluid = properties[0]
    inlet_bulk_temperature = (
        inlet_enthalpy / (fluid["cp"] * inlet_mass)
        if abs(inlet_mass) > 1.0e-30 else float("nan")
    )
    outlet_bulk_temperature = (
        outlet_enthalpy / (fluid["cp"] * outlet_mass)
        if abs(outlet_mass) > 1.0e-30 else float("nan")
    )

    return {
        "fluid_energy_J": fluid_energy,
        "solid_energy_J": solid_energy,
        "total_energy_J": fluid_energy + solid_energy,
        "fluid_volume_m3": fluid_volume,
        "solid_volume_m3": solid_volume,
        "fluid_volume_mean_temperature_K": fluid_temperature_mean,
        "solid_volume_mean_temperature_K": solid_temperature_mean,
        "inlet_mass_rate_kg_s": inlet_mass,
        "outlet_mass_rate_kg_s": outlet_mass,
        "mass_imbalance_kg_s": inlet_mass + outlet_mass,
        "inlet_enthalpy_outward_W": inlet_enthalpy,
        "outlet_enthalpy_outward_W": outlet_enthalpy,
        "net_outward_enthalpy_W": net_outward_enthalpy,
        "inlet_bulk_temperature_K": inlet_bulk_temperature,
        "outlet_bulk_temperature_K": outlet_bulk_temperature,
        "heat_bc_faces": total["boundary_faces"].get(heat_bc, 0),
        "heat_bc_area_m2": heat_area,
        "applied_power_W": applied_power,
        "fluid_tetrahedra": total["tetrahedra_by_material"].get(0, 0),
        "solid_tetrahedra": sum(
            count for material, count in total["tetrahedra_by_material"].items()
            if material != 0
        ),
    }


def nice_ticks(minimum, maximum, count):
    if not math.isfinite(minimum) or not math.isfinite(maximum):
        return []
    if maximum <= minimum:
        return [minimum]
    raw = (maximum - minimum) / float(max(1, count - 1))
    exponent = math.floor(math.log10(raw))
    fraction = raw / (10.0 ** exponent)
    if fraction <= 1.0:
        nice = 1.0
    elif fraction <= 2.0:
        nice = 2.0
    elif fraction <= 5.0:
        nice = 5.0
    else:
        nice = 10.0
    step = nice * (10.0 ** exponent)
    start = math.floor(minimum / step) * step
    end = math.ceil(maximum / step) * step
    ticks = []
    value = start
    while value <= end + 0.5 * step:
        ticks.append(value)
        value += step
    return ticks


def escape_xml(text):
    return str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def write_svg(path, title, x_label, y_label, series, width=1200, height=720):
    left, right, top, bottom = 115.0, 45.0, 70.0, 95.0
    plot_width = width - left - right
    plot_height = height - top - bottom

    all_points = []
    for item in series:
        all_points.extend((x, y) for x, y in item["points"] if math.isfinite(x) and math.isfinite(y))
    if not all_points:
        raise RuntimeError("Cannot plot empty series for {}".format(path))
    xmin = min(point[0] for point in all_points)
    xmax = max(point[0] for point in all_points)
    ymin = min(point[1] for point in all_points)
    ymax = max(point[1] for point in all_points)
    if xmax <= xmin:
        xmax = xmin + 1.0
    if ymax <= ymin:
        padding = max(abs(ymin) * 0.05, 1.0)
        ymin -= padding
        ymax += padding
    else:
        padding = 0.06 * (ymax - ymin)
        ymin -= padding
        ymax += padding

    x_ticks = nice_ticks(xmin, xmax, 7)
    y_ticks = nice_ticks(ymin, ymax, 7)
    if x_ticks:
        xmin, xmax = x_ticks[0], x_ticks[-1]
    if y_ticks:
        ymin, ymax = y_ticks[0], y_ticks[-1]

    def sx(value):
        return left + (value - xmin) * plot_width / (xmax - xmin)

    def sy(value):
        return top + (ymax - value) * plot_height / (ymax - ymin)

    palette = ("#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#17becf", "#7f7f7f")
    lines = []
    lines.append('<?xml version="1.0" encoding="UTF-8"?>')
    lines.append('<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" viewBox="0 0 {} {}">'.format(width, height, width, height))
    lines.append('<rect width="100%" height="100%" fill="white"/>')
    lines.append('<text x="{}" y="36" text-anchor="middle" font-family="Arial" font-size="24" font-weight="bold">{}</text>'.format(width / 2.0, escape_xml(title)))

    for tick in x_ticks:
        x = sx(tick)
        lines.append('<line x1="{0:.3f}" y1="{1:.3f}" x2="{0:.3f}" y2="{2:.3f}" stroke="#dddddd" stroke-width="1"/>'.format(x, top, top + plot_height))
        lines.append('<text x="{:.3f}" y="{}" text-anchor="middle" font-family="Arial" font-size="14">{:.6g}</text>'.format(x, top + plot_height + 28, tick))
    for tick in y_ticks:
        y = sy(tick)
        lines.append('<line x1="{0:.3f}" y1="{1:.3f}" x2="{2:.3f}" y2="{1:.3f}" stroke="#dddddd" stroke-width="1"/>'.format(left, y, left + plot_width))
        lines.append('<text x="{}" y="{:.3f}" text-anchor="end" dominant-baseline="middle" font-family="Arial" font-size="14">{:.6g}</text>'.format(left - 12, y, tick))

    lines.append('<rect x="{}" y="{}" width="{}" height="{}" fill="none" stroke="black" stroke-width="1.5"/>'.format(left, top, plot_width, plot_height))

    for index, item in enumerate(series):
        colour = palette[index % len(palette)]
        points = [(x, y) for x, y in item["points"] if math.isfinite(x) and math.isfinite(y)]
        if not points:
            continue
        polyline = " ".join("{:.3f},{:.3f}".format(sx(x), sy(y)) for x, y in points)
        lines.append('<polyline points="{}" fill="none" stroke="{}" stroke-width="2.5"/>'.format(polyline, colour))
        for x, y in points:
            lines.append('<circle cx="{:.3f}" cy="{:.3f}" r="3.2" fill="{}"/>'.format(sx(x), sy(y), colour))

    legend_x = left + 12
    legend_y = top + 16
    for index, item in enumerate(series):
        colour = palette[index % len(palette)]
        y = legend_y + index * 24
        lines.append('<line x1="{0}" y1="{1}" x2="{2}" y2="{1}" stroke="{3}" stroke-width="3"/>'.format(legend_x, y, legend_x + 28, colour))
        lines.append('<text x="{}" y="{}" dominant-baseline="middle" font-family="Arial" font-size="14">{}</text>'.format(legend_x + 38, y, escape_xml(item["label"])))

    lines.append('<text x="{}" y="{}" text-anchor="middle" font-family="Arial" font-size="18">{}</text>'.format(left + plot_width / 2.0, height - 30, escape_xml(x_label)))
    lines.append('<text x="28" y="{}" text-anchor="middle" font-family="Arial" font-size="18" transform="rotate(-90 28 {})">{}</text>'.format(top + plot_height / 2.0, top + plot_height / 2.0, escape_xml(y_label)))
    lines.append('</svg>')
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pvtu-dir", required=True, type=Path)
    parser.add_argument("--partition-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--dt", required=True, type=float)
    parser.add_argument("--heat-flux", required=True, type=float)
    parser.add_argument("--heat-bc", type=int, default=532)
    parser.add_argument("--inlet-bc", type=int, default=511)
    parser.add_argument("--outlet-bc", type=int, default=520)
    args = parser.parse_args()

    pvtu_dir = args.pvtu_dir.resolve()
    output_dir = args.output_dir.resolve()
    partition_root = args.partition_root.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    pvtu_files = sorted(pvtu_dir.glob("*_step_*.pvtu"), key=iteration_from_name)
    if len(pvtu_files) < 2:
        raise RuntimeError("At least two PVTU states are required")

    properties, material_file, material_copies = read_material_properties(partition_root)
    print("CBS3D CHT energy-history post-processing")
    print("PVTU directory      : {}".format(pvtu_dir))
    print("Saved states        : {}".format(len(pvtu_files)))
    print("Material file       : {}".format(material_file))
    print("Material copies     : {}".format(material_copies))
    print("Iteration dt        : {:.17g}".format(args.dt))
    print("Applied heat flux   : {:.17g} W/m^2".format(args.heat_flux))

    states = []
    for state_index, pvtu in enumerate(pvtu_files, start=1):
        iteration = iteration_from_name(pvtu)
        print("\nProcessing state {:02d}/{:02d}: {}".format(state_index, len(pvtu_files), pvtu.name))
        row = process_state(
            pvtu,
            properties,
            args.heat_bc,
            args.inlet_bc,
            args.outlet_bc,
            args.heat_flux,
        )
        row["iteration"] = iteration
        row["pseudo_time_s"] = iteration * args.dt
        states.append(row)

    initial_fluid = states[0]["fluid_energy_J"]
    initial_solid = states[0]["solid_energy_J"]
    initial_total = states[0]["total_energy_J"]
    for row in states:
        row["fluid_energy_change_J"] = row["fluid_energy_J"] - initial_fluid
        row["solid_energy_change_J"] = row["solid_energy_J"] - initial_solid
        row["total_energy_change_J"] = row["total_energy_J"] - initial_total

    intervals = []
    for index in range(1, len(states)):
        previous = states[index - 1]
        current = states[index]
        delta_time = current["pseudo_time_s"] - previous["pseudo_time_s"]
        if delta_time <= 0.0:
            raise RuntimeError("Non-positive saved-state interval")
        dfluid_dt = (current["fluid_energy_J"] - previous["fluid_energy_J"]) / delta_time
        dsolid_dt = (current["solid_energy_J"] - previous["solid_energy_J"]) / delta_time
        dtotal_dt = dfluid_dt + dsolid_dt
        applied = 0.5 * (previous["applied_power_W"] + current["applied_power_W"])
        enthalpy_out = 0.5 * (
            previous["net_outward_enthalpy_W"] + current["net_outward_enthalpy_W"]
        )
        q_interface_from_solid = applied - dsolid_dt
        q_interface_from_fluid = dfluid_dt + enthalpy_out
        global_residual = applied - enthalpy_out - dtotal_dt
        interface_difference = q_interface_from_fluid - q_interface_from_solid
        denominator = max(abs(q_interface_from_solid), abs(q_interface_from_fluid), 1.0e-30)
        intervals.append({
            "iteration_start": previous["iteration"],
            "iteration_end": current["iteration"],
            "pseudo_time_mid_s": 0.5 * (previous["pseudo_time_s"] + current["pseudo_time_s"]),
            "delta_time_s": delta_time,
            "dE_fluid_dt_W": dfluid_dt,
            "dE_solid_dt_W": dsolid_dt,
            "dE_total_dt_W": dtotal_dt,
            "applied_power_W": applied,
            "net_outward_enthalpy_W": enthalpy_out,
            "interface_power_from_solid_W": q_interface_from_solid,
            "interface_power_from_fluid_W": q_interface_from_fluid,
            "interface_power_difference_W": interface_difference,
            "interface_power_relative_difference": abs(interface_difference) / denominator,
            "global_energy_residual_W": global_residual,
            "global_energy_relative_residual": abs(global_residual) / max(abs(applied), 1.0e-30),
        })

    state_csv = output_dir / "cht_energy_states.csv"
    interval_csv = output_dir / "cht_energy_intervals.csv"
    summary_path = output_dir / "cht_energy_history_summary.txt"
    energy_svg = output_dir / "cht_energy_storage_history.svg"
    power_svg = output_dir / "cht_power_balance_history.svg"
    bulk_svg = output_dir / "cht_bulk_temperature_history.svg"

    state_fields = [
        "iteration", "pseudo_time_s",
        "fluid_energy_J", "solid_energy_J", "total_energy_J",
        "fluid_energy_change_J", "solid_energy_change_J", "total_energy_change_J",
        "fluid_volume_m3", "solid_volume_m3",
        "fluid_volume_mean_temperature_K", "solid_volume_mean_temperature_K",
        "inlet_mass_rate_kg_s", "outlet_mass_rate_kg_s", "mass_imbalance_kg_s",
        "inlet_enthalpy_outward_W", "outlet_enthalpy_outward_W",
        "net_outward_enthalpy_W", "inlet_bulk_temperature_K",
        "outlet_bulk_temperature_K", "heat_bc_faces", "heat_bc_area_m2",
        "applied_power_W", "fluid_tetrahedra", "solid_tetrahedra",
    ]
    with state_csv.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=state_fields)
        writer.writeheader()
        for row in states:
            writer.writerow({field: row[field] for field in state_fields})

    interval_fields = [
        "iteration_start", "iteration_end", "pseudo_time_mid_s", "delta_time_s",
        "dE_fluid_dt_W", "dE_solid_dt_W", "dE_total_dt_W",
        "applied_power_W", "net_outward_enthalpy_W",
        "interface_power_from_solid_W", "interface_power_from_fluid_W",
        "interface_power_difference_W", "interface_power_relative_difference",
        "global_energy_residual_W", "global_energy_relative_residual",
    ]
    with interval_csv.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=interval_fields)
        writer.writeheader()
        for row in intervals:
            writer.writerow({field: row[field] for field in interval_fields})

    write_svg(
        energy_svg,
        "CBS3D CHT thermal-energy accumulation",
        "Accumulated iteration time, s",
        "Energy change, kJ",
        [
            {"label": "Fluid", "points": [(row["pseudo_time_s"], row["fluid_energy_change_J"] / 1000.0) for row in states]},
            {"label": "Solid", "points": [(row["pseudo_time_s"], row["solid_energy_change_J"] / 1000.0) for row in states]},
            {"label": "Total", "points": [(row["pseudo_time_s"], row["total_energy_change_J"] / 1000.0) for row in states]},
        ],
    )
    write_svg(
        power_svg,
        "CBS3D CHT interval power balance",
        "Accumulated iteration time, s",
        "Power, kW",
        [
            {"label": "Applied BC 532", "points": [(row["pseudo_time_mid_s"], row["applied_power_W"] / 1000.0) for row in intervals]},
            {"label": "Solid storage", "points": [(row["pseudo_time_mid_s"], row["dE_solid_dt_W"] / 1000.0) for row in intervals]},
            {"label": "Interface from solid", "points": [(row["pseudo_time_mid_s"], row["interface_power_from_solid_W"] / 1000.0) for row in intervals]},
            {"label": "Interface from fluid", "points": [(row["pseudo_time_mid_s"], row["interface_power_from_fluid_W"] / 1000.0) for row in intervals]},
            {"label": "Net outlet enthalpy", "points": [(row["pseudo_time_mid_s"], row["net_outward_enthalpy_W"] / 1000.0) for row in intervals]},
            {"label": "Global residual", "points": [(row["pseudo_time_mid_s"], row["global_energy_residual_W"] / 1000.0) for row in intervals]},
        ],
    )
    write_svg(
        bulk_svg,
        "CBS3D coolant bulk-temperature history",
        "Accumulated iteration time, s",
        "Bulk temperature, K",
        [
            {"label": "Inlet", "points": [(row["pseudo_time_s"], row["inlet_bulk_temperature_K"]) for row in states]},
            {"label": "Outlet", "points": [(row["pseudo_time_s"], row["outlet_bulk_temperature_K"]) for row in states]},
        ],
    )

    final_state = states[-1]
    final_interval = intervals[-1]
    max_global_relative = max(row["global_energy_relative_residual"] for row in intervals)
    max_interface_relative = max(row["interface_power_relative_difference"] for row in intervals)

    lines = []
    lines.append("CBS3D CHT energy-history summary")
    lines.append("PVTU directory                    : {}".format(pvtu_dir))
    lines.append("Saved states                      : {}".format(len(states)))
    lines.append("Iteration range                   : {} to {}".format(states[0]["iteration"], final_state["iteration"]))
    lines.append("Accumulated iteration time [s]    : {:.17e}".format(final_state["pseudo_time_s"]))
    lines.append("Heat-flux BC                      : {}".format(args.heat_bc))
    lines.append("Heat-flux area [m^2]              : {:.17e}".format(final_state["heat_bc_area_m2"]))
    lines.append("Applied heat flux [W/m^2]         : {:.17e}".format(args.heat_flux))
    lines.append("Applied power [W]                 : {:.17e}".format(final_state["applied_power_W"]))
    lines.append("")
    lines.append("Final saved state")
    lines.append("  fluid mean temperature [K]      : {:.17e}".format(final_state["fluid_volume_mean_temperature_K"]))
    lines.append("  solid mean temperature [K]      : {:.17e}".format(final_state["solid_volume_mean_temperature_K"]))
    lines.append("  inlet mass rate [kg/s]          : {:.17e}".format(final_state["inlet_mass_rate_kg_s"]))
    lines.append("  outlet mass rate [kg/s]         : {:.17e}".format(final_state["outlet_mass_rate_kg_s"]))
    lines.append("  mass imbalance [kg/s]           : {:.17e}".format(final_state["mass_imbalance_kg_s"]))
    lines.append("  inlet bulk temperature [K]      : {:.17e}".format(final_state["inlet_bulk_temperature_K"]))
    lines.append("  outlet bulk temperature [K]     : {:.17e}".format(final_state["outlet_bulk_temperature_K"]))
    lines.append("  net outward enthalpy [W]        : {:.17e}".format(final_state["net_outward_enthalpy_W"]))
    lines.append("")
    lines.append("Final interval")
    for key in (
        "dE_fluid_dt_W", "dE_solid_dt_W", "dE_total_dt_W",
        "interface_power_from_solid_W", "interface_power_from_fluid_W",
        "interface_power_difference_W", "interface_power_relative_difference",
        "global_energy_residual_W", "global_energy_relative_residual",
    ):
        lines.append("  {:34s}: {:.17e}".format(key, final_interval[key]))
    lines.append("")
    lines.append("Maximum interval diagnostics")
    lines.append("  global relative residual        : {:.17e}".format(max_global_relative))
    lines.append("  interface relative difference   : {:.17e}".format(max_interface_relative))
    lines.append("")
    lines.append("Outputs")
    lines.append("  {}".format(state_csv))
    lines.append("  {}".format(interval_csv))
    lines.append("  {}".format(energy_svg))
    lines.append("  {}".format(power_svg))
    lines.append("  {}".format(bulk_svg))
    summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n" + "\n".join(lines))
    print("\nPASS: CHT energy history and diagnostic plots generated")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        print("ERROR: {}".format(error), file=sys.stderr)
        sys.exit(2)
