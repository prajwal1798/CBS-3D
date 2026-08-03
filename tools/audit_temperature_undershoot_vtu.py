#!/usr/bin/env python3
"""Locate and classify nonphysical temperature undershoots in CBS3D VTU history.

This read-only audit tracks Tmin/Tmax and sub-reference node counts through all
saved MPI snapshots.  For the final snapshot it classifies cold nodes as
fluid-only, solid-only or conformal fluid-solid interface nodes, reports
boundary IDs, velocity, coordinates, and the spatial extent of the cold region.

Only Python's standard library is required.
"""

import argparse
import math
import re
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path


def fail(message):
    raise RuntimeError(message)


def named(parent, candidates, context):
    wanted = set(candidates)
    for node in parent.findall("DataArray"):
        if node.get("Name") in wanted:
            return node
    fail("{}: missing DataArray; expected one of {}".format(context, sorted(wanted)))


def optional_named(parent, candidates):
    wanted = set(candidates)
    for node in parent.findall("DataArray"):
        if node.get("Name") in wanted:
            return node
    return None


def floats(node, context):
    try:
        return [float(value) for value in (node.text or "").split()]
    except ValueError as exc:
        raise RuntimeError("{}: invalid floating-point DataArray".format(context)) from exc


def ints(node, context):
    try:
        return [int(value) for value in (node.text or "").split()]
    except ValueError as exc:
        raise RuntimeError("{}: invalid integer DataArray".format(context)) from exc


def vector_norm(v):
    return math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])


def discover_steps(directory):
    pattern = re.compile(r"step_(\d+)_rank_\d+\.vtu$")
    counts = Counter()
    for path in directory.glob("*.vtu"):
        match = pattern.search(path.name)
        if match:
            counts[int(match.group(1))] += 1
    if not counts:
        fail("No rank-local step VTU files found in {}".format(directory))
    complete = sorted(step for step, count in counts.items() if count == max(counts.values()))
    if not complete:
        fail("No complete VTU snapshots found")
    return complete, counts


def files_for_step(directory, step):
    files = sorted(directory.glob("*step_{:08d}_rank_*.vtu".format(step)))
    if not files:
        fail("No VTU pieces found for step {}".format(step))
    return files


def read_piece(path, with_topology=False):
    root = ET.parse(path).getroot()
    piece = root.find("./UnstructuredGrid/Piece")
    if piece is None:
        fail("{}: missing UnstructuredGrid/Piece".format(path))

    points_parent = piece.find("Points")
    point_data = piece.find("PointData")
    if points_parent is None or point_data is None:
        fail("{}: missing Points or PointData".format(path))

    point_array = points_parent.find("DataArray")
    if point_array is None:
        fail("{}: missing coordinates".format(path))
    flat_points = floats(point_array, str(path))
    if len(flat_points) % 3:
        fail("{}: coordinate array length is not divisible by three".format(path))
    points = [tuple(flat_points[i:i + 3]) for i in range(0, len(flat_points), 3)]

    gids = ints(named(point_data, ("global_node_id",), str(path)), str(path))
    owned = ints(named(point_data, ("is_owned",), str(path)), str(path))
    temperature = floats(named(point_data, ("temperature",), str(path)), str(path))
    flat_velocity = floats(named(point_data, ("velocity",), str(path)), str(path))

    domain_node = optional_named(point_data, ("node_domain_kind",))
    velocity_bc_node = optional_named(point_data, ("velocity_bc_type",))
    pressure_fixed_node = optional_named(point_data, ("pressure_fixed",))

    domain_kind = ints(domain_node, str(path)) if domain_node is not None else [-1] * len(points)
    velocity_bc = ints(velocity_bc_node, str(path)) if velocity_bc_node is not None else [-1] * len(points)
    pressure_fixed = ints(pressure_fixed_node, str(path)) if pressure_fixed_node is not None else [-1] * len(points)

    n = len(points)
    if not (
        len(gids) == n
        and len(owned) == n
        and len(temperature) == n
        and len(flat_velocity) == 3 * n
        and len(domain_kind) == n
        and len(velocity_bc) == n
        and len(pressure_fixed) == n
    ):
        fail("{}: inconsistent PointData lengths".format(path))

    result = {
        "points": points,
        "gids": gids,
        "owned": owned,
        "temperature": temperature,
        "velocity": [tuple(flat_velocity[i:i + 3]) for i in range(0, len(flat_velocity), 3)],
        "domain_kind": domain_kind,
        "velocity_bc": velocity_bc,
        "pressure_fixed": pressure_fixed,
    }

    if with_topology:
        cells_parent = piece.find("Cells")
        cell_data = piece.find("CellData")
        if cells_parent is None or cell_data is None:
            fail("{}: missing Cells or CellData".format(path))

        connectivity = ints(named(cells_parent, ("connectivity",), str(path)), str(path))
        offsets = ints(named(cells_parent, ("offsets",), str(path)), str(path))
        types = ints(named(cells_parent, ("types",), str(path)), str(path))
        material_id = ints(named(cell_data, ("material_id", "mat_id"), str(path)), str(path))
        bc_id = ints(named(cell_data, ("bc_id",), str(path)), str(path))

        if not (len(offsets) == len(types) == len(material_id) == len(bc_id)):
            fail("{}: inconsistent CellData lengths".format(path))

        result.update({
            "connectivity": connectivity,
            "offsets": offsets,
            "types": types,
            "material_id": material_id,
            "bc_id": bc_id,
        })

    return result


def node_class(materials):
    if 0 in materials and any(mid != 0 for mid in materials):
        return "interface"
    if 0 in materials:
        return "fluid-only"
    if materials:
        return "solid-only"
    return "unclassified"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vtu-dir", type=Path, required=True)
    parser.add_argument("--dt", type=float, required=True)
    parser.add_argument("--reference-temperature", type=float, default=300.0)
    parser.add_argument(
        "--steps",
        default="all",
        help="Comma-separated saved steps, or 'all' to inspect every complete snapshot",
    )
    parser.add_argument("--cold-tolerance", type=float, default=0.1)
    args = parser.parse_args()

    if not args.vtu_dir.is_dir():
        fail("VTU directory does not exist: {}".format(args.vtu_dir))
    if args.dt <= 0.0 or args.cold_tolerance <= 0.0:
        fail("dt and cold-tolerance must be positive")

    available_steps, counts = discover_steps(args.vtu_dir)
    max_piece_count = max(counts.values())

    if args.steps.strip().lower() == "all":
        steps = available_steps
    else:
        steps = sorted(set(int(value.strip()) for value in args.steps.split(",") if value.strip()))
        missing = [step for step in steps if counts.get(step, 0) != max_piece_count]
        if missing:
            fail("Requested steps are absent or incomplete: {}".format(missing))

    reference = args.reference_temperature
    thresholds = [
        reference - 1.0e-9,
        reference - 0.1,
        reference - 1.0,
        reference - 5.0,
        reference - 10.0,
        reference - 20.0,
        reference - 30.0,
    ]

    print("=== Temperature undershoot history ===")
    print("VTU directory             : {}".format(args.vtu_dir))
    print("MPI pieces per snapshot   : {}".format(max_piece_count))
    print("Reference temperature     : {:.9f} K".format(reference))
    print("Saved steps               : {}".format(",".join(str(step) for step in steps)))
    print()
    print(
        "{:>10s} {:>10s} {:>13s} {:>13s} {:>10s} {:>10s} {:>10s} {:>10s}".format(
            "step", "time[s]", "Tmin[K]", "Tmax[K]", "N<Tref", "N<299.9", "N<295", "N<280"
        )
    )

    history = []
    final_records = None
    final_step = steps[-1]
    material_touch = defaultdict(set)
    boundary_touch = defaultdict(set)

    for step in steps:
        final = step == final_step
        records = {} if final else None
        minimum = None
        maximum = None
        below = Counter()

        files = files_for_step(args.vtu_dir, step)
        if len(files) != max_piece_count:
            fail("Step {} has {} pieces; expected {}".format(step, len(files), max_piece_count))

        for path in files:
            piece = read_piece(path, with_topology=final)
            n = len(piece["points"])

            for i in range(n):
                if piece["owned"][i] != 1:
                    continue

                gid = piece["gids"][i]
                temperature = piece["temperature"][i]
                record = {
                    "gid": gid,
                    "temperature": temperature,
                    "point": piece["points"][i],
                    "velocity": piece["velocity"][i],
                    "domain_kind": piece["domain_kind"][i],
                    "velocity_bc": piece["velocity_bc"][i],
                    "pressure_fixed": piece["pressure_fixed"][i],
                    "file": path.name,
                }

                if minimum is None or temperature < minimum["temperature"]:
                    minimum = record
                if maximum is None or temperature > maximum["temperature"]:
                    maximum = record

                for threshold in thresholds:
                    if temperature < threshold:
                        below[threshold] += 1

                if final:
                    if gid in records:
                        fail("Duplicate owned global node {} at step {}".format(gid, step))
                    records[gid] = record

            if final:
                start = 0
                for cell_index, end in enumerate(piece["offsets"]):
                    local_nodes = piece["connectivity"][start:end]
                    start = end
                    cell_type = piece["types"][cell_index]
                    mid = piece["material_id"][cell_index]
                    bc = piece["bc_id"][cell_index]

                    if cell_type == 10:
                        for local_node in local_nodes:
                            material_touch[piece["gids"][local_node]].add(mid)
                    elif cell_type == 5 and bc != 0:
                        for local_node in local_nodes:
                            boundary_touch[piece["gids"][local_node]].add(bc)

        if minimum is None or maximum is None:
            fail("No owned points found at step {}".format(step))

        history.append((step, minimum, maximum, below))
        print(
            "{:10d} {:10.4f} {:13.6f} {:13.6f} {:10d} {:10d} {:10d} {:10d}".format(
                step,
                step * args.dt,
                minimum["temperature"],
                maximum["temperature"],
                below[thresholds[0]],
                below[reference - 0.1],
                below[reference - 5.0],
                below[reference - 20.0],
            )
        )

        if final:
            final_records = records

    assert final_records is not None

    cold_limit = reference - args.cold_tolerance
    cold_records = [record for record in final_records.values() if record["temperature"] < cold_limit]
    cold_records.sort(key=lambda record: record["temperature"])

    class_counts = Counter()
    domain_kind_counts = Counter()
    boundary_counts = Counter()
    speeds = []

    for record in cold_records:
        gid = record["gid"]
        class_counts[node_class(material_touch.get(gid, set()))] += 1
        domain_kind_counts[record["domain_kind"]] += 1
        for bc in boundary_touch.get(gid, set()):
            boundary_counts[bc] += 1
        speeds.append(vector_norm(record["velocity"]))

    print()
    print("=== Final-snapshot cold-region classification ===")
    print("Final step                : {}".format(final_step))
    print("Final nominal time        : {:.9f} s".format(final_step * args.dt))
    print("Cold criterion            : T < {:.9f} K".format(cold_limit))
    print("Cold owned nodes          : {}".format(len(cold_records)))
    print("Material classification   : {}".format(dict(sorted(class_counts.items()))))
    print("node_domain_kind counts   : {}".format(dict(sorted(domain_kind_counts.items()))))
    print("Boundary-node incidences  : {}".format(dict(sorted(boundary_counts.items()))))

    if cold_records:
        xs = [record["point"][0] for record in cold_records]
        ys = [record["point"][1] for record in cold_records]
        zs = [record["point"][2] for record in cold_records]
        print(
            "Cold-region bounding box : x=[{:.9e},{:.9e}] y=[{:.9e},{:.9e}] z=[{:.9e},{:.9e}]".format(
                min(xs), max(xs), min(ys), max(ys), min(zs), max(zs)
            )
        )
        print(
            "Cold-node speed range    : [{:.9e},{:.9e}] m/s".format(
                min(speeds), max(speeds)
            )
        )

    final_min = history[-1][1]
    final_max = history[-1][2]

    for label, record in (("Minimum", final_min), ("Maximum", final_max)):
        gid = record["gid"]
        print()
        print("{}-temperature node".format(label))
        print("  global node ID          : {}".format(gid))
        print("  temperature             : {:.12e} K".format(record["temperature"]))
        print("  coordinate              : ({:.12e}, {:.12e}, {:.12e}) m".format(*record["point"]))
        print("  velocity                : ({:.12e}, {:.12e}, {:.12e}) m/s".format(*record["velocity"]))
        print("  velocity magnitude      : {:.12e} m/s".format(vector_norm(record["velocity"])))
        print("  material class          : {}".format(node_class(material_touch.get(gid, set()))))
        print("  incident material IDs   : {}".format(sorted(material_touch.get(gid, set()))))
        print("  incident boundary IDs   : {}".format(sorted(boundary_touch.get(gid, set()))))
        print("  node_domain_kind        : {}".format(record["domain_kind"]))
        print("  velocity_bc_type        : {}".format(record["velocity_bc"]))
        print("  pressure_fixed          : {}".format(record["pressure_fixed"]))
        print("  owner VTU               : {}".format(record["file"]))

    print()
    print("Ten coldest owned nodes at final step:")
    for index, record in enumerate(cold_records[:10], 1):
        gid = record["gid"]
        print(
            "  {:2d}: gid={} T={:.9f} K xyz=({:.6e},{:.6e},{:.6e}) "
            "|u|={:.6e} class={} BCs={}".format(
                index,
                gid,
                record["temperature"],
                record["point"][0],
                record["point"][1],
                record["point"][2],
                vector_norm(record["velocity"]),
                node_class(material_touch.get(gid, set())),
                sorted(boundary_touch.get(gid, set())),
            )
        )

    print()
    print("Interpretation rule:")
    print("  With T_initial = T_inlet = T_ref, positive heat input, and no cooling source,")
    print("  any material undershoot below T_ref violates the continuous maximum principle.")
    print("  Exact global energy conservation does not exclude paired nodal undershoot/overshoot.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print("ERROR: {}".format(exc))
        raise SystemExit(1)
