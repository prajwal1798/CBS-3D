#!/usr/bin/env python3
"""Validate CBS3D DD-5B appended-raw-binary VTU/PVTU/PVD output.

The checker understands VTK XML files whose AppendedData payload is raw binary and
therefore cannot be parsed as ordinary XML. It validates every UInt64 block
header and offset, the expected CBS3D field layout, connectivity/offset
consistency, finite Float64 fields, PVTU piece references and the PVD entry.

Compatible with Python 3.6 on Swansea Sunbird.
"""

import argparse
import math
from pathlib import Path
import re
import struct
import sys
import xml.etree.ElementTree as ET

TYPE_BYTES = {
    "Float64": 8,
    "Int64": 8,
    "Int32": 4,
    "UInt8": 1,
}

EXPECTED_ARRAYS = (
    ("__points__", "Float64", 3, "point"),
    ("connectivity", "Int32", 1, "connectivity"),
    ("offsets", "Int64", 1, "cell"),
    ("types", "UInt8", 1, "cell"),
    ("pressure", "Float64", 1, "point"),
    ("temperature", "Float64", 1, "point"),
    ("velocity", "Float64", 3, "point"),
    ("velocity_magnitude", "Float64", 1, "point"),
    ("global_node_id", "Int64", 1, "point"),
    ("owner_rank", "Int32", 1, "point"),
    ("is_owned", "UInt8", 1, "point"),
    ("vtkGhostType", "UInt8", 1, "point"),
    ("node_domain_kind", "Int32", 1, "point"),
    ("velocity_bc_type", "Int32", 1, "point"),
    ("pressure_fixed", "UInt8", 1, "point"),
    ("cell_kind", "UInt8", 1, "cell"),
    ("global_element_id", "Int64", 1, "cell"),
    ("material_id", "Int32", 1, "cell"),
    ("bc_id", "Int32", 1, "cell"),
    ("parent_global_element", "Int64", 1, "cell"),
)

ATTRIBUTE_PATTERN = re.compile(r'(\w+)="([^"]*)"')
DATA_ARRAY_PATTERN = re.compile(r"<DataArray\b([^>]*)/>")
PIECE_PATTERN = re.compile(
    r'<Piece\s+NumberOfPoints="(\d+)"\s+NumberOfCells="(\d+)">'
)


def step_tag(iteration):
    return "step_{:08d}".format(iteration)


def attributes(text):
    return dict(ATTRIBUTE_PATTERN.findall(text))


def parse_piece(path, check_finite):
    raw = path.read_bytes()
    marker = b'<AppendedData encoding="raw">'
    marker_index = raw.find(marker)
    if marker_index < 0:
        raise ValueError("{}: raw AppendedData marker is missing".format(path))

    underscore_index = raw.find(b"_", marker_index + len(marker))
    if underscore_index < 0:
        raise ValueError("{}: AppendedData underscore is missing".format(path))

    try:
        header = raw[:underscore_index].decode("ascii")
    except UnicodeDecodeError as error:
        raise ValueError("{}: VTU XML header is not ASCII: {}".format(path, error))

    if 'byte_order="LittleEndian"' not in header:
        raise ValueError("{}: byte_order is not LittleEndian".format(path))
    if 'header_type="UInt64"' not in header:
        raise ValueError("{}: header_type is not UInt64".format(path))

    piece_match = PIECE_PATTERN.search(header)
    if piece_match is None:
        raise ValueError("{}: Piece point/cell counts are missing".format(path))

    point_count = int(piece_match.group(1))
    cell_count = int(piece_match.group(2))
    payload_start = underscore_index + 1

    parsed_arrays = []
    for index, match in enumerate(DATA_ARRAY_PATTERN.finditer(header)):
        attrs = attributes(match.group(1))
        if attrs.get("format") != "appended":
            raise ValueError(
                "{}: DataArray {} is not appended".format(path, index)
            )
        if "offset" not in attrs or "type" not in attrs:
            raise ValueError(
                "{}: DataArray {} lacks type or offset".format(path, index)
            )
        name = attrs.get("Name", "__points__" if index == 0 else "")
        components = int(attrs.get("NumberOfComponents", "1"))
        parsed_arrays.append(
            {
                "name": name,
                "type": attrs["type"],
                "components": components,
                "offset": int(attrs["offset"]),
            }
        )

    if len(parsed_arrays) != len(EXPECTED_ARRAYS):
        raise ValueError(
            "{}: expected {} appended arrays, found {}".format(
                path,
                len(EXPECTED_ARRAYS),
                len(parsed_arrays),
            )
        )

    previous_end_offset = 0
    blocks = {}

    for index, (parsed, expected) in enumerate(zip(parsed_arrays, EXPECTED_ARRAYS)):
        expected_name, expected_type, expected_components, association = expected
        if parsed["name"] != expected_name:
            raise ValueError(
                "{}: array {} name is {!r}, expected {!r}".format(
                    path, index, parsed["name"], expected_name
                )
            )
        if parsed["type"] != expected_type:
            raise ValueError(
                "{}: array {} type is {}, expected {}".format(
                    path, expected_name, parsed["type"], expected_type
                )
            )
        if parsed["components"] != expected_components:
            raise ValueError(
                "{}: array {} components are {}, expected {}".format(
                    path,
                    expected_name,
                    parsed["components"],
                    expected_components,
                )
            )
        if parsed["offset"] != previous_end_offset:
            raise ValueError(
                "{}: array {} offset is {}, expected {}".format(
                    path, expected_name, parsed["offset"], previous_end_offset
                )
            )

        block_header = payload_start + parsed["offset"]
        if block_header + 8 > len(raw):
            raise ValueError(
                "{}: array {} block header is truncated".format(path, expected_name)
            )
        block_bytes = struct.unpack_from("<Q", raw, block_header)[0]
        data_start = block_header + 8
        data_end = data_start + block_bytes
        if data_end > len(raw):
            raise ValueError(
                "{}: array {} payload is truncated".format(path, expected_name)
            )

        if association == "point":
            expected_values = point_count * expected_components
            expected_bytes = expected_values * TYPE_BYTES[expected_type]
            if block_bytes != expected_bytes:
                raise ValueError(
                    "{}: array {} has {} bytes, expected {}".format(
                        path, expected_name, block_bytes, expected_bytes
                    )
                )
        elif association == "cell":
            expected_values = cell_count * expected_components
            expected_bytes = expected_values * TYPE_BYTES[expected_type]
            if block_bytes != expected_bytes:
                raise ValueError(
                    "{}: array {} has {} bytes, expected {}".format(
                        path, expected_name, block_bytes, expected_bytes
                    )
                )
        elif association == "connectivity":
            if block_bytes % TYPE_BYTES[expected_type] != 0:
                raise ValueError(
                    "{}: connectivity byte count is not Int32-aligned".format(path)
                )

        if check_finite and expected_type == "Float64":
            payload = memoryview(raw)[data_start:data_end]
            for value_index, item in enumerate(struct.iter_unpack("<d", payload)):
                if not math.isfinite(item[0]):
                    raise ValueError(
                        "{}: array {} contains non-finite value at index {}".format(
                            path, expected_name, value_index
                        )
                    )

        blocks[expected_name] = {
            "start": data_start,
            "end": data_end,
            "bytes": block_bytes,
        }
        previous_end_offset = parsed["offset"] + 8 + block_bytes

    offsets_block = blocks["offsets"]
    offsets_payload = memoryview(raw)[offsets_block["start"]:offsets_block["end"]]
    cell_offsets = [value[0] for value in struct.iter_unpack("<q", offsets_payload)]
    if len(cell_offsets) != cell_count:
        raise ValueError("{}: cell-offset count mismatch".format(path))
    if cell_offsets:
        previous = 0
        for index, value in enumerate(cell_offsets):
            if value <= previous:
                raise ValueError(
                    "{}: cell offsets are not strictly increasing at {}".format(
                        path, index
                    )
                )
            increment = value - previous
            if increment not in (3, 4):
                raise ValueError(
                    "{}: cell {} has unsupported connectivity width {}".format(
                        path, index, increment
                    )
                )
            previous = value

        connectivity_values = blocks["connectivity"]["bytes"] // 4
        if cell_offsets[-1] != connectivity_values:
            raise ValueError(
                "{}: final cell offset {} does not match {} connectivity values".format(
                    path, cell_offsets[-1], connectivity_values
                )
            )

    payload_end = payload_start + previous_end_offset
    trailer = raw[payload_end:]
    if not trailer.startswith(b"\n  </AppendedData>\n</VTKFile>"):
        raise ValueError("{}: binary payload trailer is invalid".format(path))

    return {
        "points": point_count,
        "cells": cell_count,
        "bytes": len(raw),
        "payload_bytes": previous_end_offset,
        "arrays": len(parsed_arrays),
    }


def validate_pvtu(path, case_name, iteration, ranks):
    if not path.is_file():
        raise FileNotFoundError("PVTU file does not exist: {}".format(path))
    root = ET.parse(str(path)).getroot()
    pieces = root.findall("./PUnstructuredGrid/Piece")
    if len(pieces) != ranks:
        raise ValueError(
            "{}: expected {} Piece references, found {}".format(
                path, ranks, len(pieces)
            )
        )

    expected = {
        "{}_{}_rank_{:04d}.vtu".format(case_name, step_tag(iteration), rank)
        for rank in range(ranks)
    }
    actual = {piece.attrib.get("Source", "") for piece in pieces}
    if actual != expected:
        missing = sorted(expected.difference(actual))
        extra = sorted(actual.difference(expected))
        raise ValueError(
            "{}: piece-reference mismatch; missing={} extra={}".format(
                path, missing[:5], extra[:5]
            )
        )
    return actual


def validate_pvd(path, pvtu_name, expected_time):
    if not path.is_file():
        raise FileNotFoundError("PVD file does not exist: {}".format(path))
    root = ET.parse(str(path)).getroot()
    entries = root.findall("./Collection/DataSet")
    matches = [entry for entry in entries if entry.attrib.get("file") == pvtu_name]
    if len(matches) != 1:
        raise ValueError(
            "{}: expected one DataSet for {}, found {}".format(
                path, pvtu_name, len(matches)
            )
        )
    actual_time = float(matches[0].attrib.get("timestep", "nan"))
    if not math.isfinite(actual_time):
        raise ValueError("{}: PVD timestep is non-finite".format(path))
    if expected_time is not None:
        tolerance = 1.0e-12 * max(1.0, abs(expected_time))
        if abs(actual_time - expected_time) > tolerance:
            raise ValueError(
                "{}: PVD timestep {} differs from expected {}".format(
                    path, actual_time, expected_time
                )
            )
    return actual_time


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output_root", type=Path)
    parser.add_argument("case_name")
    parser.add_argument("iteration", type=int)
    parser.add_argument("ranks", type=int)
    parser.add_argument("--expected-time", type=float)
    parser.add_argument("--skip-finite-check", action="store_true")
    args = parser.parse_args()

    if args.iteration < 0:
        raise ValueError("iteration must be non-negative")
    if args.ranks < 2:
        raise ValueError("ranks must be at least two")

    root = args.output_root.resolve()
    pvtu_name = "{}_{}.pvtu".format(args.case_name, step_tag(args.iteration))
    pvtu_path = root / pvtu_name
    piece_names = validate_pvtu(
        pvtu_path,
        args.case_name,
        args.iteration,
        args.ranks,
    )

    summaries = []
    for piece_name in sorted(piece_names):
        piece_path = root / piece_name
        if not piece_path.is_file():
            raise FileNotFoundError("Referenced VTU piece is missing: {}".format(piece_path))
        summaries.append(
            parse_piece(piece_path, check_finite=not args.skip_finite_check)
        )

    pvd_path = root / "{}_distributed.pvd".format(args.case_name)
    pvd_time = validate_pvd(pvd_path, pvtu_name, args.expected_time)

    total_bytes = sum(summary["bytes"] for summary in summaries)
    total_payload = sum(summary["payload_bytes"] for summary in summaries)
    total_points = sum(summary["points"] for summary in summaries)
    total_cells = sum(summary["cells"] for summary in summaries)

    print("DD-5B appended-binary VTU validation report")
    print("  output root              : {}".format(root))
    print("  iteration                : {}".format(args.iteration))
    print("  PVD timestep             : {:.16g}".format(pvd_time))
    print("  rank-local pieces        : {}".format(len(summaries)))
    print("  appended arrays/piece    : {}".format(summaries[0]["arrays"]))
    print("  aggregate local points   : {}".format(total_points))
    print("  aggregate local cells    : {}".format(total_cells))
    print("  aggregate VTU bytes      : {}".format(total_bytes))
    print("  aggregate binary payload : {}".format(total_payload))
    print("  finite Float64 fields    : {}".format(
        "SKIPPED" if args.skip_finite_check else "PASS"
    ))
    print()
    print("DD-5B appended-binary VTU validation: PASS")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print("ERROR: {}".format(error), file=sys.stderr)
        raise SystemExit(2)
