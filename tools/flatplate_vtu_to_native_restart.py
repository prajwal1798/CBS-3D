#!/usr/bin/env python3
"""Convert distributed CBS3D ASCII VTU output to a native restart checkpoint.

This utility exists for one specific recovery path: a completed distributed
wall-resolved flat-plate run whose VTU pieces contain the full Spalart-Allmaras
state, but for which CBS3D_CHECKPOINT_EVERY was not enabled during the run.

Unlike the generic legacy-VTU importer in RestartIO.cpp, this converter preserves
``nu_tilde``.  That is essential for a genuine continuation of an SA benchmark;
reinitialising the turbulence working variable would not be a restart.

The binary layout written here is intentionally identical to RestartIO.cpp
restart format version 1.  The next solver run still performs all native restart
checks: MPI size, global mesh sizes, rank-local node ordering/hash, case name and
turbulence/temperature flags.

Python 3.6 compatible; no third-party packages are required.
"""

from __future__ import print_function

import argparse
import math
import os
import re
import shutil
import struct
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


RESTART_VERSION = 1
RESTART_MAGIC = b"CBS3D_RESTART_1\x00"
RESTART_END_MARKER = 0xC0DEC0DEC0DEC0DE
FNV_OFFSET_BASIS = 1469598103934665603
FNV_PRIME = 1099511628211
MASK64 = (1 << 64) - 1


def fail(message):
    raise SystemExit("FATAL: " + message)


def step_tag(iteration):
    return "step_{:08d}".format(iteration)


def rank_tag(rank):
    return "rank_{:04d}".format(rank)


def parse_ascii_numbers(data_array, cast, expected, name, path):
    text = data_array.text or ""
    tokens = text.split()
    if len(tokens) != expected:
        fail(
            "{} in {} has {} values; expected {}".format(
                name, path, len(tokens), expected
            )
        )

    values = []
    for token in tokens:
        try:
            value = cast(token)
        except Exception:
            fail("{} in {} contains invalid token {!r}".format(name, path, token))
        if cast is float and not math.isfinite(value):
            fail("{} in {} contains a non-finite value".format(name, path))
        values.append(value)
    return values


def find_named_array(parent, name, path):
    for array in parent.findall("DataArray"):
        if array.get("Name") == name:
            fmt = (array.get("format") or "ascii").lower()
            if fmt != "ascii":
                fail("{} in {} is not ASCII".format(name, path))
            return array
    fail("{} is missing DataArray {}".format(path, name))


def read_piece(path, state=False):
    try:
        root = ET.parse(str(path)).getroot()
    except Exception as exc:
        fail("cannot parse {}: {}".format(path, exc))

    piece = root.find(".//Piece")
    if piece is None:
        fail("{} has no VTK Piece".format(path))

    try:
        npoin = int(piece.get("NumberOfPoints"))
    except Exception:
        fail("{} has invalid NumberOfPoints".format(path))
    if npoin < 1:
        fail("{} has no points".format(path))

    point_data = piece.find("PointData")
    cell_data = piece.find("CellData")
    if point_data is None or cell_data is None:
        fail("{} is missing PointData or CellData".format(path))

    global_ids = parse_ascii_numbers(
        find_named_array(point_data, "global_node_id", path),
        int,
        npoin,
        "global_node_id",
        path,
    )
    if len(set(global_ids)) != npoin or min(global_ids) < 1:
        fail("{} has invalid/duplicate rank-local global_node_id values".format(path))

    global_element_array = find_named_array(cell_data, "global_element_id", path)
    global_element_tokens = (global_element_array.text or "").split()
    global_elements = []
    for token in global_element_tokens:
        try:
            value = int(token)
        except Exception:
            fail("global_element_id in {} contains invalid token {!r}".format(path, token))
        if value > 0:
            global_elements.append(value)

    result = {
        "npoin": npoin,
        "global_ids": global_ids,
        "global_elements": global_elements,
    }

    if not state:
        return result

    velocity = parse_ascii_numbers(
        find_named_array(point_data, "velocity", path),
        float,
        3 * npoin,
        "velocity",
        path,
    )
    pressure = parse_ascii_numbers(
        find_named_array(point_data, "pressure", path),
        float,
        npoin,
        "pressure",
        path,
    )
    temperature = parse_ascii_numbers(
        find_named_array(point_data, "temperature", path),
        float,
        npoin,
        "temperature",
        path,
    )
    nu_tilde = parse_ascii_numbers(
        find_named_array(point_data, "nu_tilde", path),
        float,
        npoin,
        "nu_tilde",
        path,
    )

    result.update(
        {
            "velocity": velocity,
            "pressure": pressure,
            "temperature": temperature,
            "nu_tilde": nu_tilde,
        }
    )
    return result


def fnv1a_node_map(global_ids):
    value = FNV_OFFSET_BASIS
    for global_id in global_ids:
        if global_id < 0:
            fail("negative global node id cannot be hashed")
        raw = struct.pack("<Q", global_id)
        for byte in bytearray(raw):
            value ^= byte
            value = (value * FNV_PRIME) & MASK64
    return value


def write_rank_checkpoint(
    path,
    rank,
    mpi_size,
    global_npoin,
    iteration,
    physical_time,
    time_step,
    flags,
    state,
):
    npoin = state["npoin"]
    global_ids = state["global_ids"]
    velocity = state["velocity"]
    pressure = state["pressure"]
    temperature = state["temperature"]
    nu_tilde = state["nu_tilde"]
    map_hash = fnv1a_node_map(global_ids)

    with path.open("wb") as output:
        output.write(RESTART_MAGIC)
        output.write(
            struct.pack(
                "<IiiqqqddQI",
                RESTART_VERSION,
                rank,
                mpi_size,
                npoin,
                global_npoin,
                iteration,
                physical_time,
                time_step,
                map_hash,
                flags,
            )
        )

        for global_id in global_ids:
            output.write(struct.pack("<q", global_id))

        for value in velocity:
            output.write(struct.pack("<d", value))
        for value in pressure:
            output.write(struct.pack("<d", value))
        for value in temperature:
            output.write(struct.pack("<d", value))
        for value in nu_tilde:
            output.write(struct.pack("<d", value))

        output.write(struct.pack("<Q", RESTART_END_MARKER))
        output.flush()
        os.fsync(output.fileno())


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--case", default="flatplate")
    parser.add_argument("--iteration", type=int, required=True)
    parser.add_argument("--mpi-size", type=int, default=40)
    parser.add_argument("--physical-time", type=float, default=0.0)
    parser.add_argument("--time-step", type=float, default=0.0)
    parser.add_argument("--temperature-enabled", type=int, choices=(0, 1), default=0)
    parser.add_argument("--turbulence-enabled", type=int, choices=(0, 1), default=1)
    return parser.parse_args()


def main():
    args = parse_args()

    if sys.byteorder != "little":
        fail("native CBS3D restart conversion currently requires a little-endian host")
    if args.iteration < 1:
        fail("iteration must be positive")
    if args.mpi_size < 2:
        fail("distributed restart requires at least two ranks")
    if not math.isfinite(args.physical_time) or not math.isfinite(args.time_step):
        fail("physical time and time step must be finite")

    input_dir = Path(args.input_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    if not input_dir.is_dir():
        fail("input directory does not exist: {}".format(input_dir))

    tag = step_tag(args.iteration)
    pattern = re.compile(
        r"^{}_{}_rank_(\d{{4}})\.vtu$".format(re.escape(args.case), tag)
    )

    pieces = {}
    for path in input_dir.iterdir():
        match = pattern.match(path.name)
        if match:
            rank = int(match.group(1))
            pieces[rank] = path

    expected_ranks = set(range(args.mpi_size))
    if set(pieces) != expected_ranks:
        missing = sorted(expected_ranks.difference(pieces))
        extra = sorted(set(pieces).difference(expected_ranks))
        fail(
            "VTU rank set mismatch; missing={} extra={}".format(missing, extra)
        )

    all_nodes = set()
    all_elements = set()
    local_npoin = {}

    for rank in range(args.mpi_size):
        metadata = read_piece(pieces[rank], state=False)
        local_npoin[rank] = metadata["npoin"]
        all_nodes.update(metadata["global_ids"])
        all_elements.update(metadata["global_elements"])

    if not all_nodes or min(all_nodes) != 1 or max(all_nodes) != len(all_nodes):
        fail("global node ids are not contiguous 1..N")
    if not all_elements or min(all_elements) != 1 or max(all_elements) != len(all_elements):
        fail("global element ids are not contiguous 1..N")

    global_npoin = len(all_nodes)
    global_nelem = len(all_elements)

    flags = 0
    if args.temperature_enabled:
        flags |= 1
    if args.turbulence_enabled:
        flags |= 2

    tmp_dir = output_dir.with_name(output_dir.name + ".tmp")
    if tmp_dir.exists():
        shutil.rmtree(str(tmp_dir))
    tmp_dir.mkdir(parents=True)

    nu_tilde_min = None
    nu_tilde_max = None

    try:
        for rank in range(args.mpi_size):
            state = read_piece(pieces[rank], state=True)
            if state["npoin"] != local_npoin[rank]:
                fail("rank {} point count changed between conversion passes".format(rank))

            local_min = min(state["nu_tilde"])
            local_max = max(state["nu_tilde"])
            nu_tilde_min = local_min if nu_tilde_min is None else min(nu_tilde_min, local_min)
            nu_tilde_max = local_max if nu_tilde_max is None else max(nu_tilde_max, local_max)

            rank_path = tmp_dir / (
                "{}_{}_{}.cbsrst".format(args.case, tag, rank_tag(rank))
            )
            write_rank_checkpoint(
                rank_path,
                rank,
                args.mpi_size,
                global_npoin,
                args.iteration,
                args.physical_time,
                args.time_step,
                flags,
                state,
            )

        manifest = tmp_dir / "restart_manifest.txt"
        with manifest.open("w") as output:
            output.write("format_version {}\n".format(RESTART_VERSION))
            output.write("case_name {}\n".format(args.case))
            output.write("completed_iteration {}\n".format(args.iteration))
            output.write("mpi_size {}\n".format(args.mpi_size))
            output.write("global_npoin {}\n".format(global_npoin))
            output.write("global_nelem {}\n".format(global_nelem))
            output.write("physical_time {:.17g}\n".format(args.physical_time))
            output.write("time_step {:.17g}\n".format(args.time_step))
            output.flush()
            os.fsync(output.fileno())

        rank_files = list(tmp_dir.glob("*.cbsrst"))
        if len(rank_files) != args.mpi_size:
            fail("checkpoint conversion did not create all rank files")

        if output_dir.exists():
            shutil.rmtree(str(output_dir))
        tmp_dir.rename(output_dir)
    except Exception:
        if tmp_dir.exists():
            shutil.rmtree(str(tmp_dir))
        raise

    print("FLAT-PLATE VTU -> NATIVE RESTART: PASS")
    print("  source directory : {}".format(input_dir))
    print("  checkpoint       : {}".format(output_dir))
    print("  iteration        : {}".format(args.iteration))
    print("  MPI ranks        : {}".format(args.mpi_size))
    print("  global nodes     : {}".format(global_npoin))
    print("  global elements  : {}".format(global_nelem))
    print("  nu_tilde range   : [{:.12e}, {:.12e}]".format(nu_tilde_min, nu_tilde_max))
    print("  time step        : {:.12e}".format(args.time_step))


if __name__ == "__main__":
    main()
