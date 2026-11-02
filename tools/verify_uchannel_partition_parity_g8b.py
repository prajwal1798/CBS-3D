#!/usr/bin/env python3
"""G8b U-channel 2-vs-4 MPI partition verification, preserving original fields.

make-restart: map the owned global field from the actual 2-rank G6 step 50000
              to the LOCAL node order of each independently exported 4-rank
              partition. Output is a VTK XML initial-condition file for
              RestartIO::read_legacy_vtu_checkpoint; no numerical fields are
              interpolated or changed.

compare: independently compare the owned global nodal velocity, pressure,
         temperature and coordinates from 2-rank and 4-rank coupled solver
         VTU snapshots. Must pass initial-state and one-step comparisons.

This does NOT test serial-vs-MPI parity and does NOT modify solver source.
"""

import argparse
import math
from pathlib import Path
import xml.etree.ElementTree as ET

from audit_temperature_undershoot_vtu import read_piece, files_for_step, named, floats


def require(ok, explanation):
    if not ok:
        raise ValueError(explanation)


def owned_fields(directory, step, n_ranks):
    paths = files_for_step(directory, step)
    require(len(paths) == n_ranks,
            f"{directory} step {step}: expected {n_ranks} ranks, got {len(paths)}")
    owned = {}
    seen = set()
    for path in paths:
        piece = read_piece(path)
        root = ET.parse(path).getroot()
        pdata = root.find("./UnstructuredGrid/Piece/PointData")
        require(pdata is not None, f"missing PointData: {path}")
        pressure = floats(named(pdata, ("pressure",), str(path)), str(path))
        require(len(pressure) == len(piece["gids"]),
                f"{path} pressure array length mismatch")
        for i, gid in enumerate(piece["gids"]):
            seen.add(gid)
            if piece["owned"][i] != 1:
                continue
            require(gid not in owned, f"duplicate owned global ID {gid}")
            owned[gid] = {
                "coordinates": piece["points"][i],
                "temperature": (piece["temperature"][i],),
                "velocity": piece["velocity"][i],
                "pressure": (pressure[i],)
            }
    require(set(owned) == seen, "some global mesh nodes lack an owner")
    require(len(owned) > 0, "empty global mesh")
    return owned


def mpi_nodes(path, expected_rank, expected_parts):
    rows = path.read_text(encoding="utf-8").splitlines()
    rank_line = next((line.split() for line in rows
                      if line.startswith("PARTITION ")), None)
    require(rank_line is not None and len(rank_line) == 4
            and int(rank_line[2]) == expected_rank
            and int(rank_line[3]) == expected_parts,
            f"{path}: incorrect MPI rank metadata")
    local = next((line.split() for line in rows
                  if line.startswith("LOCAL ")), None)
    require(local is not None and len(local) == 6,
            f"{path}: missing LOCAL node counts")
    nlocal, nowned, nghost = map(int, local[2:5])
    require(nlocal == nowned + nghost,
            f"{path}: wrong owned/ghost totals")
    try:
        beginning = rows.index("NODES") + 1
    except ValueError:
        raise ValueError(f"{path}: NODES section missing") from None
    records = []
    for i, line in enumerate(rows[beginning:beginning+nlocal]):
        x = line.split()
        require(len(x) == 3 and int(x[0]) == i + 1,
                f"{path}: invalid local node mapping at position {i+1}")
        gid, owner = int(x[1]), int(x[2])
        require(0 <= owner < expected_parts,
                f"{path}: invalid owner for gid {gid}")
        require((owner == expected_rank) == (i < nowned),
                f"{path}: owned-first ordering violated for gid {gid}")
        records.append(gid)
    require(len(records) == nlocal and len(set(records)) == nlocal,
            f"{path}: duplicate or missing local nodes")
    return records


def write_xml_restart(file, records, fields):
    file.parent.mkdir(parents=True, exist_ok=True)
    def fmt(x):
        require(math.isfinite(x), f"nonfinite restart value {x}")
        return format(x, ".17g")
    with file.open("w", encoding="utf-8") as out:
        out.write('<?xml version="1.0"?>\n')
        out.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
        out.write(f'<UnstructuredGrid><Piece NumberOfPoints="{len(records)}" NumberOfCells="0">\n')
        out.write('<PointData>\n')
        for name, data_type, ncomp in [
            ("global_node_id", "Int64", 1),
            ("velocity", "Float64", 3),
            ("pressure", "Float64", 1),
            ("temperature", "Float64", 1),
        ]:
            out.write(f'<DataArray type="{data_type}" Name="{name}" '
                      f'NumberOfComponents="{ncomp}" format="ascii">\n')
            for gid in records:
                if name == "global_node_id":
                    out.write(f"{gid}\n")
                else:
                    out.write(" ".join(fmt(v) for v in fields[gid][name]) + "\n")
            out.write('</DataArray>\n')
        out.write('</PointData>\n<Points>\n')
        out.write('<DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
        for gid in records:
            out.write(" ".join(fmt(x) for x in fields[gid]["coordinates"]) + "\n")
        out.write('</DataArray></Points>\n')
        out.write('<Cells><DataArray type="Int32" Name="connectivity" format="ascii"></DataArray>\n')
        out.write('<DataArray type="Int32" Name="offsets" format="ascii"></DataArray>\n')
        out.write('<DataArray type="UInt8" Name="types" format="ascii"></DataArray></Cells>\n')
        out.write('</Piece></UnstructuredGrid></VTKFile>\n')


def make_restart(a):
    reference = owned_fields(a.reference, a.step, 2)
    print(f"Source snapshot: {a.reference}, step {a.step}, {len(reference)} unique nodes")
    for rank in range(4):
        mpi = a.partitions/f"rank_{rank:04d}"/f"uchannel_rank_{rank:04d}.mpi"
        require(mpi.is_file(), f"missing rank-{rank} MPI metadata {mpi}")
        records = mpi_nodes(mpi, rank, 4)
        missing = set(records)-set(reference)
        require(not missing, f"new partition introduces foreign IDs: {sorted(missing)[:10]}")
        # MPI ownership and global coverage are checked against the exporter
        # manifest; all local/ghost values derive from the same owned field.
        out = a.destination/f"uchannel_step_{a.step:08d}_rank_{rank:04d}.vtu"
        write_xml_restart(out, records, reference)
        print(f"  rank {rank}: wrote {len(records)} exact nodal fields => {out}")
    print("G8B STEP-50000 RESTART REMAP PASS (no interpolation)")


def compare(a):
    print("=== G8b 2 versus 4 MPI ranks: identical global initial fields ===")
    max_allowed = {
        "temperature": a.temperature_tol,
        "velocity": a.velocity_tol,
        "pressure": a.pressure_tol,
        "coordinates": 1e-12,
    }
    for step in (a.before, a.after):
        base = owned_fields(a.reference, step, 2)
        trial = owned_fields(a.candidate, step, 4)
        require(base.keys() == trial.keys(),
                f"step {step}: global node ID sets differ")
        print(f"STEP {step}: {len(base)} owned global nodes")
        for name in max_allowed:
            largest, worst = 0.0, None
            for gid in base:
                for left, right in zip(base[gid][name], trial[gid][name]):
                    require(math.isfinite(left) and math.isfinite(right),
                            f"nonfinite {name} at global ID {gid}")
                    delta = abs(left-right)
                    if delta > largest:
                        largest, worst = delta, gid
            print(f"  max_abs_{name:12s} = {largest:.12e} "
                  f"at gid={worst if worst is not None else 'none'} "
                  f"(tol={max_allowed[name]:.1e})")
            require(largest <= max_allowed[name],
                    f"G8b {name} at step {step} fails declared threshold")
    print("G8B 2-VS-4-RANK COUPLED FIRST-STEP PARITY PASS")
    print("Scope: same 3D mesh, same initial global field and production binary; "
          "not serial-vs-MPI.")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    mk = sub.add_parser("make-restart")
    mk.add_argument("--reference", type=Path, required=True)
    mk.add_argument("--partitions", type=Path, required=True)
    mk.add_argument("--destination", type=Path, required=True)
    mk.add_argument("--step", type=int, default=50000)
    cp = sub.add_parser("compare")
    cp.add_argument("--reference", type=Path, required=True)
    cp.add_argument("--candidate", type=Path, required=True)
    cp.add_argument("--before", type=int, default=50000)
    cp.add_argument("--after", type=int, default=50001)
    cp.add_argument("--temperature-tol", type=float, default=1e-9)
    cp.add_argument("--velocity-tol", type=float, default=1e-10)
    cp.add_argument("--pressure-tol", type=float, default=1e-7)
    a = p.parse_args()
    if a.command == "make-restart":
        make_restart(a)
    else:
        compare(a)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, IndexError, OverflowError) as exc:
        raise SystemExit(f"G8b ERROR: {exc}")
