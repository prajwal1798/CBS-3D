#!/usr/bin/env python3
"""Convert an unpartitioned Gmsh MSH4 mesh into a serial CBS3D case.

The Gmsh physical surface names/IDs are mapped to CBS3D boundary codes and the
physical volume names/IDs are mapped to material flags.  The output is the
native text format documented in ``docs/MESH_FORMAT.md``:

    <case>.plt
    <case>.bco
    <case>.material
    <case>.matprop
    <case>.par

This tool is for mesh-convergence studies on the tower.  Partitioned/MPI
meshes still use the MESH repository's ``export_rank_local_cbs.py`` workflow.
"""
from __future__ import annotations

import argparse
import shutil
from collections import Counter
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

try:
    import gmsh
except ImportError as exc:  # pragma: no cover - environment error
    raise SystemExit(
        "Install the Gmsh Python module: python -m pip install --upgrade gmsh"
    ) from exc

TET4 = 4
TRI3 = 2

VALID_BCS = {500, 501, 502, 503, 506, 507, 508, 510, 511, 520, 530, 532, 901, 902}
NAME_TO_BC = {
    "lid": 500,
    "moving_lid": 500,
    "moving_wall": 500,
    "wall": 530,
    "walls": 530,
    "fixed_wall": 530,
    "adiabatic_wall": 530,
    "inlet": 511,
    "mass_flow_inlet": 511,
    "velocity_inlet": 510,
    "outlet": 520,
    "pressure_outlet": 520,
    "interface": 901,
    "cht_interface": 901,
    "fluid_solid_interface": 901,
    "heat_flux": 532,
    "heat_flux_wall": 532,
    "heat_flux_top": 532,
}


def normalise(name: str) -> str:
    return name.strip().lower().replace("-", "_").replace(" ", "_")


def boundary_bc(physical_tag: int, name: str) -> int:
    key = normalise(name)
    if key in NAME_TO_BC:
        return NAME_TO_BC[key]
    if physical_tag in VALID_BCS:
        return physical_tag
    raise ValueError(
        f"Cannot map physical surface {physical_tag} name={name!r} to a CBS3D BC"
    )


def material_flag(physical_tag: int, name: str) -> int:
    key = normalise(name)
    if "fluid" in key or physical_tag == 10:
        return 0
    if "solid" in key or physical_tag == 20:
        return 1
    raise ValueError(
        f"Cannot map physical volume {physical_tag} name={name!r} to material 0/1"
    )


def get_elements(dim: int, entity: int, wanted_type: int) -> Tuple[List[int], List[int]]:
    types, tags_blocks, nodes_blocks = gmsh.model.mesh.getElements(dim, entity)
    tags: List[int] = []
    nodes: List[int] = []
    for etype, etags, enodes in zip(types, tags_blocks, nodes_blocks):
        if int(etype) == wanted_type:
            tags.extend(int(v) for v in etags)
            nodes.extend(int(v) for v in enodes)
    return tags, nodes


def tet_faces(conn: Sequence[int]) -> Tuple[Tuple[int, int, int], ...]:
    a, b, c, d = conn
    return (
        tuple(sorted((b, c, d))),
        tuple(sorted((a, d, c))),
        tuple(sorted((a, b, d))),
        tuple(sorted((a, c, b))),
    )


def physical_entity_map(dim: int) -> Dict[int, int]:
    result: Dict[int, int] = {}
    for pdim, ptag in gmsh.model.getPhysicalGroups(dim):
        if pdim != dim:
            continue
        for entity in gmsh.model.getEntitiesForPhysicalGroup(dim, ptag):
            if entity in result:
                raise RuntimeError(
                    f"Entity {entity} has more than one physical group at dimension {dim}"
                )
            result[entity] = int(ptag)
    return result


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Convert an unpartitioned Gmsh MSH4 mesh to a serial CBS3D case."
    )
    ap.add_argument("mesh", type=Path, help="Unpartitioned Gmsh MSH4 file")
    ap.add_argument("--case-name", required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--par", type=Path, required=True)
    ap.add_argument("--matprop", type=Path, required=True)
    ap.add_argument("--overwrite", action="store_true")
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    mesh = args.mesh.expanduser().resolve()
    par_file = args.par.expanduser().resolve()
    matprop_file = args.matprop.expanduser().resolve()
    out_dir = args.out_dir.expanduser().resolve()

    for required in (mesh, par_file, matprop_file):
        if not required.is_file():
            raise FileNotFoundError(required)

    if out_dir.exists() and any(out_dir.iterdir()):
        if not args.overwrite:
            raise FileExistsError(
                f"Output directory is not empty: {out_dir}. "
                "Use --overwrite only when deliberate."
            )
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    try:
        gmsh.open(str(mesh))

        if int(gmsh.model.getNumberOfPartitions()) > 1:
            raise RuntimeError(
                "This tool expects an unpartitioned mesh; use export_rank_local_cbs.py "
                "for partitioned/MPI meshes"
            )

        node_tags_raw, coords_raw, _ = gmsh.model.mesh.getNodes()
        if len(coords_raw) != 3 * len(node_tags_raw):
            raise RuntimeError("Global node coordinate length mismatch")
        coords = {
            int(tag): (
                float(coords_raw[3 * i]),
                float(coords_raw[3 * i + 1]),
                float(coords_raw[3 * i + 2]),
            )
            for i, tag in enumerate(node_tags_raw)
        }

        volume_physical = physical_entity_map(3)
        surface_physical = physical_entity_map(2)

        tetrahedra: List[Tuple[int, Tuple[int, int, int, int], int]] = []
        for _dim, entity in gmsh.model.getEntities(3):
            if entity not in volume_physical:
                raise RuntimeError(
                    f"Volume entity {entity} has no physical group; cannot assign a material"
                )
            ptag = volume_physical[entity]
            name = gmsh.model.getPhysicalName(3, ptag)
            mat = material_flag(ptag, name)
            tags, flat = get_elements(3, entity, TET4)
            for i, eid in enumerate(tags):
                conn = tuple(int(v) for v in flat[4 * i : 4 * i + 4])
                tetrahedra.append((eid, conn, mat))

        boundaries: List[Tuple[int, Tuple[int, int, int], int]] = []
        for _dim, entity in gmsh.model.getEntities(2):
            if entity not in surface_physical:
                continue
            ptag = surface_physical[entity]
            name = gmsh.model.getPhysicalName(2, ptag)
            bc = boundary_bc(ptag, name)
            tags, flat = get_elements(2, entity, TRI3)
            for i, tid in enumerate(tags):
                conn = tuple(int(v) for v in flat[3 * i : 3 * i + 3])
                boundaries.append((tid, conn, bc))

        if not tetrahedra:
            raise RuntimeError("No TET4 elements found")
        if not boundaries:
            raise RuntimeError("No physical TRI3 boundary faces found")

        tetrahedra.sort(key=lambda row: row[0])
        boundaries.sort(key=lambda row: row[0])

        face_to_parent: Dict[Tuple[int, int, int], int] = {}
        for lid, (_eid, conn, _mat) in enumerate(tetrahedra, start=1):
            for face in tet_faces(conn):
                if face in face_to_parent:
                    face_to_parent[face] = -1
                else:
                    face_to_parent[face] = lid

        boundary_records: List[Tuple[Tuple[int, int, int], int, int]] = []
        for _tid, conn, bc in boundaries:
            parent = face_to_parent.get(tuple(sorted(conn)), -1)
            if parent < 1:
                raise RuntimeError(
                    f"Boundary triangle {_tid} does not match one exterior tetrahedron face"
                )
            boundary_records.append((conn, parent, bc))

        used_nodes = sorted(
            {node for _eid, conn, _mat in tetrahedra for node in conn}
            | {node for conn, _parent, _bc in boundary_records for node in conn}
        )
        local_id = {global_id: local for local, global_id in enumerate(used_nodes, start=1)}

        base = out_dir / args.case_name
        with base.with_suffix(".plt").open("w", newline="\n") as f:
            f.write(f"{len(tetrahedra):10d}{len(used_nodes):10d}{len(boundary_records):10d}\n")
            for lid, (_eid, conn, _mat) in enumerate(tetrahedra, start=1):
                lconn = [local_id[n] for n in conn]
                f.write(
                    f"{lid:10d}{lconn[0]:10d}{lconn[1]:10d}{lconn[2]:10d}{lconn[3]:10d}\n"
                )
            for lid, gid in enumerate(used_nodes, start=1):
                x, y, z = coords[gid]
                f.write(f"{lid:10d}{x: .16e}{y: .16e}{z: .16e}\n")
            for conn, parent, bc in boundary_records:
                lconn = [local_id[n] for n in conn]
                f.write(
                    f"{lconn[0]:10d}{lconn[1]:10d}{lconn[2]:10d}{parent:10d}{bc:10d}\n"
                )

        unique_bcs = sorted({bc for _conn, _parent, bc in boundary_records})
        with base.with_suffix(".bco").open("w", newline="\n") as f:
            f.write(f"{len(unique_bcs):10d}{0:10d}\n")
            for bc in unique_bcs:
                f.write(f"{bc:10d}{bc:10d}\n")

        with base.with_suffix(".material").open("w", newline="\n") as f:
            for lid, (_eid, conn, mat) in enumerate(tetrahedra, start=1):
                lconn = [local_id[n] for n in conn]
                f.write(
                    f"{lid:10d}{lconn[0]:10d}{lconn[1]:10d}{lconn[2]:10d}{lconn[3]:10d}{mat:6d}\n"
                )

        shutil.copy2(par_file, base.with_suffix(".par"))
        shutil.copy2(matprop_file, base.with_suffix(".matprop"))

        bc_counts = Counter(bc for _conn, _parent, bc in boundary_records)
        mat_counts = Counter(mat for _eid, _conn, mat in tetrahedra)
        print(
            f"Wrote {base}.plt: {len(tetrahedra)} tets, {len(used_nodes)} nodes, "
            f"{len(boundary_records)} boundary faces"
        )
        print(f"  materials: {dict(sorted(mat_counts.items()))}")
        print(f"  boundary codes: {dict(sorted(bc_counts.items()))}")
        return 0
    finally:
        gmsh.finalize()


if __name__ == "__main__":
    raise SystemExit(main())
