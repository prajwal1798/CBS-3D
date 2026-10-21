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
import math
import shutil
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

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


def subtract(u: Sequence[float], v: Sequence[float]) -> Tuple[float, float, float]:
    return (u[0] - v[0], u[1] - v[1], u[2] - v[2])


def dot(u: Sequence[float], v: Sequence[float]) -> float:
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


def cross(u: Sequence[float], v: Sequence[float]) -> Tuple[float, float, float]:
    return (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )


def distance2(u: Sequence[float], v: Sequence[float]) -> float:
    dx = u[0] - v[0]
    dy = u[1] - v[1]
    dz = u[2] - v[2]
    return dx * dx + dy * dy + dz * dz


def outward_faces(
    conn: Sequence[int],
    coords: Dict[int, Tuple[float, float, float]],
) -> Tuple[Tuple[int, int, int], ...]:
    """Return the four tet faces with outward-pointing node order."""
    a, b, c, d = conn
    raw = (
        (a, b, c, d),
        (a, b, d, c),
        (a, c, d, b),
        (b, c, d, a),
    )
    result: List[Tuple[int, int, int]] = []
    for p0, p1, p2, opposite in raw:
        normal = cross(
            subtract(coords[p1], coords[p0]),
            subtract(coords[p2], coords[p0]),
        )
        if dot(normal, subtract(coords[opposite], coords[p0])) > 0.0:
            p1, p2 = p2, p1
        result.append((p0, p1, p2))
    return tuple(result)


def point_triangle_distance2(
    p: Sequence[float],
    a: Sequence[float],
    b: Sequence[float],
    c: Sequence[float],
) -> float:
    """Squared distance from p to triangle abc (Ericson, Real-Time Collision Detection)."""
    ab = subtract(b, a)
    ac = subtract(c, a)
    ap = subtract(p, a)
    d1 = dot(ab, ap)
    d2 = dot(ac, ap)
    if d1 <= 0.0 and d2 <= 0.0:
        return distance2(p, a)

    bp = subtract(p, b)
    d3 = dot(ab, bp)
    d4 = dot(ac, bp)
    if d3 >= 0.0 and d4 <= d3:
        return distance2(p, b)

    vc = d1 * d4 - d3 * d2
    if vc <= 0.0 and d1 >= 0.0 and d3 <= 0.0:
        v = d1 / (d1 - d3)
        closest = (a[0] + v * ab[0], a[1] + v * ab[1], a[2] + v * ab[2])
        return distance2(p, closest)

    cp = subtract(p, c)
    d5 = dot(ab, cp)
    d6 = dot(ac, cp)
    if d6 >= 0.0 and d5 <= d6:
        return distance2(p, c)

    vb = d5 * d2 - d1 * d6
    if vb <= 0.0 and d2 >= 0.0 and d6 <= 0.0:
        w = d2 / (d2 - d6)
        closest = (a[0] + w * ac[0], a[1] + w * ac[1], a[2] + w * ac[2])
        return distance2(p, closest)

    va = d3 * d6 - d5 * d4
    if va <= 0.0 and (d4 - d3) >= 0.0 and (d5 - d6) >= 0.0:
        w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
        closest = (
            b[0] + w * (c[0] - b[0]),
            b[1] + w * (c[1] - b[1]),
            b[2] + w * (c[2] - b[2]),
        )
        return distance2(p, closest)

    denom = 1.0 / (va + vb + vc)
    v = vb * denom
    w = vc * denom
    closest = (
        a[0] + ab[0] * v + ac[0] * w,
        a[1] + ab[1] * v + ac[1] * w,
        a[2] + ab[2] * v + ac[2] * w,
    )
    return distance2(p, closest)


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


def unit(v: Sequence[float]) -> Tuple[float, float, float]:
    length = math.sqrt(dot(v, v))
    if length <= 0.0:
        return (0.0, 0.0, 0.0)
    return (v[0] / length, v[1] / length, v[2] / length)


def classify_boundary_faces(
    tetrahedra: Sequence[Tuple[int, Tuple[int, int, int, int], int]],
    physical_triangles: Sequence[Tuple[Tuple[int, int, int], int]],
    coords: Dict[int, Tuple[float, float, float]],
) -> List[Tuple[Tuple[int, int, int], int, int]]:
    """Match volume-boundary faces to physical surface triangles.

    The structured ``Alternate`` transfinite meshes can use a different
    diagonal on the geometric surface than on the corresponding volume
    boundary.  Matching the surface triangles by node identity therefore fails
    for half of the faces.  Instead, extract the true exterior faces from the
    tetrahedra and classify each face by its nearest physical triangle
    (point-triangle distance), with the face normal used to break ties at
    surface seams.
    """
    face_owners: Dict[Tuple[int, int, int], List[int]] = defaultdict(list)
    face_nodes: Dict[Tuple[int, int, int], Tuple[int, int, int]] = {}
    for lid, (_eid, conn, _mat) in enumerate(tetrahedra, start=1):
        for face in outward_faces(conn, coords):
            key = tuple(sorted(face))
            face_owners[key].append(lid)
            face_nodes[key] = face

    exterior: List[Tuple[Tuple[int, int, int], int]] = []
    for key, owners in face_owners.items():
        if len(owners) > 2:
            raise RuntimeError(
                f"Non-manifold tet face shared by {len(owners)} elements"
            )
        if len(owners) == 1:
            exterior.append((face_nodes[key], owners[0]))

    if not exterior:
        raise RuntimeError("No exterior tetrahedron faces found")

    triangle_data: List[
        Tuple[Tuple[float, float, float], Tuple[float, float, float], Tuple[float, float, float], int, Tuple[float, float, float]]
    ] = []
    edge_lengths: List[float] = []
    for conn, bc in physical_triangles:
        p0 = coords[conn[0]]
        p1 = coords[conn[1]]
        p2 = coords[conn[2]]
        n = unit(cross(subtract(p1, p0), subtract(p2, p0)))
        triangle_data.append((p0, p1, p2, bc, n))
        edge_lengths.extend(
            (
                math.sqrt(distance2(p0, p1)),
                math.sqrt(distance2(p1, p2)),
                math.sqrt(distance2(p2, p0)),
            )
        )

    edge_lengths.sort()
    median_edge = edge_lengths[len(edge_lengths) // 2] if edge_lengths else 1.0
    all_x = [p[0] for data in triangle_data for p in data[:3]]
    all_y = [p[1] for data in triangle_data for p in data[:3]]
    all_z = [p[2] for data in triangle_data for p in data[:3]]
    diagonal = math.sqrt(
        (max(all_x) - min(all_x)) ** 2
        + (max(all_y) - min(all_y)) ** 2
        + (max(all_z) - min(all_z)) ** 2
    )
    cell = max(median_edge, diagonal * 1.0e-6, 1.0e-12)

    grid: Dict[Tuple[int, int, int], List[int]] = defaultdict(list)
    for ti, (p0, p1, p2, _bc, _n) in enumerate(triangle_data):
        lo = (
            math.floor((min(p0[0], p1[0], p2[0]) - 0.05 * cell) / cell),
            math.floor((min(p0[1], p1[1], p2[1]) - 0.05 * cell) / cell),
            math.floor((min(p0[2], p1[2], p2[2]) - 0.05 * cell) / cell),
        )
        hi = (
            math.floor((max(p0[0], p1[0], p2[0]) + 0.05 * cell) / cell),
            math.floor((max(p0[1], p1[1], p2[1]) + 0.05 * cell) / cell),
            math.floor((max(p0[2], p1[2], p2[2]) + 0.05 * cell) / cell),
        )
        for ix in range(lo[0], hi[0] + 1):
            for iy in range(lo[1], hi[1] + 1):
                for iz in range(lo[2], hi[2] + 1):
                    grid[(ix, iy, iz)].append(ti)

    tolerance2 = (max(1.0e-12, diagonal * 1.0e-8)) ** 2
    records: List[Tuple[Tuple[int, int, int], int, int]] = []
    for face, parent in exterior:
        p0 = coords[face[0]]
        p1 = coords[face[1]]
        p2 = coords[face[2]]
        centroid = (
            (p0[0] + p1[0] + p2[0]) / 3.0,
            (p0[1] + p1[1] + p2[1]) / 3.0,
            (p0[2] + p1[2] + p2[2]) / 3.0,
        )
        face_normal = unit(cross(subtract(p1, p0), subtract(p2, p0)))

        base = (
            math.floor(centroid[0] / cell),
            math.floor(centroid[1] / cell),
            math.floor(centroid[2] / cell),
        )
        candidates: List[int] = []
        for radius in (1, 2):
            candidates = []
            for ix in range(base[0] - radius, base[0] + radius + 1):
                for iy in range(base[1] - radius, base[1] + radius + 1):
                    for iz in range(base[2] - radius, base[2] + radius + 1):
                        candidates.extend(grid.get((ix, iy, iz), ()))
            if candidates:
                break
        if not candidates:
            candidates = list(range(len(triangle_data)))

        best: Optional[Tuple[float, float, int]] = None
        for ti in candidates:
            a, b, c, _bc, normal = triangle_data[ti]
            d2 = point_triangle_distance2(centroid, a, b, c)
            alignment = abs(dot(face_normal, normal))
            score = (d2, -alignment, ti)
            if best is None or score < best:
                best = score
        assert best is not None
        if best[0] > tolerance2:
            raise RuntimeError(
                f"Exterior face {face} is {math.sqrt(best[0]):.3e} from the nearest "
                "physical surface triangle; check the Gmsh physical groups"
            )
        records.append((face, parent, triangle_data[best[2]][3]))

    records.sort(key=lambda row: (row[1], tuple(sorted(row[0]))))
    return records


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

        physical_triangles: List[Tuple[Tuple[int, int, int], int]] = []
        for _dim, entity in gmsh.model.getEntities(2):
            if entity not in surface_physical:
                continue
            ptag = surface_physical[entity]
            name = gmsh.model.getPhysicalName(2, ptag)
            bc = boundary_bc(ptag, name)
            _tags, flat = get_elements(2, entity, TRI3)
            for i in range(0, len(flat), 3):
                conn = tuple(int(v) for v in flat[i : i + 3])
                physical_triangles.append((conn, bc))

        if not tetrahedra:
            raise RuntimeError("No TET4 elements found")
        if not physical_triangles:
            raise RuntimeError("No physical TRI3 boundary faces found")

        tetrahedra.sort(key=lambda row: row[0])
        boundary_records = classify_boundary_faces(
            tetrahedra, physical_triangles, coords
        )

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
