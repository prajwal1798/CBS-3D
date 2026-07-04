#!/usr/bin/env python3
"""
export_cbs3d_tet.py

Gmsh MSH 2.2 ASCII -> Swansea CBS3D tetrahedral input files.

CBS3D .plt written:
    NELEM NPOIN NBOUN
    elem_id n1 n2 n3 n4
    node_id x y z
    face_n1 face_n2 face_n3 parent_tet boundary_side_id

.bco written:
    NFLAG MASS_CHECK
    boundary_side_id solver_bc_flag

Important:
    - No material/region column is written in the .plt element block.
    - No material/region column is written in the .plt coordinate block.
    - Boundary side IDs in .plt are mapped to solver BC IDs in .bco.
    - For cube/cavity meshes, boundary side IDs are geometric face IDs 1..6.

Default cube/cavity face IDs:
    1 : xmin face
    2 : xmax face
    3 : ymin face
    4 : ymax face / moving lid by default
    5 : zmin face
    6 : zmax face

Default LidDrivenCavity3D solver BCs:
    1 -> 530
    2 -> 530
    3 -> 530
    4 -> 500
    5 -> 530
    6 -> 530

Typical LidDrivenCavity3D command:
    python export_cbs3d_tet.py LidDrivenCavity3D.msh -o LidDrivenCavity3D

Override cube-face BCs:
    python export_cbs3d_tet.py case.msh -o case --face-bc 4=500 --face-bc 1=530

Use original Gmsh physical-surface IDs instead of geometric cube face IDs:
    python export_cbs3d_tet.py case.msh -o case --side-id-mode physical --bc 301=530
"""

from __future__ import annotations

import argparse
import html
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple


VALID_SOLVER_BC_IDS = {
    500, 501, 502, 503, 506, 507, 508, 510,
    511, 520, 530, 532, 901, 902,
}

# Raw boundary IDs used for simple cube/cavity exports.
CUBE_FACE_NAMES = {
    1: "xmin",
    2: "xmax",
    3: "ymin",
    4: "ymax",
    5: "zmin",
    6: "zmax",
}

DEFAULT_CUBE_FACE_BC = {
    1: 530,
    2: 530,
    3: 530,
    4: 500,
    5: 530,
    6: 530,
}

DEFAULT_BC_BY_PHYSICAL_ID = {
    # Kept only for --side-id-mode physical.
    301: 530,
    401: 500,
    501: 530,
    502: 530,
}

DEFAULT_BC_BY_NAME = {
    "lid": 500,
    "moving_lid": 500,
    "top_lid": 500,
    "moving_wall": 500,

    "wall": 530,
    "walls": 530,
    "fixed_wall": 530,
    "stationary_wall": 530,
    "no_slip_wall": 530,
    "adiabatic_wall": 530,
    "bottom_wall": 530,
    "side_wall": 530,
    "front_wall": 530,
    "back_wall": 530,

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


@dataclass(frozen=True)
class GmshElement:
    eid: int
    etype: int
    tags: List[int]
    conn: List[int]


@dataclass(frozen=True)
class BoundaryFace:
    nodes: Tuple[int, int, int]
    parent_tet: int
    side_id: int
    gmsh_surface_id: int


def _find_section(lines: List[str], name: str) -> int:
    token = f"${name}"
    try:
        return lines.index(token)
    except ValueError as exc:
        raise ValueError(f"Required MSH2 section {token} not found") from exc


def read_msh2(path: Path) -> Tuple[Dict[Tuple[int, int], str], Dict[int, Tuple[float, float, float]], List[GmshElement]]:
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]

    i_fmt = _find_section(lines, "MeshFormat")
    fmt = lines[i_fmt + 1].split()
    if len(fmt) < 3:
        raise ValueError("Invalid $MeshFormat line")

    version = fmt[0]
    file_type = int(fmt[1])
    if version != "2.2" or file_type != 0:
        raise ValueError(
            f"Unsupported mesh format version={version}, file_type={file_type}. "
            "Regenerate using: gmsh case.geo -3 -format msh2 -o case.msh"
        )

    physical_names: Dict[Tuple[int, int], str] = {}
    if "$PhysicalNames" in lines:
        i_phys = _find_section(lines, "PhysicalNames")
        nphys = int(lines[i_phys + 1])
        for k in range(nphys):
            parts = lines[i_phys + 2 + k].split(maxsplit=2)
            if len(parts) == 3:
                dim = int(parts[0])
                tag = int(parts[1])
                name = parts[2].strip().strip('"')
                physical_names[(dim, tag)] = name

    i_nodes = _find_section(lines, "Nodes")
    nnode = int(lines[i_nodes + 1])
    nodes: Dict[int, Tuple[float, float, float]] = {}
    for k in range(nnode):
        raw = lines[i_nodes + 2 + k]
        parts = raw.split()
        if len(parts) != 4:
            raise ValueError(f"Invalid node line: {raw}")
        nid = int(parts[0])
        nodes[nid] = (float(parts[1]), float(parts[2]), float(parts[3]))

    i_elems = _find_section(lines, "Elements")
    nelem_total = int(lines[i_elems + 1])
    elements: List[GmshElement] = []
    for k in range(nelem_total):
        raw = lines[i_elems + 2 + k]
        parts = [int(x) for x in raw.split()]
        if len(parts) < 4:
            raise ValueError(f"Invalid element line: {raw}")

        eid = parts[0]
        etype = parts[1]
        ntags = parts[2]
        tags = parts[3:3 + ntags]
        conn = parts[3 + ntags:]
        if len(tags) != ntags:
            raise ValueError(f"Tag-count mismatch in element line: {raw}")
        elements.append(GmshElement(eid=eid, etype=etype, tags=tags, conn=conn))

    return physical_names, nodes, elements


def signed_tet_volume6(nodes: Dict[int, Tuple[float, float, float]], conn: Sequence[int]) -> float:
    n1, n2, n3, n4 = conn
    x1, y1, z1 = nodes[n1]
    x2, y2, z2 = nodes[n2]
    x3, y3, z3 = nodes[n3]
    x4, y4, z4 = nodes[n4]

    ax, ay, az = x2 - x1, y2 - y1, z2 - z1
    bx, by, bz = x3 - x1, y3 - y1, z3 - z1
    cx, cy, cz = x4 - x1, y4 - y1, z4 - z1

    return (
        ax * (by * cz - bz * cy)
        - ay * (bx * cz - bz * cx)
        + az * (bx * cy - by * cx)
    )


def orient_tets_positive(
    nodes: Dict[int, Tuple[float, float, float]],
    tets: List[Tuple[int, int, List[int]]],
) -> Tuple[List[Tuple[int, int, List[int]]], int, float | None]:
    fixed: List[Tuple[int, int, List[int]]] = []
    flipped = 0
    min_abs_v6: float | None = None

    for tid, region_id, conn in tets:
        v6 = signed_tet_volume6(nodes, conn)
        min_abs_v6 = abs(v6) if min_abs_v6 is None else min(min_abs_v6, abs(v6))

        if abs(v6) < 1.0e-14:
            raise ValueError(f"Degenerate tetrahedron: cbs_tet={tid}, conn={conn}, vol6={v6}")

        if v6 < 0.0:
            conn = [conn[0], conn[2], conn[1], conn[3]]
            flipped += 1
            v6_after = signed_tet_volume6(nodes, conn)
            if v6_after <= 0.0:
                raise ValueError(f"Failed to orient tetrahedron positively: cbs_tet={tid}")

        fixed.append((tid, region_id, conn))

    return fixed, flipped, min_abs_v6


def tet_face_sets(conn: Sequence[int]) -> List[Tuple[int, int, int]]:
    n1, n2, n3, n4 = conn
    return [
        (n2, n3, n4),
        (n1, n4, n3),
        (n1, n2, n4),
        (n1, n3, n2),
    ]


def bbox_for_nodes(nodes: Dict[int, Tuple[float, float, float]]) -> Tuple[float, float, float, float, float, float]:
    xs = [xyz[0] for xyz in nodes.values()]
    ys = [xyz[1] for xyz in nodes.values()]
    zs = [xyz[2] for xyz in nodes.values()]
    return min(xs), max(xs), min(ys), max(ys), min(zs), max(zs)


def cube_face_tolerance(bbox: Tuple[float, float, float, float, float, float]) -> float:
    xmin, xmax, ymin, ymax, zmin, zmax = bbox
    scale = max(xmax - xmin, ymax - ymin, zmax - zmin, 1.0)
    return 1.0e-8 * scale


def classify_cube_face(
    nodes: Dict[int, Tuple[float, float, float]],
    tri_nodes: Sequence[int],
    bbox: Tuple[float, float, float, float, float, float],
    tol: float,
) -> int:
    xmin, xmax, ymin, ymax, zmin, zmax = bbox
    coords = [nodes[n] for n in tri_nodes]

    candidates: List[int] = []

    if all(abs(x - xmin) <= tol for x, _y, _z in coords):
        candidates.append(1)
    if all(abs(x - xmax) <= tol for x, _y, _z in coords):
        candidates.append(2)
    if all(abs(y - ymin) <= tol for _x, y, _z in coords):
        candidates.append(3)
    if all(abs(y - ymax) <= tol for _x, y, _z in coords):
        candidates.append(4)
    if all(abs(z - zmin) <= tol for _x, _y, z in coords):
        candidates.append(5)
    if all(abs(z - zmax) <= tol for _x, _y, z in coords):
        candidates.append(6)

    if len(candidates) != 1:
        raise ValueError(
            "Could not uniquely classify boundary triangle into cube face 1..6. "
            f"tri_nodes={tri_nodes}, candidates={candidates}, tol={tol}. "
            "For non-cuboid geometry, rerun with --side-id-mode physical."
        )

    return candidates[0]


def build_boundary_faces(
    tets: List[Tuple[int, int, List[int]]],
    triangles: List[Tuple[int, int, int, List[int]]],
    nodes: Dict[int, Tuple[float, float, float]],
    side_id_mode: str,
    physical_to_side_id: Dict[int, int],
) -> List[BoundaryFace]:
    face_to_parent: Dict[Tuple[int, int, int], List[int]] = {}

    for tid, _region_id, conn in tets:
        for face in tet_face_sets(conn):
            key = tuple(sorted(face))
            face_to_parent.setdefault(key, []).append(tid)

    bbox = bbox_for_nodes(nodes)
    tol = cube_face_tolerance(bbox)

    boundary_faces: List[BoundaryFace] = []
    missing = []
    non_unique = []

    for tri_id, physical_id, _geom_id, tri_nodes in triangles:
        key = tuple(sorted(tri_nodes))
        parents = face_to_parent.get(key, [])
        if len(parents) == 0:
            missing.append((tri_id, tri_nodes, physical_id))
            continue
        if len(parents) > 1:
            non_unique.append((tri_id, tri_nodes, parents, physical_id))
            continue

        if side_id_mode == "cube":
            side_id = classify_cube_face(nodes, tri_nodes, bbox, tol)
        else:
            side_id = physical_to_side_id[physical_id]

        boundary_faces.append(
            BoundaryFace(
                nodes=(tri_nodes[0], tri_nodes[1], tri_nodes[2]),
                parent_tet=parents[0],
                side_id=side_id,
                gmsh_surface_id=physical_id,
            )
        )

    if missing:
        tri_id, tri_nodes, physical_id = missing[0]
        raise ValueError(
            f"{len(missing)} boundary triangles did not match any tetra face. "
            f"First unmatched: tri_id={tri_id}, nodes={tri_nodes}, physical={physical_id}"
        )
    if non_unique:
        tri_id, tri_nodes, parents, physical_id = non_unique[0]
        raise ValueError(
            f"{len(non_unique)} boundary triangles matched multiple tetra faces. "
            f"First non-unique: tri_id={tri_id}, nodes={tri_nodes}, parents={parents}, physical={physical_id}"
        )

    return boundary_faces


def parse_bc_overrides(items: List[str]) -> Dict[int, int]:
    out: Dict[int, int] = {}
    for item in items:
        if "=" not in item:
            raise ValueError(f"Invalid --bc item '{item}'. Use physical_surface_id=solver_bc_id")
        left, right = item.split("=", 1)
        physical_id = int(left.strip())
        solver_flag = int(right.strip())
        if solver_flag not in VALID_SOLVER_BC_IDS:
            raise ValueError(
                f"Invalid solver BC ID {solver_flag} in --bc {item}. "
                f"Allowed IDs: {sorted(VALID_SOLVER_BC_IDS)}"
            )
        out[physical_id] = solver_flag
    return out


def parse_face_bc_overrides(items: List[str]) -> Dict[int, int]:
    out: Dict[int, int] = {}
    for item in items:
        if "=" not in item:
            raise ValueError(f"Invalid --face-bc item '{item}'. Use face_id=solver_bc_id")
        left, right = item.split("=", 1)
        face_id = int(left.strip())
        solver_flag = int(right.strip())

        if face_id not in CUBE_FACE_NAMES:
            raise ValueError(f"Invalid cube face ID {face_id}. Valid IDs are 1,2,3,4,5,6.")
        if solver_flag not in VALID_SOLVER_BC_IDS:
            raise ValueError(
                f"Invalid solver BC ID {solver_flag} in --face-bc {item}. "
                f"Allowed IDs: {sorted(VALID_SOLVER_BC_IDS)}"
            )

        out[face_id] = solver_flag

    return out


def _normalise_name(name: str) -> str:
    return name.strip().lower().replace("-", "_").replace(" ", "_")


def solver_flag_for_surface(
    physical_names: Dict[Tuple[int, int], str],
    physical_surface_id: int,
    bc_overrides: Dict[int, int],
) -> int:
    if physical_surface_id in bc_overrides:
        return bc_overrides[physical_surface_id]

    name = physical_names.get((2, physical_surface_id), "")
    key = _normalise_name(name)
    if key in DEFAULT_BC_BY_NAME:
        return DEFAULT_BC_BY_NAME[key]

    if physical_surface_id in DEFAULT_BC_BY_PHYSICAL_ID:
        return DEFAULT_BC_BY_PHYSICAL_ID[physical_surface_id]

    if physical_surface_id in VALID_SOLVER_BC_IDS:
        return physical_surface_id

    raise ValueError(
        f"No CBS3D BC mapping found for physical surface {physical_surface_id} "
        f"name='{name}'. Add an explicit override, e.g. --bc {physical_surface_id}=530"
    )


def write_plt(
    path: Path,
    nodes: Dict[int, Tuple[float, float, float]],
    tets: List[Tuple[int, int, List[int]]],
    boundary_faces: List[BoundaryFace],
) -> None:
    with path.open("w", newline="\n") as f:
        f.write(f"{len(tets):10d}{len(nodes):10d}{len(boundary_faces):10d}\n")

        for tid, _region_id, conn in tets:
            f.write(
                f"{tid:10d}"
                f"{conn[0]:10d}"
                f"{conn[1]:10d}"
                f"{conn[2]:10d}"
                f"{conn[3]:10d}\n"
            )

        for nid in sorted(nodes):
            x, y, z = nodes[nid]
            f.write(
                f"{nid:10d}"
                f"{x: .16e}"
                f"{y: .16e}"
                f"{z: .16e}\n"
            )

        for bf in boundary_faces:
            n1, n2, n3 = bf.nodes
            f.write(
                f"{n1:10d}"
                f"{n2:10d}"
                f"{n3:10d}"
                f"{bf.parent_tet:10d}"
                f"{bf.side_id:10d}\n"
            )


def write_bco(
    path: Path,
    side_ids: List[int],
    side_to_solver_flag: Dict[int, int],
    mass_check: int = 0,
) -> None:
    with path.open("w", newline="\n") as f:
        f.write(f"{len(side_ids):10d}{mass_check:10d}\n")
        for side_id in sorted(side_ids):
            solver_flag = side_to_solver_flag[side_id]
            if solver_flag not in VALID_SOLVER_BC_IDS:
                raise ValueError(
                    f"Internal error: side_id={side_id} has invalid solver BC flag {solver_flag}"
                )
            f.write(f"{side_id:10d}{solver_flag:10d}\n")


def _write_data_array(f, name: str, values: List, vtk_type: str = "Int32", ncomp: int | None = None) -> None:
    comp_attr = f' NumberOfComponents="{ncomp}"' if ncomp else ""
    f.write(f'        <DataArray type="{vtk_type}" Name="{html.escape(name)}" format="ascii"{comp_attr}>\n')
    if values and isinstance(values[0], (tuple, list)):
        for row in values:
            f.write("          " + " ".join(str(v) for v in row) + "\n")
    else:
        for i in range(0, len(values), 12):
            f.write("          " + " ".join(str(v) for v in values[i:i + 12]) + "\n")
    f.write("        </DataArray>\n")


def write_volume_vtu(path: Path, nodes: Dict[int, Tuple[float, float, float]], tets: List[Tuple[int, int, List[int]]]) -> None:
    node_ids = sorted(nodes)
    point_index = {nid: i for i, nid in enumerate(node_ids)}
    points = [nodes[nid] for nid in node_ids]

    connectivity: List[int] = []
    offsets: List[int] = []
    types: List[int] = []
    cbs_tet_id: List[int] = []
    region_id: List[int] = []

    off = 0
    for tid, rid, conn in tets:
        connectivity.extend(point_index[n] for n in conn)
        off += 4
        offsets.append(off)
        types.append(10)
        cbs_tet_id.append(tid)
        region_id.append(rid)

    with path.open("w", newline="\n") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
        f.write('  <UnstructuredGrid>\n')
        f.write(f'    <Piece NumberOfPoints="{len(points)}" NumberOfCells="{len(tets)}">\n')
        f.write('      <Points>\n')
        _write_data_array(f, "Points", points, "Float64", ncomp=3)
        f.write('      </Points>\n')
        f.write('      <Cells>\n')
        _write_data_array(f, "connectivity", connectivity, "Int32")
        _write_data_array(f, "offsets", offsets, "Int32")
        _write_data_array(f, "types", types, "UInt8")
        f.write('      </Cells>\n')
        f.write('      <CellData Scalars="region_id">\n')
        _write_data_array(f, "cbs_tet_id", cbs_tet_id, "Int32")
        _write_data_array(f, "region_id", region_id, "Int32")
        _write_data_array(f, "gmsh_physical_volume", region_id, "Int32")
        f.write('      </CellData>\n')
        f.write('    </Piece>\n')
        f.write('  </UnstructuredGrid>\n')
        f.write('</VTKFile>\n')


def write_boundary_vtu(
    path: Path,
    nodes: Dict[int, Tuple[float, float, float]],
    boundary_faces: List[BoundaryFace],
    side_to_solver_flag: Dict[int, int],
) -> None:
    used = sorted({nid for bf in boundary_faces for nid in bf.nodes})
    point_index = {nid: i for i, nid in enumerate(used)}
    points = [nodes[nid] for nid in used]

    connectivity: List[int] = []
    offsets: List[int] = []
    types: List[int] = []
    side_id: List[int] = []
    solver_flag: List[int] = []
    gmsh_surface: List[int] = []
    parent_tet: List[int] = []

    off = 0
    for bf in boundary_faces:
        connectivity.extend(point_index[n] for n in bf.nodes)
        off += 3
        offsets.append(off)
        types.append(5)
        side_id.append(bf.side_id)
        solver_flag.append(side_to_solver_flag[bf.side_id])
        gmsh_surface.append(bf.gmsh_surface_id)
        parent_tet.append(bf.parent_tet)

    with path.open("w", newline="\n") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
        f.write('  <UnstructuredGrid>\n')
        f.write(f'    <Piece NumberOfPoints="{len(points)}" NumberOfCells="{len(boundary_faces)}">\n')
        f.write('      <Points>\n')
        _write_data_array(f, "Points", points, "Float64", ncomp=3)
        f.write('      </Points>\n')
        f.write('      <Cells>\n')
        _write_data_array(f, "connectivity", connectivity, "Int32")
        _write_data_array(f, "offsets", offsets, "Int32")
        _write_data_array(f, "types", types, "UInt8")
        f.write('      </Cells>\n')
        f.write('      <CellData Scalars="solver_flag">\n')
        _write_data_array(f, "side_id", side_id, "Int32")
        _write_data_array(f, "solver_flag", solver_flag, "Int32")
        _write_data_array(f, "gmsh_physical_surface", gmsh_surface, "Int32")
        _write_data_array(f, "parent_tet", parent_tet, "Int32")
        f.write('      </CellData>\n')
        f.write('    </Piece>\n')
        f.write('  </UnstructuredGrid>\n')
        f.write('</VTKFile>\n')


def write_report(
    path: Path,
    physical_names: Dict[Tuple[int, int], str],
    node_count: int,
    n_triangles: int,
    n_tets: int,
    n_boundary: int,
    flipped: int,
    min_abs_v6: float | None,
    side_count: Dict[int, int],
    surface_count: Dict[int, int],
    region_count: Dict[int, int],
    side_to_solver_flag: Dict[int, int],
    vtu_path: Path,
    boundary_vtu_path: Path,
) -> None:
    with path.open("w", newline="\n") as f:
        f.write("CBS3D Tetra Export Report\n")
        f.write("=========================\n\n")
        f.write(f"CBS nodes written      : {node_count}\n")
        f.write(f"Gmsh surface triangles : {n_triangles}\n")
        f.write(f"CBS tetrahedra written : {n_tets}\n")
        f.write(f"Boundary faces written : {n_boundary}\n")
        f.write(f"Tetrahedra flipped     : {flipped}\n")
        f.write(f"Minimum |6*volume|     : {min_abs_v6}\n")
        f.write(f"Volume VTU             : {vtu_path}\n")
        f.write(f"Boundary VTU           : {boundary_vtu_path}\n\n")

        f.write("Physical names\n")
        f.write("--------------\n")
        for (dim, tag), name in sorted(physical_names.items()):
            f.write(f"dim={dim:1d} tag={tag:8d} name={name}\n")

        f.write("\nCBS boundary side IDs\n")
        f.write("---------------------\n")
        for side_id in sorted(side_to_solver_flag):
            f.write(
                f"side_id={side_id:2d} "
                f"name={CUBE_FACE_NAMES.get(side_id, 'physical')} "
                f"solver_flag={side_to_solver_flag[side_id]:8d} "
                f"faces={side_count.get(side_id, 0):10d}\n"
            )

        f.write("\nGmsh physical surface counts\n")
        f.write("----------------------------\n")
        for phys in sorted(surface_count):
            name = physical_names.get((2, phys), "")
            f.write(f"physical_surface={phys:8d} faces={surface_count[phys]:10d} name={name}\n")

        f.write("\nVolume region counts\n")
        f.write("--------------------\n")
        for rid in sorted(region_count):
            name = physical_names.get((3, rid), "")
            f.write(f"physical_volume={rid:8d} tets={region_count[rid]:10d} name={name}\n")

        f.write("\nParaView audit\n")
        f.write("--------------\n")
        f.write("Open the boundary VTU and color by solver_flag or side_id.\n")
        f.write("Open the volume VTU and color by region_id.\n")


def export_cbs3d(
    msh_path: Path,
    out_base: Path,
    bc_overrides: Dict[int, int],
    face_bc_overrides: Dict[int, int],
    side_id_mode: str,
) -> Tuple[Path, Path, Path, Path, Path]:
    physical_names, gmsh_nodes, elements = read_msh2(msh_path)

    raw_triangles = [e for e in elements if e.etype == 2]
    raw_tets = [e for e in elements if e.etype == 4]

    if not gmsh_nodes:
        raise ValueError("No nodes found in .msh file")
    if not raw_triangles:
        raise ValueError("No boundary triangles found. Define Physical Surface groups in Gmsh.")
    if not raw_tets:
        raise ValueError("No tetrahedra found. Generate a 3D tetra mesh using gmsh case.geo -3 -format msh2.")

    unsupported = sorted({e.etype for e in elements if e.etype not in {1, 2, 4, 15}})
    if unsupported:
        raise ValueError(
            f"Unsupported Gmsh element types found: {unsupported}. "
            "This exporter accepts only linear line/triangle/tetra/point elements. Regenerate with -order 1."
        )

    used_gmsh_nodes = set()
    for e in raw_triangles + raw_tets:
        used_gmsh_nodes.update(e.conn)

    missing = sorted(n for n in used_gmsh_nodes if n not in gmsh_nodes)
    if missing:
        raise ValueError(f"Elements reference missing node IDs: {missing[:10]}")

    node_map = {old: new for new, old in enumerate(sorted(used_gmsh_nodes), start=1)}
    cbs_nodes = {node_map[old]: gmsh_nodes[old] for old in sorted(used_gmsh_nodes)}

    tets: List[Tuple[int, int, List[int]]] = []
    for new_tid, e in enumerate(raw_tets, start=1):
        if len(e.conn) != 4:
            raise ValueError(f"Tetra element {e.eid} does not have 4 nodes: {e.conn}")
        if not e.tags:
            raise ValueError(f"Tetra element {e.eid} has no physical volume tag")
        region_id = e.tags[0]
        conn = [node_map[n] for n in e.conn]
        tets.append((new_tid, region_id, conn))

    triangles: List[Tuple[int, int, int, List[int]]] = []
    for e in raw_triangles:
        if len(e.conn) != 3:
            raise ValueError(f"Triangle element {e.eid} does not have 3 nodes: {e.conn}")
        if not e.tags:
            raise ValueError(f"Triangle element {e.eid} has no physical surface tag")
        physical_surface_id = e.tags[0]
        geom_surface_id = e.tags[1] if len(e.tags) > 1 else -1
        conn = [node_map[n] for n in e.conn]
        triangles.append((e.eid, physical_surface_id, geom_surface_id, conn))

    surface_ids = sorted({phys for _eid, phys, _geom, _conn in triangles})
    physical_to_side_id = {phys: phys for phys in surface_ids}

    if side_id_mode == "cube":
        side_to_solver_flag = dict(DEFAULT_CUBE_FACE_BC)
        side_to_solver_flag.update(face_bc_overrides)
    else:
        side_to_solver_flag = {}
        for phys in surface_ids:
            side_id = physical_to_side_id[phys]
            side_to_solver_flag[side_id] = solver_flag_for_surface(physical_names, phys, bc_overrides)

    tets, flipped, min_abs_v6 = orient_tets_positive(cbs_nodes, tets)
    boundary_faces = build_boundary_faces(tets, triangles, cbs_nodes, side_id_mode, physical_to_side_id)

    side_count: Dict[int, int] = {}
    surface_count: Dict[int, int] = {}
    for bf in boundary_faces:
        side_count[bf.side_id] = side_count.get(bf.side_id, 0) + 1
        surface_count[bf.gmsh_surface_id] = surface_count.get(bf.gmsh_surface_id, 0) + 1

    region_count: Dict[int, int] = {}
    for _tid, region_id, _conn in tets:
        region_count[region_id] = region_count.get(region_id, 0) + 1

    if side_id_mode == "cube":
        side_ids = [1, 2, 3, 4, 5, 6]
        missing_faces = [sid for sid in side_ids if side_count.get(sid, 0) == 0]
        if missing_faces:
            raise ValueError(
                f"Cube face mode expected boundary faces on all six sides, but these side IDs are empty: {missing_faces}"
            )
    else:
        side_ids = sorted(side_count)

    plt_path = out_base.with_suffix(".plt")
    bco_path = out_base.with_suffix(".bco")
    volume_vtu_path = out_base.with_suffix(".vtu")
    boundary_vtu_path = out_base.with_name(out_base.name + "_boundary").with_suffix(".vtu")
    report_path = out_base.with_suffix(".report.txt")

    write_plt(plt_path, cbs_nodes, tets, boundary_faces)
    write_bco(bco_path, side_ids, side_to_solver_flag)
    write_volume_vtu(volume_vtu_path, cbs_nodes, tets)
    write_boundary_vtu(boundary_vtu_path, cbs_nodes, boundary_faces, side_to_solver_flag)
    write_report(
        report_path,
        physical_names,
        len(cbs_nodes),
        len(triangles),
        len(tets),
        len(boundary_faces),
        flipped,
        min_abs_v6,
        side_count,
        surface_count,
        region_count,
        side_to_solver_flag,
        volume_vtu_path,
        boundary_vtu_path,
    )

    return plt_path, bco_path, volume_vtu_path, boundary_vtu_path, report_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert Gmsh MSH2 3D tetra mesh into CBS3D .plt/.bco plus VTU audit files."
    )
    parser.add_argument("msh", type=Path, help="Input Gmsh MSH2 ASCII file")
    parser.add_argument("-o", "--out-base", type=Path, default=None, help="Output basename. Default: input mesh path without suffix.")
    parser.add_argument(
        "--side-id-mode",
        choices=("cube", "physical"),
        default="cube",
        help="cube: write raw boundary IDs as geometric cube faces 1..6. physical: use Gmsh physical surface IDs.",
    )
    parser.add_argument("--bc", action="append", default=[], help="Physical-mode mapping physical_surface_id=solver_bc_id. Can be repeated.")
    parser.add_argument("--face-bc", action="append", default=[], help="Cube-mode mapping face_id=solver_bc_id. Can be repeated.")
    args = parser.parse_args()

    msh_path = args.msh
    out_base = args.out_base if args.out_base is not None else msh_path.with_suffix("")
    bc_overrides = parse_bc_overrides(args.bc)
    face_bc_overrides = parse_face_bc_overrides(args.face_bc)

    plt, bco, vtu, boundary_vtu, report = export_cbs3d(
        msh_path,
        out_base,
        bc_overrides,
        face_bc_overrides,
        args.side_id_mode,
    )

    print("CBS3D tetra export completed successfully.")
    print(f"PLT          : {plt}")
    print(f"BCO          : {bco}")
    print(f"Volume VTU   : {vtu}")
    print(f"Boundary VTU : {boundary_vtu}")
    print(f"REPORT       : {report}")


if __name__ == "__main__":
    main()
