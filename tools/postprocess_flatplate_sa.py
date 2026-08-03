#!/usr/bin/env python3
"""
Post-process a CBS3D++_SI Spalart-Allmaras flat-plate run.

The script intentionally avoids ParaView/PyVista dependencies.  It reads the
ASCII VTU files written by Post.cpp and extracts the quantities needed for the
first verification workflow:

    residual history
    wall-boundary skin-friction coefficient Cf(x)
    u/U_inf profiles at selected x-stations
    mu_t/mu profiles at selected x-stations
    u+ versus y+ profiles at selected x-stations

The wall shear extraction is designed for the current zero-pressure-gradient
flat-plate case where the plate is BC 530 and the wall-normal direction is y.
For each wall boundary triangle, the script uses the parent tetrahedron and the
P1 velocity gradient inside that tetrahedron.  This is a code-verification
post-processing estimate; serious comparison still requires wall-normal mesh
refinement and grid-convergence checks.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

try:
    import numpy as np
except ImportError as exc:  # pragma: no cover
    raise SystemExit("This script requires numpy. Install it with: python -m pip install numpy") from exc

try:
    import matplotlib.pyplot as plt
except ImportError:  # pragma: no cover
    plt = None


@dataclass
class VtuData:
    points: np.ndarray
    cells: List[np.ndarray]
    cell_types: np.ndarray
    point_data: Dict[str, np.ndarray]
    cell_data: Dict[str, np.ndarray]


def numbers(text: Optional[str], dtype=float) -> np.ndarray:
    if text is None:
        return np.array([], dtype=dtype)
    return np.fromstring(text, sep=" ", dtype=dtype)


def read_vtu(path: Path) -> VtuData:
    tree = ET.parse(path)
    root = tree.getroot()
    piece = root.find(".//Piece")
    if piece is None:
        raise RuntimeError(f"No VTK Piece found in {path}")

    pts_array = piece.find("./Points/DataArray")
    if pts_array is None:
        raise RuntimeError(f"No point coordinate array found in {path}")
    pts = numbers(pts_array.text, float).reshape((-1, 3))

    cell_arrays = {a.attrib.get("Name", ""): a for a in piece.findall("./Cells/DataArray")}
    conn = numbers(cell_arrays["connectivity"].text, int)
    offsets = numbers(cell_arrays["offsets"].text, int)
    cell_types = numbers(cell_arrays["types"].text, int)

    cells: List[np.ndarray] = []
    start = 0
    for off in offsets:
        cells.append(conn[start:off].astype(int))
        start = int(off)

    point_data: Dict[str, np.ndarray] = {}
    for arr in piece.findall("./PointData/DataArray"):
        name = arr.attrib.get("Name")
        if not name:
            continue
        ncomp = int(arr.attrib.get("NumberOfComponents", "1"))
        vals = numbers(arr.text, float)
        if ncomp > 1:
            vals = vals.reshape((-1, ncomp))
        point_data[name] = vals

    cell_data: Dict[str, np.ndarray] = {}
    for arr in piece.findall("./CellData/DataArray"):
        name = arr.attrib.get("Name")
        if not name:
            continue
        dtype = int if arr.attrib.get("type", "").lower().startswith(("int", "uint")) else float
        cell_data[name] = numbers(arr.text, dtype)

    return VtuData(pts, cells, cell_types, point_data, cell_data)


def step_number(path: Path) -> int:
    match = re.search(r"_(\d{6,})\.vtu$", path.name)
    if not match:
        return -1
    return int(match.group(1))


def latest_vtu(output_dir: Path) -> Path:
    files = sorted(output_dir.rglob("*.vtu"), key=lambda p: (step_number(p), p.stat().st_mtime))
    if not files:
        raise SystemExit(f"No VTU files found under {output_dir}")
    return files[-1]


def latest_residual_csv(output_dir: Path) -> Optional[Path]:
    files = sorted(output_dir.rglob("*_residuals.csv"), key=lambda p: p.stat().st_mtime)
    return files[-1] if files else None


def tetra_gradient(coords: np.ndarray, values: np.ndarray) -> np.ndarray:
    # P1 tetrahedron: q(x,y,z) = a0 + a1*x + a2*y + a3*z.
    a = np.ones((4, 4), dtype=float)
    a[:, 1:4] = coords
    coeff = np.linalg.solve(a, values)
    return coeff[1:4]


def triangle_area_and_normal(coords: np.ndarray) -> Tuple[float, np.ndarray]:
    cross = np.cross(coords[1] - coords[0], coords[2] - coords[0])
    norm = float(np.linalg.norm(cross))
    if norm <= 0.0 or not math.isfinite(norm):
        return 0.0, np.zeros(3)
    return 0.5 * norm, cross / norm


def write_csv(path: Path, header: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def extract_cf(
    data: VtuData,
    rho: float,
    mu: float,
    u_inf: float,
    wall_bc: int,
    bins: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    velocity = data.point_data.get("velocity")
    if velocity is None:
        raise RuntimeError("VTU does not contain point field 'velocity'")

    bc_id = data.cell_data.get("bc_id")
    parent = data.cell_data.get("parent_element")
    if bc_id is None or parent is None:
        raise RuntimeError("VTU must contain cell fields 'bc_id' and 'parent_element'")

    xs: List[float] = []
    cfs: List[float] = []
    areas: List[float] = []

    for icell, ctype in enumerate(data.cell_types):
        if int(ctype) != 5:
            continue
        if int(bc_id[icell]) != wall_bc:
            continue

        parent_e = int(parent[icell]) - 1
        if parent_e < 0 or parent_e >= len(data.cells):
            continue

        tet = data.cells[parent_e]
        if len(tet) != 4:
            continue

        face = data.cells[icell]
        face_coords = data.points[face]
        area, normal = triangle_area_and_normal(face_coords)
        if area <= 0.0:
            continue

        tet_coords = data.points[tet]
        ux = velocity[tet, 0]
        grad_ux = tetra_gradient(tet_coords, ux)

        # For the flat plate, tau_w = mu * du/dn.  The face normal orientation is
        # not important for the magnitude, hence abs().  At the SA wall, mu_t=0,
        # so molecular mu is the correct wall viscosity.
        tau_w = mu * abs(float(np.dot(grad_ux, normal)))
        cf = 2.0 * tau_w / (rho * u_inf * u_inf)

        xs.append(float(np.mean(face_coords[:, 0])))
        cfs.append(cf)
        areas.append(area)

    if not xs:
        raise RuntimeError(f"No wall triangles found with bc_id={wall_bc}")

    x = np.asarray(xs)
    cf = np.asarray(cfs)
    area = np.asarray(areas)

    xmin = float(np.min(x))
    xmax = float(np.max(x))
    edges = np.linspace(xmin, xmax, bins + 1)

    xb: List[float] = []
    cfb: List[float] = []
    taub: List[float] = []

    for i in range(bins):
        if i == bins - 1:
            mask = (x >= edges[i]) & (x <= edges[i + 1])
        else:
            mask = (x >= edges[i]) & (x < edges[i + 1])
        if not np.any(mask):
            continue
        w = area[mask]
        cf_mean = float(np.sum(cf[mask] * w) / np.sum(w))
        xb.append(0.5 * (edges[i] + edges[i + 1]))
        cfb.append(cf_mean)
        taub.append(0.5 * rho * u_inf * u_inf * cf_mean)

    return np.asarray(xb), np.asarray(cfb), np.asarray(taub)


def station_profile(
    data: VtuData,
    x0: float,
    dx: float,
    z_window: Optional[float],
    z_mid: Optional[float],
) -> np.ndarray:
    pts = data.points
    mask = np.abs(pts[:, 0] - x0) <= dx

    if z_window is not None:
        if z_mid is None:
            z_mid = 0.5 * (float(np.min(pts[:, 2])) + float(np.max(pts[:, 2])))
        mask &= np.abs(pts[:, 2] - z_mid) <= z_window

    idx = np.where(mask)[0]
    if idx.size == 0:
        return np.empty((0, 0))

    velocity = data.point_data["velocity"]
    velocity_magnitude = data.point_data.get("velocity_magnitude")
    nu_tilde = data.point_data.get("nu_tilde")
    mu_t_over_mu = data.point_data.get("mu_t_over_mu")
    wall_distance = data.point_data.get("wall_distance")

    rows = []
    for ip in idx:
        rows.append([
            float(pts[ip, 0]),
            float(pts[ip, 1]),
            float(pts[ip, 2]),
            float(velocity[ip, 0]),
            float(velocity[ip, 1]),
            float(velocity[ip, 2]),
            float(velocity_magnitude[ip]) if velocity_magnitude is not None else math.nan,
            float(nu_tilde[ip]) if nu_tilde is not None else math.nan,
            float(mu_t_over_mu[ip]) if mu_t_over_mu is not None else math.nan,
            float(wall_distance[ip]) if wall_distance is not None else math.nan,
        ])

    arr = np.asarray(rows, dtype=float)
    arr = arr[np.argsort(arr[:, 1])]
    return arr


def interpolate_tau(x_cf: np.ndarray, tau: np.ndarray, x0: float) -> float:
    if x_cf.size == 0:
        return math.nan
    return float(np.interp(x0, x_cf, tau))


def save_plots(post_dir: Path, residual_csv: Optional[Path], cf_table: Optional[np.ndarray], profiles: Dict[float, np.ndarray], rho: float, mu: float, u_inf: float) -> None:
    if plt is None:
        print("[warn] matplotlib not available; CSV files were written but PNG plots were skipped")
        return

    post_dir.mkdir(parents=True, exist_ok=True)

    if residual_csv is not None and residual_csv.exists():
        with residual_csv.open("r", newline="") as f:
            rows = list(csv.DictReader(f))
        if rows:
            it = np.asarray([float(r["iteration"]) for r in rows])
            plt.figure(figsize=(10.5, 6.5))
            for key, label in [
                ("u_rel", "u"), ("v_rel", "v"), ("w_rel", "w"),
                ("p_rel", "p"), ("velocity_rel_max", "max velocity"),
                ("cg_relative_l2", "CG rel L2"),
            ]:
                if key in rows[0]:
                    y = np.asarray([max(float(r[key]), 1.0e-300) for r in rows])
                    plt.semilogy(it, y, label=label)
            plt.xlabel("Iteration")
            plt.ylabel("Residual")
            plt.grid(True, which="both")
            plt.legend()
            plt.tight_layout()
            plt.savefig(post_dir / "residual_history.png", dpi=200)
            plt.close()

            if "cg_iterations" in rows[0]:
                cg = np.asarray([float(r["cg_iterations"]) for r in rows])
                plt.figure(figsize=(10.5, 5.0))
                plt.plot(it, cg)
                plt.xlabel("Iteration")
                plt.ylabel("Pressure CG iterations")
                plt.grid(True)
                plt.tight_layout()
                plt.savefig(post_dir / "cg_iterations.png", dpi=200)
                plt.close()

    if cf_table is not None and cf_table.size:
        plt.figure(figsize=(10.5, 5.0))
        plt.plot(cf_table[:, 0], cf_table[:, 1], marker="o", markersize=3)
        plt.xlabel("x")
        plt.ylabel("C_f")
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(post_dir / "flatplate_cf.png", dpi=200)
        plt.close()

    for x0, arr in profiles.items():
        if arr.size == 0:
            continue
        tag = f"x{x0:.3f}".replace(".", "p")
        y = arr[:, 1]
        u = arr[:, 3]
        mut = arr[:, 8]

        plt.figure(figsize=(6.0, 6.0))
        plt.plot(u / u_inf, y, marker="o", markersize=2)
        plt.xlabel("u / U_inf")
        plt.ylabel("y")
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(post_dir / f"profile_u_{tag}.png", dpi=200)
        plt.close()

        if np.isfinite(mut).any():
            plt.figure(figsize=(6.0, 6.0))
            plt.plot(mut, y, marker="o", markersize=2)
            plt.xlabel("mu_t / mu")
            plt.ylabel("y")
            plt.grid(True)
            plt.tight_layout()
            plt.savefig(post_dir / f"profile_mut_over_mu_{tag}.png", dpi=200)
            plt.close()


def main() -> None:
    parser = argparse.ArgumentParser(description="Post-process CBS3D++_SI SA flat-plate output")
    parser.add_argument("output_dir", type=Path, help="Solver output directory containing VTU/PVD/CSV files")
    parser.add_argument("--rho", type=float, default=1.0)
    parser.add_argument("--mu", type=float, default=1.0e-6)
    parser.add_argument("--u-inf", type=float, default=1.0)
    parser.add_argument("--wall-bc", type=int, default=530)
    parser.add_argument("--stations", type=float, nargs="*", default=[0.97, 1.90])
    parser.add_argument("--station-dx", type=float, default=0.01)
    parser.add_argument("--z-window", type=float, default=None)
    parser.add_argument("--cf-bins", type=int, default=120)
    parser.add_argument("--vtu", type=Path, default=None, help="Specific VTU file. Default: latest VTU under output_dir")
    args = parser.parse_args()

    out_dir = args.output_dir
    post_dir = out_dir / "post_flatplate_sa"
    post_dir.mkdir(parents=True, exist_ok=True)

    vtu = args.vtu if args.vtu is not None else latest_vtu(out_dir)
    print(f"[info] reading {vtu}")
    data = read_vtu(vtu)

    required = ["velocity"]
    missing = [name for name in required if name not in data.point_data]
    if missing:
        raise SystemExit(f"VTU is missing required point fields: {missing}")

    if "mu_t_over_mu" not in data.point_data:
        print("[warn] VTU does not contain mu_t_over_mu. Run tools/apply_sa_postprocessing_support.py and rebuild.")
    if "wall_distance" not in data.point_data:
        print("[warn] VTU does not contain wall_distance. Run tools/apply_sa_postprocessing_support.py and rebuild.")

    x_cf, cf, tau = extract_cf(data, args.rho, args.mu, args.u_inf, args.wall_bc, args.cf_bins)
    cf_table = np.column_stack([x_cf, cf, tau])
    write_csv(post_dir / "flatplate_cf.csv", ["x", "Cf", "tau_w"], cf_table)
    print(f"[ok] wrote {post_dir / 'flatplate_cf.csv'}")

    profile_tables: Dict[float, np.ndarray] = {}
    profile_header = [
        "x", "y", "z", "u", "v", "w", "velocity_magnitude",
        "nu_tilde", "mu_t_over_mu", "wall_distance", "y_plus", "u_plus",
    ]

    for x0 in args.stations:
        arr = station_profile(data, x0, args.station_dx, args.z_window, None)
        if arr.size == 0:
            print(f"[warn] no profile nodes found near x={x0}")
            continue

        tau_w = interpolate_tau(x_cf, tau, x0)
        if math.isfinite(tau_w) and tau_w > 0.0:
            u_tau = math.sqrt(tau_w / args.rho)
            y_plus = arr[:, 1] * u_tau / (args.mu / args.rho)
            u_plus = arr[:, 3] / u_tau
        else:
            y_plus = np.full(arr.shape[0], math.nan)
            u_plus = np.full(arr.shape[0], math.nan)

        table = np.column_stack([arr, y_plus, u_plus])
        profile_tables[x0] = arr
        tag = f"x{x0:.3f}".replace(".", "p")
        write_csv(post_dir / f"profile_{tag}.csv", profile_header, table)
        print(f"[ok] wrote {post_dir / f'profile_{tag}.csv'}")

    residual_csv = latest_residual_csv(out_dir)
    save_plots(post_dir, residual_csv, cf_table, profile_tables, args.rho, args.mu, args.u_inf)
    print(f"[ok] post-processing complete: {post_dir}")


if __name__ == "__main__":
    main()
