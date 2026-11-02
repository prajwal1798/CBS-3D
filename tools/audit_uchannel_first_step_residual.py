#!/usr/bin/env python3
"""G6 read-only reconstruction of the FIRST U-channel AFC-off temperature step.

Uses independently evaluated P1 tetrahedral integrals (not the solver's
assembled RHS), the saved previous T field (step 50000), and the corrected
velocity field (step 50001). The latter is the velocity used by production
Step-4 (the restart velocity is not necessarily that velocity).

For an interior fluid node i, verify:
  T1_i - T0_i = (R_adv + R_div + R_CBS + R_diff + R_source)/H_i,
  H_i = sum_(incident e) (rhoCp_e * V_e / (4 * dt_e)).
This original U-channel run uses a communicator-uniform dt_e. BC532 contributes
nothing at interior node 5707. Original production code is NOT modified.

All inputs are existing original/replay VTU and .par/.matprop files.
"""
import argparse
import csv
import math
from pathlib import Path

from audit_temperature_undershoot_vtu import read_piece, files_for_step
from audit_cbs_explicit_timestep_vtu import (
    gradients_and_volume, read_matprop, dot,
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def param_after_comment(par, comment_prefix, index):
    lines = par.read_text(encoding="utf-8").splitlines()
    for i, line in enumerate(lines[:-1]):
        if line.strip().lower().startswith(comment_prefix.lower()):
            return float(lines[i + 1].split()[index])
    raise ValueError(f"{par}: not found: {comment_prefix}")


def find_timestep(output_dir, step):
    files = list(output_dir.glob("*distributed_residuals.csv"))
    require(len(files) == 1, f"Expected one distributed residual CSV, found {files}")
    with files[0].open(newline="") as stream:
        rows = csv.reader(stream)
        header = next(rows)
        for row in rows:
            if row and int(row[0]) == step:
                dt = float(row[2])  # iteration, pseudo-time, dt, ...
                require(math.isfinite(dt) and dt > 0, "invalid residual CSV dt")
                return dt, header[:3]
    raise ValueError(f"No global iteration {step} in {files[0]}")


def read_global_fields(output_dir, step, with_topology):
    files = files_for_step(output_dir, step)
    require(len(files) == 2, "G6 original replay should have exactly two MPI VTU pieces")
    global_nodes = {}
    pieces = []
    owners = {}
    for file in files:
        piece = read_piece(file, with_topology=with_topology)
        pieces.append((file, piece))
        for i, gid in enumerate(piece["gids"]):
            record = (piece["points"][i], piece["temperature"][i],
                      piece["velocity"][i])
            previous = global_nodes.get(gid)
            if previous is not None:
                old_coords, old_temp, old_vel = previous
                require(max(abs(a - b) for a, b in zip(old_coords, record[0])) <= 1e-10,
                        f"conflicting global-node coordinates {gid}")
                require(abs(old_temp - record[1]) <= 2e-8,
                        f"inconsistent owner/ghost temperature at {step}, gid={gid}")
                require(max(abs(a - b) for a, b in zip(old_vel, record[2])) <= 1e-9,
                        f"inconsistent owner/ghost velocity at {step}, gid={gid}")
            else:
                global_nodes[gid] = record
            if piece["owned"][i] == 1:
                require(gid not in owners, f"duplicate owner node {gid} at step {step}")
                owners[gid] = record
    require(set(global_nodes) == set(owners), "some nodes absent from owned-node records")
    return owners, pieces


def fluid_source(matprop, material):
    tokens = matprop.read_text(encoding="utf-8").split()
    n = int(tokens[0])
    for i in range(n):
        base = 1 + i * 8
        if int(tokens[base]) == material:
            return float(tokens[base + 7])
    raise ValueError(f"Material ID {material} absent")


def node_residual(args):
    out = args.output_dir
    step0, step1 = args.before, args.after
    t0, _ = read_global_fields(out, step0, with_topology=False)
    t1, pieces = read_global_fields(out, step1, with_topology=True)
    require(set(t0) == set(t1), "global node IDs changed across replay")
    require(args.node in t0, f"global node {args.node} not found")
    dt, csv_columns = find_timestep(out, step1)
    inlet_T = param_after_comment(args.par, "reference and inlet/outlet data", 4)
    qpar = param_after_comment(args.par, "alpha_sf k_ratio source_solid", 2)
    require(abs(inlet_T - args.reference_temperature) < 1e-8,
            f"unexpected inlet temperature {inlet_T}")
    require(qpar == 0.0, "G6 target audit assumes no solid source, as original U-channel")
    materials = read_matprop(args.matprop)

    components = {"advection": 0.0, "divergence": 0.0,
                  "characteristic": 0.0, "diffusion": 0.0,
                  "volumetric": 0.0}
    mass_over_dt = 0.0
    cells = []
    seen_tets = set()
    dvmax = 0.0

    for file, piece in pieces:
        start = 0
        for icell, end in enumerate(piece["offsets"]):
            local_nodes = piece["connectivity"][start:end]
            start = end
            if piece["types"][icell] != 10:
                continue
            if args.node not in [piece["gids"][i] for i in local_nodes]:
                continue
            require(len(local_nodes) == 4, "non-TET4 cell incident on target")
            ids = tuple(piece["gids"][i] for i in local_nodes)
            key = tuple(sorted(ids))
            require(key not in seen_tets,
                    f"duplicate tetrahedron attached to global node {args.node}")
            seen_tets.add(key)

            material_id = piece["material_id"][icell]
            require(material_id == 0, f"target touches nonfluid material {material_id}")
            require(material_id in materials, f"unrecognised material {material_id}")
            mat = materials[material_id]
            qvol = fluid_source(args.matprop, material_id)
            points = [piece["points"][i] for i in local_nodes]
            grad, volume = gradients_and_volume(*points)
            rhoCp = mat["rho_cp"]
            kappa = mat["k"]
            temperatures = [t0[gid][1] for gid in ids]
            velocities = [t1[gid][2] for gid in ids]
            ga = ids.index(args.node)

            gradT = tuple(sum(temperatures[a] * grad[a][d] for a in range(4))
                          for d in range(3))
            u_sum = tuple(sum(velocities[a][d] for a in range(4))
                          for d in range(3))
            u_bar = tuple(val / 4.0 for val in u_sum)
            divu = sum(dot(velocities[a], grad[a]) for a in range(4))
            tref_sum = sum(T - inlet_T for T in temperatures)
            u_weight = tuple(u_sum[d] + velocities[ga][d] for d in range(3))
            mass_over_dt += rhoCp * volume / (4.0 * dt)

            adv = -rhoCp * volume / 20.0 * dot(u_weight, gradT)
            div = -rhoCp * volume / 20.0 * divu * (
                tref_sum + (temperatures[ga] - inlet_T))
            char = -0.5 * dt * rhoCp * volume * dot(u_bar, grad[ga]) * dot(u_bar, gradT)
            diff = -kappa * volume * dot(grad[ga], gradT)
            vol = qvol * volume / 4.0
            for label, value in (("advection", adv), ("divergence", div),
                                 ("characteristic", char), ("diffusion", diff),
                                 ("volumetric", vol)):
                components[label] += value

            dvmax = max(dvmax, abs(divu))
            cells.append((ids, volume, divu, adv, div, char, diff, vol))
    require(cells, f"No incident tetrahedra for node {args.node}")
    require(mass_over_dt > 0 and math.isfinite(mass_over_dt),
            "invalid thermal capacity/time diagonal")

    T_before = t0[args.node][1]
    T_after = t1[args.node][1]
    delta_observed = T_after - T_before
    delta_terms = {key: value / mass_over_dt for key, value in components.items()}
    delta_reconstructed = math.fsum(delta_terms.values())
    error = delta_observed - delta_reconstructed

    print("=== G6 ORIGINAL U-CHANNEL FIRST-STEP THERMAL RESIDUAL ===")
    print(f"Global node ID                 : {args.node}")
    print(f"Coordinates [m]               : {t1[args.node][0]}")
    print(f"Restart step / updated step    : {step0} / {step1}")
    print(f"dt (residual CSV) [s]         : {dt:.16e}")
    print(f"Reference temperature [K]     : {inlet_T:.12f}")
    print(f"T0 / T1 [K]                  : {T_before:.12f} / {T_after:.12f}")
    print(f"Observed delta T [K]          : {delta_observed:+.12e}")
    print(f"Incident fluid tetrahedra     : {len(cells)}")
    print(f"Thermal diagonal H [W/K]      : {mass_over_dt:.12e}")
    print(f"Max element |div(u)| [1/s]    : {dvmax:.12e}")
    print("Components of first temperature increment:")
    for key, value in delta_terms.items():
        print(f"  {key:17s} : {value:+.12e} K")
    print(f"Reconstructed delta T [K]     : {delta_reconstructed:+.12e}")
    print(f"Observed - reconstructed [K]  : {error:+.12e}")
    print("Top 8 incident cells by absolute total RHS:")
    for ids, V, div, adv, corr, char, diff, src in sorted(
        cells, key=lambda item: abs(math.fsum(item[3:])), reverse=True
    )[:8]:
        print(f"  nodes={ids} vol={V:.4e} div_u={div:+.3e}"
              f" adv={adv:+.4e} div={corr:+.4e} CBS={char:+.4e}"
              f" diff={diff:+.4e} source={src:+.4e} [W]")
    tolerance = args.tolerance
    require(abs(error) <= tolerance,
            f"FAIL: reconstruction mismatch {error:.4e} K > {tolerance:.4e} K. "
            "Do not attribute individual term mechanisms yet; verify "
            "postprocess velocity, source, conductivity and timestep configuration.")
    print("G6 NODE RESIDUAL RECONSTRUCTION PASS")
    if delta_observed < -args.tolerance:
        print("G6 OBSERVED FIRST-STEP NEGATIVE UPDATE CONFIRMED")
    else:
        print("G6 node update is not negative at requested node")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--par", required=True, type=Path)
    parser.add_argument("--matprop", required=True, type=Path)
    parser.add_argument("--node", type=int, default=5707)
    parser.add_argument("--before", type=int, default=50000)
    parser.add_argument("--after", type=int, default=50001)
    parser.add_argument("--reference-temperature", type=float, default=298.15)
    parser.add_argument("--tolerance", type=float, default=1e-6,
                        help="Maximum acceptable reconstruction error [K].")
    args = parser.parse_args()
    node_residual(args)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, IndexError, ZeroDivisionError) as exc:
        raise SystemExit(f"G6 RESIDUAL AUDIT ERROR: {exc}")
