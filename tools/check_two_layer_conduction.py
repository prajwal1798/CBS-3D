#!/usr/bin/env python3
"""G4: verify conformal, two-material heat conduction against an exact solution.

This checks the existing serial/MPI VTU writer without changing production
numerics. Both materials share an interface at y=Lf. The fluid is stagnant.

T(0)=0, -k dT/dn at the top corresponds to prescribed *inward* q>0.
For 0<=y<=Lf:  T = q*y/k_fluid.
For Lf<y<=Lf+Ls: T = q*Lf/k_fluid + q*(y-Lf)/k_solid.
"""
import argparse
from collections import defaultdict, Counter
import math
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_cht_temperature_vtu import read_piece, iter_cells, triangle_area


def collect_files(target):
    target = Path(target)
    if target.is_dir():
        pvtus = sorted(target.rglob("*.pvtu"), key=lambda p: p.stat().st_mtime)
        if pvtus:
            target = pvtus[-1]
        else:
            serial = [p for p in target.rglob("*.vtu")
                      if "_rank_" not in p.name and "boundary" not in p.name]
            if not serial:
                raise ValueError("No serial VTU or MPI PVTU found in " + str(target))
            target = max(serial, key=lambda p: p.stat().st_mtime)
    if target.suffix == ".vtu":
        return [target]
    if target.suffix != ".pvtu":
        raise ValueError("Unsupported result path " + str(target))
    root = ET.parse(target).getroot()
    files = [target.parent / e.attrib["Source"]
             for e in root.findall(".//PUnstructuredGrid/Piece")]
    if not files:
        raise ValueError("No PVTU Piece entries in " + str(target))
    return files


def exact(y, lf, kf, ks, q):
    if y <= lf:
        return q * y / kf
    return q * lf / kf + q * (y - lf) / ks


def read_result(target, lf, ls, kf, ks, q, width):
    pieces = collect_files(target)
    nodes = {}
    incident = defaultdict(set)
    boundary_faces = {}
    cell_counts = Counter()
    gradients = {0: [], 1: []}
    maximum_speed = 0.0

    for path in pieces:
        piece = read_piece(path)
        root = ET.parse(path).getroot()
        pdata = root.find("./UnstructuredGrid/Piece/PointData")
        if pdata is None:
            raise ValueError("Missing PointData: " + str(path))
        vel_node = next((n for n in pdata.findall("DataArray")
                         if n.get("Name") == "velocity"), None)
        if vel_node is None:
            raise ValueError("Missing velocity PointData: " + str(path))
        flat_velocity = np.fromstring(vel_node.text or "", sep=" ")
        if flat_velocity.size != 3 * len(piece.points):
            raise ValueError("Incorrect velocity field length: " + str(path))
        velocity = flat_velocity.reshape(-1, 3)
        maximum_speed = max(maximum_speed, float(np.max(np.linalg.norm(velocity, axis=1))))

        keys = [tuple(round(float(x), 12) for x in point) for point in piece.points]
        for key, t in zip(keys, piece.temperature):
            t = float(t)
            if key in nodes and abs(nodes[key] - t) > 1.0e-9:
                raise ValueError("MPI duplicate-node temperature mismatch at " + str(key))
            nodes[key] = t

        for _cell, kind, ids, material, bc in iter_cells(piece):
            if kind == 10:
                if material not in (0, 1):
                    raise ValueError("Unexpected material ID " + str(material))
                cell_counts[material] += 1
                for n in ids:
                    incident[keys[n]].add(material)
                p = np.array([piece.points[n] for n in ids], dtype=float)
                t = np.array([piece.temperature[n] for n in ids], dtype=float)
                matrix = p[1:] - p[0]
                grad = np.linalg.solve(matrix, t[1:] - t[0])
                if min(abs(p[:, 1] - lf)) <= 1.0e-9:
                    gradients[material].append(float(grad[1]))
            elif kind == 5:
                facekey = tuple(sorted(keys[n] for n in ids))
                if facekey in boundary_faces:
                    # No duplicate physical faces across partitions are expected.
                    raise ValueError("Duplicate physical boundary triangle " + str(facekey))
                area = triangle_area(*(piece.points[n] for n in ids))
                boundary_faces[facekey] = (bc, area)

    if not nodes or not cell_counts[0] or not cell_counts[1]:
        raise ValueError("Missing nodes or one material")
    if not all(math.isfinite(x) for x in nodes.values()):
        raise ValueError("Nonfinite temperature in VTU")

    top = exact(lf + ls, lf, kf, ks, q)
    interface = exact(lf, lf, kf, ks, q)
    differences = [abs(t - exact(y, lf, kf, ks, q)) for (_, y, _), t in nodes.items()]
    err = max(differences)
    interface_nodes = [k for k in nodes
                       if incident[k] == {0, 1} and abs(k[1] - lf) < 1.0e-9]

    face_counts = Counter()
    face_areas = defaultdict(float)
    for bc, area in boundary_faces.values():
        face_counts[bc] += 1
        face_areas[bc] += area
    boundary_ids = set(face_counts)
    hot_area = face_areas[532]
    cold_area = face_areas[502]
    area_target = width * width

    print("VTU pieces                  :", len(pieces))
    print("Unique nodes                :", len(nodes))
    print("Element counts fluid/solid  :", cell_counts[0], cell_counts[1])
    print("Shared interface nodes      :", len(interface_nodes))
    print("External boundary IDs       :", sorted(boundary_ids))
    print("BC 502 area                 : {:.12e} m2".format(cold_area))
    print("BC 532 area                 : {:.12e} m2".format(hot_area))
    print("Integrated prescribed power : {:.12e} W".format(q * hot_area))
    print("Exact T(interface)          : {:.9f} K".format(interface))
    print("Exact T(top)                : {:.9f} K".format(top))
    print("Numerical Tmin/Tmax         : {:.9f} / {:.9f} K".format(min(nodes.values()), max(nodes.values())))
    print("Max nodal analytic error    : {:.12e} K ({:.5%} Ttop)".format(err, err / top))
    print("Maximum fluid speed         : {:.12e} m/s".format(maximum_speed))
    if interface_nodes:
        ierr = max(abs(nodes[k] - interface) for k in interface_nodes)
        print("Max interface temperature error: {:.12e} K".format(ierr))
    interface_flux_median = {}
    for material, conductivity in ((0, kf), (1, ks)):
        if gradients[material]:
            normalflux = [conductivity * grad for grad in gradients[material]]
            flux_median = float(np.median(normalflux))
            interface_flux_median[material] = flux_median
            print("Material {} interface-near k*dT/dy median = {:.9e} W/m2".format(
                material, flux_median))

    errors = []
    if boundary_ids != {502, 520, 530, 532}:
        errors.append("Unexpected/missing physical boundary IDs")
    if not (len(interface_nodes) >= 4):
        errors.append("No verified shared conformal interface nodes")
    if abs(hot_area - area_target) > 1.0e-3 * area_target:
        errors.append("Incorrect heated external area")
    if abs(cold_area - area_target) > 1.0e-3 * area_target:
        errors.append("Incorrect fixed-temperature external area")
    if maximum_speed > 1.0e-5:
        errors.append("Fluid does not remain effectively stagnant")
    if err > 0.005 * top:
        errors.append("Analytical temperature error exceeds 0.5% of Ttop")
    if min(nodes.values()) < -0.005 * top:
        errors.append("Unphysical temperature undershoot for no-negative-source case")
    for material in (0, 1):
        if material not in interface_flux_median:
            errors.append("Missing interface-adjacent tetrahedral gradients for material {}".format(material))
        elif abs(interface_flux_median[material] - q) > 0.05 * q:
            errors.append("Median interface heat-flux error exceeds 5% for material {}".format(material))

    for problem in errors:
        print("FAIL:", problem)
    if errors:
        return None, False
    print("G4 FIELD PASS")
    return nodes, True


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("result")
    p.add_argument("--compare", help="Optional serial result to compare with MPI")
    p.add_argument("--Lf", type=float, default=0.001)
    p.add_argument("--Ls", type=float, default=0.003)
    p.add_argument("--W", type=float, default=0.004)
    p.add_argument("--kf", type=float, default=0.606)
    p.add_argument("--ks", type=float, default=28.3)
    p.add_argument("--q", type=float, default=1000.0)
    a = p.parse_args()
    nodes, passed = read_result(a.result, a.Lf, a.Ls, a.kf, a.ks, a.q, a.W)
    if not passed:
        return 1
    if a.compare:
        baseline, other_ok = read_result(a.compare, a.Lf, a.Ls, a.kf, a.ks, a.q, a.W)
        if not other_ok:
            return 1
        if set(nodes) != set(baseline):
            print("FAIL: MPI/serial node-coordinate sets differ")
            return 1
        biggest = max(abs(nodes[k] - baseline[k]) for k in nodes)
        top = exact(a.Lf + a.Ls, a.Lf, a.kf, a.ks, a.q)
        print("MPI/serial maximum nodal difference = {:.12e} K".format(biggest))
        if biggest > 1.0e-6 * max(top, 1.0):
            print("FAIL: MPI/serial temperature parity")
            return 1
    print("G4 TWO-MATERIAL CONDUCTION PASS")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError, np.linalg.LinAlgError) as e:
        print("G4 ERROR:", e, file=sys.stderr)
        raise SystemExit(1)
