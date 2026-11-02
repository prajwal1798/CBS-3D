#!/usr/bin/env python3
"""G7: one-iteration GLOBAL CBS3D conjugate thermal-energy and mass audit.

Read-only, uses only rank-local step-50000/50001 VTUs, the same .par and
.matprop as Job 29763, and the recorded dt. The full mesh includes fluid
and solid tetrahedra and their shared conformal interface.

Exact P1 weak-form identity on the whole thermal domain:

  sum_i (C_i/dt)(T1_i-T0_i)
    = Q_532 + Q_vol - H_adv + Q_Dirichlet

where Q_Dirichlet is the implied thermal exchange from imposing fixed
temperature nodes after the explicit residual update. H_adv is the
integrated, reference-relative, outward convective enthalpy flux through
the complete boundary of the fluid region, INCLUDING any fluid/solid
interface (which should have negligible normal velocity).

The equality of these independently computed terms is an algorithmic
conservation check for one pseudo-time iteration, NOT proof of physically
accurate long-term heat transfer or global steady balance.
"""

import argparse
import math
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path

from audit_uchannel_first_step_residual import (
    read_global_fields, find_timestep, param_after_comment, fluid_source,
)
from audit_cbs_explicit_timestep_vtu import (
    read_matprop, gradients_and_volume, dot, cross, vec_sub, norm,
)


def insist(cond, msg):
    if not cond:
        raise ValueError(msg)


def area_vector(p0, p1, p2, opposite):
    # Area-weighted *outward* normal from the adjacent tetrahedron.
    v = tuple(0.5 * x for x in cross(vec_sub(p1, p0), vec_sub(p2, p0)))
    midpoint = tuple((p0[d]+p1[d]+p2[d])/3 for d in range(3))
    inside = vec_sub(opposite, midpoint)
    if dot(v, inside) > 0:
        v = tuple(-x for x in v)
    return v


def nuclear_heating_profile_from_par(path):
    """Reproduce MeshIO optional nuclear_heating block, without reading
    or modifying any solver state. Return (origin, normal, profiles).
    An absent or disabled block contributes no nuclear source.
    """
    lines = [
        line.strip() for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith(("#", "!"))
    ]
    for i, label in enumerate(lines):
        if "nuclear_heating" not in label.lower():
            continue
        insist(i+1 < len(lines), "missing nuclear_heating data line")
        controls = lines[i+1].split()
        insist(len(controls) >= 8, "invalid nuclear_heating header")
        enabled = int(controls[0])
        origin = tuple(float(x) for x in controls[1:4])
        normal = tuple(float(x) for x in controls[4:7])
        count = int(controls[7])
        insist(enabled in (0, 1) and count >= 0,
               "invalid nuclear_heating enabled/count")
        if enabled == 0:
            return origin, normal, {}
        insist(all(math.isfinite(x) for x in origin+normal),
               "nonfinite nuclear plane or normal")
        insist(norm(normal) > 0.0, "zero nuclear-heating normal")
        profiles = {}
        if count:
            insist(i+2+count < len(lines), "missing nuclear-heating profile lines")
            for record in lines[i+3:i+3+count]:
                tokens = record.split()
                insist(len(tokens) >= 3, "incomplete nuclear-heating profile")
                material = int(tokens[0])
                amplitude, length = float(tokens[1]), float(tokens[2])
                insist(material not in profiles and amplitude >= 0 and
                       length > 0 and math.isfinite(amplitude) and
                       math.isfinite(length), "invalid nuclear-heating profile")
                profiles[material] = (amplitude, length)
        return origin, normal, profiles
    return (0.0, 0.0, 0.0), (0.0, 1.0, 0.0), {}


def nuclear_element_source(mid, centroid, origin, normal, profiles):
    if mid not in profiles:
        return None
    amplitude, length = profiles[mid]
    depth = dot(vec_sub(centroid, origin), normal) / norm(normal)
    return amplitude * math.exp(-max(0.0, depth) / length)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--par", type=Path, required=True)
    ap.add_argument("--matprop", type=Path, required=True)
    ap.add_argument("--before", type=int, default=50000)
    ap.add_argument("--after", type=int, default=50001)
    ap.add_argument("--reference-temperature", type=float, default=298.15)
    ap.add_argument("--tolerance-w", type=float, default=1e-5)
    a = ap.parse_args()

    old, _ = read_global_fields(a.output_dir, a.before, with_topology=False)
    current, pieces = read_global_fields(a.output_dir, a.after, with_topology=True)
    insist(old.keys() == current.keys(), "node IDs disagree between steps")
    dt, _ = find_timestep(a.output_dir, a.after)
    materials = read_matprop(a.matprop)
    nuclear_origin, nuclear_normal, nuclear_profiles = (
        nuclear_heating_profile_from_par(a.par)
    )
    inlet = param_after_comment(a.par, "reference and inlet/outlet data", 4)
    qwall = param_after_comment(a.par, "alpha_sf k_ratio source_solid", 3)
    qsolid = param_after_comment(a.par, "alpha_sf k_ratio source_solid", 2)
    insist(abs(inlet-a.reference_temperature)<1e-8, "inlet/reference T mismatch")
    insist(abs(qsolid)<1e-12, "solid .par source unsupported in this case")
    insist(dt>0 and math.isfinite(dt), "invalid timestep")

    # Independent nodal thermal residual and capacity reconstruction.
    capacity = defaultdict(float)            # J/K
    residual = defaultdict(float)            # W
    categories = Counter()
    cell_sets = set()
    node_materials = defaultdict(set)
    node_external_bcs = defaultdict(set)
    face_adjacency = defaultdict(list)        # key -> (material, opposite_gid)
    boundary_bcs = {}                        # exact boundary faces
    total_volume = 0.0
    qvol_total = 0.0
    qvol_by_material = defaultdict(float)
    divergent_volume = 0.0

    for file, piece in pieces:
        cursor = 0
        for icell, end in enumerate(piece["offsets"]):
            loc = piece["connectivity"][cursor:end]
            cursor = end
            vtk_type = piece["types"][icell]
            mid = piece["material_id"][icell]
            bc = piece["bc_id"][icell]
            ids = tuple(piece["gids"][k] for k in loc)

            if vtk_type == 5:
                insist(len(ids)==3, "non-triangle boundary face")
                key = tuple(sorted(ids))
                insist(key not in boundary_bcs,
                       f"duplicated boundary facet {key}")
                boundary_bcs[key] = bc
                for gid in ids:
                    node_external_bcs[gid].add(bc)
                continue
            if vtk_type != 10:
                continue
            insist(len(ids)==4, "non-TET4 volume cell")
            key = tuple(sorted(ids))
            insist(key not in cell_sets, f"duplicated tetrahedron {key}")
            cell_sets.add(key)
            insist(mid in materials, f"missing material {mid}")
            mat = materials[mid]
            for gid in ids:
                node_materials[gid].add(mid)

            p = [current[gid][0] for gid in ids]
            grad, V = gradients_and_volume(*p)
            categories[mid] += 1
            total_volume += V
            rhoCp = mat["rho_cp"]
            kappa = mat["k"]
            qvol = fluid_source(a.matprop, mid)
            centroid = tuple(sum(point[d] for point in p)/4 for d in range(3))
            qnuclear = nuclear_element_source(
                mid, centroid, nuclear_origin, nuclear_normal, nuclear_profiles
            )
            if qnuclear is not None:
                insist(qvol == 0.0 and (mid == 0 or qsolid == 0.0),
                       "nuclear-heating profile conflicts with constant source")
                vol_source = qnuclear
            else:
                vol_source = qvol if mid==0 else (qsolid if qsolid else qvol)
            qvol_total += vol_source*V
            qvol_by_material[mid] += vol_source*V
            T = [old[gid][1] for gid in ids]
            velocities = [current[gid][2] for gid in ids]
            gradT = tuple(sum(T[j]*grad[j][d] for j in range(4))
                          for d in range(3))
            div_u = sum(dot(velocities[j], grad[j]) for j in range(4)) if mid==0 else 0.0
            u_sum = tuple(sum(velocities[j][d] for j in range(4))
                          for d in range(3)) if mid==0 else (0.0,0.0,0.0)
            ubar = tuple(v/4 for v in u_sum)
            theta_sum = sum(Tj-inlet for Tj in T)
            if mid==0:
                divergent_volume += div_u*V

            for j,gid in enumerate(ids):
                capacity[gid] += rhoCp*V/4
                value = -kappa*V*dot(grad[j],gradT) + vol_source*V/4
                if mid==0:
                    uw = tuple(u_sum[d]+velocities[j][d] for d in range(3))
                    value -= rhoCp*V/20*(
                        dot(uw,gradT)+div_u*(theta_sum+(T[j]-inlet))
                    )
                    value -= 0.5*dt*rhoCp*V*dot(ubar,grad[j])*dot(ubar,gradT)
                residual[gid] += value

            for omit in range(4):
                face_ids = tuple(ids[j] for j in range(4) if j!=omit)
                face_adjacency[tuple(sorted(face_ids))].append((mid,ids[omit]))

    insist(set(capacity)==set(old), "unattached/extra thermal nodes")
    insist(all(v>0 for v in capacity.values()), "nonpositive heat capacity")

    # Include only prescribed external heat-flux cells (BC 532).
    qheat = 0.0
    fixed = set()
    fixed_bcs = {501,502,503,508,510,511}
    bc_faces = Counter()
    mass_bc = defaultdict(float)
    enthalpy_bc = defaultdict(float)
    fluid_external_enthalpy = 0.0
    fluid_interface_enthalpy = 0.0
    total_fluid_boundary_enthalpy = 0.0
    total_fluid_boundary_mass = 0.0

    for face, bc in boundary_bcs.items():
        bc_faces[bc] += 1
        insist(face in face_adjacency, f"orphan marked boundary face {face}")
        insist(len(face_adjacency[face])==1,
               f"marked face is not external to full mesh: {face}")
        mid, opposite = face_adjacency[face][0]
        p = [current[gid][0] for gid in face]
        av = area_vector(*p, current[opposite][0])
        if bc==532:
            qheat += qwall*norm(av)
            for gid in face:
                residual[gid] += qwall*norm(av)/3.0
        if bc in fixed_bcs:
            fixed.update(face)
        if mid==0:
            vn = [dot(current[gid][2],av) for gid in face]
            theta = [old[gid][1]-inlet for gid in face]
            mflux = materials[mid]["rho"]*sum(vn)/3
            hflux = materials[mid]["rho_cp"]/12*(
                sum(vn)*sum(theta)+sum(v*t for v,t in zip(vn,theta))
            )
            mass_bc[bc] += mflux
            enthalpy_bc[bc] += hflux

    # Reconstruct all outward fluid-boundary face integrals from actual tets,
    # not from externally labelled physical boundaries alone. The fluid-solid
    # interface has exactly one adjoining FLUID tet plus one SOLID tet.
    for face, adj in face_adjacency.items():
        insist(1<=len(adj)<=2, f"nonmanifold tetrahedron face {face}")
        fluid_adjacencies = [z for z in adj if z[0]==0]
        if len(fluid_adjacencies)!=1:
            continue
        _mid, opposite = fluid_adjacencies[0]
        p = [current[gid][0] for gid in face]
        av = area_vector(*p,current[opposite][0])
        vn = [dot(current[gid][2],av) for gid in face]
        theta = [old[gid][1]-inlet for gid in face]
        mflux = materials[0]["rho"]*sum(vn)/3
        hflux = materials[0]["rho_cp"]/12*(
            sum(vn)*sum(theta)+sum(v*t for v,t in zip(vn,theta))
        )
        total_fluid_boundary_enthalpy += hflux
        total_fluid_boundary_mass += mflux
        if len(adj)==2:
            insist(adj[0][0]!=adj[1][0], f"unrecognized fluid face {face}")
            fluid_interface_enthalpy += hflux
        else:
            fluid_external_enthalpy += hflux

    storage = sum(capacity[gid] * (current[gid][1]-old[gid][1])/dt
                  for gid in old)
    constrained_exchange = sum(
        capacity[gid]*(current[gid][1]-old[gid][1])/dt - residual[gid]
        for gid in fixed
    )
    interior_update_error = max(
        abs(capacity[gid]*(current[gid][1]-old[gid][1])/dt-residual[gid])
        for gid in old if gid not in fixed
    )

    # Independent exact boundary integral vs integrated volume conservative
    # transport, including any fluid/solid internal interface advective flux.
    adv_volume = 0.0
    discrete_r_sum = sum(residual.values())
    transport_from_volume = qheat + qvol_total - discrete_r_sum
    adv_face_mismatch = total_fluid_boundary_enthalpy - transport_from_volume
    predicted_storage = qheat + qvol_total - total_fluid_boundary_enthalpy + constrained_exchange
    closure = storage - predicted_storage
    scale = max(abs(qheat), abs(storage), abs(total_fluid_boundary_enthalpy),1.0)

    print("=== G7 FIRST-STEP GLOBAL THERMAL AND MASS BALANCE ===")
    print(f"Steps (before / after)       : {a.before} / {a.after}")
    print(f"Global dt [s]               : {dt:.14e}")
    print(f"Unique nodes / tetrahedra   : {len(old)} / {len(cell_sets)}")
    print(f"Fluid/solid tetra counts    : {categories[0]} / {sum(v for k,v in categories.items() if k!=0)}")
    print(f"External boundary IDs       : {dict(sorted(bc_faces.items()))}")
    print(f"Fixed temperature nodes     : {len(fixed)}")
    print(f"Prescribed Q_532 [W]         : {qheat:+.12e}")
    print(f"Volumetric source [W]       : {qvol_total:+.12e}")
    print(f"Nuclear profile materials   : {sorted(nuclear_profiles)}")
    for mid in sorted(qvol_by_material):
        print(f"  material {mid} volumetric Q : {qvol_by_material[mid]:+.12e} W")
    print(f"Fluid boundary enthalpy [W] : {total_fluid_boundary_enthalpy:+.12e}")
    print(f"  external only [W]         : {fluid_external_enthalpy:+.12e}")
    print(f"  solid interface [W]       : {fluid_interface_enthalpy:+.12e}")
    print(f"Fluid net mass flux [kg/s]  : {total_fluid_boundary_mass:+.12e}")
    print(f"Fluid integrated divu [m3/s]: {divergent_volume:+.12e}")
    for bc in sorted(mass_bc):
        print(f"  BC {bc}: mdot={mass_bc[bc]:+.12e} kg/s,"
              f" relative Q={enthalpy_bc[bc]:+.12e} W")
    print(f"Fixed-T constraint exchange : {constrained_exchange:+.12e} W")
    print(f"Thermal storage/dt [W]      : {storage:+.12e}")
    print(f"Predicted storage/dt [W]    : {predicted_storage:+.12e}")
    print(f"Global closure error [W]    : {closure:+.12e}")
    print(f"Relative closure error      : {abs(closure)/scale:.12e}")
    print(f"Adv volume/surface err [W]  : {adv_face_mismatch:+.12e}")
    print(f"Worst free-node residual err : {interior_update_error:.12e} W")
    print("=== NODE-LEVEL UPDATE MISMATCH DIAGNOSIS (no gate relaxation) ===")
    mismatches = []
    classification = defaultdict(lambda: [0,0.0,0.0,0.0])
    for gid in old:
        error_i = capacity[gid]*(current[gid][1]-old[gid][1])/dt - residual[gid]
        mats = node_materials[gid]
        if 0 in mats and len(mats)>1:
            label = "interface"
        elif 0 in mats:
            label = "fluid-only"
        else:
            label = "solid-only"
        if gid in fixed:
            label = "FIXED-" + label
        group = classification[label]
        group[0] += 1
        group[1] += error_i
        group[2] += abs(error_i)
        group[3] = max(group[3], abs(error_i))
        if gid not in fixed:
            mismatches.append((abs(error_i),gid,error_i,label))
    for label, (count,total,absolute,maximum) in sorted(classification.items()):
        print(f"  {label:18s} N={count:6d} sum_err={total:+.12e} W"
              f" abs_sum={absolute:.12e} W max_abs={maximum:.12e} W")
    print("  Worst 20 FREE nodes (error = observed C*deltaT/dt - reconstructed RHS):")
    for _abs,gid,err,label in sorted(mismatches, reverse=True)[:20]:
        print(f"    gid={gid:6d} error={err:+.12e} W"
              f" class={label} BCs={sorted(node_external_bcs.get(gid,set()))}"
              f" T0={old[gid][1]:.9f} T1={current[gid][1]:.9f}"
              f" R={residual[gid]:+.9e} W")
    print("Note: this is ONE pseudo-time update, not a steady-state"
          " Qout=Qin verification.")
    insist(abs(closure)<a.tolerance_w,
           "global thermal energy balance does not close")
    insist(abs(adv_face_mismatch)<a.tolerance_w,
           "volume/surface divergence-theorem check failed")
    insist(interior_update_error < a.tolerance_w,
           "unconstrained node update inconsistent with reconstructed RHS")
    print("G7 FIRST-STEP GLOBAL THERMAL CONSERVATION PASS")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, IndexError, ZeroDivisionError) as exc:
        raise SystemExit(f"G7 AUDIT ERROR: {exc}")
