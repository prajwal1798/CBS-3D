# Heated panel with a U-channel: generic CHT unit cell

`HeatedPanel_UChannel_CHT.geo` is a parametric conjugate-heat-transfer cell:
a solid panel heated on one face and cooled by a square channel with two
parallel passes joined by a 180-degree bend. It is the common building block
of actively cooled heated walls (fusion first-wall and blanket panels, cold
plates, turbine-blade internal cooling) without copying any one design.

## Why this geometry

| Feature | Reason |
|---|---|
| Thin heated wall (3 mm) between the heated face and the coolant | Sets the peak solid temperature, as in real cooled walls; short conduction paths keep the solid's thermal time constant small (fast steady convergence) |
| Square channel, two passes, 180-degree bend (mean radius 1.25 a) | Dean-type secondary flow and bend heat-transfer enhancement; same family as published square-duct bend benchmarks |
| Unheated entrance (heating starts at x = 5a) | Separates hydrodynamic and thermal development |
| Ports on one face, unit-cell width 2p | One closed U-loop per cell; side faces are planes of symmetry of a periodic panel (adiabatic) |
| Conformal fluid/solid volumes (shared interface nodes) | Exact flux continuity in CBS3D's monolithic CHT coupling |

## Default dimensions (metres, all `DefineConstant`, change with `-setnumber`)

| Parameter | Value | Meaning |
|---|---|---|
| `a` | 0.010 | channel side (= hydraulic diameter) |
| `p` | 0.025 | leg centre spacing; bend mean radius p/2 |
| `Lleg` | 0.250 | straight leg length (25 a) |
| `tFront` / `tBack` | 0.003 / 0.007 | heated wall / back wall thickness |
| `xHeat` | 0.050 | heated face starts here (0 = heat the whole face) |
| `mEnd` | 0.005 | solid beyond the bend |

Derived: panel 272.5 x 50 x 20 mm, fluid volume 5.3927e-5 m3.

## Boundary and volume IDs

`10` fluid, `20` solid, `511` inlet (y = -p/2 leg), `520` outlet (y = +p/2 leg),
`532` heated face, `530` all other external solid faces (adiabatic).

## Mesh

Linear tetrahedra only. Fluid and solid are graded separately (Gmsh
`Restrict` fields): 0.15 mm at the interface, 0.45 mm fluid core, 0.5 mm in
the heated wall, 1.2 mm solid bulk (all times `hScale`, or `-clscale`).

| `hScale` | Tetrahedra | Use |
|---|---|---|
| 2 | 2.06 M (measured; worst minSICN 0.29, worst gamma 0.17) | quick look, laminar trials |
| 1.5 | about 4.9 M (estimate) | laminar production |
| 1 | about 16 M (estimate, 8x hScale 2) | finest; 128-640 MPI ranks |

Checked at `hScale` 2: conversion to CBS passes `mesh_summary.py`, and the
fluid and solid volumes match the exact values to 4e-7.

```bash
# laptop: quick look
gmsh cases/heated_panel_uchannel/HeatedPanel_UChannel_CHT.geo -3 -setnumber hScale 2 -o panel.msh
# Sunbird2: mesh -> CBS -> partitions on a compute node (case .par/.matprop still to be set)
sbatch $T/build_case.slurm $G PAR MATPROP $CBS3D_SCRATCH/cases/panel --clscale 1.5 --parts 256
```
