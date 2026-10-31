# Manifold-fed cooled panel: parametric CHT benchmark geometry

`ManifoldPanel_CHT.geo`: a layered solid panel (back wall, coolant layer,
front wall, optional armour) heated on its front face and cooled by N
parallel square channels fed by an inlet header and drained by an outlet
header, U-type (ports on the same side) or Z-type (opposite sides).
All lengths are ratios of the channel side `a`, so the benchmark is defined
by dimensionless groups; `a` sets the scale.

What it adds over a single-channel cell: flow distribution between parallel
channels, lateral heat spreading between channels, layered solids, and a size
that needs a full HPC allocation.

| Parameter | Default | Meaning |
|---|---|---|
| `a` | 0.010 m | channel side (hydraulic diameter) |
| `nCh` | 8 | channels |
| `pR`, `LR` | 2, 30 | pitch / a, channel length / a |
| `wHR`, `hHR` | 2, 1 | header width / a, header height / a |
| `mR` | 1 | solid margin / a |
| `tFrontR`, `tBackR`, `tArmR` | 0.3, 0.7, 0.2 | front, back, armour thickness / a (armour 0 = none) |
| `manifold` | 0 | 0 U-type, 1 Z-type |

IDs: volumes `10` fluid, `20` structure, `21` armour (`solid_armour`; the
converter maps it to solid material 1 until multi-solid materials exist);
surfaces `511` inlet, `520` outlet, `532` heated face, `530` other external
faces (adiabatic). All interfaces are conformal.

Mesh: linear tetrahedra; fluid and solid graded separately (Restrict fields);
sizes AND grading distances scale with `hScale`, so the mesh is self-similar
and element quality does not depend on resolution.

| `hScale` | Tetrahedra | Quality |
|---|---|---|
| 3 | 10.2 M (measured: fluid 3.87 M, structure 5.99 M, armour 0.37 M) | worst gamma 0.21, worst minSICN 0.34, none below 0.2 |
| 2 | about 35 M (estimate, x (3/2)^3) | same grading |
| 1 | about 275 M (estimate) | same grading |

Default panel: 360 x 180 x 22 mm, fluid volume 3.06e-4 m3, total 1.4256e-3 m3
(conversion checked: mesh_summary PASS, total volume exact).
Large meshes: mesh, convert and partition on a Sunbird2 `highmem` node.
