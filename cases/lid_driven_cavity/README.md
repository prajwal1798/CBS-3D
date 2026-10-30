# Lid-driven cavity cases

| Case | Geometry | Use |
|---|---|---|
| `ldc3d` | unit cube, no-slip side walls, lid y = 1 moving in +x | 3D cavity (side walls weaken the mid-plane vortex; not comparable 1:1 with 2D data) |
| `ldc_slab2d` | 1 x 1 x 2h slab, symmetry (BC 506) on z = 0 and z = 2h | 2D-equivalent cavity: compare directly with Ghia et al. (1982) |

Both: Re = 100, energy off, fixed dt = 1e-3 (dt = 2e-3 is unstable on
these meshes), steady-state stop when the largest relative velocity change
per step is below 1e-7.

Verification of `ldc_slab2d` (30 Oct 2026, 2 MPI ranks, steady state reached
in ~13,300 steps; profiles at mid-depth, `tools/plot_ldc_centerlines.py`):

| Mesh | dt | u min | v max | v min | RMS vs Ghia (u / v) |
|---|---|---|---|---|---|
| h = 1/64 | 1e-3 | -0.1979 | 0.1673 | -0.2287 | 0.0089 / 0.0073 |
| h = 0.0125 | 1e-3 | -0.2000 | 0.1689 | -0.2316 | 0.0084 / 0.0058 |
| h = 0.0125 | 5e-4 | -0.2000 | 0.1689 | -0.2317 | 0.0085 / 0.0057 |
| Ghia et al. (1982), 129 x 129 | | -0.2109 | 0.1753 | -0.2453 | |

The time step does not change the steady answer; the remaining few-percent
gap shrinks with mesh refinement (spatial error of linear tetrahedra).

```bash
bash tools/cbs_case.sh cases/lid_driven_cavity/ldc_slab2d.geo \
    cases/lid_driven_cavity/ldc_slab2d.par cases/lid_driven_cavity/ldc_slab2d.matprop \
    $CBS3D_SCRATCH/cases/ldc_slab2d --parts 8
```
