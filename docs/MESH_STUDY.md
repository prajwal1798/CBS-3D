# Flat-Plate Mesh Study

## Geometry family

Use the `.geo` family in `C:\PhD\MESH\examples\3D\flat_plate`:

| Level | Source | Tetrahedra | Nodes | Boundary faces |
| --- | --- | ---: | ---: | ---: |
| corrected | `flat_plate_sa_3d_corrected.geo` | 284,400 | 64,320 | 66,548 |
| refined | `flat_plate_sa_3d_refined.geo` | 1,235,520 | 234,900 | 114,448 |
| quality | `flat_plate_sa_3d_quality.geo` | about 6.48M | not generated locally | - |

All three use the same domain: `x[-0.5, 2]`, `y[0, 1]`, `z[0, 0.02]`,
volume `5.000000000000e-02`, and BCs `506` symmetry, `511` inlet, `520`
outlet, `530` plate wall.  The refined mesh was generated with the Gmsh Python
module and converted with:

```powershell
python tools\gmsh_to_cbs.py `
    C:\Users\prajw\AppData\Local\Temp\opencode\flat_plate_refined.msh `
    --case-name flat_plate_refined `
    --out-dir C:\Users\prajw\AppData\Local\Temp\opencode\refined_alt_case2 `
    --par C:\PhD\MESH\examples\3D\flat_plate\flat_plate_sa_3d_corrected.par `
    --matprop C:\PhD\MESH\examples\3D\flat_plate\flat_plate_sa_3d_corrected.matprop
```

The converter handles the `Alternate` transfinite surface triangulation; see
[MESH_FORMAT.md](MESH_FORMAT.md).

> Do **not** use the legacy `coarse.msh` or `fine.msh` as refinement levels.
> The coarse mesh starts at `x=-1/3`, and the fine mesh spans `z=1`; neither
> matches the canonical thin domain.

## Current PETSc/AMG smoke results

Both runs use the standard SA model, PETSc GAMG, `dt=1.0e-5`, and 20 OpenMP
threads.

| Mesh | Cells | Iterations | CG iterations | Wall clock |
| --- | ---: | ---: | --- | ---: |
| corrected | 284,400 | 3 | 275 / 396 / 287 | 5.9 s |
| refined | 1,235,520 | 3 | 356 / 429 / 280 | 35.1 s |

The refined run ended at `RelU=7.31e-3`, `RelV=1.62e-1`, `RelW=3.64e-1`,
`RelP=1.49e-1`.  These are three-iteration startup values, not a converged
solution.  The last pressure solve had an initial true residual of
`2.36e-9`, a final true residual of `1.75e-12`, and a maximum absolute
residual of `7.35e-14`.

Profiler timings for the refined run were about 6.4-9.5 s per pressure solve
and 3.36 s per SA assembly.  The observed working set was about 1 GB.  The
corrected run is the committed `cases\flat_plate_petsc` smoke case; the refined
run is a scratch case.

### Timestep test

A corrected-mesh run with `dt=1.0e-4` (10x the baseline) failed in its first
iteration:

```text
SpalartAllmarasAssembly::updateEddyViscosity - nu_tilde/nu ceiling exceeded
at element 54037, ratio=2727929.341793
```

The baseline `dt=1.0e-5` is stable on the corrected mesh, but it is **not**
stable on the refined mesh.  The refined 200-iteration fixed-step run failed
after the first completed iteration at element 267301, ratio
`7.38887811475e11`.  That element is the first cell off the symmetry plane just
upstream of the leading edge (`x ~ -0.0116`, `y = 0` to `7.98e-6`).  Its minimum
altitude is `7.98e-6`.

The SA transport equation is advanced explicitly, so the solver's diffusive
limit is

```text
dt_diff = h^2 / (2 D_SA),    D_SA = (nu + nu_tilde) / sigma.
```

At the freestream SA value `nu_tilde = 3 nu = 3e-6`, with `nu = 1e-6` and
`sigma = 2/3`, this gives `dt_diff ~ 5.3e-6` for the refined cell, below the
fixed `1.0e-5`.  The corresponding corrected-mesh cell has `h = 2.12e-5` and
`dt_diff ~ 3.76e-5`, so `1.0e-5` is stable there.  The corrected `dt=1.0e-4`
failure is likewise above its local limit.  This is a timestep-stability
failure, not evidence of a model-form error.

Local pseudo-time stepping (`ilots=1`, `dtfixed=0`, with a positive `dtfix`
placeholder) removes the fixed-step restriction.  A refined diagnostic
completed 10 iterations in `1m 32s`, with `dt` settling at `2.56e-6`, no SA
ceiling, `CG=195`, and a final true pressure residual of about `5e-13`.

A 200-iteration refined LTS run completed in `20m 52s` with final `dt=2.04e-6`,
`RelU=1.78e-4`, `RelV=2.01e-3`, `RelW=3.76e-3`, `RelP=2.22e-3`, last `CG=77`,
final true pressure residual `1.14e-13`, and maximum absolute residual
`5.46e-15`.

The matching corrected-mesh LTS run completed 200 iterations in `2m 40s` with
final `dt=7.64e-6`, `RelU=1.92e-4`, `RelV=2.31e-3`, `RelW=5.32e-3`,
`RelP=2.36e-3`, last `CG=67`, final true pressure residual `2.43e-13`, and
maximum absolute residual `2.33e-14`.  The corrected step is limited mainly by
advection (about `7.6e-6`), while the refined step is limited by SA diffusion
and settles near `2.0e-6`.  At 200 iterations the two meshes have very similar
residual magnitudes.

## What is established

- The 1.235M refined mesh is valid and runs through the full 3D SA/PETSc path.
- Fixed `dt=1.0e-5` is stable on the corrected mesh but outside the explicit
  SA diffusive limit on the refined mesh.  Refined runs should use LTS or a
  substantially smaller fixed step.
- Refined and corrected LTS are both stable for 200 iterations and reach
  similar residual magnitudes.
- The refined LTS 200-iteration run takes about `7.8x` the corrected LTS wall
  clock while increasing the cell count by 4.35x.  SA assembly is about 62% of
  the refined iteration cost versus 38% for the corrected mesh.
- The refined mesh costs about six times the corrected wall clock for the same
  three fixed-step iterations, while increasing the cell count by 4.35x.
- The pressure absolute residuals remain small; the reported `CGrel` is
  normalized by an already tiny initial residual (see
  [LOCAL_TOWER.md](LOCAL_TOWER.md)).

This is **not** a mesh-convergence result.  A convergence study needs a common
converged state (or a common physical time) on each level, plus an integral
quantity such as plate drag or skin friction.  The 200-iteration corrected and
refined runs are still far from the steady tolerance.

## Next steps

1. Run corrected and refined LTS to the same physical time or steady tolerance
   and compare an integral quantity.
2. Use `ilots=1` local pseudo-time stepping for the refined mesh.  A 10x larger
   fixed step (`dt=1.0e-4`) is unstable and cannot be used as-is.
3. Run the quality mesh on Sunbird; it is likely too large for the 16 GB tower.
4. Keep the legacy `coarse.msh`/`fine.msh` out of the convergence claim.
