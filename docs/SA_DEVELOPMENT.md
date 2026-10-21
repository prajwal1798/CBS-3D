# Spalart-Allmaras Development Status

## Model selection

The parameter file controls the turbulence model with:

```
Spalart-Allmaras controls: turbulence_on turbulence_model turbulent_thermal_diffusivity_on
1  0  0
```

| Value | Model |
| --- | --- |
| `turbulence_on = 0` | Laminar / no SA transport |
| `turbulence_model = 0` | Liu/Nithiarasu SA-noft2 |
| `turbulence_model = 1` | NASA SA-neg |

The SA-neg selection was added by local commit `1e3663d` and the production
branch by `e4038b7`.  These commits were not on `origin` and are preserved in
the master.

## Production sources

| File | Role |
| --- | --- |
| `src/assembly/SpalartAllmarasAssemblyLiu.cpp` | P1/TET4 SA assembly, source terms, eddy viscosity |
| `src/turbulence/SpalartAllmaras.cpp` | SA constants and pointwise kernels |
| `include/cbs/turbulence/SpalartAllmaras.hpp` | Public constants and helpers |
| `src/turbulence/WallDistance.cpp` | Wall distance and BVH |
| `src/turbulence/TurbulencePreprocess.cpp` | SA preprocessing and inlet/wall state |
| `src/boundary/TurbulenceBoundary.cpp` | SA boundary conditions |
| `src/timestep/TimeStep.cpp` | SA contribution to the local timestep |
| `src/turbulence/WallTreatment.cpp` | Spalding wall law |
| `src/turbulence/WallModelCoupling.cpp` | Momentum coupling of the wall model |
| `src/turbulence/CHTWallModelCoupling.cpp` | Conformal CHT wall treatment |

The historical `SpalartAllmarasAssembly.cpp` and
`SpalartAllmarasAssemblyVerified.cpp` are **not** compiled.  They are retained
only for A/B audit.  The production model is the Liu assembly.

## Liu SA-noft2 formulation

The assembly follows Chun-Bin Liu (2005), Chapters 3-4 and Appendix B:

- conservative convection `div(u nu_tilde)`;
- exact P1/TET4 CBS characteristic term `K_u_nu`;
- consistent P1 mass-weighted production/destruction;
- fully turbulent SA-noft2 form;
- destruction in the explicit Step-4 residual.

## NASA SA-neg branch

For negative `nu_tilde`, the negative branch uses:

```
chi = nu_tilde / nu
fv1 = 0
fn  = (cn1 + chi^3) / (cn1 - chi^3)
```

The production term is replaced by the negative-branch source, and the
destruction is modified so that negative states decay rather than grow.  The
eddy viscosity is zero for negative `nu_tilde`, so the momentum equation sees
only molecular viscosity.

The `sa_nu_tilde_floor` control is applied only to the standard model.  The
SA-neg model must be allowed to carry a negative state; applying the positive
floor would destroy the branch.

## Verification

| Test | Coverage |
| --- | --- |
| `sa_neg_negative_branch` | NASA SA-neg negative-state residual, preservation, zero eddy viscosity |
| `liu_nithiarasu_element` | Liu momentum and SA element residuals, including production/destruction |
| `wall_distance_exactness` | Wall-distance and BVH correctness |
| `wall_treatment_nasa_reference` | Spalding inversion against NASA reference values |
| `wall_treatment_reynolds_sweep` | Wall-law behaviour over a Reynolds sweep |
| `turbulence_wall_topology` | Wall-face topology and model-node inventory |

The retired `SAFiniteElementQuadratureVerification` test was removed because it
checked the old centroid/four-point source quadrature, not the production
Liu source.  The Liu element test covers the actual production term.

## Local flat-plate run

`cases\flat_plate_petsc` uses the 284,400-tetrahedron flat plate with SA on
(model 0), PETSc GAMG pressure, and three CBS iterations.  The run completes
without divergence:

```
Pressure solver : PETSc KSP + AMG
free stream      : U=1, T=300, Re=1.0e6
dt              : 1.0e-5
CG iterations   : 275, 396, 287
CG rel L2       : 5.54e-4, 1.59e-5, 6.22e-5
wall clock      : 5.9 s
```

The native CG path on the same mesh stalls at 5,000 iterations with a true
relative residual of about `1.87e-4` (canonical `U=1` free stream), which is why
the PETSc/AMG path is the default for this mesh.

### 200-iteration stability run

The same case was also run to `ntime=200` in a scratch directory with file
output disabled.  It completed in 2 min 31 s on 20 OpenMP threads.  The final
residuals were:

| Quantity | Value |
| --- | ---: |
| `RelU` | 2.06e-4 |
| `RelV` | 2.42e-3 |
| `RelW` | 5.56e-3 |
| `RelP` | 2.41e-3 |
| CG iterations (last solve) | 67 |
| CG absolute L2 (last solve) | 2.81e-13 |
| stop reason | maximum iteration count |

The run shows that the SA transport remains stable over a longer pseudo-time
window, but it is **not** converged: the velocity residuals are still falling
and the canonical case uses 30,000 iterations.  The last pressure solve reports
a large relative CG residual (`3.62e-2`) because it is normalized by the current
initial residual, which is already only `7.76e-12`; the final absolute residual
is `2.81e-13`.

### SA-neg full-case regression

The same 200-iteration case was rerun with `turbulence_model = 1` (NASA SA-neg).
The final residuals were `RelU=2.056e-4`, `RelV=2.404e-3`, `RelW=5.508e-3` and
`RelP=2.418e-3`, within 1% of the standard SA result at every sampled
iteration.  This is the expected result when `nu_tilde` stays positive: the
SA-neg branch is inactive and the run checks that model selection and the
production assembly are wired through the full solver.  The negative branch
itself is exercised by `sa_neg_negative_branch`.

## Current limitations

- The 284k flat plate is a coarse near-wall mesh; the first-cell `y+` is not
  suitable for a quantitative drag comparison.
- The 8-part case is partitioned and validated, but the local MSYS2 PETSc is
  MPIUNI, so eight-rank execution still requires Sunbird or a local MPI PETSc
  build.
- The packaged PETSc has no hypre, so the local runs use GAMG.  BoomerAMG
  comparisons require rebuilding PETSc with `--with-hypre=1`.
- The three-iteration case is a stability smoke test; the 200-iteration run is
  stable but remains far from the steady convergence tolerance.
- A 10x larger pseudo-time step (`dt=1.0e-4`) exceeded the SA `nu_tilde/nu`
  ceiling in its first iteration, so the baseline `dt=1.0e-5` is retained.

## Next steps

1. Extend the flat plate to true steady convergence (thousands of iterations,
   or a larger pseudo-time step/steady acceleration) and retain the full
   residual history.
2. Compare standard SA and SA-neg on the same case, including a forced
   negative-`nu_tilde` region or a rough-wall/inlet perturbation.
3. Perform a mesh-convergence study with the Gmsh `coarse`, corrected and
   refined flat-plate meshes.
4. Validate against the NASA TMR flat-plate skin-friction and velocity
   profiles once the near-wall resolution is adequate.
5. Rebuild PETSc with hypre and compare BoomerAMG against GAMG iteration
   counts and wall time.
6. Build a real MPI PETSc on the tower or use Sunbird for the 8-part and
   40-part cases.
