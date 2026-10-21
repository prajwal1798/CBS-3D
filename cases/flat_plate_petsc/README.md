# Flat Plate - PETSc / AMG

Full 284,400-tetrahedron flat plate with Spalart-Allmaras standard model and
the PETSc pressure solver.

- Pressure solver: PETSc KSP + AMG (`Pressure Field Solver = 3`)
- Preconditioner: GAMG (the MSYS2 PETSc has no hypre)
- Timestep: `1.0e-5`
- Iterations: 3 (stability smoke test, not a converged solution)
- Measured: 275/396/287 CG iterations, 5.9 s wall clock on 20 OpenMP threads

```powershell
$env:PATH = "C:\msys64\ucrt64\bin;$env:PATH"
$env:OMP_NUM_THREADS = "20"
Set-Location cases\flat_plate_petsc
..\..\build\local-msys2-petsc\cbs3dpp_si.exe flat_plate_petsc
```

The run writes `residual.csv`, `run_petsc.log` and the solver output files.
Mesh bytes are not tracked; see `docs/MESH_FORMAT.md` to regenerate them.

For a longer stability run, set `ntime` in the timestep line to `200` (or more).
A 200-iteration run on the tower took 2 min 31 s and ended at `RelU=2.06e-4`,
`RelV=2.42e-3`, `RelW=5.56e-3` and `RelP=2.41e-3`.  It remains far from steady
convergence; see `docs\SA_DEVELOPMENT.md` for the interpretation.
