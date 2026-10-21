# Flat Plate - Native CG Baseline

Full 284,400-tetrahedron flat plate with Spalart-Allmaras standard model.

- Mesh: `flat_plate_local.plt` (hard link to the canonical corrected mesh)
- Pressure solver: native conjugate gradient (`Pressure Field Solver = 1`)
- Timestep: `1.0e-5`
- Purpose: baseline for comparison with PETSc/AMG

The native CG solve does not converge on this mesh.  It reaches the 5,000
iteration cap with a true relative residual around `1.55e-4`.  Use
`../flat_plate_petsc` for production runs.

```powershell
Set-Location cases\flat_plate_local
..\..\build\local-serial-openmp\Release\cbs3dpp_si.exe flat_plate_local
```
