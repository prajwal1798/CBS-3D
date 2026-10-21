# CBS3D++_SI

Three-dimensional tetrahedral incompressible-flow solver using the
characteristic-based split (CBS) method, with Spalart-Allmaras turbulence,
wall modelling and conjugate heat transfer.

This tree is the consolidated master at `C:\PhD\CFD\CBS3D_MASTER`.
It is based on the local `integration/sa-mpi` branch and includes the two
SA-neg commits that were never pushed (`e4038b7`, `1e3663d`).

## Quick start

```powershell
cmake -S . -B build\local-serial-openmp -G "Visual Studio 18 2026" -A x64 `
      -DCBS3D_ENABLE_OPENMP=ON -DCBS3D_ENABLE_PETSC=OFF
cmake --build build\local-serial-openmp --config Release --parallel 20

Set-Location cases\flat_plate_local
..\..\build\local-serial-openmp\Release\cbs3dpp_si.exe flat_plate_local
```

The native CG solver is not adequate for the 284k-cell flat plate.  Use the
PETSc/AMG build for that case:

```powershell
cmake --preset local-msys2-petsc
cmake --build --preset build-local-msys2-petsc

$env:PATH = "C:\msys64\ucrt64\bin;$env:PATH"
$env:OMP_NUM_THREADS = "20"
Set-Location cases\flat_plate_petsc
..\..\build\local-msys2-petsc\cbs3dpp_si.exe flat_plate_petsc
```

## Documentation

- [docs/MASTER.md](docs/MASTER.md) — branch consolidation, build matrix,
  test suite, cleanup record.
- [docs/LOCAL_TOWER.md](docs/LOCAL_TOWER.md) — 20-core Windows workflow,
  MSYS2 PETSc, measured flat-plate timings.
- [docs/MESH_FORMAT.md](docs/MESH_FORMAT.md) — native `.plt/.bco/.material`
  format, boundary codes, Gmsh partition/export workflow.
- [docs/SA_DEVELOPMENT.md](docs/SA_DEVELOPMENT.md) — SA/SA-neg status and
  next validation steps.

## Tests

```powershell
cmake --preset tests-msvc-release
cmake --build --preset build-tests-msvc-release
ctest --preset test-tests-msvc-release
```

## Layout

```
apps/            serial and parallel entry points
include/cbs/     public headers
src/             solver, assemblies, turbulence, I/O, linear algebra
tests/           verification and regression programs
cases/           local flat-plate cases
tools/           mesh validation and diagnostics
jobs/            Sunbird Slurm jobs
docs/            master, mesh, tower and SA documentation
```
