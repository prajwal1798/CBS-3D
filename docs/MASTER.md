# CBS3D Master

`C:\PhD\CFD\CBS3D_MASTER` is the single consolidated development tree for the
CBS-3D solver.  It is based on the authoritative `integration/sa-mpi` branch at
commit `1e3663d`, which is two local commits ahead of `origin`:

| Commit | Purpose |
| --- | --- |
| `e4038b7` | NASA SA-neg production branch |
| `1e3663d` | Enable SA-neg model selection in the parameter parser |

The master branch is `integration/master`.  The remote is
`https://github.com/prajwal1798/CBS-3D.git`.

## What is integrated

`integration/sa-mpi` already contains the work from:

- `feature/spalart-allmaras`
- `feature/domain-decomposition-sunbird`
- `feature/massflow-discrete-normalisation`
- `feature/temperature-afc`
- `integration/production-consolidation`

The older `feature/step4-validation`, `feature/sa-wall-treatment-foundation`
and `feature/sa-wall-topology-audit` branches are not merged because their
functionality is superseded by the integrated temperature-AFC and wall-model
implementations.  Their history remains available in the original repository.

The master adds:

- a working MSVC OpenMP build (`/openmp:llvm`, required by the OpenMP 3.1
  atomic-write pragmas in the Liu/SA kernels);
- a local MSYS2 UCRT64 GCC + PETSc build;
- a single CMake test suite covering the turbulence, wall-treatment, CHT and
  Spalart-Allmaras verification programs;
- a Gmsh-free mesh validator;
- local flat-plate cases produced from `C:\PhD\MESH\examples\3D`.

## Build configurations

### Native Windows serial / OpenMP

```powershell
cmake -S . -B build\local-serial-openmp -G "Visual Studio 18 2026" -A x64 `
      -DCBS3D_ENABLE_OPENMP=ON -DCBS3D_ENABLE_MPI=OFF -DCBS3D_ENABLE_PETSC=OFF
cmake --build build\local-serial-openmp --config Release --parallel 20
```

Executable: `build\local-serial-openmp\Release\cbs3dpp_si.exe`.

### Local PETSc / AMG (MSYS2 UCRT64)

The tower has no system PETSc.  MSYS2 UCRT64 provides a native GCC + PETSc
package, which is the no-administrator route used here.

```powershell
C:\msys64\usr\bin\bash.exe -lc "pacman -S --needed ucrt64/mingw-w64-ucrt-x86_64-gcc ucrt64/mingw-w64-ucrt-x86_64-cmake ucrt64/mingw-w64-ucrt-x86_64-ninja ucrt64/mingw-w64-ucrt-x86_64-petsc"

cmake --preset local-msys2-petsc
cmake --build --preset build-local-msys2-petsc
```

The preset sets `PETSC_DIR=C:/msys64/ucrt64`, `PETSC_ARCH=dto` and prepends
`C:/msys64/ucrt64/bin` to `PATH`.  The `dto` package is the double-precision,
OpenMP, MPIUNI build.  It provides PETSc GAMG; it is **not** an MPI build and
does not include hypre.

### MPI / Sunbird

The production distributed configuration remains:

```bash
cmake -S . -B build/mpi-release -DCBS3D_ENABLE_MPI=ON \
      -DCBS3D_ENABLE_PETSC=ON -DCBS3D_ENABLE_OPENMP=OFF
```

`PETSC_DIR` must point to an MPI-enabled PETSc with hypre if BoomerAMG is
desired.  This path is unchanged and is validated on Sunbird, not on the
20-core Windows tower.

## Tests

```powershell
cmake --preset tests-msvc-release
cmake --build --preset build-tests-msvc-release
ctest --preset test-tests-msvc-release
```

The suite builds:

| CTest name | What it checks |
| --- | --- |
| `step4_solid_volumetric_source` | Step-4 solid volumetric source |
| `wall_distance_exactness` | Exact wall-distance/BVH behaviour |
| `wall_treatment` | Spalding wall-law inversion |
| `wall_treatment_nasa_reference` | NASA reference wall-law values |
| `wall_treatment_reynolds_sweep` | Reynolds-number sweep |
| `thermal_wall_treatment` | Kader thermal wall treatment |
| `turbulence_wall_topology` | Wall-face topology and distance |
| `wall_model_coupling` | Momentum coupling of the wall model |
| `cht_wall_model_coupling` | Conformal CHT wall treatment |
| `sa_neg_negative_branch` | NASA SA-neg branch residual |
| `liu_nithiarasu_element` | Liu/Nithiarasu momentum and SA element residuals |

The old `SAFiniteElementQuadratureVerification` test was removed.  It checked
the retired centroid/four-point quadrature assembly, not the production
Liu assembly that now uses the consistent P1 mass-weighted source.  The
`liu_nithiarasu_element` test covers the production source term directly.

The two-rank MPI regressions are opt-in with `-DCBS3D_BUILD_MPI_TESTS=ON` and
a real MPI implementation.  The MSYS2 PETSc package cannot build them because
it is an MPIUNI (single-process) build.

## Mesh and cases

The native mesh workflow is documented in [MESH_FORMAT.md](MESH_FORMAT.md).
The local tower workflow is documented in [LOCAL_TOWER.md](LOCAL_TOWER.md).

Quick validation of any native case:

```powershell
python tools\mesh_summary.py cases\flat_plate_local\flat_plate_local --json cases\flat_plate_local\mesh_summary.json
```

The 284,400-tetrahedron flat plate validates as:

- 64,320 nodes
- 66,548 boundary faces
- total volume `5.000000000000e-02`
- BC counts `506=64640`, `511=474`, `520=474`, `530=960`
- result `PASS`

## Spalart-Allmaras status

See [SA_DEVELOPMENT.md](SA_DEVELOPMENT.md).  The production model is the
Liu/Nithiarasu SA-noft2 formulation; `turbulence_model=1` selects NASA SA-neg.
Both are covered by element-level verification.  The remaining work is
case-level validation on the flat plate and a mesh-convergence study, not
assembly implementation.

## Cleanup decisions

Removed from the master because they were dead, duplicated or generated:

- `src/main.cpp` — the application entry point is `apps/serial/main.cpp` or
  `apps/parallel/main.cpp`;
- the standalone `Makefile` — it referenced the dead entry point and the
  retired assemblies, and CMake is the authoritative build;
- seven one-shot `tools/apply_*.py` patch scripts;
- `tools/SolverProfiler.cpp` — duplicate of `src/utils/SolverProfiler.cpp`;
- `tools/CHT_UTube.msh` and `tools/.cht_patch_trigger` — generated artefacts;
- `tests/liu/CMakeLists.txt` — the Liu test is now a normal target in the
  master CMake suite;
- `tests/SAFiniteElementQuadratureVerification.cpp` — superseded by the Liu
  element test.

The historical experimental assemblies (`MomentumAssembly.cpp`,
`SpalartAllmarasAssembly.cpp`, `SpalartAllmarasAssemblyVerified.cpp`) are
kept in the tree for A/B audit but are not compiled.  The production sources
are the `*Liu.cpp` assemblies.
