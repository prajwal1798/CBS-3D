# CBS3D production consolidation: branch audit and integration plan

Date: 2026-08-06

Working branch: `integration/production-consolidation`

Base branch: `feature/massflow-discrete-normalisation`

## 1. Branch topology

The repository currently contains six long-lived branches:

- `main`
- `feature/spalart-allmaras`
- `feature/domain-decomposition-sunbird`
- `feature/temperature-afc`
- `feature/massflow-discrete-normalisation`
- `feature/step4-validation`

The principal development lineage is:

```text
feature/spalart-allmaras
    -> feature/domain-decomposition-sunbird
        -> feature/temperature-afc
            -> feature/massflow-discrete-normalisation
```

`feature/massflow-discrete-normalisation` is therefore the most complete
validated laminar production lineage. It contains the SA source modules,
domain-decomposition infrastructure, distributed PETSc pressure system,
restart support, temperature AFC wrapper and the corrected discrete BC511
normalisation.

`feature/step4-validation` diverged from the same development family and
contains valuable Step-4/AFC tests and heat-transfer post-processing work, but
it does not contain the final BC511 correction and is not suitable as the
production base without selective integration.

`main` contains public documentation, CI, contribution templates and the
initial examples documentation, but its solver source is older than the
production branches.

## 2. Current capability audit

| Capability | Current status | Production status |
|---|---|---|
| Serial laminar | Implemented | Validated baseline |
| Serial OpenMP | Implemented through OpenMP build option | Requires formal regression matrix |
| Serial PETSc pressure solve | Implemented through `solver_opt=3` and PETSc build option | Requires dedicated regression and packaging |
| MPI distributed laminar | Implemented with PETSc persistent pressure system | Validated on Sunbird |
| Hybrid MPI + OpenMP | Build combination exists | Not yet production-validated; current benchmark used one OpenMP thread per MPI rank |
| Serial Spalart-Allmaras | Source path exists | Development/validation incomplete |
| Distributed Spalart-Allmaras | Not integrated into the production distributed loop | Unsupported |
| CHT temperature AFC | Production wrapper exists | Validated for the current CHT benchmark |
| Discrete BC511 mass flow | Implemented for serial and distributed preprocessing | Validated on Sunbird, 40 ranks |
| Restart/checkpoint | Implemented for distributed production wrapper | Needs consolidation and formal tests |

## 3. Important technical findings

### 3.1 One codebase should be retained

Serial and parallel calculations should not be maintained as duplicated solver
source trees. Duplicate numerical kernels will diverge and reproduce the
current branch/folder problem.

The production repository should contain:

- one shared numerical core;
- a thin serial application target;
- a thin MPI application target;
- separate build presets and build directories;
- separate serial and parallel launch scripts.

This gives visibly separate serial and parallel applications without
maintaining two copies of momentum, pressure, energy, boundary and turbulence
kernels.

### 3.2 The current distributed production loop is laminar only

The distributed production loop executes momentum, pressure, velocity
correction and energy, but does not advance the Spalart-Allmaras equation.
Until distributed SA assembly, halo exchange, constraints, convergence and
output are implemented, an MPI run with `turbulence_on=1` must fail explicitly
rather than silently execute a laminar calculation.

### 3.3 The production wrapper structure must be simplified

The current restart/AFC implementation compiles the distributed production
loop by including `.cpp` translation units under macro substitutions. This
preserved a validated loop during development, but it is not an appropriate
long-term production architecture.

Restart iteration mapping, checkpointing and AFC selection should become
normal classes/policies called by one compiled production loop.

### 3.4 Step-4 branches must be reconciled selectively

The repository currently contains both `TemperatureAFC` and `ThermalAfc`
development approaches. The validated production implementation and the
standalone Step-4 test suite must be compared and reduced to one implementation
with one naming convention and one test matrix.

### 3.5 Output units must be explicit

CBS pressure is stored internally as kinematic pressure. Output and
post-processing must state this convention explicitly or write a second
physical-pressure field in pascals when dimensional mode is active. The recent
straight-pipe post-processing required `p_Pa = rho * p_kinematic`.

## 4. Target production layout

```text
CBS-3D/
├── apps/
│   ├── serial/
│   │   └── main.cpp
│   └── parallel/
│       └── main.cpp
├── include/cbs/
├── src/
│   ├── assembly/
│   ├── boundary/
│   ├── core/
│   ├── io/
│   ├── linalg/
│   ├── parallel/
│   ├── preprocess/
│   ├── solver/
│   ├── timestep/
│   ├── turbulence/
│   └── utils/
├── examples/
│   ├── Blanket_Laminar/
│   ├── Blanket_Turbulent/
│   ├── Pipe_Flow/
│   └── Flat_Plate_Turbulent/
├── tests/
│   ├── serial/
│   ├── parallel/
│   ├── pressure/
│   ├── thermal/
│   └── turbulence/
├── cmake/
├── jobs/
├── tools/
├── docs/
├── CMakeLists.txt
└── CMakePresets.json
```

The two application folders are intentionally separate. The numerical source
is intentionally shared.

## 5. Planned executable/build matrix

| Preset | Target | MPI | OpenMP | PETSc | Intended use |
|---|---|---:|---:|---:|---|
| `serial-release` | `cbs3d_serial` | OFF | OFF | OFF | deterministic serial reference |
| `serial-openmp-release` | `cbs3d_serial` | OFF | ON | OFF | workstation multicore |
| `serial-petsc-release` | `cbs3d_serial` | OFF | optional | ON | PETSc-preconditioned serial pressure solve |
| `mpi-release` | `cbs3d_mpi` | ON | OFF | ON | production distributed run |
| `hybrid-release` | `cbs3d_mpi` | ON | ON | ON | MPI + OpenMP after scaling validation |
| `sa-serial-dev` | `cbs3d_serial` | OFF | ON | selectable | SA development |
| `sa-mpi-dev` | `cbs3d_mpi` | ON | ON | ON | future distributed SA development |

## 6. Integration order

1. Preserve `feature/massflow-discrete-normalisation` as the immutable validated
   source of the corrected 1M benchmark.
2. Merge documentation, CI and repository governance files from `main` without
   replacing validated production source.
3. Add CMake presets and produce separate serial/parallel application targets
   over the shared core.
4. Remove `.cpp`-inclusion wrappers and replace them with normal restart/AFC
   interfaces.
5. Import the Step-4 test suite selectively and reconcile the AFC
   implementations.
6. Add explicit rejection for distributed SA until the distributed SA path is
   complete.
7. Build the examples hierarchy and add manifest files documenting every input,
   mesh, material, boundary condition, expected output and validation metric.
8. Establish CI for serial GCC/Clang/MSVC builds and a Sunbird regression script
   for MPI/PETSc.
9. Validate all build modes against a common reference-result manifest.
10. Only after validation, merge the consolidation branch to `main` and archive
    obsolete branches.

## 7. Scratch-directory cleanup policy

No source or run directory should be deleted before the consolidated branch is
built and validated.

The cleanup sequence is:

1. inventory every scratch directory;
2. identify its Git branch/commit and whether it contains uncommitted files;
3. preserve benchmark inputs, final outputs, logs and manifests;
4. remove reproducible build directories first;
5. archive obsolete source clones after commit verification;
6. remove duplicate run directories only after confirming which one is the
   authoritative result;
7. retain one production source clone and one controlled worktree per active
   development branch.

The visible scratch directory is a deployment/workspace problem, not the
canonical source history. GitHub must become the authoritative master, while
`/scratch` contains only builds, examples, submitted runs and retained results.
