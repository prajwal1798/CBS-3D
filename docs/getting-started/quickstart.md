# Quick start

This guide shows the minimum portable workflow. Site-specific paths, Slurm settings and private development cases are intentionally excluded.

## 1. Clone

```bash
git clone https://github.com/prajwal1798/CBS-3D.git
cd CBS-3D
```

## 2. Build the reference solver

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCBS3D_ENABLE_OPENMP=ON \
  -DCBS3D_ENABLE_MPI=OFF \
  -DCBS3D_ENABLE_PETSC=OFF

cmake --build build -j
```

The executable is normally located at:

```text
build/cbs3dpp_si
```

Multi-configuration generators may place it under `build/Release/`.

## 3. Prepare a case

The solver argument is the **case base name**, not an individual filename. For:

```text
cases/cavity/cavity.plt
cases/cavity/cavity.bco
cases/cavity/cavity.par
```

run with:

```text
cases/cavity/cavity
```

Material-enabled cases also require:

```text
cases/cavity/cavity.material
cases/cavity/cavity.matprop
```

All files must share the same base name. See [case files](../user-guide/case-files.md).

## 4. Run

```bash
./build/cbs3dpp_si cases/cavity/cavity
```

On Windows with a Visual Studio generator:

```powershell
.\build\Release\cbs3dpp_si.exe cases\cavity\cavity
```

## 5. Check successful startup

Confirm that the startup report shows the intended:

```text
case name
mesh dimensions
simulation mode
material mode
boundary classification
pressure solver
time-step controls
```

Stop the run when any reported count or solver mode disagrees with the intended case.

## 6. Inspect output

Typical output contains:

```text
<case>.pvd
<case>_step_XXXXXXXX.vtu
<case>_residuals.csv
```

Open the `.pvd` file in ParaView and inspect velocity, pressure, temperature and any turbulence/debug fields written by the active build.

## PETSc pressure path

Build:

```bash
export PETSC_DIR=/path/to/petsc

make USE_PETSC=1 USE_MPI=0 USE_OPENMP=1 \
  PETSC_DIR="$PETSC_DIR" -j
```

Run:

```bash
./cbs3dpp_si_petsc cases/cavity/cavity \
  -pc_type hypre \
  -pc_hypre_type boomeramg \
  -ksp_converged_reason
```

The `.par` pressure-solver selection and executable capability must agree.

## MPI development path

The MPI target currently supports distributed-development infrastructure. It must not be used as evidence of a fully distributed CBS solver until the acceptance gates in [parallel status](../parallel/status.md) pass.

## Public example status

A compact, redistributable first-run case is required before the first formal release. Large research meshes and legacy-protected benchmark packages should not be used as the only onboarding route.
