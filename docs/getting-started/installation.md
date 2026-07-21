# Installation

## Supported build configurations

CBS3D++ supports four principal capability combinations:

| Build | OpenMP | MPI | PETSc | Intended use |
|---|---:|---:|---:|---|
| Reference | optional | no | no | Baseline solver and regression |
| PETSc/OpenMP | yes | no | yes | Current large-case pressure-solver path |
| MPI bootstrap | optional | yes | no | Distributed infrastructure development |
| MPI/PETSc development | optional | yes | yes | Integration work; not yet a released distributed CFD path |

## Requirements

### Core

- C++20 compiler
- CMake 3.16 or newer, or GNU Make
- Python 3 for repository utilities

### Optional

- OpenMP for shared-memory element assembly
- MPI for distributed-development builds
- PETSc for KSP pressure solution
- HYPRE when BoomerAMG is selected through PETSc
- Gmsh for mesh generation
- ParaView for VTU/PVD post-processing

## Toolchain compatibility

PETSc, MPI and the application must use an ABI-compatible compiler and MPI family. A configuration can compile successfully and still fail during linking or MPI initialisation when PETSc was built against a different MPI implementation.

Record at minimum:

```text
compiler version
MPI implementation and version
PETSc version and configuration
PETSC_DIR
PETSC_ARCH, when used
CMake or Make command
Git commit
```

## CMake builds

### Serial/OpenMP

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCBS3D_ENABLE_OPENMP=ON \
  -DCBS3D_ENABLE_MPI=OFF \
  -DCBS3D_ENABLE_PETSC=OFF

cmake --build build -j
```

The CMake target is:

```text
cbs3dpp_si
```

### MPI development build

```bash
cmake -S . -B build-mpi \
  -DCMAKE_BUILD_TYPE=Release \
  -DCBS3D_ENABLE_OPENMP=ON \
  -DCBS3D_ENABLE_MPI=ON \
  -DCBS3D_ENABLE_PETSC=OFF

cmake --build build-mpi -j
```

### PETSc build

Set the PETSc environment first:

```bash
export PETSC_DIR=/path/to/petsc
# export PETSC_ARCH=arch-name   # only when the PETSc installation uses it
```

Then configure:

```bash
cmake -S . -B build-petsc \
  -DCMAKE_BUILD_TYPE=Release \
  -DCBS3D_ENABLE_OPENMP=ON \
  -DCBS3D_ENABLE_MPI=OFF \
  -DCBS3D_ENABLE_PETSC=ON

cmake --build build-petsc -j
```

## Make builds

The Makefile uses separate object directories for incompatible build variants.

### Serial/OpenMP

```bash
make USE_OPENMP=1 USE_MPI=0 USE_PETSC=0 -j
```

Executable:

```text
cbs3dpp_si
```

### PETSc/OpenMP

```bash
make USE_OPENMP=1 USE_MPI=0 USE_PETSC=1 \
  PETSC_DIR="$PETSC_DIR" -j
```

Executable:

```text
cbs3dpp_si_petsc
```

### MPI bootstrap

```bash
make USE_OPENMP=1 USE_MPI=1 USE_PETSC=0 \
  MPI_CXX=mpicxx -j
```

Executable:

```text
cbs3dpp_si_mpi
```

### MPI/PETSc development

```bash
make USE_OPENMP=1 USE_MPI=1 USE_PETSC=1 \
  MPI_CXX=mpicxx \
  PETSC_DIR="$PETSC_DIR" -j
```

Executable:

```text
cbs3dpp_si_mpi_petsc
```

## Verify the binary

On Linux:

```bash
ldd ./cbs3dpp_si_petsc | grep -Ei 'petsc|hypre|mpi|gomp'
```

For a PETSc/HYPRE/OpenMP build, the output should resolve the corresponding libraries without `not found` entries.

## Clean rebuilds

When changing compiler, MPI, PETSc, OpenMP flags or build mode, use a clean build directory or run:

```bash
make clean
```

Do not reuse object files produced by an incompatible toolchain.

## Next step

Proceed to the [quick-start guide](quickstart.md), then consult the [case-file reference](../user-guide/case-files.md).
