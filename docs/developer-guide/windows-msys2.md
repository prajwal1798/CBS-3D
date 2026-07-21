# Windows and MSYS2 workflow

This workflow uses one shell consistently. Do not mix path syntaxes from PowerShell, WSL, Visual Studio terminals and MSYS2 within the same command sequence.

## Recommended terminal

Use:

```text
MSYS2 UCRT64
```

Verify tools:

```bash
which make
which g++
which cmake
which python

g++ --version
cmake --version
python --version
```

## Clone and enter the repository

```bash
git clone https://github.com/prajwal1798/CBS-3D.git
cd CBS-3D
```

## Make build

```bash
make clean
make USE_OPENMP=1 USE_MPI=0 USE_PETSC=0 -j4
```

Expected executable:

```text
cbs3dpp_si.exe
```

## CMake build

```bash
cmake -S . -B build-ucrt64 \
  -G "MSYS Makefiles" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCBS3D_ENABLE_OPENMP=ON \
  -DCBS3D_ENABLE_MPI=OFF \
  -DCBS3D_ENABLE_PETSC=OFF

cmake --build build-ucrt64 -j4
```

Do not reuse a build directory created by a different generator.

## Gmsh

When Gmsh is not in `PATH`, define its location once:

```bash
GMSH_EXE=/c/path/to/gmsh.exe
"$GMSH_EXE" -version
```

Generate a version-2 mesh when required by the active exporter:

```bash
"$GMSH_EXE" \
  input/<case>.geo \
  -3 \
  -format msh2 \
  -o input/<case>.msh
```

## Mesh conversion

The expected workflow is:

```text
.geo -> .msh -> .plt/.bco/.material/.matprop
```

Use the actual repository or controlled external exporter path. Do not document a guessed script location.

Typical pattern:

```bash
python /path/to/export.py \
  input/<case>.msh \
  -o input/<case> \
  --side-id-mode physical
```

The material writer consumes the `.msh` file, not `.plt`:

```bash
python /path/to/material_writer.py \
  input/<case>.msh \
  -o input/<case> \
  --fluid-name Fluid \
  --fluid-rho <rho> \
  --fluid-cp <cp> \
  --fluid-k <k> \
  --fluid-mu <mu>
```

## Run from the case directory

```bash
cd input
../cbs3dpp_si.exe <case>
```

Pass only the case base name.

Correct:

```bash
../cbs3dpp_si.exe cavity
```

Incorrect:

```bash
../cbs3dpp_si.exe cavity.par
../cbs3dpp_si.exe cavity.plt
```

Whether `input/cavity` is correct depends on the current working directory. Prefer entering the case directory and using the base name to remove ambiguity.

## Output

Inspect:

```text
output/<case>/<case>.pvd
output/<case>/<case>_residuals.csv
output/<case>/<case>_step_XXXXXXXX.vtu
```

## Path rules

- MSYS2 path: `/c/Users/...`
- Windows path: `C:\Users\...`
- Quote paths containing spaces.
- Use forward slashes in Bash commands.
- Do not paste PowerShell continuation characters into Bash.

## Safe mesh regeneration

Do not overwrite the accepted mesh until the exporter audit passes. Write a test mesh first:

```bash
"$GMSH_EXE" input/<case>.geo -3 -format msh2 \
  -o input/<case>_test.msh
```

Verify topology, boundary tags and material regions before replacing the controlled input.

## Common mistakes

| Mistake | Correction |
|---|---|
| Running `gmsh` when it is not in `PATH` | Use the verified full executable path |
| Supplying `.plt` to the material writer | Supply `.msh` |
| Passing a file extension to the solver | Pass the case base name |
| Mixing MSVC objects with UCRT64 objects | Clean and rebuild with one toolchain |
| Reusing a Visual Studio CMake directory with MSYS Makefiles | Use a new build directory |
| Overwriting accepted case files during exporter testing | Write `_test` outputs first |

## PETSc on Windows

PETSc toolchain compatibility on native Windows requires additional control and is not covered by this basic workflow. Use the verified project environment rather than combining arbitrary MSYS2, MS-MPI and prebuilt PETSc packages.
