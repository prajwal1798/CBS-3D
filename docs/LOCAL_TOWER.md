# Local Tower Workflow

## Hardware

| Item | Value |
| --- | --- |
| CPU | Intel Core Ultra 7 265F |
| Cores / logical processors | 20 / 20 |
| RAM | 16 GB |
| OS | Windows 11 |
| Compiler | Visual Studio 18 Community, MSVC 19.51 |
| CMake | 4.3.4 |
| Native PETSc | MSYS2 UCRT64 GCC 16.2 + PETSc 3.25.1 |

The tower is a shared-memory machine.  There is no multi-node requirement for
the current development loop.  The 8-part flat-plate mesh is valuable for
testing the partition format and the MPI path, but it cannot be run as eight
ranks with the MSYS2 PETSc package because that package is an MPIUNI build.

## Native MSVC build

```powershell
cmake -S . -B build\local-serial-openmp -G "Visual Studio 18 2026" -A x64 `
      -DCBS3D_ENABLE_OPENMP=ON -DCBS3D_ENABLE_MPI=OFF -DCBS3D_ENABLE_PETSC=OFF
cmake --build build\local-serial-openmp --config Release --parallel 20
```

This is the fastest way to compile the solver.  It is the build used for the
wall-treatment and SA element tests.

## Local PETSc build

MSYS2 was installed at `C:\msys64`.  The required packages are:

```powershell
C:\msys64\usr\bin\bash.exe -lc "pacman -S --needed --noconfirm ucrt64/mingw-w64-ucrt-x86_64-gcc ucrt64/mingw-w64-ucrt-x86_64-cmake ucrt64/mingw-w64-ucrt-x86_64-ninja ucrt64/mingw-w64-ucrt-x86_64-petsc ucrt64/mingw-w64-ucrt-x86_64-hypre"
```

The CMake preset carries the environment:

```powershell
cmake --preset local-msys2-petsc
cmake --build --preset build-local-msys2-petsc
```

Equivalent manual configuration:

```powershell
$env:PATH = "C:\msys64\ucrt64\bin;$env:PATH"
$env:PETSC_DIR = "C:/msys64/ucrt64"
$env:PETSC_ARCH = "dto"

cmake -S . -B build\local-msys2-petsc -G Ninja `
      -DCMAKE_BUILD_TYPE=Release `
      -DCMAKE_C_COMPILER=C:/msys64/ucrt64/bin/gcc.exe `
      -DCMAKE_CXX_COMPILER=C:/msys64/ucrt64/bin/g++.exe `
      -DCBS3D_ENABLE_OPENMP=ON `
      -DCBS3D_ENABLE_MPI=OFF `
      -DCBS3D_ENABLE_PETSC=ON `
      -DCBS3D_BUILD_TESTS=ON
```

`dto` means double precision, MPIUNI, OpenMP.  PETSc is installed under
`C:\msys64\ucrt64\lib\petsc\dto`; the CMake logic understands this packaged
layout as well as a normal `${PETSC_DIR}/${PETSC_ARCH}` install.

### PETSc packaging quirk

The MSYS2 PETSc is compiled without C++ support, so its `petscconf.h` omits
`PETSC_FUNCTION_NAME_CXX` and `PETSC_CXX_RESTRICT`.  `include/cbs/linalg/PetscHeaders.hpp`
defines them before including PETSc and is the only PETSc include point in the
solver.  Do not include `<petscksp.h>` directly.

### hypre / BoomerAMG

The packaged PETSc does **not** define `PETSC_HAVE_HYPRE`.  The solver therefore
takes the `PCGAMG` fallback branch.  hypre is installed as a separate library,
but PETSc must be rebuilt with it to expose `PCHYPRE`:

```bash
./configure --with-hypre=1 --with-hypre-dir=/ucrt64 \
            --with-debugging=0 --with-openmp=1 \
            COPTIMFLAGS="-O3" CXXOPTIMFLAGS="-O3"
make -j20
```

This is optional.  GAMG already gives a working 284k-cell solve.

## Running the flat plate

The PETSc smoke case is `cases\flat_plate_petsc`:

```powershell
$env:PATH = "C:\msys64\ucrt64\bin;$env:PATH"
$env:OMP_NUM_THREADS = "20"
Set-Location cases\flat_plate_petsc
..\..\build\local-msys2-petsc\cbs3dpp_si.exe flat_plate_petsc
```

The case uses the full 284,400-tetrahedron mesh and three CBS iterations.
Measured on the tower:

| Iteration | Pressure CG iterations | CG relative L2 | Pressure step | SA step | Total |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1 | 256 | 1.15e-03 | 1.35 s | 0.23 s | 1.60 s |
| 2 | 450 | 2.91e-06 | 1.96 s | 0.23 s | 2.22 s |
| 3 | 295 | 2.91e-06 | 1.33 s | 0.24 s | 1.60 s |

The wall-clock for the three-iteration run was 5.4-5.8 s.  The pressure solve is
about 82-88% of the iteration; SA transport is about 10-16%.

The native CG baseline is `cases\flat_plate_cg`.  It uses the same mesh and
controls with `Pressure Field Solver = 1`.

### Important `.par` detail

The timestep line is:

```
ntime transient_on dtfixed dtfix iwrite
```

`dtfixed` is a flag (`-1`, `0`, `1`), not the timestep.  The timestep is
`dtfix`.  A value of `100000` in the `dtfix` column produces `dt=1.0e+05` and
makes the pressure system unnecessarily hard; the flat-plate cases use
`1.0e-5`.

## Threading

Set `OMP_NUM_THREADS` explicitly.  The default on Windows may be the number of
logical processors, but being explicit makes runs reproducible:

```powershell
$env:OMP_NUM_THREADS = "20"
```

The assembly colouring is the main OpenMP consumer.  The PETSc GAMG setup and
solve are partly serial; increasing threads beyond the useful assembly width
will not scale linearly.  Measure with the per-section profiler printed at the
end of each iteration.

## Tests

The MSVC test suite is the reliable way to execute the verification binaries on
this machine:

```powershell
cmake --preset tests-msvc-release
cmake --build --preset build-tests-msvc-release
ctest --preset test-tests-msvc-release
```

The MSYS2 build also has the tests, but Windows code-integrity policy can block
freshly built unsigned executables.  When that happens, run the MSVC equivalent
or rebuild; it is an operating-system policy, not a solver failure.

## What is still 20-core-only

- The 8-part `cases\flat_plate_p8` mesh is validated but needs a real MPI PETSc
  build to run.
- The Sunbird `jobs/sunbird` scripts remain the production path for MPI/hypre
  scaling.
- A local hypre-enabled PETSc rebuild is the next step if BoomerAMG is needed
  for comparison with GAMG.
