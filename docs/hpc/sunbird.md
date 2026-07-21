# Sunbird HPC workflow

This page defines the public, sanitised Sunbird workflow. Usernames, project accounts and private filesystem paths must not be committed. Replace placeholders locally.

## Scope

The current accepted large-case execution model is:

```text
one solver process
OpenMP threaded element loops
PETSc/HYPRE pressure solve when enabled
Slurm compute allocation
scratch working directory
```

MPI mesh-distribution development is separate from the accepted production CFD path.

## Directory policy

Use high-performance scratch for builds and numerical runs:

```text
/scratch/<username>/CBS3D++_SI
```

Do not treat scratch as archival storage. Retain important source revisions, run records and validation outputs in backed-up storage or an approved research-data archive.

Do not run production cases or large PETSc builds directly on a login node.

## Environment discovery

Start clean:

```bash
module purge
module list
```

Discover live modules rather than copying stale names:

```bash
module spider GCC
module spider OpenMPI
module spider PETSc
module spider CMake
```

Select one ABI-compatible compiler/MPI/PETSc family. Record the exact loaded modules:

```bash
module list 2>&1 | tee logs/modules.txt
which g++
which mpicxx
mpicxx --showme 2>/dev/null || mpicxx -show
```

## Controlled environment file

Create an untracked site-local file such as:

```text
sunbird_modules.sh
```

Example structure:

```bash
#!/usr/bin/env bash
set -euo pipefail

module purge
module load <compiler-module>
module load <mpi-module>
module load <petsc-module-or-local-environment>

export PETSC_DIR=<petsc-prefix>
# export PETSC_ARCH=<petsc-arch>  # only when required
```

Do not commit personal paths or credentials.

## Preflight

Before building:

```bash
source ./sunbird_modules.sh

g++ --version
mpicxx --version
cmake --version
make --version

printf 'PETSC_DIR=%s\n' "${PETSC_DIR:-unset}"
printf 'PETSC_ARCH=%s\n' "${PETSC_ARCH:-unset}"
```

Verify PETSc headers and libraries:

```bash
test -f "$PETSC_DIR/include/petsc.h"
find "$PETSC_DIR" -name 'libpetsc*' -print
```

When HYPRE is required, confirm it was enabled in PETSc and resolves at runtime.

## Build the current PETSc/OpenMP path

```bash
make clean

make USE_PETSC=1 USE_MPI=0 USE_OPENMP=1 \
  SERIAL_CXX=mpicxx \
  PETSC_DIR="$PETSC_DIR" \
  -j 8
```

Verify:

```bash
ls -lh cbs3dpp_si_petsc
ldd ./cbs3dpp_si_petsc | grep -Ei 'petsc|hypre|mpi|gomp'
```

No library should resolve to `not found`.

## Slurm template

```bash
#!/usr/bin/env bash
#SBATCH --job-name=cbs3d
#SBATCH --account=<project-account>
#SBATCH --partition=<compute-partition>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err

set -euo pipefail

cd /scratch/<username>/CBS3D++_SI
source ./sunbird_modules.sh

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OMP_PROC_BIND=close
export OMP_PLACES=cores

CASE=<case-base-name>

cd input
srun --ntasks=1 \
  --cpus-per-task="${SLURM_CPUS_PER_TASK}" \
  --cpu-bind=cores \
  ../cbs3dpp_si_petsc "$CASE" \
  -pc_type hypre \
  -pc_hypre_type boomeramg \
  -ksp_converged_reason
```

Use the live site-recommended `--mpi` option when explicitly required. Do not guess it from another cluster.

## Mandatory run record

Every accepted build/run should retain:

```text
date and time
Git commit
input checksums
loaded modules
compiler version
MPI version
PETSc configuration
build command
executable checksum
case name
Slurm job ID
account and partition
nodes, tasks and threads
wall time and MaxRSS
stdout and stderr
residual history
VTU/PVD output
acceptance result
```

Suggested commands:

```bash
git rev-parse HEAD
sha256sum cbs3dpp_si_petsc
sha256sum input/<case>.plt input/<case>.bco input/<case>.par
sacct -j "$SLURM_JOB_ID" \
  --format=JobID,JobName,Partition,AllocCPUS,Elapsed,State,ExitCode,MaxRSS
```

## Distributed-mesh development

The DMPlex/partition acceptance executable is intended to validate:

- distributed mesh creation;
- ownership and ghost points;
- boundary labels;
- material labels;
- global topology totals.

Passing these checks does not mean that CBS Steps 1–4 are distributed. Follow the gates in [parallel status](../parallel/status.md).

## Common failures

| Symptom | Likely cause | Response |
|---|---|---|
| `mpicxx: command not found` | MPI module missing | Purge modules, rediscover and load one compatible family |
| `PETSC_DIR is not set` | PETSc environment incomplete | Inspect module output or controlled local prefix |
| `petsc.h` missing | Wrong prefix or `PETSC_ARCH` | Locate `petscconf.h` and correct the installation layout |
| `libpetsc.so` missing at runtime | Runtime library path/rpath problem | Re-source environment and inspect `ldd` |
| MPI initialisation failure | PETSc and executable use different MPI libraries | Compare wrapper and PETSc configuration; rebuild consistently |
| Unknown partitioner | PETSc lacks the requested package | Use an available diagnostic partitioner or rebuild PETSc |
| Job timeout | Wall-time request too short | Inspect progress and resubmit with justified resources |

## Security and publication

Never commit:

- personal cluster usernames;
- private SSH hosts or keys;
- account/project identifiers not intended for publication;
- home-directory software paths;
- proprietary case data;
- credentials or access tokens.
