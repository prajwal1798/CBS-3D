# Contributing to CBS3D++

CBS3D++ is research software. A change is acceptable only when its numerical meaning, implementation impact and verification evidence are clear.

## Before opening a change

- Search existing issues and pull requests.
- Identify whether the change is numerical, infrastructural, documentary or case-specific.
- Keep unrelated changes in separate branches.
- Do not commit private cluster paths, credentials, proprietary meshes or restricted legacy material.

## Branches

Use a descriptive branch name:

```text
feature/<capability>
fix/<problem>
docs/<topic>
validation/<case>
parallel/<phase>
```

## Build requirements

At minimum, verify the serial/OpenMP configuration:

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCBS3D_ENABLE_OPENMP=ON \
  -DCBS3D_ENABLE_MPI=OFF \
  -DCBS3D_ENABLE_PETSC=OFF

cmake --build build -j
```

When the change affects MPI or PETSc, test the corresponding build with a compatible toolchain.

## Numerical pull requests

A numerical change must document:

1. physical or mathematical purpose;
2. governing equation and source/reference;
3. discrete finite-element form;
4. coefficients and units;
5. input and output arrays;
6. boundary-condition behaviour;
7. material-domain behaviour;
8. OpenMP race-safety strategy;
9. MPI ownership/communication implication;
10. regression or verification evidence.

Do not describe a numerical change only as a refactor when operation ordering, boundary sequencing, solver tolerance, accumulation or floating-point reductions change.

## Required checks by change type

| Change | Minimum evidence |
|---|---|
| Documentation only | MkDocs strict build |
| Build-system change | Clean configure/build for affected variants |
| Input/parser change | Accepted old case plus new parser test |
| Assembly change | Single-step array comparison and final regression |
| Pressure solver | Operator/RHS comparison and convergence record |
| Boundary condition | Minimal case exercising the boundary |
| Multiphysics | Material/interface audit and conservation check |
| Turbulence | Model-term diagnostics and benchmark comparison |
| MPI | One-rank equivalence plus multi-rank ownership checks |
| Performance | Correctness first, then measured profile/scaling data |

## Validation data

A validation contribution should include:

```text
case description
reference citation
mesh statistics
input checksums
solver configuration
extraction script
comparison script
acceptance tolerances
summary table
representative plot
```

Do not commit large binary results without a documented storage policy. Prefer scripts and compact extracted data; archive large datasets separately and link them by DOI or stable record.

## Coding conventions

- C++20.
- Preserve module boundaries under `include/cbs/` and `src/`.
- Use existing solver types and containers consistently.
- Preserve one-based scientific indexing unless a complete migration is approved.
- Keep numerical kernels explicit and auditable.
- Document units and indexing assumptions.
- Compile with warnings enabled.
- Avoid hidden global state and undocumented side effects.
- Keep rank ownership explicit in distributed code.

## Documentation

Update documentation in the same pull request when changing:

- equations or algorithm stages;
- capability status;
- build options;
- executable names;
- case-file grammar;
- boundary IDs or semantics;
- material conventions;
- validation results;
- parallel acceptance status.

The root README remains a landing page. Detailed material belongs under `docs/`.

Build documentation locally:

```bash
python -m pip install -r docs/requirements.txt
mkdocs build --strict
```

## Pull-request description

Use this structure:

```text
Purpose

Numerical/software change

Files and interfaces affected

Validation and regression

Parallel implications

Known limitations
```

## Legacy material

The legacy CBSFlow source and documentation may be subject to Zeta Computational Resources copyright and redistribution restrictions. Do not add translated or copied legacy material unless the project owner has documented permission and the repository licence is compatible.

## Review standard

A pull request can be merged when:

- its scope is coherent;
- the build checks pass;
- numerical claims are supported;
- capability status remains accurate;
- restricted data are absent;
- documentation is updated;
- unresolved risks are stated explicitly.
