# Legacy regression protocol

## Purpose

Legacy regression determines whether CBS3D++ reproduces the accepted behaviour of the corresponding Fortran CBS solver. It is a redevelopment-control activity, not a substitute for independent physical validation.

## Legal boundary

Use only legacy executables, source, documentation and case data that the project is authorised to access. Zeta Computational Resources states that CBSFlow source is copyrighted and may not be redistributed without written permission.

Regression artefacts stored publicly must therefore be limited to:

- newly created CBS3D++ inputs;
- numerical data whose redistribution is authorised;
- scripts written for this project;
- summary norms and plots derived lawfully;
- links and citations to protected source material.

Do not commit protected legacy source or copied manuals merely to make a test self-contained.

## Comparison hierarchy

### Level 1 — input interpretation

Confirm identical interpretation of:

```text
NELEM, NPOIN and NBOUN
connectivity
coordinates
boundary mappings
material IDs and properties
initial conditions
solver options
time-step controls
```

A result comparison is invalid when the two executables did not solve the same case.

### Level 2 — preprocessing quantities

Compare:

- tetrahedral determinants and volumes;
- shape-function gradients;
- lumped nodal mass;
- pressure stiffness coefficients;
- boundary parent faces and normals;
- fluid, solid and interface masks;
- pressure-reference nodes.

### Level 3 — single-step quantities

For a small deterministic mesh, export and compare after each CBS stage:

```text
Step 1 predictor velocity and RHS
Step 2 pressure RHS and pressure
Step 3 corrected velocity
Step 4 temperature/scalar RHS and field
```

This localises differences more effectively than comparing only a final converged solution.

### Level 4 — iterative history

Compare:

- time step;
- pressure iterations and residuals;
- velocity/pressure/temperature convergence histories;
- global conservation quantities;
- output sampling times.

### Level 5 — final fields and engineering quantities

Compare:

- velocity components;
- pressure;
- temperature;
- turbulence variables when applicable;
- mass-flow balance;
- drag/lift or pressure drop;
- wall heat flux and energy balance.

## Controlled test procedure

1. Select a case supported by both solvers.
2. Freeze and checksum all input files.
3. Record both executable versions and build environments.
4. Disable optional algorithmic differences where possible.
5. Run one iteration and compare stage-level arrays.
6. Run a short fixed number of iterations.
7. Run to a common convergence criterion.
8. Extract fields at identical nodes or map them conservatively.
9. Compute predefined norms.
10. Archive the report with the CBS3D++ commit.

## Suggested tolerance classes

Tolerance values must be tailored to the quantity and algorithm, but the report should distinguish:

| Class | Meaning |
|---|---|
| Bitwise | Identical binary values; normally limited to tightly controlled serial tests |
| Round-off | Differences consistent with operation ordering and compiler arithmetic |
| Solver-equivalent | Different iterative paths converge to the same accepted algebraic solution |
| Discretisation-equivalent | Engineering fields agree within the predefined numerical tolerance |

Do not use a loose final-field tolerance to conceal a sign, indexing or boundary-condition error visible in stage-level arrays.

## Recommended report format

```text
Case:
CBS3D++ commit:
Legacy executable/source identifier:
Input checksums:
Compiler/toolchain:
Solver modes:
Reference pressure treatment:
Time-step policy:
Pressure tolerance:
Iterations compared:

Mesh interpretation: PASS/FAIL
Preprocessing comparison: PASS/FAIL
Step 1: PASS/FAIL
Step 2: PASS/FAIL
Step 3: PASS/FAIL
Step 4: PASS/FAIL/NA
Final field norms: table
Conservation: table
Known intentional differences:
Conclusion:
```

## Intentional differences register

A modernisation may deliberately change:

- storage layout;
- linear-solver implementation;
- operation ordering;
- output format;
- parallel decomposition;
- defensive input validation.

It must not silently change the mathematical formulation. Every intentional numerical difference requires:

1. a written reason;
2. the affected equation or algorithm;
3. an independent verification case;
4. an update to the project lineage documentation.
