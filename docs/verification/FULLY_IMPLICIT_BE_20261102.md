# CBS3D P1/TET4: fully implicit backward-Euler flow + thermal CHT (experimental)

**Branch:** `feat/fully-implicit-be-20261102` (based on `1338bff`).

**Status:** implementation committed; NOT compiled or validated on Sunbird2. Do not use for production forecasts until all gates below pass.

## Mathematical model

Physical time is `dt` in seconds. Let `u^n` be the velocity from the previous accepted physical step. At Picard iterate `k`, solve for corrected velocity `u^(k+1)`, predictor `u*^(k+1)`, and pressure `p^(k+1)`:

1. Predictor: `(M_L/dt) u* + A_u(w^k;dt) u = (M_L/dt)u^n`.
2. Pressure: `dt H p - B_vol u* = -B_boundary u^n`, with the sign convention used in the existing CBS Step-2 assembly.
3. Correction: `(M_L/dt)(u-u*) + G p = 0`.

For the steady fixed-point mode only, Step 1 instead uses
`(M_L/tau)(u*-u)+A_u(w^k;tau)u=0`, and the pressure boundary term multiplies corrected velocity.

The flow solve uses PETSc FGMRES/PCFIELDSPLIT(Schur), strong nodal boundary conditions, and Picard underrelaxation. The physical BE flow step is followed by one BE energy solve at the **same** `dt` using the existing graph-viscosity low-order CHT discretisation. Constant properties, prescribed inlet mass flow, and no thermal feedback permit this one-way coupling; it is **not** a monolithic four-field thermo-fluid Newton method.

Only P1 tetrahedra, laminar flow, and the existing distributed constraints are targeted.

## Reproducible initial test

On Sunbird2, in a shell with the CBS3D environment loaded:

```bash
cd "$CBS3D_CODE"
git status --short
git fetch origin feat/fully-implicit-be-20261102
git switch feat/fully-implicit-be-20261102
git pull --ff-only
bash jobs/sunbird/templates/build_cbs3d.sh
bash jobs/sunbird/templates/blanket.sh be-smoke
bash jobs/sunbird/templates/blanket.sh status
```

The `be-smoke` command submits five physical BE flow/thermal steps on the existing 40-rank partition set. It does not modify or restart the existing production run.

**First validate PETSc/MPI runtime.** Sunbird2 job 30020 failed in `MatCreateAIJ` with PETSc error 98. In PETSc 3.26, this is `PETSC_ERR_MPI` (a general MPI error), **not** a matrix preallocation error. Check the full PETSc stack trace, matching `libmpi` dependencies in the application and `libpetsc`, `mpirun`/Slurm PMI integration, and that PETSc was initialized collectively. A successful compilation alone does not validate this.

## Acceptance gates

1. **MPI:** clean 2-rank and 40-rank matrix creation and KSP solve; no rank hangs, no MPI mismatch. Check `ldd` and `PETSC_OPTIONS=-on_error_attach_debugger` only in a controlled debug test, never in a production batch.
2. **Algebra:** report `||b-Ax||/||b||` for each Picard iteration. Require KSP convergence **and** the independently recomputed residual. Note: finite-Picard change tolerances alone do not prove the nonlinear residual is small.
3. **Continuity:** independently reassemble weak continuity and evaluate nodal and volume-norm measures. Avoid relying only on `DivRMS`.
4. **Temporal order:** use manufactured unsteady Stokes flow and thermal conduction solutions with dt, dt/2, dt/4 and fixed fine mesh; demonstrate first-order convergence where spatial error is negligible.
5. **Conservation:** perform integrated fluid mass balance, fluid momentum residual, total thermal energy budget, and conjugate heat-flux balance. Verify correct inlet/outlet/nuclear heating and BC 532 terms.
6. **Physical bounds:** test thermal minimum and maximum at each step, **without clipping**. A previous explicit-flow/BE-CHT run violated the optional 298.15 K audit floor on step 1; diagnose the material/source/Dirichlet/operator contribution before claiming physical admissibility.
7. **MPI independence:** compare fields and balances at 2, 4, and 40 ranks on identical meshes; verify global timestep, ownership, and halo reductions.
8. **Efficiency:** measure matrix setup/preallocation, Schur preconditioning, KSP time, Picard count, and memory per rank. The current prototype rebuilds matrices/PCs every physical time step and is not an optimized 50k-step production solver.

## Explicit limitations

- MPI/PETSc build and validation cannot be performed from the GitHub connector.
- The current predictor/pressure/correction blocks inherit the existing CBS stabilisation and equal-order formulation, not a newly derived inf-sup-stable mixed method.
- The energy method is backward Euler but spatially low-order graph-viscosity stabilised; do not call it high-order.
- Time-step selection in the transient production loop still uses the existing communicator-wide timestep policy, and a BE-specific timestep-independence audit is required.
- Large-scale time integration is too expensive unless matrix graph and preconditioner reuse is implemented.
