# CBS3D++_SI — MPI Stage 2

## Distributed ownership and weak-form lumped-operator assembly

This package is the first **algorithmic** domain-decomposition change to the
CBS3D++_SI solver. It replaces the Stage 1 metadata-only path with a genuine
finite-element assembly distributed over MPI ranks.

It does **not** yet advance CBS Steps 1–4. The distributed pressure operator and
global Conjugate Gradient method are the next stage and must be validated before
the complete solver loop is enabled.

---

## Mathematical basis

The semi-implicit CBS intermediate-momentum step contains the finite-element
mass matrix

\[
M_u\,\Delta U^* = \Delta t\,R_u,
\qquad
M_u = \int_\Omega N^T N\,d\Omega.
\]

This is the Step 1 matrix form and mass-matrix definition appearing in the CBS
formulation, including Eq. (3.49)–(3.50) in *A General Algorithm for the CBS
Scheme*.

For a four-node linear tetrahedron with volume \(V_e\),

\[
M_{aa}^{(e)}=\frac{V_e}{10},\qquad
M_{ab}^{(e)}=\frac{V_e}{20},\qquad
M_{L,a}^{(e)}=\frac{V_e}{4}.
\]

For conjugate heat transfer, the row-summed thermal capacitance is

\[
C_{L,a}^{(e)}=(\rho c_p)_e\frac{V_e}{4}.
\]

Under domain decomposition, the global nodal sum is partitioned between ranks:

\[
g_i=\sum_{e\ni i}g_i^{(e)}
   =\sum_r\sum_{e\in\Omega_r,\,e\ni i}g_i^{(e)}.
\]

Therefore, every shared-node assembly follows exactly this sequence:

1. Each rank scatters contributions from its owned tetrahedra.
2. Ghost-node contributions are sent to the unique owner.
3. The owner adds all neighbour contributions.
4. The owner broadcasts the completed value to all ghosts.

The partition boundary is not a physical boundary and receives no wall,
adiabatic, heat-flux or interface condition.

---

## Main structural changes

### `DistributedLayout`

Centralises and validates:

- local-to-global node and element numbering;
- one unique owner for every nodal unknown;
- owned and ghost node lists;
- neighbour SEND and RECV maps;
- runtime MPI rank and partition consistency.

### `HaloExchange`

Provides two deliberately separate operations:

- `accumulateToOwners(...)`: ghost-to-owner addition after local assembly;
- `updateGhosts(...)`: owner-to-ghost overwrite after an owner update.

`synchroniseAssembledSum(...)` performs both operations in the correct order.

### Split mass preprocessing

`Preprocess::massMatrix()` is split into:

- `assembleMassMatrixContributions()`;
- `finaliseMassMatrix()`.

This prevents the solver from inverting an incomplete partition-local diagonal
at shared nodes.

### New persistent nodal coefficient

`CBSStateSI::Cdiag_real` stores the assembled thermal-capacitance diagonal before
its inverse is formed in `elcoe2p`.

---

## Files to copy into the existing project

Copy the package over the existing project root while preserving its directory
structure. It contains only new or changed files.

```bash
cd /scratch/s.2337862/CBS3D++_SI
cp -a /path/to/CBS3D_MPI_Stage2/. .
```

Do not copy the package into a nested subdirectory of the solver.

---

## Build on Sunbird

```bash
cd /scratch/s.2337862/CBS3D++_SI
chmod +x build_stage2_sunbird.sh
./build_stage2_sunbird.sh
```

The build uses:

- GNU 12.1.0;
- OpenMPI 4.1.6;
- CMake 3.31.10;
- OpenMP enabled;
- PETSc disabled.

---

## Run on two physical nodes

```bash
sbatch run_stage2_40_2nodes.slurm
```

After completion:

```bash
sacct -j <JOB_ID> \
  --format=JobID,JobName,NNodes,NTasks,State,ExitCode,Elapsed

cat cbs_stage2_2n_<JOB_ID>.out
cat cbs_stage2_2n_<JOB_ID>.err
```

The decisive output is:

```text
owner-to-ghost field update        : PASS
ghost-to-owner FE accumulation     : PASS
distributed lumped operators       : PASS
RESULT                              : PASS
```

The following integral errors should be near round-off:

```text
relative volume error
relative thermal error
```

---

## Local deterministic tests

```bash
./tests/run_unit_tests.sh
```

Expected:

```text
partition layout test: PASS
single-tetra mass test: PASS
```

---

## Validation already performed before packaging

- Complete serial/OpenMP project build with CMake: PASS.
- Complete serial project build with the supplied Makefile: PASS.
- MPI-enabled source syntax and full link check using an MPI interface stub:
  PASS.
- One-tetrahedron analytical mass/capacitance test: PASS.
- Two-partition ownership-layout test: PASS.
- Full 2,383,581-tetrahedron blanket mesh serial preprocessing test: PASS.
  The integrated volume was `3.0000000000000001e-03` from both elements and
  lumped nodes; the relative discrepancies were `2.67e-17` for volume and
  `5.27e-17` for thermal capacitance.

A real two-node runtime test must still be executed on Sunbird because this
container has no MPI runtime or access to the exported 40-partition mesh.

---

## Next coding stage

The next package should implement a matrix-free distributed pressure operator
and global Conjugate Gradient method:

1. update pressure-vector ghosts;
2. apply local element pressure-Laplacian terms;
3. accumulate shared-node matrix-vector contributions to owners;
4. compute dot products and residual norms over owned pressure nodes only;
5. use `MPI_Allreduce` for every global CG scalar;
6. enforce pressure-outlet constraints only on owner ranks.

Only after this pressure stage matches the one-rank solution should CBS Steps 1,
3 and 4 be enabled in the distributed solver loop.
