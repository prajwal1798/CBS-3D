# CBS3D++_SI Distributed-Memory MPI Migration Audit

## 1. Scope of this audit

This audit was produced from the complete legacy archive:

- `MPIMainGG07SI.f90`
- all files under `Program/`
- the legacy Makefile and run scripts
- the serial fake-MPI compatibility files
- runtime option files

The purpose is not to copy the old MPI implementation. The purpose is to identify the numerical communication semantics that must be preserved while replacing the old manual communication layer with a modern PETSc-native distributed architecture.

---

## 2. Important qualification about the archive

The driver declares:

```fortran
parameter (ndim=2,nep=ndim+1,nsid=ndim+1,nsidp=ndim)
```

in `MPIMainGG07SI.f90`, even though the driver comments describe linear tetrahedral elements and many routines contain explicit 3D branches.

Therefore, this archive is a solver framework capable of 2D and 3D operation, but the saved driver snapshot is configured for a 2D run. The 3D branches in:

- `getgeom123linear.f90`
- `getface123linear.f90`
- `pstiff.f90`
- `multi4nPdiag.f90`
- the CBS element routines

remain directly relevant to CBS3D++_SI.

The modern implementation must use the 3D tetrahedral branches only.

---

## 3. Legacy execution order

The driver sequence is:

1. Initialise MPI and obtain rank and communicator size.
2. Define custom contiguous MPI datatypes for vectors and geometry records.
3. Rank zero reads runtime options.
4. Broadcast runtime options to all ranks.
5. Each rank opens its own partitioned mesh and communication files.
6. Each rank reads its local elements, nodes, boundary faces and global ids.
7. Each rank reads its shared-node communication map.
8. Compute local geometry, faces and boundary sets.
9. Assemble and synchronise characteristic lengths.
10. Assemble and synchronise lumped mass.
11. Assemble and synchronise the pressure diagonal.
12. Select one unique rank for every duplicated shared node.
13. Enter the CBS loop:
    - global time-step reduction
    - Step 1 local element assembly
    - shared-node residual accumulation
    - boundary conditions
    - Step 2 pressure RHS assembly
    - shared-node RHS accumulation
    - distributed Pressure CG
    - Step 3 local element assembly
    - shared-node residual accumulation
    - boundary conditions
    - global convergence reductions
    - rank-local output
14. Finalise MPI.

The central calls in the driver occur in this order:

```text
findlen / findlen_e
getmassmat
pstiff
converge_poin_alloc*
alotim5
step1gg
Step2Semi
step_3gg
cnvrg_chck
outputG1
```

---

## 4. Legacy partition data model

### 4.1 The solver does not partition the mesh at runtime

Each rank opens pre-generated files with a rank suffix:

```text
<case>_1.plt
<case>_2.plt
...
<case>_1.com
<case>_2.com
...
```

This is controlled by `input_file_open.f90` and `filemaker2.f90`.

The archive does not contain the partition-generator source. Therefore, the exact graph-partitioning algorithm used to create the old rank-local files is not available and should not be reconstructed.

### 4.2 Rank-local mesh

For MPI runs, `input6G.f90` reads:

```fortran
local element connectivity + global element id
local node coordinates      + global node id
local boundary faces
```

The mappings are:

```fortran
ge(local_element) = global_element
gp(local_node)     = global_node
```

Each rank stores only its local elements, but interface nodes are duplicated on every rank whose elements use them.

### 4.3 Shared-node communication file

Each rank reads one `.com` file:

```fortran
numcommproc, ntotcomm

for each rank j:
    ncomms(j)
    comms(...)
```

The compressed arrays are:

```text
ncomms(j)   number of local nodes shared with rank j
scomms(j)   displacement of rank-j data inside comms
comms(k)    local node number to send to or receive from a peer rank
```

The old code uses the same `ncomms` and `scomms` arrays as both the send and receive descriptions in `MPI_ALLTOALLV`.

This requires pairwise communication lists to be symmetric:

- equal counts between each pair of ranks;
- corresponding shared nodes in the same pairwise order.

That hidden ordering requirement is fragile and must not be reproduced.

---

## 5. Shared-node numerical semantics

The old solver is based on:

```text
owned local elements
duplicated shared nodes
local element assembly
additive accumulation of shared-node contributions
```

For a nodal residual assembled from element contributions:

1. Each rank assembles contributions from its own elements.
2. For every shared local node, it sends its partial value to all neighbour ranks that also contain the node.
3. `MPI_ALLTOALLV` exchanges the partial values.
4. Each rank adds the received values to its local copy.

The pattern appears in:

- `getmassmatGG.f90`
- `pstiff.f90`
- `step1ggMPI.f90`
- `Step2Semi.f90`
- `step3ggMPI.f90`
- `multi4nPdiag.f90`
- `multi4nPdiag_vec.f90`
- the artificial-diffusion routines

For a scalar quantity:

```fortran
send(k) = local_value(comms(k))
MPI_ALLTOALLV(...)
local_value(comms(k)) += receive(k)
```

For velocity vectors, the same operation is performed with a contiguous vector MPI datatype.

This is an additive finite-element assembly operation, not merely a ghost-value copy.

---

## 6. Geometry and time-step reductions

Not every shared-node operation is additive.

### Characteristic length

`findlen123.f90` uses:

```text
minimum across all elements touching a shared node
maximum for aspect-ratio diagnostics
```

The shared-node operations are:

```fortran
alen(shared_node)       = min(local and received values)
alen_ratio(shared_node) = max(local and received values)
```

### Global time step

`alotim5.f90` computes a local minimum stable time step and then uses:

```fortran
MPI_ALLREDUCE(..., MPI_MIN, ...)
```

The global time step is therefore the minimum over all ranks.

These different reduction operators must be represented explicitly in the modern communication API:

- `ADD_VALUES` for assembled residuals and masses;
- `MPI_MIN` for characteristic lengths and global time step;
- `MPI_MAX` for diagnostics;
- `INSERT_VALUES` or broadcast for synchronised solution fields.

---

## 7. Legacy Pressure CG

### 7.1 Pressure operator storage

Each rank stores pressure coefficients only for its local elements:

```text
gstif   local tetrahedral off-diagonal coefficients
pdiag   assembled nodal diagonal
```

`pstiff.f90` assembles the local diagonal and then additively synchronises it across shared nodes. Every duplicated shared node therefore receives the complete global diagonal value.

### 7.2 Distributed matrix-vector multiplication

`multi4nPdiag.f90` performs:

1. local element off-diagonal multiplication;
2. additive shared-node accumulation with `MPI_ALLTOALLV`;
3. addition of the already-globalised nodal diagonal.

For a tetrahedron:

```text
y1 += A12 x2 + A13 x3 + A14 x4
y2 += A12 x1 + A23 x3 + A24 x4
y3 += A13 x1 + A23 x2 + A34 x4
y4 += A14 x1 + A24 x2 + A34 x3
```

This produces the complete `A*x` value on every local copy of a shared node.

### 7.3 Unique-node ownership for dot products

The pressure vector is duplicated on shared nodes, so CG dot products cannot sum every local copy.

The routines `converge_point_alloc*.f90` create:

```text
nresp   list of local nodes counted by this rank
nperrc  number of owned/countable local nodes
```

The basic method assigns each duplicated node to one rank, effectively favouring the lowest participating rank. Later methods rebalance shared-node ownership to reduce dot-product work imbalance.

CG computes local dot products over `nresp` only, followed by:

```fortran
MPI_ALLREDUCE(..., MPI_SUM, ...)
```

for:

- `(r,z)`
- `(p,A*p)`
- the new `(r,z)`

This ensures each global node is counted once.

### 7.4 Shared pressure-vector consistency

All copies of a shared pressure value remain equal because:

- every rank begins with the same duplicated value;
- every rank receives the same complete matrix-vector result for that node;
- every CG scalar coefficient is globally reduced;
- every local copy applies the same vector update.

The old solver therefore avoids a separate pressure-vector halo broadcast after every CG update.

### 7.5 Legacy convergence defect to avoid

In `conjugate_gradient.f90`, the relative-L2 convergence branch forms a local sum of squared residuals but combines it with:

```fortran
MPI_ALLREDUCE(..., MPI_MAX, ...)
```

rather than `MPI_SUM`.

That does not form the global Euclidean norm. It appears to be a legacy defect or an undocumented maximum-patch norm. PETSc KSP norms must be treated as the authoritative distributed pressure convergence measure.

---

## 8. Step 1 and Step 3

### Step 1

`step1ggMPI.f90`:

1. assembles momentum residuals from local elements;
2. additively synchronises the nodal vector residual;
3. applies the global inverse lumped mass;
4. updates all local copies of velocity.

### Step 3

`step3ggMPI.f90`:

1. assembles the pressure-gradient correction from local elements;
2. additively synchronises the nodal correction residual;
3. applies the global inverse lumped mass;
4. updates all local velocity copies.

The modern implementation must preserve the sequence:

```text
local element assembly
    -> global additive nodal assembly
    -> nodal scaling/update
    -> owner-to-ghost synchronisation
```

---

## 9. Convergence and ownership

`cnvrg_chck.f90` evaluates residuals only on nodes in `nresp` and then uses global reductions.

The ownership-balancing routines are not part of the physical domain decomposition. They exist only because the old solver stores duplicated vectors and needs a unique set of nodes for:

- CG inner products;
- residual norms;
- convergence diagnostics.

A PETSc distributed vector already assigns one global DOF owner. Therefore, the modern solver does not need:

```text
converge_poin_alloc
converge_poin_alloc2
converge_poin_alloc3
converge_poin_alloc4
nresp
nperrc
```

Their numerical purpose is replaced by PETSc ownership ranges.

---

## 10. Legacy output

The old solver does not gather a global solution onto rank zero.

Each rank opens and writes its own:

```text
VALTEMP
TPLOT
RESTART
debug
```

files with a rank suffix. `gp` is written so that post-processing can recover global node numbers.

For CBS3D++_SI, the scalable replacements should be:

- rank-local VTU pieces plus one rank-zero PVTU file; or
- parallel HDF5/XDMF through PETSc viewers.

Gathering the complete multi-million-element solution to rank zero should not be the default output path.

---

## 11. What must be preserved

The following legacy semantics are correct and important:

1. Spatial cell decomposition.
2. Fluid and solid cells may exist on the same rank.
3. One owner per global node/DOF.
4. Shared local representations for element kernels.
5. Additive accumulation of finite-element nodal contributions.
6. Global minimum for the stable time step.
7. Global sums for pressure-CG dot products.
8. Global sums for convergence norms.
9. Correct distribution of boundary faces and material identifiers.
10. Parallel output without gathering the entire solution.

---

## 12. What must not be copied

The modern implementation should not reproduce:

1. Pre-generated rank-specific `.plt` and `.com` files as the primary workflow.
2. Manual `comms/ncomms/scomms` arrays passed through numerical routines.
3. Pairwise-order-dependent `MPI_ALLTOALLV`.
4. Manual shared-node ownership balancing for dot products.
5. Custom distributed Pressure CG.
6. Legacy `mpif.h`.
7. Fake MPI stubs.
8. Repeated raw MPI communication embedded directly in element kernels.
9. One output file per rank without a collection descriptor.
10. Any legacy reduction that uses `MPI_MAX` where a global sum/norm is required.

---

## 13. Recommended modern architecture

### 13.1 Mesh and partitioning

Use PETSc `DMPlex` as the distributed unstructured-mesh layer:

```text
existing global tetrahedral mesh
    -> DMPlex creation
    -> attach material and BC labels
    -> PetscPartitioner
    -> ParMETIS or PT-Scotch
    -> DMPlexDistribute
```

Attach before distribution:

```text
material_id   on tetrahedral cells
bc_id         on boundary faces
global ids    for validation and output
```

Use a spatial partition across the complete CHT mesh. Never partition fluid and solid as separate MPI domains.

### 13.2 Staged integration with the current C++ solver

Do not immediately rewrite every numerical kernel around raw DMPlex traversal.

After distribution:

1. Extract each rank's owned cells and local/ghost nodes.
2. Populate compact rank-local CBS arrays.
3. Preserve the existing P1 tetrahedral element loops.
4. Use PETSc communication objects for local/global field assembly.

This keeps the validated CBS equations intact while modernising the distributed data layer.

### 13.3 Field layouts

Use PETSc sections with physically active DOFs.

Velocity:

```text
3 DOFs on fluid-connected nodes
interface velocity DOFs constrained to zero
no DOFs on solid-only nodes
```

Pressure:

```text
1 DOF on fluid-connected nodes
no DOF on solid-only nodes
prescribed-pressure DOFs constrained
```

Temperature:

```text
1 DOF on every thermal node
fluid, solid and interface included
```

Separate velocity, pressure and temperature sections/vectors are recommended for the first implementation because the CBS steps solve them independently.

### 13.4 Communication

Use:

```text
DMLocalToGlobal ADD_VALUES
```

for:

- lumped mass;
- pressure diagonal;
- Step 1 momentum residual;
- Step 2 pressure RHS;
- Step 3 correction residual;
- Step 4 thermal residual.

Use:

```text
DMGlobalToLocal INSERT_VALUES
```

after owner updates to refresh ghost values.

Use `PetscSF` reductions with the appropriate operator for:

- minimum characteristic length;
- maximum diagnostic quantities;
- specialised label or boundary metadata exchange.

### 13.5 Distributed pressure solve

Replace the current serial PETSc objects:

```text
PETSC_COMM_SELF
MatCreateSeqAIJ
VecCreateSeq
```

with distributed objects on the solver communicator:

```text
PETSC_COMM_WORLD
distributed Mat
distributed Vec
distributed KSP
```

Initial baseline:

```text
KSPCG
PCHYPRE + BoomerAMG
PCGAMG fallback
```

PETSc must own:

- global pressure numbering;
- off-process matrix insertion;
- matrix assembly;
- parallel matrix-vector products;
- halo communication;
- global dot products;
- pressure convergence norms;
- AMG setup and application.

This removes the legacy custom Pressure CG and is the primary route to significant acceleration.

---

## 14. CHT-aware load balancing

Equal element counts are not equal computational work.

Approximate costs:

```text
fluid cell:
    Step 1 momentum
    Step 2 pressure
    Step 3 correction
    Step 4 energy

solid cell:
    Step 4 energy only
```

The partition should therefore support cell weights. A starting model is:

```text
fluid tetrahedron    weight 4
solid tetrahedron    weight 1
interface-heavy cell optional additional weight
```

Weights must later be calibrated from measured per-cell timings.

The first partition validation must report per rank:

```text
owned fluid cells
owned solid cells
owned pressure DOFs
owned thermal DOFs
shared nodes
neighbour ranks
estimated work weight
```

---

## 15. Recommended implementation phases

### Phase A: distributed mesh only

- Create `PetscMeshDistributor`.
- Convert the existing global mesh to DMPlex.
- Attach material and BC labels.
- Partition and distribute.
- Extract rank-local cell/node maps.
- Verify global counts and label preservation.
- No numerical solver changes yet.

### Phase B: distributed field communication

- Create velocity, pressure and temperature sections.
- Create local and global vectors.
- Implement additive local-to-global residual assembly.
- Implement owner-to-ghost field updates.
- Validate against serial mass, pressure diagonal and residual arrays.

### Phase C: distributed pressure solve

- Assemble distributed pressure Mat/Vec.
- Use KSPCG + BoomerAMG.
- Validate pressure against the serial solver.
- Profile setup, MatMult, PCApply and global reductions.

### Phase D: distributed CBS steps

- Step 1 local assembly and global accumulation.
- Step 3 local assembly and global accumulation.
- Step 4 fluid-solid thermal assembly.
- Global time-step and convergence reductions.

### Phase E: scalable output and performance

- PVTU or parallel HDF5 output.
- Weighted partitioning.
- Hybrid MPI/OpenMP tuning.
- Strong and weak scaling on Sunbird.

---

## 16. Required validation gates

Do not proceed to the next phase unless the current phase passes.

### Mesh gate

- sum of owned cells equals global cell count;
- every global cell has exactly one owner;
- every global pressure DOF has exactly one owner;
- material and boundary labels are unchanged;
- interface faces remain conformal.

### Assembly gate

For a small case, compare serial and distributed global vectors:

```text
lumped momentum mass
thermal capacitance
pressure diagonal
Step 1 residual
Step 2 RHS
Step 3 residual
Step 4 residual
```

Required discrepancy should be near floating-point round-off, allowing for changed summation order.

### Pressure gate

Compare:

```text
pressure field
KSP residual
velocity correction
mass conservation
```

### Full-solver gate

Compare:

```text
residual history
velocity field
pressure field
temperature field
boundary fluxes
VTU output
```

### Scaling gate

Measure:

```text
1 rank x OpenMP threads
2 ranks
4 ranks
8 ranks
16 ranks
40 ranks
multiple nodes
```

The distributed pressure stage must show meaningful scaling because it is the dominant cost.

---

## 17. Immediate next development action

The archive provides enough information to stop requesting additional legacy MPI files.

The first new module should be:

```text
include/cbs/parallel/PetscMeshDistributor.hpp
src/parallel/PetscMeshDistributor.cpp
```

Its first release should only:

1. create a DMPlex from the existing global tetrahedral mesh;
2. attach cell material labels and boundary-face BC labels;
3. partition with PETSc;
4. distribute the mesh;
5. report owned and shared entities per rank;
6. verify global counts.

It should not yet alter CBS Step 1, Step 2, Step 3 or Step 4.
