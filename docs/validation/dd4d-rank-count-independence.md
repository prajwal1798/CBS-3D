# DD-4D rank-count independence validation

**Status:** PASS  
**Validation date:** 23 July 2026  
**Development branch:** `feature/domain-decomposition-sunbird`

## Purpose

DD-4D verifies that the distributed CBS3D++_SI solution is independent of MPI rank count and node placement within the accuracy of the configured inexact PETSc pressure solve.

The validation compares the same 2,383,581-tetrahedron conjugate heat-transfer case using:

- 2 MPI ranks on one compute node; and
- 40 MPI ranks distributed over two compute nodes, with 20 ranks per node and one CPU core per rank.

The comparison is performed directly on every unique owned global node rather than from screenshots or aggregate norms alone.

## Mesh and partition

| Quantity | Value |
|---|---:|
| Global tetrahedra | 2,383,581 |
| Global nodes | 407,458 |
| Global boundary faces | 28,328 |
| Two-rank pieces | 2 |
| Forty-rank pieces | 40 |

All 40 rank-local partitions passed the metadata and file-completeness audit. The partition format uses a one-based partition identifier and zero-based MPI rank:

```text
partition_id = mpi_rank + 1
```

The legacy 40-rank export retained original Gmsh tetrahedron tags. The solver collectively verified that these formed a unique contiguous interval and normalised the IDs by the detected offset of 28,328 to the dense range `1..2,383,581`.

## Distributed execution

The 40-rank run used two compute nodes:

```text
nodes               : 2
MPI ranks           : 40
MPI ranks per node  : 20
CPU cores per rank  : 1
OpenMP threads/rank : 1
```

The production loop completed normally:

```text
MPI ranks                  : 40
completed iterations       : 2
pressure matrix builds     : 1
AMG hierarchy builds       : 1
return code                : 0
```

The distributed pressure matrix and AMG hierarchy were assembled once and reused across both iterations.

## Field-level comparison

The final state at iteration 2 was reconstructed from owner copies only and matched by `global_node_id` across all 407,458 nodes.

| Field | Relative L2 difference | Maximum absolute difference | Allowed maximum | Result |
|---|---:|---:|---:|:---:|
| `u` | 6.36365719e-7 | 1.01465584e-6 | 2.77042598e-6 | PASS |
| `v` | 8.91123400e-5 | 5.77165325e-7 | 2.77042598e-6 | PASS |
| `w` | 1.41131699e-6 | 8.14193697e-7 | 2.77042598e-6 | PASS |
| Velocity vector | 7.77776064e-7 | 1.06406501e-6 | 2.77042598e-6 | PASS |
| Pressure | 4.48044747e-6 | 1.62101891e-3 | 2.59568754e-3 | PASS |
| Temperature | 0.0 | 0.0 | 3.10013345e-6 | PASS |

The velocity-vector norm is the primary velocity acceptance metric. Component-wise relative errors remain diagnostic because a physically small component can exhibit a large relative value while its absolute discrepancy remains negligible.

The maximum physical residual-history scalar difference was:

```text
4.95523901e-4  (approximately 0.0496%)
```

## PETSc convergence

The Krylov histories were:

```text
2 ranks  : CG iterations [11, 9]
40 ranks : CG iterations [12, 10]
```

True relative residuals were:

```text
2 ranks  : [1.078764663438848e-9, 4.492329639702683e-6]
40 ranks : [1.390811652231012e-9, 9.220242791689309e-6]
```

The small difference in iteration count is expected because BoomerAMG constructs a decomposition-dependent hierarchy. Both rank counts satisfy the configured pressure-solver tolerance and produce equivalent physical fields.

## Numerical invariants

The principal printed solution diagnostics agree across rank counts:

| Quantity | Iteration 1 | Iteration 2 |
|---|---:|---:|
| `DivRMS` | 1.6228e-9 | 1.0225e-9 |
| `Umax` | 0.0706133 | 0.0770426 |
| `dUmax` | 0.072782 | 0.0101353 |
| `Tmax` | 300.007 | 300.013 |

## Performance observation

For the two-iteration validation case:

```text
2-rank mean iteration time  : approximately 1.52031 s
40-rank mean iteration time : approximately 0.110237 s
speedup                     : approximately 13.79x
strong-scaling efficiency   : approximately 68.9% relative to 2 ranks
```

The 40-rank maximum setup time was approximately 7.19 s. For a two-iteration test, setup dominates; for long production runs, the persistent setup cost is amortised.

## Accepted distributed capabilities

DD-4D establishes the following capabilities for the current development branch:

- unique owned-element decomposition;
- owned/ghost nodal storage;
- forward owner-to-ghost field exchange;
- reverse-add shared-node residual assembly;
- owner-only nodal updates and strong boundary conditions;
- distributed CBS Steps 1–4;
- distributed PETSc pressure assembly and solve;
- persistent matrix, KSP and AMG objects;
- distributed convergence monitoring and residual CSV output;
- rank-local VTU with PVTU/PVD aggregation;
- cross-node MPI execution; and
- rank-count-independent primary solution fields.

## Conclusion

```text
CBS3D++_SI distributed-memory numerical implementation : PASS
DD-4D rank-count independence                          : PASS
Core domain decomposition                              : COMPLETE
```

## Remaining production-hardening work

The next stage is DD-5:

1. long-run stability and steady-state stopping;
2. appended binary/compressed distributed output;
3. rank-filtered parallel logging;
4. formal strong-scaling instrumentation; and
5. distributed restart/checkpoint validation.
