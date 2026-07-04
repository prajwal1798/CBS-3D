# CBS3D++_SI — Complete Source, Mathematics and Parallelisation Audit

**Audit scope:** active July 2026 `src/` and `include/` trees, `CMakeLists.txt`, `Makefile`, the supplied `blanket` input set, and the four supplied CBS references.  
**Explicit exclusion:** the older duplicate `src/src/` tree.  
**Code changes made:** none.  
**Build check:** the active serial/OpenMP source tree was configured and compiled cleanly as C++20 with warnings enabled.

---

## 1. Audit status legend

| Mark | Meaning |
|---|---|
| **CC** | Confirmed correct for the stated formulation and active path. |
| **S** | Suspicious, incomplete, fragile, or a deliberate simplification that requires documentation/validation. |
| **I** | Confirmed incorrect, internally inconsistent, or an advertised path that is not actually implemented. |
| **D** | Dormant/inactive code or configuration; not currently affecting the blanket run. |
| **MPI** | A required change or communication dependency for distributed execution. |

A finding marked **S** is not automatically a mathematical error. Several are legitimate reduced forms of the general CBS algorithm, but they must be stated explicitly because they differ from the full equations in the references.

---

## 2. Sources examined

### 2.1 Active code and build files

- `CMakeLists.txt`
- `Makefile`
- all headers under `include/cbs/`
- all active implementation files under `src/`
- the `blanket.plt`, `.bco`, `.par`, `.material`, and `.matprop` files

### 2.2 CBS references

1. O. C. Zienkiewicz and co-authors, **A general algorithm for compressible and incompressible flows — the characteristic-based split algorithm**, especially Split A, Steps 1–4 and boundary conditions.
2. P. Nithiarasu and C.-B. Liu, **An artificial-compressibility-based CBS scheme for steady and unsteady turbulent incompressible flows**.
3. P. Nithiarasu, **An efficient artificial-compressibility scheme based on the CBS method for incompressible flows**.
4. P. Nithiarasu, R. Codina and O. C. Zienkiewicz, **The Characteristic-Based Split scheme — a unified approach to fluid dynamics**.

The active solver is a **semi-implicit incompressible pressure-Poisson CBS solver**, not the fully explicit CBS-AC algorithm. The two AC papers are therefore used mainly to assess time-step, residual and dual-time concepts; their artificial-compressibility equations must not be inserted into the present pressure-Poisson formulation without a deliberate method change.

---

## 3. Executive conclusion

### 3.1 What is already sound

- **CC:** The active solver has a coherent four-stage CBS structure: pressure-free momentum predictor, pressure solve, velocity correction, and scalar/temperature transport.
- **CC:** P1 tetrahedral shape gradients, determinant/volume relations, lumped mass, pressure stiffness, boundary-face matching, outward normals and the conformal CHT interface treatment are internally consistent.
- **CC:** The dimensional CHT diffusion, heat capacity, volumetric source and prescribed heat-flux terms are dimensionally correct on the active path.
- **CC:** Element colouring correctly prevents shared-node OpenMP scatter races in the coloured kernels.
- **CC:** The native pressure method is a conventional matrix-free preconditioned conjugate-gradient solver with a Jacobi option.
- **CC:** The current blanket input is topologically consistent: 2,383,581 unique tetrahedra, 407,458 unique nodes, 28,328 valid physical boundary triangles, matching material connectivity, exact material totals and exact BC totals.

### 3.2 What blocks a trustworthy fully parallel release

1. **I — The executable is not distributed.** MPI currently launches ranks, but rank 0 alone reads, advances and writes the complete mesh.
2. **I — PETSc is sequential.** It uses `PETSC_COMM_SELF`, sequential matrix/vector types, and rank-0-only invocation.
3. **I — The advertised transient mode is not implemented.** Physical-time history arrays and real-time hooks exist, but no valid BDF/dual-time update is assembled.
4. **I — PETSc can reuse a stale pressure matrix when `dt` changes.** The cached-system test does not include the time-scaled matrix coefficients.
5. **I — Optional Step-2 pressure time-step correction is dimensionally invalid in dimensional mode and changes `dt` after Step 1, producing a mixed-time-step iteration.
6. **I — BC 504 is implemented in boundary/preprocessing code but rejected by the input reader.
7. **I — `bad_detj` in the OpenMP Step-1 assembly has a data race.
8. **S — The momentum viscous term is a componentwise Laplacian, not the complete deviatoric-stress operator written in the references.
9. **S — The active Step 3 omits the characteristic higher-order pressure contribution shown in the general CBS formulation.
10. **S — Several input controls are accepted but ignored or rejected later, including explicit CBS, banded pressure, restart, `theta1`, buoyancy controls and several stabilisation controls.

### 3.3 Overall mathematical judgement

The **active fixed-time-step, steady, laminar, incompressible, constant-property CHT path is coherent and has already been benchmarked**, but it should be described as a **simplified semi-implicit projection/CBS specialization**, not as a literal implementation of every term in the most general Split-A equations. No evidence was found of a fundamental sign error in the active pressure projection or CHT diffusion/flux terms. The principal risks are conditional/dormant paths, transient claims, PETSc caching, boundary-contract inconsistencies, and the absence of distributed ownership/communication.

---

## 4. Complete file and call graph

```text
main
 ├─ optional MPI_Init_thread(MPI_THREAD_FUNNELED)
 ├─ broadcastCaseName
 ├─ Solver(case)
 ├─ Solver::setMpiContext
 └─ rank 0 only: Solver::run
      ├─ Solver::initialise
      │   ├─ MeshIO::readAll
      │   │   ├─ CaseFiles construction
      │   │   ├─ MeshIO::readSizes
      │   │   ├─ CBSStateSI::set_problem_sizes
      │   │   ├─ MeshIO::readMeshFile
      │   │   ├─ MeshIO::readBoundaryFile
      │   │   ├─ MeshIO::readParameterFile
      │   │   ├─ MeshIO::readMaterialFiles / defaults
      │   │   └─ MeshIO::initialiseFields
      │   ├─ Preprocess::validateBoundaryFlags
      │   ├─ Preprocess::shapeFunctionDerivatives
      │   ├─ Preprocess::assignBoundaryFaceNumbers
      │   ├─ Preprocess::getNormals
      │   ├─ Preprocess::massMatrix
      │   ├─ Preprocess::classifyFaceEdges
      │   ├─ Preprocess::elementSize
      │   ├─ Preprocess::wallDetermination
      │   ├─ Preprocess::computeMassFlowInletVelocity
      │   ├─ Preprocess::initialiseVelocityMagnitude
      │   ├─ Preprocess::detectPressureBoundaryNodes
      │   ├─ Coloring::build
      │   └─ Boundary::{applyTemperature,applyVelocity,applyPressure}
      ├─ Solver::preparePressureSystem
      │   ├─ PressureAssembly::buildElementPressureTerms
      │   └─ PressureAssembly::buildGlobalPressureTerms
      ├─ initial TimeStep::{computeTimeStep,updateLhsDiagonal}
      └─ iteration loop: Solver::advanceOneStep
          ├─ copy current fields to iteration-old fields
          ├─ TimeStep::computeTimeStep
          ├─ TimeStep::updateLhsDiagonal
          ├─ Steps::step1
          │   ├─ MomentumAssembly::assembleStep1Rhs
          │   ├─ nodal predictor update
          │   └─ Boundary velocity/symmetry/backflow package
          ├─ optional TimeStep::applyStep2PressureTimeStepCorrection
          ├─ optional second TimeStep::updateLhsDiagonal
          ├─ Steps::step2
          │   ├─ PressureAssembly::assembleStep2Rhs
          │   ├─ native ConjugateGradient::solvePressure
          │   │   └─ MatrixVectorCalc::pressureMultiply
          │   └─ or PetscPressureSolver::solvePressure
          ├─ Steps::step3
          │   ├─ element pressure-gradient correction assembly
          │   ├─ nodal velocity correction
          │   └─ Boundary velocity/symmetry/backflow package
          ├─ Steps::step4
          │   ├─ EnergyAssembly::assembleStep4Rhs
          │   ├─ nodal temperature update
          │   └─ Boundary::applyTemperature
          ├─ Solver::updateVelocityMagnitude
          ├─ Convergence::evaluate
          └─ Post::{writeResidualRow,printProgressLine,writeSolution}
```

### Distributed implication

The eventual MPI solver must preserve this numerical order. Communication is required **between**, not inside, the mathematical steps:

```text
forward halo old field
→ owned-element assembly
→ reverse additive halo of nodal contributions
→ owner nodal update / strong BC
→ forward halo updated field
```

---

## 5. State variables, storage and indexing audit

### 5.1 Fundamental types and arrays

| Item | Audit |
|---|---|
| `Int = std::int32_t` | **CC** for the current mesh. **S/MPI:** global IDs above 2,147,483,647 would overflow; PETSc may be configured with 64-bit indices. Decide one global-index type before general release. |
| `Real = double` | **CC** for CFD/FEM operations. |
| `Array1D/2D/3D` | **CC:** deliberately 1-based and contiguous; index zero is unused. Column-major-style offsets reproduce the legacy scientific array ordering. |
| Debug bounds checks | **CC**, but disabled in release by design. Input validation must therefore remain strict. |

### 5.2 Principal state groups

| State group | Meaning and active use | Audit |
|---|---|---|
| `intma`, `coord`, `iside`, `flag_list` | Element connectivity, coordinates, boundary records and BCO mapping | **CC** |
| `unkno`, `pres`, `temperature` | Current velocity, pressure and temperature | **CC** |
| `unkn1`, `pres1`, `temperature1` | Field at the beginning of the current CBS iteration | **CC** |
| `unknn1`, `unknn2`, `tempert1`, `tempert2` | Intended physical-time histories | **D/I:** allocated and initialized but not advanced into a valid transient discretisation |
| `rhs`, `rhs1` | Momentum/pressure-correction vector storage and thermal residual | **CC**, but the same arrays are reused step by step and therefore require strict ordering |
| `dNkdx`, `detJ` | Constant P1 tetrahedral shape gradients and determinant | **CC** |
| `M_diag`, `elcoe_e` | Lumped nodal and element mass quantities | **CC** |
| `elcoe2` | Inverse momentum mass/time diagonal | **CC** on fluid nodes |
| `elcoe2p` | Inverse thermal-capacitance/time diagonal | **CC** on full CHT domain |
| `pdiagE`, `gstifE` | Unscaled element pressure stiffness | **CC** |
| `pdiag`, `gstif` | Time-scaled global/element compact pressure operator | **CC** for native serial/OpenMP path |
| material arrays | Element material ID and constant element properties | **CC** for one fluid ID `0` plus one-or-more solid IDs; **S:** phase name is ignored |
| colouring arrays | Full and pressure-only colour groupings | **CC**; pressure-only grouping is not used by the matrix-vector kernel |
| MPI rank/size fields | Context only | **D:** no local ownership, ghost or communication state in the active solver |
| banded arrays | Legacy/inactive pressure path | **I/D:** dimensions do not represent a real bandwidth and path is rejected at runtime |

### 5.3 Time-level convention

The active steady/pseudo-time convention is:

```text
begin iteration:
    unkn1       = unkno
    pres1       = pres
    temperature1= temperature

Step 1: unkno becomes u*
Step 2: pres becomes p^(n+1)
Step 3: unkno becomes u^(n+1)
Step 4: temperature becomes T^(n+1)
```

**CC:** This is internally consistent for the active pseudo-time iteration.  
**I for true transient:** the physical-history arrays do not participate in a BDF or dual-time scheme.

### 5.4 Material-domain convention

```text
mat_elem(e) == 0 : fluid
mat_elem(e) != 0 : solid
```

- **CC** for the supplied blanket case.
- **S:** the `.matprop` phase/name field does not control phase classification. A second fluid material with ID other than zero would be treated as solid.
- **MPI:** every rank needs material IDs for owned elements and material-domain node masks after reverse/forward ownership reconciliation.

---

## 6. File-by-file and function-by-function audit

## 6.1 Build system

### `CMakeLists.txt`

**Role:** selects C++20 sources and optional OpenMP, MPI and PETSc capabilities.

| Block | Inputs/outputs/state | Audit |
|---|---|---|
| project/options | User CMake options; generates one executable | **CC** for serial/OpenMP. |
| source list | Active `.cpp` files | **CC:** excludes duplicate `src/src`. **MPI:** partition and halo modules are absent. |
| warning/optimisation flags | Compiler flags | **CC**; clean build completed. |
| OpenMP | `find_package(OpenMP)` and target link | **CC**. |
| MPI | `find_package(MPI CXX)` | **CC** for bootstrap only; no distributed source path. |
| PETSc | manually finds include path and one `libpetsc` | **S/I:** fragile because PETSc can require transitive libraries and ABI-consistent MPI/compiler flags. It also permits PETSc without MPI even though distributed PETSc is the target. |

**Required MPI change:** use PETSc's exported CMake/package configuration or its generated variables/pkg-config, require the same MPI compiler wrapper, and add partition/halo sources.

### `Makefile`

**Role:** alternative build with isolated variant object directories.

- **CC:** separate build directories prevent stale object reuse across serial/MPI/PETSc/OpenMP variants.
- **CC:** dependency files and compiler wrappers are conventional.
- **S:** PETSc is linked as only `-lpetsc`; this may omit PETSc's required transitive libraries.
- **MPI:** `run-mpi` currently launches a bootstrap executable, not a distributed numerical solver.

---

## 6.2 Program entry and driver

### `src/main.cpp`

#### `broadcastCaseName(std::string&, int)`

- **Role:** broadcasts rank-0 case name.
- **Input:** rank-0 string, communicator rank.
- **Output/modified state:** resizes and fills the string on all non-root ranks.
- **Mathematics/CBS:** none.
- **OpenMP:** none; called before threaded kernels.
- **MPI:** collective and correct for the bootstrap.
- **Finding:** **CC**.

#### `main(int, char**)`

- **Role:** initializes optional MPI, obtains the case name, constructs the solver and controls finalisation.
- **Input:** command line or interactive case name.
- **Output:** process exit code and solver files from rank 0.
- **Modified state:** global MPI state; solver instance.
- **OpenMP:** requests `MPI_THREAD_FUNNELED`, consistent with MPI calls made by the main thread only.
- **PETSc/CG:** PETSc shutdown is rank-0-only; initialization is deferred into the PETSc pressure solver.
- **Finding:**
  - **CC:** MPI thread-level check and errors-return handler.
  - **I/MPI:** only rank 0 calls `Solver::run`; all other ranks wait.
  - **I/PETSc:** distributed PETSc cannot be initialized and used only on rank 0 while MPI world contains all ranks. PETSc initialization/finalization must be a program-level all-participating-ranks responsibility.

### `include/cbs/solver/Solver.hpp` and `src/solver/Solver.cpp`

#### constructor / `state()` / `setMpiContext()`

- **Role:** stores case name, exposes state, records rank/size.
- **Finding:** **CC**, but rank/size are informational only.

#### `Solver::run()`

- **Role:** complete steady/pseudo-time loop.
- **Inputs:** initialized case through `case_name_`.
- **Outputs:** final state and output files.
- **Modified state:** all numerical fields, residuals, timing and iteration counters.
- **CBS:** enforces Step 1 → Step 2 → Step 3 → Step 4 order.
- **Finding:** **CC** for the active steady mode; **I** if `transient_on=1` is interpreted as a physical transient calculation.
- **MPI:** iteration decisions, stop reasons and output triggers must be globally identical; convergence uses reductions and only rank 0 writes metadata.

#### `Solver::initialise()`

- **Role:** input, geometry, mass, BC classification, colouring and initial strong BCs.
- **Finding:** call order is **CC**. Pressure fixed nodes are detected before the pressure operator is used. Element colouring occurs after material IDs are known.
- **Performance:** material/interface masks are reconstructed in later modules instead of being cached once.
- **MPI:** replace full-file read with rank-local partition read and reconcile shared nodal quantities.

#### `Solver::preparePressureSystem()`

- **Role:** computes fixed geometry pressure stiffness before iterations.
- **Finding:** **CC** when mesh/material topology is static.
- **MPI:** build local element contributions and distributed row adjacency; geometry terms remain reusable, but time scaling may change.

#### `Solver::advanceOneStep(Int)`

- **Role:** copies iteration-old fields, computes timestep/diagonals, executes all four steps and convergence.
- **Finding:**
  - **CC:** normal fixed-`dt` blanket path.
  - **I conditional:** when Step-2 pressure timestep correction is enabled, Step 1 uses the original `dt` but Steps 2–4 use a newly altered `dt`; one CBS iteration then has inconsistent temporal coefficients.
  - **I transient:** no physical-history advancement or inner pseudo-time convergence structure.
- **MPI:** forward/reverse halo boundaries must be inserted at precise points; see Section 13.

#### `Solver::updateVelocityMagnitude()`

- **Role:** derives nodal speed for diagnostics/timestep.
- **Finding:** **CC**, apart from initialization floor discussed under preprocessing.
- **MPI:** owner nodes update; forward halo before any element calculation that reads speed.

#### `steadyStateReached()` / `transientEndTimeReached()`

- **Role:** stopping controls.
- **Finding:** steady test delegates correctly. Physical end-time test alone does not make the scheme transient; **I** as part of the currently advertised transient path.
- **MPI:** root-independent global Boolean decision required.

---

## 6.3 Core headers

### `Types.hpp`

- **Role:** numerical aliases.
- **Finding:** **CC**, with 32-bit global-index limitation noted above.

### `CaseFiles.hpp`

- **Role:** deterministic construction of case filenames.
- **Finding:** **CC** for monolithic cases.
- **MPI:** add rank-partition naming without changing numerical modules.

### `Array.hpp`

#### `Array1D`, `Array2D`, `Array3D`

- **Role:** simple 1-based scientific arrays.
- **Inputs:** positive dimensions and indices.
- **Outputs/state:** owned contiguous `std::vector` storage.
- **Finding:** **CC**. The plain design matches the requested student-readable style.
- **MPI:** arrays should be resized to local owned+ghost counts rather than global counts; do not replace them with opaque frameworks unnecessarily.

### `RunConfig.hpp`

- **Role:** all mesh constants, solver controls, physical parameters and BC IDs.
- **Finding by group:**

| Control | Finding |
|---|---|
| P1 tetra dimensions/constants | **CC** |
| `theta[1]`/theta2 forced to 1 | **S:** only fully implicit pressure weighting is supported. |
| `theta[0]`/theta1 | **D/I:** parsed and checked but not used in active assembly. |
| explicit `cbs_scheme=0` | **I contract:** accepted, rejected by `Steps`. |
| `solver_opt=2` banded | **I contract:** accepted, rejected by Step 2. |
| restart option | **I contract:** accepted; no restart restoration path. |
| `alpha_sf`, `k_ratio`, `source_solid`, `art_diff`, `rem_deltp`, parts of `beta1`/`p_ref` | **D:** no active numerical effect or incomplete effect. |
| Rayleigh/Richardson/body-force controls | **D/I:** no buoyancy force assembled in Step 1. |
| material/dimensional/mass-flow controls | **CC** for supplied blanket case. |

### `CBSStateSI.hpp`

#### `initialise_local_topology()`

- **Role:** defines tetrahedron local face-node mapping and six compact off-diagonal pairs.
- **Finding:** **CC** and consistent with pressure compact storage.

#### `set_problem_sizes(...)`

- **Role:** validates counts and allocates every field.
- **Finding:**
  - **CC:** active arrays correctly dimensioned for monolithic P1 tetra solver.
  - **S:** very high peak memory because all global fields, histories and inactive solver workspaces are allocated on one rank.
  - **I/D:** banded work arrays have placeholder dimensions and cannot support the advertised banded solve.
  - **MPI:** replace global sizes by `n_owned + n_ghost` for nodal fields and `n_owned_elements` for assembly; preserve local/global maps separately.

---

## 6.4 Input and mesh I/O

### `MeshIO::readAll()`

- **Role:** complete input orchestration.
- **Finding:** **CC** ordering: sizes → allocation → mesh/BC/parameters/materials → field initialization.

### `MeshIO::readSizes()`

- **Role:** reads solver option and mesh/BC headers before allocation.
- **Finding:**
  - **CC:** basic positivity/range checks.
  - **S:** accepts modes that later fail (`solver_opt=2`, explicit CBS).
  - **MPI:** monolithic global header is unsuitable for independent rank-local input unless metadata supplies global and local counts.

### `MeshIO::readMeshFile()`

- **Role:** loads tetra connectivity, coordinates and boundary records.
- **Geometry:** no mathematical transformation; assumes P1 tetrahedra.
- **Finding:**
  - **CC:** node/connectivity/range validation and boundary-parent validation.
  - **S:** does not explicitly detect duplicate or missing element/node IDs; duplicate entries can overwrite allocated rows. The supplied blanket input independently passed uniqueness/completeness checks.
  - **MPI:** each rank reads only its owned cells plus required ghost-node coordinates and physical boundary faces owned by that partition.

### `MeshIO::readBoundaryFile()`

- **Role:** maps raw side IDs to solver BC IDs.
- **Finding:**
  - **CC:** all supplied blanket IDs 511/520/530/532 map correctly.
  - **I:** solver BC 504 is documented and implemented elsewhere but omitted from the reader's accepted solver-BC set, so it cannot be used through normal input.
  - **S:** validation logic is duplicated across MeshIO and Preprocess, increasing divergence risk.

### `MeshIO::readParameterFile()`

- **Role:** parses `.par` controls.
- **Finding:**
  - **CC:** supplied blanket file parses into the expected fixed-dt, steady, dimensional CHT setup.
  - **S:** optional trailing groups are positional. A partially old file can shift later mandatory meanings rather than fail at a named key.
  - **S/I:** accepts controls that are inactive or rejected later.
  - **S:** insufficient early validation for all denominators, such as strictly positive Reynolds/Prandtl/property combinations in every branch.
  - **S:** `nusselt_diameter` may be required even when Nusselt calculation is not implemented/enabled.

### material readers/defaults

- **Role:** read element material mapping and constant properties or supply a default fluid.
- **Dimensional consistency:** `rho`, `cp`, `k`, `mu`, `Q` units are correct for the dimensional Step 1/Step 4 paths.
- **Finding:**
  - **CC:** current blanket material connectivity and counts match the mesh exactly.
  - **S:** phase/name is read but classification is hard-coded by element material ID zero versus nonzero.
  - **S:** no explicit duplicate/missing element-ID set audit inside the C++ reader, although line count/connectivity checks catch many failures.
  - **MPI:** map material records by global element ID during partition export; avoid each rank scanning the global file.

### field initialization

- **Role:** initializes velocity, pressure, temperature, histories and transport coefficients.
- **Finding:**
  - **CC:** temperature and velocity initial values for current steady blanket run.
  - **S:** fluid starts at `p0=101325` while outlet is immediately set to gauge pressure zero; `p_ref` is not used. Output semantics alternate between absolute-looking initial pressure and gauge pressure.
  - **I transient:** history arrays are zeroed rather than initialized/rotated into a valid physical-time method.
  - **S:** `ani` remains a legacy global diffusion scale and does not represent all material diffusivities in dimensional CHT.

---

## 6.5 Input-set audit: `blanket`

| Quantity | Verified value | Audit |
|---|---:|---|
| tetrahedra | 2,383,581 | **CC** unique and valid |
| nodes | 407,458 | **CC** unique and all referenced IDs valid |
| physical boundary triangles | 28,328 | **CC** unique; parent tetra and face-node subsets valid |
| fluid elements | 757,851 | **CC** |
| solid elements | 1,625,730 | **CC** |
| BC 511 inlet | 757 | **CC** |
| BC 520 outlet | 753 | **CC** |
| BC 530 adiabatic wall | 21,938 | **CC** |
| BC 532 heat-flux wall | 4,880 | **CC** |

Current properties:

```text
Water:       rho=997 kg/m^3, cp=4182 J/(kg K), k=0.606 W/(m K), mu=8.9e-4 Pa s
EUROFER-like solid: rho=7850 kg/m^3, cp=500 J/(kg K), k=28.3 W/(m K)
Applied heat flux: 5.0e5 W/m^2
Mass flow: 0.0157 kg/s
Initial/inlet temperature: 300 K
Fixed dt: 1.0e-5 s
```

**Exposure of dormant defects:** the current case has `transient_on=0`, fixed `dt`, `step2_check=0`, native CG and no BC 504. Therefore the transient, PETSc-stale-matrix, Step-2-dt and BC504 defects are not exercised by the present two-iteration blanket run.

---

## 6.6 Preprocessing and geometry

### `Preprocess::validateBoundaryFlags()`

- **Role:** validates mapped solver BC IDs.
- **Finding:** **S:** useful defensive check, but duplicates and disagrees with MeshIO on BC 504. Consolidate one authoritative table.

### `Preprocess::shapeFunctionDerivatives()`

- **Role:** computes affine P1 tetrahedron Jacobian determinant and constant global shape gradients.
- **Equation:** for each tetrahedron, `grad N_a = J^{-T} grad_hat N_a`, `V = det(J)/6`.
- **Finding:** **CC**. Positive orientation is enforced by rejecting nonpositive/nonfinite determinants.
- **OpenMP:** safe because each iteration writes a distinct element.
- **MPI:** owned-element local calculation only; no communication.

### `Preprocess::assignBoundaryFaceNumbers()`

- **Role:** matches each triangular boundary record to one of the four parent-tetra local faces using sorted node triplets.
- **Finding:** **CC** and robust to face-node order.
- **MPI:** partition exporter should preserve parent owned element and face; no cross-rank face ownership duplication.

### `Preprocess::getNormals()`

- **Role:** computes outward area-weighted face normals and element-face normal storage.
- **Equation:** cross product gives `2A n`; sign is corrected against opposite tetra node.
- **Finding:** **CC**.
- **MPI:** local physical faces only; shared internal partition faces are not physical BC faces.

### `Preprocess::massMatrix()`

- **Role:** assembles element/nodal lumped mass and thermal capacitance.
- **Equations:**
  - momentum lump `V/4` per P1 tetra node;
  - thermal lump `rho cp V/4`;
  - consistent P1 mass coefficients retained in element form.
- **Finding:** **CC**.
- **MPI:** shared nodal mass/capacity contributions require reverse-add to owners, followed by forward distribution.

### `Preprocess::classifyFaceEdges()`

- **Role:** marks interior, prescribed-velocity and other exterior faces for convection/diffusion boundary treatment.
- **Finding:** **CC** for active BC taxonomy; dependent on the duplicated BC table.

### `Preprocess::elementSize()`

- **Role:** computes tetrahedral length based on minimum altitude.
- **Equation:** `h_a = 3V/A_opposite`, then element minimum.
- **Finding:** **CC**, matching the characteristic-length definition used in the supplied CBS-AC references.
- **MPI:** local elements; nodal minima need owner reduction if nodal local stepping is retained.

### `Preprocess::wallDetermination()`

- **Role:** marks external no-slip walls and fluid-solid interface nodes; accumulates wall normals.
- **Finding:** **CC** for conformal CHT: interface velocity is zero and interface temperature remains shared/unconstrained.
- **S:** material node masks are recomputed in multiple later modules instead of stored once.
- **MPI:** node classification must be reconciled across all incident partitions; a node is interface/no-slip if any global incident cell/face establishes it.

### `Preprocess::computeMassFlowInletVelocity()`

- **Role:** computes inlet area and uniform speed from mass flow.
- **Equation:** `U = mdot/(rho A)`.
- **Finding:** **CC** dimensionally.
- **MPI:** inlet area is a global sum over owned inlet faces (`MPI_Allreduce`).

### `Preprocess::initialiseVelocityMagnitude()`

- **Role:** computes nodal speed.
- **Finding:** **S/minor I:** adds/floors with approximately `1e-16`, so an exactly stationary node can report a nonzero speed of about `1e-8`. Use a zero value and protect only divisions.

### `Preprocess::detectPressureBoundaryNodes()`

- **Role:** collects prescribed-pressure nodes, using outlet nodes or a fallback reference node.
- **Finding:** **CC** for a pressure-Poisson gauge constraint.
- **S:** reference-node choice and `p_ref` semantics should be explicit; multiple ranks must agree on one global fallback if no outlet exists.
- **MPI:** each rank identifies local candidates; globally select/own constrained global IDs.

### `Coloring::build()`

- **Role:** greedy node-conflict element colouring; also prepares a fluid pressure-only subset.
- **Finding:**
  - **CC:** no two elements of one colour share a node; coloured scatter is race-free.
  - **D/performance:** pressure-only colour arrays are built but the active pressure matrix-vector multiply loops the full colours and skips solid elements.
  - **S/minor:** an unused local `used` object remains.
- **MPI:** colour only owned local elements. Partition sharing is handled by reverse-add communication, not by global colouring.

---

## 6.7 Boundary conditions

### `Boundary::applyTemperature()`

- **Role:** applies strong nodal temperature values and leaves natural flux/interface conditions for assembly.
- **Inputs:** boundary-face IDs, configured inlet/wall temperatures, material topology.
- **Modified state:** `temperature` at constrained nodes.
- **CHT:**
  - BC 511 applies inlet temperature;
  - BC 532 is not strongly constrained because its heat flux is assembled in Step 4;
  - BC 901 is a conformal shared-node interface and correctly receives no artificial interface source.
- **Finding:** **CC** for current blanket CHT.
- **S:** at nodes shared by conflicting strong thermal BC faces, last-write/priority behavior depends on traversal order. Such corner conflicts should be detected once during preprocessing.
- **MPI:** owner must resolve each constrained global node once and forward the value to ghosts.

### `Boundary::applyVelocity()`

- **Role:** applies inlet, benchmark, no-slip and interface velocity conditions.
- **Modified state:** `unkno`.
- **Finding:**
  - **CC:** solid-touch/interface velocity is finally forced to zero.
  - **S:** BC 511 computes direction from each face normal and writes shared nodal values face by face. On a curved/non-planar inlet, a shared node can be overwritten by the last face. The supplied inlet is effectively planar, so current exposure is low.
  - **S:** benchmark BC 507/508 formulas are geometry-specific legacy behavior, not general inlet profiles.
  - **Performance:** reconstructs material-touch masks repeatedly.
- **MPI:** reduce/resolve boundary markers and inlet direction at owner nodes; do not allow ghost-face traversal order to determine values.

### `Boundary::applyPressure()`

- **Role:** applies fixed pressure/gauge nodes.
- **Finding:** **CC** for current outlet gauge pressure.
- **S:** `p_ref` is not used to define absolute output pressure.
- **MPI:** only owner rows impose the constraint; ghost values are synchronized afterward.

### `Boundary::applySymmetry()`

- **Role:** removes the normal component of velocity at slip/symmetry faces.
- **Equation:** `u <- u - (u·n)n`.
- **Finding:** **CC** on a single planar symmetry surface.
- **S:** a node lying on two non-coplanar symmetry faces is projected sequentially; the final velocity can depend on face order and need not satisfy both normals. Assemble an orthonormal constraint basis per node instead.
- **MPI:** owner constructs the complete set of incident symmetry normals across partitions.

### `Boundary::applyOutletBackflowControl()`

- **Role:** clips velocity entering through a nominal outlet.
- **Finding:** **S:** useful stabilization but not a neutral Navier–Stokes BC. It can modify mass conservation and should be optional, reported and measured.
- **MPI:** owner-face/node decisions must be deterministic; global outlet flow should be monitored.

### Boundary application order

The active package applies strong velocity, symmetry and outlet corrections after Steps 1 and 3. This agrees with the practical CBS recommendation to impose correct velocity values on solid boundaries in predictor/correction stages. However, duplicated passes and final solid-zero enforcement make the actual priority implicit. **S:** preprocess a node-level BC priority map and reject incompatible combinations.

---

## 6.8 Time stepping and diagonals

### `TimeStep::computeTimeStep()`

- **Role:** computes fixed or stability-limited pseudo-time steps.
- **Equations:** `dt_adv ~ h/|u|`, `dt_diff ~ h²/(2D)`, safety-scaled minimum.
- **Finding:**
  - **CC:** fixed-dt blanket path.
  - **CC/S:** dimensional diffusivity choice uses `nu=mu/rho` for fluid and `alpha=k/(rho cp)` for thermal control; conservative maximum selection is reasonable.
  - **S:** `htype=2/3` are presented as alternatives but do not form a complete streamline/SUPG directional length implementation.
  - **MPI:** global-time mode requires `MPI_Allreduce(MIN)`; local values require owner reconciliation.

### `TimeStep::computeNodalLocalTimeStep()`

- **Role:** derives nodal minima/legacy local-step data.
- **Finding:**
  - **S:** uses legacy global `ani` rather than each incident material's true dimensional diffusivity.
  - **D:** several nodal arrays do not drive the active element operator.
  - **S:** denominator safety should fail on invalid properties instead of silently producing fallback behavior.

### `TimeStep::computeGlobalRealTimeStep()`

- **Role:** establishes global/real-time value and caps.
- **Finding:** **CC** as a pseudo-time utility; **I** if interpreted as a complete physical transient integration.

### `TimeStep::applyStep2PressureTimeStepCorrection()`

- **Role:** optional pressure-based reduction of the timestep before Step 2.
- **Finding:** **I** for two independent reasons:
  1. The candidate proportional to `h*|p|/|u|` is not dimensionally a time in SI dimensional mode.
  2. It is invoked after Step 1. Step 1 has already used the old `dt`, while pressure/velocity/thermal diagonals are recomputed with the new `dt`.
- **Current blanket exposure:** none (`step2_check=0`).
- **Required action before MPI:** disable/remove until a derived dimensionally consistent formula is placed before every step of an iteration.

### `TimeStep::updateLhsDiagonal()`

- **Role:** forms inverse momentum and thermal diagonals and time-scales pressure stiffness.
- **Equations:**
  - `elcoe2_i = dt_i / M_i` on fluid velocity nodes;
  - `elcoe2p_i = dt_i / (sum_e rho_e cp_e V_e/4)`;
  - pressure `A = sum_e dt_e H_e`.
- **Finding:** **CC** for active fixed-dt path.
- **MPI:** reverse-add nodal mass/capacity/diagonal coefficients, then owner invert; distributed pressure entries are assembled with global row ownership.

### `TimeStep::updateRealTimeTerms()`

- **Role:** declared hook for physical-time/BDF terms.
- **Finding:** **I:** empty. Together with unused history arrays, this proves the physical transient mode is not implemented.

---

## 7. CBS mathematics audit against the references

## 7.1 Formulation actually implemented

The active algorithm is best written as a simplified semi-implicit projection/CBS method:

### Step 1 — pressure-free predictor

```text
M_L/dt (u* - u^n)
  = R_convection(u^n)
  + R_characteristic(u^n)
  + R_viscous(u^n)
```

### Step 2 — pressure projection

```text
A_p p^(n+1) = b_div(u*)
A_p = sum_fluid dt_e ∫ grad(N)^T grad(N) dV
```

with strong gauge/outlet pressure rows.

### Step 3 — correction

```text
M_L/dt (u^(n+1)-u*) = -∫ N^T grad(p^(n+1)) dV
```

### Step 4 — dimensional CHT

```text
rho cp (T^(n+1)-T^n)/dt
  = -rho cp u^(n+1)·grad(T^n)
    + div(k grad(T^n)) + Q
    + characteristic convection stabilization
```

This is structurally consistent with Split A and the four-stage sequence in the references. It does **not** implement every higher-order residual term of the most general published equations.

---

## 7.2 Step 1: `MomentumAssembly`

### Internal helpers

The file contains local helpers for index/geometry validation, material checks, gradient access, nodal interpolation, convection, characteristic and diffusion assembly. They read `unkn1`, geometry and material properties and scatter into `rhs`.

### `MomentumAssembly::assembleStep1Rhs(CBSStateSI&)`

- **Role:** complete fluid-element predictor residual.
- **Inputs:** iteration-old velocity `unkn1`, shape gradients, determinant, element timestep, density/viscosity, face classification/normals and colouring.
- **Outputs:** `rhs(1:3,1:npoin)`.
- **Modified state:** RHS only, apart from an error flag.
- **Fluid/solid:** fluid elements only; solids skipped.
- **OpenMP:** colours processed serially; elements within one colour in parallel.
- **Boundary handling:** convection and natural viscous face terms are included according to `fedge`.
- **Future MPI:** read forward-halo `u^n`; assemble owned elements; reverse-add shared RHS to owners.

#### Galerkin convection

- **CC:** the code performs exact P1 tetrahedral integration for the chosen nodal interpolation of the nonlinear product and includes the corresponding surface contribution.
- **Dimensional check:** acceleration residual has consistent momentum-mass scaling.

#### Characteristic correction

- **CC/S:** the `dt/2` streamline-like term is a valid characteristic stabilization of convection.
- **S:** the general CBS residual shown in the references can include viscous and body/source residual parts. The code retains a simplified convection-only characteristic residual.
- **Interpretation:** acceptable as a documented reduced method; do not claim the full general residual without adding/validating the omitted terms.

#### Viscous diffusion

- **S:** code uses componentwise `nu ∇²u`. The references write the divergence of the full deviatoric stress
  `mu(grad u + grad u^T - 2/3 div(u) I)`.
- These are equivalent for constant viscosity and exactly divergence-free velocity in the final continuum solution, but the explicit predictor is not exactly divergence-free and variable viscosity would break equivalence.
- **No confirmed sign error** was found in the implemented Laplacian weak form.

#### Natural momentum boundary term

- **S:** the implementation reconstructs an interior normal derivative on natural faces rather than receiving a clearly prescribed total/deviatoric traction. This is not a transparent realization of the CBS traction treatment in the references, especially at pressure outlets.
- Current blanket walls are strongly no-slip and the outlet configuration has been validated empirically, but the code should state the exact outlet traction assumption.

#### Body force/buoyancy

- **I relative to advertised controls:** no gravity, Boussinesq or other body-force term is assembled despite stored Rayleigh/Richardson-related controls. Forced convection is the actual active model.

#### OpenMP defect

- **I:** a shared `bool bad_detj` can be written by multiple threads without reduction/atomic protection. All writers store `true`, but this is still a C++ data race and undefined behavior.

### `MomentumAssembly::applyRealTimeMomentumTerm()`

- **Finding:** **I/D:** placeholder only; no physical-time contribution.

### Step-1 correspondence to the papers

| Published component | Active code |
|---|---|
| pressure removed in Split A | **CC** |
| Galerkin convection | **CC** |
| characteristic `dt/2` stabilization | **CC**, simplified residual |
| full deviatoric viscous stress | **S**, replaced by component Laplacian |
| body/source terms | **I/D**, absent |
| prescribed traction treatment | **S**, implicit/reconstructed assumption |

---

## 7.3 Step 2: `PressureAssembly`

### `buildElementPressureTerms()`

- **Role:** computes compact unscaled P1 pressure stiffness for fluid tetrahedra.
- **Equation:** `H_ab^e = V_e grad(N_a)·grad(N_b)`.
- **Outputs:** four diagonals and six unique off-diagonals per element.
- **Finding:** **CC**; symmetry and positive-semidefinite element form are correct.
- **OpenMP:** element-private writes are safe.
- **MPI:** owned fluid elements only.

### `buildGlobalPressureTerms()`

- **Role:** initializes/accumulates global pressure structures before timestep scaling.
- **Finding:** **CC** for monolithic state; name is somewhat misleading because active time scaling is completed in `TimeStep`.
- **MPI:** becomes local row preallocation/assembly rather than one global array.

### `assembleStep2Rhs()`

- **Role:** assembles weak divergence of predictor velocity plus normal boundary contribution.
- **Inputs:** predictor `unkno`, geometry, face normals/classification, local dt assumptions.
- **Output:** scalar pressure RHS (stored in the pressure work vector).
- **Finding:** **CC** for the implemented theta1=theta2=1 absolute-pressure projection.
- **S:** comments sometimes call the unknown a pressure increment/correction, but the solver stores and solves the newest absolute/gauge pressure directly.
- **S:** `theta1` is not used and `theta2` is fixed to one, so the general weighted CBS Step 2 in the references is not available.
- **OpenMP:** active RHS assembly is serial.
- **MPI:** reverse-add nodal RHS; constrained owner rows; distributed global solve.

### `assembleRhs()`

- **Role:** compatibility wrapper.
- **Finding:** **CC/D**.

### Step-2 correspondence

- The pressure Laplacian and divergence projection are correct for the simplified semi-implicit method.
- The full general Split-A equation includes theta weighting and pressure-history/stabilization details not exposed by this code. This is **S**, not a detected sign error.

---

## 7.4 Step 3: `Steps::step3SemiImplicit()`

- **Role:** element pressure-gradient correction over fluid elements.
- **Input:** newest `pres`, P1 shape gradients, determinant, inverse mass/time diagonal.
- **Output/modified state:** velocity `unkno`; temporary vector `rhs`.
- **Equation:** each local node receives `-(V/4) grad(p)`, then multiplies by `dt/M_L`.
- **Finding:**
  - **CC:** first-order projection sign, `V/4` integration and fluid-only domain are correct.
  - **CC:** strong velocity BCs are reapplied after correction.
  - **S:** the general CBS Step 3 in the supplied references contains a characteristic higher-order pressure term (often written with `dt²/2 P p^n`). It is absent here.
  - **S:** assembly is serial although naturally colour-parallel.
- **Dimensional check:** `V grad(p)` divided by density-weighted/lumped momentum scaling must be interpreted consistently with how pressure and density enter the active equation. The current dimensional implementation has produced physically scaled validated results, but the code comments should explicitly state whether the stored pressure operator is divided by constant density or uses a kinematic pressure convention. The present source mixes absolute-pressure-looking initialization with projection/gauge use; this deserves a unit test.
- **MPI:** forward halo pressure; owned-element correction; reverse-add nodal correction; owner update/BC; forward halo velocity.

---

## 7.5 Step 4: `EnergyAssembly`

### `EnergyAssembly::assembleStep4Rhs()`

- **Role:** complete explicit dimensional CHT residual.
- **Inputs:** iteration-old `temperature1`, corrected velocity, geometry, timestep, material properties, heat-source and heat-flux data.
- **Output:** `rhs1`.
- **Modified state:** thermal RHS only.
- **Fluid/solid:** convection/stabilization only in fluid; diffusion/source in all elements.
- **OpenMP:** volume assembly is serial; nodal update in `Steps` is parallel.
- **Future MPI:** forward halo `T^n` and corrected `u`; owned-element assembly; reverse-add RHS/capacity; owner update; forward halo temperature.

#### Thermal capacity

- **CC:** lumped `rho cp V/4` and update `dt/(rho cp mass)` are dimensionally correct.

#### Fluid convection

- **CC:** exact P1 integration of `rho cp N_a (u·grad T)` for linearly interpolated velocity is represented by the tetrahedral `V/20` coefficients.

#### Characteristic/SUPG-style term

- **CC/S:** `dt/2 rho cp V (ubar·grad N_a)(ubar·grad T)` is a conventional characteristic streamline stabilization.
- **S:** it is a simplified convection-residual form rather than the complete Step-4 residual in the general compressible-energy equation.

#### Diffusion

- **CC:** `-k V grad(N_a)·grad(T)` is the correct weak internal residual for both fluid and solid.

#### Source

- **CC:** `Q V/4` for constant volumetric source.
- **D:** a separate `source_solid` control is inactive; the actual source is the per-material `Q` from `.matprop`.

#### Heat flux BC 532

- **CC:** positive inward heat flux contributes `q A/3` to each triangular face node, consistent with `∫ N_a q dA`.

#### Conformal interface BC 901

- **CC:** shared temperature nodes impose temperature continuity; summing fluid and solid diffusion contributions gives weak flux balance naturally. No duplicate interface heat source is required.

#### Property model

- **S:** properties are constant by material. Temperature-dependent `mu`, `k`, `cp` are not implemented.

### `EnergyAssembly::applyRealTimeEnergyTerm()`

- **Finding:** **I/D:** placeholder only.

### `Steps::step4Energy()`

- **Role:** applies `T^{n+1}=T^n+elcoe2p*rhs1` and temperature BCs.
- **Finding:** **CC** for pseudo-time/explicit Step 4; **I** as a true-transient update because history terms are absent.

---

## 8. Linear algebra and pressure-solver audit

## 8.1 `MatrixVectorCalc`

### `pressureMultiply()`

- **Role:** matrix-free multiplication by the compact pressure operator.
- **Inputs:** nodal vector, `pdiag`, element `gstif`, connectivity, fluid mask, fixed-pressure markers and colours.
- **Output:** nodal product.
- **Finding:**
  - **CC:** diagonal and six symmetric off-diagonal contributions are assembled correctly.
  - **CC:** fixed-pressure rows return identity behavior for the native correction-space CG.
  - **CC:** element colouring prevents shared-node races.
  - **D/performance:** traverses all element colours and tests material inside, despite pressure-only colour lists already existing.
- **MPI:** each MatVec needs forward halo of trial vector, owned-element multiplication, reverse-add output and global reductions in CG.

### `pressureMultiplyVector()`

- **Role:** vector/multiple-component variant.
- **Finding:** **D/S:** not used by the active pressure path and does not reproduce all scalar fixed-row behavior. Do not use it for constrained pressure without completing its BC semantics.

## 8.2 Native `ConjugateGradient`

### `solvePressure(CBSStateSI&)`

- **Role:** solves the symmetric pressure projection system.
- **Inputs:** compact pressure operator, pressure RHS, old/current pressure, fixed nodes and tolerances.
- **Outputs:** updated `pres` and convergence result.
- **Modified state:** CG work vectors and pressure.
- **Algorithm:** standard preconditioned CG with optional Jacobi.
- **Finding:**
  - **CC:** conventional residual, search direction, alpha/beta and Jacobi formulas.
  - **CC:** fixed-pressure values are handled using an initial residual and zero correction directions on constrained DOFs; free/fixed coupling is retained consistently.
  - **S:** CG work arrays/masks are reconstructed every time step, adding significant overhead.
  - **S:** breakdown test uses `abs(pAp)`; for an SPD system, `pAp <= 0` should be reported explicitly as loss of positive definiteness rather than accepting a negative value of sufficient magnitude.
  - **S:** max-norm reduction is serial.
  - **Performance:** Jacobi is mathematically valid but weak for the 3D Poisson-like system, explaining the large iteration count.
- **OpenMP:** vector loops/dot products and matrix-vector product are threaded safely.
- **MPI:** local dot sums require `MPI_Allreduce(SUM)`; max residual requires `MAX`; pressure MatVec needs halo exchange. A hand-built distributed native CG is possible but PETSc is the safer primary path.

## 8.3 `PetscPressureSolver`

### cached system construction/helpers

- **Role:** compress active fluid-connected pressure nodes into PETSc indices, assemble a sequential sparse matrix, set constraints, create KSP and preconditioner.
- **Current communicator:** `PETSC_COMM_SELF`.
- **Current object types:** sequential AIJ matrix/vector.
- **Finding:**
  - **CC:** active-DOF compression and symmetric treatment of fixed-pressure columns/rows are conceptually correct for the sequential matrix.
  - **I:** this is not a distributed solver.
  - **I:** lazy PETSc initialization inside a rank-0-only numerical call is incompatible with the intended MPI-world solve.
  - **I:** cache validity tests topology/counts but not the actual time-scaled matrix coefficients. If `dt` changes, `pdiag/gstif` change while the cached PETSc matrix/AMG hierarchy can remain stale.
  - **S:** fixed-node cache comparison appears count-based rather than a complete global-ID/value signature.
  - **S:** fixed estimate of about 80 nonzeros per row is safe but memory-heavy and not derived from local adjacency.

### `PetscPressureSolver::solvePressure()`

- **Role:** loads RHS/initial guess, calls KSP and copies pressure back.
- **Krylov:** CG.
- **Preconditioner:** HYPRE BoomerAMG when available, otherwise PETSc GAMG.
- **Finding:**
  - **CC:** AMG is an appropriate class of preconditioner for this Poisson-like problem and should be far stronger than Jacobi.
  - **S:** CG mathematically requires an SPD operator and an effectively symmetric positive preconditioner. BoomerAMG defaults are not guaranteed to be symmetric under all relaxation/coarsening choices. Configure symmetric relaxation/coarse solve explicitly, verify `KSPCG` compatibility, or use a flexible/non-CG Krylov method if the selected AMG application is nonsymmetric.
  - **S:** input controls named for the native CG do not fully govern PETSc because command-line PETSc options can override KSP/PC behavior. This is useful but must be reported.
  - **CC:** recomputed true algebraic residual is valuable.
  - **I/MPI:** result is copied only within rank 0's full vector.

### `shutdown()`

- **Role:** destroys cached PETSc objects.
- **Finding:** **CC** sequentially; **MPI:** destruction/finalization must occur collectively/consistently after all ranks finish using PETSc.

### Required distributed PETSc form

```text
communicator     = MPI communicator used by the solver
matrix           = MatMPIAIJ
vectors          = VecMPI
row ownership    = active pressure DOFs owned by each rank
global numbering = contiguous PETSc numbering mapped from global mesh-node IDs
assembly         = owned element contributions, PETSc off-process insertion allowed
constraints      = symmetric row/column elimination or explicit RHS elimination
KSP              = KSPCG only after SPD verification
PC               = HYPRE BoomerAMG with symmetric options, or a compatible flexible KSP
```

For fixed `dt`, matrix values and hierarchy may be reused. For changing `dt`, matrix values must be updated and `KSPSetOperators` called; `SAME_NONZERO_PATTERN` can preserve sparsity while preconditioner reuse is controlled intentionally.

## 8.4 `BandedMatrix` / `BandedGaussianSolver`

- **Role:** legacy small/banded direct solver utility.
- **Finding:**
  - **CC:** storage/addressing and no-pivot Gaussian procedure are conventional for a correctly sized banded system.
  - **S:** `BandedMatrix::add` silently ignores entries outside the declared band, which can conceal a wrong bandwidth.
  - **I/D:** active CBS3D Step 2 rejects the banded option; state storage is not ready for the real 3D bandwidth.
  - **MPI:** no role in the distributed large system.

## 8.5 `CSRMatrix`

- **Role:** generic triplet-to-CSR construction and serial multiplication.
- **Finding:** **CC** sorting/merging and CSR MatVec; **D:** not used by active pressure solver.
- **S:** writable internal arrays can invalidate CSR invariants if misused.
- **MPI:** not a replacement for PETSc's distributed matrix unless a full distributed CSR layer is deliberately built.

---

## 9. OpenMP safety and performance audit

| Kernel/path | Safety | Finding |
|---|---|---|
| shape derivatives | element-private writes | **CC** |
| pressure element stiffness | element-private writes | **CC** |
| coloured momentum scatter | no shared node within colour | **CC**, except shared `bad_detj` flag race |
| pressure MatVec | coloured scatter | **CC** |
| native CG vector/dot loops | OpenMP reductions/independent indices | **CC** |
| Step 3 pressure correction | serial | safe but unthreaded |
| Step 4 volume assembly | serial | safe but unthreaded |
| nodal Step 1/4 updates | independent nodes | **CC** |
| boundary loops | serial | safe; priority/order concerns remain |
| convergence | partly serial/threaded | safe; repeated mask construction costly |

### Confirmed OpenMP defect

`MomentumAssembly` writes a shared `bad_detj` Boolean from parallel iterations without `reduction(||:bad_detj)`, atomic, or thread-local reduction. This must be fixed before claiming a race-free OpenMP implementation.

### Recommended hybrid policy

First validate **40 MPI ranks × 1 OpenMP thread**. Only after bitwise/tolerance agreement should ranks-per-node and threads-per-rank be tuned. MPI calls remain outside worker regions under `MPI_THREAD_FUNNELED`; local element colouring remains the simple race-control method.

---

## 10. Convergence and transient audit

### `Convergence::evaluate()`

- **Role:** builds relative change, current L2 and RHS-scale diagnostics for velocity, pressure and temperature.
- **Domain masks:** velocity fluid-only excluding solid-touch/interface nodes; pressure fluid-connected; temperature full CHT domain.
- **Finding:**
  - **CC:** relative change form and domain separation are reasonable.
  - **CC:** pressure includes interface fluid-connected nodes, matching pressure assembly support.
  - **S:** velocity/thermal “RHS residuals” are reconstructed from update divided by inverse diagonal. They are useful update-equivalent diagnostics, not full PDE residuals.
  - **S:** pressure relative norm depends on pressure gauge/offset; fixed outlet gauge makes current case usable.
  - **S:** invalid/zero inverse weights can fall back rather than hard-fail, potentially hiding preprocessing defects.
  - **Performance:** material masks are rebuilt every iteration.
- **MPI:** every squared norm is local sum followed by global `Allreduce(SUM)`; maximums use `MAX`.

### residual accessors

- `velocityResidual()`: maximum of three component relative residuals — **CC**.
- `pressureResidual()`: configured pressure measure — **CC/S**, dependent on gauge and selected criterion.
- `temperatureResidual()`: relative thermal residual — **CC**.

### `steadyStateReached()`

- **Role:** applies configured tolerances/minimum iteration.
- **Finding:** **CC** when at least one criterion is enabled.
- **S:** a configuration with all tests disabled can report convergence after the minimum iteration. Validate that at least one stopping criterion is active.
- **MPI:** one globally identical Boolean.

### Physical transient model

The supplied AC references recover true transience by adding physical-time terms and using pseudo-time convergence inside each real time step. The active code has:

- history arrays,
- `transient_on`, `delt`, `rtime`,
- empty real-time assembly hooks,
- no BDF term,
- no inner pseudo-time convergence loop per physical step,
- no history rotation.

Therefore **`transient_on=1` is confirmed not to implement a valid physical transient CBS/dual-time method**. It must be disabled in input validation until implemented and separately validated.

---

## 11. Output, diagnostics and profiling

## 11.1 `Post`

### banner/stage/progress/final-summary functions

- **Role:** human-readable run diagnostics.
- **Finding:** **CC**; ANSI formatting can clutter batch logs but is cosmetic.
- **MPI:** root-only console reporting; aggregate timing/statistics first.

### `initialiseRunOutputs()` / `writeResidualRow()`

- **Role:** create residual CSV and append diagnostics.
- **Finding:** **CC** serially.
- **MPI:** rank 0 writes one global row after reductions.

### `writeSolution()` and VTU/PVD helpers

- **Role:** ASCII VTK unstructured-grid output and time collection.
- **Finding:**
  - **CC:** connectivity/cell types and core point/cell fields are structurally valid.
  - **I diagnostic:** a `safe_value` helper writes zero for NaN/Inf. This hides numerical failure in visualization instead of stopping/reporting it.
  - **S/performance:** reconstructs material node masks for each output.
  - **S/scalability:** one giant ASCII VTU from 2.4 million tetrahedra is serial, slow and memory/I/O heavy.
  - **S:** steady iterations can all carry the same physical PVD time (`rtime=0`), making the collection ambiguous.
  - **D/I:** Tecplot/Nusselt controls are present but no complete implementation exists.
  - **S:** pressure output is gauge-like but `p_ref` is not added; initial `p0=101325` can look absolute.
- **MPI:** one VTU piece per rank containing owned cells, plus a root-written PVTU and PVD. Ghost points may be duplicated between pieces, but owned cells must not be duplicated. Binary appended VTU or parallel HDF5 should be considered later.

## 11.2 `SolverProfiler`

- **Role:** RAII timing by solver section and pressure iteration statistics.
- **Finding:** **CC**.
- **D:** a separate Step-2 pressure-solve section exists but current instrumentation may combine RHS/solve timing.
- **MPI:** report rank maximum (critical path), average and imbalance, not rank-0 time alone.

---

## 12. Final defect and risk register

### 12.1 Critical/blocking before distributed release

| ID | Status | Defect | Current blanket exposure | Required resolution |
|---|---|---|---|---|
| C-01 | **I** | MPI ranks do not execute distributed Steps 1–4 | Always; current MPI is audit/bootstrap only | Rank-local state, ownership, halos, reductions and distributed output |
| C-02 | **I** | PETSc uses `PETSC_COMM_SELF` and sequential objects | Avoided because native CG selected | `PETSC_COMM_WORLD`/solver communicator, MatMPIAIJ, VecMPI |
| C-03 | **I** | Physical transient path is not implemented | Avoided (`transient_on=0`) | Disable option or implement dual-time/BDF completely |
| C-04 | **I** | PETSc cached matrix becomes stale when timestep coefficients change | Avoided by fixed dt | Include coefficient/version in cache; update/reassemble values and KSP operators |
| C-05 | **I** | PETSc initialization is rank-0 lazy rather than all-participating-rank program lifecycle | Avoided by current native path | Initialize/finalize PETSc from `main` on all solver ranks |
| C-06 | **I** | Step-2 timestep correction is dimensionally wrong and creates mixed dt within one iteration | Avoided (`step2_check=0`) | Disable/remove; rederive and apply before Step 1 |

### 12.2 Major mathematical/software issues

| ID | Status | Finding | Resolution |
|---|---|---|---|
| M-01 | **I** | BC 504 is implemented but rejected by MeshIO | Single authoritative BC table; include 504 or remove dead implementation |
| M-02 | **I/S** | `theta1` ignored; theta2 fixed to one | Document fixed formulation or implement weighted equations consistently |
| M-03 | **I** | OpenMP race on `bad_detj` | OpenMP Boolean reduction or prevalidate determinant serially |
| M-04 | **S** | Componentwise Laplacian replaces full deviatoric stress | Document assumption; add full tensor form if required, then benchmark |
| M-05 | **S** | Higher-order characteristic pressure term absent in Step 3 | Confirm intended reduced CBS formulation against legacy code; document or implement/validate |
| M-06 | **I/D** | Buoyancy/body-force controls exist but Step 1 has no force term | Reject unsupported mode or implement body-force residual |
| M-07 | **S** | Natural momentum/traction BC is not transparently the reference formulation | State exact outlet traction model; add traction benchmark |
| M-08 | **S** | Material phase is inferred only from ID zero/nonzero | Store explicit phase enum from material data |
| M-09 | **I contract** | explicit CBS, banded solver and restart are accepted but not implemented | Fail during input parsing with precise message or implement |
| M-10 | **S** | Absolute/gauge pressure semantics and `p_ref` are inconsistent | Define stored pressure and output pressure clearly; validate units |
| M-11 | **I diagnostic** | NaN/Inf output is silently replaced by zero | Throw or write NaN and terminate with node/field context |
| M-12 | **S** | KSPCG + default AMG symmetry is not guaranteed | Symmetric AMG settings/SPD verification or compatible flexible Krylov |
| M-13 | **I/MPI** | no halo exchange/global reductions in numerical code | Implement explicit forward/reverse communication and collectives |

### 12.3 Medium/conditional findings

| ID | Finding |
|---|---|
| R-01 | `htype=2/3` are partial characteristic-length alternatives, not a complete streamline formulation. |
| R-02 | `rem_deltp`, parts of beta controls and several config fields are dormant. |
| R-03 | curved mass-flow inlet nodal direction can depend on face iteration order. |
| R-04 | multi-plane symmetry nodes can depend on projection order. |
| R-05 | outlet backflow clipping changes the physical/mass balance and should be optional/reported. |
| R-06 | duplicate/missing mesh IDs are not fully audited by C++ input code. |
| R-07 | native CG should reject nonpositive `pAp`, not only near-zero magnitude. |
| R-08 | convergence configuration can be vacuous if all criteria are disabled. |
| R-09 | pressure-only colouring is generated but not used. |
| R-10 | Step 3 and Step 4 volume kernels remain serial. |
| R-11 | material/interface masks are repeatedly rebuilt. |
| R-12 | serial ASCII output is not scalable. |
| R-13 | banded add silently drops out-of-band entries. |
| R-14 | initial velocity magnitude floor makes rest nonzero in diagnostics. |

---

## 13. MPI integration plan that preserves the validated solver

The MPI work should not rewrite the CBS mathematics. It should replace the global storage/assembly mechanics while preserving the exact step order and element formulas.

## Phase 0 — freeze the reference

1. Keep the present serial/OpenMP executable unchanged as a reference target.
2. Record for small deterministic cases:
   - nodal `u,v,w,p,T` after each individual CBS step;
   - pressure RHS and matrix action;
   - residual history;
   - material/BC counts and integral mass/energy checks.
3. Add a one-rank MPI mode that must match the reference within roundoff before adding more ranks.

**Reason:** distributed reverse-add ordering can change floating-point summation. Validation must use both tight norms and selected exact/topological invariants.

## Phase 1 — rank-local mesh/state model

Add simple, explicit metadata:

```text
PartitionMetadata
    rank, nranks
    n_owned_nodes, n_ghost_nodes, n_local_nodes
    n_owned_elements
    local_to_global_node
    global_to_local_node
    node_owner_rank
    neighbour ranks
    send/receive global-node lists
```

Rules:

- Every tetrahedron has exactly one owning rank.
- A rank stores all nodes of its owned tetrahedra.
- A local node is owned by exactly one rank; other copies are ghosts.
- Numerical arrays use local indices `1..n_owned+n_ghost`.
- Only owners perform final nodal updates and strong BC resolution.
- Owned cells, not ghost cells, contribute to global element counts and output.

This remains readable and close to the current array-based style.

## Phase 2 — rank-local input

- Read `rank_XXXX` mesh/material/BC metadata generated by the verified 40-partition pipeline.
- Preserve original global node and element IDs for PETSc numbering, audit and output.
- Validate collectively:
  - sum of owned tetrahedra = 2,383,581;
  - unique owned nodes = 407,458;
  - sum of owned physical faces = 28,328;
  - material and BC totals match the verified values.

## Phase 3 — two explicit halo operations

Implement one simple `HaloExchange` class with two operations:

1. **Forward exchange:** owner value → all ghosts. Used for `u`, `p`, `T` and derived fields read by element loops.
2. **Reverse additive exchange:** ghost-partition contribution → owner sum. Used for assembled nodal RHS, lumped mass, thermal capacity and diagonals.

Do not hide these calls inside generic operator overloading. Put them visibly between numerical stages.

## Phase 4 — distributed preprocessing

- Shape gradients, determinants, element lengths and owned-element colouring are local.
- Reverse-add nodal lumped mass and thermal capacity.
- Reconcile wall/interface/BC node masks by owner reduction.
- Compute mass-flow inlet area from owned inlet faces and `MPI_Allreduce(SUM)`.
- Select a global fallback pressure reference deterministically if no pressure outlet exists.

## Phase 5 — distributed Steps 1, 3 and 4

### Step 1

```text
forward u^n
zero local RHS
assemble owned fluid elements using local colours
reverse-add momentum RHS
owner update u*
apply owner BCs
forward u*
```

### Step 3

```text
forward p^(n+1)
assemble owned fluid-element pressure-gradient correction
reverse-add correction RHS
owner update u^(n+1)
apply owner BCs/symmetry/backflow
forward u^(n+1)
```

### Step 4

```text
forward T^n and corrected u
assemble owned fluid/solid elements and owned heat-flux faces
reverse-add thermal RHS
owner update T^(n+1)
apply owner thermal BCs
forward T^(n+1)
```

First reproduce the active serial formulas exactly. Mathematical improvements such as the full deviatoric stress or omitted characteristic pressure term must be separate, independently validated changes.

## Phase 6 — distributed PETSc Step 2

1. Initialize PETSc in `main` on every solver rank.
2. Build a global pressure DOF numbering for owned active fluid-connected nodes.
3. Create `MatMPIAIJ` with row ownership matching node ownership.
4. Derive diagonal/off-diagonal preallocation from local element adjacency.
5. Insert owned-element `dt_e H_e` contributions.
6. Apply fixed-pressure constraints symmetrically.
7. Create `VecMPI` RHS/solution.
8. Use `KSPCG` only after verifying operator and PC symmetry.
9. Prefer HYPRE BoomerAMG with explicitly symmetric relaxation/coarse settings; retain GAMG fallback.
10. Expose `-ksp_*` and `-pc_hypre_boomeramg_*` options in the Slurm command.
11. Fixed dt: reuse values/hierarchy only after checking coefficient invariance.
12. Variable dt: update matrix values and KSP operators every required iteration; reuse only the nonzero pattern unless tested otherwise.

## Phase 7 — global convergence and timestep reductions

Use explicit collectives:

```text
sum of squares / dot products : MPI_Allreduce SUM
maximum absolute residual     : MPI_Allreduce MAX
global stable timestep        : MPI_Allreduce MIN
converged Boolean             : MPI_Allreduce logical AND
```

The native distributed CG would need these reductions each iteration. PETSc already manages its internal Krylov collectives.

## Phase 8 — parallel output

- Each rank writes one VTU piece containing owned tetrahedra and the required local points.
- Rank 0 writes the PVTU descriptor and PVD series.
- Do not output partition-interface faces as physical boundary faces.
- Never silently zero nonfinite values.
- After correctness, move from huge ASCII files to appended-binary VTK or parallel HDF5/XDMF if necessary.

## Phase 9 — validation ladder

| Test | Required result |
|---|---|
| 1 MPI rank, tiny cavity | serial/OpenMP agreement after every step |
| 2 ranks, two-tetra mini mesh | exact ownership/halo and pressure projection checks |
| 2–4 ranks, small CHT interface | temperature/flux continuity and energy balance |
| 40 ranks, blanket, two iterations | global counts exact; field differences within roundoff tolerance |
| pressure solver comparison | native CG and PETSc converge to the same algebraic solution |
| scaling | report max/average rank timings and communication fraction |

## Phase 10 — only after parity: optional mathematical improvements

Treat each of the following as a separate branch/test, not part of basic MPI conversion:

- full deviatoric-stress tensor;
- omitted characteristic pressure correction term;
- body force/buoyancy;
- true physical transient/dual-time stepping;
- variable properties;
- general curved inlet constraint handling.

---

## 14. Recommended implementation order

1. Correct the non-distributed confirmed defects that can invalidate tests:
   - reject transient mode;
   - disable Step-2 pressure timestep correction;
   - repair BC 504 contract;
   - fix `bad_detj` OpenMP race;
   - stop hiding NaN/Inf.
2. Add rank-local state and reader without changing Steps.
3. Add forward/reverse halo exchange and distributed preprocessing.
4. Parallelize Steps 1, 3 and 4 around the exact existing formulas.
5. Move PETSc lifecycle to `main` and implement distributed Step 2.
6. Add global convergence/timestep reductions.
7. Add PVTU/PVD parallel output.
8. Validate 1, 2, 4 and 40 ranks.
9. Profile and tune MPI/OpenMP layout.
10. Only then revisit mathematical simplifications.

---

## 15. Public file/function coverage checklist

The following active public functions/classes were individually traced. Anonymous helper routines were audited with their owning function and are summarized in the relevant sections above.

| File | Functions/classes | Status summary |
|---|---|---|
| `Array.hpp` | `Array1D/2D/3D` constructors, resize, fill, accessors, offset/check functions | **CC** |
| `CBSStateSI.hpp` | `initialise_local_topology`, `set_problem_sizes` | **CC/S/MPI** |
| `MeshIO` | `readAll`, `readSizes`, `readMeshFile`, `readBoundaryFile`, `readParameterFile`, `readMaterialFile`, `readMaterialPropertyFile`, `initialiseFields` | mixed; BC504 and unsupported-mode contracts flagged |
| `Preprocess` | boundary validation, shape gradients, face matching, normals, mass, face classification, size, walls, mass-flow velocity, velocity magnitude, pressure nodes | geometry mostly **CC**; minor/parallel issues flagged |
| `Boundary` | temperature, velocity, pressure, symmetry, outlet backflow | active blanket path **CC**; corner/curved-face/order risks flagged |
| `TimeStep` | timestep, nodal local timestep, global/real timestep, Step-2 correction, LHS diagonal, real-time terms | fixed-dt **CC**; Step-2 correction and transient **I** |
| `MomentumAssembly` | Step-1 RHS, real-time hook | simplified CBS **S**; OpenMP race **I** |
| `PressureAssembly` | element/global terms, Step-2 RHS, wrapper | projection path **CC/S** |
| `EnergyAssembly` | Step-4 RHS, real-time hook | dimensional CHT **CC**; true transient absent |
| `Steps` | selectors, four semi-implicit steps, explicit rejection | active sequence **CC**; accepted unsupported modes **I contract** |
| `ConjugateGradient` | pressure solve and helper reductions/preconditioning/constraints | **CC/S**, serial-global only |
| `MatrixVectorCalc` | scalar and vector pressure multiply | scalar **CC**; vector constrained semantics incomplete/dormant |
| `PetscPressureSolver` | solve, shutdown and cached-system helpers | sequential only; stale-cache and lifecycle defects |
| `BandedMatrix` | all storage/access methods | utility **CC/S**, active path absent |
| `BandedGaussianSolver` | `solve` | conventional no-pivot utility, inactive |
| `CSRMatrix` | clear/size/triplets/accessors/MatVec | **CC**, inactive |
| `Convergence` | evaluate, three residual accessors, steady test | **CC/S**, needs global reductions |
| `Post` | banner, summaries, stage output, residual, scheduling, VTU/PVD, progress, final summary, live plot | serially functional; nonfinite/output scalability defects |
| `SolverProfiler` | timers, iteration reset, data setters, report, names/total | **CC**, future distributed aggregation |
| `Solver` | constructor/state/context/run/initialize/prepare/advance/update/stop tests | order **CC**; distributed/transient limitations |
| `main.cpp` | broadcast and `main` | MPI bootstrap **CC**; numerical distribution/PETSc lifecycle absent |
| `CMakeLists.txt`, `Makefile` | all build branches | serial/OpenMP **CC**; PETSc linkage and MPI source set need revision |

---

## 16. Final release judgement

### Present serial/OpenMP steady blanket solver

**Conditionally releasable as a validated steady, laminar, incompressible, constant-property CHT solver**, provided its formulation is described accurately and unsupported modes are disabled rather than advertised.

### Present MPI executable

**Not a parallel solver.** It is a valid MPI bootstrap and mesh/halo audit environment only.

### Present PETSc path

**Useful sequential prototype, not safe as the final distributed pressure solver.** Fixed-dt sequential use can work, but changing-dt cache reuse and MPI lifecycle/communicator design must be corrected.

### Mathematical integrity

No fundamental sign or geometry error was identified in the active pressure projection, CHT diffusion, source or heat-flux terms. The principal mathematical qualifications are the simplified viscous stress, reduced characteristic residuals, omitted higher-order Step-3 pressure term, absent body force and absent real-transient formulation.

### Immediate decision

Proceed to MPI integration only after freezing the current reference and correcting/disabling the confirmed conditional defects. Preserve the existing readable array-and-loop style. Do not mix the explicit CBS-AC equations into the semi-implicit pressure-Poisson code during parallelisation.

