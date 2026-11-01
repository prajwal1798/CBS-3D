# CBS3D handover: state of the code, the blanket study and the implicit-solver plan

Written 1 November 2026 for the next coding agent (OpenCode). Read this whole file before touching the code. Owner: P. Bharadwaj (EngD, Swansea University, supervisor Prof. P. Nithiarasu). Repository: <https://github.com/prajwal1798/CBS-3D>.

---

## 0. TL;DR

- **Code**: CBS3D++_SI, a C++ finite element solver (P1 tetrahedra) for incompressible flow + conjugate heat transfer (CHT), Characteristic-Based Split (CBS) scheme, MPI domain decomposition, PETSc/hypre for the pressure equation.
- **Branch to use**: `claude/blanket` on GitHub. It is complete only after the two patches in section 1.2 have been applied on SUNBIRD2 and pushed.
- **Current study**: laminar CHT of a fusion breeder-blanket slice (water + EUROFER97, 3.2 M tets). The steady flow never converges (residual plateau). Best current explanation: genuinely unsteady separated flow behind the outer bend of channel 1 (evidence in section 7).
- **Main bottleneck**: the scheme is explicit except the pressure, so steady CHT needs 10^5+ iterations (heat creeps through the steel one element per iteration).
- **First implicit piece exists**: `SteadyEnergySolver` (implicit steady conjugate energy on a frozen velocity, PETSc GMRES, with a bounded AFC option). It is verified on small cases, not yet run on the blanket.
- **Next goal**: an implicit/accelerated solver built on PETSc (KSP/SNES/TS). Roadmap in section 10.

---

## 1. Repository, branches and what is on GitHub

### 1.1 Branch map

| Where | Branch / commit | Content |
|---|---|---|
| GitHub `origin/main` | as of 21 Oct 2026 | older; **27 commits behind** the work below |
| GitHub `claude/blanket` | `f50de6b` (pushed 1 Nov 2026 from SUNBIRD2) | everything up to *running time averages* (see 1.3) |
| SUNBIRD2 `$CBS3D_CODE` | `claude/blanket` | same + whatever patches have been applied there |
| Claude sandbox (not reachable) | `workflow/linux-sunbird-case-pipeline` | identical content; hashes differ because SUNBIRD2 got the commits via `git am` |
| Windows tower `C:\PhD\CFD\CBS3D_MASTER` | local `main` `ea7cd29`, 63 commits ahead of `origin/main`, **never pushed** | the owner's consolidated master. Its relation to `claude/blanket` was never reconciled. **Ask the owner before merging anything.** |

GitHub also has `pr/1..4` refs, plus older local branches (`build/linux-petsc-mpi`, `integrate`, `style/classic-coding-rulebook`).

### 1.2 Not yet on GitHub (do this first)

Two commits exist only as patch files in the owner's SUNBIRD2 home directory:

| Patch file | sha256 (prefix) | Commit (sandbox) | What |
|---|---|---|---|
| `~/cbs3d-afc-fluid-only.patch` | `7605a146a39b…` | `b389412` | temperature AFC limits fluid elements only; self-check runs with AFC |
| `~/cbs3d-steady-energy.patch` | `8b819ee5b2eb…` | `2bb5f84` | implicit steady energy solver |

Apply them in order, then push:

```bash
source ~/cbs3d_env.sh; T=$CBS3D_CODE/jobs/sunbird/templates
bash $T/blanket.sh update ~/cbs3d-afc-fluid-only.patch     # git am + build + self-check job
bash $T/blanket.sh update ~/cbs3d-steady-energy.patch
cd $CBS3D_CODE && git log --oneline -3 && git push origin claude/blanket
```

GitHub push from SUNBIRD2 uses a fine-grained personal access token (repository CBS-3D, Contents read/write) entered as the password; `git config --global credential.helper 'cache --timeout=28800'`.

### 1.3 Commit history of this work (oldest first, on top of `origin/main`)

| Commit (sandbox) | Summary |
|---|---|
| `4501fa6` … `1cfba6f` | Linux build with PETSc; Linux/Sunbird case pipeline `.geo` → MPI rank files; scratch inventory tool |
| `719791a` … `685877d` | SUNBIRD2 environment, PETSc/CBS3D build scripts, 128-core job templates, any-N submission, `--mem-per-cpu=2G` |
| `7711331` … `d4f1e5b` | lid-driven-cavity centreline plots vs Ghia et al. (1982); residual plots; 2D-equivalent cavity verified |
| `c68fc84` | partition ladder + automatic core picking |
| `338f102` | partition exporter: full BC table in every rank's `.bco` (interior ranks failed before) |
| `99aefbe` … `1c1dc0f` | heated-panel U-channel unit cell, manifold panel, blanket slice geometry (dimensioned) |
| `b02cdeb` | **nuclear heating source** (exponential depth profile), blanket L0/L1/L2 cases, `blanket.sh`, self-check slab case, `energy_balance.py` |
| `08f28d5` | self-check: one MPI rank per physical core (`--hint=nomultithread --cpu-bind=cores`) |
| `5f47d76` | `blanket.sh status/check` read MPI progress lines |
| `597ef02` | **energy equation**: correct (diffusive) sign of the CBS streamline term; conservative advection form |
| `b4bd2bd` | blanket L0/L1/L2 use local time stepping (`ilots = 1`) |
| `4809fcf` | **LTS monotone refresh** (`CBS3D_LTS_MONOTONE=1`) |
| `378c952` | **MPI outlet backflow control** (was serial-only) |
| `ebb3614` | **running time averages** for transient MPI runs; stage L1T; VTU-as-initial-condition restart |
| `b389412` | AFC fluid-only (patch, 1.2) |
| `2bb5f84` | implicit steady energy solver (patch, 1.2) |

---

## 2. Machines and environment

### SUNBIRD2 (Supercomputing Wales, Swansea)

- Login: `login.sunbird.swansea.ac.uk` (sb-login-01/02). The owner connects with PuTTY from Windows; VPN is not needed; MFA is required. Account (Slurm) `scw7008`, user `s.2337862`.
- Environment file `~/cbs3d_env.sh`. **Source it in every new login**, and also set `T`:
  ```bash
  source ~/cbs3d_env.sh; T=$CBS3D_CODE/jobs/sunbird/templates
  ```
  It defines `CBS3D_CODE=/pixhome/s.2337862/CBS3D_SB2/code`, `CBS3D_SCRATCH=/scratch/scw7008/s.2337862`, `CBS3D_MPI_EXE`, `CBS3D_SERIAL_EXE`, `CBS3D_SRUN_MPI`, the modules, and a Python venv with gmsh.
- Nodes: **128 physical cores** (Slurm shows 256 hyperthreads). MPI jobs must use `--hint=nomultithread --cpu-bind=cores`. Without them, 4 ranks were bound to 2 hyperthreads and ran about 1000× slower.
- Blanket case folder: `$CBS3D_SCRATCH/cases/blanket_slice` (`serial/`, `partitions_40/80/128/256/`, `run_blanket_<JOBID>/`).
- Mesh: `$CBS3D_SCRATCH/meshes/slice_h5.msh` (made by the owner on the Windows tower).

### Windows tower, laptop, local builds

- Windows tower: about 20 cores, gmsh, WSL. The owner meshes here. Do not run meshing on the owner's behalf unless asked.
- Linux laptop (Ubuntu): blanket master files at `/home/victuslnx/PhD/MESH/gmsh/cases/blanket`.

### Local build for development (any Linux)

```bash
PETSC_DIR=/usr/lib/petsc cmake -S . -B b-mpi -DCMAKE_BUILD_TYPE=Release -DCBS3D_ENABLE_MPI=ON -DCBS3D_ENABLE_PETSC=ON
make -C b-mpi -j cbs3d_parallel          # -> b-mpi/cbs3dpp_si (MPI executable)
cmake -S . -B b-tests -DCMAKE_BUILD_TYPE=Release && make -C b-tests -j && (cd b-tests && ctest)   # serial + 11 unit tests
mpirun -np 2 b-mpi/cbs3dpp_si <case> <partition_root>          # MPI run
```

On SUNBIRD2: `bash $T/build_cbs3d.sh` (called by `blanket.sh update`).

---

## 3. Code architecture (what lives where)

```
include/cbs/...      headers (same tree as src/)
src/core             CBSStateSI (the global state: arrays, cfg = RunConfig)
src/io               MeshIO (.plt/.bco/.material/.matprop/.par), Post/DistributedPost (VTU/PVTU/PVD),
                     RestartIO (native checkpoints + legacy VTU import), TimeAverage
src/preprocess       Preprocess (geometry, masks, BC state, normals), NuclearHeating
src/boundary         Boundary (velocity/pressure/temperature BCs, outlet backflow), DistributedSymmetry
src/timestep         TimeStep (stability dt, LTS, elcoe2/elcoe2p lumped inverse mass/dt)
src/assembly         MomentumAssembly (Step 1), PressureAssembly (Step 2), VelocityCorrectionAssembly (Step 3),
                     EnergyAssembly (Step 4, + elementOperator), TemperatureAFC (limiter)
src/linalg           ConjugateGradient (serial), PetscPersistentDistributedPressureSystem (MPI pressure)
src/parallel         HaloExchange (owner/ghost broadcast, reverse sums, min/or reductions), Coloring
src/solver           Solver, Steps (serial), SolverDistributed* (MPI), SteadyEnergySolver, Convergence
src/turbulence       Spalart-Allmaras (+ wall models)
tools/               python: gmsh_to_cbs.py, export_rank_local_cbs.py, energy_balance.py, check_nuclear_slab.py,
                     plot_residuals_pub.py, plot_ldc_centerlines.py, ...
jobs/sunbird/templates/   blanket.sh, run_mpi.slurm, submit_auto.sh, submit_mpi.sh, prepare_partitions.sh,
                     selfcheck.sh, selfcheck_slab.slurm, build_case.slurm, build_cbs3d.sh, ...
cases/               blanket_slice, nuclear_heating_slab, cht_uchannel_check, lid_driven_cavity, flat_plate_*, ...
docs/                LINUX_SUNBIRD_RUNBOOK.md (§6.4a/6.4b), user-guide/case-files.md, MESH_FORMAT.md, ...
```

### Traps in the architecture

1. **The MPI production loop is compiled through macro wrappers.** `SolverDistributedProductionLoopRestart.cpp` does
   `#define DistributedPost RestartDistributedPost`, `#define TimeStep RestartTimeStep` and `#include "SolverDistributedProductionLoop.cpp"`. `SolverDistributedProductionLoopRestartAfc.cpp` wraps that again. Only the AFC wrapper is in `CMakeLists.txt`.
   - Edit `SolverDistributedProductionLoop.cpp` for loop changes.
   - Never name a new class `TimeStep…` or `DistributedPost…`.
   - In that file, `DistributedPost::writeSolution(s, case, k)` writes **global** iteration `iiter_total + k`.
2. Arrays are 1-based (`Array1D`/`Array2D`, position 0 unused). `intma(a, ie)` is local node `a = 1..4` of element `ie`. Element gradients are `dNkdx` via helper `grad(s, ie, dim, a)`. `detJ = 6 V`.
3. MPI data model:
   - Each rank holds its own tets plus **ghost nodes** (no ghost cells). `owned_nodes` lists the owner copies.
   - Assemble element contributions into local nodal arrays, then call `HaloExchange::sumGhostContributionsToOwners`, update the owners, then `broadcastOwnedToGhosts`.
   - MPI reproduces serial to round-off (1e-13) on the slab.
4. `std::cout` is silenced during MPI set-up (ScopedStdoutSilence). Messages that must appear there use `std::fputs` to stdout (see NuclearHeating).
5. Material: `mat_elem(ie) == 0` means fluid; anything else is solid. Node masks: `node_material_mask` (touches fluid = 1, touches solid = 2); interface nodes are no-slip.

---

## 4. Numerical method (as implemented)

### CBS, semi-implicit form (`cbs_scheme = 1`). Only the pressure is implicit.

| Step | What | Treatment | Code |
|---|---|---|---|
| 1 | intermediate velocity u* (convection, viscous diffusion, characteristic stabilisation −Δt/2 (u·∇N)(u·∇u)) | **explicit**, lumped mass: `unkno = unkn1 + rhs·elcoe2` | MomentumAssembly |
| 2 | pressure Poisson `Δt K p = b` | **implicit**: PETSc KSPCG + hypre BoomerAMG (`PetscPersistentDistributedPressureSystem`; operator weighted by element Δt, rebuilt when the LTS field is refreshed) | PressureAssembly, SolverDistributedStep2 |
| 3 | velocity correction | explicit: `unkno += rhs·elcoe2` | VelocityCorrectionAssembly, SolverDistributedStep3 |
| 4 | temperature (water + steel, monolithic CHT on conformal mesh) | **explicit**: `T = T1 + rhs1·elcoe2p` | EnergyAssembly |

- **Time step**: `dt = csafm · min(h/|u|, h²/(2D))` per element. `h` is the minimum altitude (3V/A_max). Fluid D = max(ν, α); solid has only the diffusion limit.
  - Global time step = minimum over all elements.
  - **LTS** (`ilots = 1`, steady only) uses the element's own dt. The frozen field is margin × current (margin 0.8, `CBS3D_LTS_REFRESH_MARGIN`). With `CBS3D_LTS_MONOTONE=1` a refresh may only lower values; this stopped a refresh-every-2.5-iterations limit cycle.
- **Energy (after the 597ef02 fixes)**:
  - Conservative advection `div(ρc_p u (T − T_ref))` with `T_ref` = inlet T.
  - The characteristic streamline term has the **diffusive** sign `−Δt/2 ρc_p (u·∇N_a)(u·∇T)`. It is already the integrated-by-parts form, so `+` is anti-diffusive. An external review once argued `+`; that argument is wrong for P1.
- **AFC** (`CBS3D_TEMPERATURE_AFC=1`, MPI only, `TemperatureAFC.cpp`): Zalesak-type flux limiter on the explicit Step-4 update. Without it the water temperature under/overshoots by tens of kelvin (cell Péclet number in the hundreds, Pr ≈ 6). Since `b389412` it limits **fluid elements only**; limiting the steel conduction gave 3.3 K error on the slab. `CBS3D_AFC_SOLID=1` restores the old behaviour.
- **Outlet backflow** (BC 520): inward normal velocity is clipped. The MPI path uses a rank-reconciled nodal outlet normal (`node_outlet_normal`) inside `Boundary::applyOwnedVelocityConstraints`.
- **Steady stop**: RelU, RelV, RelW < `l2norm_vel` and RelT < `l2norm_temp` after `steady_min_iterations`. These are per-iteration relative changes (Convergence.cpp), so a small dt can fake convergence. Always check the energy balance.
- Restart: native checkpoints every `CBS3D_CHECKPOINT_EVERY` iterations (`output/<case>/restart/step_N/`, same rank count). Legacy VTU import (`CBS3D_RESTART_FORMAT=legacy_vtu`) is possible too.

### Boundary-condition IDs (`docs/Boundary Condition IDs.txt`)

500 adiabatic + prescribed velocity · 501 T=1 no-slip · 502 T=0 no-slip · 503 legacy inlet u=1 · 506 symmetry · 507/508 legacy parabolic inlets · 510 velocity inlet + T · **511 mass-flow inlet + T** · **520 pressure outlet** · **530 no-slip adiabatic wall** · **532 no-slip heat-flux wall** · 901 CHT interface · 902 heat-flux marker only.

---

## 5. Case files and how they are produced

Full reference: `docs/user-guide/case-files.md` and `docs/MESH_FORMAT.md`. The owner will also send examples. Summary:

| File | Content |
|---|---|
| `<case>.plt` | mesh: header, tetrahedra (connectivity + material), nodes, boundary faces (3 nodes + BC id) |
| `<case>.bco` | boundary-ID mapping table (every rank needs the full table) |
| `<case>.material` | element → material id |
| `<case>.matprop` | material table: `id name fluid/solid rho cp k mu Qvol` (example below) |
| `<case>.par` | run controls: labelled lines, each followed by a value line; `#` lines are comments |
| `<case>_rank_XXXX.*` + `.mpi` | rank-local case files + SEND/RECV owner/ghost maps (made by the exporter) |

Pipeline:

```bash
gmsh model.geo -3 -format msh2 -o model.msh                     # owner meshes (Windows tower)
python3 tools/gmsh_to_cbs.py model.msh --case-name blanket --out-dir serial \
        --par blanket_L1.par --matprop blanket.matprop --overwrite   # serial CBS case
python3 tools/export_rank_local_cbs.py serial/blanket --parts 256 --out-dir partitions_256 --gmsh gmsh
# on SUNBIRD2 all of this is:  bash $T/blanket.sh setup $CBS3D_SCRATCH/meshes/slice_h5.msh
```

Exporter rules:
- METIS (via `gmsh -part N`) assigns tets to ranks.
- A node belongs to the lowest rank that has a tet on it.
- Local numbering: owned nodes first, then ghosts, each in ascending global id.

`blanket.matprop`:
```
2
0 Water_298K fluid 997.0 4182.0 0.606 0.00089 0.0
1 RAFM_EUROFER97 solid 7850.0 500.0 28.3 0.0 0.0
```

Key `.par` lines for the blanket (dimensional, SI; full files in `cases/blanket_slice/`):
```
ntime transient_on dtfixed dtfix iwrite                  -> 1000000  0  0  1.0e-3  1000000
CBS3D timestep controls: cbs_scheme ilots htype ...      -> 1  1  1  0  0  1  0  2.0  1.0e-6  1.0e20   (ilots=1: LTS)
Output Controls - write_output ... end_rtime             -> 0  0.0  0.0  0.0       (end_rtime used when transient_on=1)
Steady State Convergence Tolerances                      -> 1  1.0e-8  1  1.0e-9
alpha_sf k_ratio source_solid heat_flux_bc               -> 1.0  1.0  0.0  1.0e5   (BC 532 flux, W/m2)
Dimensional controls: dimensional_mode material_properties_enabled mass_flow_inlet_enabled -> 1 1 1
Reference and inlet/outlet data: p_ref model_depth inlet_mass_flow_rate inlet_density inlet_temperature outlet_pressure_gauge
                                                         -> 0.0  1.0  0.050  997.0  298.15  0.0
CBS3D output/monitor controls: ... vtu_output_every_iterations ... steady_min_iterations -> 1 10 500 0 1 50000 0.0 0 20000
# optional trailing block (L2):
Nuclear heating: nuclear_heating_on x0 y0 z0 nx ny nz nprofiles
1  0.0  0.28  0.0  0.0  -1.0  0.0  1
Nuclear heating profiles: material_id A_W_per_m3 lambda_m
1  1.6e6  0.15
```

Nuclear heating: Q = A·exp(−d/λ) at element centroids, where d is the distance from the plane through (x0, y0, z0) along n. The solver prints `HEAT INPUT …` and `HEAT_INPUT_TOTAL_W = …` at start-up.

---

## 6. Environment variables (runtime switches)

| Variable | Meaning |
|---|---|
| `CBS3D_TEMPERATURE_AFC=1` | energy flux limiter (MPI); blanket.sh sets it |
| `CBS3D_AFC_SOLID=1` | also limit solid elements (old behaviour, not recommended) |
| `CBS3D_AFC_MIN_TEMPERATURE`, `CBS3D_AFC_LOG_EVERY` | AFC floor / logging |
| `CBS3D_LTS_MONOTONE=1`, `CBS3D_LTS_REFRESH_MARGIN` | LTS refresh rules |
| `CBS3D_CHECKPOINT_EVERY`, `CBS3D_CHECKPOINT_ROOT` | native checkpoints |
| `CBS3D_RESTART=1`, `CBS3D_RESTART_ROOT`, `CBS3D_RESTART_FORMAT=native\|legacy_vtu`, `CBS3D_RESTART_ITERATION` | restart (same rank count) |
| `CBS3D_RESTART_AS_INITIAL_CONDITION=1` | allow a VTU to start a new transient run at t = 0 |
| `CBS3D_RESTART_VELOCITY_ARRAY=mean_velocity` | import the time-averaged velocity |
| `CBS3D_CONTINUATION_ITERATIONS=N` | length of a restarted segment |
| `CBS3D_AVERAGE_START_TIME` / `CBS3D_AVERAGE_START_ITERATION` | running time averages (transient, global dt) |
| `CBS3D_STEADY_ENERGY=1` (+ `_TOL`, `_MAXIT`, `_RELAX`) | implicit steady energy instead of the loop |
| `CBS3D_PAR`, `CBS3D_MATPROP` | master .par/.matprop copied to every rank by run_mpi.slurm |
| `CBS3D_VERBOSE`, `CBS3D_SA_WALL_TREATMENT`, `CBS3D_CHT_WALL_TREATMENT`, `CBS3D_WALL_B`, `CBS3D_WALL_KAPPA`, `CBS3D_MPI_PARTITION` | other |
| `PETSC_OPTIONS` | PETSc options; prefix `energy_` for the steady energy KSP |

`blanket.sh` sub-commands (SUNBIRD2):
- `update PATCH|BUNDLE`: git am / fetch, then build, then self-check.
- `setup MSH`: convert the mesh, then build the partition ladder 40/80/128/256.
- `smoke L1`: short start-up test.
- `run L0|L1|L1T|L2 [--sets N --time T]`: production run.
- `energy RUN_DIR L1|L2`: implicit steady energy on that run's last flow field.
- `check [RUN_DIR]`: mass and energy balance plus residual plot.
- `status`: your jobs and their latest progress.

---

## 7. The blanket study

### Geometry, mesh, loads

- Geometry: `cases/blanket_slice/BlanketSlice_CHT.geo`.
  - Tapered EUROFER plate: back width 350 mm, height 280 mm, sides at 85°, front fillet R6; slice thickness 30 mm.
  - Four nested U-channels, square cross-sections 15/10/10/5 mm, fillet radii 10/5/5/3.5 mm.
  - Inlets (x < 0) and outlets (x > 0) on the back plane y = 0; the first wall (heated) is at y = 0.28 m.
  - Each channel is the outer region minus the outer region shifted in by w, both with the **same fillet radius**. So the width across the bend diagonal is about √2·w (channel 1: 21 mm against 15 mm straight), a built-in expansion that promotes separation.
- Mesh h5: 3,212,185 tets (1,345,009 water, 1,867,176 steel), 544,316 nodes. Boundary faces: 511 = 1144, 520 = 1158, 530 = 49,316, 532 = 2454. Volume 2.7338e-3 m³. BC 532 area 9.231e-3 m², so 1e5 W/m² gives 923.1 W.
- Water: 997 / 4182 / 0.606 / 8.9e-4 (constant, 298 K). EUROFER97: 7850 / 500 / 28.3. Inlet 0.050 kg/s total at 298.15 K.
  - Measured split over the four inlets: 0.0256 / 0.0110 / 0.0110 / 0.0025 kg/s.
  - Re ≈ 1870 / 1250 / 1250 / 620 (channel 1 is the outermost, nearest the first wall).
- Stages:
  - L0: flow only.
  - L1: + wall flux 1e5 W/m².
  - L2: + nuclear heating 1.6e6·exp(−d/0.15) W/m³ in the steel.
  - L1T: L1 in physical time.
- Partitions 40/80/128/256 all pass checks. 256 ranks ≈ 0.017 s/iteration.
- Research basis for the loads (EU-DEMO WCLL): first-wall flux about 0.5 MW/m² nominal; EUROFER nuclear heating 7–9 W/cm³ at the first wall. Laminar water can carry only about 0.1 MW/m², so the laminar study uses 1e5 W/m² and 298 K water. Realistic loads need turbulent flow (about 2 m/s, RANS/SA).

### Run history on SUNBIRD2

| Run | Settings | Result |
|---|---|---|
| 29482 | L1, global dt 4.96e-4 s, CP_SCALE 4, old energy code, 290k it | plateau RelU ~1.7e-3, RelW ~7e-3; Umax 0.30–0.33; mass 0.08 %; Q_out/Q_in 1.075 |
| 29510 | L1, LTS, new energy code | LTS field refreshed every ~2.5 iterations; exact 12,500-iteration limit cycle; cancelled |
| 29513 | L1, LTS + monotone + AFC (old AFC on solids), CP_SCALE 1, 256 ranks, 117.5k it | refreshes stopped at ~340; **flat plateau** RelU ≈ 3e-3, RelW ≈ 1.1e-2, RelT ≈ 7e-5; Umax 0.297 (jumps to 0.303); Tmax 407–413 K; mass balance 3e-4; **Q_out/Q_in = 1.041** at 100k; cancelled. VTUs at 50k and 100k kept |
| (submitted 1 Nov) | L0T: flow-only transient from the 29513 100k field (`legacy_vtu`, initial condition), averages from 10 s, end 50 s, 256 ranks | check its status; its last VTU carries `mean_velocity` |
| (suggested) | L0 at 0.025 kg/s (half flow), LTS, 128 ranks | not yet confirmed as run: the decisive physics test |

### Evidence for unsteady flow

`diffsnap.py` (appendix A) compared the 29513 snapshots at 50k and 100k iterations:
- 6 % of fluid nodes changed by more than 5 % of Umax.
- 22.5k of the 39.6k changed nodes lie within 40 mm of the heated wall (y = 240–280 mm).
- The largest changes are at x ≈ 113–123 mm, y ≈ 232–268 mm: the **outlet-side top bend of channel 1**. There |du| (up to 0.18 m/s) exceeds the local |u| (0.07–0.10), so a recirculation zone is moving. ΔT up to 60 K near the heated wall.

The pattern is smooth, localised and convected. That is the signature of physical separation, not of a numerical checkerboard. Not yet proven: the half-flow run and/or an OpenFOAM `pimpleFoam` cross-check would settle it.

### External review (another agent, on the old `main` ea7cd29): verified verdicts

- Correct: energy is one-way coupled, so the energy fixes cannot cure RelU/RelW. Equal fillets widen the bends. 5° leg tilt. **LTS projection inconsistency**: the pressure operator is weighted by element dt, but the Step-3 correction uses the nodal harmonic `elcoe2`. This is real in the code; its impact is unproven (the U-channel converges to 1e-13 with it). **Missing MPI outlet backflow**: real, fixed in 378c952.
- Wrong: "energy stabilisation sign + is correct" (it is anti-diffusive for P1). "Channel 4 has the highest velocity" (it has the lowest flow). Its bend-width formula √2·w − 2R(√2−1) has a sign error.
- Valid only for the old global-dt run: under-stabilisation from the global dt (CBS stabilisation about 40× below SUPG at core cells). LTS runs give each element its own dt.

---

## 8. Verification cases (re-run after every change)

| Case | Where | Reference / result |
|---|---|---|
| Conduction slab with nuclear heating | `cases/nuclear_heating_slab`, `tools/check_nuclear_slab.py`, `selfcheck_slab.slurm` | exact T(L) = 513.83 K; serial error 0.21 K (4.2e-4); MPI = serial to 2e-13 K; heat input 16.918 W (rel. err 8e-5). Explicit: 89,434 iterations. Implicit steady energy: 16 GMRES its, 3 ms, same field to 3e-4 K |
| CHT U-channel | `cases/cht_uchannel_check` (60×50×10 mm block, 5 mm U-channel, ṁ 1.25e-3 kg/s, flux 1e5 + nuclear) | explicit LTS + AFC: steady at ~13.6k iterations, Q_out/Q_in 0.9933, Tmin 298.15, Tmax 408.6 K |
| Lid-driven cavity | `cases/lid_driven_cavity`, `tools/plot_ldc_centerlines.py` | Ghia et al. (1982) centrelines |
| Flat plate SA | `cases/flat_plate_*` | NASA TMR |

U-channel implicit steady energy (frozen converged flow):
- **No AFC**: 63 GMRES its, 0.08 s, explicit-residual check 3.5e-12. T range 214.9–506.2 K: the unlimited scheme undershoots the inlet by 83 K.
- **AFC**: 167 defect corrections, 2.9 s, Tmin = 298.15, Tmax 389.7 K, Q_out/Q_in 0.9983.

The explicit-AFC and implicit-AFC peak temperatures differ by 19 K on this coarse mesh. The explicit limiter's bounds depend on the local dt, so its steady state is not unique. This means the near-wall thermal layer is under-resolved; see section 9.

---

## 9. Issues found and fixed (chronological, with the lesson)

1. `build_case` wiped the case dir and its own log. Fix: submit from the parent directory with `--output` in the parent.
2. The smoke-test awk matched comment lines. Fix: skip `#` lines.
3. `pkill -f` killed the agent's own shell. Use `pkill -x`.
4. Heat-input report invisible in MPI (stdout silenced). Fix: print with `fputs`.
5. Oversubscription / bad binding made runs about 1000× slower. Fix: `--hint=nomultithread --cpu-bind=cores --distribution=block:block`.
6. `energy_balance.py`: gap-based port grouping split ports, and absolute-T enthalpy let mass error swamp heat. Fix: union-find on shared face nodes; enthalpy relative to the mixed-mean inlet T.
7. Energy: anti-diffusive streamline sign plus non-conservative advection gave Q_out/Q_in 1.236 on the U-channel. Fix: diffusive sign plus conservative form gives 0.993.
8. LTS refreshed every ~2.5 iterations, so the residuals cycled with period 12,500. Fix: monotone refresh.
9. The MPI path lacked outlet backflow control, so serial and MPI behaved differently at outlets. Fix: 378c952 (36 % → 0.3 % serial/MPI difference on a backflow test).
10. AFC limited conduction in solids, giving 3.3 K error on the slab; found because `blanket.sh` exported AFC into the self-check. Fix: fluid only; self-check now sets production env explicitly.
11. SUNBIRD2 user pitfalls:
    - Typing `<JOBID>` literally.
    - `rm -rf` of a running job's folder.
    - `tail -1` with several files (use `tail -n 1`).
    - Forgetting to source the env or set `T` in a new login.
    - WSL + VPN breaks networking; use PuTTY without VPN.
12. Pushing from the Claude sandbox is impossible (403). Code moved as git bundles/patches pasted into PuTTY as gzip+base64 heredocs with a sha256 check.
13. Classical BoomerAMG breaks down (`DIVERGED_BREAKDOWN`) on the non-symmetric energy operator. Use GMRES + ASM/ILU(1) (default) or BoomerAMG with AIR restriction (`-energy_pc_hypre_boomeramg_restriction_type 1`).
14. The steady AFC defect correction cycles with relaxation 1.0. Use 0.5 (the default).

---

## 10. The goal now: an implicit, accelerated solver (PETSc)

### Why

Steps 1 and 4 are explicit:
- Convection limits Δt ≤ h/|u|.
- Diffusion limits Δt ≤ h²/(2α); in the steel this means ~(L/h)² iterations for heat to cross the plate.
- The global Δt is set by the smallest wall cell.

The PETSc/hypre AMG work accelerated only the pressure solve inside each iteration, not the iteration count. The flow is one-way coupled to T (constant properties, Ri = 0), which allows the flow and the energy to be split exactly.

### What exists: `SteadyEnergySolver` (template for everything else)

`src/solver/SteadyEnergySolver.cpp`, `include/cbs/solver/SteadyEnergySolver.hpp`, `EnergyAssembly::elementOperator`.

- **DOF map**: one PETSc row per global node; owned nodes contiguous (`MPI_Exscan`); ghosts get the owner's DOF via `broadcastOwnedToGhosts`.
- **Matrix**: element operator `M_e` (the exact linear part of the explicit residual), `MatSetValues(ADD_VALUES)` into MPIAIJ. Dirichlet rows: skip element rows and put 1 on the diagonal (ADD_VALUES; never mix with INSERT). Constant part = explicit residual at T = 0, so all production terms are reused.
- **Solver**: KSPGMRES (restart 200, rtol 1e-12) + PCASM/ILU(1); prefix `energy_`; `KSPSetFromOptions`.
- **Self-check**: re-assemble the original explicit residual at the solution; it must be ~1e-12. **Keep this pattern for every new implicit operator.**
- **AFC (bounded)**:
  - Low-order matrix = M + graph viscosity d_ab = max(0, −M_ab, −M_ba) on fluid elements.
  - Antidiffusive fluxes f_ab = d_ab (T_a − T_b), Zalesak limiter with Q_i = q_i (T_i^max/min − T_i), q_i = Σ max(0, l_ib).
  - Defect correction `A_L T^{k+1} = b + g(T^k)` with one preconditioner and relaxation 0.5.
- **Use on SUNBIRD2**: `bash $T/blanket.sh energy <run_dir> L1` (or `L2`). It picks `mean_velocity` when present and uses the same rank count as that run.

### Roadmap (choose with the owner; deadlines are tight)

1. **Now (hours)**: run steady energy on the 29513 field and on the L0T time-averaged flow (L1, L2). Check `HEAT_INPUT_TOTAL_W` against Q_out of the mean fields (`energy_balance.py --mean`).
2. **Implicit transient energy** (BDF1/BDF2): `(C/Δt + A) T^{n+1} = C/Δt T^n + b`, with the lumped C (`elcoe2p⁻¹`) already available. Removes the h²/(2α) limit in the steel. Same matrix code as above, plus a mass term.
3. **Implicit momentum** (Route A, "implicit CBS"):
   - Picard-linearised convection (frozen u^n) in Step 1, solved with GMRES + ASM/ILU or BoomerAMG per component (or block).
   - Keep Steps 2/3 (already PETSc).
   - Steady flow: pseudo-time with a growing CFL (SER), or PETSc `TSPSEUDO`.
   - Fix the LTS projection inconsistency (element dt in Step 2 vs nodal `elcoe2` in Step 3) while there.
4. **Fully coupled Newton–Krylov** (Route B):
   - Equal-order P1/P1 with PSPG/SUPG. The τ can be CBS-derived: Nithiarasu, Bevan & Murali (2013) show AC-CBS with dual time stepping is structurally PSPG.
   - PETSc `SNES` (Newton or Picard plus line search).
   - Preconditioner `PCFIELDSPLIT` (Schur): AMG on the velocity block; LSC or PCD Schur approximation for steady flow. SIMPLE-type (`selfp`) only for small Δt (PETSc developers' advice; Shadid & Tuminaro 2012).
   - `TS` BDF2 for large-Δt transients (needed for FNO training data anyway); `TSPSEUDO` for steady.
5. **Mesh**: the water thermal boundary layer on the heated walls is about 0.8–1.4 mm (Lévêque estimate) against 0.75 mm wall cells, so only 1–2 cells. Add graded near-wall layers (first cell about 0.1 mm), then a mesh study (h3/h4/h5) on wall T and Nu. The implicit energy solver makes each mesh cheap.
6. **Partitioning** (speed only, not convergence): METIS multi-constraint weights so fluid and solid elements balance per rank.
7. If the flow is physically unsteady: either run transient and time-average, use a lower Re for a steady benchmark, use constant-width bends (R_outer = R_inner + w; the owner decides the geometry), or move to RANS (SA is in the code) for realistic loads.

Do not claim methodological novelty without a literature review (owner's rule).

---

## 11. Working rules from the owner

- **Develop only through commits on GitHub** (`claude/blanket` → PR to `main` when stable). Check all branches for the latest state before starting.
- **C++ style**: classic, readable student-research C++ (no fancy modern idioms). Comments in ALL CAPS with banner blocks:
  ```cpp
  // ================================
  // WHAT THIS BLOCK DOES
  // ================================
  ```
  Shell/Slurm scripts use the same ALL-CAPS banner headers (`# FILE / PURPOSE / USAGE / OUTPUT`).
- Every change gets a verification against a case with a known answer (section 8) **before** it touches the blanket. Put the numbers in the commit message.
- Formal, humble tone; plain language plus simple diagrams when explaining. The owner meshes on the Windows tower; do not mesh for him unless asked.
- Never store or repeat MFA secrets or tokens.
- The owner reports fortnightly to the supervisor with Beamer slides; results should be presentable (tables with numbers, residual plots via `tools/plot_residuals_pub.py`).

---

## 12. How SUNBIRD2 was operated with the owner (copy this workflow)

The agent never has a shell on SUNBIRD2. **The owner types every command** in PuTTY and pastes the output back. Every instruction must therefore be:
- exact and copy-pasteable;
- complete, including `source ~/cbs3d_env.sh; T=…` whenever a new login is possible;
- free of placeholders the owner must edit (no `<JOBID>`; derive it with a command instead);
- explicit about what the output must show for success and what to send back.

### 12.1 Start of every session

```bash
source ~/cbs3d_env.sh; T=$CBS3D_CODE/jobs/sunbird/templates
cd $CBS3D_CODE && git branch --show-current && git log --oneline -3     # which code is on the cluster
squeue -u $USER                                                         # what is running
```

A prompt without `(cbs)` in front means the env is not loaded: `$CBS3D_CODE` is empty and `cd` goes to `~`.

### 12.2 Moving code to the cluster

1. **Preferred**: the agent commits and pushes to GitHub. On SUNBIRD2:
   ```bash
   cd $CBS3D_CODE && git fetch origin && git checkout claude/blanket && git pull --ff-only
   bash $T/build_cbs3d.sh && bash $T/selfcheck.sh
   ```
   `bash $T/blanket.sh update` without arguments does pull + build + self-check.
2. **When the agent cannot push** (this happened all along): a gzip+base64 heredoc with a checksum, pasted into PuTTY:
   ```bash
   # AGENT SIDE
   git format-patch --stdout BASE..HEAD > x.patch; sha256sum x.patch
   { echo "cd ~ && base64 -d > x.patch.gz <<'EOF'"; gzip -9c x.patch | base64 -w 76; echo EOF;
     echo "gunzip -f x.patch.gz && sha256sum x.patch"; } > paste_x.txt
   # OWNER SIDE: paste paste_x.txt, compare the printed sha256, then
   bash $T/blanket.sh update ~/x.patch          # git am --3way + build + self-check
   ```
   - Test the paste file locally (`HOME=tmpdir bash paste_x.txt`).
   - Test that the patch applies to the exact cluster base (`git am` in a clean clone).
   - Compare tree hashes (`git rev-parse HEAD^{tree}`).
   - Files for the owner can also go by WinSCP into `~`.
3. **Self-check result**: `selfcheck.sh` submits from the current directory, so the report is `cbs3d_selfcheck_<JOBID>.out` there and the work files are in `$CBS3D_SCRATCH/selfcheck/<JOBID>/`. Find the newest report with:
   ```bash
   F=$(ls -t $(find ~ $CBS3D_SCRATCH -maxdepth 3 -name 'cbs3d_selfcheck_*.out' 2>/dev/null) | head -1); cat $F
   ```
   The last line must be `SELF-CHECK PASS`.

### 12.3 Running

```bash
cd $CBS3D_SCRATCH/cases/blanket_slice
bash $T/free_cores.sh                                  # idle physical cores per node
bash $T/blanket.sh smoke L1                            # 300 iterations on 40 ranks (start-up check)
bash $T/blanket.sh run L1 --sets 256 --time 24:00:00   # production; --sets pins the partition set
bash $T/blanket.sh status                              # jobs + latest progress line of each run
```

- A run writes `run_blanket_<JOBID>/` with `solver.out`, `solver.err`, `run_record.txt` (settings, git commit, exe sha256), `output/blanket/*.pvtu`, `*_residuals.csv` and `restart/step_N/`. The Slurm log `cbs3d_mpi_<JOBID>.out` is in the case folder.
- Custom `.par` files: write them with `sed`/`awk` from the repo `.par`, show the changed lines with `grep -A1` before submitting, then call `submit_auto.sh $PWD blanket --sets N --time T --par $PWD/file.par --matprop $CBS3D_CODE/cases/blanket_slice/blanket.matprop`.
- Exports (AFC, LTS, average, restart) must be set **before** `submit_auto.sh`; sbatch inherits the environment. `unset` restart variables afterwards so a later run does not restart by accident.

Live output (the owner asked for this often):
```bash
J=$(squeue -u $USER -h -o %i | tail -1)
tail -f $CBS3D_SCRATCH/cases/blanket_slice/run_blanket_$J/solver.out | grep --line-buffered Iteration
```

Checks and stopping:
```bash
bash $T/blanket.sh check $CBS3D_SCRATCH/cases/blanket_slice/run_blanket_$J   # balances + residuals.png
sacct -j $J --format=JobID,State,Elapsed,MaxRSS
scancel $J                                                                   # only after the wanted VTU exists
```

- Use **absolute paths** for run folders: the runs live in `$CBS3D_SCRATCH/cases/blanket_slice/`, not in `$CBS3D_SCRATCH`.
- Restart (same rank count): `CBS3D_RESTART=1 CBS3D_RESTART_ROOT=<run>/output/blanket/restart/step_N`, plus `--sets <same N>`.
- Delete only finished or cancelled runs: `rm -rf run_blanket_<id>` after `squeue` confirms it is not running.

### 12.4 Costs to keep in mind

- 256 ranks (2 nodes): 0.017 s per iteration on the 3.2 M mesh. 100k iterations take about 30 minutes.
- Each VTU is about 0.3 GB (ASCII), so write them every 20k–50k iterations.
- Partition export takes minutes; the conversion job about 10 minutes.
- The owner pays in wall time and patience. **Never start a long run that cannot answer a specific question.** Cancel runs whose answer is already known.

---

## 13. Rules of engagement for the next agent (the owner's explicit request)

The owner has a hard deadline and limited credit. Work in a **disciplined, fast, verifiable** way.

1. **Understand before you change.** Read the files you will touch, and the call sites, before editing. State in one line what you will change and why.
2. **One concern per commit.** Small commits, each building and passing tests. The commit message states the problem, the fix and the **measured** check (numbers, not adjectives).
3. **Verify every numerical change on a case with a known answer before the blanket**: slab (exact), U-channel (balance 0.993), cavity (Ghia). Also compare serial and MPI. New implicit operators must reproduce the explicit residual at their solution (the `SteadyEnergySolver` check).
4. **Do not break the explicit solver.** New methods sit behind switches (env var or `.par` flag); the default behaviour stays bit-identical (check with a before/after run).
5. **Instructions to the owner must be exact.**
   - Full commands in order, with environment set-up included.
   - The expected output for success, and what to send back if it differs.
   - No placeholders; no "adjust as needed".
   - One block per step, short.
   - If something depends on an earlier result, say so and stop there.
6. **Be accurate and say what is unknown.** Do not present guesses as facts. When the evidence is incomplete, give the specific test that settles it and its cost. Check documentation (PETSc manual pages) instead of guessing option names.
7. **Be quick.** No long preambles or repeated summaries. Prefer the smallest change that answers the question. Do not spend cluster hours on runs that cannot change a decision.
8. **Respect the owner's style**: classic C++, ALL-CAPS banner comments, plain-language explanations with simple diagrams, formal and humble tone.
9. **Keep GitHub the single source of truth.** Push every verified commit to `claude/blanket`. Keep `HANDOVER.md` (this file) up to date at the end of every working session: what changed, what was run, what is open.

---

## Appendix A: diagnostic scripts (not in the repo; copies in `$CBS3D_SCRATCH` on SUNBIRD2)

`diffsnap.py`: where does the flow change between two snapshots?
```python
import os, sys
sys.path.insert(0, os.path.join(os.environ["CBS3D_CODE"], "tools"))
import numpy as np
import energy_balance as eb
A, B = sys.argv[1], sys.argv[2]
pa, _ = eb.read_result(A)
pb, _ = eb.read_result(B)
P = np.vstack([p[0] for p in pa])
du = np.concatenate([np.linalg.norm(b[2] - a[2], axis=1) for a, b in zip(pa, pb)])
dT = np.concatenate([np.abs(b[3] - a[3]) for a, b in zip(pa, pb)])
U = np.concatenate([np.linalg.norm(b[2], axis=1) for b in pb])
umax = U.max()
print("Umax %.4f m/s   max|du| %.4f m/s   max|dT| %.3f K" % (umax, du.max(), dT.max()))
for f in (0.01, 0.05, 0.10):
    print("nodes with |du| > %4.0f%% of Umax: %d of %d" % (100 * f, (du > f * umax).sum(), du.size))
idx = np.argsort(du)[::-1][:15]
print("largest |du| (mm):   x       y       z     |du|    |u|")
for i in idx:
    print("            %8.1f %7.1f %6.1f  %.4f  %.4f" % (P[i,0]*1e3, P[i,1]*1e3, P[i,2]*1e3, du[i], U[i]))
hot = du > 0.05 * umax
if hot.any():
    print("region with |du| > 5%% Umax: x %.1f..%.1f  y %.1f..%.1f mm" % (
        P[hot,0].min()*1e3, P[hot,0].max()*1e3, P[hot,1].min()*1e3, P[hot,1].max()*1e3))
    h, e = np.histogram(P[hot,1]*1e3, bins=14, range=(0, 280))
    print("count by depth y (mm):", " ".join("%d-%d:%d" % (e[k], e[k+1], h[k]) for k in range(14) if h[k]))
```

`backflow.py`: inflow through the outlets (blanket: inlets x < 0, outlets x > 0, ports on y = 0).
```python
import os, sys
sys.path.insert(0, os.path.join(os.environ["CBS3D_CODE"], "tools"))
import numpy as np
import energy_balance as eb
pieces, path = eb.read_result(sys.argv[1])
P = np.vstack([p[0] for p in pieces]); U = np.vstack([p[2] for p in pieces])
on = (np.abs(P[:, 1]) < 1e-7) & (np.linalg.norm(U, axis=1) > 0)
for name, side in (("OUTLETS (x>0)", P[:, 0] > 0), ("INLETS  (x<0)", P[:, 0] < 0)):
    sel = on & side
    un = -U[sel, 1]                      # OUTWARD NORMAL OF y = 0 IS -y
    if un.size == 0:
        continue
    print("%s: %d wetted nodes, %d with u.n<0 (into domain), %d with u.n>0; u.n range %.4f .. %.4f m/s"
          % (name, sel.sum(), (un < 0).sum(), (un > 0).sum(), un.min(), un.max()))
```

## Appendix B: references

- P. Nithiarasu, R.L.T. Bevan, K. Murali (2013). An artificial compressibility based fractional step method … similarity with a monolithic method. *Computational Mechanics*, doi:10.1007/s00466-012-0719-5. <https://cronfa.swan.ac.uk/Record/cronfa20247>
- J.N. Shadid, R.S. Tuminaro (2012). Stabilization and scalable block preconditioning for the Navier–Stokes equations. *J. Comput. Phys.* 231(2) 345–363.
- petsc-users thread "Scalable Solver for Incompressible Flow" (May 2023): <https://lists.mcs.anl.gov/pipermail/petsc-users/2023-May/048513.html>
- petsc-users thread on BoomerAMG AIR restriction (April 2023): <https://lists.mcs.anl.gov/pipermail/petsc-users/2023-April/048347.html>
- Zienkiewicz, Taylor & Nithiarasu, *The Finite Element Method for Fluid Dynamics* (CBS scheme).
- Ghia, Ghia & Shin (1982), lid-driven cavity benchmark.
