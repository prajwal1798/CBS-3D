# CBS3D handover (v2): code, mathematics, SUNBIRD2 workflow, runs and next steps

Written 5 November 2026, 05:00 UK time, by the Claude agent that worked with the owner from 29 October to 5 November 2026. It replaces the v1 handover of 1 November (commit `395507d`).

- **Owner**: P. Bharadwaj, EngD/PhD, Swansea University (Zienkiewicz Institute); supervisor Prof. P. Nithiarasu.
- **Repository**: <https://github.com/prajwal1798/CBS-3D>.

**Read sections 0–3 before doing anything.** Sections 4–14 are reference material. Whenever this document and the code disagree, the code wins; check it.

---

## 0. TL;DR: where we are, at 5 Nov 2026 05:00

1. **What the code is.** CBS3D++_SI is a C++ P1-tetrahedron finite-element solver for incompressible flow with conjugate heat transfer (water + EUROFER97 steel).
   - It uses the Characteristic-Based Split (CBS) scheme, semi-implicit: only the pressure is implicit (PETSc CG + hypre BoomerAMG). Momentum (Step 1) and energy (Step 4) are explicit.
   - It runs with MPI domain decomposition on SUNBIRD2.
2. **What the study is.** A laminar CHT model of a fusion breeder-blanket slice, a WCLL-like slice with four nested U-channels. Two meshes exist:
   - **h5**: 3.21 M tets.
   - **h4**: 4.89 M tets.
3. **What has been produced** (job IDs in §9.4):
   - A **converged h4 L2 steady result**, run **30102**. L2 means first-wall flux 1e5 W/m² plus exponential nuclear heating.
   - A **600 s transient with a plasma pulse**, run **31157**, restarted from 30102. It produced probe traces every 0.01 s, a regular-grid FNO dataset every 0.1 s, and a VTU every 5 s.
4. **Open problem 1: the steady and transient answers disagree.** At full power, the transient's equilibrium is **~40 K hotter at the first-wall centre** than the steady solution 30102 (corner +13 K). Energy is conserved in both:
   - steady: Q_out/Q_in = 1.000;
   - transient (550–600 s): ⟨Q_out⟩ = 2549 ± 18 W against Q_in − dE_s/dt = 2530 W.

   The leading hypothesis is the CBS stabilisation τ ∝ Δt: the steady solve used local time steps (LTS), the transient used one small global Δt. **Diagnostic job 31524** (steady energy with a global Δt, stage `L2G`) was submitted at 03:50 on 5 Nov; **its result has not been seen yet** (§9.7).
5. **Open problem 2: FNO data.** One pulse trajectory is not enough to train a generalisable FNO. The plan is an implicit transient energy solver on frozen flows, giving many cheap trajectories (§12).
6. **Git state** (§2):
   - GitHub `claude/blanket` = `09221c1`.
   - SUNBIRD2 = `39e2bdb` (6 commits ahead, not pushed).
   - The agent sandbox has 2 more commits (`09e83e1`, `c15288e`), delivered to the owner as `cbs3d-v3.bundle`. Whether he has applied it is unknown.
   - **This file is commit `HEAD` of the v4 bundle.**
   - The agent could never push (403). The owner pushes from SUNBIRD2.

---

## 1. The owner and how to work with him (read this)

### 1.1 Rules he set (non-negotiable)

- **Develop only through GitHub commits.** Check *all* branches for the latest state before starting.
  ```bash
  git fetch --all
  git branch -a
  git log --oneline -5 origin/<branch>
  ```
- **C++ style**: classic, readable research C++, with no fancy modern idioms. Comments go in ALL CAPS in banner blocks:
  ```cpp
  // ================================
  // WHAT THIS BLOCK DOES
  // ================================
  ```
  Python and shell files use the same banner style for headers (`# FILE / PURPOSE / USAGE / OUTPUT / STATUS`).
- **Tone**: formal, humble, plain language; simple diagrams when explaining. He reports fortnightly to his supervisor with Beamer slides, so results must be presentable.
- **No novelty claims** without a literature review.
- **Never run meshing for him** unless he asks ("Just dont run the mesh on your own"). He meshes on the Windows tower.
- **Never store, repeat or ask for passwords, MFA codes or tokens.** He once pasted his SCW password into the chat and was advised to change it. Do not refer to it again.
- **Instructions must be exact and copy-pasteable**:
  - No placeholders such as `<JOBID>`; derive values with commands.
  - Include `source ~/cbs3d_env.sh; T=$CBS3D_CODE/jobs/sunbird/templates` whenever a new login is possible.
  - Say what success looks like and what to send back.

  The previous agent broke the placeholder rule twice, and he was rightly annoyed.
- He wants **fast, disciplined, accurate** answers, not long preambles. He is often under deadline pressure and sometimes frustrated: stay calm, own mistakes, fix them.

### 1.2 How he interacts with the cluster

**The agent never has a shell on SUNBIRD2.** He types or pastes every command and pastes the output back, sometimes as screenshots or attached `.txt` files. Design every step so that one paste does one thing, with a clear expected output.

He uses two computers; details in §3.3.
- **Windows tower** (home/office). Works:
  - PuTTY (right-click pastes);
  - WinSCP or FileZilla for files;
  - ParaView **6.2.0** GUI, but its `pvpython.exe` is broken (§10.4).
- **Ubuntu 26.04 laptop** (HP Victus, US keyboard layout; do not change the layout). Works:
  - **`plink -t -ssh s.2337862@login.sunbird.swansea.ac.uk`** in the normal terminal (copy/paste works there);
  - **FileZilla** (SFTP, Site Manager, *Logon Type: Interactive*, limit simultaneous connections to 1);
  - ParaView 6.1.1 (Kitware binary).

  Does **not** work: OpenSSH `ssh` and `pscp` from the terminal (§3.3). PuTTY's own window works, but pasting into it from Linux is awkward.

---

## 2. Git state, exactly

### 2.1 Branches

| Where | Ref | Commit | Notes |
|---|---|---|---|
| GitHub | `main` | `5153f99` (21 Oct) | old consolidated master; `claude/blanket` contains all of it (0 commits behind, 33 ahead) |
| GitHub | `claude/blanket` | `09221c1` (2 Nov) | last push by the owner |
| GitHub | `feat/fully-implicit-be-20261102` | `5a4a112` | **other agent**: FullyImplicitFlowSolver (§13.3) |
| GitHub | `audit/*-20261101/02` (7 branches) | various | **other agent**: AFC and implicit-energy experiments, CBS literature audit; not merged; not reviewed in depth |
| SUNBIRD2 `$CBS3D_CODE` | `claude/blanket` | **`39e2bdb`** | `09221c1` + 6 commits below (applied 4 Nov from `cbs3d-pulse-v2.bundle` by fast-forward); built and self-checked (job 31156 PASS) |
| Agent sandbox | `claude/blanket-figs` | `c15288e` → + this handover | 2 more commits; delivered as `cbs3d-v3.bundle` (v4 adds this file) |

Commits on top of `09221c1` (oldest first):

| Commit | What | Verified by |
|---|---|---|
| `76debe8` | `tools/coherent_structures.py` (ParaView figures: Q, λ2, λci, helicity; temperature/velocity fields) | synthetic 2-piece duct in ParaView 6.1.1 (§10.3) |
| `fb18a92` | **HeatPulse** f(t) + **PointSampler** (probes + regular grid) | U-channel MPI tests (§7.5) |
| `4b09d7c` | stage `L2P`, `pulse_600s.txt`, `blanket_probes.txt`, `blanket.sh pulse/pulse-resume` | mock dry run; probe positions checked on a gmsh mesh |
| `39e2bdb` | `tools/heat_transfer_coefficient.py` (consistent wall flux, h, Nu) | U-channel: ∫q dA = 76.245 W vs Q_in 76.246 W |
| `09e83e1` | `blanket.sh update` fetches a bundle into a side ref and fast-forwards (git refuses to fetch into the checked-out branch) | bash -n; the manual equivalent worked on SUNBIRD2 |
| `c15288e` | `tools/probe_vtu.py` (probe values from saved VTUs) | reproduces the in-solver probes exactly (U-channel) |
| (v4) | this `HANDOVER.md` | |

### 2.2 Getting code onto SUNBIRD2

- The agent sandbox **cannot push** (403: repository not in the session's authorised set). If your environment can push, push `claude/blanket` and have the owner `git pull --ff-only`.
- Otherwise the agent makes a **git bundle**; the owner uploads it with FileZilla/WinSCP to `/scratch/scw7008/s.2337862` and runs:
  ```bash
  source ~/cbs3d_env.sh; T=$CBS3D_CODE/jobs/sunbird/templates
  cd $CBS3D_CODE && git fetch $CBS3D_SCRATCH/cbs3d-v4.bundle +refs/heads/claude/blanket-figs:refs/remotes/bundle/latest && git merge --ff-only bundle/latest && git log --oneline -3
  bash $T/build_cbs3d.sh && bash $T/selfcheck.sh      # only if C++ changed; self-check report: $CBS3D_CODE/cbs3d_selfcheck_<job>.out, last line SELF-CHECK PASS
  cd $CBS3D_CODE && git push origin claude/blanket     # owner's fine-grained PAT as password; credential.helper cache
  ```
  - Bundles are made with `git bundle create X.bundle origin/claude/blanket..claude/blanket-figs`. Their prerequisite is `09221c1`, which the cluster has.
  - After `09e83e1` is on the cluster, `bash $T/blanket.sh update $CBS3D_SCRATCH/X.bundle` does all of the above, including build and self-check.
- **Untracked file on the cluster**: `$CBS3D_CODE/cases/blanket_slice/blanket_L2G.par`, created with `sed` for job 31524. Either commit an equivalent file to the repo or delete it after the test.

---

## 3. Machines, access and environment

### 3.1 SUNBIRD2 (Supercomputing Wales)

- **Login**:
  - Host `login.sunbird.swansea.ac.uk` (round-robin `sb-login-01/02`; `sunbird.swansea.ac.uk` resolves to 137.44.57.130).
  - User `s.2337862`, Slurm account `scw7008`.
  - MFA: password, then verification code. The server asks both in **one keyboard-interactive round with two prompts**.
- **Home**: `/pixhome/s.2337862` (FileZilla shows it on connect). **Scratch**: `/scratch/scw7008/s.2337862`, GPFS, 15 TB quota.
- **Environment** (every new login; a prompt without `(cbs)` means it is not loaded):
  ```bash
  source ~/cbs3d_env.sh; T=$CBS3D_CODE/jobs/sunbird/templates
  ```
  - `CBS3D_CODE=/pixhome/s.2337862/CBS3D_SB2/code`
  - `CBS3D_SCRATCH=/scratch/scw7008/s.2337862`
  - executables `build/mpi-release/cbs3d_parallel` and `build/serial-petsc-release/cbs3d_serial`
  - OpenMPI 5.0.10 (GCC 13.5);
  - PETSc 3.26 at `/pixhome/s.2337862/software/petsc-sb2` (hypre 3.2.0, ParMETIS);
  - conda env `(cbs)` with numpy, scipy (assumed, needed by the HTC tool), matplotlib and gmsh.
- **Nodes**: 128 physical cores (Slurm shows 256 hyperthreads). MPI jobs use `--hint=nomultithread --cpu-bind=cores --distribution=block:block`; without these, runs were about 1000× slower.
- **Cost**: h4 on 256 ranks takes ≈ 0.021–0.024 s per iteration. A 600 s transient took 9 h 11 min (job 31157).
- **Disk**:
  - one h4 VTU snapshot ≈ 0.9 GB (ASCII, 256 pieces of 3.3–3.8 MB); one h5 snapshot ≈ 0.75 GB;
  - run 31157 output ≈ 120 × 0.9 GB VTU + ~71 GB FNO grid + checkpoints.

### 3.2 Windows tower

- About 20 cores; gmsh; WSL (do not use WSL with the VPN, which breaks networking).
- Local results:
  - `C:\PhD\CFD\runs\blanket\`: run 30102 step 159055 files, `blanket_distributed.pvd`, `blanket_30102.pvd` (one-step) and `Contour_Plots\`.
  - Other copies under `C:\PhD\CFD\runs\`.
- ParaView 6.2.0 at `C:\Program Files\ParaView 6.2.0`. The GUI works. `pvpython.exe` exits with **-1073741515 (0xC0000135, DLL not found)**; the VC++ 2015-2022 redistributable is already installed (setup said a newer version exists). **Workaround**: run scripts from the GUI's *View → Python Shell* (§10.4), or reinstall ParaView (6.1.1 .msi is the tested version).

### 3.3 Ubuntu laptop (HP Victus 15, Ubuntu 26.04)

- **Keyboard**: the physical keyboard is US/India, so the X11 layout must stay **`us`**. A previous agent wrongly switched it to `gb` and back; it is `us pc105` now.
- **OpenSSH `ssh` fails to log in to SUNBIRD2**: the server accepts the password and then closes the connection ("Connection closed by 137.44.57.130 port 22"), with or without a post-quantum KEX, and with only keyboard-interactive enabled. `plink`, PuTTY and FileZilla (all PuTTY-based) work. Probable cause: the server's combined two-prompt keyboard-interactive exchange. Do not spend more time on it:
  - log in with `plink -t -ssh s.2337862@login.sunbird.swansea.ac.uk` and press **Enter** at "Access granted. Press Return to begin session";
  - transfer files with FileZilla.
- `~/.ssh/config` has a `Host sunbird` block with ControlMaster settings, unused.
- **ParaView 6.1.1**: Kitware tarball, installed under `~/opt/ParaView-6.1.1-MPI-Linux-Python3.12-x86_64` with `~/.local/bin/paraview` linked (as instructed; verify). `pvpython` is in the same `bin`.
- **Local data**, `~/PhD/CFD/runs/`:
  - `run30082/`: step 159054 of the h4 flow run, 257 files, `blanket_30082.pvd`.
  - `blanket_step_00196836*`: **unidentified** run, pieces 2.3 MB. Possibly the earlier L2 heat-up transient from 30101 (159055 + 37781 steps), but this is not verified. Do not use it without identifying it.
  - `pulse_now/`: run 31157 snapshot at t ≈ 35 s plus probes up to t = 39 s.

---

## 4. Code architecture

```
include/cbs/... mirrors src/...
src/core         CBSStateSI (global state; cfg = RunConfig), Array1D/Array2D (1-based)
src/io           MeshIO (.plt/.bco/.material/.matprop/.par), Post, DistributedPost (VTU/PVTU/PVD),
                 RestartIO (native checkpoints + legacy VTU import), TimeAverage, PointSampler (NEW)
src/preprocess   Preprocess (geometry, masks, BC state, outlet normals), NuclearHeating
src/boundary     Boundary (velocity/pressure/temperature BCs, MPI outlet backflow), DistributedSymmetry
src/timestep     TimeStep (stability dt, LTS, elcoe2/elcoe2p = lumped inverse mass x dt)
src/assembly     MomentumAssembly (Step 1), PressureAssembly (2), VelocityCorrectionAssembly (3),
                 EnergyAssembly (4; elementOperator for the implicit solver), TemperatureAFC, HeatPulse (NEW)
src/linalg       ConjugateGradient (serial), PetscPersistentDistributedPressureSystem (MPI pressure)
src/parallel     HaloExchange (owner/ghost broadcast, reverse sums), Coloring
src/solver       Solver, Steps (serial); SolverDistributed* (MPI); SteadyEnergySolver; Convergence
src/turbulence   Spalart-Allmaras (+ wall models), not used in the blanket study
apps/parallel/main.cpp   usage: cbs3dpp_si <case> <partition_root>   (executable name on SUNBIRD2: cbs3d_parallel)
tools/           python (see section 10)
jobs/sunbird/templates/  blanket.sh, run_mpi.slurm, submit_auto.sh, selfcheck.sh, selfcheck_slab.slurm,
                 build_cbs3d.sh, build_case.slurm, prepare_partitions.sh, free_cores.sh, ...
cases/blanket_slice/     BlanketSlice_CHT.geo, blanket.matprop, blanket_L0/L1/L1T/L2/L2P.par,
                         pulse_600s.txt, blanket_probes.txt
docs/            LINUX_SUNBIRD_RUNBOOK.md, user-guide/case-files.md, MESH_FORMAT.md, Boundary Condition IDs.txt
```

**Traps:**
1. **The MPI production loop is compiled through macro wrappers.**
   - `SolverDistributedProductionLoopRestart.cpp` `#define`s `DistributedPost` to `RestartDistributedPost` and `TimeStep` to `RestartTimeStep`, then `#include`s `SolverDistributedProductionLoop.cpp`.
   - `SolverDistributedProductionLoopRestartAfc.cpp` wraps that again, and only that file is in CMake.
   - Edit `SolverDistributedProductionLoop.cpp`. Never name a new class `TimeStep…` or `DistributedPost…`.
   - `DistributedPost::writeSolution(s, case, k)` writes **global** iteration `iiter_total + k`.
2. **Indexing.** Arrays are 1-based.
   - `intma(a, ie)` (a = 1..4); `coord(d, node)`.
   - Element gradients: `dNkdx((ie-1)*ndim*nep + (dim-1)*nep + a)`.
   - `detJ = 6V`. `mat_elem(ie) == 0` means fluid.
3. **MPI data model.**
   - Elements are disjoint across ranks (local `nelem` = owned elements). Nodes on partition boundaries are duplicated as **ghosts**, with `vtkGhostType = 1` on non-owned copies in the VTU.
   - Assemble into local nodal arrays, then `HaloExchange::sumGhostContributionsToOwners`, update the owners, then `broadcastOwnedToGhosts`.
   - MPI equals serial to round-off on the slab.
4. `std::cout` is silenced during MPI set-up (ScopedStdoutSilence); use `std::fputs` for start-up messages there.
5. Every distributed VTU piece contains its ghost nodes. **Merge pieces (ParaView *Clean to Grid*) before computing gradients**; otherwise seams appear (§10.3).

---

## 5. Mathematics, as implemented

### 5.1 Governing equations (dimensional SI, constant properties, no buoyancy)

Fluid (water, material 0):

$$\nabla\cdot\mathbf u = 0,\qquad \rho\Big(\frac{\partial \mathbf u}{\partial t}+\mathbf u\cdot\nabla\mathbf u\Big) = -\nabla p + \mu\nabla^2\mathbf u$$

Energy in both materials (monolithic CHT on a conforming mesh), with Q the volumetric source:

$$\rho c_p\Big(\frac{\partial T}{\partial t} + \nabla\cdot(\mathbf u (T - T_{ref}))\Big) = \nabla\cdot(k\nabla T) + Q$$

- **u = 0 in the steel.** `T_ref` = inlet temperature (conservative advection form).
- **One-way coupling.** The flow does not depend on T (constant properties, Ri = 0). Hence:
  - **T is linear** in the sources, the boundary heat flux and T_in, for a given velocity history;
  - the flow can be computed once and reused for many thermal problems.

### 5.2 CBS (`cbs_scheme = 1`, semi-implicit; only the pressure is implicit)

| Step | Equation (schematic) | Treatment |
|---|---|---|
| 1 | $\Delta\mathbf u^* = -\Delta t[\mathbf u\cdot\nabla\mathbf u - \nu\nabla^2\mathbf u] + \frac{\Delta t^2}{2}\mathbf u\cdot\nabla(\mathbf u\cdot\nabla\mathbf u)$ | explicit, lumped mass |
| 2 | $\Delta t\,\nabla^2 p = \nabla\cdot \mathbf u^*$ | **implicit**: PETSc KSPCG + hypre BoomerAMG; the operator is weighted by the element Δt |
| 3 | $\mathbf u^{n+1} = \mathbf u^* - \Delta t\nabla p$ | explicit |
| 4 | $C_L (T^{n+1}-T^n)/\Delta t = R(T^n)$ | explicit, lumped mass `elcoe2p` |

**Energy residual** (Galerkin P1, with the characteristic stabilisation in its integrated-by-parts **diffusive** sign):

$$R_a = -\int N_a\nabla\cdot(\rho c_p\mathbf u (T-T_{ref})) - \int k\nabla N_a\cdot\nabla T - \frac{\Delta t_e}{2}\int \rho c_p(\mathbf u\cdot\nabla N_a)(\mathbf u\cdot\nabla T) + \int N_a Q + \int_{\Gamma_{532}} N_a q''$$

In code: `EnergyAssembly`, element operator `elementOperator`, with τ factor `0.5 * delte(ie) * rho_cp_e(ie) * V`.

**Time step:**
- Per element: $\Delta t_e = c_{safm}\min(h/|\mathbf u|,\ h^2/(2D))$, with h = minimum altitude = 3V/A_max. Fluid D = max(ν, α); solid: diffusion limit only.
- Global Δt = min over elements (`ilots = 0`; required for physical-time runs).
- **LTS** (`ilots = 1`, steady): each element uses its own Δt_e. Frozen-field refresh with margin 0.8. `CBS3D_LTS_MONOTONE=1` lets a refresh only lower Δt, which stopped a limit cycle.

**Stabilisation depends on Δt** (CBS ≈ GLS/SUPG with τ = Δt/2; Codina & Zienkiewicz 2002):
- With LTS, τ_e ≈ c·h/(2|u|), roughly SUPG-like.
- With a global Δt, τ_e = Δt_min/2 ≪ h_e/(2|u_e|) in the large core cells, so the scheme is nearly Galerkin there, with AFC providing the bounds.
- **Steady solutions therefore depend on how Δt is chosen.** This is the prime suspect for the 40 K discrepancy (§9.7).

**AFC** (`CBS3D_TEMPERATURE_AFC=1`; `blanket.sh` sets it):
- A Zalesak flux limiter on the explicit Step-4 update.
- **Fluid elements only** (since `6077df2`; `CBS3D_AFC_SOLID=1` restores the old behaviour, which gave a 3.3 K error on the slab).
- Without AFC the water under/overshoots by tens of K: the cell Péclet number is in the hundreds and Pr ≈ 6.

**Outlet backflow** (BC 520): the inward normal velocity is removed, using a rank-reconciled nodal outlet normal.

**Steady stop:** RelU/V/W < tol and RelT < tol after `steady_min_iterations`. These are per-iteration relative changes, so a small Δt fakes convergence. **Always check the energy balance instead.**

### 5.3 Implicit steady energy (`SteadyEnergySolver`, `CBS3D_STEADY_ENERGY=1`)

For a frozen velocity the explicit residual is affine: R(T) = M T + c + q. The solver finds the steady state directly:

$$-M T = c + q \ \text{(free nodes)},\qquad T = T_D \ \text{(Dirichlet nodes, identity rows)}$$

- **DOF map**: owned nodes contiguous (`MPI_Exscan`); MPIAIJ; `ADD_VALUES` only.
- **Solver**: KSPGMRES (restart 200, rtol 1e-12) + PCASM/ILU(1), options prefix `energy_`. Classical BoomerAMG breaks down (DIVERGED_BREAKDOWN) on this non-symmetric operator.
- **AFC version**:
  - low-order operator with graph viscosity $d_{ab} = \max(0,-M_{ab},-M_{ba})$ on fluid elements;
  - limited antidiffusive fluxes via defect correction $A_L T^{k+1} = b + g(T^k)$, relaxation 0.5 (1.0 cycles);
  - env `CBS3D_STEADY_ENERGY_TOL` (1e-6 K), `_MAXIT` (300), `_RELAX` (0.5).
- **Self-check**: the original explicit residual is re-assembled at the solution and must be ≈ 1e-12.
- **Δt dependence**: τ uses `delte` as set by the set-up `TimeStep::computeTimeStep` with the **par's ilots**. Runs 30039–30102 used `ilots = 1`, i.e. LTS τ.
- **Use**: `blanket.sh energy RUN_DIR L1|L2|L2G`. It takes the last `.pvtu` of RUN_DIR, `mean_velocity` if present, and the same rank count, and writes one VTU at iteration (restart + 1).

### 5.4 Plasma pulse (NEW, `HeatPulse`)

$$q''(\mathbf x,t) = f(t)\,q''(\mathbf x),\qquad Q(\mathbf x,t) = f(t)\,Q(\mathbf x)$$

- f is piecewise linear from the table in `CBS3D_HEAT_PULSE_FILE` (`t f` rows, `#` comments), constant outside the table range. Optional period `CBS3D_HEAT_PULSE_PERIOD`.
- t is the solver's `rtime`: 0 for a VTU-initial-condition start; restored from native checkpoints, so a resumed run stays in phase.
- Implemented as `s.heat_source_scale`, applied in:
  - `EnergyAssembly` (volumetric source and the serial BC 532 term);
  - the production loop (`heat_flux_load`);
  - `SteadyEnergySolver`.
- Updated every iteration after `rtime += dt`. Logged in the probe CSV (`heat_scale`) and the grid `index.csv`.

### 5.5 Nuclear heating

$$Q = A\,e^{-d/\lambda},\qquad d = \max(0,(\mathbf x_c - \mathbf x_0)\cdot\hat n)$$

- Evaluated at element centroids, in the material listed in the profile.
- Blanket: A = 1.6e6 W/m³, λ = 0.15 m, x₀ = (0, 0.28, 0) (first wall), n = (0, −1, 0), steel only.
- Start-up log: `HEAT INPUT … HEAT_INPUT_TOTAL_W = 2.637232e+03` (h4: nuclear 1714.065 W + BC 532 flux 923.167 W, flux area 9.2317e-3 m²).

### 5.6 Point sampling (NEW, `PointSampler`)

- For a point p in tet e with barycentric weights λ_a: $\phi(p) = \sum_a\lambda_a\phi_a$ for T, u, v, w, p.
- Vorticity is the element curl, constant in a P1 tet: $\omega = \nabla\times\mathbf u$ with $\partial u_k/\partial x_j = \sum_a (\partial N_a/\partial x_j)\,u_{k,a}$.
- **Owner**: the lowest rank whose owned element contains p (`MPI_Allreduce` MIN). Values are summed onto rank 0. On a material interface the first-found element is used (fluid or solid), so keep probes off interfaces.
- **Grid**: cell-centred, $x_i = x_0 + (i+½)(x_1-x_0)/n_x$. Output is float32 `[field][k][j][i]` with a material mask (−1 outside the mesh).

### 5.7 Wall heat flux, h and Nu (NEW, `tools/heat_transfer_coefficient.py`)

Interface = faces shared by a fluid tet and a steel tet; n points from the steel into the fluid.

**Consistent (residual) flux** at an interface node a, from the steel's own discrete equation:

$$F_a = \sum_{e\in\text{steel}}\Big[f\,Q_e\frac{V_e}{4} + f\,\frac{q'' A_{532}}{3}\Big|_a - k_s V_e\nabla N_a\cdot\nabla T\Big],\qquad q_w = F_a/A_a$$

- $A_a$ is the lumped interface area.
- At steady state $\sum F_a = Q_{in}$ exactly. U-channel: 76.245 W against 76.246 W. Element gradients give only 59.9 W on the steel side and 6.1 W on the coolant side, so do not use element gradients for the wall flux.
- **Storage caveat**: for transient snapshots F omits the steel's heat storage, so it is valid for (quasi-)steady states only.

From the flux:
- $h_{in} = q_w/(T_w - T_{in})$.
- $h_b = q_w/(T_w - T_b)$, where T_b is a local mixing-cup temperature: u·e- and volume-weighted coolant mean in a slab |(x − x_f)·e| < 0.15 D_h within 1.6 D_h of the face, in the same channel (connected fluid region).
- $Nu = hD_h/k_f$, with D_h = 15/10/10/5 mm by channel (ordered by decreasing max y).

### 5.8 Coherent-structure criteria (`tools/coherent_structures.py`)

With J = ∇u, S = (J + Jᵀ)/2 and W = (J − Jᵀ)/2:

| Criterion | Definition | Reference |
|---|---|---|
| Q | $Q = ½(\|W\|^2-\|S\|^2) > 0$ | Hunt, Wray & Moin 1988 |
| λ₂ | middle eigenvalue of $S^2+W^2$, < 0 in a vortex | Jeong & Hussain 1995 |
| Swirling strength λ_ci | \|Im\| of the complex eigenvalue pair of J | Zhou et al. 1999 |
| Normalised helicity | $H_n = \mathbf u\cdot\omega/(\|\mathbf u\|\|\omega\|)$; its sign separates the cells of a Dean-vortex pair | Levy, Degani & Seginer 1990 |

Iso-levels are a fraction f of the 99.5th percentile of the positive values in the fluid: Q = f·Q_ref, λ₂ = −f·l2_ref, λ_ci = √f·lci_ref.

### 5.9 Global energy budget and analytical outlet temperature

- **Coolant capacity rate**: ṁc_p = 0.050 × 4182 = 209.1 W/K.
- **Steady mixed outlet**: T_exit = T_in + Q_in/209.1:

  | Stage | Q_in | T_exit |
  |---|---|---|
  | L1 | 923.2 W | **302.565 K** |
  | L2 | 2637.2 W | **310.762 K** |

- **Transient budget**: $f(t)\,Q_{in} = \langle Q_{out}\rangle + dE_s/dt\ (+dE_f/dt \approx 0)$, with $E_s = \rho c_p\sum_e V_e\bar T_e$ (`tools/solid_energy.py`).
- **Per-snapshot outflow is noisy**: Q_out = Σ ṁ c_p (T_mix − T_ref) over the outlet ports varies by ±60 W because of the unsteady flow, so average it over many snapshots.

### 5.10 Thermal time scales

- Conduction estimate for the 10 mm first wall: L²/α, with α_EUROFER = 28.3/(7850·500) = 7.21e-6 m²/s, gives 14 s.
- **Measured** from the pulse dwell, where the excess over T_in decays by a factor of about 0.45–0.6 per 30 s:

  | Region | Time constant |
  |---|---|
  | First wall | ≈ 40–60 s |
  | Corner | ≈ 60 s |
  | Steel between channels 1 and 2 | ≈ 90–100 s |
  | Bulk steel (mean T still rising at 600 s) | several hundred s |

- **Flow**: unsteady, with probe speed fluctuations of 0.01–0.03 m/s at about 0.1–0.15 m/s; the dominant slow motion is ≈ 2 s in the channel-1 bends.

---

## 6. Case files and pipeline

- **Mesh**: the owner runs gmsh on the Windows tower → `.msh` → `tools/gmsh_to_cbs.py` → serial CBS case → `tools/export_rank_local_cbs.py --parts N` → rank-local files plus `.mpi` maps. On SUNBIRD2: `bash $T/blanket.sh setup $CBS3D_SCRATCH/meshes/slice_h4.msh`.
- **Folders**:
  - h5 case: `$CBS3D_SCRATCH/cases/blanket_slice`;
  - h4 case: `$CBS3D_SCRATCH/cases/blanket_slice_h4`. Select it with `BLANKET_DIR=$CBS3D_SCRATCH/cases/blanket_slice_h4` in front of `blanket.sh`.
- **Partition ladder**: `partitions_40/80/128/256`. Restarts need the **same rank count**; all h4 production uses 256.
- **`.par` format**: labelled lines, each followed by a value line; `#` lines are comments. Key lines:
  - `ntime transient_on dtfixed dtfix iwrite`;
  - `CBS3D timestep controls: cbs_scheme ilots htype step2_check rem_deltp beta_opt dtfix_end csafm2 epsilon1 deltr` → `1 1 1 0 0 1 0 2.0 1.0e-6 1.0e20` (L2; `1 0 …` in L1T/L2P/L2G);
  - `Output Controls … end_rtime`;
  - `alpha_sf k_ratio source_solid heat_flux_bc` → `1.0 1.0 0.0 1.0e5`;
  - mass-flow inlet 0.050 kg/s at 298.15 K;
  - `CBS3D output/monitor controls: residual_log_enabled residual_log_every console_log_every live vtu_output_enabled vtu_output_every_iterations vtu_output_every_sim_time write_boundary_debug_arrays steady_min_iterations` (L2P: `1 10 500 0 1 0 5.0 0 1000000`);
  - optional nuclear block (see L2).
- **`blanket.matprop`**: water 997 / 4182 / 0.606 / 8.9e-4; EUROFER97 7850 / 500 / 28.3.
- **BC IDs**: 511 mass-flow inlet + T; 520 pressure outlet; 530 no-slip adiabatic wall; 532 no-slip heat-flux wall (the first wall). Full list in `docs/Boundary Condition IDs.txt`.

---

## 7. Runtime switches and `blanket.sh`

### 7.1 Environment variables

| Variable | Meaning |
|---|---|
| `CBS3D_TEMPERATURE_AFC=1`, `CBS3D_AFC_SOLID=1` | energy limiter; solids too (not recommended) |
| `CBS3D_LTS_MONOTONE=1`, `CBS3D_LTS_REFRESH_MARGIN` | LTS refresh rules |
| `CBS3D_CHECKPOINT_EVERY`, `CBS3D_CHECKPOINT_ROOT` | native checkpoints → `output/<case>/restart/step_N/` |
| `CBS3D_RESTART=1`, `CBS3D_RESTART_FORMAT=native\|legacy_vtu`, `CBS3D_RESTART_ROOT`, `CBS3D_RESTART_ITERATION` | restart (same rank count). Native: ROOT = the `step_N` directory |
| `CBS3D_RESTART_AS_INITIAL_CONDITION=1` | a VTU starts a new run at rtime = 0 |
| `CBS3D_RESTART_VELOCITY_ARRAY=mean_velocity` | import the time-averaged velocity as `velocity` |
| `CBS3D_AVERAGE_START_TIME` / `_ITERATION` | running averages (mean/rms velocity, mean pressure, mean/rms/max T, u′T′ covariance; rank file version 2) |
| `CBS3D_STEADY_ENERGY=1` (+ `_TOL`, `_MAXIT`, `_RELAX`) | implicit steady energy instead of the loop |
| **`CBS3D_HEAT_PULSE_FILE`, `CBS3D_HEAT_PULSE_PERIOD`** | plasma pulse f(t) (NEW) |
| **`CBS3D_PROBE_FILE`, `CBS3D_PROBE_DT` (0.01 s), `CBS3D_PROBE_EVERY` (100, steady)** | probes → `output/<case>/<case>_probes.csv` (NEW) |
| **`CBS3D_GRID_DIMS="nx ny nz"`, `CBS3D_GRID_BOUNDS`, `CBS3D_GRID_DT` (0.1 s), `CBS3D_GRID_EVERY`** | regular grid → `output/<case>/fno_grid/` (NEW) |
| `CBS3D_PAR`, `CBS3D_MATPROP` | master files copied to every rank by `run_mpi.slurm` |
| `PETSC_OPTIONS` | e.g. `-energy_ksp_monitor` |

### 7.2 `blanket.sh` (`jobs/sunbird/templates/blanket.sh`)

| Command | What |
|---|---|
| `update [PATCH\|BUNDLE]` | git am / bundle fetch, then build, then self-check. Before `09e83e1` a bundle failed with "refusing to fetch into branch checked out"; use the manual fetch + `merge --ff-only` (§2.2) |
| `setup MSH` | convert the mesh, then partition 40/80/128/256 |
| `smoke L1` | start-up test |
| `run L0\|L1\|L1T\|L2 [--sets N --time T]` | production |
| `energy RUN_DIR L1\|L2\|L2G` | implicit steady energy on RUN_DIR's last flow field |
| **`pulse RUN_DIR`** | stage L2P from RUN_DIR's last VTU at t = 0. Defaults: pulse `pulse_600s.txt`, probes `blanket_probes.txt` every 0.01 s, grid `176 140 15` every 0.1 s, checkpoints every 50000, 24 h limit |
| **`pulse-resume RUN_DIR`** | continue from the newest `restart/step_*` |
| `check [RUN_DIR]` | mass/energy balance (`energy_balance.py`) plus residual plot |
| `status` | jobs plus latest progress |

`run_blanket_<JOBID>/` contains `solver.out`, `solver.err`, `run_record.txt`, `partitions/` and `output/blanket/` (VTUs, `blanket_distributed.pvd`, `*_residuals.csv`, `restart/`, and for pulse runs `blanket_probes.csv` and `fno_grid/`).

---

## 8. Verification cases (re-run after every change)

| Case | Result |
|---|---|
| Nuclear-heating conduction slab (`cases/nuclear_heating_slab`, `selfcheck.sh`) | exact T(L) = 513.83 K; serial error 0.21 K; MPI = serial to 2e-13; heat input 16.917 W (rel. err 8e-5). **Self-check job 31156 PASS** on `39e2bdb` |
| CHT U-channel (`cases/cht_uchannel_check`, 5 mm channel) | explicit LTS + AFC: steady ~13.6k iterations, Q_out/Q_in 0.993, Tmax 408.6 K. Implicit AFC: Tmax 389.7 K, Q_out/Q_in 0.998 |
| HeatPulse + PointSampler (U-channel, transient, 2 ranks, 0.3 s) | f = 0.5 gives exactly half the ΔT of f = 1 everywhere (1.5e-5 K on 16 K, float32); f = 0 leaves T = T_in exactly; velocity unchanged; probes equal VTU element interpolation to 1e-9; works with AFC on |
| HTC tool (U-channel steady) | consistent flux 76.245 W vs Q_in 76.246 W; off-interface steel residual 8e-4 W |
| Lid-driven cavity | Ghia et al. 1982 (2D). **Open**: 3D validation against Ku, Hirsh & Taylor 1987 (JCP 70:439–462) and Albensoeder & Kuhlmann 2005 (JCP 206:536–558; citation from memory, verify); needs the owner's data tables; `tools/plot_ldc_centerlines.py --ref-u/--ref-v` |

---

## 9. The blanket study: results and run registry

### 9.1 Geometry and loads (`cases/blanket_slice/BlanketSlice_CHT.geo`)

- **Coordinates**: x = 0 on the centre line, y = 0 on the back face, first wall at y = 0.28 m, z ∈ [0, 0.030] m (slice thickness 30 mm; channel mid-plane z = 0.015).
- **Plate**: EUROFER plate, back width 350 mm, sides at 85° (front half-width ≈ 150.5 mm), front fillet R6.
- **Channels**: four nested U-channels with square sections 15/10/10/5 mm, centred in z.
  - Front-run top edges at y = 0.270 / 0.220 / 0.150 / 0.100 m.
  - Legs start (outer edge on the back face) at drawing x = 20/65/105/135 mm from the left corner.
  - Bend fillets 10/5/5/3.5 mm. Equal outer and inner fillets make the bend diagonal ≈ √2·w, a built-in expansion.
  - **Inlets are the left legs (x < 0), outlets the right legs**; all ports are on y = 0.
- **Flow split** (measured): 0.0254 / 0.0110 / 0.0110 / 0.0026 kg/s; Re ≈ 1870 / 1250 / 1250 / 620.
- **Stages**:
  - L0: flow only;
  - L1: + BC 532 flux 1e5 W/m²;
  - L2: + nuclear heating;
  - L1T / L2P: physical time (L2P with the pulse);
  - L2G: L2 with a global Δt (diagnostic only).
- **Realistic DEMO loads** (≈ 0.5 MW/m² first-wall flux, 7–9 W/cm³ nuclear) need turbulent flow (≈ 2 m/s, RANS); this study is laminar on purpose.

### 9.2 Meshes

| Mesh | Tets (water + steel) | Nodes | BC faces 511 / 520 / 530 / 532 | Case dir |
|---|---|---|---|---|
| h5 | 3,212,185 (1,345,009 + 1,867,176) | 544,316 | 1144 / 1158 / 49316 / 2454 | `cases/blanket_slice` |
| h4 | 4,890,411 (2,186,830 + 2,703,581) | 823,234 | 1552 / 1576 / 54974 / 3772 | `cases/blanket_slice_h4` |

- Steel volume (h4): 2.452903e-3 m³.
- The water thermal boundary layer on heated walls is ≈ 0.8–1.4 mm against ≈ 0.75 mm wall cells, so it is under-resolved. **A third, near-wall refined mesh is needed** (first cell ≈ 0.1 mm, 8–10 layers on the heated channel walls). Ask the owner whether an "h3" exists and its tet count. He meshes; do not mesh for him.

### 9.3 Physics conclusion on steadiness

- The flow at Re ≈ 1870 is **unsteady**: separated flow behind the channel-1 bends, a slow ≈ 2 s motion, and plateaued residuals (RelU ≈ 2–3e-3).
- **Method adopted**: transient flow (global Δt), then time-average the velocity, then solve the implicit steady energy on the mean flow.
- RelU ~ 2e-3 in a transient is **not** a convergence criterion. Judge by statistical stationarity of the means (`meanconv.py`), energy closure and mesh independence.

### 9.4 Run registry (SUNBIRD2, `$CBS3D_SCRATCH/cases/…`)

| Job / run dir | Mesh | What | Key result |
|---|---|---|---|
| 29482, 29510, 29513 | h5 | early steady L1 attempts (old energy code; LTS limit cycle; LTS + monotone + AFC) | plateaus; 29513: Tmax 407–413 K, Q_out/Q_in 1.041; VTUs at 50k/100k |
| 29535, 29537, 29541, 29574, 29576 | h5 | earlier/other runs | **contents not verified**; avoid |
| **30029** | h5 | L0T flow, 50 s, averages 10–50 s | max⟨u⟩ 0.3058 → 0.3055 → 0.3054 (converging); pointwise change 5.4 %; 6 snapshots (4.5 GB) |
| 30037 / 30038 / 30039 | h5 | L1 steady energy on the 40 s / 31 s averages / 1500 AFC iterations | Tmax 405.31 / 405.13 / 405.31 K; Q_out/Q_in 1.00000; outlets ch1..4 306.21 / 299.38 / 298.26 / 298.19 K |
| 30040 | h5 | L2 steady energy | **Tmax 429.27 K**; Q_out/Q_in 1.00003; Q_in 2637.15 W; outlets 314.01 / 310.18 / 304.54 / 307.58 K. Files: step 180001 (256 pieces) |
| **30082** | h4 | L0T flow from rest, 70 s, averages 30–70 s (dt ≈ 4.4e-4, 0.021 s/step) | max⟨u⟩ 0.2959 (3.1 % below h5); change 4.5 %; mass imbalance 1.7e-4. **Final step 159054** (instantaneous `velocity` + `mean_velocity`, `rms_velocity`); 9 snapshots (9.3 GB) |
| 30083, 30084 | h4 | energy runs mistakenly started on step 0 of 30082 | cancelled; folders removed |
| 30101 | h4 | L1 steady energy on the 30082 mean flow | Tmax 411.78 K (+6.5 K vs h5); outlets 306.29 / 299.29 / 298.26 / 298.19 K |
| **30102** | h4 | **L2 steady energy** on the 30082 mean flow (`ilots = 1`, LTS τ) | **Tmax 437.06 K** (+7.8 K vs h5); Q balance exact; outlets 314.11 / 310.07 / 304.51 / 307.35 K. **Step 159055** = converged result (`velocity` = mean flow); step 159054 in the same folder = unsolved initial state, do not use |
| 31156 | — | self-check of `39e2bdb` | PASS (heat in rel. err 7.97e-5) |
| **31157** | h4 | **L2P pulse**, 600 s, from 30102 step 159055 at t = 0, AFC on, global Δt ≈ 4.34e-4 s | COMPLETED 4 Nov 22:28 after 9 h 10 min; 1,367,004 iterations (global iterations 159055 → **1526059**); stop reason "transient end time reached". 23 probes every 0.01 s; FNO grid 176×140×15 every 0.1 s (~6001 samples); VTU every 5 s (~121); checkpoints every 50k (~27). Final Tmax 450.26 K |
| **31524** | h4 | **L2G diagnostic**: steady energy on the 30082 mean flow with `ilots = 0` (global Δt τ); par `blanket_L2G.par` (untracked) | **pending**: compare with `probe_vtu.py` (§9.7) |

**Mesh independence, h5 → h4** (steady, LTS τ):
- Outlet temperatures agree to ≤ 0.2 K.
- Tmax changes by +6.5 K (L1) and +7.8 K (L2), about 6 % of the temperature rise, so **Tmax is not mesh-converged**.

### 9.5 Pulse run 31157: probe history (every 30 s)

```
   t[s]     f   FW_surf    FW_mid    corner    S_c1c2    C1_out
    0.0  1.00     378.6     361.1     411.8     321.9     313.8     <- 30102 steady field (initial condition)
   30.0  1.00     398.9     383.2     413.0     325.5     309.9     <- still heating at full power!
   60.0  0.00     363.7     360.1     384.3     325.2     309.1
   90.0  0.00     327.4     325.9     346.9     319.4     304.8
  120.0  0.00     311.9     311.1     327.3     313.5     302.1
  150.0  0.00     305.0     304.7     315.6     308.9     300.9
  180.0  0.00     301.9     301.7     308.8     305.6     299.7
  210.0  0.00     300.4     300.3     304.7     303.4     299.4
  240.0  0.00     299.6     299.5     302.3     301.9     298.9
  270.0  1.00     344.0     331.4     333.5     305.3     301.0
  300.0  1.00     384.0     369.1     372.4     313.5     305.6
  360.0  1.00     408.8     392.7     405.9     326.7     309.6
  420.0  1.00     414.0     397.5     417.4     333.4     310.5
  480.0  1.00     415.0     398.7     422.0     336.5     312.5
  540.0  1.00     416.5     399.5     424.6     338.0     311.6
  600.0  1.00     416.7     400.3     424.7     338.8     312.7
```

- Probe names are as in `blanket_probes.txt`. `C1_out` is the channel-1 outlet centre, which fluctuates by ±1–2 K.
- **During the first 0–39 s**:
  - the coolant outlets dropped 3–5 K within ≈ 4 s, then recovered slowly;
  - the flow regained its unsteadiness within ≈ 2 s;
  - speed after the first bend of channel 1 varies between 0.03 and 0.15 m/s.

  Plot: `pulse_probes_first_39s.png` (agent-made; regenerate from the CSV).

### 9.6 Energy budget of 31157 (550–600 s)

| Item | Value |
|---|---|
| Q_in (f = 1) | 2637.2 W |
| Steel storage dE_s/dt (`solid_energy.py`, 550 → 600 s: ΔE = 5362.7 J; mean T 323.04 → 323.60 K) | 107.3 W |
| Expected ⟨Q_out⟩ | 2529.9 W |
| Measured ⟨Q_out⟩, 11 snapshots 550–600 s (range 2469–2672 W, σ ≈ 59 W) | **2549 ± 18 W** |
| Closure | **0.7 % of Q_in** (within one standard error) |

At t = 600 s alone, Q_out/Q_in = 0.936. That is storage plus instantaneous fluctuation, not an error.

Port temperatures at 600 s:

| Port | ṁ [kg/s] | Mixed T [K] |
|---|---|---|
| 8o (ch1) | 0.02544 | 313.01 |
| 7o (ch2) | 0.01097 | 309.62 |
| 6o (ch3) | 0.01103 | 304.13 |
| 5o (ch4) | 0.00258 | 306.10 |

Mass imbalance 4.7e-4.

### 9.7 OPEN PROBLEM 1: why is the transient ~40 K hotter at the first wall?

Full-power equilibrium, transient (31157 at 600 s) against steady 30102:

| Probe | Transient | Steady | Difference |
|---|---|---|---|
| FW_mid_C | 400.3 | 361.1 | **+39 K** |
| FW_surf_C | 416.7 | 378.6 | **+38 K** |
| FW_corner_R | 424.7 | 411.8 | +13 K |
| S_c1c2_C | 338.8 | 321.9 | +17 K, still rising |
| Global Tmax | 450.3 | 437.1 | |

Both conserve energy, so the coolant carries the heat, but the transient discretisation transfers it from the wall to the coolant less effectively.

**Hypotheses:**
1. **(Leading) Stabilisation τ ∝ Δt.** 30102's steady solve used LTS element Δt (large τ in big coolant cells, more numerical transport); the transient used the global Δt_min (small τ, nearly Galerkin + AFC). CBS steady states depend on Δt (§5.2).
2. **Mean flow against instantaneous flow.** 30102 used the 30–70 s *time-averaged* velocity; the transient uses the instantaneous unsteady flow. The mean of T under a fluctuating u is not T under the mean u (u′T′ transport). This would be physics, not an error.

**Diagnostic: job 31524** (stage L2G = L2 with `ilots 0`, same mean flow as 30102). When it is finished, compare at the probes:
```bash
source ~/cbs3d_env.sh; R=$CBS3D_SCRATCH/cases/blanket_slice_h4/run_blanket_31157
python3 $CBS3D_CODE/tools/probe_vtu.py $CBS3D_SCRATCH/cases/blanket_slice_h4/run_blanket_30102/output/blanket/blanket_step_00159055.pvtu $CBS3D_SCRATCH/cases/blanket_slice_h4/run_blanket_31524/output/blanket/blanket_step_00159055.pvtu $R/output/blanket/blanket_step_01526059.pvtu
```
(`probe_vtu.py` needs `c15288e` on the cluster: bundle v3 or v4.)

- **31524 ≈ 31157** (FW_mid ≈ 400 K): τ is the cause. Then:
  - (a) make the steady and transient τ definitions consistent, for example τ from h/|u| and h²/α, independent of the run's Δt;
  - (b) quantify the τ-sensitivity as a discretisation error;
  - (c) the refined third mesh should shrink it. **Do not quote first-wall temperatures to better than ±20 K until this is resolved.**
- **31524 ≈ 30102** (≈ 361 K): the flow unsteadiness (u′T′) is the cause, which is a physical result. Confirm it with the u′T′ covariance accumulated by `TimeAverage` (version-2 rank files) in a transient window.
- Also check that 31524 printed a sane `t_max`, `converged`, and `residual_check` ≈ 1e-12.

**Sign argument (added 2026-11-06).** Unsteady transport ⟨u′T′⟩ is, on the whole, down-gradient: it acts as extra diffusion and *improves* wall-to-coolant heat transfer. A steady solve on the mean flow omits it and should therefore be *hotter* at the wall than the transient, yet 30102 is ~40 K *cooler*. Hypothesis 2 has the wrong sign to explain the gap. The trend with numerical diffusion has the right sign: Tmax 429.3 K (h5, LTS τ) < 437.1 K (h4, LTS τ) < 450.3 K (h4, transient, small τ); every step that removes streamline diffusion makes the wall hotter. The CBS energy stabilisation in `EnergyAssembly` is −(Δt_e/2)∫ρc_p(u·∇N_a)(u·∇T): advective part only, not the full residual, so it does not vanish at steady state. With LTS, Δt_e ≈ c·h/|u|, which gives a streamwise artificial conductivity ≈ c·Pe_h/2 times k_f (cell Péclet in the hundreds). Expect 31524 ≈ 31157.

### 9.8 Model-validity audit (2026-11-06): what the paper must state or fix

Physics (the model is the forced-convection, constant-property, single-phase limit):
- **Buoyancy is neglected but not negligible at these loads.** Wall-to-bulk ΔT ≈ 90 K at U ≈ 0.11 m/s gives Gr = gβΔT D_h³/ν² ≈ 1–2e6, Ra ≈ 1e7 (laminar) and Ri = Gr/Re² ≈ 0.3–0.6 in channel 1 (0.2–0.4 in channels 2–3, 0.1–0.2 in channel 4; β = 2.6–6e-4 1/K). That is the mixed-convection range (Ri ≳ 0.1; Metais & Eckert 1964; Morcos & Bergles 1975 for laminar heated tubes). Ri < 0.1 in channel 1 needs wall-to-bulk ΔT ≲ 15–30 K, i.e. loads ≈ 4× lower. (An earlier version of this note said Gr ≈ 1e9 and Ri ≈ 10²–10³: an arithmetic error by a factor of 1000.) DEMO flows are turbulent at ~2 m/s, so Ri ≪ 1 there; that justification does not carry over to this laminar case.
- **Adding buoyancy (assessed 2026-11-06)**: Ra/Ri are read from the .par but no body force is assembled (`docs/CBS3Dpp_SI_Full_Solver_Audit.md`). A Boussinesq term −β(T − T_ref)g, lagged explicitly in Step 1 (`MomentumAssembly::assembleStep1Rhs`, fluid elements), with a modified pressure p − ρ₀g·x, is ~1 day of code. Validate it on the differentially heated cubic cavity (BC 501/502; Fusegi et al. 1991, Tric et al. 2000). It ends one-way coupling: `CBS3D_STEADY_ENERGY` (frozen flow) must refuse to run with gravity on, and the linear-in-load FNO ensemble no longer holds.
- **Constant properties**: μ of water falls ≈ 4× from 298 K to ≈ 410 K (wall); the Sieder–Tate factor (μ_b/μ_w)^0.14 ≈ 1.2.
- **Boiling**: coolant-wall temperatures of ≈ 390–410 K exceed T_sat at 1 bar. State a system pressure (≥ ~4 bar) or no single-phase claim holds.
- **Defence and fix**: T is linear in the loads for a given flow, so every result is exactly the small-load forced-convection limit scaled up. Report θ = (T − T_in)·k/(q″·L) and Nu, not absolute kelvin. Otherwise lower the loads (ΔT_wall ≲ 10 K), or add Boussinesq buoyancy + μ(T) (this breaks one-way coupling and the linearity used for the FNO ensemble).

Numerics:
- Δt-dependent, non-residual energy stabilisation (§9.7). The steady result carries an O(h) streamline-diffusion error.
- **CBS equal-order (P1/P1) pressure stabilisation ∝ Δt** (Codina 2001). With the global Δt ≈ 4.4e-4 s against h/|u| ≈ 0.02 s in core cells, it is ~40× weaker than an h/|u|-scaled value. The momentum Step-1 stabilisation is equally small (cell Re ≈ 250). Check mid-plane p and element div(u) for node-to-node oscillation. Show that the unsteadiness (St, amplitude) is independent of mesh and Δt before calling it physical.
- **Mesh**: h5 → h4 moved Tmax by +7.8 K (6 % of the rise above inlet); the thermal boundary layer has 1–2 cells. Wall T, h and Nu are not mesh-converged until h3 (GCI, Celik et al. 2008).
- **BCs**: the uniform inlet velocity in every port imposes the flow split (area-weighted), so it is not computed. Outlet backflow is clipped (BC 520); measure the fraction of outlet faces that would back-flow.
- **Pulse start**: 31157 began from the steady state, not the transient equilibrium. The +13 K Tmax drift is not pulse physics (see §15, item 6).

What is sound: monolithic conforming CHT; conservative advective form (energy closes 0.7 % over 550–600 s and exactly when steady); AFC keeps T bounded; one-way coupling is consistent with the stated assumptions; physical-time stepping is explicit with CFL ≪ 1; unsteadiness at Re ≈ 1870 is plausible (Pothérat & Zhang, sharp 180° bend: unsteady from Re ≈ 800).

---

## 10. Tools, post-processing and visualisation

### 10.1 Python tools (`tools/`; numpy-based; run on SUNBIRD2 inside `(cbs)` or locally)

| Tool | Use |
|---|---|
| `energy_balance.py RESULT --plane y=0 [--log solver.out] [--mean] [--storage W]` | ports by union-find, ṁ and mixed T per port, Q_out, mass balance. `blanket.sh check` wraps it |
| `solid_energy.py R1 R2 …` | steel energy E_s per snapshot (ΔE/Δt = storage) |
| `heat_transfer_coefficient.py RESULT [--log solver.out] [--par …] [--heat-scale f]` | consistent q_w, h_in, h_b, Nu; writes `*_htc_interface.vtu`, `*_htc_channels.csv`, `*_htc_profiles.png/.csv`. **Not yet run on 30102**; it must print `consistent / input ≈ 1.000` |
| `probe_vtu.py R1 [R2 …] [--probes …]` | T at the probe points from saved results |
| `coherent_structures.py INPUT.pvd --what structures\|fields\|all` (pvpython) | §10.3 |
| `plot_residuals_pub.py`, `plot_ldc_centerlines.py`, `check_nuclear_slab.py`, `gmsh_to_cbs.py`, `export_rank_local_cbs.py` | as named |

Scratch helpers in `$CBS3D_SCRATCH` (not in the repo):
- `meanconv.py`: compares `mean_velocity` between snapshots;
- `diffsnap.py`, `backflow.py`: v1 handover, appendix A.

### 10.2 Getting data to the owner's computer

- **Single snapshot**: link the `.pvtu` plus its 256 `_rank_` pieces into a folder, write a **one-entry `.pvd`**, and `tar -ch` it into `$CBS3D_SCRATCH/transfer/<name>.tar`. The owner downloads it with FileZilla/WinSCP and runs `tar -xf`. The block that does this (newest snapshot of 31157) is in Appendix A.
- **One folder per run, always.** Different runs share file names (`blanket_step_00100000_rank_0000.vtu`), and mixing them silently corrupts data. ParaView's "File for piece N cannot be read" came from a `.pvd` listing missing steps; fix it by writing a one-step `.pvd`.
- **Which file is which**: in an energy-run folder the **last** step is the solution and the first is the unsolved start. Iteration numbers differ between meshes (h5 50 s ≈ 180,000 iterations; h4 70 s ≈ 159,054), so they are not comparable.

### 10.3 `coherent_structures.py` (ParaView 6.1 pvpython)

```bash
pvpython --force-offscreen-rendering coherent_structures.py ~/PhD/CFD/runs/run30082/blanket_30082.pvd
pvpython --force-offscreen-rendering coherent_structures.py <run30102 one-step .pvd> --what all --velocity-label "time-averaged velocity" --sweep
```

- **Pipeline**: PassArrays → **CleantoGrid (merges partition pieces)** → Threshold `material_id` = 0 → Gradient (velocity; vorticity, Q) → ProgrammableFilter (numpy λ₂, λ_ci, H_n, with a check of the gradient storage convention) → contours and slices.
- **Output**: `figures/*.png` (white background, Viridis/Inferno/Cool-to-Warm, MathText labels, parallel projection for top views), `coherent_structures.pvsm` and `summary.json`.
- **Tested** on a synthetic two-piece duct in ParaView 6.1.1:
  - Q within 6 % of the hand value;
  - VTK Q equals numpy Q to 1e-15;
  - without merging, Q at a partition seam is wrong by 10 % of Q_ref.
- **Not yet run on real data.**
- **ParaView 6 differences**: presets are named `Viridis`/`Inferno`; the scalar-bar format uses std::format (`{:.4g}`); MathText cannot contain `|` (use `\vert`).

### 10.4 ParaView on the owner's machines

- **Windows**: `pvpython.exe` is broken (§3.2). In the GUI, open *View → Python Shell* and paste one line:
  ```python
  import sys, os; sys.argv = ["cs", r"C:\PhD\CFD\runs\blanket\blanket_30102.pvd", "--what", "all", "--velocity-label", "time-averaged velocity", "--sweep"]; exec(open(os.path.join(os.path.expanduser("~"), "Downloads", "coherent_structures.py")).read(), {"__name__": "__main__"})
  ```
  This has not been confirmed to work by the owner yet.
- **GUI pipeline**:
  1. Clean to Grid.
  2. Threshold `material_id` 0 (water) / 1 (steel).
  3. Gradient on `velocity` (advanced: Compute Vorticity, Compute QCriterion).
  4. Contour on QCriterion at ~1–2 % of its max (real data: Q range ≈ −1.6e4..1.2e4 s⁻², level 100–300).
  5. Slice normal to the channel after a bend, coloured by the streamwise vorticity (diverging map) for Dean vortices.
- **Fixed colour ranges** when comparing times: 298–450 K.
- **Touchpad navigation**: two-finger scroll zooms. Set *Edit → Settings → Camera → Left Button + Shift = Pan* (gear icon for advanced options).
- **Owner's earlier `Contour_Plots`** (made in the GUI from the `C:\PhD\CFD\runs\blanket` folder): physically plausible (pressure drop, bend stagnation, Dean-type tubes behind bends, hot first-wall edges), but not publication-ready:
  - the white "vein" network is partition faces (no Clean to Grid), and many small Q blobs are seam artefacts;
  - jet colour map on black; perspective in the "top" views (trapezoid);
  - missing units; pressure units unverified (values ≈ 0.07 suggest kinematic p/ρ; check the code before labelling).

---

## 11. Data locations (cluster)

```
$CBS3D_SCRATCH/cases/blanket_slice/run_blanket_{30029,30037,30038,30039,30040,...}/output/blanket/   h5
$CBS3D_SCRATCH/cases/blanket_slice_h4/run_blanket_30082/output/blanket/   h4 flow (step 159054 = final, with averages)
$CBS3D_SCRATCH/cases/blanket_slice_h4/run_blanket_30101/output/blanket/   h4 L1 steady (step 159055)
$CBS3D_SCRATCH/cases/blanket_slice_h4/run_blanket_30102/output/blanket/   h4 L2 steady (step 159055)
$CBS3D_SCRATCH/cases/blanket_slice_h4/run_blanket_31157/output/blanket/   pulse run:
      blanket_step_*.pvtu (+_rank_*.vtu), blanket_distributed.pvd (times), blanket_probes.csv (~70 MB),
      fno_grid/{grid.json, mask.i8, index.csv, sample_<iteration>.f32}, restart/step_*/, blanket_residuals.csv
$CBS3D_SCRATCH/cases/blanket_slice_h4/run_blanket_31524/   L2G diagnostic
$CBS3D_SCRATCH/transfer/                                   tarballs for download
$CBS3D_SCRATCH/{meanconv.py, diffsnap.py, ...}             helpers
```

**FNO grid layout** (31157):
- dims (nx, ny, nz) = (176, 140, 15) over the mesh box x ∈ [−0.175, 0.175], y ∈ [0, 0.28], z ∈ [0, 0.03] (≈ 2 mm spacing; verify the box in `grid.json`).
- 369,600 points × 8 fields × float32 = 11.8 MB per sample.
- Fields `T, u, v, w, p, omega_x, omega_y, omega_z`. Mask int8: −1 outside the mesh, 0 water, 1 steel.

Loader:
```python
import numpy as np, json, csv
d = "run_blanket_31157/output/blanket/fno_grid/"
g = json.load(open(d + "grid.json")); nx, ny, nz = g["dims"]
rows = list(csv.DictReader(open(d + "index.csv")))            # file, iteration, time, heat_scale
mask = np.fromfile(d + "mask.i8", dtype="i1").reshape(nz, ny, nx)
x = np.fromfile(d + rows[k]["file"], dtype="<f4").reshape(8, nz, ny, nx)   # x[0] = T
```

`blanket_probes.csv` columns: `iteration,time,heat_scale,<probe>_{T,u,v,w,p,wx,wy,wz}` for 23 probes:
- steel: `FW_surf_C, FW_mid_C, FW_mid_L, FW_mid_R, FW_corner_L, FW_corner_R, FW_mid_C_zface, S_c1c2_C, S_back_C`;
- coolant: `C1_in, C1_leg_L, C1_run_L, C1_run_C, C1_run_C_wall` (1 mm below the heated wall), `C1_run_R, C1_leg_R, C1_out, C2_run_C, C2_out, C3_run_C, C3_out, C4_run_C, C4_out`.

All points were checked against a gmsh mesh of the `.geo` for the right material; coolant points are centred in their channel to ±0.1 mm.

**Housekeeping**: after the analysis, keep only the newest `restart/step_*` of 31157 (≈ 27 checkpoints exist). Ask the owner before deleting anything.

---

## 12. FNO / neural-operator data: assessment and plan

**Is run 31157 enough? No.**
- It is one trajectory of one scenario. Its 6001 grid samples are consecutive and strongly correlated.
- An FNO trained on it learns to step that run forward; it will not generalise to other pulse shapes, power levels or flow rates.
- FNO papers (Li et al., ICLR 2021) train on hundreds to about a thousand trajectories.

**The physics is linear in the heating.** With constant properties and no buoyancy, T responds linearly to f(t), q'' and Q for a given flow. A linear reduced-order model (POD/DMD with inputs, or a transfer function) could reproduce the pulse response almost exactly, so reviewers will ask why an FNO is needed. Nonlinearity must enter through:
- **flow rate**, the strongest and available now (the flow field changes nonlinearly with Re);
- temperature-dependent μ(T) of water and k(T) of EUROFER (not implemented);
- geometry variation.

**Plan:**
1. Keep 31157 as the **high-fidelity test trajectory** (full CBS, unsteady flow). Use it to quantify the frozen-flow error.
2. Implement an **implicit transient energy solver** on a frozen velocity, extending `SteadyEnergySolver`:
   $$\Big(\frac{C_L}{\Delta t} - M\Big) T^{n+1} = \frac{C_L}{\Delta t}T^n + f(t^{n+1})\,(c_Q + q) + c_0$$
   - BDF1 first, then BDF2. Same matrix for constant Δt: set up the preconditioner once and warm-start.
   - Δt = 0.1–0.5 s instead of 4.4e-4 s, so a 600 s trajectory costs minutes.
   - **Verify** against 31157 with the same pulse, using the 30082 mean flow, and against the slab analytic transient.
   - **Use the same τ definition as the transient** (open problem 1).
3. **Ensemble**: 3–5 flow rates (one 70 s flow run each, ≈ 8 h at 256 ranks) × 50–200 randomised power histories (ramps, dwells, multiple pulses, levels). Store on the same regular grid, plus probes.
4. **FNO task**: input = (T field at t, future f(t…t+H), flow-rate parameter or mean-flow fields), output = T at t+Δt…t+H. Test on held-out pulse shapes, flow rates, and the full CBS trajectory.
5. **Grid caveat**: 2 mm gives 2–5 points across the first wall and channel 4. That is adequate for steel temperature (the design quantity), coarse for coolant velocity.

---

## 13. Other agents' work (reviewed)

1. **Commits by another agent on `claude/blanket`** (`f2ffbcc`…`09221c1`): self-check clears the inherited production environment; executable-age check; **u′T′ covariance** in `TimeAverage` (rank file version 2); `tools/solid_energy.py`; `energy_balance.py --storage`. All good.
2. **`audit/*` branches** (1–2 Nov): AFC direction/globalised-defect/pseudotime/residual-merit experiments, implicit-energy AFC residual, CBS literature verification. Not merged; read them before redoing similar work.
3. **`feat/fully-implicit-be-20261102`** (FullyImplicitFlowSolver: u, u*, p unknowns; FGMRES + fieldsplit Schur a11; Picard relaxation 0.7). Review findings:
   - steady τ = global min Δt (under-stabilised);
   - backward-Euler τ = physical Δt (over-stabilises at large Δt);
   - not verified against a known solution;
   - its implicit high-order AFC failed. Its low-order operator gives Tmax 375.8 K on the U-channel, against 389.7 K (`SteadyEnergySolver` AFC) and 408.6 K (explicit);
   - it extended graph viscosity to solids, which is wrong for conduction (cf. the 3.3 K slab error).

   Do not merge it without verification.

---

## 14. Open tasks, prioritised

1. **Read the result of job 31524** and run the `probe_vtu.py` comparison (§9.7). Decide on the stabilisation fix. *(Hours.)*
2. **Make the steady and transient stabilisation consistent**, behind a switch, with the default bit-identical. Verify on the U-channel and the slab, then recompute L2 on h4.
3. **Third mesh**, near-wall refined (owner meshes), then the mesh study on wall T, Tmax and Nu (h5/h4/h3), with paper-style plots:
   - T through the first wall at Tmax and along the channel-1 heated wall;
   - velocity profiles across channel 1 before and after the bend;
   - outlet T bars;
   - Tmax against h.
4. **Run `heat_transfer_coefficient.py` on 30102** (and on 31524 or the new consistent result). Check `consistent / input ≈ 1.000`, then the profiles plot.
5. **Full pulse analysis**: download `blanket_probes.csv` of 31157 and plot all 600 s (power, steel probes, outlets, flow). Extract time constants. Make VTU figures at t = 30, 120, 240, 330, 600 s with a fixed 298–450 K range.
6. **Coherent-structure figures** on 30082 step 159054 (instantaneous and mean) with `coherent_structures.py`. On Linux: `~/opt/ParaView-6.1.1-MPI-Linux-Python3.12-x86_64/bin/pvpython --force-offscreen-rendering …`.
7. **Implicit transient energy solver** plus the FNO ensemble (§12).
8. **3D lid-driven-cavity validation** (needs the owner's reference tables).
9. **Housekeeping**: push `claude/blanket` (owner, from SUNBIRD2); commit or remove `blanket_L2G.par`; prune checkpoints; keep this file current.

---

## 15. Issues log (new since v1; v1 items still apply)

1. `blanket.sh update BUNDLE` failed when `claude/blanket` was checked out ("refusing to fetch into branch"). Fixed in `09e83e1`; manual fetch + `merge --ff-only` meanwhile.
2. The OpenSSH login failure on the Ubuntu laptop (§3.3) cost hours. Use plink/FileZilla.
   - Keyboard layout was a red herring: the laptop is US.
   - Eduroam vs visitor Wi-Fi made no difference (the visitor network has an outside IP and still reached Sunbird).
3. The steady energy run of 30102 started on the time-averaged flow. Its step 159054 (initial state) was mistaken for a result; always use the last step of an energy run.
4. ParaView errors came from mixed folders and `.pvd` files listing missing steps. One folder per run, one-step `.pvd` (§10.2).
5. The owner expected visible changes at t = 35 s: power was still 83 % and the steel responds in 40–100 s. Explain timing before showing snapshots, and use fixed colour ranges.
6. The pulse run did **not** start in thermal equilibrium: the steady → transient discretisation offset (§9.7). Lesson: before a pulse study, run the transient at full power until its own equilibrium, then start the pulse (for example a second run from the end of 31157 with `CBS3D_HEAT_PULSE_PERIOD`).
7. The agent sandbox could not push to GitHub (403) for the whole period; `add_repo` needs the owner to link GitHub in claude.ai settings.
8. Windows `pvpython` is broken (§3.2).

---

## 16. Rules of engagement for the next agent

1. **Understand before changing**: read the files and their call sites. State in one line what you will change and why.
2. **One concern per commit.** Each commit builds, passes tests, and states the measured verification in its message.
3. **Verify every numerical change** on a case with a known answer (slab, U-channel, cavity) **before** the blanket; compare serial and MPI. New implicit operators must reproduce the explicit residual at their solution.
4. **Do not break the explicit solver**: new methods sit behind switches, and the default stays bit-identical.
5. **Exact instructions to the owner** (§1.1). One block per step. Include the expected output and what to send back. Stop where a decision depends on a result.
6. **Accuracy over speed of talk**: say what is unknown, give the test that settles it and its cost. Check PETSc and ParaView documentation instead of guessing option names.
7. **Do not waste cluster hours**: never start a long run that cannot answer a specific question.
8. **Respect the owner's style**: classic C++, ALL-CAPS banners, formal and humble tone, plain-language explanations.
9. **Keep GitHub the source of truth**, and update this file at the end of every session.

---

## Appendix A: command blocks the owner uses (paste into PuTTY/plink)

Session start:
```bash
source ~/cbs3d_env.sh; T=$CBS3D_CODE/jobs/sunbird/templates
cd $CBS3D_CODE && git branch --show-current && git log --oneline -3; squeue -u $USER
```

Live watch of a pulse run (refresh 30 s; Ctrl+C stops watching only):
```bash
P=$(ls -1t $CBS3D_SCRATCH/cases/blanket_slice_h4/run_blanket_*/output/blanket/blanket_probes.csv | head -1); R=${P%/output/blanket/blanket_probes.csv}; S=$(stat -c %Y "$(ls -1tr $R/output/blanket/*.pvtu | head -1)")
while true; do clear; t=$(tail -1 $P | cut -d, -f2); f=$(tail -1 $P | cut -d, -f3); e=$(( $(date +%s) - S ))
awk -v t=$t -v f=$f -v e=$e -v r=$R 'BEGIN{eta=(t>0)?(600-t)*e/t:0; printf "RUN   %s\nSIM   %.1f / 600 s  (%.1f %%)    POWER f = %.2f\nWALL  %.1f h elapsed, about %.1f h to go\n", r, t, 100*t/600, f, e/3600, eta/3600}'
awk -F, 'NR==FNR{for(i=1;i<=NF;i++)c[$i]=i;next}{printf "TEMP  first wall: surface %.1f K, mid-depth %.1f K, corner %.1f K | outlet ch1 %.1f K, ch4 %.1f K\n",$c["FW_surf_C_T"],$c["FW_mid_C_T"],$c["FW_corner_R_T"],$c["C1_out_T"],$c["C4_out_T"]}' <(head -1 $P) <(tail -1 $P)
grep "^Iteration" $R/solver.out | tail -1 | cut -c1-150
echo "FILES $(ls $R/output/blanket/*.pvtu | wc -l) VTU snapshots, $(( $(wc -l < $R/output/blanket/fno_grid/index.csv) - 1 )) FNO samples"; squeue -u $USER -h; sleep 30; done
```

Probe table every 30 s:
```bash
P=$CBS3D_SCRATCH/cases/blanket_slice_h4/run_blanket_31157/output/blanket/blanket_probes.csv
awk -F, 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; printf "%7s %5s %9s %9s %9s %9s %9s\n","t[s]","f","FW_surf","FW_mid","corner","S_c1c2","C1_out"; next} {s=int($2+1e-6); if (s%30==0 && !(s in d)) {d[s]=1; printf "%7.1f %5.2f %9.1f %9.1f %9.1f %9.1f %9.1f\n",$2,$3,$c["FW_surf_C_T"],$c["FW_mid_C_T"],$c["FW_corner_R_T"],$c["S_c1c2_C_T"],$c["C1_out_T"]}}' $P
```

Pack the newest snapshot of a run, with a one-step `.pvd` and the probes, for download:
```bash
R=$CBS3D_SCRATCH/cases/blanket_slice_h4/run_blanket_31157; O=$R/output/blanket
F=$(grep -o 'file="[^"]*"' $O/blanket_distributed.pvd | cut -d'"' -f2 | tail -1); TS=$(grep "file=\"$F\"" $O/blanket_distributed.pvd | sed 's/.*timestep="\([^"]*\)".*/\1/')
echo "SNAPSHOT: $F at t = $TS s, pieces $(ls $O/${F%.pvtu}_rank_*.vtu | wc -l) of $(grep -c 'Piece Source' $O/$F)"
D=$CBS3D_SCRATCH/transfer/pulse_now; rm -rf $D; mkdir -p $D; ln -s $O/$F $O/${F%.pvtu}_rank_*.vtu $O/blanket_probes.csv $D/
printf '<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n  <Collection>\n    <DataSet timestep="%s" group="" part="0" file="%s"/>\n  </Collection>\n</VTKFile>\n' "$TS" "$F" > $D/pulse_now.pvd
tar -chf $CBS3D_SCRATCH/transfer/pulse_now.tar -C $CBS3D_SCRATCH/transfer pulse_now; ls -lh $CBS3D_SCRATCH/transfer/pulse_now.tar
```
(To pick a specific time, replace `tail -1` with an `awk -F'"' '$2>195 && $2<205'` filter on the `DataSet` lines.)

Energy budget over a window (outflow per snapshot, then storage):
```bash
cd $R/output/blanket; for f in $(grep -o 'file="[^"]*"' blanket_distributed.pvd | cut -d'"' -f2 | tail -11); do printf "%s  " $f; python3 $CBS3D_CODE/tools/energy_balance.py $f --plane y=0 2>/dev/null | grep -o "Q_out = [0-9.e+]*"; done
python3 $CBS3D_CODE/tools/solid_energy.py $(grep -o 'file="[^"]*"' blanket_distributed.pvd | cut -d'"' -f2 | tail -11 | sed -n '1p;$p')
```

Resume an interrupted pulse run:
```bash
BLANKET_DIR=$CBS3D_SCRATCH/cases/blanket_slice_h4 bash $T/blanket.sh pulse-resume $CBS3D_SCRATCH/cases/blanket_slice_h4/run_blanket_31157
```

---

## Appendix B: references (verify those marked *)

| Topic | Reference |
|---|---|
| CBS method | Zienkiewicz, Taylor & Nithiarasu, *The Finite Element Method for Fluid Dynamics* |
| AC-CBS ≈ PSPG | Nithiarasu, Bevan & Murali (2013), *Comput. Mech.*, doi:10.1007/s00466-012-0719-5 |
| CBS ≈ GLS, τ ↔ Δt | Codina & Zienkiewicz (2002), *Commun. Numer. Meth. Eng.* 18:99–112 |
| Characteristic Galerkin / CBS stabilisation | Codina, Coppola-Owen, Nithiarasu & Liu (2006), *IJNME* 66:1672–1689 |
| Discrete maximum principle (review) | Barrenechea, John & Knobloch (2024), *SIAM Review* 66(1) |
| SOLD methods | John & Knobloch (2007/2008), *CMAME** |
| Algebraic flux correction | Kuzmin et al., *Flux-Corrected Transport* (Springer) |
| | Giles, Oxford NA report 95/12* |
| Consistent boundary fluxes | Hughes, Engel, Mazzei & Larson (2000), *JCP* 163:467–488* |
| Block preconditioning | Shadid & Tuminaro (2012), *JCP* 231:345–363; petsc-users threads (May 2023 incompressible solvers; April 2023 BoomerAMG AIR) |
| Vortex criteria | Hunt, Wray & Moin (1988), CTR-S88; Jeong & Hussain (1995), *JFM* 285:69–94; Zhou, Adrian, Balachandar & Kendall (1999), *JFM* 387:353–396; Levy, Degani & Seginer (1990), *AIAA J.* 28:1347–1352 |
| Cavity benchmarks | Ghia, Ghia & Shin (1982); Ku, Hirsh & Taylor (1987), *JCP* 70:439–462; Albensoeder & Kuhlmann (2005), *JCP* 206:536–558* |
| Neural operators | Li et al. (2021), Fourier Neural Operator, ICLR; Lu et al. (2021), DeepONet, *Nat. Mach. Intell.* |
