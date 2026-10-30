# CBS3D on Linux and Sunbird: a new 3D case from scratch

This is the Linux replacement for the PowerShell steps in
`LOCAL_TOWER_CBS3D_GMSH_MANUAL.md`. Every command below was run end to end on
Ubuntu (GCC, OpenMPI, PETSc with HYPRE, Gmsh 4.15). The `.par` reference
(section 8 of that manual) is still the authority for what each control means.

## 1. The whole pipeline

```text
  laptop or Sunbird                                       files
  ─────────────────                                       ─────
  1. write geometry ..................................... NAME.geo
  2. gmsh NAME.geo -3 -format msh2 -o NAME.msh .......... NAME.msh
  3. tools/gmsh_to_cbs.py ............................... NAME.plt .bco .material
                                                          + your NAME.par .matprop
  4. tools/mesh_summary.py  -> RESULT: PASS
  ── serial / OpenMP / PETSc ─────────────────────────────────────────────────
  5. cd serial && cbs3dpp_si NAME ....................... output/NAME/*.vtu .pvd .csv
  ── MPI (N ranks) ───────────────────────────────────────────────────────────
  5. tools/export_rank_local_cbs.py --parts N ........... partitions_N/rank_XXXX/
                                                          NAME_rank_XXXX.{plt,bco,
                                                          material,par,matprop,mpi}
  6. tools/dd4d_check_partition.py -> PASS
  7. mpirun -np N cbs3d_parallel NAME partitions_N ...... output/NAME/*.pvtu, per-rank .vtu
```

Steps 1-6 are one command:

```bash
bash tools/cbs_case.sh NAME.geo NAME.par NAME.matprop OUT_DIR --parts N \
    --python ~/venvs/phd/bin/python
```

## 2. What every file is

| File | Made by | Read by | Contents |
|---|---|---|---|
| `NAME.geo` | you | Gmsh | Geometry, mesh sizes, **physical groups** (the BC codes and materials) |
| `NAME.msh` | Gmsh | converter | Mesh, MSH 2.2 |
| `NAME.plt` | converter | solver | `nelem npoin nboun`, tets, node coordinates, boundary faces `n1 n2 n3 parent bc` |
| `NAME.bco` | converter | solver | List of BC codes that occur |
| `NAME.material` | converter | solver | One material id per tet (0 = fluid, 1 = solid) |
| `NAME.matprop` | you | solver | `count`, then `id name phase rho cp k mu Qvol` per material |
| `NAME.par` | you | solver | All solver controls, positional label/value pairs (54 lines) |
| `NAME_rank_XXXX.mpi` | partition exporter | MPI solver | Local-to-global node/element ids, owner of each node, halo SEND/RECV lists |
| `partition_manifest.json`, `partition_counts.csv`, `element_partition.txt` | partition exporter | you | Load balance and checks |
| `output/NAME/NAME_residuals.csv` | solver | `tools/residual_plot.py` | One row per iteration |
| `output/NAME/NAME.pvd` + `.vtu` (serial) or `.pvtu` + per-rank `.vtu` (MPI) | solver | ParaView | Fields |

CBS3D reads no `.dat` file. The Slurm `.slurm` files are only job wrappers.

### Physical groups the converter understands

Boundary surfaces: the numeric tag **is** the CBS BC code.

| Code | Meaning | Code | Meaning |
|---|---|---|---|
| 500 | Wall with prescribed velocity (moving lid) | 510 | Velocity inlet |
| 501, 502 | Constant-temperature wall | 511 | Mass-flow inlet |
| 503, 507, 508 | Other legacy codes (see `Boundary Condition IDs.txt`) | 520 | Pressure outlet |
| 506 | Symmetry | 530 | No-slip adiabatic wall |
| 532 | Heat-flux wall (`heat_flux_bc` in `.par`) | 901, 902 | CHT interface variants |

Volumes: a name containing `fluid` (or tag 10) becomes material 0, `solid`
(or tag 20) becomes material 1. A conformal fluid/solid interface needs no
surface group: the material change across the shared faces is enough
(`cases/fw_straight` works this way).

```text
Physical Surface("inlet_511", 511)  = {6};
Physical Surface("outlet_520", 520) = {4};
Physical Surface("walls_530", 530)  = {1, 2, 3, 5};
Physical Volume("fluid", 1)         = {1};
Mesh.MshFileVersion = 2.2;
```

## 3. Laptop: once

The setup script already built `build/linux-serial` (OpenMP + PETSc) and
`build/linux-mpi` (MPI + OpenMP + PETSc); both write `cbs3dpp_si`. Check:

```bash
ls ~/PhD/CFD/CBS-3D/build/linux-serial/cbs3dpp_si ~/PhD/CFD/CBS-3D/build/linux-mpi/cbs3dpp_si
~/venvs/phd/bin/gmsh --version          # 4.15.x
```

**Use the Gmsh from the Python environment, not Ubuntu's `gmsh` package.**
Ubuntu's build has no Netgen, so `Mesh.Algorithm3D = 4` (Frontal) prints
`Error : No elements in volume 1` and still exits with status 0: you get a
surface-only mesh. `cbs_case.sh` picks the venv Gmsh automatically when you pass
`--python ~/venvs/phd/bin/python` and stops on any Gmsh `Error` line.

## 4. Laptop: smoke test (2 minutes)

```bash
cd ~/PhD/CFD/CBS-3D
phd                                             # activate ~/venvs/phd
bash tools/cbs_case.sh cases/box/box.geo cases/box/box.par cases/box/box.matprop \
    ~/PhD/runs/box --parts 4 --python ~/venvs/phd/bin/python
```

Expected: `1089 tets, 291 nodes, 400 boundary faces`, BC counts
`511=68, 520=66, 530=266`, two `RESULT: PASS`, `DD-4D partition validation: PASS`.

Serial (native CG, `solver_opt = 1`):

```bash
cd ~/PhD/runs/box/serial
OMP_NUM_THREADS=8 ~/PhD/CFD/CBS-3D/build/linux-serial/cbs3dpp_si box
ls output/box        # box.pvd  box_residuals.csv  box_step_00000000.vtu  box_step_00000005.vtu
```

PETSc: change line 2 of `box.par` from `1` to `3` and run the same command.
PETSc options go in an environment variable, because the serial executable
accepts only the case name:

```bash
PETSC_OPTIONS="-pc_type gamg" OMP_NUM_THREADS=8 ~/PhD/CFD/CBS-3D/build/linux-serial/cbs3dpp_si box
```

MPI, 4 ranks:

```bash
mkdir -p ~/PhD/runs/box/run_p4 && cd ~/PhD/runs/box/run_p4
OMP_NUM_THREADS=1 mpirun -np 4 ~/PhD/CFD/CBS-3D/build/linux-mpi/cbs3dpp_si box ../partitions_4
ls output/box        # box_distributed.pvd, .pvtu, per-rank .vtu, box_distributed_residuals.csv
```

Native CG, PETSc and MPI all end with `RelU=5.214855e-03 RelV=1.202367e-01
RelW=1.732350e-01 RelP=4.308202e-01` after 5 iterations. Compare two rank
counts field by field with

```bash
python tools/dd4d_compare_rank_outputs.py run_p2/output/box run_p4/output/box --case box
```

## 5. Laptop: your own 3D case

1. Write `NAME.geo` (copy `cases/box/box.geo` or `cases/fw_straight/mesh/FW_SquareChannel_CHT.geo`).
   Give every boundary surface a physical group with a BC code and every volume
   a `fluid`/`solid` group. Check it by eye: `~/venvs/phd/bin/gmsh NAME.geo`.
2. Write `NAME.matprop` (`id name phase rho cp k mu Qvol`). Layout example,
   helium values illustrative only (take yours from the retained
   `fw_straight.matprop`):
   ```text
   2
   0 helium fluid 6.7 5193 0.27 3.0e-5 0.0
   1 EUROFER solid 7850 500 28.3 0.0 0.0
   ```
3. Copy a `.par` that is close (`cases/box/box.par` laminar, the flat-plate
   `.par` for SA) and change only what you need. The lines you touch most:

   | Line | Block | Typical change |
   |---|---|---|
   | 2 | `solver_opt` | `1` native CG (small meshes), `3` PETSc AMG (anything over ~100k tets, all MPI runs) |
   | 6 | energy | `1` for heat transfer / CHT |
   | 8 | free stream `Ux Uy Uz P T` | must agree with the inlet velocity |
   | 10 | `ntime transient_on dtfixed dtfix iwrite` | iteration count |
   | 16 | `Re Pr Ra Ri` | |
   | 32 | `alpha_sf k_ratio source_solid heat_flux_bc` | CHT source and BC-532 flux |
   | 34-38 | dimensional mode, `p_ref ... mdot ...`, inlet velocity | dimensional runs |
   | 44 | output/monitor | `vtu_output_every_iterations`, `steady_min_iterations` |
   | 48 | `turbulence_on turbulence_model ...` | `1 0 x` SA-noft2, `1 1 x` SA-neg |

4. Build the case, coarse first:
   ```bash
   bash tools/cbs_case.sh NAME.geo NAME.par NAME.matprop ~/PhD/runs/NAME_coarse \
       --clscale 4 --parts 8 --python ~/venvs/phd/bin/python
   ```
   `--clscale 4` multiplies every mesh size by 4 (about 64x fewer tets) so you
   can test the whole chain in minutes. Drop it for the real mesh.
5. Run a few iterations serially, open `output/NAME/NAME.pvd` in ParaView,
   then the MPI case.

Laptop limits (8 cores, 16 GB): about 1-2 M tets with PETSc. Larger meshes
belong on Sunbird.

## 6. Sunbird2

SA2C Sunbird2 (from October 2026): login `login.sunbird.swansea.ac.uk` with
password + MFA code; project `scw7008`; nodes of 2 x 64 cores (128 cores,
256 hyperthreads) and 773 GB; partitions `compute` (3 days), `ondemand`
(12 h), `highmem`. Home (`/pixhome/$USER`) was migrated from old Sunbird;
scratch was not, and is now per project: `/scratch/scw7008/$USER`.

### 6.1 Which mode

| Mode | Executable | When |
|---|---|---|
| One process, 128 OpenMP threads, PETSc AMG | `build/serial-petsc-release/cbs3d_serial` | Up to a few M tets; no partitioning needed |
| MPI, 128 ranks per node | `build/mpi-release/cbs3d_parallel` | Production; strong scaling; anything large |

### 6.2 Once: PETSc, Gmsh, environment, build

Sunbird2 has no PETSc module, so build it once in `$HOME/software`
(20-40 min, login node, because configure downloads HYPRE and METIS):

```bash
mkdir -p ~/CBS3D_SB2 && cd ~/CBS3D_SB2
git clone https://github.com/prajwal1798/CBS-3D.git code
T=~/CBS3D_SB2/code/jobs/sunbird/templates

mkdir -p ~/software && cd ~/software
nohup bash $T/build_petsc_sunbird2.sh > build_petsc_sb2.log 2>&1 &
tail -f build_petsc_sb2.log                   # ends with "PETSC BUILD DONE"

module load python/3.14.6                     # Python venv with the gmsh module
python3 -m venv ~/venvs/cbs && ~/venvs/cbs/bin/pip install gmsh
module unload python/3.14.6

cp $T/sunbird_env.sh.example ~/cbs3d_env.sh   # the one site-specific file
mkdir -p /scratch/scw7008/$USER
bash $T/build_cbs3d.sh                        # OK serial-petsc-release, OK mpi-release
```

Check with the box case (seconds, login node is fine):

```bash
source ~/cbs3d_env.sh && cd $CBS3D_CODE
bash tools/cbs_case.sh cases/box/box.geo cases/box/box.par cases/box/box.matprop \
    $CBS3D_SCRATCH/runs/box --parts 4 --python ~/venvs/cbs/bin/python --overwrite
cd $CBS3D_SCRATCH/runs/box/serial && $CBS3D_SERIAL_EXE box | grep "Final RelU"
mkdir -p ../run_p4 && cd ../run_p4
mpirun -np 4 $CBS3D_MPI_EXE box ../partitions_4 | grep "Iteration 5/5"
```

Both give `RelU=5.214855e-03 RelV=1.202367e-01 RelW=1.732350e-01
RelP=4.308202e-01` (checked 29 Oct 2026 with GCC 13.5, OpenMPI 5.0.10,
PETSc 3.26).

### 6.3 Every case

Gmsh, the converter and the exporter all run on Sunbird2, so a case can be
built there directly. Small meshes on the login node:

```bash
source ~/cbs3d_env.sh && cd $CBS3D_CODE
bash tools/cbs_case.sh GEO PAR MATPROP $CBS3D_SCRATCH/cases/NAME \
    --python ~/venvs/cbs/bin/python --threads 8
```

Large meshes (above ~0.5 M tets) on a compute node, mesh to partitions in one job:

```bash
T=$CBS3D_CODE/jobs/sunbird/templates
sbatch $T/build_case.slurm GEO PAR MATPROP $CBS3D_SCRATCH/cases/NAME --parts 128
```

Or only new partitions of an existing serial case:

```bash
cd $CBS3D_SCRATCH/cases/NAME
sbatch $CBS3D_CODE/jobs/sunbird/templates/export_partitions.slurm serial/NAME 128 partitions_128
```

(Or build it on the laptop and `rsync` the `serial/` or `partitions_N/`
folder up.)

Submit from the case folder:

```bash
T=$CBS3D_CODE/jobs/sunbird/templates
sbatch $T/run_mpi.slurm partitions_128 NAME                          # 128 ranks, 1 node
sbatch --nodes=2 --ntasks=256 $T/run_mpi.slurm partitions_256 NAME   # 256 ranks, 2 nodes
sbatch $T/run_openmp_petsc.slurm serial NAME                         # 1 process x 128 threads
```

The templates request `--hint=nomultithread` (one rank per physical core),
`--exclusive` and `--mem-per-cpu=2G` (256 GB or more per node; the
2.3 M-tet cavity needs about 11 GB in total). Do not use `--mem=0`: it asks
for more than the per-CPU memory limit and the job pends forever with
`(MaxMemPerLimit)`.

#### Any number of ranks, whole nodes or scattered free cores

The rank count is fixed when the mesh is partitioned (N partitions = N
ranks). `submit_mpi.sh` makes `partitions_N` if it does not exist yet (a
partition job) and submits the run behind it:

```bash
T=$CBS3D_CODE/jobs/sunbird/templates
bash $T/free_cores.sh                                    # idle cores per node, right now
bash $T/submit_mpi.sh $CBS3D_SCRATCH/cases/NAME NAME 256           # 2 whole nodes
bash $T/submit_mpi.sh $CBS3D_SCRATCH/cases/NAME NAME 300 any       # 300 cores anywhere
bash $T/submit_mpi.sh $CBS3D_SCRATCH/cases/NAME NAME 256 nodes 24:00:00
```

- `nodes` (default): ceil(N/128) exclusive nodes. Fastest per rank; use it
  whenever whole nodes are idle.
- `any`: `run_mpi_anycores.slurm`, no `--nodes`, no `--exclusive`,
  2 GB per core. Slurm fills idle cores on partly used nodes (20 here, 30
  there). Ranks then share memory bandwidth with other jobs, and the slowest
  rank sets the pace, so expect less speed per rank than on whole nodes.
Prepare several sets once, then let the script pick the largest set that
fits the idle cores at submission time:

```bash
bash $T/prepare_partitions.sh $CBS3D_SCRATCH/cases/NAME NAME 40 80 128 256 384
bash $T/submit_auto.sh $CBS3D_SCRATCH/cases/NAME NAME --dry-run      # see the choice
bash $T/submit_auto.sh $CBS3D_SCRATCH/cases/NAME NAME --time 24:00:00
```

Whole exclusive nodes are used only when N fills them (N = k x 128) and k
nodes are idle; otherwise the ranks go on scattered idle cores. Long steady
runs should write checkpoints (`export CBS3D_CHECKPOINT_EVERY=10000` before
submitting); a restart must use the same N as the checkpoint.

- More ranks help only while each rank keeps enough work: below roughly
  10-20 k tets per rank, communication dominates. For the 2.3 M-tet cavity
  that is about 128-256 ranks; a finer mesh (`lc = 0.010`, ~4.5 M tets)
  justifies more.

Watch and fetch:

```bash
squeue -u $USER
tail -f run_NAME_<jobid>/solver.out
sacct -j <jobid> --format=JobID,Elapsed,State,MaxRSS
rsync -av sunbird:/scratch/scw7008/<user>/cases/NAME/run_NAME_<jobid>/output/ ~/PhD/runs/NAME/sunbird/
```

`run_mpi.slurm` hard-links the partitions into a private run folder, gives each
rank its own `.par`, records the git commit and checksums, and launches
`srun --mpi=pmix ... cbs3d_parallel NAME <run>/partitions`.

### 6.4 Pitfalls found on Sunbird2

- `~/cbs3d_env.sh` starts with `module purge`; the venv Python needs the
  `python/3.14.6` module loaded, so the env file loads it.
- The pip `gmsh` command is a `#!/usr/bin/env python3` script: the venv must
  be activated (the env file does this), not just its `python` used.
- The OpenMPI 5.0.10 module also loads `cuda/13.0.3`; this is harmless.
- `make check` for PETSc runs from the source tree (`~/software/petsc-src`),
  not from the install prefix.

## 7. How the partitioned case is built

`tools/export_rank_local_cbs.py` reproduces the original tower exporter. On
`cases/flat_plate_p8`, given the same element-to-rank map, it writes
`.bco`, `.material`, `.par`, `.matprop` and the `.mpi` topology byte-identical
to the shipped files; `.plt` differs only in the last digit of 4 coordinates
per rank (the originals were re-read from Gmsh's binary file). Rules:

- METIS (through `gmsh -part N`) assigns each tet to a rank; rank = Gmsh partition id - 1.
- A node belongs to the lowest rank that has a tet on it; local numbering is
  owned nodes, then ghost nodes, each ascending by global id.
- Each rank keeps its own tets only (no ghost cells) and the boundary faces whose parent tet it owns.
- `SEND` to rank q: owned nodes that q holds as ghosts. `RECV` from q: ghosts owned by q.
- Global ids in `.mpi` are the serial-case ids, so serial and MPI results map onto each other.

## 8. Pitfalls found while writing this

| Symptom | Cause | Fix |
|---|---|---|
| `.msh` has only triangles, Gmsh exit code 0 | Ubuntu Gmsh lacks Netgen (`Algorithm3D = 4`) | Use `~/venvs/phd/bin/gmsh`; `cbs_case.sh` now stops on Gmsh errors |
| `gmsh in.msh -part N -o out.msh` writes nothing | Gmsh 4.12+ needs `-save` when re-saving a loaded mesh | The exporter tries `-save`, then without (old Sunbird 4.4) |
| `.geo` + `-save` gives an empty mesh | `-save` writes before meshing | Never combine `-3` with `-save` |
| `ERROR: Usage: cbs3d_serial <case_name>` | PETSc flags on the command line (old `docs/hpc/sunbird.md` template) | `PETSC_OPTIONS="..."` |
| `mesh_summary.py: cannot parse id and coordinates` on CHT meshes | A negative y or z touches the previous number (`0.0e+00-7.5e-03`) | Parser fixed; the solver always read these correctly |
| `mesh_summary.py` fails on Sunbird | `from __future__ import annotations` needs Python 3.7 | Removed; it runs on 3.6 again |
| `export_rank_local_cbs.py not found` | The tower exporter lived outside the repo | Now `tools/export_rank_local_cbs.py` |
| `make USE_PETSC=1 ...` fails | `docs/hpc/sunbird.md` still describes a Makefile build | Use the CMake presets (6.2) |
| Gmsh on Sunbird rejects the first-wall `.geo` | `gmsh/4.4.0` module has no OpenCASCADE | Use `FW_SquareChannel_CHT_gmsh44.geo`, or mesh on the laptop |
| Rank-local case run serially | It holds one partition only | Always `mpirun -np N` with all ranks |
