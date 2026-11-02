# CBS3D++ numerical verification gates — 2026-11-02

## Rules

- Baseline: b8f7d8d08e829a746eeac34d5a8654e7f2e09dd2.
- This branch adds verification tests and documentation only. Do not modify production numerical operators.
- AFC OFF and experimental implicit energy OFF until baseline is verified.
- A test passes only after a real executable produces recorded results; formulas in comments are not evidence.
- Separate (i) consistency with PDE, (ii) code correctness, (iii) solution accuracy, and (iv) physical validation.
- Diagnose the first failure on the smallest problem before trying an alternative formulation.

## Checkpoint register

| Gate | Required verification | Pass criterion | Status |
| --- | --- | --- | --- |
| G0 | Baseline provenance: commit, binary SHA256, mesh, boundary data, materials, timestep flags | Reproducible case snapshot | Commit and Job 29706 identified; full case archive pending |
| G1a | Actual Step-4 P1-tet FEM volume and surface terms vs independently evaluated quadrature | All nodal convection, conduction, characteristic, heating, Neumann contributions agree | PASS: Sunbird2 Step-4 operator verification CTest, 2026-11-02 |
| G1b | Exact frozen-velocity element linear operator versus finite-difference residual derivative | All 16 entries agree | PASS: Sunbird2 Step-4 operator verification CTest, 2026-11-02 |
| G1c | Lumped thermal mass and timestep diagonal | Code vs direct P1-tet reference with positive entries | PASS: Sunbird2 CTest 3/3 (2026-11-02), shared-node and material masks |
| G2 | Actual inlet/outlet/wall/interface boundary inventory and flux sign | Every applied BC and flux accounted for | PASS (structural): Job 29706 has only BC 511/520/530/532, 648 BC532 faces = 0.0006 m^2, 60 W applied, zero duplicate faces, zero MPI duplicate mismatch; full global energy closure remains G7 |
| G3 | Solid-only diffusion | Uniform/linear patch tests; analytical heated slab; mesh refinement | Analytical slab and serial/MPI PASS, Job 29742; systematic mesh refinement still pending |
| G4 | Two-material conduction | Shared interface temperature, weak flux continuity, integrated balance | PASS: Job 29749 reached thermal steady state, max analytic error 3.681e-5 K, two-sided flux 999.978/999.984 W/m², serial/MPI 1.11e-15 K, qA=0.016 W |
| G5 | Prescribed-flow scalar transport | Uniform T, analytical convection-diffusion and manufactured solutions | G5a PASS (29754): 3-grid steady Pe=3, bound-preserving exact solution; G5b implementation PASS (29756): production mass/RHS/matrix/first-step consistency, but unlimited update fails DMP on synthetic meshes; full coupled/transient MMS pending |
| G6 | Discrete positivity / DMP | Inspect assembled update matrix and *first* bound violation | CONFIRMED non-DMP on actual U-channel: Job 29763 node 5707 has 5 negative neighbour links (total -0.310272 K) vs 8 positive links (+0.211759 K); first-step undershoot -0.098513 K reproduced to ~3e-15 K. Galerkin advection dominates negative links. |
| G7 | Integrated energy and mass | Sources + boundary heat = advection + storage within stated tolerance | G7a PASS: actual U-channel first-step (50000→50001) thermal balance 2.37e-11 W abs / 3.17e-13 relative after including nuclear Q=16.246255 W. Max unconstrained nodal discrepancy 5.60e-12 W. Net mass imbalance 0.503%; long-run steady energy/mass validation outstanding. |
| G8 | Serial/MPI and restart parity | Independent solutions agree at declared tolerances | G8a PASS (29786); G8b PASS (29788) one-step 2-vs-4 ranks; G8c PASS (29789) 100-step 2-vs-4 ranks max dT 7.33e-12 K, du 2.52e-13 m/s, dp 1.88e-13 Pa; serial comparison pending. |
| G9 | Time/iteration convergence | Stationary residual, dt-variation and stopping robustness | Pending |
| G10 | Grid convergence of quantities of interest | Tmax, interface flux, wall gradient and outlet enthalpy | Pending |
| G11 | Additional stabilisation only if needed | Isolated failure, specific literature choice, accuracy + bounds verified | G11a one-step AFC PASS (Job29832) after single-call source correction: AFC call=1; Tmin OFF 298.051487K vs ON 298.15K; 94 cold nodes ->0; pressure/velocity unchanged; correction balance 3.23e-13. Long-run accuracy/conservation pending. |
| G12 | Implicit acceleration only after verification | Same discrete result, lower runtime, robust nonlinear residual | G12d isolated stationary low-order PETSc baseline PASS (Job29859): 49 GMRES, 0.0356358 s, residual Linf 5.58048e-11, Tmin=298.15K and Tmax=375.849K over 9613 nodes. This is graph-viscosity low-order scheme, not high-order AFC; G12b/c nonlinear high-order attempts failed. Energy and high-order QoI accuracy not yet verified. |

## Literature anchors

- Chun-Bin Liu (2005), Swansea PhD thesis, chapter 4 and appendix B: CBS flow formulation and P1 tet matrices.
- Codina and Zienkiewicz (2002), doi:10.1002/cnm.470: characteristic time-step stabilisation relationship; not a positivity proof.
- Codina et al. (2006), doi:10.1002/nme.1697: CBS versus SGS numerical comparisons; not direct thermal CHT validation.
- Giles (1997), Oxford report: stability of fully and loosely coupled diffusion interfaces.
- Barrenechea, John and Knobloch (2023), SIAM Review: DMP, time-dependent positivity, mesh constraints, AFC.
- Explicit thermal FEM volume, source and flux equations follow the P1 weak form; any custom divergence or stabilisation terms need their own derivation and reference mapping.

## G1 execution on SUNBIRD2

~~~bash
source ~/cbs3d_env.sh
cd "$CBS3D_CODE"
git fetch origin
git switch audit/cbs-literature-verification-20261102
cmake -S . -B build/verification -DCBS3D_BUILD_SOURCE_VALIDATION=ON -DCBS3D_BUILD_TESTS=ON -DCBS3D_ENABLE_MPI=OFF -DCBS3D_ENABLE_PETSC=OFF
cmake --build build/verification --target cbs3d_step4_operator_verification cbs3d_source_validation cbs3d_lumped_thermal_mass_verification -j 8
(cd build/verification && ctest --output-on-failure -R '^step4_')
~~~

Expected: three CTests pass. G1a/G1b already passed on Sunbird2 (2/2 tests, 2026-11-02). G1c PASSED on Sunbird2 on 2026-11-02 (3/3 CTests). It checks the production inverse thermal diagonal and fluid-only momentum/pressure masks at a shared solid-fluid node. Full interface heat-flux continuity remains G4.

If a test fails, halt and diagnose the specific element contribution. Do not patch production code before an independent reproduction.

## Change acceptance criteria

1. State failing input, governing equation, predicted value and observed value.
2. Identify first affected node/element and local operator terms.
3. Reproduce independently on a smaller problem.
4. Make a single targeted change on a separate branch; record before/after results.
5. Re-run applicable gates, especially serial/MPI and analytical benchmarks.

## G2 audit-tool compatibility

The serial VTU writer outputs CellData `mat_id`; DistributedPost.cpp writes `material_id` for MPI. `tools/audit_cht_temperature_vtu.py` previously accepted only `mat_id` and therefore failed before examining the Job 29706 result. The audit now accepts both and checks conflicting dual fields. Synthetic regression tests: `python3 -m unittest discover -s tests -p test_audit_cht_temperature_vtu.py -v`. No thermal numerical operator was modified. Do not mark G2 passed until the actual boundary-face mapping, power and final field are inspected.

## G2 evidence (2026-11-02, Sunbird2)

- Read-only VTU parser: 5/5 regression checks passed.
- Job 29706 BC faces/area: 511: 44 / 2.5e-5 m2; 520: 44 / 2.5e-5 m2; 530: 8124 / 7.55e-3 m2; 532: 648 / 6.0e-4 m2.
- Dirichlet boundary IDs prescribing zero or one temperature did not appear in this output. Reference/free-stream and inlet temperatures are 298.15 K.
- BC532 q=100000 W/m2 produces 60 W input. No duplicate boundary triangles, no mismatched temperatures on duplicated MPI nodes.
- Global Tmin=235.4509136 K, Tmax=440.3715889 K, still unexplained. Not proof of a bad BC. Convection/diffusion numerical checks and independent energy closure remain pending.
- Heat-load/capacitance coefficient of variation 0.40297 is a diagnostic, not a failed physical condition. The 836.7166 J relative stored energy is not comparable to 60 J from arbitrary local pseudo-time.
- The pre-existing read-only tools/audit_temperature_undershoot_vtu.py locates and classifies cold owned nodes; its time and threshold column labels were clarified. Next: run this on Job 29706 output; no production modification.

## G2 cold-node localisation (2026-11-02)

- Job 29706 restart snapshot at global iteration 50,000: Tmin 298.15 K, Tmax 398.241306 K; no nodal temperatures below 298.15 K.
- After 50,000 AFC-off iterations (global iteration 100,000): Tmin 235.4509136183 K, Tmax 440.3715888721 K; 180 nodes below Tref; 176 below 298.05 K.
- All 176 cold nodes are fluid-only. Two touch outlet BC520; remaining cold nodes have no boundary-face incidence (classification histogram has BC520:2). Min node global ID 5792, xyz=(0.0461112664048,0.03861369882182,0.003859505728298) m, interior fluid, |u|=0.02508956 m/s. Max node 608, solid heated wall adjacent BC530+532.
- These are empirical diagnostics, NOT proof of a defect in CBS characteristic stabilisation or thermal boundary assembly. Next verify the pure conduction slab and then prescribed-flow scalar transport without flow coupling. The earlier solver log *does* show Tmin=298.051 K at local iteration 1/50000; therefore at least one node crosses 298.15 K during the first AFC-disabled update. That first-step node's ID and term-by-term residual are not available from the current VTU checkpoints.

## G3 isolated non-AFC benchmark

- Runner: jobs/sunbird/templates/verify_slab_g3_no_afc.slurm (new; not executed as of writing).
- Analytic stationary EUROFER slab solution: T(y)= (A/k)[lambda*y-lambda^2(exp(-(L-y)/lambda)-exp(-L/lambda))]+q*y/k.
- Parameters L=.10 m, W=.01 m, k=28.3 W/(m K), A=1.6e6 W/m3, lambda=.05 m, q=1e5 W/m2.
- Comparison: tools/check_nuclear_slab.py; relative maximum nodal error criterion 1%, power relative error 0.1%, and serial/4 MPI nodal parity 1e-6*abs(T(L)). Both jobs must explicitly have AFC=0 and implicit steady energy unset.
- This benchmark tests a disconnected solid conduction bar with a passive fluid pressure reference. It does NOT test the fluid-solid interface, which remains G4, nor mesh convergence, which remains open.

## G3 measured evidence (Job 29742, 2026-11-02)

- sacct COMPLETED 0:0, 20 s; git checkout aae523d. Thermal AFC=0, local-timestep monotone=0, no restart or implicit steady energy.
- Serial executable SHA256 2bc843a610a18d9e166cbb0b859c0b72d67e427ef01c162e938aef42730c393a.
- MPI executable SHA256 96ca7134f6affc2cd17cdfb6b14ba70bbda9f243d7041ba196d70a00bc1c16fa.
- Exact T(L)=513.8283 K. Max abs nodal error=0.21482 K (4.181e-4 of T(L)), same for serial and 4 MPI ranks.
- Printed Q=16.91597 W; exact Q=16.91732 W; relative error 7.97e-5.
- Serial/MPI maximum nodal difference=1.137e-13 K, matched step 89434.
- Result: G3-NO-AFC PASS. This verifies the particular stationary homogeneous-solid problem but NOT mesh-convergence order, fluid-solid interface heat-flux continuity, or time accuracy.

## G4 two-material conduction verification design

- Input mesh: cases/two_layer_conduction/cht2.geo; generated using OpenCASCADE BooleanFragments for conformal shared fluid/solid interface. The internal material interface must not become a boundary face.
- Materials: water k_f=0.606 W/(m K), EUROFER k_s=28.3 W/(m K). Stagnant fluid as mathematical conduction region, no advection or nuclear source.
- Geometry: W=0.004 m, L_f=0.001 m, L_s=0.003 m. Cold fluid face BC502 T=0 (numerical temperature-offset test; not intended as physically realistic absolute water temperature); heated solid end BC532 has inward q=1000 W/m2. BC520 on a fluid side provides pressure reference; remaining exterior insulated BC530.
- Exact stationary solution T(y)=q*y/k_f for y<=L_f, T(y)=q*L_f/k_f+q*(y-L_f)/k_s above interface.
- Exact integrated input Q=qW^2=0.016 W. Nodal max error threshold 0.5% of exact top temperature, zero intended bulk flow, reference area within 0.1%, surface pressure BC inventory and shared interface-node connectivity required. Median reconstructed k*dT/dy at interface-adjacent tetrahedra must be within 5% of prescribed flux on both sides.
- Runner: jobs/sunbird/templates/verify_two_layer_g4_no_afc.slurm, AFC disabled and no implicit energy solver. Checker: tools/check_two_layer_conduction.py. This is a new unexecuted verification case; do not claim pass until logs show it.
- Note: interface-near k*dT/dy medians are diagnostic; pointwise element fluxes need not be continuous across a linear FEM interface even when the global weak heat balance is correct.

## G4 measured status — Job 29745 (2026-11-02)

- sacct: FAILED 1:0 because the **verification runner** flagged analytical comparisons; its serial/MPI solver substeps completed 0:0 according to supplied output.
- Tested source checkout 923ca68; AFC=0, LTS_MONOTONE=0; same pre-existing serial/MPI executable hashes as G3.
- Geometry checks: 337 unique nodes, 354 fluid tet, 763 solid tet, 58 shared interface nodes; BC502 and BC532 each 1.6e-5 m²; four intended BC types only; q*A=0.016 W exactly.
- Analytic T_int=1.650165017 K, T_top=1.756172084 K. Measured Tmax=1.334267022 K; max analytic nodal error=0.422242384 K (24.04334%). Interface reconstructed k*dT/dy medians: fluid=751.164742 W/m², solid=811.731062 W/m² vs 1000 W/m².
- Both serial and 4-rank outputs give the same printed values; exact maximum difference has not been reported by the verifier because it stops at failed analytic tests. Do not claim quantified serial/MPI parity until measured.
- CRITICAL UNRESOLVED QUESTION: did simulations stop at ntime=300000 without thermal steady convergence? Examine both solver.out iteration and stop reason, RelT and temperature history, actual effective timestep and residual CSV before changing any solver operator. This experiment is not yet evidence of defective interface assembly. Steady state analytic profile is inapplicable to an unconverged transient/pseudotransient iterate.
- Separate issue: median cellwise interface-near reconstructed fluxes are not the same as the integrated weak interface-flux balance. Ultimately check both after convergence.

## G4 convergence diagnosis and rerun controls (2026-11-02)

- Job 29745 serial and MPI final stop reason: **maximum iteration count reached** at 300,000. G4 analytic steady-state comparison at that checkpoint is NOT valid as a test of stationary interface physics.
- RelT (300,000)=1.520100200195e-6, target 1e-10; dt reported 0.000107266491. Flow residuals and velocity zero.
- MPI Tmax still increases: 1.23549 K at 255,000 -> 1.33427 K at 300,000. RelT decreases from 2.03055e-6 to 1.52010e-6, without a plateau; expected stationary value 1.756172 K.
- Trial estimate from relative residual 299,600 to 300,000 gives a characteristic decay of order 160k iterations; this is a local fit, NOT an extrapolation guarantee.
- The exact binary/matprop/mesh and CFL settings remain unchanged. Separate case input cases/two_layer_conduction/cht2_convergence.par raises only ntime to 2,500,000, retains 1e-10 temperature convergence threshold. G4 runner now points to that .par and requires explicit 'steady-state convergence reached' in BOTH solver logs before interpreting numerical temperature/flux accuracy.
- No changes to fluid/solid FEM, interface coupling, Neumann integration, pressure equation, or AFC implementation. G4 not passed; require rerun results.

## G4 rerun intermediate evidence — Job 29749 (2026-11-02)

- User-reported serial stage has terminated with **steady-state convergence reached**, final RelT=9.999981e-11 (target 1e-10). The analytic output error and interface flux check are not yet available in the main G4 summary.
- Four-rank MPI stage still running at iteration 1,535,000/2,500,000; RelT=3.57e-9, Tmax=1.75486 K, Umax=0, T_min=0. The exact stationary top temperature is 1.756172084 K: top-temperature difference at this intermediate step is ~0.001312 K (0.0747%).
- MPI residual is decaying consistently; no claim of complete G4 PASS until MPI reaches thermal steady convergence and both numerical fields satisfy nodal temperature and two-sided interface heat-flux checks.

## G4 verified evidence — Job 29749 (2026-11-02)

- Job 29749: Slurm COMPLETED 0:0, 00:10:08; checkout 4b58dad, with AFC=0, LTS monotone=0, restart=0, implicit steady solver=0.
- Existing serial binary SHA256 2bc843a610a18d9e166cbb0b859c0b72d67e427ef01c162e938aef42730c393a; MPI binary SHA256 96ca7134f6affc2cd17cdfb6b14ba70bbda9f243d7041ba196d70a00bc1c16fa.
- Both runs reached 'steady-state convergence reached' (serial final RelT=9.999981e-11; MPI stop criterion enforced by runner; final MPI residual value not pasted).
- 337 unique nodes, 354 fluid and 763 solid tetrahedra, 58 conformal shared interface nodes. BC502 and BC532 each have area 1.6e-5 m². Exact applied power qA=0.016 W; both solver outputs print 0.016 W.
- Analytical temperatures: T_interface=1.650165017 K, T_top=1.756172084 K; numerical Tmax=1.756135305 K. Max nodal error 3.680762630642e-5 K (0.00210% of analytic top). Max interface T error 3.581702319710e-5 K.
- Reconstructed median interface-adjacent conductive k*dT/dy flux: fluid=999.9783086 W/m²; solid=999.9835883 W/m²; exact=1000 W/m². Zero fluid speed.
- Maximum serial/MPI nodal difference 1.110223024625e-15 K. Runner reports **G4-CONFORMAL-NO-AFC PASS**.
- Prior Job 29745 stopped at 300k iterations without steady convergence, explaining its premature comparison; no numerical production code change made.
- Scope: G4 verifies the stationary, stagnant two-layer conformal interface for this mesh and flux. It does not verify transient accuracy, mesh-convergence order, or advection; proceed to G5.

## G5a prescribed-velocity scalar-transport verification (prepared; execution pending)

- C++ source: `tests/PrescribedVelocityThermalTransportVerification.cpp`; target: `cbs3d_g5_prescribed_transport_verification`; CTest: `g5_prescribed_velocity_transport`.
- Uses **production** `EnergyAssembly::elementOperator` and `EnergyAssembly::assembleStep4Rhs` on conformal structured P1 tetrahedral meshes, **not** a rewritten element operator. Its small test-only dense partial-pivot linear solver computes the steady assembled temperature field.
- Mathematical PDE: rhoCp U dT/dx = k d²T/dx², `U=(0.3,0,0)`, `rhoCp=1`, `k=0.1`, `x in [0,1]`, `y,z in [0,0.25]`; Dirichlet T(0)=1, T(1)=2, lateral homogeneous Neumann.
- Exact stationary solution `T(x)=1+expm1(3*x)/expm1(3)` (Peclet 3); test imposed tau=dt/2 using `dt=0.3*h_x/U` to retain CBS characteristic diffusion. The finite-element approximation is expected to approach the physical continuum result under refinement with O(h) streamline-diffusion truncation, no claim of second-order rate.
- Three meshes `nx=8,16,32`, with `ny=nz=2`. Expected reference max nodal errors based on independent Python matrix check ~0.01615, 0.00899, 0.00474 K; these are **predictions**, not validated C++ outputs. Acceptance: max error reduces by factors >1.3 and >1.4, finest max error <0.0075 K, nodal T remains [1,2], and production residual and uniform-field invariance check < specified tolerances.
- Test imposes temperature Dirichlet rows algebraically, independent of production BC handling; G2 already addressed BC mapping, but G5a is not an end-to-end `.bco` run and cannot establish the correctness of real production solver time stepping, velocity update or parallel halo.
- Job runner: `jobs/sunbird/templates/verify_g5_prescribed_transport.slurm`, which compiles only test target with MPI/PETSc disabled, runs CTest verbose and prints `G5A ISOLATED TRANSPORT PASS` only after all checks succeed. This job has **not yet executed** on Sunbird2.
- If G5a passes, proceed to G5b: imposed flow, first-step CBS explicit positivity and time-step sensitivity (and eventually a fully coupled fluid calculation). If it fails, inspect matrix/RHS and temporal coefficient before altering production numerics.

## G5a observed verification results (Job 29754, 2026-11-02)

- Slurm COMPLETED 0:0 in 4 seconds. Source checkout 11f1b530a1815ddc983aada448cd61dd29fa33f2. GCC 13.5.0, MPI=OFF, PETSc=OFF.
- Source tested includes actual `EnergyAssembly::elementOperator` and `assembleStep4Rhs` with imposed constant `U=(0.3,0,0)`, Pe=3, x in [0,1], Dirichlet temperature [1,2], insulating lateral surfaces, uniform rhocp=1 and conductivity k=0.1.
- Measured maximum analytic nodal errors for nx 8/16/32: 0.0161499349760 K, 0.00898798833336 K, 0.00474007524202 K; error ratios 1.79683532922 and 1.89616997082.
- Maximum steady residual across three grids <= 1.05926552252e-16, uniform-state residual 0, Tmin=1 and Tmax=2. All G5a CTest assertions PASS; job printed `G5A ISOLATED TRANSPORT PASS`.
- Scope: the exact steady operator, not the explicit iterative update, was exercised; it does not show DMP for arbitrary fields or for the U-channel velocity field.

## G5b first explicit-update positivity audit (prepared; execution pending)

- C++: `tests/G5ExplicitThermalUpdateAudit.cpp`; CTest `g5_explicit_update_monotonicity_audit`; runner `jobs/sunbird/templates/verify_g5_explicit_update.slurm`.
- Calls existing production `TimeStep::updateLhsDiagonal`, `EnergyAssembly::elementOperator` and `EnergyAssembly::assembleStep4Rhs`. All test mesh, material data and prescribed velocities are isolated, and neither solver numerical code nor production binaries are changed.
- Multitet P1 fluid meshes nx=8/16 with 2 x 2 transverse hexahedral subdivisions (6 tetrahedra per subdivision), U=0.3, rhocp=1, k=0.1 or 0.001 -> Pe_L=3 or 300. Courant_x=0.1,0.3,0.6, using dt=Courant/(nx U). Uniform reference temperature Tref=298.15 K.
- At each combination, independently checks production reciprocal lumped heat capacity/time diagonal against sum_e rhoCp V/(4dt), verifies uniform-temperature invariance for nonzero divergence-free velocity, assembles the linearized explicit one-step matrix B=I+diag(elcoe2p)*dR/dT, and counts negative offdiagonal entries among unconstrained interior nodes.
- If a negative B_ij exists, use a bounded initial temperature hat T_j=Tref+25 K, all other nodes Tref. The actual production residual and explicit update must reproduce the matrix-predicted first-step undershoot to 1e-9 K. This is a **discrete maximum-principle (DMP) counterexample** for the tested unlimited stencil, but is NOT a proof of long-time instability or of identical causation in U-channel Job 29706.
- Test PASS signifies mass/RHS/matrix/update consistency only; observed negative couplings or overshoots are reported diagnostics, not grounds to mark the numerical test executable as failing. Positivity preservation is assessed separately under G6.
- If G5b exhibits negative offdiagonals, follow with input-specific Job 29706 first-step stencil/mesh and timestep diagnostics; do not directly change CBS stabilization without identifying their origin.

## G5b numerical evidence: Job 29756 (2026-11-02)

- Slurm completed 0:0 in 2 s, source commit 50a164aa671ccfbbc9b8226a297688fe2d3049bd, GCC 13.5.0, MPI/PETSc off. 12/12 test combinations passed diagnostic consistency checks. NO production numerics were changed.
- All cases (nx=8,16; Pe_L=3 or 300; Courant_x=0.1,0.3,0.6): independent lumped thermal diagonal error zero; uniform-field residual zero; direct production assembled residual and B=I+diag(elcoe2p)*dR/dT predictions agree exactly for one-step hat perturbations; no erroneous matrix/residual inconsistency found.
- Full explicit update B has NEGATIVE interior offdiagonals. Counts for nx=8 Pe=3: 124 at all three Courant choices; nx=8 Pe=300: 262,262,208; nx=16 Pe=3: 284; nx=16 Pe=300: 590,590,464. Worst B_ij=-0.025, -0.075, -0.150 for Courant_x=0.1,0.3,0.6 respectively (for both Pe and nx shown).
- A *bounded* hat perturbation T_j=323.15 K and all other nodes=298.15 K yielded T_i^{new}=297.525,296.275,294.400 K at those Courant values respectively (undershoot below Tref of 0.625,1.875,3.750 K). This is a direct counterexample to unconditional DMP of the **unlimited explicit** thermal scheme on tested meshes. It is NOT evidence of long-time unbounded instability or proof of the first-step mechanism at U-channel Job 29706.
- This result separates **G5b implementation consistency PASS** from **G6 DMP FAIL (specific synthetic configurations)**. No proposed production numerical correction until the original U-channel first-step stencil and actual effective dt field are captured and analyzed.
- Next: inspect the Job 29706 run_record, flow_50000.par (ilots, dtfixed, csafm, htype), restart metadata, AFC switches, and iteration-1 residual extrema. Before reproducing its exact first step, preserve original run inputs, compile/binary hashes and partitioning; do not assume VTU snapshot velocity is the corrected velocity within the next CBS iteration.

## G6 original U-channel replay provenance and next experiment (2026-11-02)

- The U-channel Job 29706 run_record.txt was located after correctly sourcing ~/cbs3d_env.sh. Case root /scratch/scw7008/s.2337862/cases/uchannel_afc_audit; original running directory run_uchannel_29706. First 50k AFC-on restart snapshot came from run_uchannel_29580/output/uchannel, 50000 iterations completed. Next resumed global step is 50001.
- The master flow_50000.par specifies ntime=50000, transient_on=0, dtfixed=0, csafm=0.50, cbs_scheme=1, **ilots=0**, htype=1, art_diff=0. This is GLOBAL UNIFORM pseudo-time, not local LTS; export CBS3D_LTS_MONOTONE=1 with margin 0.8 is INACTIVE in this configuration.
- Run record: CBS3D_TEMPERATURE_AFC=0, CBS3D_STEADY_ENERGY=0, CBS3D_RESTART=1, format=legacy_vtu, iteration=50000, CBS3D_RESTART_AS_INITIAL_CONDITION=0, partition root partitions_2. Rank0 original .par SHA256 706b64a8e641a53294c393c23dba3c086cd104fd2961541af7cfe3d01839e6ac; matprop SHA256 928bc80463893a6d80c19eb9b5df27e5553ccd8b9821473935acf7e12641a3fb; binary SHA256 96ca7134f6affc2cd17cdfb6b14ba70a00bc1c16fa.
- Source code confirms the legacy VTU import copies the saved velocity, pressure, and temperature into the solver, initializes iteration histories to those values and then advances at global restart+1. First-step VTU can be recorded as 50001.
- A **not-yet-executed** G6 Slurm job is committed: jobs/sunbird/templates/replay_uchannel_first_step_g6.slurm. It requires all three original SHA256s, builds a NEW master par containing only 1 local step plus output every step, stages a private copy using run_mpi.slurm, imports the same original AFC-on step50000 VTU, runs same production binary with AFC disabled and global dt, then audits steps 50000 and 50001 with tools/audit_temperature_undershoot_vtu.py.
- The checker identifies the first new sub-reference node, its global ID, location, velocity and material adjacency. This is the necessary entry point before numerical term decomposition. The first-step velocity after Step 3 may differ from the restart snapshot; avoid attributing a particular stencil without the corrected velocity and dt analysis.
- **No production numerical source or original case inputs altered**; G6 replay must be run and output inspected before declaring exact cause.

## G6 original U-channel first-step replay — observed Job 29763 (2026-11-02)

- The one-step replay used the original, hash-pinned U-channel binary 96ca7134f6affc2cd17cdfb6b14ba70a00bc1c16fa, original Job 29706 rank0 input .par SHA 706b64a8e641a53294c393c23dba3c086cd104fd2961541af7cfe3d01839e6ac and .matprop SHA 928bc80463893a6d80c19eb9b5df27e5553ccd8b9821473935acf7e12641a3fb.
- Restart imported from original AFC-on Job 29580 at step 50000; no production code changes; private par changes only ntime=1 and output/monitor cadence. Global timestep (ilots=0) in first iteration = 0.000458758 s.
- Solver completed first local iteration 1/1, printing **Tmin=298.051 K**, Tmax=398.241 K, maximum velocity=0.232705 m/s, RelT=8.94382e-6, velocity correction max=0.00337992 m/s, DivRMS=3.62957e-9, negligible velocity residuals. This reproduces the previously reported immediate AFC-disabled temperature undershoot relative to 298.15 K (~0.099 K).
- Provided snippet shows the Slurm job still RUNNING at 4 seconds, while solver's srun step 29763.0 COMPLETED 0:0; postprocessor had printed '=== Cold-node comparison ===' but not completed its output. Slurm sacct encountered a network connection error to sb-wlm-01; no basis to interpret this as solver numerical failure. **Do not mark full runner completed until final postprocessing and status are available.**
- Next: inspect G6 job batch output after completion for first cold global node ID, coordinates, material and boundary touches, compare step 50000 and 50001. Then determine the controlling thermal residual and signed element contributions using actual corrected first-step velocity and global dt. Existing VTU output may expose only the end-of-step velocity, so avoid silently equating restart velocity with the corrected field.

## G6 observed node-level first-step audit (Job 29763; provided output)

- Existing two-rank VTU auditor successfully inspected both initial global step 50000 and first resumed AFC-disabled step 50001 and printed `G6 FIRST-STEP REPLAY COMPLETE`. The only reported error was a separate Slurm accounting-service networking failure (sacct to sb-wlm-01); numerical solver and auditor output completed. No final sacct status was provided.
- Unique owned-node initial Tmin/Tmax=298.150000/398.241306 K; first-step Tmin/Tmax=298.051487/398.241349 K. 94 nodes are below reference after 1 step, of which 21 are below 298.140000K. These 21 are *all* fluid-only domain-kind 1, with no physical boundary incidence; velocities range 0.03715839 to 0.20936081 m/s.
- Global coldest ID=5707 at x=0.04286716590163, y=0.03644056582626, z=0.004375m, T=298.0514872962 K; post-step velocity=(0.189359135962,-0.04367619653701,0.01826773500022)m/s, speed=0.1951876088747m/s. Node belongs to fluid material 0 only, no external boundary, no pressure fixation and no prescribed velocity BC.
- The other nine coldest IDs are 5708,5728,5733,5819,5800,5764,5634,5717,5825. Because 94 includes sub-0.01 K drops, **do not claim all 94 nodes are boundary-free**; classification was run on the threshold-selected 21.
- Read-only independent exact first-step residual reconstruction tool now committed at `tools/audit_uchannel_first_step_residual.py` (NOT YET EXECUTED). It combines T at 50000, corrected u at 50001, unique tetrahedra adjacent to node 5707 and exact CSV dt. Decomposes nodal thermal update into Galerkin advective, discrete div(u) conservative correction, characteristic CBS, thermal diffusion and fluid volumetric source; checks production RHS/diagonal-consistent update against measured first-step T difference to within 1e-6 K.
- If the reconstruction mismatches, inspect source-property values, config and post-step velocity consistency; **do not infer a particular numerical defect** until matched. If it matches, isolate which signed term dominates and confirm negative off-diagonal mechanism with the original stencil. No production numerics changed.

## G6 exact first-step source decomposition — verified output (2026-11-02)

- User ran `python3 tools/audit_uchannel_first_step_residual.py` on the actual Job 29763 replay, using original Job 29706 material properties, case par, step50000 initial temperature and corrected step50001 velocity. Target: global node 5707, interior fluid-only, no thermal boundary faces.
- Exact global one-step dt from CSV 4.5875770122873882e-4 s, inlet/reference temperature 298.15 K; target initial T 298.149999999982 K, post-step T 298.051487296235 K; observed ΔT = -0.09851270374679 K.
- Node incident to **22 fluid tetrahedra**, production thermal diagonal independently reconstructed H=30.00608461324 W/K; maximum incident-element |div(u)|=57.63623879709 s^-1.
- Independently reconstructed per-term nodal ΔT: Galerkin advection **-0.09983648768179 K**, relative-enthalpy divergence correction **-0.004152433214213 K**, CBS characteristic term **+0.003669035464927 K**, conduction **+0.001807181684282 K**, volumetric source 0. Net -0.09851270374680 K. Observed minus reconstructed = +3.302913498260e-15 K; PASS.
- Interpretation: the Galerkin transport term dominates the first negative increment. Divergence correction makes it more negative at the target; characteristic stabilization and diffusion are positive yet insufficient to preserve its reference lower bound in this actual one-step case. **Do not conflate dominating signed residual contribution with proof of a negative off-diagonal DMP coefficient**; actual-mesh stencil reconstruction is next.
- `tools/audit_uchannel_first_step_residual.py` now includes an optional in-place printout of original-step first-row matrix/stencil coefficients (B_ij = δ_ij + (1/H)dR_i/dT_j), and identifies signed negative-link weighted contributions. These checks are committed but **have not been executed on Sunbird2 yet**. Their goal is to establish the exact original-node DMP mechanism without numerical production changes.
- Caution: a large P1 elementwise divergence is not an assessment of global mass-conservation error by itself. It warrants separate velocity/continuity investigation; do not attribute the temperature issue to it without evidence.

## G6 actual-mesh stencil counterexample (Sunbird2 output, 2026-11-02)

- An independent reconstruction of the exact first-step explicit thermal matrix row on Job 29763 **PASSED**: for node 5707, stencil contains 14 nodes (13 neighbours), with 5 negative and 8 positive off-diagonal B_ij. Diagonal B_ii=0.9982969979740.
- Weighted contributions sum negative-link B_ij*(Tj-Tref)=-0.3102716641480 K, positive-link sum=+0.2117589604012 K, diagonal -1.776169030824e-11 K, yielding T_after-Tref=-0.09851270376459 K; direct saved field=-0.09851270376458 K. Difference -2.983724378680e-15 K.
- Dominant negative links: gid5621 B=-0.008027575798, theta_j=15.48508836K, effect=-0.1243077205K; gid1134 B=-0.003504110872, theta=30.66024604K, effect=-0.1074369015K; gid5769 B=-0.006957484972, theta=6.953538381K, effect=-0.04837913879K; gid5657 B=-0.0008799167342, theta=34.26222298K, effect=-0.03014790335K; gid5708 B=-0.004509991945 but theta ~0 K so negligible first-step effect.
- In these negative coefficients, the separately reconstructed Galerkin advective term is predominantly negative. Thus a local **positivity/DMP violation of the unlimited explicit CBS thermal advection operator** is directly verified on the production U-channel field, not merely on a synthetic mesh.
- Scope: this demonstrates the instantaneous cause of the first-step 0.098513K undershoot. Whether it explains the eventual 62.7K cold excursion at step 100000 needs further convergence/long-run evidence. The observed large local P1 div(u) remains a separate continuity accuracy question. No production solver edits made.
- Next G7: verify one-step integrated global thermal energy and inlet/outlet mass flux using actual fields; a DMP failure need not violate conservative global balance. Be careful: legacy tools/energy_balance.py is designed specifically for ports on a common plane (blanket slice y=0) and a steady state, so it must **not** be applied blindly to the U-channel transient/pseudo-time step.

## G7 conservation experiment — prepared, not yet executed

- Script `tools/audit_uchannel_global_thermal_balance_g7.py` works from the **same saved** two-rank VTUs as Job 29763, steps 50000/50001, with original `.par` and `.matprop`. No solver edits, no re-running long transient.
- Independently assembles P1-tetra node heat capacities C_i=rhoCp V/4 and Step-4 thermal contributions, and integrates all external BC532 heat fluxes, material volumetric heat sources and **exact P1 normal-temperature flux** across the closed fluid-region boundary. Fluid/solid interface flux is reported separately (ideally negligible due to no-slip).
- For the pseudo-time update, the identity checked is `sum C_i*(T1_i-T0_i)/dt = Q532+Qvol-Hadv+QDirichlet`, where `QDirichlet = sum_(fixed T nodes)[C_i*DeltaT_i/dt - assembled R_i]` accounts for temperature nodes overwritten by strong Dirichlet constraints after thermal advancement. Ignoring this exchange could cause an apparent false energy imbalance.
- Also checks equality of finite-element volume-integrated conservative advection against exact fluid surface enthalpy integral, and worst unconstrained node update against its assembled residual. Convergence tolerance for absolute integrated W mismatch is 1e-5 W; **cannot be marked passed until Sunbird2 reports actual output**.
- This verifies discrete conservation for one explicit *pseudo-time* step and separately reports net mass boundary flux/volume velocity divergence. It does not establish physical transient energy conservation, global steady heat transfer or ultimate outlet enthalpy at 100,000 iterations.
- The older `tools/energy_balance.py` assumes all inlet/outlet ports lie on a common plane (blanket y=0) and is therefore not a valid default tool for arbitrary U-channel 3D inlet/outlet geometry; do not invoke it blindly.

## G7 provisional audit output (2026-11-02): not a valid solver failure verdict

- User executed `tools/audit_uchannel_global_thermal_balance_g7.py` on Job 29763 restart VTU step50000 and first postrestart step50001 (dt=4.58757701228739e-4 s). Parsed 9613 unique nodes, 45567 tetrahedra (4880 fluid/40687 solid); BC counts 511:44, 520:44, 530:8124, 532:648; 31 fixed-T nodes.
- Read-only audit reports 60.000000 W BC532 source, 0.0 volumetric source, 74.84550869219 W outward relative enthalpy flux, 0 fluid-solid interface advective enthalpy, -0.2939321204950 W Dirichlet reaction, actual measured field storage/dt +1.106171546681 W, predicted by these reconstructed terms -15.13944081269 W; apparent difference +16.24561235937 W (21.7% of largest flux).
- Fluid net mass flux=-6.289028047455e-6 kg/s; fluid integrated P1 div(u)=-6.307951903164e-9 m3/s (rho~997): same independent check. Inlet mdot=-1.250000e-3 kg/s; outlet=+1.243710971953e-3 kg/s (approx 0.503% imbalance relative to inlet).
- **Crucially, largest free-node update mismatch** `C_i*(T1-T0)/dt - R_i` is **0.01178141406118 W**, which is orders of magnitude beyond numerical roundoff; thus the independent global residual calculation does not yet reproduce all production nodal updates. Do not label the 16.2456 W discrepancy a proven solver conservation failure. On the other hand, separately computed volume/face advective enthalpy fluxes match to 2.84217e-13 W, so the checker geometry/divergence theorem is consistent with its own reconstructed field.
- G6 reconstruction at node5707 had matched to 3e-15 K, but that validates only this interior fluid node; it cannot establish that solid, interface, turbulent-thermal or constrained updates are captured throughout the mesh.
- A follow-up DIAGNOSTIC revision of the same read-only script prints per-node residual mismatch grouped by fluid-only, solid-only, interface and Dirichlet node, with the 20 largest free-node discrepancies and BC incidences; pending user execution. Do not loosen acceptance tolerance or change production FEM yet. Distinguish possible postprocessing field/coefficients differences from real MPI assembly/update issues based on node classifications.

## G7 audit correction: optional nuclear-heating profile (2026-11-02, execution pending)

- Localization of free-node residual discrepancies showed fluid-only 273 nodes agree (sum +2.351186399081e-11 W; max 5.6046e-12 W), while solid-only 8110 nodes had sum +14.98599373775 W and interface 1199 nodes had sum +1.259618621597 W; fixed-node sums +0.1768666808623 W (fluid) and -0.4707988013573 W (interface). Together these account for the global 16.24561235937 W apparent residual gap.
- Source inspection revealed `src/preprocess/NuclearHeating.cpp::applyProfiles`, which writes the optional `nuclear_heating_on` exponential per-material volumetric heating source into `Qvol_e` at tetrahedron **centroids** before energy assembly. Its profile is `A exp(-max(0, ((x-x0)·n)/|n|)/lambda)`. `src/io/MeshIO.cpp` reads these profiles from the optional `.par` tail. The original G7 audit used **only constant .matprop volumetric Q**, so it omitted this possible source from both aggregate energy and nodal RHS.
- A 16.245612 W omitted source would explain the observed global discrepancy, and its predominant solid-only/interface distribution strongly supports this hypothesis. **This is not yet proven:** need inspect actual U-channel par profile enabled and run the corrected checker; do not assert energy conservation PASS without independent updated output.
- Updated read-only `tools/audit_uchannel_global_thermal_balance_g7.py` parses optional nuclear-heating profiles, evaluates exact same element-centroid law, rejects conflicting constant and profile sources, includes source in nodal residuals/storage identity, and reports volumetric source by material. No numerical solver source modified.
- Run checker on Job29763 saved VTU after `git pull --ff-only`, preserving all original case files. Required evidence: nuclear heating per material, worst free-node residual, whole-domain storage closure, and surface/volume transport mismatch, preferably all <= configured 1e-5W acceptance in closure and local residual checks.

## G7 first-step GLOBAL energy conservation VERIFIED (2026-11-02)

- Verified by direct user execution on Job 29763 restart steps 50000 to 50001; actual audit output after fixing optional nuclear-heating source. .par nuclear heating = 1, origin=(0,0.05,0), normal=(0,-1,0), one profile on solid material 1 with A=1.6e6 W/m3 and decay length lambda=0.02m. Source integrated using exact production element-centroid law.
- Mesh: 9613 unique nodes / 45567 tetrahedra (4880 fluid / 40687 solid), 31 fixed-temperature nodes, external faces BC511=44, BC520=44, BC530=8124, BC532=648; pseudo-dt = 4.58757701228739e-4 seconds.
- External BC532 Q=+60W, nuclear volumetric Q=+16.24625521637W (all solid), outward coolant relative enthalpy=+74.84550869219W, fixed-temperature exchange=-0.2945749775253W; observed storage/dt=+1.106171546681W, independently predicted storage/dt=+1.106171546657W. **Global closure=+2.374322960463e-11 W, relative 3.172298514568e-13: PASS.** Fluid-solid normal advective enthalpy 0; volume/surface advection mismatch=-2.131628207280e-13W.
- Free-node max |C ΔT/dt - R|=5.604600117337e-12 W. Per class free-node sum error: fluid-only +2.351186399081e-11W (273 nodes), interface +6.676151945929e-12W (1199 nodes), solid-only -6.237777267416e-12W (8110 nodes). All demonstrate exact local update consistency within roundoff. Fixed-node reaction intentionally nonzero due to imposed Dirichlet BC and is accounted for in balance.
- Net fluid mass flux=-6.289028047455e-6 kg/s with inlet magnitude 1.25e-3 kg/s, giving **0.5031% inlet mass imbalance**. Integrated divu=-6.307951903164e-9 m3/s. Mass balance has NOT been asserted PASS at a strict criterion; separate G7b continuation/continuity study required. At this pseudo-time step, outlet heat flux > wall heating because stored thermal energy from existing initial field and volumetric heating are present; this is not a steady thermodynamic balance.
- Crucial: DMP failure (G6) and exact conservative global energy transport (G7a) **coexist**. The initial 16.2456W apparent mismatch was an audit-source omission (nuclear heating), NOT a production conservation defect. No production numerical code modified.

## G8a native-versus-legacy two-rank restart parity (prepared; unexecuted)

- G6 replay automatically writes a native checkpoint at **step 50000** on initial legacy VTU import, beneath `$CBS3D_SCRATCH/verification_g6/29763/run_uchannel_29763/output/uchannel/restart/step_00050000`. This was confirmed from the production RestartIO code; the target path must exist and is checked before the job proceeds.
- Audit runner `jobs/sunbird/templates/verify_g8_native_restart_parity.slurm` replays exactly one additional pseudo-iteration from the native checkpoint, using the **same original 2-rank partition, production binary SHA256, and copied one-step .par/material inputs as G6**. No other source or original case modifications. It checks all original SHA256 input files before running, and compares original G6 .par against the recorded staged one.
- Parity checker `tools/compare_cbs_restart_parity_g8.py` compares **every global owned-node** pressure, velocity and temperature at steps 50000 and 50001 against the legacy-VTU G6 replay, by global node ID, with coordinate identity. T toler=1e-9K, velocity tol=1e-10m/s, pressure tol=1e-7Pa, coordinate tol=1e-12m. Reports max errors and offending global IDs. JOB NOT YET EXECUTED.
- Scope: verifies native-checkpoint and legacy-VTU restart path equivalence on the same MPI partition and original solver binary, not serial/MPI or alternative MPI rank-count invariance. Continue G8b after G8a PASS; do not mark whole G8 PASS from G8a alone.

## G8a verified result — Job 29786 (2026-11-02)

- Slurm job 29786 COMPLETED 0:0 (batch, extern and solver step all 0:0) in 3 seconds. Original binary .par and .matprop SHA256 checks match Job 29706. Same original 2-rank partition, one restart continuation step 50000 -> 50001. Original G6 legacy-VTU run Job 29763 compared against native-checkpoint restart Job 29786.
- Independent owned global-node-ID comparison of all 9613 nodes at BOTH steps: maximum absolute temperature difference = 0 K; velocity difference = 0 m/s; pressure difference = 0 Pa; coordinates difference = 0 m. G8A NATIVE / LEGACY RESTART PARITY PASS. This checks two restart formats at fixed partitioning, not rank-count independence.
- `sacct: ... Network is unreachable` appeared on stderr after completed result from the `sacct` command embedded at end of `jobs/sunbird/templates/run_mpi.slurm`; it concerns Slurm accounting connectivity to sb-wlm-01, port 6819, **not** CBS3D solver or MPI. Login-node sacct successfully reported all jobs COMPLETED, exit 0:0. The staging script intentionally tolerates this accounting failure using `|| true`. No rerun required.
- Next G8b: examine rank-count / serial consistency with explicitly reconstructed partition/restart controls, preserving exact solver binary and case, without changing numerical production code.

## G8b independent MPI rank-count 2-versus-4 first-step audit (prepared; not yet executed)

- G8b test files: `tools/verify_uchannel_partition_parity_g8b.py`, `jobs/sunbird/templates/verify_g8_partition_count_2v4.slurm`. Exact same original production binary SHA256 96ca7134..., .par SHA256 706b64..., .matprop SHA256 928bc804... as first-step Job29763, no production source edits.
- Job independently exports the **same full 45567 TET4, 9613-node, 8860-physical-triangle serial U-channel mesh** to 4 ranks using existing Gmsh/partition exporter. The job discovers the original unpartitioned mesh by exact header and matching .bco/.material; fails safely if zero or multiple candidates. Optional `CBS3D_G8B_SERIAL_BASE=/path/to/original/uchannel` overrides discovery.
- G8b remaps 2-rank G6 step50000 VTU global owned fields T/u/p to 4-rank rank-local node ordering from each new .mpi NODES mapping, with 17-digit precision, **no interpolation**. RestartIO legacy-VTU importer accepts the minimal VTK XML fields, validates expected rank-local global_node_id ordering, and sets same completed iteration=50000, so 4-rank one-step output is 50001. Original 2-rank first-step Job29763 is comparison reference.
- Independent comparator reads actual two-rank and four-rank solver VTU at steps 50000 and 50001, verifies same 9613 globally owned IDs, and reports maximum per-field differences and worst ID. PREDECLARED maximum absolute tolerances: T 1e-9K, velocity 1e-10m/s, pressure 1e-7Pa, coordinate 1e-12m. Any discrepancy above tolerance fails and must be diagnosed rather than silently relaxed. Test requires exactly 2 vs 4 piece counts.
- This tests MPI decomposition sensitivity on the identical physically hot restart field and actual coupled CBS time step; G8b is not proof of serial/MPI parity or longer-run convergence. Independent G8c serial consistency remains open. Unexecuted awaiting Sunbird2 Slurm job results.


G8b Job 29788: Slurm COMPLETED 0:0, 6s; independent four-rank partition, identical step-50000 global field; step-50001 maximum temperature difference 0 K, velocity component difference 1.040834e-17 m/s, pressure difference 0 Pa. Original binary SHA256 matched. G8b first-step parity PASS; multi-step and serial remain open.
## Directed development pivot after G8c (2026-11-02)

- Job 29789 COMPLETED 0:0. Same full U-channel mesh and identical initial step50000 global fields, 2 versus 4 MPI ranks, 100 coupled iterations through 50100. At step50100, across all 9613 global nodes max abs deltaT=7.332801033044e-12 K (gid909), max abs velocity component=2.519685848856e-13 m/s (gid5617), max abs p=1.883094374877e-13 Pa (gid1326), max abs coordinate=0 m. All predeclared tolerances PASS. This is not serial/MPI parity.
- **Correct conceptual partition:** TemperatureAFC.cpp limits Step-4 temperature convection/diffusion flux and is not in momentum/continuity equation. Pressure/velocity CBS update controls mass conservation; G7 showed relative inlet/outlet mass-flow deficit 0.503% at step50001. AFC cannot be credited with correcting that. Thermal first-step global energy conservation PASS independent of AFC on/off; the checked G7a test used AFC OFF.
- **Concrete code-level bug:** CMake builds `src/solver/SolverDistributedProductionLoopRestartAfc.cpp`, which defines `EnergyAssembly` macro as `RestartAfcEnergyAssembly` and includes `SolverDistributedProductionLoopRestart.cpp`, which in turn includes `SolverDistributedProductionLoop.cpp`. Thus the wrapper's `RestartAfcEnergyAssembly::assembleStep4Rhs` already called `TemperatureAFC::limitDistributedResidual(s)` and the main explicit loop called it a second time. A second nonlinear graph-viscosity limiter pass acts on already-limited owner residual; this differs from one-call AFC residual used in the PETSc steady solver, jeopardizing formal boundedness and consistency.
- Minimal correction committed on audit branch: **removed the AFC invocation inside wrapper**, leaving the established main production Step4 call as the single owner after CHT assembly; no other numerical code changed. **NOT YET COMPILED OR VALIDATED ON SUNBIRD2.** Existing original executable remains intact, its SHA pinned.
- G11a test `jobs/sunbird/templates/verify_afc_one_call_uchannel_g11.slurm` builds modified source in isolated scratch, runs exactly one U-channel AFC-on step from original step50000 restart and checks runtime AFC invocation count=1, all 9613 owned-node initial fields equal, final pressure/velocity match AFC-off baseline, temperature lower bound >=298.15K. Independent comparator `tools/audit_afc_uchannel_first_step_g11.py`. If run fails, keep outcome visible; do not automatically change AFC settings or assert its effectiveness.
- **Existing implicit implementation:** `src/solver/SteadyEnergySolver.cpp` solves `-M T=c+q` at fixed velocity with PETSc GMRES+ASM/ILU(1) by default. With AFC=1 it solves the *single-call* nonlinear bounded residual via defect correction and line search, fixed low-order preconditioning matrix, optional Anderson acceleration. `CBS3D_STEADY_ENERGY=1` enters solver instead of coupled iterations. This is a **steady frozen-velocity CHT solve**, not a fully coupled implicit pressure/momentum/thermal nonlinear system.
- G12a benchmark `jobs/sunbird/templates/benchmark_implicit_cht_uchannel_g12.slurm`: isolated MPI rebuild, same hot U-channel step50000 restart, runs both AFC=0 and AFC=1 implicit steady frozen-flow PETSc solves, records KSP counts, nonlinear defects, T range and PETSc solver seconds. **UNEXECUTED.** Must assess convergence, conservation, boundedness and acceleration vs valid equivalent explicit steady criterion. Do not claim it is faster until those metrics are measured.


## G11a/G12a Sunbird2 job results (2026-11-02)
Job29832 G11a COMPLETED 0:0. Isolated corrected binary SHA256 6a565fb90209b0f8aa15a5b3c8ab4919cfff2f8e79d10d9a70321ae5c7de0f0d. AFC called exactly once: 20853 edges, 7460 limited, alpha_min=0, Tlow_min=298.15K, Tafc_min=298.15K, Tafc_max=398.072K; flux correction balance 3.23297e-13. All 9613 owned global nodes identical initial T,velocity,pressure,coordinates. At first thermal iteration AFC OFF Tmin=298.051487296235 K, 94 cold nodes; AFC ON Tmin=298.149999999976 K, zero nodes below 298.15K at defined comparison tolerance; node5707 AFC OFF=298.051487296235K vs AFC ON=298.149999999982K. Velocity/pressure unchanged at step50001. G11a bounds and one-call PASS, not long-run convergence or global bounded-CHT accuracy.
Job29833 G12a FAILED 1:0 because AFC implicit variant failed (Slurm step .1 cancelled). Unlimited Galerkin implicit variant SOLVED: 94 GMRES iterations; linear relative residual 2.74207e-10; reconstructed unlimited explicit residual relative 9.77124e-09; temperature range 233.482..444.333K (below physical 298.15K); reported frozen-flow thermal solve 0.0490794s (not whole-job runtime). AFL nonlinear solver failure reason not available in submitted outer job log: request actual benchmark_implicit_cht_g12/29833/afc/run_uchannel_29833/solver.err and ending solver.out. Do NOT claim AFC implicit validated or that the 0.049s unlimited algebraic solve is a usable CHT solution. No new production numerics changes pending diagnosis.

## G12a post-mortem and G12b controlled parameter experiment (2026-11-02)
Job 29833 implicit AFC solver.err: 'Implicit steady energy did not converge: AFC iterations=150, last max|dT|=5.42467e-05 K, nonlinear relative residual=4.75107e-05, total Krylov iterations=450'. Both MPI ranks report same controlled solver nonconvergence. The subsequent Slurm SIGKILL concerns termination of the failed step, not a demonstrated out-of-memory failure.
At iteration 1 PETSc defect-correction trial alpha=1 achieved L2 pseudo-time residual 0.00266321 from initial 0.00266331; explicit fallback alpha=1 achieved 0.00266306, so the line search correctly preferred explicit fallback. In 150 iterations all reported methods remained 'explicit fallback', 3 KSP iterations each; nonlinear relative raw residual fell only from 4.81992e-05 to 4.75107e-05. Both PETSc and explicit directions were finite, satisfying Armijo and max-rate guard on iteration 1. Thus this is **nonlinear convergence-rate failure**, not PETSc KSP divergence or a line-search rejection at iteration 1.
Working hypothesis: the fixed pseudo-time shift gamma=1 used as A_low + gamma C/dt makes the low-order preconditioning nearly the explicit time step at very small dt=4.5876e-4 s, preventing useful large stable stationary corrections. This must be tested rather than asserted. G12b job 'jobs/sunbird/templates/diagnose_implicit_afc_shift_g12b.slurm' uses the **same hash-pinned corrected executable from Job29833** (6a565fb9...), same U-channel step50000 restart, no numerical source changes. Five arms: (gamma=1, Anderson=3), (gamma=0.01, Anderson=0), (gamma=0.0001, Anderson=0), (gamma=0, Anderson=0), (gamma=0.01, Anderson=3). Identical 150 iterations, AFC relax=0.5, update tolerance=1e-6 K, nonlinear raw residual tolerance=1e-6, same PETSc GMRES/ASM/ILU configuration. Record convergence, relative residual and absolute temperatures for every successful arm; mark failed arms honestly. This is parameter diagnosis, **not** an automatic justification to relax tolerance or alter scheme.
Any successful implicit AFC candidate still needs physical Tmin bound, total energy, accuracy and wall/outlet quantities checked. Unbounded Galerkin implicit result cannot be accepted simply because its linear solver was fast.

## G12b Job 29843 partial sweep diagnosis (2026-11-02)
Job29843 batch COMPLETED 0:0, one srun step FAILED 1:0. Output shows ONLY shift1_anderson3 arm was run; despite misleading '0 / 5' message the other four gamma/Anderson configurations did not execute. Measured first arm gamma=1, Anderson m=3: after 150 AFC iterations and 450 KSP iterations, nonlinear relative residual 4.55944e-5 (required 1e-6), no converged steady thermal field. This is not evidence that gamma=0.01, 0.0001, or 0 failed.
Found shell orchestration defect: while-read loop fed through stdin here-string, and nested Slurm srun can consume remaining stdin, skipping other arms. Corrected verification runner to iterate over an in-memory bash array (no stdin stream to consume), explicitly redirect launcher stdin from /dev/null, count actual executed arms, and return nonzero status if no arm converges with physical temperature floor. No production C++ numerical changes; all original configurations retained. New G12b job required for completed sweep.

## G12b complete sensitivity and G12c algorithmic experiment (2026-11-02)
Job29844 completed all five shift/Anderson arms; none converged to nonlinear tolerance 1e-6. gamma=1/m3 final raw relative residual 4.55944e-5 (150 nonlinear iterations); gamma=0.01,0.0001,0 with m=0 all final 4.75107e-5 and selected explicit fallback at every logged iteration; gamma=0.01/m3 stalled with Armijo failure near raw relative 4.6e-5. Lower gamma increased KSP work (gamma0 49 vs gamma0.01 10 Krylov per iteration) without improving nonlinear progress. Parameter tuning exhausted.
Experimental development branch implementation: source src/solver/SteadyEnergySolver.cpp includes optional CBS3D_STEADY_ENERGY_LOW_ORDER_BOOTSTRAP=1. Assemble same AFC low-order matrix, remove pseudo-time shift from a temporary matrix, solve stationary -(M+D)T=c+q with PETSc subject to unchanged Dirichlet rows; use that field as the starting point for the existing exact production AFC nonlinear residual/line search. The original default remains unchanged. Research hypothesis: a global low-order stationary initial guess may circumvent tiny explicit pseudo-time steps without sacrificing the bounded nonlinear target. No convergence or monotonicity claims until tested.
G12c job jobs/sunbird/templates/benchmark_implicit_afc_bootstrap_g12c.slurm builds candidate in scratch, runs gamma0 Anderson0/3 with identical original U-channel hot restart, and verifies solver status and independent final Tmin>=298.15K on 9613 owned nodes. This is NOT a production release; test may fail and requires CHT energy/accuracy validation after any success.

## G12c Job 29852 outcome and G12d solver baseline (2026-11-02)
G12c COMPLETED as failed Slurm Job29852 exit1:0. Corrected candidate rebuilt (SHA256 1c293135aaf1d866663ea9d2fb1a97b386f5e728265cc16f025b3374bf8c0431). A steady low-order PETSc solve with 49 KSP iterations produced bounded Tmin=298.15K, Tmax=375.849K. However, subsequent nonlinear AFC solver for both m=0 and m=3 FAILED: no PETSc/explicit line-search direction reduced L2 scaled residual. At stall m0 L2=0.238511K and raw relative residual=0.048046; m3 L2=0.172194K and raw relative residual=0.045432. Changing initialization/Anderson did not resolve the nonlinear issue. Do not claim that an implicit AFC solution exists or is accurate.
Engineering pivot: Rather than repeatedly tuning the non-smooth explicit-timestep AFC stationary fixed-point map, isolate a reliable stationary graph-viscosity finite-element baseline, then design a separately verified, conservative, high-resolution bounded *steady* transport correction (a mathematically distinct scheme from the explicit Zalesak limiter). CBS3D_STEADY_ENERGY_LOW_ORDER_ONLY=1 is a new **opt-in experimental** mode requiring CBS3D_TEMPERATURE_AFC=1 and CBS3D_STEADY_ENERGY=1. It solves A_low T=(c+q) by one PETSc GMRES+ASM/ILU call and checks ||b-A_low T||_inf/||b||_inf. It bypasses the untrusted nonlinear steady AFC corrections. Not a production method: graph diffusion can materially degrade gradients/Tmax; no energy-budget or accuracy validation yet. Existing production modes and original binary unchanged.
G12d isolated compilation+trial job: jobs/sunbird/templates/benchmark_implicit_low_order_g12d.slurm; uses original 2-rank U-channel hot restart and SHA256-checked par/matprop; checks assembled low-order residual and independent VTU Tmin>=298.15K across 9613 nodes. **NOT EXECUTED**. Success is only low-order bounded baseline readiness, not high-order CHT success.

## G12d completed standalone low-order implicit PETSc U-channel baseline (2026-11-02)
Job29859 COMPLETED 0:0; MPI solver Slurm step 29859.0 COMPLETED in 2 sec. Candidate scratch binary SHA256 7902003bada55429b3541d00e36ec7c33fd779cace6e6814fbfe934e93e25c67. Source commit 6ebf7822b08db7f6d0a20e7a8c225b0d802adde4.
With CBS3D_STEADY_ENERGY_LOW_ORDER_ONLY=1, AFC switch=1 solely to enable fluid-element graph viscosity, no nonlinear AFC corrections: GMRES 49 iterations; independently checked ||b-A_low*T||inf/||b||inf=5.58048e-11, accepted; solver time 0.0356358 seconds (not whole job); across 9613 owned global nodes Tmin=298.15K and Tmax=375.849353817137K, so first case benchmark bound PASS. Result physically bounded **for this case**, not a theorem for all meshes, because solid elements are not graph-limited by default. Nonlinear AFC code has NOT become usable.
Next development: implement and benchmark a conservative high-resolution *steady* algebraic flux correction/monolithic convex limiting formulation with consistent treatment of the physical volumetric nuclear heat source and BC532 surface heat load. A posteriori clipping of temperatures would destroy conservation and should not be used. Keep graph-viscosity baseline as an opt-in reference, quantify changes in Tmax, wall gradients, CHT interface flux, outlet enthalpy and integrated energy closure before accepting a new scheme. Research basis includes Knobloch, Kuzmin and Jha, JCP 518 (2024), article 113305, DOI 10.1016/j.jcp.2024.113305; source-compatible, well-balanced stationary convex limiting rather than time-dependent Zalesak FCT repeated to stationarity.
Because new source commits 742a17d and 4560130 change labels only, Job29859 executable hash and raw measured results remain associated with its original commit 6ebf782. Rebuilding would yield a different binary but should not alter numerical algorithm.