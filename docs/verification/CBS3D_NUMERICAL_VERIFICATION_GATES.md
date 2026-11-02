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
| G6 | Discrete positivity / DMP | Inspect assembled update matrix and *first* bound violation | FAIL (unconditional DMP, synthetic meshes, Job 29756): negative offdiagonals and first-step 0.625/1.875/3.750 K undershoots; actual U-channel first-step mechanism still unresolved |
| G7 | Integrated energy and mass | Sources + boundary heat = advection + storage within stated tolerance | Existing balance script; pending for U-channel |
| G8 | Serial/MPI and restart parity | Independent solutions agree at declared tolerances | Pending |
| G9 | Time/iteration convergence | Stationary residual, dt-variation and stopping robustness | Pending |
| G10 | Grid convergence of quantities of interest | Tmax, interface flux, wall gradient and outlet enthalpy | Pending |
| G11 | Additional stabilisation only if needed | Isolated failure, specific literature choice, accuracy + bounds verified | Deferred |
| G12 | Implicit acceleration only after verification | Same discrete result, lower runtime, robust nonlinear residual | Deferred |

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
