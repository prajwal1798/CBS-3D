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
| G1a | Actual Step-4 P1-tet FEM volume and surface terms vs independently evaluated quadrature | All nodal convection, conduction, characteristic, heating, Neumann contributions agree | Added C++ test; not yet executed |
| G1b | Exact frozen-velocity element linear operator versus finite-difference residual derivative | All 16 entries agree | Added C++ test; not yet executed |
| G1c | Lumped thermal mass and timestep diagonal | Code vs direct P1-tet reference with positive entries | PASS: Sunbird2 CTest 3/3 (2026-11-02), shared-node and material masks |
| G2 | Actual inlet/outlet/wall/interface boundary inventory and flux sign | Every applied BC and flux accounted for | PASS (structural): Job 29706 has only BC 511/520/530/532, 648 BC532 faces = 0.0006 m^2, 60 W applied, zero duplicate faces, zero MPI duplicate mismatch; full global energy closure remains G7 |
| G3 | Solid-only diffusion | Uniform/linear patch tests; analytical heated slab; mesh refinement | Analytical slab and serial/MPI PASS, Job 29742; systematic mesh refinement still pending |
| G4 | Two-material conduction | Shared interface temperature, weak flux continuity, integrated balance | Conformal two-layer analytic test and serial/4-rank runner added; execution pending |
| G5 | Prescribed-flow scalar transport | Uniform T, analytical convection-diffusion and manufactured solutions | Uniform T check added; remaining pending |
| G6 | Discrete positivity / DMP | Inspect assembled update matrix and *first* bound violation | U-channel AFC-off run: 176 cold nodes, ALL fluid-only, zero interface/solid; root mechanism and first crossing unresolved |
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
- Exact integrated input Q=qW^2=0.016 W. Nodal max error threshold 0.5% of exact top temperature, zero intended bulk flow, reference area within 0.1%, surface pressure BC inventory and shared interface-node connectivity required.
- Runner: jobs/sunbird/templates/verify_two_layer_g4_no_afc.slurm, AFC disabled and no implicit energy solver. Checker: tools/check_two_layer_conduction.py. This is a new unexecuted verification case; do not claim pass until logs show it.
- Note: interface-near k*dT/dy medians are diagnostic; pointwise element fluxes need not be continuous across a linear FEM interface even when the global weak heat balance is correct.
