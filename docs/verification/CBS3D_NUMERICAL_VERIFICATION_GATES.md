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
| G1c | Lumped thermal mass and timestep diagonal | Code vs direct P1-tet reference with positive entries | Pending |
| G2 | Actual inlet/outlet/wall/interface boundary inventory and flux sign | Every applied BC and flux accounted for | Pending |
| G3 | Solid-only diffusion | Uniform and linear patch tests; analytical heated slab; mesh refinement | Existing slab check available; re-evaluation pending |
| G4 | Two-material conduction | Shared interface temperature, weak flux continuity, integrated balance | Pending |
| G5 | Prescribed-flow scalar transport | Uniform T, analytical convection-diffusion and manufactured solutions | Uniform T check added; remaining pending |
| G6 | Discrete positivity / DMP | Inspect assembled update matrix and *first* bound violation | Pending; never infer a cause from final Tmin alone |
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
cmake --build build/verification --target cbs3d_step4_operator_verification cbs3d_source_validation -j 8
(cd build/verification && ctest --output-on-failure -R '^step4_')
~~~

Expected: the existing source test and the new P1 tetrahedral operator test both pass. The new executable prints five individual PASS lines, then STEP-4 OPERATOR VERIFICATION PASS.

If a test fails, halt and diagnose the specific element contribution. Do not patch production code before an independent reproduction.

## Change acceptance criteria

1. State failing input, governing equation, predicted value and observed value.
2. Identify first affected node/element and local operator terms.
3. Reproduce independently on a smaller problem.
4. Make a single targeted change on a separate branch; record before/after results.
5. Re-run applicable gates, especially serial/MPI and analytical benchmarks.
