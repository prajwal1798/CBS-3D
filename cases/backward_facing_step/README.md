# Backward-facing step, turbulent, Re = 3025 (validation of Spalart–Allmaras)

Reproduces Nithiarasu & Liu, *CMAME* 195 (2006) 2961–2982, Sec. 8.2 / Fig. 4(b),
and C.-B. Liu (2005) thesis Sec. 8.1 (3D: the 2D section extruded by unit span).
Experiment: Denham, Briard & Patrick, *J. Phys. E* 8 (1975) 681–683.

## Set-up (non-dimensional, step height h = 1)

| Item | Value | Source |
|---|---|---|
| Domain | inlet at x = 0 (4h upstream of the step), step at x = 4, outlet at x = 40; channel 2h upstream, 3h downstream; span z = 0–1 | Liu Fig. 7.6 / paper Fig. 2 |
| Re | U₀h/ν = 3025, U₀ = mean upstream velocity (Denham: 202 mm/s, h = 15 mm, water) | Denham |
| Fluid | ρ = 1, μ = 1/3025 = 3.3057851e-4 (`bfs3d_liu.matprop`) | |
| Inlet (507) | measured u(y) of Denham's first station, `reference/inlet_denham_station0.txt` (mean 0.958), via `CBS3D_BFS_INLET_PROFILE`; v = w = 0 | paper: "inlet profile obtained from Denham et al." |
| SA inlet | ν̃ = 0.05 → `sa_inlet_ratio` = 0.05 × 3025 = 151.25 | paper Sec. 8.2 |
| SA model | `turbulence_model 0` (Liu–Nithiarasu SA-noft2) | |
| Walls (530) | no-slip, ν̃ = 0 | |
| Outlet (520) | p = 0 | |
| Side faces (506) | symmetry (2D-equivalent; Liu: 2D and 3D results identical) | |
| Solve | steady, local time stepping (`ilots 1`), stop at RelU < 1e-7 and SA residual < 1e-6 | |

Without `CBS3D_BFS_INLET_PROFILE`, BC 507 is the legacy CBS2D++_SI polynomial
(peak 1.38 at y = 1.77, mean 1.03). That is **not** Denham's profile; set the
variable for this case.

## Reference data (`reference/`)

Extracted exactly from the vector paths of the paper's Fig. 4(b) by
`tools/digitise_liu2006_bfs.py` (no pixel reading):

| File | Content |
|---|---|
| `liu2006_fig4b_sa.csv` | Liu's SA solution, stations x/h = 0, 1, 2, 4, 6, 8 |
| `liu2006_fig4b_exp.csv` | Denham's data as plotted by Liu |
| `inlet_denham_station0.txt` | BC 507 table: Denham's first station, u = 0 at y = 1 and 3 |

Check: ∫u dy of Liu's SA curves is 1.932–1.955 at all six stations (±0.6 %),
as mass conservation requires.

Stations: Denham measured at x = −20, 0, 15, 30, 60, 90, 120 mm (h = 15 mm).
Liu's profiles are x/h = 0, 1, 2, 4, 6, 8; his plot offsets 0, 2, 4, 6, 8, 12
are display shifts. Evidence: the deepest reversal (−0.235) is in his third
profile and Denham's is at 30 mm (2h); 15 mm is weaker (−0.12); 60 and 90 mm are
≈ 0 at the wall.

## Mesh

Gmsh with the background size field of `bfs3d_liu.geo`: recirculation box,
step lip, and wall refinement by one Box field per planar wall (lc_wall = 0.03h
at the wall, growing linearly to lc_far over d_wall = 1.0h, a size gradient
of ~0.2).

| `--clscale` | Tetrahedra | First-element height (measured) | Use |
|---|---|---|---|
| 1.75 | 281,457 | 0.056–0.075 h (y⁺ ≈ 4–5) | Liu-sized (297,054) |
| 1.25 | — | ~0.04 h | mesh study |
| 1.0 | — | ~0.03 h (y⁺ ≈ 2) | mesh study |

History (do not repeat): a Gmsh `Distance` field with `Sampling = 200` on the
40h-long walls places its samples 0.2h apart, so it overestimates the
distance and under-refines the walls; a steep gradient (d_wall = 0.3, gradient
0.7) made the first element 2–3.5 times thicker than lc_wall. Job 33994
(338,050 tets) used that earlier version.

Wall refinement is required. The measured inlet profile reaches u = 0.35 at
0.016h from the wall; with the wall nodes at u = 0, a 0.25 node spacing loses
~7 % of the inlet flux (a 2-rank trial without refinement: ∫u dy = 1.55
against 1.917).

## Commands (SUNBIRD2)

```bash
source ~/cbs3d_env.sh
cd $CBS3D_CODE
bash tools/cbs_case.sh cases/backward_facing_step/bfs3d_liu.geo \
    cases/backward_facing_step/bfs3d_liu.par cases/backward_facing_step/bfs3d_liu.matprop \
    $CBS3D_SCRATCH/cases/bfs3d_liu_c175 --parts 16 --clscale 1.75
export CBS3D_BFS_INLET_PROFILE=$CBS3D_CODE/cases/backward_facing_step/reference/inlet_denham_station0.txt
bash jobs/sunbird/templates/submit_auto.sh $CBS3D_SCRATCH/cases/bfs3d_liu_c175 bfs3d_liu --sets 16 --time 12:00:00
```

The solver log must show `BC 507 INLET : TABLE ... (17 points, ...)`.

## Comparison

```bash
python3 $CBS3D_CODE/tools/bfs_profiles.py $CBS3D_SCRATCH/cases/bfs3d_liu_c175 --plot profiles.png --wall-law wall_law.png
```

(The folder argument takes the newest `.pvtu` under it.) Prints, per station,
the RMS and maximum difference from Liu's SA curve, the RMS difference from
Denham's data, the section flux ∫u dy (Liu: 1.93–1.96), and the reattachment
length x_r/h. `--wall-law` adds u⁺ against y⁺ at attached stations (x/h = −2,
12, 20, 30) with u_τ from the exact P1 wall gradient, and prints u_τ, Re_τ and
the first-element y⁺. Liu's SA curves reverse at 2h and are ≈ 0 at 4h, so his x_r is
2h–4h; Denham reports about 6h. `--plot FILE.png` writes the Fig. 4(b) layout.

Verified on a 2-rank `--clscale 2.5` trial (35k tetrahedra): the log shows the
inlet table, ν̃ = 0.05 at every interior inlet node, steady convergence
(RelU 1e-16) after 5,000 iterations, ∫u dy = 1.87–1.90 at all stations.
