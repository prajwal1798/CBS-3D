# CBS3D fully implicit physical backward Euler flow

## Scope

This path advances the existing laminar P1/TET4 CBS discretisation in physical
time. It does not replace CBS by a different mixed finite-element formulation.

The monolithic unknown at each fluid node is

[
  x=[u,v,w,u^*,v^*,w^*,p]^T.
]

For the present constant-property blanket case there is no temperature-to-flow
feedback. Therefore one physical timestep is block triangular:

1. fully implicit BE flow,
2. fully implicit BE energy using the converged (u^{n+1}).

This is exact for the present one-way coupling. It is not a claim of two-way
thermofluid coupling.

## Discrete flow equations

Let (M_L) be the lumped momentum mass matrix, (H) the existing CBS pressure
Laplacian, (B) the Step-2 weak-divergence operator, and (G) the Step-3
pressure-gradient operator. Let (R_m(u)) be exactly the production Step-1
momentum residual from `MomentumAssembly.cpp`.

At Picard iteration (k), freeze the advecting velocity at (w^k) and write

[
  A_u(w^k)u^{n+1} = -R_m(u^{n+1};w^k).
]

The physical backward-Euler CBS system is

[
  A_u(w^k)u^{n+1}
  + rac{M_L}{Delta t}u^{*,n+1}
  = rac{M_L}{Delta t}u^n,
]

[
  Delta t,H p^{n+1} - B u^{*,n+1}
  = b_{Gamma}(u^n),
]

[
  rac{M_L}{Delta t}
  left(u^{n+1}-u^{*,n+1}ight)
  + Gp^{n+1}=0.
]

Adding the first and third equations gives

[
  rac{M_L}{Delta t}
  left(u^{n+1}-u^night)
  + A_u(w^k)u^{n+1}
  + Gp^{n+1}=0,
]

which is the required first-order backward-Euler momentum balance.

The Step-2 boundary flux intentionally uses (u^n), matching
`PressureAssembly::assembleStep2Rhs`.

## Exact Step-1 spatial operator retained

The implicit operator is a Picard linearisation of the production CBS residual,
not a substitute stress formulation.

For a TET4 element, with (V=det J/6), frozen nodal advecting velocity (w_b),
(ar w=rac14sum_b w_b), and
(q_b=
abla N_bcdot w_b), the volume contribution to
(A_u=-R_m) is componentwise

[
 A_{ab} =
 -rac{V}{4}
abla N_acdot w_b
 +rac{Delta t}{2}V
   (ar wcdot
abla N_a)q_b
 +
u V
abla N_acdot
abla N_b.
]

The existing conservative convective face term, CBS characteristic face term,
and componentwise natural-diffusion face term are retained with the same
`fcon`, `fdif`, `annxf`, and face-node conventions as
`MomentumAssembly.cpp`.

## Nonlinear and linear solution

Picard/Oseen iteration is used at every physical timestep. Defaults:

- relative velocity change: (10^{-7}),
- relative pressure change: (10^{-6}),
- divergence RMS: (10^{-7},mathrm{s}^{-1}),
- maximum Picard iterations: 100,
- initial relaxation: 0.7.

Each Picard system is solved with PETSc FGMRES and PCFIELDSPLIT/Schur.
The pressure subsolve defaults to HYPRE BoomerAMG.

## Energy equation

After flow converges,

[
  rac{C_L}{Delta t}
  (T^{n+1}-T^n)
  + A_T(u^{n+1})T^{n+1}
  = q.
]

The current physical-BE thermal path uses the low-order graph-viscosity operator
for all materials and does not clip the resulting temperature.

## Activation

Use the transient case with both

`CBS3D_FULLY_IMPLICIT_FLOW=1`

and

`CBS3D_IMPLICIT_TRANSIENT_ENERGY=1`.

For the blanket test helper:

`bash jobs/sunbird/templates/blanket.sh fullybe 1`

Then increase to 5, 20, etc. only after the one-step invariants pass.

## Verification gates

Do not start a long production run until all of these pass:

1. MPI/PETSc matrix creation and one Picard solve complete on 40 ranks.
2. Coupled linear residual is finite and within tolerance.
3. Picard (RelU), (RelP), and divergence converge.
4. Mass-flow inlet/outlet closure is consistent with the existing CBS path.
5. One-step temperature remains above the 298.15 K physical floor without
   clipping.
6. BE temporal refinement ((Delta t,Delta t/2,Delta t/4)) shows first-order
   convergence for transient observables.
7. Small-(Delta t) BE results agree with the established transient CBS result.
8. Only after correctness: increase physical (Delta t) and optimise persistent
   PETSc matrix/KSP reuse across timesteps.

## Current limitations

- laminar flow only,
- P1/TET4 only,
- distributed MPI + PETSc only,
- no distributed symmetry constraints in the monolithic flow path,
- temperature-to-flow feedback is not implemented,
- the current implementation rebuilds the monolithic PETSc objects per physical
  timestep; persistence is a performance optimisation, not a change to the BE
  formulation.
