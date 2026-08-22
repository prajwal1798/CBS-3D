# SA Wall-Treatment Mathematical Foundation

## Status

This document defines the first mathematically isolated stage of a wall-modelled Spalart-Allmaras development for CBS3D++_SI.

**This stage does not change production BC 530/532/901 behaviour.** Existing wall-resolved calculations therefore remain unchanged. The wall-law kernel is introduced and verified independently before any CBS Step-1 or boundary-state coupling is enabled.

The immediate objective is to remove the requirement that a production RANS mesh resolve the viscous sublayer at `y+ ~= 1`, while retaining the present P1/TET4 CBS architecture.

---

## 1. Current CBS3D momentum boundary state

The live `integration/sa-mpi` implementation currently has two relevant properties:

1. `MomentumAssembly::step1_diffusion()` uses a scalar kinematic-viscosity Laplacian,

   ```text
   R^nu_(a,i)
     = - integral_Omega nu grad(N_a).grad(u_i) dOmega
       + integral_Gamma N_a nu grad(u_i).n dGamma.
   ```

2. BC 530, BC 532 and the conformal CHT interface are presently classified as strong no-slip velocity boundaries. The persistent boundary state therefore sets all three velocity components to zero at those nodes.

A wall model must **not** simply add an empirical wall shear to this existing state. Doing so would simultaneously impose strong zero tangential velocity and a modeled tangential traction, and the present natural diffusion flux could also be counted a second time.

---

## 2. Velocity decomposition at a wall

For a wall unit normal `n` and an off-wall sampled velocity `u_s`, define

```text
u_n = u_s . n
u_t = u_s - (u_s . n) n
U_t = |u_t|
t_hat = u_t / U_t.
```

The wall model must be invariant to the orientation of the geometric normal. Replacing `n` by `-n` therefore leaves `u_t`, `U_t` and the predicted wall shear unchanged.

For the future wall-modelled CBS boundary state, impermeability and tangential momentum treatment are separate conditions:

```text
normal velocity       : u . n = 0
tangential momentum   : supplied by a wall traction law
```

The present wall-resolved mode remains

```text
u_n = 0
u_t = 0
```

strongly.

---

## 3. Continuous Spalding law used for the algebraic foundation

The first kernel uses the smooth-wall Spalding relation because it is continuous through the viscous, buffer and logarithmic regions and therefore avoids a piecewise switch at a prescribed `y+`.

Define

```text
u+ = U_t / u_tau
y+ = u_tau y / nu
u_tau = sqrt(tau_w / rho).
```

The implemented relation is

```text
y+
  = u+
  + exp(-kappa B)
    [ exp(kappa u+)
      - sum_(m=0)^4 (kappa u+)^m / m! ].
```

Default smooth-wall constants are

```text
kappa = 0.41
B     = 5.2
```

and remain explicit options so later model-consistent/automatic-wall-treatment studies can change them without altering the numerical kernel.

At small `u+`, the exponential remainder begins at fifth order, so

```text
y+ -> u+.
```

Consequently the modeled shear has the correct viscous asymptote:

```text
tau_w -> mu U_t / y.
```

At larger wall units the exponential term produces the logarithmic-wall behaviour.

### 3.1 Root variable

Rather than solving directly for `u_tau`, use the exact identity

```text
Re_y = U_t y / nu = u+ y+.
```

The scalar residual is therefore

```text
F(u+) = u+ y+(u+) - Re_y.
```

For non-negative `u+`, Spalding gives

```text
y+(u+) >= u+,
```

so

```text
F(0) = -Re_y < 0,
F(sqrt(Re_y)) >= 0.
```

Thus the root is rigorously bracketed by

```text
0 <= u+ <= sqrt(Re_y).
```

CBS3D uses a safeguarded Newton iteration inside this bracket. If a Newton proposal leaves the bracket, the method falls back to bisection. This retains deterministic convergence while avoiding the cost of pure bisection in the logarithmic regime.

The derivative used by Newton is analytic:

```text
F'(u+) = y+ + u+ dy+/du+,
```

with

```text
dy+/du+
  = 1
    + exp(-kappa B) kappa
      [ exp(kappa u+)
        - sum_(m=0)^3 (kappa u+)^m / m! ].
```

The exponential remainders are evaluated with a series near zero to avoid catastrophic cancellation.

---

## 4. Wall shear returned by the kernel

Once `u+` is known,

```text
u_tau = U_t / u+,
tau_w = rho u_tau^2.
```

The traction exerted by a stationary wall on the fluid is defined to oppose the tangential sampled velocity:

```text
t_w = -tau_w t_hat.
```

Therefore

```text
t_w . u_t <= 0,
```

which is an explicit dissipation invariant checked by the unit test.

The present CBS Step-1 diffusion equation is written in kinematic form. The corresponding traction quantity required by that residual is therefore

```text
g_w = t_w / rho = -u_tau^2 t_hat.
```

Its units are `m^2/s^2`, identical to `nu du/dn` in the existing natural diffusion term.

---

## 5. Required future CBS weak-form coupling

For a wall-modelled TRI3 boundary face `Gamma_w`, the existing reconstructed natural diffusion flux on that face must be **replaced**, not augmented, by the wall-model flux:

```text
integral_Gamma_w N_a nu grad(u_i).n dGamma

          becomes

integral_Gamma_w N_a g_w,i dGamma.
```

If the face traction is treated as constant over one linear triangular face of area `A_f`,

```text
R_wall_(a,i) = (A_f / 3) g_w,i,
```

for each of its three face nodes.

The sign must be checked against the final Step-1 residual convention using a synthetic single-face test before production activation.

### Non-negotiable rule

Do not simultaneously:

```text
1. strongly set u_t = 0,
2. retain the reconstructed natural viscous face flux, and
3. add the wall-model traction.
```

That would over-constrain and/or double-count the wall momentum transfer.

---

## 6. Future nodal boundary semantics

The current persistent velocity states distinguish free, prescribed, moving-wall and no-slip nodes. Wall-model mode will require an additional concept representing **impermeability without strong tangential no-slip**.

The existing BC506 symmetry projector is useful precedent because it removes velocity components in the span of incident face normals. A separate wall-model projector/state should be implemented rather than silently reusing BC506 semantics, because wall-model faces also carry non-zero modeled tangential stress.

At wall intersections, all incident wall normals must be respected. A single averaged normal is insufficient.

For conformal CHT interface nodes, the velocity field is inactive in solid elements, so the wall-modelled fluid momentum semantics must be derived explicitly before relaxing the present material-interface zero-velocity priority. Temperature continuity remains a separate shared-node condition.

---

## 7. Future wall-face inventory

Before Step-1 coupling, preprocessing should create one owned record per physical fluid wall face containing at least

```text
boundary face id
BC id
face nodes
adjacent owned fluid tetrahedron
opposite off-wall node
face area
unit normal
face centroid
sample distance
thermal wall type
```

For a TET4 wall face the opposite tetrahedral node provides a natural first off-wall velocity sample. The production implementation must nevertheless verify its face-normal distance and corner behaviour rather than assuming tetrahedral altitude is always the appropriate wall-model height.

Partition interfaces are never physical turbulence walls.

---

## 8. MPI assembly rule

A physical wall face contribution is assembled exactly once by the rank that owns the adjacent fluid tetrahedron / physical boundary face.

The wall residual is additive, so shared-node contributions follow the existing CBS distributed rule:

```text
local owned-face assembly
    -> rank-local nodal residual
    -> reverse SUM to node owners
    -> owner velocity update
    -> wall normal constraint
    -> forward COPY owner state to ghosts.
```

No wall shear may be duplicated merely because a shared wall node has ghost copies.

---

## 9. Thermal wall treatment is deliberately not included in this stage

The final first-wall CHT application also requires a model for unresolved thermal resistance. The existing

```text
k_eff = k + rho cp nu_t / Pr_t
```

is a bulk RANS closure and is not, by itself, a thermal wall model when the first off-wall point lies in the logarithmic region.

Thermal wall treatment will therefore be developed only after the momentum wall treatment passes flat-plate/channel verification. Its face heat flux must replace the appropriate unresolved near-wall thermal transport rather than being added on top of an already counted conductive flux.

---

## 10. QCR is also a later stage

QCR2000 changes the turbulent stress tensor. It must not be represented by another scalar `mu_eff` value. Step 1 will need a full turbulent-stress divergence/weak-form refactor after the wall-model baseline is validated.

---

## 11. Verification gates before production activation

The wall-model mode is not acceptable for the breeder-blanket first wall until all of the following have passed:

1. Algebraic wall-law unit tests.
2. Single TRI3 face traction sign/conservation test.
3. Laminar regression with wall treatment disabled.
4. Existing wall-resolved SA flat-plate regression.
5. Wall-modelled flat plate over multiple first-point `y+` values.
6. Pressure-drop and skin-friction comparison against reference data.
7. MPI rank-independence of wall traction and velocity fields.
8. LTS versus global-pseudo-time equality of the converged steady solution.
9. Heated-wall thermal wall-treatment verification.
10. CHT energy balance.
11. Square-duct QCR verification.
12. Helium/Eurofer first-wall sensitivity study.

---

## 12. References used for this foundation

- D. B. Spalding, “A Single Formula for the Law of the Wall,” *Journal of Applied Mechanics*, 28(3), 455-458, 1961. DOI: `10.1115/1.3641728`.
- T. Knopp, T. Alrutz and D. Schwamborn, “A grid and flow adaptive wall-function method for RANS turbulence modelling,” *Journal of Computational Physics*, 220(1), 19-40, 2006.
- J.-R. Carlson, V. N. Vatsa and J. White, “Validation of a Node-Centered Wall Function Model for the Unstructured Flow Code FUN3D,” AIAA 2015-2758 / NASA NTRS `20160005970`, 2015. This work is particularly relevant because it couples an SA wall-function treatment to a three-dimensional unstructured solver and weakly applies wall momentum/energy fluxes.

The present CBS3D kernel is intentionally only the algebraic foundation. It does not claim to reproduce the complete Knopp/FUN3D model-consistent SA treatment yet.
