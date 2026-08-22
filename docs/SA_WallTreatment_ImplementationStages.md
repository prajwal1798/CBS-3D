# SA Wall-Treatment Implementation Stages

The production integration is intentionally staged so a defect can be localized to one numerical layer.

## Stage 0: algebraic foundation

Current branch scope.

- continuous smooth-wall Spalding kernel;
- robust wall-unit root solve;
- tangential/dissipative traction construction;
- algebraic verification and wide-`Re_y` sweep;
- no production boundary semantics changed.

## Stage 1: wall-face inventory only

Add preprocessing records for physical fluid wall faces. Do not alter Step 1 yet.

Required audits:

- each physical wall face is represented exactly once globally;
- adjacent tetrahedron is fluid and locally owned;
- TRI3 area and normal are finite;
- opposite/sample node is valid;
- partition interfaces never enter the inventory;
- P1 versus P40 inventory counts are identical globally.

## Stage 2: impermeability state only

Introduce an explicit wall-model velocity state that constrains the span of incident wall normals while leaving tangential velocity unconstrained.

Required tests:

- one planar wall removes one normal component;
- two intersecting walls remove two independent normal components;
- corner rank-three case is deterministic;
- normal orientation does not matter;
- resolved-wall mode remains bitwise/round-off equivalent to the current implementation.

## Stage 3: synthetic TRI3 traction

Before calling the wall law from Step 1, apply a prescribed constant tangential traction on a single boundary face.

Required checks:

```text
sum_a R_face,a = A_f g_w
R_face,a = A_f g_w / 3  for a linear TRI3 and constant traction
```

The resulting sign must oppose a prescribed tangential velocity in a one-element/patch test.

## Stage 4: wall-law Step-1 coupling

Replace the current reconstructed natural viscous boundary flux on wall-model faces by the wall-model traction. Do not add both.

Required checks:

- no strong tangential zero remains on modelled walls;
- normal impermeability remains active;
- modeled wall work is non-positive;
- total wall force equals the assembled face-force sum;
- serial and distributed residuals agree.

## Stage 5: SA-consistent near-wall scalar treatment

The present `nu_tilde=0` wall condition belongs to wall-resolved SA. The wall-function formulation for the SA transported variable must be selected from a verified model-consistent treatment and implemented separately rather than guessed from the momentum law.

This stage requires comparison against published SA wall-function reference results.

## Stage 6: flat-plate verification

Use the existing wall-resolved SA reference and deliberately vary first-point location.

Minimum comparison set:

```text
y+ target: 1, 5, 10, 30, 50, 100
Cf(x)
u+(y+)
pressure distribution
nu_tilde
mu_t/mu
```

A wall model is acceptable only if outer-flow and skin-friction predictions become substantially less sensitive to first-point placement than wall-resolved SA on under-resolved meshes.

## Stage 7: thermal wall treatment

Derive and verify the thermal wall flux in the energy weak form. The wall-law heat flux must replace the unresolved near-wall flux contribution consistently and must not double-count the existing volume conductivity.

Verification precedes CHT:

```text
heated channel / flat plate
q_w
T_w
T+
Nu
bulk-temperature rise
global energy balance
```

## Stage 8: QCR2000

Refactor Step 1 to accept a full turbulent stress tensor and then implement QCR2000. Do not encode QCR as a scalar viscosity multiplier.

Validate on a square-duct secondary-flow benchmark.

## Stage 9: breeder-blanket first wall

Only after Stages 0-8 pass:

- helium/Eurofer CHT;
- realistic manifold/channel flow distribution;
- target `y+` sensitivity;
- mesh-quality sensitivity;
- LTS/global-pseudo-time steady equivalence;
- P2/P40 or broader MPI rank independence;
- global mass and energy conservation.
