# Public example cases

This directory is reserved for compact, redistributable CBS3D++ examples suitable for continuous integration and new-user onboarding.

## Required contents of each example

```text
<case>.plt
<case>.bco
<case>.par
<case>.material      when required
<case>.matprop       when required
README.md
reference/           compact accepted values or checksums
scripts/             extraction/comparison scripts when needed
```

## Acceptance requirements

An example must:

- be legally redistributable;
- run in practical CI time;
- document physics and boundary conditions;
- identify the intended executable/build mode;
- provide expected startup totals;
- define a numerical acceptance check;
- avoid machine-specific paths;
- avoid protected legacy source or case archives.

## Planned initial examples

1. Small tetrahedral pressure/operator smoke test.
2. Laminar lid-driven cavity.
3. Analytical channel-flow verification.
4. Thermal diffusion or advection-diffusion verification.
5. Small conformal fluid-solid heat-transfer case.

Large research meshes should be stored in a research-data archive and referenced by stable identifier rather than committed directly to Git.
