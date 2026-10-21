# Example cases

This directory contains compact, redistributable CBS3D++ examples suitable for
continuous integration and new-user onboarding, together with the production
case hierarchy used by the research workflow.

## Production hierarchy

```text
examples/
  Blanket_Laminar/
  Blanket_Turbulent/
  Pipe_Flow/
  Flat_Plate_Turbulent/
```

Each promoted example uses the following internal layout:

```text
<Example>/
  case/              serial/global CBS input files
  mesh/              source mesh and generation/export scripts
  materials/         material-property provenance
  partitions/        generated MPI partitions or a generation manifest
  jobs/              portable launch templates
  postprocessing/    extraction and plotting scripts
  reference/         accepted numerical values and checksums
  manifest.json      machine-readable case definition
  README.md          physics, BCs, execution and validation
```

## Required contents

At minimum, a runnable example contains:

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
- document geometry, mesh provenance and boundary conditions;
- identify the intended executable and build mode;
- provide expected element, node and boundary counts;
- define a numerical acceptance check;
- avoid machine-specific paths;
- avoid protected legacy source or case archives.

Every promoted example must also document material identifiers and units,
expected serial and MPI launch modes, expected mass/momentum/energy checks,
accepted reference quantities and tolerances, and whether it is laminar,
thermal or turbulent.

## Repository-size policy

Small redistributable inputs may be committed directly. Large meshes, partition
folders and result fields should be retained in an external research-data
archive or Sunbird project storage and referenced by checksum and stable
location from `manifest.json`.

The current scratch folders will be inventoried before any files are moved or
deleted. Only authoritative inputs, scripts and compact reference outputs will
be promoted into this hierarchy.
