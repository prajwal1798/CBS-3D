# CBS3D production examples

This directory is the canonical home for reproducible solver cases. Source code
must not be copied into an example. Each case uses one of the production
executables built from the shared numerical core.

## Production hierarchy

```text
examples/
├── Blanket_Laminar/
├── Blanket_Turbulent/
├── Pipe_Flow/
└── Flat_Plate_Turbulent/
```

Each example will use the following internal layout:

```text
<Example>/
├── case/              serial/global CBS input files
├── mesh/              source mesh and generation/export scripts
├── materials/         material-property provenance
├── partitions/        generated MPI partitions or a generation manifest
├── jobs/              portable launch templates
├── postprocessing/    extraction and plotting scripts
├── reference/         accepted numerical values and checksums
├── manifest.json      machine-readable case definition
└── README.md          physics, BCs, execution and validation
```

## Required acceptance information

Every promoted example must document:

- geometry and mesh provenance;
- element/node/boundary counts;
- material identifiers and units;
- boundary-condition identifiers;
- expected serial and MPI launch modes;
- expected mass, momentum and energy checks;
- accepted reference quantities and tolerances;
- whether the case is laminar, thermal or turbulent;
- whether the case is production-validated or development-only.

## Repository-size policy

Small redistributable inputs may be committed directly. Large meshes,
partition folders and result fields should be retained in an external research
data archive or Sunbird project storage and referenced by checksum and stable
location from `manifest.json`.

The current scratch folders will be inventoried before any files are moved or
deleted. Only authoritative inputs, scripts and compact reference outputs will
be promoted into this hierarchy.
