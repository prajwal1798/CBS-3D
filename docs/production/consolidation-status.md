# Production consolidation status

Updated: 2026-08-06

Branch: `integration/production-consolidation`

## Completed in Phase 1

- created dedicated serial application entry point: `apps/serial/main.cpp`;
- created dedicated parallel application entry point: `apps/parallel/main.cpp`;
- retained one shared numerical source list for both applications;
- added an explicit runtime rejection for distributed Spalart-Allmaras cases;
- introduced CMake target names `cbs3d_serial` and `cbs3d_parallel`;
- retained optional legacy executable output name `cbs3dpp_si` for existing jobs;
- added serial, OpenMP, serial-PETSc, MPI, hybrid and SA-development presets;
- added Linux serial/OpenMP and Windows MSVC CI builds;
- added a one-job Sunbird build matrix for all production configurations;
- added a non-destructive Sunbird workspace inventory tool;
- established the requested examples hierarchy;
- documented the validated straight-pipe mass-flow benchmark.

## Continuous-integration result

GitHub Actions run `31074698723` completed successfully for:

- Linux deterministic serial build;
- Linux OpenMP build;
- Windows MSVC Release serial build.

This verifies the new application split and CMake structure for the non-PETSc
build modes.

## Validation still required

The following must be validated on Sunbird before this phase is promoted:

1. serial PETSc build and one-case regression;
2. MPI/PETSc build and one-iteration 40-rank regression;
3. corrected BC511 mass-flow verification;
4. serial versus MPI field comparison;
5. hybrid MPI/OpenMP build and rank/thread scaling smoke test;
6. confirmation that a turbulent MPI parameter file terminates with the new
   unsupported-feature diagnostic.

## Next integration phase

After the Sunbird build matrix passes:

1. import the Step-4 validation tests selectively;
2. reconcile `TemperatureAFC` and `ThermalAfc` into one implementation;
3. remove `.cpp`-inclusion production wrappers;
4. import public documentation and repository governance files from `main`;
5. inventory scratch directories and identify authoritative example inputs;
6. move reproducible builds and retained results into the controlled Sunbird
   workspace hierarchy;
7. archive or remove redundant source clones only after commit and content
   verification.
