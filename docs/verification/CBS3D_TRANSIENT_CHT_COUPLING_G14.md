# G14 — transient segregated CBS–CHT coupling

**Status (2026-11-02): implemented in the production MPI loop; not yet executed on Sunbird2.**

## Research basis

The implementation follows the standard segregated ordering of the Characteristic-Based Split (CBS) method: momentum predictor, pressure correction, velocity correction, then solution of additional scalar equations such as temperature. This ordering is described by Nithiarasu, Codina & Zienkiewicz, *The Characteristic-Based Split (CBS) scheme—a unified approach to fluid dynamics*, Int. J. Numer. Meth. Engng. 66 (2006) 1514–1546, DOI **10.1002/nme.1698**. Conjugate heat-transfer use of the CBS framework is also reported by Nithiarasu, *A unified fractional step method for compressible and incompressible flows, heat transfer and incompressible solid mechanics*, Int. J. Numer. Meth. Heat & Fluid Flow 18 (2008) 111–130, DOI **10.1108/09615530810846284**.

For transient CHT, the thermal equation is advanced by backward Euler after the corrected velocity is available:

[
\left(\frac{C}{\Delta t}+A_L(\mathbf u^{n+1})\right)T^{n+1}
=
\frac{C}{\Delta t}T^n+b^{n+1}.
]

This is a standard loosely/segregated coupled treatment for the present one-way thermal coupling (constant properties, no buoyancy/temperature feedback). The current fluid and solid regions use one conforming continuous P1 temperature field, so temperature continuity is inherent at the shared interface and the assembled diffusion operator enforces the weak heat-flux balance across it.

## Production algorithm

For every physical CBS timestep:

1. Save (\mathbf u^n,p^n,T^n).
2. Compute the communicator-wide physical timestep (\Delta t).
3. CBS Step 1: momentum predictor.
4. CBS Step 2: pressure Poisson solve.
5. CBS Step 3: velocity correction, yielding (\mathbf u^{n+1}).
6. **Implicit Step 4:** assemble the CHT operator using (\mathbf u^{n+1}) and solve one backward-Euler temperature step for (T^{n+1}).
7. Continue with time averaging, residual logging, output and restart/checkpoint handling.

The standalone frozen-flow G13 path remains available and unchanged in purpose. G14 activates the same physical BE operator from inside the real production time loop.

## Controls

Required for the first coupled implementation:

```bash
export CBS3D_IMPLICIT_TRANSIENT_ENERGY=1
export CBS3D_TEMPERATURE_AFC=1
export CBS3D_STEADY_ENERGY_LOW_ORDER_ONLY=1
```

The case must have `temp_calc = 1` and `transient_on = 1`.

The thermal timestep is **not** taken from `CBS3D_IMPLICIT_ENERGY_BE_DT` in coupled mode. It is the same `global_dt` used by the CBS flow step. This prevents flow and energy from drifting onto different time levels.

The present coupled path uses the validated low-order graph-viscosity operator only. The nonlinear high-order AFC limiter is not claimed as an implicit transient method.

## Important state-history correction

The G13 residual helper temporarily assigns `temperature1 = temperature` because the legacy explicit energy assembly reads `temperature1`. In a coupled transient run this would overwrite the production history (T^n), corrupting the reported temperature increment and potentially restart/convergence diagnostics.

G14 now preserves and restores `temperature1` around residual evaluation. In coupled mode the BE solver leaves `temperature1 = T^n` until the production diagnostics have consumed it.

## PETSc behaviour

The current first implementation rebuilds the thermal matrix/KSP object for each coupled physical timestep because the advection operator changes with (\mathbf u^{n+1}). This is correct but not yet performance-optimal.

The next optimisation, once the coupled run is physically confirmed, is to retain the sparse matrix/KSP objects and update only matrix values. PETSc explicitly supports repeated solves with the same nonzero structure; `MatZeroEntries()` retains the sparse structure, and a changed operator passed to KSP causes the preconditioner to be rebuilt unless reuse is explicitly requested.

## First execution gate

Do **not** create another family of toy verification studies. The first execution should only establish that the integrated production path is sound:

- build succeeds;
- physical transient run reaches several coupled timesteps;
- pressure and thermal KSP solves converge;
- `temperature1` history is preserved (nonzero physical `RelT`);
- no NaN/Inf or temperature-floor violation;
- restart/checkpoint contains the coupled (T^{n+1});
- then move directly to the intended large case and measure memory, wall time and thermal KSP scaling.

After that, the physically necessary validation remains energy conservation, CHT interface flux balance, timestep sensitivity and grid sensitivity. Algebraic residual convergence alone is not a physical validation.
