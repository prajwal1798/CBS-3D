# Physical-time implicit thermal CHT: G13a development note

**Status (2026-11-02): IMPLEMENTED ON AUDIT BRANCH, NOT BUILT OR RUN ON SUNBIRD2.** Do not use for production data until measured convergence, energy budgets, time-step and grid sensitivity are accepted. The original executable and U-channel inputs are unchanged.

## Scope and equations

For an imported **frozen** 3D velocity field, fixed thermal material properties and continuous P1 tetrahedra, the new *opt-in* thermal backward-Euler path solves

\[
  \left(\frac{C}{\Delta t}+A_L\right)T^{n+1}
  = b+\frac{C}{\Delta t}T^n,\qquad
  T_i^{n+1}=T_{D,i}\quad(i\in\Gamma_D).
\]

- \(C\): diagonal **physical lumped heat capacity**, \(C_i=\sum_{K\ni i} (\rho c_p)_K |K|/4\), assembled in MPI from contributions on every local element and summed to owners.
- \(A_L\): *steady* graph-viscosity-stabilized energy matrix, assembled using the **same** \`EnergyAssembly::elementOperator\` as CBS3D. The fluid elements receive conservative symmetric graph viscosity \(d_{ij}=\max(0,-M_{ij},-M_{ji})\). Solid diffusion remains standard FEM, i.e. not generally M-matrix on arbitrary tetrahedra.
- \(b=c+q\): constant source/reference-temperature contribution obtained from the exact existing element thermal RHS evaluated at \(T=0\), plus the *same* prescribed face-flux load (\`BC 532\`) used by production CBS.
- In the opt-in physical-time mode, `EnergyAssembly::add_fluid_convection_stabilisation` explicitly accepts the scoped zero pseudo-time and returns without adding characteristic diffusion; outside this mode, non-positive element timestep remains an error.
- The CBS characteristic thermal stabilization term \(\frac12\Delta t_\text{pseudo}\rho c_p V(\mathbf u\cdot\nabla N_i)(\mathbf u\cdot\nabla N_j)\) is **disabled only for physical-time BE assembly**. Physical \(\Delta t\) enters solely through \(C/\Delta t\). This removes the dependency of the physical transient spatial operator on the solver's unrelated local pseudo-time controls.
- Fixed thermal boundary values occupy identity PETSc rows, and the mass term is added only for free nodes.
- PETSc GMRES + ASM/ILU(1) is initialized for the single fixed matrix; the RHS is updated and KSP reused at each physical-time step.
- No temperature clipping, no arbitrary field rescaling. Any optional physical floor violation raises an error.

This is **not** a fully coupled implicit Navier–Stokes + thermal CHT formulation; velocity and coefficients do not evolve with temperature. It is a consistent implicit discretization of the *frozen-velocity, graph-viscosity-stabilized* semi-discrete heat equation. For a digital twin, evolving flow, material nonlinearities and the energy advection operator will need further work.

### Relevant controls

- \`CBS3D_STEADY_ENERGY=1\` (production entry point; requires PETSc/MPI)
- \`CBS3D_TEMPERATURE_AFC=1\` + \`CBS3D_STEADY_ENERGY_LOW_ORDER_ONLY=1\` (enables low-order fluid graph viscosity; **does not run the nonlinear limiter**)
- \`CBS3D_STEADY_ENERGY_PSEUDOTIME_SHIFT=0\` (mandatory for physical BE, avoids double time mass)
- \`CBS3D_IMPLICIT_ENERGY_BE_DT=0.01\` (positive physical timestep in seconds)
- \`CBS3D_IMPLICIT_ENERGY_BE_STEPS=10\` (number of frozen-flow steps)
- \`CBS3D_IMPLICIT_ENERGY_BE_MIN_TEMPERATURE=298.15\` (optional *audit floor*; zero clipping)

The PVD/VTU final physical time is **relative to the imported restart temperature**, \(t_\text{final}=N\Delta t\); it must not inherit the CBS pseudo-time or be interpreted as the elapsed time of the original steady solver restart. Only the **last** temperature snapshot is exported by this first experimental path. Full time-series output and changing fluid states are future requirements for surrogate training.

## G13a first experiment

\`jobs/sunbird/templates/benchmark_physical_be_cht_g13a.slurm\` compiles a separate candidate binary in scratch (no overwrite of original). It uses the verified 2-rank original U-channel restart at global step 50000 and SHA256-checked input \`.par\`/\`.matprop\`. Both test arms evolve to exactly **0.1 physical seconds**:

| Arm | Physical dt | BE steps | Final physical time |
|---|---:|---:|---:|
| Coarser | 0.01 s | 10 | 0.1 s |
| Finer | 0.001 s | 100 | 0.1 s |

Acceptance is limited to: successful clean MPI/PETSc solve for **both** arms, independent matrix-residual tolerance on each BE step (\(<10^{-6}\)), independent owned-node final VTU temperature check (9613 nodes), **every intermediate step** satisfying \(T_\min\ge298.15\) K, and measured maximum temperature difference \(\|T_{\Delta t=0.01}-T_{\Delta t=0.001}\|_\infty\) at equal final physical time. These are *necessary, not sufficient*, for transient accuracy.

## Subsequent physical/accuracy gates (NOT completed)

1. **Thermal energy budget:** measure \(\frac{U^{n+1}-U^n}{\Delta t}\), integrated \(Q_\mathrm{nuclear}+Q_{532}\), inlet/outlet advective enthalpy, wall/interface conductive flux and boundary reaction flux. Report signed residual in W and its relative tolerance; do not confuse linear \(\|A T-b\|\) with the full physical budget.
2. **Conjugate-interface continuity:** continuity of FEM temperature and signed balance of heat flux across solid/fluid internal faces, with mesh refinement. Source power in original case: \(\sim16.246\) W volumetric + 60 W BC 532 = \(\sim76.246\) W whole-domain input (flow exits carry enthalpy).
3. **Spatial/time accuracy:** compare against verified analytical conduction/advective benchmarks and \(\Delta t\) refinements (two time levels alone do not establish temporal order). Check hot-spot \(T_\max\), interface heat flux, outlet enthalpy and total energy, not only boundedness.
4. **High-order bounded thermal transport:** implement well-balanced *steady-source-aware* monolithic convex limiting / conservative antidiffusive edge corrections with a demonstrably convergent nonlinear solver. The **transient Zalesak** AFC reused at steady state failed G12b/G12c. Avoid presenting G13a's low-order solution as high-fidelity or as a successful high-order implicit AFC scheme.
5. **Coupled CHT:** integrate thermal update with CBS pressure/velocity and temperature-dependent \(\rho,\mu,k,c_p\) and possibly buoyancy; validate nonlinear coupling and conservation. This is distinct from frozen-flow energy BE.

## References for high-order follow-on

- Knobloch, Kuzmin & Jha, *Well-balanced convex limiting for finite element discretizations of steady convection-diffusion-reaction equations*, JCP 518, 113305 (2024), DOI 10.1016/j.jcp.2024.113305. Treat steady source terms as part of the flux-limiting/equilibrium design, rather than clipping source-driven extrema.
- Jha et al., *An assessment of solvers for algebraically stabilized discretizations of convection-diffusion-reaction equations*, J. Numer. Math. 31 (2023) 79–103, DOI 10.1515/jnma-2021-0123.

Neither paper is implemented wholesale here; this document records a forward design direction, not a claim of reproducing the published algorithm.
