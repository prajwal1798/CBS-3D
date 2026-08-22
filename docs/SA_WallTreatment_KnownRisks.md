# SA Wall-Treatment Known Risks

The following points are intentionally unresolved and must not be hidden by the foundation implementation.

1. **Boundary-node semantics:** current BC 530/532 and material-interface logic strongly zero all velocity components. Production wall modelling requires a separately derived normal-only constraint.
2. **Natural-flux replacement:** current Step-1 diffusion reconstructs a natural viscous face flux. The wall-model traction must replace that contribution on modelled faces, not be added to it.
3. **CHT interface nodes:** conformal fluid/solid nodes share storage. Relaxing tangential velocity for a fluid wall model must not accidentally introduce solid convection or alter temperature continuity.
4. **Sampling point:** a TET4 wall face has a natural opposite node, but the appropriate wall-model sampling distance on skewed elements and at corners still requires verification.
5. **SA transported-variable boundary condition:** the current `nu_tilde=0` condition is wall-resolved SA behaviour. A model-consistent high-Re SA scalar treatment has not yet been selected or implemented.
6. **Thermal resistance:** turbulent `k_eff` alone is not a high-Re thermal wall treatment.
7. **QCR:** the current scalar `mu_eff` diffusion assembly cannot represent QCR2000 tensor stresses.
8. **LTS:** wall-model activation changes near-wall residual stiffness. Stable convergence and equality with the global-pseudo-time steady state must be re-established.
9. **MPI:** physical wall faces must be assembled once globally; shared nodal wall residuals must reverse-sum to owners before updates.
10. **Geometry:** a wall function reduces wall-normal resolution requirements but does not excuse poor tetrahedral quality, abrupt growth or under-resolution of separation/secondary-flow structures.

These are development gates, not optional post-processing checks.
