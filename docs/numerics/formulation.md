# Mathematical formulation and code architecture

This page maps the CBS3D++ solver through four levels:

$$\boxed{\text{strong PDE}\;\longrightarrow\;\text{weak finite-element form}\;\longrightarrow\;\text{element algebra}\;\longrightarrow\;\text{C++ arrays and updates}}$$

The objective is not only to state the CBS equations. Every major term is connected to its test function, tetrahedral coefficient, residual array and update routine.

> **Implementation boundary.** CBS3D++ currently uses a semi-implicit incompressible CBS pressure-projection formulation on four-node linear tetrahedra. General CBS notation is shown for context, but the boxed **implemented equations** correspond to the active C++ source.

## 1. Scope, indices and stored variables

The active formulation contains:

- three-dimensional incompressible flow;
- P1 tetrahedral finite elements;
- element-wise material properties;
- lumped momentum and thermal mass/time diagonals;
- CBS momentum prediction, pressure solution and velocity correction;
- optional fluid-solid temperature transport;
- optional Spalart-Allmaras transport;
- OpenMP-coloured element assembly;
- native matrix-free PCG or PETSc pressure solution.

The present iteration is steady or pseudo-transient. A validated physical-time BDF or dual-time path is not yet released.

### 1.1 Indices

$$i,j,k\in\{1,2,3\},\qquad x_1=x,\quad x_2=y,\quad x_3=z.$$

Local tetrahedral nodes are $a,b\in\{1,2,3,4\}$, elements are indexed by $e$, and global mesh nodes by $I,J$.

### 1.2 Conservative theory variable and code variable

CBS derivations commonly introduce momentum per unit volume,

$$U_i=\rho u_i.$$

The active arrays `unkn1` and `unkno` store velocity components $u_i$, not a separate conservative momentum field. Code-aligned equations below therefore use velocity form after division by density where required.

| Mathematical quantity | C++ representation |
|---|---|
| coordinates $x_i$ | `coord(i,I)` |
| local-to-global map $I=I_e(a)$ | `intma(a,e)` |
| $N_{a,i}^{(e)}=\partial N_a/\partial x_i$ | flattened `dNkdx` |
| $\det J_e=6V_e$ | `detJ(e)` |
| previous velocity $u_i^n$ | `unkn1(i,I)` |
| predicted/corrected velocity | `unkno(i,I)` |
| pressure $p$ | `pres(I)` |
| previous temperature $T^n$ | `temperature1(I)` |
| updated temperature | `temperature(I)` |
| Step-1 or Step-3 residual | `rhs(i,I)` |
| Step-2 or Step-4 residual | `rhs1(I)` |
| local element time scale $\Delta t_e$ | `delte(e)` |
| $D_u^{-1}$ | `elcoe2(I)` |
| $D_T^{-1}$ | `elcoe2p(I)` |

`rhs` and `rhs1` are work arrays. Their mathematical meaning changes with the CBS stage, so documentation must identify the active step whenever either name is used.

## 2. Governing equations

### 2.1 Continuity

For incompressible flow,

$$\boxed{\frac{\partial u_i}{\partial x_i}=0}\qquad\Longleftrightarrow\qquad\boxed{\nabla\!\cdot\mathbf{u}=0}.$$

### 2.2 Momentum

The general incompressible Newtonian momentum equation is

$$\rho\left(\frac{\partial u_i}{\partial t}+u_j\frac{\partial u_i}{\partial x_j}\right)=-\frac{\partial p}{\partial x_i}+\frac{\partial\tau_{ij}}{\partial x_j}+\rho b_i.$$

The complete Newtonian deviatoric stress is

$$\tau_{ij}=\mu\left(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}-\frac{2}{3}\delta_{ij}\frac{\partial u_k}{\partial x_k}\right).$$

The active `MomentumAssembly.cpp` diffusion kernel does **not** assemble this complete tensor. It currently contracts $\nabla N_a$ with $\nabla u_i$, which corresponds to a component-wise Laplacian:

$$\boxed{\left.\frac{\partial u_i}{\partial t}\right|_{\mathrm{diff}}=\nu\nabla^2u_i}.$$

Body force belongs to the governing equation but is not assembled in the current Step-1 core.

### 2.3 Temperature

For a fluid element,

$$\boxed{\rho c_p\left(\frac{\partial T}{\partial t}+u_i\frac{\partial T}{\partial x_i}\right)=\frac{\partial}{\partial x_i}\left(k\frac{\partial T}{\partial x_i}\right)+Q}.$$

For a solid element, $\mathbf{u}=\mathbf{0}$ and only transient/pseudo-transient conduction remains:

$$\rho c_p\frac{\partial T}{\partial t}=\nabla\!\cdot(k\nabla T)+Q.$$

At a conformal fluid-solid interface, both materials share one nodal temperature degree of freedom. Diffusion is assembled from the adjacent fluid and solid elements; no duplicate interface temperature is introduced.

## 3. P1 tetrahedral finite elements

### 3.1 Interpolation

For a four-node tetrahedron $\Omega_e$,

$$u_i^h(\mathbf{x})=\sum_{a=1}^{4}N_a(\mathbf{x})u_{i,a},\qquad p^h(\mathbf{x})=\sum_{a=1}^{4}N_a(\mathbf{x})p_a,\qquad T^h(\mathbf{x})=\sum_{a=1}^{4}N_a(\mathbf{x})T_a.$$

The affine physical map is

$$\mathbf{x}=\mathbf{x}_1+J_e\begin{bmatrix}\xi&\eta&\zeta\end{bmatrix}^{T},\qquad J_e=\begin{bmatrix}x_2-x_1&x_3-x_1&x_4-x_1\\y_2-y_1&y_3-y_1&y_4-y_1\\z_2-z_1&z_3-z_1&z_4-z_1\end{bmatrix}.$$

The positive-orientation convention gives

$$\boxed{V_e=\frac{\det J_e}{6}}.$$

### 3.2 Constant gradients and flattened storage

For P1 tetrahedra,

$$N_{a,i}^{(e)}\equiv\frac{\partial N_a}{\partial x_i}=\text{constant in }\Omega_e.$$

Hence

$$\left.\frac{\partial u_i^h}{\partial x_j}\right|_e=\sum_{a=1}^{4}u_{i,a}N_{a,j}^{(e)},\qquad\left.\frac{\partial p^h}{\partial x_i}\right|_e=\sum_{a=1}^{4}p_aN_{a,i}^{(e)}.$$

The code flattens the logical object $N_{a,i}^{(e)}$ using

$$\mathcal{I}(e,i,a)=(e-1)(3)(4)+(i-1)(4)+a.$$

`dNkdx_index(...)` computes this address; the derivative values are stored in `s.dNkdx`.

### 3.3 Exact integration identities

The principal P1 tetrahedral identities are

$$\int_{\Omega_e}N_a\,d\Omega=\frac{V_e}{4}=\frac{\det J_e}{24}.$$

$$\int_{\Omega_e}N_a^2\,d\Omega=\frac{V_e}{10},\qquad\int_{\Omega_e}N_aN_b\,d\Omega=\frac{V_e}{20}\quad(a\ne b).$$

For a triangular face $\Gamma_f$,

$$\int_{\Gamma_f}N_a\,d\Gamma=\frac{A_f}{3},\qquad\int_{\Gamma_f}N_a^2\,d\Gamma=\frac{A_f}{6},\qquad\int_{\Gamma_f}N_aN_b\,d\Gamma=\frac{A_f}{12}\quad(a\ne b).$$

These explain the code constants:

| Code expression | FEM meaning |
|---|---|
| `detJ * fcon[1]` | $\det J_e/24=V_e/4$ |
| `fcon[2]` | $1/12$ for consistent triangular-face interpolation |
| `detJ * fdif[1]` | $\det J_e/6=V_e$ |
| `fdif[2]` | $1/3$ for $\int_{\Gamma_f}N_a\,d\Gamma$ |

### 3.4 Consistent mass, lumping and time diagonals

The consistent scalar element mass matrix is

$$M^{(e)}=\frac{V_e}{20}\begin{bmatrix}2&1&1&1\\1&2&1&1\\1&1&2&1\\1&1&1&2\end{bmatrix}.$$

Row-sum lumping gives

$$\boxed{M_L^{(e)}=\frac{V_e}{4}I_4}.$$

Preprocessing stores $V_e/4$ in `elcoe_e(e,a)`. The field updates do not use a general $M^{-1}$; they use time-dependent diagonal inverses:

$$[D_u]_{II}=\sum_{e\ni I,\;e\in\Omega_f}\frac{V_e/4}{\Delta t_e},\qquad\boxed{\texttt{elcoe2}(I)=[D_u]_{II}^{-1}}.$$

$$[D_T]_{II}=\sum_{e\ni I}\frac{(\rho c_p)_eV_e/4}{\Delta t_e},\qquad\boxed{\texttt{elcoe2p}(I)=[D_T]_{II}^{-1}}.$$

Only fluid elements contribute to $D_u$; fluid and solid elements contribute to $D_T$.

## 4. Solver information flow

```text
mesh: coord, intma, iside
        │
        ▼
geometry preprocessing
(detJ, dNkdx, face normals, face classes)
        │
        ▼
mass and pressure preprocessing
(elcoe_e, pdiagE, gstifE)
        │
        ▼
TimeStep::updateLhsDiagonal
(delte, elcoe2, elcoe2p, pdiag, gstif)
        │
        ├── Step 1: MomentumAssembly       → u*
        ├── Step 2: PressureAssembly + KSP → p
        ├── Step 3: Steps.cpp correction   → u^(n+1)
        ├── optional SA transport           → turbulent properties
        └── Step 4: EnergyAssembly          → T^(n+1)
```

At the start of an iteration, `unkn1` represents $\mathbf{u}^n$ and `temperature1` represents $T^n$. After Step 1, `unkno` represents $\mathbf{u}^*$; after Step 3, the same array represents $\mathbf{u}^{n+1}$.

## 5. CBS Step 1: pressure-free momentum predictor

### 5.1 General CBS interpretation

After division by density, define the non-pressure operator

$$R_i(\mathbf{u})=-\frac{\partial(u_ju_i)}{\partial x_j}+\nu\nabla^2u_i+b_i.$$

A characteristic expansion has the structure

$$\Delta u_i^*=\Delta t\,R_i^n-\frac{\Delta t^2}{2}\frac{\partial}{\partial x_k}\left(u_kR_i^n\right)+\cdots.$$

The active implementation retains Galerkin convection, the convective characteristic correction and component-wise viscous diffusion. It omits body force, artificial diffusion and characteristic derivatives of viscous/body-force terms.

### 5.2 Implemented nodal equation

The Step-1 element residuals assemble into

$$\boxed{D_u\left(\mathbf{u}^*-\mathbf{u}^n\right)=\mathbf{r}_m},\qquad\mathbf{r}_m=\mathbf{r}_{\mathrm{conv}}+\mathbf{r}_{\mathrm{char}}+\mathbf{r}_{\mathrm{diff}}.$$

The nodal update is

$$\boxed{\mathbf{u}^*=\mathbf{u}^n+D_u^{-1}\mathbf{r}_m}.$$

```cpp
s.unkno(idim, ip) =
    s.unkn1(idim, ip)
    + s.rhs(idim, ip) * s.elcoe2(ip);
```

### 5.3 Galerkin convection

At local node $b$, define the momentum-flux product

$$F_{ij,b}=u_{j,b}u_{i,b}.$$

The strong Step-1 convective contribution is $-\partial F_{ij}/\partial x_j$. Multiplication by $N_a$ and integration by parts gives

$$r_{\mathrm{conv},i,a}^{(e)}=\int_{\Omega_e}N_{a,j}F_{ij}\,d\Omega-\int_{\Gamma_e}N_aF_{ij}n_j\,d\Gamma.$$

With $F_{ij}^h=N_bF_{ij,b}$ and constant $N_{a,j}$,

$$\boxed{r_{\mathrm{conv},i,a}^{(e),\Omega}=N_{a,j}^{(e)}\frac{V_e}{4}\sum_{b=1}^{4}F_{ij,b}}.$$

Theory-to-code correspondence:

| Mathematical operation | C++ object |
|---|---|
| gather $u_{i,b}^n$ | `lunkno[i][b]` from `unkn1` |
| form $F_{ij,b}$ | `lunk[i][j][b]` |
| compute $\sum_bF_{ij,b}$ | `lunksum[i][j]` |
| $N_{a,j}^{(e)}$ | `grad(s,ie,j,a)` |
| $V_e/4$ | `volume_factor = detJ * fcon[1]` |
| local residual | `lrhs[i][a]` |
| global assembly | `rhs(i,I) += lrhs[i][a]` |

The face term uses the consistent $A_f/12$ coefficient through `annxf` and `fcon[2]`.

### 5.4 Characteristic correction

The implemented volume contribution is

$$\boxed{r_{\mathrm{char},i,a}^{(e),\Omega}=-\frac{\Delta t_e}{2}\int_{\Omega_e}N_{a,k}\,\bar u_k\,\frac{\partial F_{ij}}{\partial x_j}\,d\Omega}.$$

The element-average velocity is

$$\bar u_k=\frac{1}{4}\sum_{a=1}^{4}u_{k,a}.$$

Because the adopted gradients and averages are constant in a P1 element,

$$r_{\mathrm{char},i,a}^{(e),\Omega}=-\frac{\Delta t_e}{2}V_eN_{a,k}^{(e)}\bar u_k\frac{\partial F_{ij}}{\partial x_j}.$$

The code factor is therefore

$$\texttt{ldelte}=\frac{1}{2}\Delta t_e\det J_e\left(\frac{1}{6}\right)=\frac{1}{2}\Delta t_eV_e.$$

This maps to `umean`, `lunksum`, `lunksumk`, `ldelte` and `step1_characteristic_correction(...)`.

### 5.5 Viscous diffusion

The constant element velocity gradient is

$$\left.\frac{\partial u_i}{\partial x_j}\right|_e=\sum_{b=1}^{4}u_{i,b}N_{b,j}^{(e)}.$$

The implemented weak volume residual is

$$\boxed{r_{\mathrm{diff},i,a}^{(e),\Omega}=-\nu_eV_eN_{a,j}^{(e)}\left.\frac{\partial u_i}{\partial x_j}\right|_e}.$$

The diffusivity is

$$\nu_e=\begin{cases}\mu_e/\rho_e,&\text{dimensional material mode},\\\texttt{ani},&\text{non-dimensional mode}.\end{cases}$$

This maps to `dNuidxj`, `momentum_diffusivity(...)`, `fdif[1]` and `step1_diffusion(...)`.

### 5.6 Step-1 code contract

| Item | Meaning |
|---|---|
| routine | `MomentumAssembly::assembleStep1Rhs` |
| domain | fluid elements only: `mat_elem(e) == 0` |
| reads | `unkn1`, `intma`, `dNkdx`, `detJ`, `delte`, face geometry and material data |
| writes | `rhs` |
| update routine | `Steps::step1SemiImplicit` |
| post-update action | reapply symmetry, velocity/wall, backflow and solid-zero constraints |

Element colouring is an execution strategy for race-free scatter; it does not change the finite-element operator.

## 6. CBS Step 2: pressure equation

### 6.1 Projection equation

The pressure field is chosen so that the Step-3 corrected velocity satisfies continuity. The active algebraic system is

$$\boxed{A_p\mathbf{p}=\mathbf{b}_p}.$$

The local time scale is carried by $A_p$, while $\mathbf{b}_p$ is the raw weak divergence of the predictor.

### 6.2 Element pressure stiffness

The pressure stiffness matrix is

$$H_{ab}^{(e)}=\int_{\Omega_e}\nabla N_a\!\cdot\nabla N_b\,d\Omega.$$

Since the P1 gradients are constant,

$$\boxed{H_{ab}^{(e)}=V_eN_{a,i}^{(e)}N_{b,i}^{(e)}}.$$

`PressureAssembly::buildElementPressureTerms()` stores the four diagonal coefficients in `pdiagE` and the six unique off-diagonal coefficients in `gstifE`, ordered as $(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)$.

The time-scaled element operator is

$$\boxed{A_p^{(e)}=\Delta t_eH^{(e)}}.$$

`TimeStep::updateLhsDiagonal()` forms the active `pdiag` and `gstif` arrays from `pdiagE` and `gstifE`.

### 6.3 Pressure right-hand side

The element volume term is

$$b_{p,a}^{(e),\Omega}=\int_{\Omega_e}\nabla N_a\!\cdot\mathbf{u}^*\,d\Omega.$$

With linear velocity interpolation,

$$\boxed{b_{p,a}^{(e),\Omega}=N_{a,i}^{(e)}\frac{V_e}{4}\sum_{b=1}^{4}u_{i,b}^*}.$$

This maps to `u_sum`, `vol4`, `lrhs` and the assembled `rhs1` in `PressureAssembly::assembleStep2Rhs()`.

> **Time-step convention.** `rhs1` is not divided by `dtreal`. The time factor is already embedded in $A_p^{(e)}=\Delta t_eH^{(e)}$ and in the Step-3 diagonal update. Moving this factor independently changes the numerical method.

### 6.4 Pressure space and constraints

A pressure node is free when it belongs to at least one fluid element and is not prescribed:

$$\operatorname{free}(I)=\operatorname{active}(I)\land\neg\operatorname{fixed}(I).$$

Solid-only nodes are excluded from the pressure space. Prescribed pressure nodes are imposed algebraically. A reference pressure or equivalent null-space treatment is required for an all-Neumann pressure boundary configuration.

### 6.5 Matrix-free operator

The native solver applies

$$\boxed{\mathbf{y}=A_p\mathbf{x}}$$

without constructing a global CSR matrix. Element off-diagonal actions use `gstif`; the assembled diagonal action uses `pdiag`. The compact pair order must be identical in pressure assembly and `MatrixVectorCalc::pressureMultiply()`.

### 6.6 Step-2 theory-to-code map

| Theory object | C++ representation |
|---|---|
| $H_{aa}^{(e)}$ | `pdiagE` |
| six $H_{ab}^{(e)}$, $a<b$ | `gstifE` |
| $A_p^{(e)}=\Delta t_eH^{(e)}$ | `gstif` and assembled `pdiag` |
| $\mathbf{b}_p$ | `rhs1` during Step 2 |
| $A_p\mathbf{x}$ | matrix-free pressure multiplication |
| pressure solution | `pres` |
| solver selection | native PCG or PETSc in `Steps::step2SemiImplicit` |

## 7. Native preconditioned conjugate gradient

Given an initial pressure estimate $\mathbf{p}_0$,

$$\mathbf{r}_0=\mathbf{b}_p-A_p\mathbf{p}_0,\qquad M\mathbf{z}_0=\mathbf{r}_0,\qquad\mathbf{d}_0=\mathbf{z}_0.$$

At iteration $k$,

$$\alpha_k=\frac{\mathbf{r}_k^T\mathbf{z}_k}{\mathbf{d}_k^TA_p\mathbf{d}_k},\qquad\mathbf{p}_{k+1}=\mathbf{p}_k+\alpha_k\mathbf{d}_k.$$

$$\mathbf{r}_{k+1}=\mathbf{r}_k-\alpha_kA_p\mathbf{d}_k,\qquad\beta_k=\frac{\mathbf{r}_{k+1}^T\mathbf{z}_{k+1}}{\mathbf{r}_k^T\mathbf{z}_k},\qquad\mathbf{d}_{k+1}=\mathbf{z}_{k+1}+\beta_k\mathbf{d}_k.$$

For Jacobi preconditioning,

$$M=\operatorname{diag}(A_p),\qquad z_I=\frac{r_I}{(A_p)_{II}}.$$

The diagonal is exactly `pdiag(I)`. Dot products, norms and updates operate only on free pressure nodes. Before pressure constraints, the Laplacian is positive semidefinite; removing the constant null mode gives the SPD free problem required by CG.

## 8. CBS Step 3: velocity correction

### 8.1 Element pressure gradient

For each fluid tetrahedron,

$$\boxed{\left.\frac{\partial p}{\partial x_i}\right|_e=\sum_{a=1}^{4}p_aN_{a,i}^{(e)}}.$$

The constant gradient is accumulated into `grad_pres[idim]`.

### 8.2 Local residual and nodal update

The P1 nodal pressure-gradient contribution is

$$\boxed{r_{p,i,a}^{(e)}=-\frac{V_e}{4}\left.\frac{\partial p}{\partial x_i}\right|_e}.$$

It is scattered to `rhs(i,I)`, after which

$$\boxed{\mathbf{u}^{n+1}=\mathbf{u}^*+D_u^{-1}\mathbf{r}_p}.$$

The required time factor is already contained in `elcoe2`; Step 3 must not multiply by a second explicit $\Delta t$.

```cpp
s.unkno(idim, ip) +=
    s.rhs(idim, ip) * s.elcoe2(ip);
```

Velocity boundary constraints are reapplied after correction because the algebraic update modifies every free nodal component.

## 9. CBS Step 4: temperature and conjugate heat transfer

### 9.1 Nodal equation

The thermal residual is decomposed as

$$\mathbf{r}_T=\mathbf{r}_{\mathrm{conv}}+\mathbf{r}_{\mathrm{stab}}+\mathbf{r}_{\mathrm{diff}}+\mathbf{r}_{\mathrm{source}}+\mathbf{r}_{\mathrm{flux}}.$$

The update is

$$\boxed{T^{n+1}=T^n+D_T^{-1}\mathbf{r}_T}.$$

`temperature1` supplies $T^n$, `rhs1` stores $\mathbf{r}_T$, and `elcoe2p` stores $D_T^{-1}$.

### 9.2 Temperature gradient

$$\boxed{\left.\frac{\partial T}{\partial x_i}\right|_e=\sum_{a=1}^{4}T_a^nN_{a,i}^{(e)}}.$$

This is computed by `compute_temperature_gradient(...)` and stored temporarily as `dTdx`, `dTdy` and `dTdz`.

### 9.3 Fluid convection

The weak thermal-convection term is

$$r_{\mathrm{conv},a}^{(e)}=-(\rho c_p)_e\int_{\Omega_e}N_a\,\mathbf{u}^{n+1}\!\cdot\nabla T^n\,d\Omega.$$

Using the consistent P1 mass integral,

$$\boxed{r_{\mathrm{conv},a}^{(e)}=-(\rho c_p)_e\frac{V_e}{20}\left(\sum_{b=1}^{4}\mathbf{u}_b^{n+1}+\mathbf{u}_a^{n+1}\right)\!\cdot\nabla T^n}.$$

This is implemented by `add_fluid_convection(...)` and uses the corrected velocity in `unkno`.

### 9.4 Thermal characteristic stabilisation

The scalar streamline correction is

$$\boxed{r_{\mathrm{stab},a}^{(e)}=\frac{\Delta t_e}{2}(\rho c_p)_eV_e\left(\bar{\mathbf{u}}\!\cdot\nabla N_a\right)\left(\bar{\mathbf{u}}\!\cdot\nabla T^n\right)}.$$

It is assembled only in fluid elements by `add_fluid_convection_stabilisation(...)`.

### 9.5 Diffusion, source and prescribed heat flux

The volume diffusion residual is

$$\boxed{r_{\mathrm{diff},a}^{(e)}=-k_eV_e\nabla N_a\!\cdot\nabla T^n}.$$

It is assembled in both fluid and solid. For a constant volumetric source,

$$\boxed{r_{\mathrm{source},a}^{(e)}=Q_e\frac{V_e}{4}}.$$

For an inward prescribed heat flux $q''$ on a triangular face,

$$\boxed{r_{\mathrm{flux},a}^{(f)}=q''\frac{A_f}{3}}.$$

The conformal fluid-solid interface receives no separate imposed flux term; shared temperature degrees of freedom and diffusion on both sides provide the discrete coupling.

### 9.6 Step-4 theory-to-code map

| Theory object | C++ representation |
|---|---|
| $T^n$ | `temperature1` |
| $\nabla T^n$ | `dTdx`, `dTdy`, `dTdz` |
| corrected $\mathbf{u}^{n+1}$ | `unkno` after Step 3 |
| $(\rho c_p)_e$ | `rho_cp_e` |
| $k_e$ or $k_{\mathrm{eff},e}$ | `k_e`, `k_eff_e` |
| $Q_e$ | `Qvol_e` |
| $\mathbf{r}_T$ | `rhs1` during Step 4 |
| $D_T^{-1}$ | `elcoe2p` |
| $T^{n+1}$ | `temperature` |

## 10. Optional Spalart-Allmaras transport

The transported SA working variable is $\tilde\nu$, not the eddy viscosity itself. The constitutive conversion is

$$\chi=\frac{\tilde\nu}{\nu},\qquad f_{v1}=\frac{\chi^3}{\chi^3+c_{v1}^3},\qquad\nu_t=\tilde\nu f_{v1},\qquad\mu_t=\rho\nu_t.$$

The standard transport architecture is

$$\frac{\partial\tilde\nu}{\partial t}+u_j\frac{\partial\tilde\nu}{\partial x_j}=\text{production}-\text{destruction}+\text{diffusion}+\text{nonlinear gradient source}.$$

The current coupling order is

```text
Step 3 corrected velocity
        ↓
SA residual and nodal update
        ↓
SA wall/inlet conditions
        ↓
nu_t, mu_t, mu_eff and optional k_eff refresh
        ↓
Step 4 temperature assembly
```

SA is an additional transported scalar, not a fifth CBS split. Wall distance, near-wall resolution, inlet $\tilde\nu$, source stiffness and SA/SA-neg behaviour require independent validation.

## 11. Residual ownership by stage

| Stage | Mathematical residual | Work array | Resulting field |
|---|---|---|---|
| Step 1 | $\mathbf{r}_m$ | `rhs` | predictor `unkno = u*` |
| Step 2 | $\mathbf{b}_p$ | `rhs1` | `pres` |
| Step 3 | $\mathbf{r}_p$ | `rhs` | corrected `unkno = u^(n+1)` |
| SA | $r_{SA}$ | SA-specific storage | `nu_tilde` and effective properties |
| Step 4 | $\mathbf{r}_T$ | `rhs1` | `temperature` |

The reuse is safe only because each assembly routine clears the relevant work array before accumulation.

## 12. Verification checklist

### Geometry and P1 operators

1. Reproduce $V_e$, $N_{a,i}^{(e)}$ and $V_e/4$ for one analytical tetrahedron.
2. Confirm $\sum_aN_a=1$ and $\sum_a\nabla N_a=\mathbf{0}$.
3. Confirm positive `detJ` and consistent local-node orientation.
4. Compare flattened `dNkdx` lookup with an explicit $(e,i,a)$ table.

### Step 1

1. Uniform velocity must give zero interior velocity gradient and zero viscous volume residual.
2. Reproduce one-element Galerkin convection and characteristic correction by hand.
3. Verify `elcoe2` is inverse mass-over-time at the instant of the predictor update.
4. Confirm essential velocity conditions after Step 1.
5. Keep SA disabled during baseline laminar verification.

### Pressure and Step 3

1. Confirm $H^{(e)}$ is symmetric and has zero row sum before constraints.
2. Compare compact matrix-free $A_p\mathbf{x}$ with an explicitly assembled small matrix.
3. Verify the six off-diagonal pair order in both assembly and multiplication.
4. Confirm one pressure reference removes the constant null mode.
5. Check $\mathbf{x}^TA_p\mathbf{x}>0$ for non-zero free vectors after constraints.
6. Confirm divergence decreases after Step 3.

### Temperature and multiphysics

1. Uniform temperature must give zero gradient, convection and diffusion.
2. A linear temperature field must reproduce the exact constant P1 gradient.
3. Solid elements must have diffusion but no convection.
4. Interface nodes must share one temperature degree of freedom.
5. Verify $q''A_f/3$ at every node of a heat-flux face.
6. Check global thermal energy balance.

## 13. Source-file map

| Numerical responsibility | Primary source |
|---|---|
| Step-1 residual | `MomentumAssembly.cpp` |
| pressure stiffness and RHS | `PressureAssembly.cpp` |
| Step order and Step-3 correction | `Steps.cpp` |
| time-scaled diagonals/operator | `TimeStep.cpp` |
| native pressure multiplication | `MatrixVectorCalc.cpp` |
| native PCG | `ConjugateGradient.cpp` |
| PETSc pressure option | `PetscPressureSolver.cpp` |
| temperature residual | `EnergyAssembly.cpp` |
| SA transport | `SpalartAllmarasAssembly.cpp` |

## 14. References

1. P. Nithiarasu, R. Codina and O. C. Zienkiewicz, “The Characteristic-Based Split scheme—a unified approach to fluid dynamics,” *International Journal for Numerical Methods in Engineering*, 66, 1514–1546, 2006.
2. O. C. Zienkiewicz, R. L. Taylor and P. Nithiarasu, *The Finite Element Method for Fluid Dynamics*.
3. P. R. Spalart and S. R. Allmaras, “A One-Equation Turbulence Model for Aerodynamic Flows,” 1994.
4. CBS3D++ implementation files listed above.
