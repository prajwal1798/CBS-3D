# Mathematical formulation and code architecture

This page connects the **continuous equations**, the **weak finite-element form**, the **element algebra**, and the exact **CBS3D++ C++ storage/update path**.

The governing rule is:

$$
\boxed{
\text{strong PDE}
\;\longrightarrow\;
\text{weak form}
\;\longrightarrow\;
\text{element residual}
\;\longrightarrow\;
\text{assembled arrays and nodal update}
}
$$

A numerical operation is not considered fully documented until all four levels can be identified.

!!! important "Implemented formulation"
    CBS3D++ currently implements a **semi-implicit incompressible pressure-projection CBS specialisation** on four-node linear tetrahedra. General CBS equations are included to establish the derivation, but every section distinguishes the general theory from the active C++ implementation.

## 1. Scope and conventions

The active solver path supports:

- three-dimensional incompressible flow;
- P1 tetrahedral finite elements;
- element-wise constant material properties;
- lumped nodal mass/time diagonals;
- CBS momentum prediction, pressure projection and velocity correction;
- optional temperature transport in fluid and solid regions;
- native matrix-free PCG or the current PETSc pressure path;
- OpenMP-coloured element assembly.

True physical-time BDF or dual-time integration is not yet a released capability. The active iteration is steady or pseudo-transient.

### 1.1 Indices

Cartesian indices use

$$
 i,j,k\in\{1,2,3\},
 \qquad
 x_1=x,\;x_2=y,\;x_3=z.
$$

Local tetrahedral nodes use $a,b\in\{1,2,3,4\}$, elements use $e$, and global nodes use $I,J$.

### 1.2 Theory variables versus stored variables

The conservative CBS literature often introduces momentum per unit volume

$$
U_i=\rho u_i.
$$

The active C++ arrays `unkn1` and `unkno` store **velocity components**, not a separate conservative-momentum field. Consequently, code-aligned equations below are written in velocity form after the momentum equation has been divided by density where required.

| Mathematical meaning | C++ storage |
|---|---|
| coordinates $x_i$ | `coord(i,I)` |
| local-to-global map $I=I_e(a)$ | `intma(a,e)` |
| shape gradients $N_{a,i}^{(e)}$ | flattened `dNkdx` |
| Jacobian determinant $\det J_e$ | `detJ(e)` |
| previous velocity $u_i^n$ | `unkn1(i,I)` |
| predictor/corrected velocity | `unkno(i,I)` |
| pressure $p$ | `pres(I)` |
| previous temperature $T^n$ | `temperature1(I)` |
| updated temperature | `temperature(I)` |
| Step-1 or Step-3 vector residual | `rhs(i,I)` |
| Step-2 pressure or Step-4 thermal residual | `rhs1(I)` |
| local element time scale $\Delta t_e$ | `delte(e)` |
| inverse fluid mass/time diagonal $D_u^{-1}$ | `elcoe2(I)` |
| inverse thermal capacitance/time diagonal $D_T^{-1}$ | `elcoe2p(I)` |

The work arrays `rhs` and `rhs1` are deliberately reused. Their mathematical interpretation therefore depends on the active CBS step.

## 2. Governing equations

### 2.1 Incompressible continuity

$$
\boxed{
\frac{\partial u_i}{\partial x_i}=0
}
$$

or, in vector form,

$$
\nabla\!\cdot\mathbf{u}=0.
$$

### 2.2 Momentum

The incompressible Newtonian momentum equation is

$$
\rho
\left(
\frac{\partial u_i}{\partial t}
+u_j\frac{\partial u_i}{\partial x_j}
\right)
=
-\frac{\partial p}{\partial x_i}
+\frac{\partial \tau_{ij}}{\partial x_j}
+\rho b_i.
$$

The complete Newtonian deviatoric stress is

$$
\tau_{ij}
=
\mu
\left(
\frac{\partial u_i}{\partial x_j}
+
\frac{\partial u_j}{\partial x_i}
-
\frac{2}{3}\delta_{ij}
\frac{\partial u_k}{\partial x_k}
\right).
$$

For an exactly divergence-free field, the final isotropic term vanishes. However, the Step-1 predictor is not pointwise divergence-free before pressure correction.

!!! warning "Active diffusion operator"
    `MomentumAssembly.cpp` currently assembles a **component-wise velocity Laplacian**,

    $$
    \nu\,\nabla^2 u_i,
    $$

    through $\nabla N_a\cdot\nabla u_i$. It does not currently assemble the complete symmetric-gradient/deviatoric-stress operator. The documentation follows the active code rather than silently substituting the more general stress form.

Body force is part of the governing equation but is not assembled in the current Step-1 core.

### 2.3 Temperature equation

For fluid elements,

$$
\boxed{
\rho c_p
\left(
\frac{\partial T}{\partial t}
+u_i\frac{\partial T}{\partial x_i}
\right)
=
\frac{\partial}{\partial x_i}
\left(
 k\frac{\partial T}{\partial x_i}
\right)
+Q
}
$$

For solid elements, $\mathbf{u}=\mathbf{0}$ and the equation reduces to transient/pseudo-transient heat conduction:

$$
\rho c_p\frac{\partial T}{\partial t}
=
\nabla\!\cdot(k\nabla T)+Q.
$$

A conformal fluid-solid interface shares the same nodal temperature unknown. Diffusion is assembled from both adjacent materials; no duplicate interface temperature is introduced.

## 3. P1 tetrahedral finite-element architecture

### 3.1 Interpolation

For one four-node tetrahedron $\Omega_e$,

$$
 u_i^h(\mathbf{x})
 =\sum_{a=1}^{4}N_a(\mathbf{x})u_{i,a},
 \qquad
 p^h(\mathbf{x})
 =\sum_{a=1}^{4}N_a(\mathbf{x})p_a,
$$

$$
 T^h(\mathbf{x})
 =\sum_{a=1}^{4}N_a(\mathbf{x})T_a.
$$

The physical map is

$$
\mathbf{x}
=
\mathbf{x}_1
+J_e
\begin{bmatrix}
\xi\\[2pt]\eta\\[2pt]\zeta
\end{bmatrix},
$$

with

$$
J_e=
\begin{bmatrix}
 x_2-x_1 & x_3-x_1 & x_4-x_1\\
 y_2-y_1 & y_3-y_1 & y_4-y_1\\
 z_2-z_1 & z_3-z_1 & z_4-z_1
\end{bmatrix}.
$$

The element volume is

$$
\boxed{
V_e=\frac{\det J_e}{6}
}
$$

for the positive-orientation convention enforced by preprocessing.

### 3.2 Constant shape gradients

For a P1 tetrahedron,

$$
N_{a,i}^{(e)}
\equiv
\frac{\partial N_a}{\partial x_i}
=
\text{constant in }\Omega_e.
$$

Therefore,

$$
\frac{\partial u_i^h}{\partial x_j}
=
\sum_{a=1}^{4}u_{i,a}N_{a,j}^{(e)},
\qquad
\frac{\partial p^h}{\partial x_i}
=
\sum_{a=1}^{4}p_aN_{a,i}^{(e)},
$$

and both gradients are constant inside one element.

The flattened storage address is

$$
\mathcal{I}(e,i,a)
=
(e-1)(3)(4)+(i-1)(4)+a.
$$

This is implemented by `dNkdx_index(...)`; the values themselves reside in `s.dNkdx`.

### 3.3 Exact integration identities

The constants used throughout the code follow from exact P1 integration:

$$
\int_{\Omega_e}N_a\,d\Omega
=\frac{V_e}{4}
=\frac{\det J_e}{24},
$$

$$
\int_{\Omega_e}N_a^2\,d\Omega
=\frac{V_e}{10},
\qquad
\int_{\Omega_e}N_aN_b\,d\Omega
=\frac{V_e}{20}
\quad(a\ne b),
$$

$$
\int_{\Gamma_f}N_a\,d\Gamma
=\frac{A_f}{3},
$$

$$
\int_{\Gamma_f}N_a^2\,d\Gamma
=\frac{A_f}{6},
\qquad
\int_{\Gamma_f}N_aN_b\,d\Gamma
=\frac{A_f}{12}
\quad(a\ne b).
$$

These identities explain the principal code factors:

| Code quantity | Mathematical value |
|---|---:|
| `detJ * fcon[1]` | $\det J_e/24=V_e/4$ |
| `fcon[2]` | $1/12$ for consistent triangular-face integration |
| `detJ * fdif[1]` | $\det J_e/6=V_e$ |
| `fdif[2]` | $1/3$ for $\int_{\Gamma_f}N_a\,d\Gamma$ |

### 3.4 Consistent and lumped mass

The consistent scalar element mass matrix is

$$
M^{(e)}
=
\frac{V_e}{20}
\begin{bmatrix}
2&1&1&1\\
1&2&1&1\\
1&1&2&1\\
1&1&1&2
\end{bmatrix}.
$$

Row-sum lumping gives

$$
\boxed{
M_L^{(e)}=\frac{V_e}{4}I_4
}
$$

and preprocessing stores

$$
\texttt{elcoe\_e}(e,a)=\frac{V_e}{4}.
$$

The active momentum diagonal is not merely $M_L^{-1}$. The local element time scale is already incorporated:

$$
[D_u]_{II}
=
\sum_{e\ni I,\;e\in\Omega_f}
\frac{V_e/4}{\Delta t_e},
$$

$$
\boxed{
\texttt{elcoe2}(I)=[D_u]_{II}^{-1}
}
$$

Only fluid elements contribute. For temperature,

$$
[D_T]_{II}
=
\sum_{e\ni I}
\frac{(\rho c_p)_eV_e/4}{\Delta t_e},
$$

$$
\boxed{
\texttt{elcoe2p}(I)=[D_T]_{II}^{-1}
}
$$

and both fluid and solid elements contribute.

## 4. Global CBS data flow

```text
coord, intma, iside
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

At the beginning of one CBS iteration,

$$
\texttt{unkn1}=\mathbf{u}^n,
\qquad
\texttt{temperature1}=T^n.
$$

After Step 1, `unkno` means $\mathbf{u}^*$. After Step 3, the same array means $\mathbf{u}^{n+1}$.

## 5. CBS Step 1: pressure-free momentum predictor

### 5.1 General CBS interpretation

Let the non-pressure momentum operator be

$$
R_i(\mathbf{u})
=
-\frac{\partial (u_ju_i)}{\partial x_j}
+\nu\nabla^2u_i
+b_i.
$$

A characteristic expansion gives a pressure-free increment of the form

$$
\Delta u_i^*
=
\Delta t\,R_i^n
-
\frac{\Delta t^2}{2}
\frac{\partial}{\partial x_k}
\left(u_kR_i^n\right)
+\cdots.
$$

The active implementation retains:

- Galerkin convection;
- the convective characteristic correction;
- component-wise viscous diffusion.

It does not currently retain body force, artificial diffusion, or characteristic derivatives of viscous/body-force terms.

### 5.2 Code-level nodal equation

The assembled Step-1 equation is

$$
\boxed{
D_u
\left(
\mathbf{u}^*-\mathbf{u}^n
\right)
=
\mathbf{r}_m
}
$$

with

$$
\mathbf{r}_m
=
\mathbf{r}_{\mathrm{conv}}
+
\mathbf{r}_{\mathrm{char}}
+
\mathbf{r}_{\mathrm{diff}}.
$$

The nodal update in `Steps::step1SemiImplicit()` is

$$
\boxed{
\mathbf{u}^*
=
\mathbf{u}^n
+D_u^{-1}\mathbf{r}_m
}
$$

which maps directly to

```cpp
s.unkno(idim, ip) =
    s.unkn1(idim, ip)
    + s.rhs(idim, ip) * s.elcoe2(ip);
```

### 5.3 Galerkin convection

Define the nodal momentum-flux product

$$
F_{ij,a}=u_{j,a}u_{i,a}.
$$

The strong convective contribution to the Step-1 right-hand side is

$$
-\frac{\partial F_{ij}}{\partial x_j}.
$$

Multiplying by test function $N_a$ and integrating by parts gives

$$
r_{\mathrm{conv},i,a}^{(e)}
=
\int_{\Omega_e}
N_{a,j}F_{ij}\,d\Omega
-
\int_{\Gamma_e}
N_aF_{ij}n_j\,d\Gamma.
$$

With $F_{ij}^h=N_bF_{ij,b}$ and constant $N_{a,j}$,

$$
\boxed{
r_{\mathrm{conv},i,a}^{(e),\Omega}
=
N_{a,j}^{(e)}
\frac{V_e}{4}
\sum_{b=1}^{4}F_{ij,b}
}
$$

This is the exact meaning of `lunksum[i][j]` and the code line

```cpp
lrhs[i][a] +=
    grad(s, ie, j, a)
    * volume_factor
    * lunksum[i][j];
```

where `volume_factor = detJ/24 = V_e/4`.

For a triangular face, consistent interpolation produces the $A_f/12$ factor used with `annxf` and `fcon[2]`.

### 5.4 Characteristic correction

The implemented volume correction is

$$
\boxed{
r_{\mathrm{char},i,a}^{(e),\Omega}
=
-\frac{\Delta t_e}{2}
\int_{\Omega_e}
N_{a,k}\,\bar u_k
\frac{\partial F_{ij}}{\partial x_j}
\,d\Omega
}
$$

where

$$
\bar u_k
=
\frac{1}{4}
\sum_{a=1}^{4}u_{k,a}.
$$

Using constant P1 gradients and the element-average velocity,

$$
r_{\mathrm{char},i,a}^{(e),\Omega}
=
-\frac{\Delta t_e}{2}
V_e
N_{a,k}^{(e)}\bar u_k
\frac{\partial F_{ij}}{\partial x_j}.
$$

The code factor

$$
\texttt{ldelte}
=
\frac{1}{2}\Delta t_e\det J_e\left(\frac{1}{6}\right)
=
\frac{1}{2}\Delta t_eV_e
$$

is therefore not arbitrary; it is the exact constant-element integral.

### 5.5 Viscous diffusion

The active velocity gradient is

$$
\left.\frac{\partial u_i}{\partial x_j}\right|_e
=
\sum_{b=1}^{4}u_{i,b}N_{b,j}^{(e)}.
$$

The implemented weak volume term is

$$
\boxed{
r_{\mathrm{diff},i,a}^{(e),\Omega}
=
-\nu_eV_e
N_{a,j}^{(e)}
\frac{\partial u_i}{\partial x_j}
}
$$

with

$$
\nu_e=
\begin{cases}
\mu_e/\rho_e, & \text{dimensional material mode},\\[4pt]
\texttt{ani}, & \text{non-dimensional mode}.
\end{cases}
$$

This maps to `dNuidxj`, `momentum_diffusivity(...)`, `fdif[1]`, and `step1_diffusion(...)`.

### 5.6 Step-1 theory-to-code map

| Theory object | Element/C++ representation |
|---|---|
| $u_i^n$ | gathered from `unkn1(i,I)` into `lunkno[i][a]` |
| $F_{ij,a}=u_{j,a}u_{i,a}$ | `lunk[i][j][a]` |
| $r_{\mathrm{conv}}^{(e)}$ | `step1_convective_galerkin(...)` |
| $r_{\mathrm{char}}^{(e)}$ | `step1_characteristic_correction(...)` |
| $r_{\mathrm{diff}}^{(e)}$ | `step1_diffusion(...)` |
| local residual $r_{m,i,a}^{(e)}$ | `lrhs[i][a]` |
| assembled $r_{m,i,I}$ | `rhs(i,I)` |
| $D_u^{-1}$ | `elcoe2(I)` |
| $u_i^*$ | `unkno(i,I)` after Step 1 |

Element colouring changes only execution order and race avoidance; it does not alter the finite-element assembly operator.

## 6. CBS Step 2: pressure equation

### 6.1 From continuity to a Poisson system

The pressure correction is chosen so that the corrected velocity satisfies continuity. In the active algebraic convention, the code solves

$$
\boxed{
A_p\mathbf{p}=\mathbf{b}_p
}
$$

where the local time factor is stored in the pressure operator rather than dividing the right-hand side by a global time step.

### 6.2 Element pressure stiffness

The pressure stiffness matrix is

$$
H_{ab}^{(e)}
=
\int_{\Omega_e}
\nabla N_a\cdot\nabla N_b
\,d\Omega.
$$

Because the gradients are constant,

$$
\boxed{
H_{ab}^{(e)}
=
V_e
N_{a,i}^{(e)}N_{b,i}^{(e)}
}
$$

`PressureAssembly::buildElementPressureTerms()` stores:

- the four diagonal entries in `pdiagE`;
- the six unique off-diagonal entries in `gstifE`;
- the pair order $(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)$.

The time-scaled element operator is

$$
\boxed{
A_p^{(e)}
=
\Delta t_eH^{(e)}
}
$$

and `TimeStep::updateLhsDiagonal()` forms `pdiag` and `gstif` from `pdiagE` and `gstifE`.

### 6.3 Pressure right-hand side

The element volume contribution is the weak divergence of the predictor:

$$
b_{p,a}^{(e),\Omega}
=
\int_{\Omega_e}
\nabla N_a\cdot\mathbf{u}^*
\,d\Omega.
$$

With linear velocity interpolation,

$$
\boxed{
b_{p,a}^{(e),\Omega}
=
N_{a,i}^{(e)}
\frac{V_e}{4}
\sum_{b=1}^{4}u_{i,b}^*
}
$$

This maps directly to `u_sum`, `vol4`, `lrhs`, and finally `rhs1` in `PressureAssembly::assembleStep2Rhs()`.

!!! important "Time-step convention"
    The pressure RHS is the **raw weak divergence**. It is not divided by `dtreal` in `assembleStep2Rhs()`. The time scale is already carried by $A_p^{(e)}=\Delta t_eH^{(e)}$ and by the Step-3 mass/time diagonal. Moving the factor to the RHS without changing the other stages would change the method.

### 6.4 Pressure constraints

Pressure is active at a node if the node belongs to at least one fluid element. A pressure node is free when it is active and not prescribed:

$$
\mathrm{free}(I)
=
\mathrm{active}(I)
\land
\neg\mathrm{fixed}(I).
$$

Solid-only nodes are excluded from the pressure space. Prescribed pressure nodes are enforced algebraically. A pressure reference or equivalent null-space treatment is required when all pressure boundaries are Neumann.

### 6.5 Matrix-free multiplication

For the native solver, a global CSR matrix is not required. The off-diagonal element contributions are scattered, followed by the assembled diagonal action:

$$
\boxed{
\mathbf{y}=A_p\mathbf{x}
}
$$

with `gstif` providing the six compact off-diagonal pairs and `pdiag` providing the global diagonal. The same pair ordering must be used by assembly and `MatrixVectorCalc::pressureMultiply()`.

### 6.6 Step-2 theory-to-code map

| Theory object | C++ representation |
|---|---|
| $H_{aa}^{(e)}$ | `pdiagE` |
| six $H_{ab}^{(e)}$, $a<b$ | `gstifE` |
| $A_p^{(e)}=\Delta t_eH^{(e)}$ | `gstif` plus assembled `pdiag` |
| $\mathbf{b}_p$ | `rhs1` during Step 2 |
| pressure-active/fixed masks | pressure-space classification in solver state |
| $A_p\mathbf{x}$ | matrix-free pressure multiplication |
| $\mathbf{p}$ | `pres` |

## 7. Native preconditioned conjugate gradient

For an initial pressure guess $\mathbf{p}_0$,

$$
\mathbf{r}_0
=
\mathbf{b}_p-A_p\mathbf{p}_0,
$$

$$
M\mathbf{z}_0=\mathbf{r}_0,
\qquad
\mathbf{d}_0=\mathbf{z}_0.
$$

At iteration $k$,

$$
\alpha_k
=
\frac{\mathbf{r}_k^T\mathbf{z}_k}
     {\mathbf{d}_k^TA_p\mathbf{d}_k},
$$

$$
\mathbf{p}_{k+1}
=
\mathbf{p}_k+\alpha_k\mathbf{d}_k,
$$

$$
\mathbf{r}_{k+1}
=
\mathbf{r}_k-\alpha_kA_p\mathbf{d}_k,
$$

$$
\beta_k
=
\frac{\mathbf{r}_{k+1}^T\mathbf{z}_{k+1}}
     {\mathbf{r}_k^T\mathbf{z}_k},
$$

$$
\mathbf{d}_{k+1}
=
\mathbf{z}_{k+1}+\beta_k\mathbf{d}_k.
$$

For Jacobi preconditioning,

$$
M=\operatorname{diag}(A_p),
\qquad
z_I=\frac{r_I}{(A_p)_{II}}.
$$

The diagonal is exactly `pdiag(I)`. Dot products and norms must include only free pressure nodes.

## 8. CBS Step 3: velocity correction

### 8.1 Element pressure gradient

For each fluid tetrahedron,

$$
\boxed{
\left.\frac{\partial p}{\partial x_i}\right|_e
=
\sum_{a=1}^{4}p_aN_{a,i}^{(e)}
}
$$

The gradient is stored temporarily in `grad_pres[idim]`.

### 8.2 Local pressure residual

The weak nodal pressure-gradient contribution is

$$
r_{p,i,a}^{(e)}
=
-\int_{\Omega_e}
N_a\frac{\partial p}{\partial x_i}
\,d\Omega.
$$

Because $\nabla p$ is constant,

$$
\boxed{
r_{p,i,a}^{(e)}
=
-\frac{V_e}{4}
\left.\frac{\partial p}{\partial x_i}\right|_e
}
$$

which is implemented as

```cpp
s.rhs(idim, ip) -= vol4 * grad_pres[idim];
```

### 8.3 Nodal velocity correction

The code-aligned correction is

$$
\boxed{
\mathbf{u}^{n+1}
=
\mathbf{u}^*
+D_u^{-1}\mathbf{r}_p
}
$$

or

```cpp
s.unkno(idim, ip) +=
    s.rhs(idim, ip) * s.elcoe2(ip);
```

Do not append another explicit $\Delta t$ to this code equation: the required time scaling is already embedded in `elcoe2` through $D_u^{-1}$.

After correction, symmetry, prescribed velocity, wall conditions, outlet backflow control and zero solid velocity are reapplied. This is part of the algebraic enforcement of essential boundary conditions.

## 9. Step 4: temperature equation

### 9.1 Code-level nodal equation

The thermal update is

$$
\boxed{
D_T
\left(
T^{n+1}-T^n
\right)
=
r_T
}
$$

with

$$
r_T
=
r_{\mathrm{conv}}
+r_{\mathrm{stab}}
+r_{\mathrm{diff}}
+r_{\mathrm{source}}
+r_{\mathrm{flux}}.
$$

The nodal update is

$$
\boxed{
T^{n+1}
=
T^n+D_T^{-1}r_T
}
$$

and maps to

```cpp
s.temperature(ip) =
    s.temperature1(ip)
    + s.rhs1(ip) * s.elcoe2p(ip);
```

### 9.2 Element temperature gradient

$$
\boxed{
T_{,i}^{(e)}
=
\sum_{a=1}^{4}T_a^nN_{a,i}^{(e)}
}
$$

This is calculated from `temperature1` by `compute_temperature_gradient(...)`.

### 9.3 Fluid convection

The weak convective residual is

$$
r_{\mathrm{conv},a}^{(e)}
=
-(\rho c_p)_e
\int_{\Omega_e}
N_a\,u_iT_{,i}
\,d\Omega.
$$

Since $T_{,i}$ is constant and velocity is linearly interpolated,

$$
\boxed{
r_{\mathrm{conv},a}^{(e)}
=
-(\rho c_p)_e
\frac{V_e}{20}
\left(
\sum_{b=1}^{4}\mathbf{u}_b
+\mathbf{u}_a
\right)
\cdot\nabla T
}
$$

The corrected Step-3 velocity from `unkno` is used. This term is assembled only in fluid elements.

### 9.4 Thermal characteristic stabilisation

The scalar characteristic/SUPG-style term is

$$
\boxed{
r_{\mathrm{stab},a}^{(e)}
=
\frac{\Delta t_e}{2}
(\rho c_p)_eV_e
\left(
\bar{\mathbf{u}}\cdot\nabla N_a
\right)
\left(
\bar{\mathbf{u}}\cdot\nabla T
\right)
}
$$

where

$$
\bar{\mathbf{u}}
=
\frac{1}{4}
\sum_{a=1}^{4}\mathbf{u}_a.
$$

It is the scalar analogue of the Step-1 characteristic correction.

### 9.5 Thermal diffusion

After integration by parts, the element volume contribution is

$$
\boxed{
r_{\mathrm{diff},a}^{(e)}
=
-k_eV_e\nabla N_a\cdot\nabla T
}
$$

It is assembled in both fluid and solid elements. When turbulent heat transfer is enabled in a fluid element,

$$
k_{\mathrm{eff},e}
=
k_e
+(\rho c_p)_e\frac{\nu_{t,e}}{Pr_t}.
$$

The turbulent contribution is not applied to solid elements.

### 9.6 Volumetric source and surface heat flux

For constant volumetric generation $Q_e$,

$$
\boxed{
r_{\mathrm{source},a}^{(e)}
=
Q_e\frac{V_e}{4}
}
$$

For a prescribed inward-positive heat flux $q''$ on a triangular face,

$$
\boxed{
r_{\mathrm{flux},a}^{(f)}
=
q''\frac{A_f}{3}
}
$$

for each node on that face.

### 9.7 Step-4 theory-to-code map

| Theory object | C++ representation |
|---|---|
| $T^n$ | `temperature1` |
| $\nabla T$ | local `dTdx`, `dTdy`, `dTdz` |
| corrected $\mathbf{u}^{n+1}$ | `unkno` after Step 3 |
| $r_{\mathrm{conv}}$ | `add_fluid_convection(...)` |
| $r_{\mathrm{stab}}$ | `add_fluid_convection_stabilisation(...)` |
| $r_{\mathrm{diff}}$ | thermal diffusion kernel |
| $r_{\mathrm{source}}$ | volumetric-source kernel |
| $r_{\mathrm{flux}}$ | prescribed heat-flux kernel |
| assembled $r_T$ | `rhs1` during Step 4 |
| $D_T^{-1}$ | `elcoe2p` |
| $T^{n+1}$ | `temperature` |

## 10. Optional Spalart-Allmaras transport

The transported SA working variable is $\tilde\nu$, not the eddy viscosity itself. The standard model architecture is

$$
\frac{\partial\tilde\nu}{\partial t}
+u_j\frac{\partial\tilde\nu}{\partial x_j}
=
\underbrace{c_{b1}(1-f_{t2})\tilde S\tilde\nu}_{\text{production}}
-
\underbrace{
\left(
 c_{w1}f_w-\frac{c_{b1}}{\kappa^2}f_{t2}
\right)
\left(\frac{\tilde\nu}{d}\right)^2
}_{\text{destruction}}
+
\underbrace{
\frac{1}{\sigma}
\left[
\nabla\!\cdot
\left((\nu+\tilde\nu)\nabla\tilde\nu\right)
+c_{b2}|\nabla\tilde\nu|^2
\right]
}_{\text{diffusion and nonlinear source}}.
$$

The eddy viscosity is recovered from

$$
\nu_t=\tilde\nu f_{v1},
\qquad
f_{v1}=\frac{\chi^3}{\chi^3+c_{v1}^3},
\qquad
\chi=\frac{\tilde\nu}{\nu}.
$$

SA is an additional transported scalar executed after Step 3 and before Step 4; it is not a fifth CBS split. Its wall distance, near-wall resolution and source-term treatment require separate validation.

## 11. Boundary and material-domain semantics

The current material test is

$$
\texttt{mat\_elem}(e)=0
\quad\Longrightarrow\quad
\text{fluid},
$$

$$
\texttt{mat\_elem}(e)\ne0
\quad\Longrightarrow\quad
\text{non-fluid/solid in the current cases}.
$$

Consequences:

- momentum and pressure are assembled only in fluid elements;
- thermal diffusion and volumetric sources are assembled in fluid and solid;
- velocity is set to zero on every node of a non-fluid element;
- a conformal fluid-solid interface therefore receives no-slip velocity and one shared temperature degree of freedom;
- solid-only nodes are excluded from the pressure space.

Velocity constraints are reapplied after Steps 1 and 3 in this order:

1. symmetry projection;
2. prescribed velocity and wall conditions;
3. outlet backflow control;
4. zero velocity in the solid material.

Changing this order is a numerical change, not cosmetic refactoring.

## 12. General CBS versus the active implementation

| General formulation feature | Active CBS3D++ status |
|---|---|
| pressure-free predictor | implemented |
| Galerkin convection | implemented |
| second-order convective characteristic correction | implemented |
| complete characteristic derivative of all residual terms | not implemented |
| body force | not assembled in current Step 1 |
| complete Newtonian deviatoric-stress operator | not implemented; component Laplacian active |
| pressure-Poisson projection | implemented |
| higher-order characteristic pressure term in Step 3 | not implemented |
| true physical-time BDF/dual-time term | not implemented |
| temperature equation in fluid and solid | implemented development path |
| SA transported scalar | active development/validation |
| distributed Steps 1-4 | under development |

## 13. Verification checklist

### Geometry and interpolation

1. Verify $V_e=\det J_e/6$ on an analytical tetrahedron.
2. Verify $\sum_aN_a=1$ and $\sum_a\nabla N_a=\mathbf{0}$.
3. Verify `dNkdx_index(e,i,a)` addresses the intended derivative.
4. Verify a linear scalar field produces the exact constant element gradient.

### Step 1

1. Uniform velocity must give zero velocity gradient.
2. Verify one-element Galerkin convection against a hand calculation.
3. Verify characteristic signs and the $\Delta t_eV_e/2$ factor.
4. Verify `elcoe2` is the inverse mass/time diagonal at the instant of update.
5. Verify velocity constraints after the predictor.

### Pressure

1. Verify $H^{(e)}$ is symmetric.
2. Verify each unconstrained element stiffness row sums to zero.
3. Compare matrix-free $A_p\mathbf{x}$ with an explicitly assembled small matrix.
4. Verify the compact six-pair ordering in assembly and multiplication.
5. Verify one pressure reference removes the constant null space.
6. Verify solid-only nodes are excluded.

### Step 3

1. Verify the element pressure gradient for a linear pressure field.
2. Verify $r_{p,a}^{(e)}=-(V_e/4)\nabla p$.
3. Verify divergence is reduced after correction.
4. Verify no extra explicit time factor is applied outside `elcoe2`.

### Energy

1. Uniform temperature must produce zero gradient, convection and diffusion.
2. Solid elements must contain diffusion but no convection.
3. Verify $Q_eV_e/4$ for volumetric generation.
4. Verify $q''A_f/3$ for a prescribed triangular-face flux.
5. Verify one shared temperature field across a conformal interface.

## 14. Primary implementation files

| Numerical responsibility | Active source |
|---|---|
| geometry, shape gradients and lumped mass | `src/preprocess/Preprocess.cpp` |
| local time scales and time-dependent diagonals | `src/timestep/TimeStep.cpp` |
| Step-1 residual | `src/assembly/MomentumAssembly.cpp` |
| pressure stiffness and RHS | `src/assembly/PressureAssembly.cpp` |
| native pressure operator | `src/linalg/MatrixVectorCalc.cpp` |
| native PCG | `src/linalg/ConjugateGradient.cpp` |
| Steps 1-4 orchestration and Step-3 correction | `src/solver/Steps.cpp` |
| temperature residual | `src/assembly/EnergyAssembly.cpp` |
| SA residual and property coupling | `src/assembly/SpalartAllmarasAssembly.cpp` |

## 15. References

1. P. Nithiarasu, R. Codina and O. C. Zienkiewicz, “The Characteristic-Based Split scheme—a unified approach to fluid dynamics,” *International Journal for Numerical Methods in Engineering*, 66, 1514-1546, 2006.
2. O. C. Zienkiewicz, R. L. Taylor and P. Nithiarasu, *The Finite Element Method for Fluid Dynamics*, Elsevier.
3. P. R. Spalart and S. R. Allmaras, “A One-Equation Turbulence Model for Aerodynamic Flows,” 1994.
4. The active CBS3D++ source files listed above. Where a historical derivation and the current source disagree, the current source defines the documented implementation and the discrepancy must be resolved through code review and regression testing.
