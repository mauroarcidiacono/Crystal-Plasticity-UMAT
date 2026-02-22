# Variable Dictionary

This file documents the variables, arrays, flags, and constants used in the provided UMAT and its helper subroutines.  
Dimensions are given as `(rows, cols)` or `(length)`. All reals are double precision.

---

## 0) Conventions

* **Stress/strain ordering (Voigt):**  
  Normal: `1→11`, `2→22`, `3→33`; Shear: `4→12`, `5→13`, `6→23`.  
  UMAT passes `NTENS = NDI + NSHR`, where `NSHR ∈ {0,1,2,3}`.

* **Sets and systems:**  
  `NSET` = number of slip **sets** (≤ 3).  
  `NSLIP(i)` = number of slip systems in set `i`.  
  `NSLPTL` = total slip systems across all sets.

* **CFIX / CFIXA / CFIXB comments:** author modifications to track **cumulative** slip per system and overall.

---

## 1) UMAT interface: global input/output

| Name               | Dim.             | Meaning                                                        | Notes                              |
| ------------------ | ---------------- | -------------------------------------------------------------- | ---------------------------------- |
| `STRESS`           | `(NTENS)`        | Cauchy stress vector (Voigt)                                   | In/out, updated in place           |
| `STATEV`           | `(NSTATV)`       | Solution-dependent state variables                             | Layout detailed in §7              |
| `DDSDDE`           | `(NTENS, NTENS)` | Algorithmic Jacobian (∂Δσ/∂Δε)                                 | In global axes; may be unsymmetric |
| `SSE, SPD, SCD`    | scalars          | Specific strain energy, plastic dissipation, creep dissipation | Not actively set in this code      |
| `RPL`              | scalar           | Volumetric heat generation                                     | Not used here                      |
| `DDSDDT`           | `(NTENS)`        | ∂Δσ/∂T                                                         | Not used here                      |
| `DRPLDE`           | `(NTENS)`        | ∂RPL/∂ε                                                        | Not used here                      |
| `DRPLDT`           | scalar           | ∂RPL/∂T                                                        | Not used here                      |
| `STRAN`            | `(NTENS)`        | Total nominal/logarithmic strain (see code comments)           | Input                              |
| `DSTRAN`           | `(NTENS)`        | Strain increment                                               | Input                              |
| `TIME`             | `(2)`            | Step time; `TIME(1)` total, `TIME(2)` current step             | Input                              |
| `DTIME`            | scalar           | Time increment Δt                                              | Input                              |
| `TEMP, DTEMP`      | scalars          | Temperature and increment                                      | Not used here                      |
| `PREDEF, DPRED`    | `(1)`            | Predefined fields and increments                               | Not used here                      |
| `CMNAME`           | char\*8          | Material name                                                  | Input                              |
| `NDI, NSHR, NTENS` | ints             | Normal, shear, and total components                            | Input                              |
| `NSTATV`           | int              | Number of state variables                                      | Input                              |
| `PROPS`            | `(NPROPS)`       | Material constants                                             | Layout in §8                       |
| `NPROPS`           | int              | Number of material constants                                   | Input                              |
| `COORDS`           | `(3)`            | Spatial coordinates                                            | Not used here                      |
| `DROT`             | `(3,3)`          | Incremental rotation tensor                                    | Used if `NLGEOM=1`                 |
| `PNEWDT`           | scalar           | Suggested new time increment                                   | Not set here                       |
| `CELENT`           | scalar           | Characteristic element length                                  | Not used here                      |
| `DFGRD0, DFGRD1`   | `(3,3)`          | Deformation gradient (start, end)                              | Not used directly                  |
| `NOEL, NPT`        | ints             | Element and integration point ids                              | Input                              |
| `LAYER, KSPT`      | ints             | Layer and section point                                        | Not used                           |
| `KSTEP, KINC`      | ints             | Step and increment counters                                    | Not used                           |

---

## 2) Global parameters, flags, and common scalars

| Name                | Dim.      | Meaning                                            | Notes                                                      |
| ------------------- | --------- | -------------------------------------------------- | ---------------------------------------------------------- |
| `ND`                | PARAMETER | Upper bound for slip systems (here 150)            | Sizes many work arrays                                     |
| `NSET`              | int       | Number of slip **sets**                            | From `PROPS(25)`                                           |
| `NSLIP`             | `(3)`     | Systems per set                                    | Filled in initialisation                                   |
| `NSLPTL`            | int       | Total slip systems                                 | Computed                                                   |
| `THETA`             | scalar    | Implicit integration parameter ∈ [0,1]             | `PROPS(145)`; 0 explicit, 0.5 typical, 1 fully implicit    |
| `NLGEOM`            | int       | Geometry flag                                      | 0 small strain; else finite rotation/strain (`PROPS(146)`) |
| `ITRATN`            | int       | Iteration flag                                     | 0 no NR iteration; else use iteration (`PROPS(153)`)       |
| `ITRMAX`            | int       | Max NR iterations                                  | `NINT(PROPS(154))`                                         |
| `GAMERR`            | scalar    | Slip-increment absolute tolerance                  | `PROPS(155)`                                               |
| `NITRTN`            | int       | Current iteration counter                          | Starts at −1, then 0,1,…                                   |
| `CHECK`             | scalar    | Switch between hardening models                    | `0` → hyper-secant; otherwise Bassani (see code)           |
| `DEV`               | scalar    | Dilatational strain increment `∑ DSTRAN(i), i≤NDI` | Used in Jaumann term                                       |
| `GSHEAR, E11, E12`  | scalars   | Isotropic elastic moduli intermediates             | Built from `PROPS(1:2)`                                    |
| `J1,J2,J3,I1,I2,I3` | ints      | Sorting helpers for Miller indices                 | In `SLIPSYS`                                               |
| `ID, IDNOR, IDDIR`  | ints      | Running indices into arrays                        | Context specific                                           |

---

## 3) Elasticity and orientation

| Name         | Dim.    | Meaning                                           | Notes                                     |
| ------------ | ------- | ------------------------------------------------- | ----------------------------------------- |
| `DLOCAL`     | `(6,6)` | Elastic stiffness in **crystal** axes             | Built from `PROPS(1:21)` by symmetry case |
| `D`          | `(6,6)` | Elastic stiffness in **global** axes              | `D = ROTD * DLOCAL * ROTDᵀ`               |
| `ROTATE`     | `(3,3)` | Direction cosines of crystal axes in global frame | From `ROTATION`                           |
| `ROTD`       | `(6,6)` | 6×6 transformation from crystal→global            | Built from `ROTATE`                       |
| `TERM, TRM0` | `(3,3)` | Temporary matrices for `DSPIN` solve              | Used when `NLGEOM≠0`                      |
| `ITRM`       | `(3)`   | Pivot indices (LU) for 3×3 solves                 | For `TERM/TRM0`                           |

---

## 4) Slip system geometry and kinematics

| Name     | Dim.              | Meaning                                       | Notes                                        |
| -------- | ----------------- | --------------------------------------------- | -------------------------------------------- |
| `SLPDIR` | `(3, ND)`         | Slip direction unit vectors                   | In **global** axes                           |
| `SLPNOR` | `(3, ND)`         | Slip plane normal unit vectors                | In **global** axes                           |
| `SLPDEF` | `(6, ND)`         | Schmid tensors in Voigt                       | Built from `SLPDIR, SLPNOR`                  |
| `SLPSPN` | `(3, ND)`         | Slip spin tensors (skew parts)                | Used only if `NLGEOM≠0`                      |
| `DSPDIR` | `(3, ND)`         | Increment of slip directions                  | Updated if `NLGEOM≠0`                        |
| `DSPNOR` | `(3, ND)`         | Increment of slip plane normals               | Updated if `NLGEOM≠0`                        |
| `DDEMSD` | `(6, ND)`         | `D : SLPDEF` plus spin–stress terms           | See code; spin terms only if finite rotation |
| `DGAMMA` | `(ND)`            | Slip increments Δγ                            | Solved via LU                                |
| `DTAUSP` | `(ND)`            | Increment of resolved shear stress per system | From `DDEMSD·DELATS`                         |
| `DGSLIP` | `(ND)`            | Increment of current strength per system      | `H * Δγ`                                     |
| `FSLIP`  | `(ND)`            | Slip rates `γ̇` at start of step              | From `F(x)`                                  |
| `DFDXSP` | `(ND)`            | dF/dx at start of step                        | `x = τ/G`                                    |
| `H`      | `(ND, ND)`        | Self/latent hardening moduli                  | From `LATENTHARDEN`                          |
| `DDGDDE` | `(ND,6 or NTENS)` | ∂Δγ / ∂Δε                                     | Used to build Jacobian                       |

---

## 5) Stress, strain, rotation, and Jacobian helpers

| Name     | Dim.      | Meaning                                          | Notes                              |
| -------- | --------- | ------------------------------------------------ | ---------------------------------- |
| `DELATS` | `(6)`     | Lattice-stretching strain increment              | Δεᵉ = Δε − Σ SLPDEF·Δγ             |
| `DVGRAD` | `(3,3)`   | Incremental velocity gradient for stretch + spin | Only if `NLGEOM≠0`                 |
| `DSPIN`  | `(3)`     | Material spin increment components (12,31,23)    | From `DROT` via 3×3 solves         |
| `DSTRES` | `(NTENS)` | Corotational stress increment                    | Includes Jaumann and plastic parts |
| `DSOLD`  | `(NTENS)` | Stored Δσ from previous NR iterate               | For iteration                      |

---

## 6) Newton–Raphson iteration work arrays

| Name     | Dim.             | Meaning                                       | Notes                     |
| -------- | ---------------- | --------------------------------------------- | ------------------------- |
| `WORKST` | `(ND, ND)`       | NR system matrix                              | LU-factorised by `LUDCMP` |
| `INDX`   | `(ND)`           | Pivot indices (LU)                            | For `WORKST`              |
| `FSLIP1` | `(ND)`           | Saved slip rates from no-iteration pass       | For residual              |
| `STRES1` | `(NTENS)`        | Saved stress (no-iteration)                   | Fallback if NR fails      |
| `GAMMA1` | `(ND)`           | Saved slip γ (no-iteration)                   |                           |
| `TAUSP1` | `(ND)`           | Saved resolved shear stress (no-iteration)    |                           |
| `GSLP1`  | `(ND)`           | Saved current strengths (no-iteration)        |                           |
| `SPNOR1` | `(3, ND)`        | Saved normals (no-iteration)                  |                           |
| `SPDIR1` | `(3, ND)`        | Saved directions (no-iteration)               |                           |
| `DDSDE1` | `(NTENS, NTENS)` | Saved Jacobian (no-iteration)                 |                           |
| `DGAMOD` | `(ND)`           | Stored Δγ from previous NR iterate            |                           |
| `DTAUOD` | `(ND)`           | Stored Δτ from previous NR iterate            |                           |
| `DGSPOD` | `(ND)`           | Stored ΔG (strength) from previous NR iterate |                           |
| `DSPNRO` | `(3, ND)`        | Stored Δn from previous NR iterate            |                           |
| `DSPDRO` | `(3, ND)`        | Stored Δs from previous NR iterate            |                           |
| `DHDGDG` | `(ND, ND)`       | ∑ ∂H/∂γ · Δγ accumulator                      | For NR matrix augmentation |

---

## 7) `STATEV` layout (by ranges)

Let `NS = NSLPTL` be total slip systems.

| Range (1-based)     | Meaning                                              |
| ------------------- | ---------------------------------------------------- |
| `1 … NS`            | Current strength `G` per slip system                 |
| `NS+1 … 2NS`        | Shear strain `γ` per slip system                     |
| `2NS+1 … 3NS`       | Resolved shear stress `τ` per slip system            |
| `3NS+1 … 6NS`       | Components of slip plane normals (current)           |
| `6NS+1 … 9NS`       | Components of slip directions (current)              |
| `9NS+1 … 10NS`      | **Cumulative absolute shear** per slip system (`∑Δγ`) |
| `10NS+1`            | **Total cumulative absolute shear** over all systems |
| `10NS+2 … NSTATV-4` | Free room for extra model parameters (optional)      |
| `NSTATV-3`          | `NSLIP(1)`                                           |
| `NSTATV-2`          | `NSLIP(2)`                                           |
| `NSTATV-1`          | `NSLIP(3)`                                           |
| `NSTATV`            | `NSLPTL` (total systems)                             |


**Constraint:**  
`NSTATV ≥ 10*NSLPTL + 5` (or `NPARMT + 10*NSLPTL + 5` if extra parameters used).

---

## 8) `PROPS` layout (material constants)

Only non-zero segments matter for a given model.

### 8.1 Elasticity (`PROPS(1:21)`)

| Symmetry         | Required entries                                  |
| ---------------- | ------------------------------------------------- |
| Isotropic        | `E=PROPS(1)`, `ν=PROPS(2)`                        |
| Cubic            | `c11=PROPS(1)`, `c12=PROPS(2)`, `c44=PROPS(3)`    |
| Orthotropic      | `PROPS(1:9)` map to ABAQUS orthotropic constants  |
| Full anisotropic | `PROPS(1:21)` map to ABAQUS anisotropic constants |

### 8.2 Slip set descriptor and example system

| Indexes                   | Meaning                                                           |
| ------------------------- | ----------------------------------------------------------------- |
| `PROPS(25)`               | Number of slip sets `NSET` (real but integer valued)              |
| For each set `i=1..3`     |                                                                   |
| `PROPS(33..35 + 8*(i-1))` | Typical plane normal `(p q r)` in set `i` (local crystal indices) |
| `PROPS(36..38 + 8*(i-1))` | Typical slip direction `[l m n]` in set `i` (local indices)       |

### 8.3 Orientation (`PROPS(57:72)`)

Two non-parallel vectors specified in **local** and **global** frames:

| Indexes  | Meaning         |
| -------- | --------------- |
| `57..59` | Local vector 1  |
| `60..62` | Global vector 1 |
| `65..67` | Local vector 2  |
| `68..70` | Global vector 2 |

### 8.4 Rate sensitivity (`PROPS(73:96)`)

Per set (8 values each; only first two used by `F/DFDX` here):

| Entry per set | Meaning                         |
| ------------- | ------------------------------- |
| `PROP(1)`     | Rate exponent `n`               |
| `PROP(2)`     | Coefficient `γ̇₀`               |
| others        | Reserved in this implementation |

### 8.5 Hardening parameters (`PROPS(97:144)`)

Per set (16 values each), used by both hyper-secant **or** Bassani laws:

| Common entries | Meaning                                   |
| -------------- | ----------------------------------------- |
| `PROP(1)`      | `H0` initial hardening modulus            |
| `PROP(2)`      | `TAUs` (hyper-secant) or `TAUI` (Bassani) |
| `PROP(3)`      | `TAU0` initial CRSS                       |
| `PROP(9)`      | `Q` latent/self ratio (within set)        |
| `PROP(10)`     | `Q1` cross-set latent/self ratio          |

**Bassani-specific additional entries (per set):**

| Entry     | Meaning                                    |
| --------- | ------------------------------------------ |
| `PROP(4)` | `Hs` easy-glide hardening modulus          |
| `PROP(5)` | `Γ₀` peak interaction slip within set      |
| `PROP(6)` | `Γ₀` peak interaction slip across sets     |
| `PROP(7)` | `F_ab` magnitude of within-set interaction |
| `PROP(8)` | `F_ab` magnitude of cross-set interaction  |

### 8.6 Time integration and geometry (`PROPS(145:152)`)

| Index | Meaning                                                       |
| ----- | ------------------------------------------------------------- |
| `145` | `THETA` implicitness (0–1)                                    |
| `146` | `NLGEOM` switch (0 small strain, else finite rotation/strain) |

### 8.7 Iteration controls (`PROPS(153:160)`)

| Index | Meaning                                      |
| ----- | -------------------------------------------- |
| `153` | `ITRATN` switch (0 off, else on)             |
| `154` | `ITRMAX` maximum NR iterations               |
| `155` | `GAMERR` absolute tolerance on slip residual |

---

## 9) Subroutines: interface variables

### 9.1 `ROTATION(PROP, ROTATE)`

| Name             | Dim.    | Meaning                                                        |
| ---------------- | ------- | -------------------------------------------------------------- |
| `PROP`           | `(16)`  | Orientation descriptor (two vectors in local & global)         |
| `ROTATE`         | `(3,3)` | Direction cosines crystal→global                               |
| `TERM1, TERM2`   | `(3,3)` | Temporary matrices (local and its inverse; then global)        |
| `INDX`           | `(3)`   | Pivot indices for 3×3 LU                                       |
| `ANGLE1, ANGLE2` | scalars | Angles between vectors (local vs global) for consistency check |

### 9.2 `CROSS(A, B, C, ANGLE)`

| Name    | Dim.    | Meaning                                |
| ------- | ------- | -------------------------------------- |
| `A, B`  | `(3)`   | Input vectors                          |
| `C`     | `(3,3)` | Columns: unit `A`, unit `B`, and `A×B` |
| `ANGLE` | scalar  | Angle between `A` and `B` (radians)    |

### 9.3 `SLIPSYS(ISPDIR, ISPNOR, NSLIP, SLPDIR, SLPNOR, ROTATE)`

| Name             | Dim.     | Meaning                                                        |
| ---------------- | -------- | -------------------------------------------------------------- |
| `ISPDIR, ISPNOR` | `(3)`    | Integer Miller indices of exemplar direction and plane (local) |
| `NSLIP`          | int      | Number of independent systems in set                           |
| `SLPDIR, SLPNOR` | `(3,50)` | Slip directions and plane normals (global)                     |
| `ROTATE`         | `(3,3)`  | Direction cosines                                              |
| `IWKDIR, IWKNOR` | `(3,24)` | Enumerated candidate directions/planes (integer, local)        |
| `TERM`           | `(3)`    | Temporary for transforms                                       |

### 9.4 `LINE(I1,I2,I3,IARRAY)` and `LINE1(J1,J2,IARRAY,ID)`

Helpers to enumerate `{lmn}` or `{0mn}` integer triplets.

### 9.5 `GSLPINIT(GSLIP0, NSLIP, NSLPTL, NSET, PROP)` and `GSLP0(...)`

| Name     | Dim.         | Meaning                           |
| -------- | ------------ | --------------------------------- |
| `GSLIP0` | `(NSLPTL)`   | Initial strengths per system      |
| `PROP`   | `(16, NSET)` | Per-set initialisation parameters |
| `GSLP0`  | function     | Returns `TAU0` here (`PROP(3)`)   |

### 9.6 `STRAINRATE(GAMMA, TAUSLP, GSLIP, NSLIP, FSLIP, DFDXSP, PROP)`

| Name     | Dim.      | Meaning                                 |
| -------- | --------- | --------------------------------------- |
| `GAMMA`  | `(NSLIP)` | Slip at start of step (not used by `F`) |
| `TAUSLP` | `(NSLIP)` | Resolved shear stress per system        |
| `GSLIP`  | `(NSLIP)` | Current strength per system             |
| `FSLIP`  | `(NSLIP)` | Slip rates `γ̇ = F(τ/G)`                |
| `DFDXSP` | `(NSLIP)` | Derivative dF/d(τ/G)                    |
| `PROP`   | `(8)`     | Rate parameters (uses 1–2 here)         |

**Functions:**
```fortran
F(X, PROP) = PROP(2) * |X|^{PROP(1)} * sign(X)
DFDX(X, PROP) = PROP(1) * PROP(2) * |X|^{PROP(1)−1}
```

### 9.7 `LATENTHARDEN(GAMMA, TAUSLP, GSLIP, GMSLTL, GAMTOL, NSLIP, NSLPTL, NSET, H, PROP, ND)`

| Name     | Dim.           | Meaning                                                                 |
| -------- | -------------- | ----------------------------------------------------------------------- |
| `GAMMA`  | `(NSLPTL)`     | Slip at start                                                           |
| `TAUSLP` | `(NSLPTL)`     | Resolved shear stresses                                                 |
| `GSLIP`  | `(NSLPTL)`     | Current strengths                                                       |
| `GMSLTL` | `(NSLPTL)`     | **Per-system cumulative γ** (**CFIX**)                                  |
| `GAMTOL` | scalar         | **Total cumulative γ** (**CFIX**)                                       |
| `NSLIP`  | `(NSET)`       | Systems per set                                                         |
| `H`      | `(ND, NSLPTL)` | Filled with self/latent moduli                                          |
| `PROP`   | `(16, NSET)`   | Per-set hardening parameters                                            |
| `CHECK`  | scalar         | Switch: 0 hyper-secant, else Bassani (decided inside from `PROP(4..8)`) |

Self and latent are computed via `HSELF` and `HLATNT`.

### 9.8 `HSELF(...)` and `HLATNT(...)`

* **Hyper-secant**: functions of `GAMTOL` with `Q, Q1` multipliers for latent.  
* **Bassani**: functions of per-system cumulative `GMSLTL`, interaction magnitudes `FAB`, and `tanh(·/Γ₀)` aggregation.

### 9.9 `ITERATION(GAMMA, TAUSLP, GSLIP, GMSLTL, GAMTOL, NSLPTL, NSET, NSLIP, ND, PROP, DGAMOD, DHDGDG)`

Builds the matrix:
`DHDGDG = Σ ∂H/∂γ_k · |Δγ_k(old)|`
via `DHSELF/DHLATN`.

### 9.10 `DHSELF(...)` and `DHLATN(...)`

Derivatives of the self/latent moduli with respect to each `γ_k`, for both hardening laws.

### 9.11 `LUDCMP(A, N, NP, INDX, D)` and `LUBKSB(A, N, NP, INDX, B)`

Standard LU factorisation with partial pivoting and back substitution.

| Name    | Dim.       | Meaning                              |
| ------- | ---------- | ------------------------------------ |
| `A`     | `(NP, NP)` | Matrix to factor / triangular system |
| `N, NP` | ints       | Problem size and leading dim         |
| `INDX`  | `(N)`      | Pivot indices                        |
| `D`     | scalar     | Parity of number of row interchanges |
| `B`     | `(N)`      | Right-hand side / solution vector    |

---

## 10) Key algorithmic relations (for quick reference)

* **Schmid tensor per system `k`:**  
  `SLPDEF(1..6,k)` assembled from `s_k ⊗ m_k` in Voigt.

* **Spin terms (finite rotation only):**  
  `SLPSPN(:,k)` contribute to `DDEMSD` with `± SLPSPN * STRESS`.

* **Slip system linear solve:**  
  Build:
```fortran
WORKST = I + θΔt [ (∂F/∂x)/G · (D:SLPDEF)
                 + (∂F/∂x)·(x/G) · H · sign(γ̇)
                 + (∂F/∂x)·(x/G) · DHDGDG ]
```
and solve for `DGAMMA`.

* **Elastic-plastic split:**  
  `DELATS = Δε − Σ_k SLPDEF_k Δγ_k`.

* **Stress update (Jaumann form when `NLGEOM≠0`):**  
  `Δσ = D : Δε − Σ_k (D:SLPDEF_k) Δγ_k − σ·DEV` (with proper Voigt handling).

* **Jacobian:**  
  Elastic part from `D`, plastic part subtracts `DDEMSD · DDGDDE`.

---

## 11) Minimal checklist to wire a new case

* Set `PROPS(25)` and per-set plane/direction exemplars.  
* Provide elasticity block (isotropic/cubic/etc.).  
* Provide rate sensitivity per set (`n, γ̇₀`).  
* Provide hardening block per set (hyper-secant vs Bassani).  
* Choose `THETA`, `NLGEOM`, and iteration controls.  
* Ensure `NSTATV ≥ 10*NSLPTL + 5` and DEPVAR matches §7 layout.

---
