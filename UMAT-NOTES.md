# UMAT Notes

Huang’s UMAT contains several intricate sections. The version included in this repository is a modified form of Huang’s original implementation.  
These notes serve as a guide to walk through `core.for` and to assist in understanding the UMAT’s structure and logic.

## Variable Dictionary

**DELATS**

```math
\Delta\boldsymbol{\varepsilon}^{e}
= \Delta\boldsymbol{\varepsilon}
- \sum_{\alpha}\Delta\gamma^{(\alpha)}\,\boldsymbol{\mu}^{(\alpha)}
```

**DVGRAD**

```math
\Delta\mathbf{L}^{*}
= \Delta\boldsymbol{\varepsilon}^{e}
+ \boldsymbol{\Omega}\,\Delta t
- \sum_{\alpha}\Delta\gamma^{(\alpha)}\,\boldsymbol{\omega}^{(\alpha)}
```

**DDEMSD**

```math
S_{ij}
= L_{ijkl}\,\mu^{(\alpha)}_{kl}
+ \omega^{(\alpha)}_{ik}\,\sigma_{jk}
+ \omega^{(\alpha)}_{jk}\,\sigma_{ik}
```

**DDGDDE**

```math
\frac{\partial R^{(\alpha)}}{\partial \Delta\varepsilon}
= -\theta\,\Delta t\,
\frac{\partial\dot{\gamma}^{(\alpha)}}{\partial\tau}\,
\mathbf{S}^{(\alpha)}
```

## core.for
### Lines 122 - 227
In this block of code, the `ROTATION` subroutine is called. Its main purpose is to compute the rotation matrix by using direction cosines defined in the local crystal coordinate system together with a pair of vectors that represent the crystal orientation in the global frame.

The intention of Huang was to allow an input based on the orientation in the crystal coordinates and also where the basis is not necessarily of unit length in global space. This is an elegant idea but it is much more pragmatic to keep pure rotation separated from space stretching. This is why in this implementation his strategy to solve for cubic crystals was left untouched while the HCP crystal structure option was added as an external subroutine which avoids the need of modifying the `ROTATION` subroutine.

As a side note, Huang included a check for angle relationship after inverting the matrix. Why not check first? You could argue that you are reusing TERM1 for less memory usage. However, you could still use TERM2 and keep a similar implementation. Checking first would save a little computation in cases where the angles are different enough.

> **Note**: The implementation of LU decomposition and the solver that Huang used seems to originate from a book called Numerical Recipes. 

`LUDCMP` uses scaled partial row pivoting. The L and U matrices are stored in the same matrix in this UMAT.

IMPORTANT: in the `ROTATION` subroutine Huang reuses TERM1. Initially, TERM1 is the change of basis matrix from global to local coordinates. The inverse of that matrix goes from the crystal coordinates to the global ones. This is what Huang's notes referred to when he mentions that modifying the `ROTATION` and `SLIPSYS` subroutines would be needed for non-cubic crystals. That part could be changed to adjust for non-orthogonal and different dimension basis, as mentioned earlier (HCP, for example). However, this implementation does it directly with the change of basis matrix in the noncubiccrystal subroutine. TERM2 ends up being the inverse of TERM1 up to that moment (local -> global). Then, Huang REUSES TERM1! This is confusing but saves memory. TERM1 is the new linear map that generates the rotation of the vectors. This is why the `ROTATE` matrix ends up being `TERM1 * TERM2`. 

The elastic tensor rotation from local elastic matrix to a global one does not require the aspect ratio consideration for the case of non-cubic crystals. The reason is that any anisotropy or stretch due to lattice aspect ratio is already baked into `DLOCAL`.

### Lines 227 - 318

There are six important configuration parameters that are introduced in the UMAT in this section. 

1- `NSET`: this is the total number of sets of slip systems. Slip systems are grouped in sets in order to apply physical properties to each. In this particular UMAT, some examples are the hardening modulus, the stage I stress and the initial yield stress. In the case of HCP, the sets could be the basal, prismatic and first-order pyramidal ones.

2- `THETA`: this is a time integration parameter that tells you where, within the current time increment, you evaluate your shear strain rate for the slip system alpha. If it is 0, you take the start of the step (explicit forward Euler), if it is 1 you take the end of the step (implicit backward Euler) and if it is between 0 and 1 you take a weighted average between start and end (generalised trapezoidal rule). The slip strains are obtained via integration over a finite time increment. Expression (3.1.2) in Huang's technical report is a numerical quadrature for that integral. In rate-dependent crystal plasticity, the slip rate depends on the current stress and hardening, which in turn depend on the slip strain. A value between 0.5 and 1 is recommended because this makes the scheme unconditionally stable for stiff rate equations like crystal plasticity. Low theta makes the solver more explicit and less accurate. Higher values makes it more implicit and more stable. If `THETA` is bigger than zero it usually means therefore that iteration is required. 

3- `NLGEOM`: if it is zero, small deformation theory is used, otherwise, the theory of finite rotation and finite strain is used. Activating `NLGEOM` in Abaqus tells the solver to include large-deformation kinematics and update the geometry as the analysis progresses. If inactive, Abaqus basically “freezes” the geometry at the initial configuration and pretends all your deformation is infinitesimally small. The post-processor takes the computed nodal displacements and adds them to the original coordinates to draw a deformed shape. That is why the ODB will still look bent or twisted, but the numbers you see (stresses, strains) are from a linearised kinematic assumption. Internally, Abaqus uses the initial shape for element shape functions, Jacobian determinants and the B matrix. From a user's perspective, it could be argued that `NLGEOM` should now be active by default; however, this convention originates from a time when such computations were computationally expensive. In this UMAT, if `NLGEOM` is not active in Abaqus FEA while the internal formulation employs the theory of finite rotations and strains, a mismatch arises. Consequently, numerical errors may occur, and the simulation may fail to converge. 
`NLGEOM` set to zero uses the small-strain theory that leads to a linearised form of the update of the slip plane directions and normals. In Huang's technical report, Equations (3.3.1a) and (3.3.1b) end up simply as $s_{\text{new}} \approx (\mathbf{I} + \boldsymbol{\Omega})\ s_{\text{old}}$ (the same concept applies to $m$ since the exponential rotation operator is linearised under the assumption of infinitesimal rotations). In this context, finite deformation theory means that the kinematics (the way in which stresses, slip systems' orientations, and rotations are updated) come from the exact finite‐strain continuum mechanics relations, not small‐strain approximations.

4- `ITRATN`: if this is not active, that is, it is zero, the iterative Newton-Raphson solver is not used. 

5- `ITRMAX`: maximum number of iterations.

6- `GAMERR`: gamma error (used for tolerance convergence).

> **Note**: A set of mathematically inconsistent parameters could be picked. Be careful with the parameters you select.

The number of iterations parameter `NITRTN` is initialised with -1.

The increment of spin associated with the material element (`DSPIN`) is also commonly referred to as $\Delta \Omega$. For a finite rotation case, Huang used the exact Rodrigues formula. 

$$(I + \Delta R^{\mathsf{T}}) \ \Delta \Omega = \Delta R^{\mathsf{T}} - I$$

$$\Delta \Omega = (I + \Delta R^{\mathsf{T}})^{-1} \ (\Delta R^{\mathsf{T}} - I)$$

For the case of small deformation theory, the following approximations are used:

$$\Delta \Omega \approx \Delta R - I$$

$$\delta \Omega \approx \frac{\delta R - \delta R^{\mathsf{T}}}{2}$$

### Lines 318 - 383

> **Note**: From 318 to 449 the initial state is developed. This lines only run one time at the start of the simulation.

In this section the slip direction and normal unit vectors are defined and stored in `SLPDIR` and `SLPNOR`.  

The `SLIPSYS` subroutine constructs all crystallographically equivalent slip systems for a given cubic slip family, starting from one representative slip direction and one slip plane normal provided as integer Miller indices. It first classifies the input direction and plane normal into standard cubic families (e.g., ⟨100⟩, ⟨110⟩, ⟨111⟩ for directions; {100}, {110}, {111} for planes) and, using helper routines `LINE` and `LINE1`, generates all permutations and sign variations allowed by cubic symmetry while avoiding duplicates that differ only by reversal of both direction and normal. It then filters the combinations to keep only those that satisfy orthogonality between slip direction and slip plane normal. The resulting sets are normalised to unit vectors and finally rotated from the crystal frame into the global coordinate frame using the provided lattice orientation matrix. The output lists, `SLPDIR` and `SLPNOR`, contain the slip directions and plane normals for every valid slip system in that family, and the total count `NSLIP` is returned.

The `NONCUBICCRYSTAL` subroutine defines slip systems for noncubic crystals, with its current implementation tailored to the HCP structure of α-Ti. Instead of generating equivalent systems via cubic symmetry, it explicitly specifies the slip plane normals and slip directions for each system in the crystal's native, nonorthogonal coordinate frame. A change-of-basis matrix, built from the HCP lattice parameters, is applied to transform these vectors into an orthogonal Cartesian frame. The transformed directions and normals are then rotated into the global frame using the provided lattice orientation matrix. Each vector is normalised to unit length to ensure correct geometric representation. The subroutine concludes by defining the number of systems per slip family (basal, prismatic, pyramidal) and setting the total slip system count, returning fully initialised arrays SLPDIR and SLPNOR for use in subsequent UMAT calculations.

Huang computed the symmetric part of the slip deformation tensor (Schmid tensor $\mu$) because only the symmetric component of the dyadic product between slip direction and slip plane normal contributes to the plastic strain rate. The antisymmetric part represents a rigid-body rotation, which does not produce deformation, so it is discarded in the constitutive update.

$$
\mu_{ij}^{(\alpha)} = \frac{1}{2} \left( s_i^{(\alpha)}\ m_j^{(\alpha)} + s_j^{(\alpha)}\ m_i^{(\alpha)} \right)
$$

This is stored in `SLPDEF` where each column contains the components $(\mu_{11}, \mu_{22}, \mu_{33}, 2\mu_{12}, 2\mu_{13}, 2\mu_{23})$.

### Lines 383 - 449

From line 383 to line 404 the number of slip systems per set, the slip directions and slip normals are saved as state dependent variables. 

Afterwards, the `GSLPINIT` subroutine sets the initial slip resistance (current strength) for every slip system in the model. The initialisation is isotropic within each set of slip systems, meaning all systems in a set start with the same strength. The routine loops over each set and slip system, calling `GSLP0` to get the starting value. In the present implementation, `GSLP0` simply returns the value of `PROP(3)` for the corresponding set, which is the initial critical resolved shear stress (local yield stress). This value represents the starting resistance to slip for that system. Huang’s structure allows easy modification of `GSLP0` for alternative initialisation laws, but here it is kept minimal.

Between lines 408 and 418, the shear strain of all slip systems, as well as the total cumulative shear strain on each slip system, are initialised to zero.

In lines 418-433, the initial resolved shear stresses are computed as $\tau^{(\alpha)} = \mathbf{P} : \boldsymbol{\sigma}$ according to Schmid's law. It is worth noting that this section clearly shows what would occur if plane strain or plane stress conditions were applied, specifically which components contributing to the shear stress in each slip system would be omitted.

In the final block of code, the kinematic hardening variables (usually referred as $X$, which represents a backstress) and their increments are initialised to zero as well. Three groups of variables were included in the code to allow the implementation of Chaboche's model. 

### Lines 449 - 557

The current stress state retrieval starts here. 

This section starts by retrieving from the state dependent variables the total number of slip systems and the number of slip systems per set. Then, it also retrieves the slip directions and normals. After that, it calculates the symmetric Schmid matrix per slip system and stores it in a Schmid factors matrix. This is stored in `SLPDEF` where each column contains the components $(\mu_{11}, \mu_{22}, \mu_{33}, 2\mu_{12}, 2\mu_{13}, 2\mu_{23})$.

If finite strain theory is active (`NLGEOM` is input as PROPS(146)), the antisymmetric part of the Schmid tensor for each slip system is computed (`SLPSPN`).

$$
\omega_{ij}^{(\alpha)} = \frac{1}{2} \left( m_i^{(\alpha)}\ n_j^{(\alpha)} - m_j^{(\alpha)}\ n_i^{(\alpha)} \right)
$$

From line 503 to 557 the double dot product of the elasticity matrix (6x6 in Voigt notation) and the Schmid factors matrix is done. The way in which the product is computed in terms of indices seems wrong since the D indices appear to be reversed. However, elasticity matrices are symmetric to their main diagonal and therefore the product is actually correct. Then, why did Huang did this? The reason lies in Fortran's memory layout, where data is stored contiguously following the columns of the matrices. The way in which the multiplication is done gives better cache performance because each successive K is right next to the last in memory. As a result, the inner loop is vectorisable and the code runs much faster.

The rest is the effect that is caused by the rotation, which only contributes when `NLGEOM` is not zero.

### Lines 557 - 614

This part of the code can be divided in three main sections. The `STRAINRATE` subroutine call to compute the slipping rates and the function derivative, the `LATENTHARDEN` call to calculate the H matrix that contains the self- and latent- hardening moduli, and the Newton-Raphson linear system assembly for slip increments.

In the first part, the `STRAINRATE` subroutine is called per set of slip systems. Within the subroutine, computations are done to solve for the the Hutchinson’s power-law viscoplasticity model (1976). As a result, a function value and a function derivative value are calculated for each slip system. The function value is the slipping rate according to the constitutive formulation. The derivative is calculated to be used later as part of the Newton-Raphson scheme that is used to obtain the increment of shear strain per slip system.

$$
x_\alpha = \frac{\tau_\alpha - \chi_\alpha}{g_\alpha}
$$

$$
F_\alpha(x) = \dot\gamma_0 \ \text{sgn}(x)\ |x|^n
$$

$$
\frac{dF_\alpha}{dx} = \dot\gamma_0 \  n\ |x|^{n-1}
$$

In the second part, the `LATENTHARDEN` subroutine is called. The self-hardening modulus represents how slip system *i* hardens itself as it shears. The latent-hardening modulus represents how slip system *j* makes slip system *i* harder. 

Just as a reference, the equivalence used to calculate the $\mathrm{sech}(x)$ in the code is the following:

$$
\mathrm{sech}(x) = \frac{2\ e^{-x}}{1 + e^{-2x}}
$$

The third part is particularly dense and to understand it, it is necessary to consider the Taylor expansion and the Newton-Raphson (NR) formulation for a system of nonlinear equations. In many textbooks, however, only the scalar form of the NR method is presented, whereas in this context the full system formulation must be considered. Furthermore, Huang uses `TERMX` where `X` is a number to enumerate expressions that are used to build the `WORKST` matrix. The `WORKST` matrix is the Jacobian of the NR scheme. "WORK" was commonly used to name working arrays is many old Fortran codes. The abbreviation "ST" following "WORK" could have multiple interpretations; however, it is reasonable to assume that it refers to "system", which is consistent with the context of the Newton-Raphson formulation.

In the third part of these lines, the `WORKST` matrix is built and then it is decomposed using the `LUDCMP` subroutine. The precise definition of each `TERM` can be found by expanding the residual that results from linearising Equation (3.4.1) of Huang's technical report around a guess using the Taylor-expansion.

### Lines 614 - 654

In this section, the right-hand side of the Taylor expansion or the Newton-Raphson formulation is built. The Taylor expansion approximation is used when `NITRTN.EQ.0` is true. 

Afterwards, the LUBKSB subroutine is called, where the lower-upper backward substitution is performed to calculate the new correction factor $\Delta\gamma_{\text{corr}}^{(\alpha)}$.

The correction is added to the previous guess (which is stored in `DGAMOD`).

### Lines 654 - 706

This section of Huang's UMAT performs state-variable updates after each Newton-Raphson iteration of the crystal plasticity algorithm. At this stage, the algorithm has already solved the linearised system that provides the new estimate of the slip increments $\Delta\gamma^{(\alpha)}\$ for all slip systems. All the material state variables (shear strain, slip resistance, and kinematic backstress) are updated to reflect the most recent iteration of the solution.

The key to understanding this section is realising that the updates are performed while looping, not at the end of the Newton–Raphson process. Each iteration immediately overwrites the previous estimates of the internal variables with the newly corrected values. This might seem unnecessary at first, but it is a deliberate design choice aimed at saving memory and keeping the algorithm compact. Instead of storing multiple temporary copies of all the state variables for every iteration, this code aims at simply updating them in place. This means the UMAT always holds the latest estimate of the material state without allocating extra arrays or buffers, which was an important consideration given the limited memory resources available at the time the code was written.

### Lines 706 - 767

In this part of the code, the increment of strain associated with lattice stretching is computed first. Abaqus FEA passes strain increments in the Jaumann corotational frame. This is calculated as:

$$
\Delta\varepsilon^{e} = \Delta\varepsilon - \Delta\varepsilon^{p}
$$

or, in the form used in the UMAT:

$$
\texttt{DELATS} = \texttt{DSTRAN} - \sum_{\alpha} \boldsymbol{\mu}^{(\alpha)} \ \Delta\gamma^{(\alpha)}
$$

If the corotational frame was not used, `DSTRAN` would include rigid-body rotation components. Then, subtracting only the plastic part would not yield the pure elastic (stretching) contribution, because the rotation-induced strains would contaminate it. But since the Jaumann formulation already filters those out, the subtraction is physically meaningful: it directly isolates the elastic lattice stretching in the corotating axes (exactly what is needed to update stress and the lattice orientation correctly).

From line 725 to 767 the lattice velocity gradient increment is computed by using `DELATS`:

$$
\underbrace{\texttt{DELATS}}_{\texttt{DSTRAN} - \sum_\alpha \mu^{(\alpha)}\ \Delta\gamma^{(\alpha)}\ \text{(symmetric part)}}\
+\
\underbrace{\texttt{DSPIN}}_{\Omega\ \Delta t\ \text{(macroscopic spin)}}\
-\
\underbrace{\sum_\alpha \texttt{DGAMMA(α)}\ \texttt{SLPSPN(α)}}_{\omega^{(\alpha)}\ \Delta\gamma^{(\alpha)}\ \text{(slip-system spin)}}
$$

### Lines 767 - 866

In the first group of computations, the increment of resolved shear stress in a slip system is calculated and the updated resolved shear stress per slip system is stored in state variables. From line 767 to line 775, Equation (3.2.3) of Huang's technical report is computed.

From line 782 to line 827 Equation (3.2.4) is computed and the increment of corotational stress in calculated. The stress is updated from 827 to 832. 

From line 832 to 866 Equations (3.3.2a) and (3.3.2b) are computed and the updates of the normal and slipping direction are updated. 

### Lines 866 - 961

In this section the material Jacobian matrix is calculated. 

### Lines 961 - 1088

This section of the code is for the NR iterative solution. If the UMAT was configured to provide a solution without iterations, i.e., a solution using Taylor expansion, the branches is this section will be unused.

In the first part, the initial solution is saved in alternative arrays. In this way, even if convergence doesn't occur with the desired tolerance, a solution is provided.

Later, the increments are stored to keep track of the accumulated increments of the variables utilised during the loop. 

Between lines 1026 and 1041, whether the iteration solution converges is analysed. `TERM1` is simply the effective resolved shear stress driving slip in a particular slip system. The variable `RESIDU` is:

$$
\texttt{RESIDU} = R^{(\alpha)} = \Delta t \left[ \theta\ \dot{\gamma}_{\text{new}} + (1 - \theta)\ \dot{\gamma}_{\text{old}} \right] - \Delta\gamma_{\text{current}}
$$

and equivalently:

$$
R^{(\alpha)}=
\underbrace{\theta\ \Delta t\ F(X)}_{\text{new slip increment}}\
+\
\underbrace{(1-\theta)\ \Delta t\ F_{\text{old}}}_{\text{old slip increment}}\
-\
\underbrace{\Delta\gamma^{(\alpha)}}_{\text{computed increment}}
$$

The call to the `ITERATION` subroutine is to compute and update the array `DHDGDG`, a part of the Jacobian of the NR scheme. 

If the solution does not converge, the last block of code saves the solution that results from the first step of the iterative solver.

### Lines 1088 - 1134

The first part is the calculation of the aggregated cumulative shear strains on all slip systems, and the cumulative strain in individual slip systems. 

The effective plastic slip parameter is a scalar measure of equivalent plastic deformation and it is calculated in the second block of code. 
