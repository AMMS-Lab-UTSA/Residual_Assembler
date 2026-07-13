# C3D8 residual — tangent & patch-test verification

Abaqus-**independent** correctness gate for `residual_core/c3d8_residual.py`
(Steps 4 & 7). Run:

```bash
python tests/c3d8_tangent_fd_check/run_checks.py     # table + last_report.txt, exit≠0 on fail
python residual_core/c3d8_residual.py                # same checks, verbose
```

None of these need Abaqus. They are pure math checks of the element machinery
(`B`, `detJ`, integration-point ordering, `∫Bᵀσ dV`, and the tangent mapping for
the linear case). All conventions are pinned in `residual_core/CONTRACT.md`
sections 1–2 (Voigt order `(11,22,33,12,13,23)`, engineering shear, node-major
24-DOF ordering, Abaqus C3D8 node & IP ordering).

## What each check proves

### Check 1 — Divergence-theorem patch test  *(exact, analytic)*
For a **constant** stress field `σ0`, `∇·σ0 = 0`, so the divergence theorem gives
exactly

```
∫_Ω Bᵀ σ0 dV  =  ∫_∂Ω Nᵀ (σ0·n) dA
```

The RHS is computed **independently** of the `B`-matrix path: the six C3D8 faces
are integrated as 4-node bilinear quads with 2×2 Gauss, the outward area vector
`n·dA = dx/du × dx/dv` is oriented away from the element centroid, and the
traction is `t = σ0·n`. The LHS is `element_internal_force_small_strain`.

- Runs on the **real** element 1 of `Compression111.inp`
  (connectivity `37,38,44,43,1,2,8,7`; the 8 coordinate triples are hardcoded in
  `COMPRESSION111_ELEM1_XE` and were verified against the `.inp`), **and** on a
  deliberately **distorted** hex (random ±0.6 nodal perturbation) so the result
  is not an artifact of a rectangular cell.
- Two stress states: a random symmetric tensor and a uniaxial `σ33=250`.
- Asserts relative error `< 1e-8` **and** self-equilibrium `Σ_nodes F ≈ 0`.

This is a rigorous, Abaqus-free proof that `∫Bᵀσ dV` — hence the `B` matrix,
`detJ`, quadrature weights, and node ordering — is assembled correctly. Measured
relative error is `~5e-16` (machine precision) on both geometries.

### Check 2 — Uniform uniaxial sanity  *(exact, analytic)*
Unit cube, `σ = diag(0,0,S33)`. Asserts each `+z`-face node carries `+S33·A/4`,
each `−z`-face node `−S33·A/4`, and all in-plane (x,y) components are zero.
Measured error `~1e-14`.

### Check 3a — Linear tangent finite-difference  *(exact machinery check)*
With a constant isotropic linear-elastic `D(E,ν)`, define the residual
`r(U) = element_internal_force_small_strain(Xe, σ(U))` with per-IP
`σ(U) = D·(B0·U)`. The analytic tangent is `K = Σ B0ᵀ D B0 detJ0 w`
(via `element_tangent(..., mode='small')`). Central-difference each of the 24
DOFs and compare `dr/dU` to `K`.

Report has **absolute** error, **relative (Frobenius)** error, and **max
single-entry** error; asserts relative `< 1e-6`. This validates the
**DDSDDE → element-tangent** assembly mapping for the *linear* case: the `Bᵀ D B`
contraction, the Voigt/engineering-shear pairing between `D` and `B`, and the
per-IP integration. Measured relative error `~2e-16`.

### Check 3b — Finite-strain force finite-difference  *(consistency check)*
Central-difference `d/dU` of `element_internal_force_finite_strain` at a **fixed**
Cauchy-stress field and compare to the **exact analytic Jacobian-at-fixed-σ**
`force_tangent_fixed_sigma`, i.e. the honest derivative of the actual code path

```
f_a,i(U) = ∫ (∂N_a/∂x_j) σ_ji dv ,   x = X+U
⇒  dF_a,i/dU_b,k = ∫ [ (σ g_a)_i g_b,k − (σ g_b)_i g_a,k ] dv ,   g_a = ∇N_a
```

accounting for the current-config spatial gradient `∂N/∂x` **and** the current
volume `detJ` both varying with `U`. Run at `U=0` and at a genuinely deformed
baseline (`U_amp=0.5`). Asserts relative `< 1e-6`; measured `~5e-11`.

**What 3b proves:** the finite-strain internal-force routine is smooth in `U` and
its current-configuration kinematics (spatial `B`, `detJ`) are differentiated
correctly — the assembly plumbing of the Abaqus-matching path is self-consistent.

**What 3b does NOT prove:** the *physical* finite-strain element tangent. This
fixed-σ Jacobian is **not** the standard geometric/initial-stress stiffness
returned by `element_tangent` (that is `δ_ik · g_a·σ·g_b`, `kron(G,I₃)`); the two
differ by exactly the objective-rate terms. A correct finite-strain tangent also
needs the **material** derivative `dσ/dU`.

## Known open item (to be closed later against Abaqus)

The error-prone piece **not** verified here is mapping the **finite-strain UMAT
`DDSDDE`** (ngrilli/Oxford CP, a Cauchy/objective-rate tangent) into the
finite-strain element material tangent `𝔻`, together with the conventionally-split
geometric term, so that `element_tangent(mode='finite')` reproduces Abaqus'
`AMATRX` for an `nlgeom=YES` C3D8. That requires an actual Abaqus run
(`S`, `RF`, `DDSDDE`) and is deliberately deferred — Abaqus is not installed in
this environment. Checks 1–3a are exact and Abaqus-free; 3b validates the
finite-strain *force* kinematics only.

## Files
- `run_checks.py` — runs all checks, prints the table, writes `last_report.txt`,
  exits nonzero on failure.
- `last_report.txt` — machine-generated numeric report (regenerated each run).
