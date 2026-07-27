# The OTI Execution Algorithm

This is the method, exactly as implemented in `resasm_user/oti_global.py`
(`solve_python` for Paths A/C, `solve_executable` for the black-box Path B).

Nothing here is aspirational — every step below corresponds to code.

---

## The sequence

```
 1. Load the converged real solution  u.
 2. Load the real parameter values    a.
 3. Build the OTI parameter vector    a*   (seed each a_i with imaginary direction i).
 4. Initialize u*  with real part u and ZERO imaginary coefficients.
 5. For p = 1 .. order:
      a. Evaluate  R* = R(u*, a*)                     <- the user's residual, in OTI algebra
      b. Extract the p-th order coefficients into R^(p)
      c. Form     rhs^(p) = -R^(p)
      d. Solve    T U^(p) = rhs^(p)                   <- free DOFs only
      e. Inject   U^(p) into the imaginary directions of u*
 6. Export U^(p), R^(p) and the reports.
```

> **Step 5e is the one that is easy to get wrong and easy to omit.**
> For `order > 1`, the lower-order solution sensitivities **must** be injected into
> `u*` before the next order is computed. The residual depends on the parameters
> both directly and through `u(a)`; without the injection, the order-`p` coefficient
> extracted from `R*` is not the right-hand side of the order-`p` system, and the
> higher-order sensitivities are wrong. This is why the orders are computed **in
> sequence** and cannot be parallelised over `p`.

---

## Step by step, with the code

### 1–2. Load `u` and `a`
`u` comes from `solution.file` (a `.npy`) or `solution.values`; `a` is the
`parameters:` map. The key order of `parameters` **is** the basis order: parameter
`i` (1-based) owns imaginary direction `i`. That mapping is written to
`private/parameter_map.json` and echoed in the black-box request as
`seed_directions`.

### 3. Build the OTI parameter vector `a*`
Each real parameter is promoted to an OTI number carrying its own imaginary
direction:

```
    a*_i = a_i + 1·e_i
```

```python
ctx = OtiContext(num_bases=m, order=order)
seeded = {n: ctx.seed(float(parameters[n]), i + 1) for i, n in enumerate(names)}
```

`num_bases = m` (one per parameter) and `order` = the truncation order. This is
the "define the OTI number system" + "perturb the variable of interest" step.

### 4. Initialize `u*`
The solution vector starts as a **pure real** OTI vector — the converged solution
in the real part, nothing in the imaginary directions (we do not yet know the
sensitivities; that is what we are solving for):

```python
u_star = [ctx.scalar(float(u[j])) for j in range(ndof)]
```

### 5a. Evaluate the residual in OTI algebra
The user's residual is called **once per order**, with OTI arguments:

```python
R_oti = list(residual_fn(u_star, seeded, state, time))
```

This is the entire coupling to the user's model. The residual must be written in
generic arithmetic so the OTI type survives (Path A). `state` and `time` are passed
through unchanged.

### 5b. Extract the order-`p` coefficients
The order-`p` directions are enumerated (multiset combinations of the `m` bases),
and the coefficient of each direction is read off every component of `R*`:

```python
dirs = ctx.order_directions(p)                 # N^(p) directions
R_p  = np.zeros((ndof, len(dirs)))
for col, d in enumerate(dirs):
    for j in range(ndof):
        R_p[j, col] = ctx.coeff(R_oti[j], d["exponents"])
```

Each column of `R^(p)` is one derivative direction — e.g. at order 2 with
parameters `(k, f)`, the columns are `∂²/∂k²`, `∂²/∂k∂f`, `∂²/∂f²`. The
column ↔ direction mapping is the **direction map**.

> **⚠ These are OTI *coefficients*, not derivatives.** The code calls
> `ctx.coeff(...)`, not `ctx.deriv(...)`. For a direction with exponents `κ`, the
> true partial derivative is `Π_i (κ_i !) · coefficient` — the **recovery factor**.
> At order 1 every factor is `1`, so `U^(1)` really is `∂u/∂a`. From order 2 on,
> *repeated* directions must be scaled (`∂²u/∂k² = 2! · U^(2)[:, col(k²)]`), while
> mixed directions (`∂²u/∂k∂f`, all `κ_i = 1`) need no scaling.
>
> **You do not have to apply it yourself.** Every order exports both conventions —
> `U_coefficients` / `U_derivatives` (and `residual_*` / `rhs_*`), together with
> `recovery_factors`, `direction_exponents`, and
> `private/direction_map_order<p>.json`. The **public report quotes recovered
> derivatives**. Raw coefficients are kept for traceability. See
> [output_objects.md](output_objects.md).

### 5c–5d. Form the RHS and solve
The right-hand side is the **negative** of the extracted coefficients, and the
system is solved on the free DOFs only, with the tangent `T` evaluated once at the
converged real solution:

```python
U_p = np.zeros((ndof, len(dirs)))
U_p[free, :] = np.linalg.solve(Tff, -R_p[free, :])     # T_ff U^(p)_f = -R^(p)_f
```

Prescribed DOFs keep `U^(p) = 0` — their values are imposed, not solved, so they
have no solution sensitivity. `Tff` is factorised from the same `T` at every order.

### 5e. Inject `U^(p)` back into `u*`
The freshly solved order-`p` sensitivities are written into the corresponding
imaginary directions of `u*`, so that the next evaluation of the residual sees a
solution vector that carries the correct derivative information:

```python
for col, d in enumerate(dirs):
    for j in range(ndof):
        u_star[j] = ctx.set_coeff(u_star[j], d["exponents"], U_p[j, col])
```

This closes the loop. On the next iteration (`p+1`), step 5a re-evaluates the
residual with this updated `u*`, and the coefficient extracted at order `p+1` is
now the correct right-hand side.

### 6. Export
`R^(p)` → `private/rhs_order<p>.npz` (both `residual` = `R^(p)` and `rhs` = `-R^(p)`),
`U^(p)` → `private/solution_sensitivities_order<p>.npz` (key `U`), plus the tangent,
the maps, and the public report. See [output_objects.md](output_objects.md).

---

## The black-box variant (Path B)

`solve_executable` runs the **same** mathematics, but the split of labour changes:
the user's executable evaluates and extracts, the framework solves.

Per order `p`:

1. The framework writes a `request.json` containing `u`, `parameters`,
   `seed_directions`, `basis_count`, `truncation_order`, `order`, the
   `direction_map` for order `p`, and **`u_star_coefficients`** — the previously
   solved `U^(1) … U^(p-1)`.
2. The user's executable rebuilds `u*` from `u` + `u_star_coefficients`
   (**this is how step 5e is honoured across the process boundary**), evaluates its
   own residual, extracts the order-`p` coefficients, and writes
   `R_order_<p>` (npz) or `residual_coefficients_by_order` (json), plus optionally
   the `tangent`.
3. The framework solves `T U^(p) = -R^(p)` on the free DOFs and stores `U^(p)` into
   `u_star_coefficients` for the next order.

The tangent may be supplied once as a file (`tangent.type: file`) or returned in
the response (`tangent.type: response`).

> The framework never sees the user's model — only the request it writes and the
> response it reads. That is the privacy guarantee of Path B.
>
> It also means the framework never receives the **real residual vector**, so it
> cannot compute a residual norm and cannot verify convergence. The public report
> says so explicitly rather than printing a fabricated zero.

---

## Cost

| quantity | count |
|---|---|
| nonlinear solves performed by Residual_Assembler | **0** (the user brings a converged `u`) |
| residual evaluations | one per order (`order` total) |
| tangent factorisations | **1** (`T` is formed once and reused at every order) |
| linear solves | one per order, each with `N^(p)` right-hand sides |

Finite differences, by contrast, pay for a full nonlinear re-solve per parameter
per order. The residual method pays once for `T` and then only for back-substitutions.

---

## Failure modes this algorithm is sensitive to

- **The residual drops the OTI type** (Path A): a `float()` cast, or a numpy op that
  coerces to real, silently zeroes the imaginary directions → sensitivities come back
  as zeros or garbage. Write the residual in generic arithmetic.
- **`u` is not converged**: the whole derivation assumes `R(u,a) = 0`.
- **`T` is wrong**: nothing in the implemented verification catches this — the RHS
  finite-difference check reuses the same `T` on both sides. See
  [verification_contract.md](verification_contract.md).
- **`T_ff` singular**: limit point, bifurcation, or a missing constraint.
- **Order > 1 without injection** (only relevant if you reimplement Path B yourself):
  ignoring `u_star_coefficients` in your executable gives silently wrong
  higher-order results.

---

Related: [problem_setting.md](problem_setting.md) ·
[residual_provider_contract.md](residual_provider_contract.md) ·
[verification_contract.md](verification_contract.md)
