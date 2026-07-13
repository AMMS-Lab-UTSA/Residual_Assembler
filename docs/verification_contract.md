# Verification Contract

What Residual_Assembler checks, what each check actually proves, and — stated
plainly — **which levels are not implemented today**.

> This ladder is about the **residual-sensitivity execution contract** (the user
> layer, `resasm check` / `resasm run`). It is a different ladder from
> `residual_core/docs/verification_strategy.md`, which grades the *finite-element
> backends* inside the framework. Do not conflate them.

---

## The ladder at a glance

| level | check | proves | status |
|---|---|---|---|
| **0** | config loads | the `resasm.yml` is well-formed and complete | ✅ implemented (`resasm check`) |
| **1** | solution vector loads; `ndof` consistent | `u` is readable, and `problem.unknowns` (if given) agrees with it | ✅ implemented |
| **2** | residual evaluates; shape is `(ndof,)` | the provider runs and returns the right shape | ✅ implemented |
| **3** | real residual norm on free DOFs | the supplied `u` really is converged | ✅ implemented **only when the residual is exposed** (Paths A/C). **Impossible in black-box mode** — see below. |
| **4** | tangent shape + free/prescribed partition | `T` is `(ndof, ndof)` and the partition is usable | ✅ implemented |
| **5** | tangent finite-difference check | that `T` really is `∂R/∂u` | ❌ **NOT implemented** |
| **6** | RHS finite-difference check: `dR/dp` at fixed `u` | the generated `R^(p)` and the solve that consumes it | ✅ implemented (`validation.rhs_finite_difference_check`) |
| **7** | sensitivity solve residual: `T U^(p) + R^(p) ≈ 0` | the linear solve was performed accurately (conditioning) | ❌ **NOT implemented** |
| **8** | full solution finite-difference check (re-solve at perturbed parameters) | the end-to-end sensitivity, including `T` | ❌ **NOT implemented** — reserved config field only |

Levels 0–4 run inside `resasm check` (before you commit to a job). Level 6 runs
inside `resasm run` when enabled.

---

## Level 0 — config loads
The `resasm.yml` parses, every required key is present, and the residual/tangent
types are supported. Errors are actionable: they name the missing key and print the
YAML block to add.

## Level 1 — solution vector loads, `ndof` consistent
`solution.file` (`.npy`) or `solution.values` is read and flattened to `(ndof,)`.
`problem.unknowns` is **optional**: if omitted it is inferred from this vector; if
supplied it must match, and a mismatch is a **hard error**, never a silent
truncation.

## Level 2 — residual evaluates, shape `(ndof,)`
The provider is invoked once at the real `(u, a)` and the returned vector's length
is checked against `ndof`. A provider that raises, or returns the wrong length, is
reported with the exact shapes.

## Level 3 — real residual norm (only if the residual is exposed)
`‖R(u,a)‖` is computed on the **free** DOFs and compared to a tolerance; a large
norm means the supplied `u` is not actually converged, and the sensitivities that
follow would be meaningless.

> **This level is impossible in black-box mode.** A Path-B executable returns only
> the perturbed RHS coefficients (and optionally the tangent) — it never hands back
> the real residual vector. There is therefore nothing to take a norm of.
> The tool does **not** guess, and does **not** print `0.000e+00`. It records
> `residual_free_norm: null`, `residual_free_norm_available: false`, prints
> `n/a — black-box did not expose real residual`, and states explicitly:
> **"Equilibrium was NOT verified by this run."**

## Level 4 — tangent shape and partition
`T` is `(ndof, ndof)`, and the free/prescribed partition from `constraints:` yields
a `T_ff` block that can be factorised. Wrong shapes are reported with the expected
dimensions.

## Level 5 — tangent finite-difference check — ❌ NOT IMPLEMENTED
*What it would do:* finite-difference `∂R/∂u` column by column and compare to the
supplied `T`.
*Why it matters:* this is the **only** check that would catch a wrong tangent. Level
6 cannot (see below), and level 8 is not implemented either. **Today, an incorrect
`T` is not detected by any implemented check.**
*(The FE side of the framework does have a tangent-FD check —
`residual_core/core/verification.py::finite_difference_tangent` — but it validates
framework backends, not a user-supplied `T` in the user layer.)*

## Level 6 — RHS finite-difference check — ✅ IMPLEMENTED
Enable with:

```yaml
validation:
  rhs_finite_difference_check: true
```

**What it does.** For each parameter `a_i`, central-difference the residual with
respect to `a_i` **at the fixed solution `u`**:

```
    dR/da_i  ≈  (R(u, a_i+h) - R(u, a_i-h)) / 2h
```

then solve the **same** linear system with the **same** tangent `T` the
hypercomplex path used, and compare the result to the order-1 sensitivity `U^(1)`:

```
    u1_fd = -T_ff⁻¹ (dR/da_i)_f      vs      U^(1)[:, i]
```

**What it proves.** The generated right-hand side `R^(1)` (i.e. the OTI extraction)
and the linear solve that consumes it.

**What it does NOT prove — read this.**
- It does **not** re-solve the nonlinear problem at the perturbed parameters.
- It **cannot detect an error in the tangent**: the same (possibly wrong) `T` appears
  on both sides of the comparison and cancels out.
- It says nothing about whether `u` is converged.

> ### It is not a full solution finite-difference validation
> Do not describe it as one, and do not read a green result as end-to-end validation
> of the sensitivities. It is a **residual-derivative (RHS)** check. The name in the
> code, the config key, the JSON key, and the public report all say
> `rhs_finite_difference_check` for exactly this reason.
> *(The legacy config key `validation.finite_difference` is still accepted as a
> deprecated alias — it meant the same thing under a misleading name.)*

## Level 7 — sensitivity solve residual — ❌ NOT IMPLEMENTED
*What it would do:* after solving, evaluate `‖T U^(p) + R^(p)‖` and confirm it is
near zero.
*What it would prove:* that the linear solve itself was accurate — it would catch an
ill-conditioned or near-singular `T_ff` (a limit point, a bifurcation, a missing
constraint). It is cheap; it is simply not wired up yet.

## Level 8 — full solution finite-difference check — ❌ NOT IMPLEMENTED
*What it would do:* re-run the **user's own nonlinear solver** at `a_i ± h`,
difference the two converged solutions, and compare to `U^(1)`.
*What it would prove:* the sensitivity **end to end**, including the correctness of
`T` — the one thing no implemented level covers.
*Why it is not implemented:* it requires the user to expose a **re-solve callback**,
which the current provider contract does not include.

A **reserved** config field exists so the intent is visible and so setting it is not
silently ignored:

```yaml
validation:
  solution_finite_difference_solver: ...   # RESERVED — NOT IMPLEMENTED
```

If you set it, the run performs no such check and says so explicitly in
`validation_summary.json` (`notes`) and in `summary.md`.

---

## Where the results land

| artifact | contents |
|---|---|
| `resasm check` (stdout) | levels 0–4, one `[ok]` / `[warn]` / `[fail]` line each, stopping at the first blocker |
| `public/validation_summary.json` | `status`, `residual_free_norm` (+ `residual_free_norm_available`, and `residual_free_norm_reason` when null), `rhs_finite_difference_check`, `solution_finite_difference_check: null` (not implemented), `notes` |
| `public/summary.md` | the same, human-readable, with the RHS-check caveat spelled out inline |
| `private/validation_full.json` | the full validation record |

---

## Honest summary of coverage

**What is verified today:** the config, the solution vector, the provider's shape
contract, the tangent's shape, the convergence of `u` (Paths A/C only), and the
generated RHS + the solve that consumes it.

**What is NOT verified today:** the tangent itself (level 5), the conditioning of
the sensitivity solve (level 7), and the end-to-end sensitivity against a real
re-solve (level 8). If your `T` is wrong, every implemented check can still pass.

---

Related: [oti_execution_algorithm.md](oti_execution_algorithm.md) ·
[privacy_contract.md](privacy_contract.md) ·
[simple_config_contract.md](simple_config_contract.md)
