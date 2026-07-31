# Problem Setting — the mathematics

## The nonlinear problem

The user's model is a nonlinear residual problem:

```
    R(u, a, q, t) = 0
```

| symbol | meaning |
|---|---|
| `u` | unknown solution vector, shape `(ndof,)` — displacements, temperatures, whatever the user's model solves for |
| `a` | parameters, shape `(m,)` — the quantities we want derivatives with respect to (moduli, rates, loads, …) |
| `q` | internal / history state — plastic strain, damage, hardening variables, … |
| `t` | time / increment information (`time`, `dtime`) |

Residual_Assembler makes **no assumption** about what is inside `R`. It only
requires that the user can evaluate it. `q` and `t` are passed through to the
user's evaluator untouched.

## The converged real solution

The user solves their own problem first, by whatever means they already use
(Newton, arc-length, a commercial code). They arrive at a converged solution `u`
at the real parameter values `a`:

```
    R(u, a) ≈ 0
```

This vector is an **input** to Residual_Assembler, not an output. The tool never
performs the nonlinear solve.

> "Converged" means the free-DOF residual norm is small. If the user's residual is
> exposed (Path A/C), the tool measures `‖R(u,a)‖` on the free DOFs and warns if it
> is not small. If the residual is **not** exposed (black-box, Path B), the tool
> **cannot** check this and says so — it does not assume convergence.

## The sensitivity equation

The converged solution is an implicit function of the parameters, `u = u(a)`.
Differentiating the equilibrium condition `R(u(a), a) = 0` with respect to `a`
and collecting terms order by order gives the **sensitivity equation**:

```
    T U^(p) = -R^(p)                     for p = 1, 2, …, order
```

| symbol | definition | shape |
|---|---|---|
| `T`     | `∂R/∂u` — the tangent/Jacobian evaluated at the converged real solution | `(ndof, ndof)` |
| `R^(p)` | p-th order OTI residual coefficient matrix — one column per order-`p` derivative direction | `(ndof, N^(p))` |
| `U^(p)` | p-th order solution derivative coefficient matrix — the answer | `(ndof, N^(p))` |

`N^(p)` is the number of distinct order-`p` directions over `m` parameters
(multiset combinations): at order 1 there are `m` of them (`∂/∂a_1 … ∂/∂a_m`); at
order 2 there are `m(m+1)/2` (the `∂²/∂a_i∂a_j` with `i ≤ j`); and so on. The
mapping from a column to the derivative it holds is the **direction map**.

### ⚠ Coefficients are not derivatives (orders ≥ 2)

`R^(p)` and `U^(p)` are, as their names say, **OTI coefficient** matrices — not
partial derivatives. For a direction with exponents `κ = (κ_1 … κ_m)` (where
`Σκ_i = p`), the true partial derivative is recovered by multiplying by the
**recovery factor**

```
    ∂^p u / ∂a_1^κ1 … ∂a_m^κm   =   recovery_factor · U^(p)[:, col]
    recovery_factor = Π_i (κ_i !)
```

| order | direction | exponents | recovery factor | coefficient == derivative? |
|---|---|---|---|---|
| 1 | `∂/∂k` | `(1,0)` | `1! = 1` | **yes** |
| 2 | `∂²/∂k∂f` | `(1,1)` | `1!·1! = 1` | **yes** |
| 2 | `∂²/∂k²` | `(2,0)` | `2! = 2` | **no — multiply by 2** |
| 3 | `∂³/∂k³` | `(3,0)` | `3! = 6` | **no — multiply by 6** |

So at **order 1 the coefficients *are* the sensitivities** `∂u/∂a_i` (every factor
is 1), which is why the public `parameter_ranking.csv` can legitimately label its
column `‖du/da‖`. From **order 2 onward, pure (repeated) directions must be scaled**
by `Π κ_i!` before they are read as derivatives; mixed directions with all `κ_i = 1`
need no scaling.

### How the output handles this (you do not have to)

Nothing is left for you to infer. Every order exports **both** conventions plus the
map that connects them:

| artifact | contents |
|---|---|
| `private/solution_sensitivities_order<p>.npz` | `U_coefficients`, `U_derivatives`, `recovery_factors`, `direction_exponents` (+ legacy `U` = coefficients) |
| `private/rhs_order<p>.npz` | `residual_coefficients`, `rhs_coefficients`, `residual_derivatives`, `rhs_derivatives`, `recovery_factors`, `direction_exponents` (+ legacy `residual`, `rhs` = coefficients) |
| `private/direction_map_order<p>.json` | per column: `index`, `exponents`, `label` (e.g. `d2/dk2`), `recovery_factor`, `parameter_names` |

**The public report always quotes recovered derivatives**
(`public/sensitivity_norms.csv` carries the `recovery_factor` used for each
direction). Order-1 numbers are identical under either convention.

### Why the tangent is the *same* `T` at every order

`T` is evaluated once, at the converged real solution, and reused for every order.
This is the point of the residual method: the expensive object (a factorised
tangent) is formed once, and every additional order costs only a new right-hand
side and a back-substitution — not a new nonlinear solve. It is also why the same
`T` appears on the left of every order-`p` system above.

### Boundary conditions

The system is solved only on the **free** DOFs:

```
    T_ff U^(p)_f = -R^(p)_f  ,        U^(p) = 0 on prescribed DOFs
```

Prescribed DOFs have no solution sensitivity (their value is imposed, not solved),
so their rows are removed and their entries of `U^(p)` are left at zero. The
partition comes from `constraints.free` / `constraints.prescribed` in `resasm.yml`;
if neither is given, every DOF is free.

## Where OTI enters

`R^(p)` is the only genuinely new object, and it is what hypercomplex algebra
gives us for free.

The real solve is done **first**. Then the *already converged* residual is
**re-evaluated with OTI-perturbed parameters**: each parameter `a_i` is promoted to
an OTI number carrying an imaginary direction, the residual is evaluated once in
that algebra, and the derivative coefficients are read straight off the imaginary
components — to machine accuracy, with no step size to choose, and with every order
up to the truncation order available from the single evaluation.

```
    a_i        ->   a*_i = a_i + 1·ε_i                (perturb: seed direction i)
    u          ->   u*   = u  + 0                     (real part u, zero imaginary)
    R(u*, a*)  ->   R*                                (evaluate in OTI algebra)
    R^(p)[:,c] =  coefficient of R* in direction c    (extract order-p coefficients)
```

This is the same Hypercomplex Taylor Series Expansion (HYPAD/OTI) idea used to
extract material Jacobians from a UMAT — but applied one level up, to the **global
residual**, so that what falls out is the parameter sensitivity of the **solution**
rather than of the stress.

### The coupling that makes order > 1 work

For `p = 1`, `u*` can carry a zero imaginary part: the first-order residual
coefficient `R^(1)` is the *explicit* parameter derivative `∂R/∂a` at fixed `u`.

For `p ≥ 2` this is no longer enough. The residual depends on the parameters both
directly and through `u(a)`, so the order-`p` coefficient of `R` involves the
lower-order solution sensitivities. The algorithm therefore **injects the solved
`U^(1)`, …, `U^(p-1)` back into the imaginary directions of `u*`** before evaluating
the residual for order `p`. Once that is done, the order-`p` coefficient extracted
from `R*` is exactly the right-hand side of the order-`p` system.

> **Consequence:** the orders must be computed in sequence, `p = 1, 2, …`. They
> cannot be computed independently or in parallel over `p`.

The exact loop is in [oti_execution_algorithm.md](oti_execution_algorithm.md);
it is implemented in `resasm_user/oti_global.py`.

## What this buys, and what it costs

| | finite differences | residual method (OTI) |
|---|---|---|
| nonlinear solves | one per parameter per order (re-solve) | **none** (the real solve is already done) |
| cost per extra order | another full sweep | one more RHS + back-substitution with the *same* `T` |
| step size | must be chosen; truncation vs round-off | **none** |
| accuracy | limited by the step | machine accuracy |
| user's code | unchanged | must be evaluable in generic arithmetic (Path A) or expose OTI/derivative coefficients (Paths B/C) |

## Assumptions and limits (stated, not hidden)

1. **`u` is converged.** The equation `T U^(p) = -R^(p)` follows from `R(u,a) = 0`.
   If `u` is not converged, the sensitivities are wrong, and in black-box mode the
   tool cannot detect it.
2. **`T` is the correct tangent at `u`.** The tool does not verify it. The
   implemented RHS finite-difference check reuses the same `T` on both sides, so an
   incorrect tangent passes it. See [verification_contract.md](verification_contract.md).
3. **`T_ff` is nonsingular.** A singular/ill-conditioned free-free block (a
   bifurcation, a limit point, a missing constraint) makes the linear solve
   meaningless.
4. **History dependence must be handled by the user's provider.** `q` is passed
   through; if the model is path-dependent, the state supplied must be the state at
   the converged solution.
5. **The residual must survive the OTI algebra (Path A).** Casting to `float`,
   or any operation that discards the imaginary directions, silently destroys the
   derivatives.

---

Next: [oti_execution_algorithm.md](oti_execution_algorithm.md) —
the exact sequence, step by step.
