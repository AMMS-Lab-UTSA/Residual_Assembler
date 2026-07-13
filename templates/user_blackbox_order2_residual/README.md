# Black-box residual template — ORDER 2 (coefficients, not derivatives)

Use this when you need **order ≥ 2** sensitivities from a private solver.
For order-1-only jobs the simpler `user_blackbox_residual` template is enough.

```
resasm init --template blackbox-order2 --out my_case
cd my_case
resasm check resasm.yml
resasm run   resasm.yml
```

## The rule this template exists to enforce

> **The arrays you return in `residual_coefficients_by_order` / `R_order_<p>` are
> OTI Taylor COEFFICIENTS — not partial derivatives.**
>
> ```
> coefficient = derivative / Π_i (κ_i !)
> ```
>
> For direction `[2,0]` (`d²/dk²`) return `(1/2!) · d²R/dk²`, **not** `d²R/dk²`.
> The framework multiplies by the recovery factor itself when it exports
> `*_derivatives` and writes the public report. Pre-multiplying makes your second
> derivatives **2× too large** (6× at order 3).

At order 1 every factor is 1, so the distinction is invisible — which is exactly
why it is easy to get wrong the first time you go to order 2.

Full contract: [`docs/blackbox_order2_contract.md`](../../docs/blackbox_order2_contract.md).

## How this template gets it right

It does **not** hand-derive coefficients. It evaluates the residual once in a tiny
**truncated Taylor algebra** (`class T2` in `my_solver.py`) and reads the
coefficients straight off — the same thing the framework does internally with
OTILib. You write the model in ordinary arithmetic and the algebra carries the
derivatives.

You edit only:

```python
def residual(u, k, f):
    return [k * k * u[0] ** 3 - f]      # R(u,k,f) = k² u³ − f

def dR_du(u, k, f):
    return [[3.0 * k * k * u[0] ** 2]]  # T = ∂R/∂u
```

Use only `+ - * **`. **No `float()` casts and no numpy ufuncs** on the parameters
or on `u` — either would collapse the Taylor number back to a plain float and
silently destroy every derivative.

## The model, and why not the plain spring

`R = k·u³ − f` (the order-1 template) has **all second parameter-derivatives equal
to zero**, so its order-2 residual coefficients are a vector of zeros — useless as
a reference. This template uses

```
R(u, k, f) = k² u³ − f
```

so the order-2 residual coefficients are visibly non-zero.

## Reference numbers (k = 2, f = 32 ⇒ u = 2, T = ∂R/∂u = 48)

Exact solution `u(k,f) = (f/k²)^{1/3}`.

| what | `k²` | `k·f` | `f²` |
|---|---|---|---|
| `R^(2)` coefficient *(what your solver returns)* | `−40/3` | `1/3` | `1/96` |
| `U^(2)` coefficient | `5/18` | `−1/144` | `−1/4608` |
| recovery factor | `2! = 2` | `1` | `2! = 2` |
| **recovered derivative** | **`5/9`** | `−1/144` | `−1/2304` |
| analytic `∂²u/∂…` | `5/9` ✓ | `−1/144` ✓ | `−1/2304` ✓ |

Note `U^(2)` coefficient for `k²` is `5/18 ≈ 0.2778` while the true `∂²u/∂k²` is
`5/9 ≈ 0.5556`. **They differ by exactly 2.** Read `U_derivatives`, not
`U_coefficients` (nor the legacy `U` key).

Order 1: `∂u/∂k = −2/3`, `∂u/∂f = 1/48`.

Pinned by `tests/framework/test_blackbox_order2_coefficients.py`, which also
asserts that a solver returning *derivatives* instead of coefficients produces a
**wrong** answer — i.e. the test really does catch the mistake.

## Protocol recap

The framework calls your command once **per order** `p = 1, 2, …`:

- `request.json` gives `u`, `parameters`, `seed_directions`, `order` (= the `p`
  being asked for **now**), `direction_map` (the exponent vector for each column
  you must fill), and **`u_star_coefficients`** — the `U^(1) … U^(p-1)` the
  framework has already solved.
- You rebuild `u*` from `u` + `u_star_coefficients` (its order-`p` part is **zero**
  — that is what is being solved for), evaluate the residual, and return the
  order-`p` **coefficients** as `R_order_<p>`, plus the `tangent` (once).
- The framework solves `T U^(p) = −R^(p)` and hands `U^(p)` back to you for the
  next order.

Order-1 columns are in basis order (column `i` = parameter `i` of
`seed_directions`). Order-`p` columns follow `direction_map[str(p)]`.

Compiled (C++/Fortran) providers may write `response.json` instead of
`response.npz` — see `docs/residual_provider_contract.md`.
