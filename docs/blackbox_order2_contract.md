# Black-box contract at order ≥ 2 — coefficients, not derivatives

> ## ⚠ THE RULE
>
> **`response.residual_coefficients_by_order` (and the `R_order_<p>` arrays) must
> contain OTI *Taylor coefficients*, NOT partial derivatives.**
>
> For a direction with exponents `κ = (κ_1 … κ_m)`:
>
> ```
>     coefficient = derivative / Π_i (κ_i !)
> ```
>
> For direction `[2, 0]` (i.e. `d²/dk²`) your response must carry
>
> ```
>     (1/2!) · d²R/dk²          ✅
> ```
>
> and **not**
>
> ```
>     d²R/dk²                   ❌  (2× too large)
> ```
>
> **The framework applies the recovery factor later**, when it writes
> `*_derivatives` and the public report. If you pre-multiply, it is applied twice.

**At order 1 every recovery factor is 1**, so a provider that returns derivatives
is indistinguishable from a correct one. The mistake only appears at order ≥ 2 —
and it appears as a *plausible-looking* number, not a crash. That is why this page
exists.

| order | direction | `κ` | recovery factor | coefficient vs derivative |
|---|---|---|---|---|
| 1 | `d/dk` | `(1,0)` | `1! = 1` | identical — mistake invisible |
| 2 | `d2/dk_df` | `(1,1)` | `1!·1! = 1` | identical |
| 2 | `d2/dk2` | `(2,0)` | `2! = 2` | **differ by 2×** |
| 3 | `d3/dk3` | `(3,0)` | `3! = 6` | **differ by 6×** |

---

## Why it is defined this way

The framework's own Python path extracts with `ctx.coeff` (the raw imaginary
coefficient), not `ctx.deriv`. A black-box provider stands in for exactly that
extraction, so it must speak the same convention. The recovery factor is applied
once, centrally, in `resasm_user/output.py` — for both paths — so that
`U_derivatives` means the same thing no matter where `R^(p)` came from.

The field is *named* `residual_coefficients_by_order` for this reason.

---

## What the framework asks you for, per order

The framework calls your command **once per order** `p = 1, 2, …, order`.

`request.json` (the fields that matter here):

| field | meaning |
|---|---|
| `u` | the converged real solution, length `ndof` |
| `parameters` | `{name: real value}` |
| `seed_directions` | `{name: 1-based basis index}` — **defines the basis order** |
| `order` | the order `p` being asked for **right now** (not the total) |
| `direction_map` | `{"<p>": [[κ…], …]}` — one exponent vector per **column** you must fill, in column order |
| `u_star_coefficients` | `{"1": U^(1), …, "<p-1>": U^(p-1)}` — the lower-order sensitivities the framework has **already solved** |

You must return, in the response:

| key | meaning |
|---|---|
| `R_order_<p>` (npz) or `residual_coefficients_by_order["<p>"]` (json) | `(ndof, N^(p))` — the order-`p` **Taylor coefficients** |
| `tangent` | `(ndof, ndof)` — `∂R/∂u` at the converged real solution (at least once; or supply it via `tangent.file`) |
| `diagnostics` | optional |

### Building `u*` (this is the part people miss)

`u*` is **not** just the real solution. It carries the imaginary coefficients the
framework has already solved:

```
    u*  =  u  +  Σ_{q=1}^{p-1} Σ_κ  U^(q)[:, κ] · e^κ      (order-p part = ZERO)
```

The order-`p` coefficients of `u*` are **zero** — they are precisely what the
framework is about to solve for. Feed `u*` and the seeded parameters
`a_i* = a_i + 1·e_i` through your residual, and read off the order-`p` coefficient.

That the framework then solves `T U^(p) = −R^(p)` is *exactly* the statement that
the order-`p` coefficient of `R(u*(a*), a*)` must vanish.

**Column ordering.** Order-1 columns are in basis order (column `i` = parameter `i`
of `seed_directions`). Order-`p` columns follow `direction_map["<p>"]`, which is
`itertools.combinations_with_replacement` over the basis indices — the same
enumeration the framework uses.

---

## The safe way to comply

**Do not hand-derive the coefficients.** Write your residual once in ordinary
arithmetic and evaluate it in a **truncated Taylor algebra**, then read the
coefficients straight off. That is what the framework does internally (with
OTILib), and what the reference template does with ~60 lines of local code:

```
resasm init --template blackbox-order2 --out my_case
```

→ [`templates/user_blackbox_order2_residual/`](../templates/user_blackbox_order2_residual/)

Its `class T2` is a minimal order-2 multivariate Taylor number. Multiplication is
a plain coefficient convolution:

```
    (A·B)_κ  =  Σ_{α+β=κ}  A_α · B_β
```

which yields coefficients directly. **No factorials appear anywhere in the
provider** — that is the tell-tale of a correct implementation.

> If you find yourself typing `math.factorial(...)` while building the response,
> stop: you are almost certainly about to convert coefficients into derivatives.

---

## Worked reference (the shipped template)

The plain spring `R = k·u³ − f` is **useless** as an order-2 reference: all of its
second parameter-derivatives are zero, so `R^(2)` is a zero vector and the bug
cannot show. The template therefore uses a residual that is nonlinear in a
parameter:

```
    R(u, k, f) = k² u³ − f          ⇒     u(k, f) = (f / k²)^{1/3}
```

At **k = 2, f = 32 ⇒ u = 2**, `T = ∂R/∂u = 3k²u² = 48`:

| direction | recovery | `R^(2)` **coefficient** *(you return this)* | `U^(2)` coefficient | **recovered derivative** | analytic `∂²u/∂…` |
|---|---|---|---|---|---|
| `k²` | `2! = 2` | `−40/3` | `5/18` | **`5/9`** | `5/9` ✓ |
| `k·f` | `1` | `1/3` | `−1/144` | `−1/144` | `−1/144` ✓ |
| `f²` | `2! = 2` | `1/96` | `−1/4608` | `−1/2304` | `−1/2304` ✓ |

Order 1: `R^(1) = [32, −1]` ⇒ `∂u/∂k = −2/3`, `∂u/∂f = 1/48`.

Note the headline number: the `k²` **coefficient** is `5/18 ≈ 0.278`, while the
true `∂²u/∂k²` is `5/9 ≈ 0.556`. **A provider that returns derivatives makes the
framework report `10/9` — exactly 2× wrong**, while its order-1 output still looks
perfect.

---

## How to check yourself

1. Run the reference template and compare to the table above.
2. Grep your provider for `factorial` — there should be none.
3. Sanity check one repeated direction by hand:
   `coefficient(k²) · 2! ` must equal your `d²R/dk²`.
4. `tests/framework/test_blackbox_order2_coefficients.py` pins all of this,
   **including a deliberately-wrong provider that must fail** — so the guard is
   itself guarded.

## Where the framework puts the result

- `private/rhs_order<p>.npz` — `residual_coefficients` / `rhs_coefficients` **and**
  `residual_derivatives` / `rhs_derivatives`, plus `recovery_factors`.
- `private/solution_sensitivities_order<p>.npz` — `U_coefficients` **and**
  `U_derivatives`.
- `private/direction_map_order<p>.json` — what every column means.
- `public/` — **recovered derivatives** only.

See [output_objects.md](output_objects.md) ·
[residual_provider_contract.md](residual_provider_contract.md) ·
[oti_execution_algorithm.md](oti_execution_algorithm.md)
