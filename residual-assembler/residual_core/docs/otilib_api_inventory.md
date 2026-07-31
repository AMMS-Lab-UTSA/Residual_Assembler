# OTILib API Inventory

Concrete API of the genuine OTILib (<https://github.com/mauriaristi/otilib.git>,
branch `master`, **GPLv3**), recorded so `otilib_adapter.py` uses the *real* names
— not guesses. Source references are paths inside that repository.

## Python module names

| Module | Implementation | Notes |
|---|---|---|
| `pyoti.sparse` | dynamic-sparse | Most versatile; **the adapter targets this**. `import pyoti.sparse as oti`. |
| `pyoti.static` | static-dense | Fastest; fixed precompiled scalar families (below). |
| `pyoti` | package root | Re-exports; adapter falls back to `.sparse`. |
| `pyoti.core` / `pyoti.real` | helpers | real-only creators, direction utilities. |

Import is available after a source build + `conda develop .` (see
[otilib_integration.md](otilib_integration.md)).

## Scalar families

- **Sparse** (`pyoti.sparse`): dynamic scalar type `sotinum` — arbitrary bases `m`
  and truncation order `nt` (up to the build limit). Real value + sparse imaginary
  coefficients.
- **Static** (`pyoti.static`): precompiled families named `onummXnY` = `X` bases,
  order `Y`. Found in the repo: `X ∈ {1,2,3,10}`, `Y ∈ {1,2,3,4,10,20,30}`
  (e.g. `onumm10n1`, `onumm2n3`, `onumm1n20`). C types `onummXnY_t`.

## C header names (for the C/Fortran seams)

- Core scalar: `include/oti/sparse/scalar/base.h` (`sotinum_t`, `soti_*`).
- Dense: `include/oti/dense/scalar/base.h` (`otinum_t`, `oti_*`).
- Static families: `include/oti/static/onummXnY/scalar/base.h`
  (`onummXnY_t`, `onummXnY_*`).
- Cython decls: `include/pyoti/c_otilib/...`.

## Creation / seeding functions (Python, `pyoti.sparse`)

| Call | Meaning |
|---|---|
| `oti.number(re, order=nt)` | real-valued OTI scalar (all imaginary coeffs zero) |
| `oti.zero(order=nt)` / `oti.one(order=nt)` | 0 / 1 scalars |
| `oti.e(imdir, order=nt)` | scalar = 1 along imaginary direction `imdir`, truncation `nt` |
| `re + oti.e(i, order=nt)` | **seed**: perturb parameter `i` (1-based basis) |

C equivalents: `soti_createReal(re, order, dhl)`, `soti_set_item(1.0, idx, ord, &num, dhl)`.

## Coefficient access functions (Python)

| Call | Returns |
|---|---|
| `x.real` (property, get **and** set) | real coefficient |
| `x.get_im(humdir)` | **raw** imaginary coefficient along `humdir` |
| `x.get_deriv(humdir)` | true derivative = coefficient × recovery factor (∏ κ_i!) |
| `x.extract_im(humdir)` / `x.extract_deriv(humdir)` | same, returned as an OTI number |
| `x.get_order_im(p)` | OTI number holding only order-`p` directions |

C equivalents: `soti_get_item` / `soti_get_im` (raw), `soti_get_deriv`
(raw × `dhelp_get_deriv_factor`), `soti_set_im_r`, `soti_set_deriv_r`.

> The framework extracts **raw** coefficients with `get_im` (to assemble
> `R^(p)` and solve `T U^(p) = -R^(p)`), and reports true derivatives by applying
> the recovery factor itself. `get_deriv` is available as a cross-check
> (`OtiContext.deriv`).

## Supported basis / order combinations

- **Sparse**: arbitrary `m`, `nt` (subject to the build's `gendata` limits). `nt`
  is passed per number via `order=` on `number`/`e`.
- **Static**: only the precompiled `onummXnY` combos listed above.

## Direction map ↔ OTILib coefficient storage

Our canonical direction is an **exponent multi-index** `κ` of length `m`
(`κ_i` = differentiation order in parameter `i`). OTILib's "human-friendly"
direction (`humdir`) is the **flat, 1-based basis-index list with multiplicity**:

```
κ = (2, 0, 1)   ⟶   humdir = [1, 1, 3]      # ε₁² ε₃
κ = (1, 1, 0)   ⟶   humdir = [1, 2]         # ε₁ ε₂
κ = (0, 0, 3)   ⟶   humdir = [3, 3, 3]      # ε₃³
```

`OtiContext._exponents_to_index_list(κ)` performs this map; the same list is used
by `oti.e(humdir, order=nt)` (build the direction) and `x.get_im(humdir)`
(read the coefficient). Verified against the repo tutorials:
`oti.e([2,2])` = ε₂², `f.get_deriv([1,1])` = ∂²f/∂x₁², `f.get_im(1)` = ε₁ coeff.

Example (`examples/python/basic_highorder.py`):

```python
import pyoti.sparse as oti
x = 3.5 + oti.e(1, order=2)
y = 0.5 + oti.e(2, order=2)
f = oti.sin(x * y)
f.real                 # value
f.get_im(1)            # ε₁ coefficient (== ∂f/∂x for order 1)
f.get_deriv([1,1])     # ∂²f/∂x²
f.get_deriv([1,2])     # ∂²f/∂x∂y
```
