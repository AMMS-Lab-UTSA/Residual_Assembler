# Partner Quickstart — OTILib (production sensitivity backend)

OTILib is the **production** hypercomplex backend for arbitrary-order,
multi-parameter residual sensitivity. (The earlier Dual1 path was only a
first-order smoke test.) You compile your residual **once**, generic in the scalar
type, and run it in two modes — real and OTI sensitivity — without sharing code
with us.

## 1. Install OTILib (external, GPLv3)

OTILib (genuine `pyoti` OTI numbers) builds from source:

- **Repo**: <https://github.com/mauriaristi/otilib.git> (branch `master`)
- **License**: **GPLv3** — install it **locally**; do not vendor it into your
  proprietary tree, and note that redistributing a binary linked against it
  triggers GPLv3 obligations.
- **Build**: `conda env create -f environment.yml && conda activate pyoti`, then
  `mkdir build && cd build && cmake .. && make && make gendata`, then
  `conda develop .`. **Windows: WSL only.**

**Do not** `pip install pyoti` — that PyPI name is squatted by an unrelated
library. Confirm `import pyoti.sparse as oti` works (Python partners) or that the
C/Fortran headers under `include/oti` are available (C++/Fortran partners).

Your proprietary residual code stays entirely local. The kit only needs the OTI
residual coefficients (or the private sensitivity package) — never your source,
mesh, or material model.

## 2. Write your residual generic in the scalar type

C++:

```cpp
#include "otilib_scalar.hpp"     // ../include

template <class Scalar>
void eval_residual(const Scalar* U, const Scalar* params,
                   std::size_t ndof, Scalar* R) {
  // ordinary arithmetic -> works for double AND an OTILib scalar
  Scalar u = U[0];
  R[0] = params[0] * u * u * u - params[1];   // k u^3 - f
}
```

- **Real mode**: instantiate with `double`.
- **OTI mode**: seed each parameter with `a_i + e_i` (all parameters at once) and
  instantiate with your OTILib scalar. Read `R^(p)` coefficients via
  `oti_coeff(R[j], exponents_to_index_list(kappa))`.

Template: [../templates/cpp_otilib_residual](../templates/cpp_otilib_residual).
Fortran: [../templates/fortran_otilib_residual](../templates/fortran_otilib_residual).

## 3. Seeding and extraction (the OTI conventions)

- Seed: `a_i* = a_i + e_i` — **all** selected parameters, simultaneously; one OTI
  evaluation carries every direction.
- Configure `m` = number of parameters, `nt` = truncation order = max order `q`.
- Extract `R^(p)`: for each order-`p` exponent multi-index `kappa`, read the
  coefficient of `R*` along `kappa`. The true partial derivative is
  `(prod_i kappa_i!) * coeff`.

## 4. Order-by-order loop (why lower orders matter)

For `p > 1`, inject the solved lower-order `U^(k)` (k<p) into `u*` **before**
evaluating order `p`. The kit / our framework does this automatically when you use
a Python provider; in a black-box executable you must implement the loop yourself
(or return all `R^(p)` for a supplied `u*` that already carries the lower orders).

## 5. Black-box OTI mode (no linking)

Your executable can handle OTI internally and return coefficient arrays. The
request carries the OTI configuration:

```json
{ "order": 2,
  "solution": [ ... ],
  "parameters": { "k": 2.0, "f": 16.0 },
  "seed_directions": { "k": 1 },
  "basis_count": 1,
  "truncation_order": 2,
  "direction_map": { "1": [[1]], "2": [[2]] },   // order -> list of exponent vectors
  "time": [0,0], "dtime": 0 }
```

Response:

```json
{ "status": "ok",
  "residual_coefficients_by_order": { "1": [[...]], "2": [[...]] },
  "tangent": [[...]],           // optional
  "diagnostics": { ... } }
```

Your model/code stays entirely behind the executable. See
[residual_provider_contract.md](residual_provider_contract.md).

## 6. Validate locally

Cross-check first-order OTI sensitivities against finite differences, and (if you
have a closed form) higher orders against analytics. See
[validation_checklist.md](validation_checklist.md).
