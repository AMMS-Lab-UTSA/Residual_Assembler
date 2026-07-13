# Template — C++ OTILib residual (compile once, run real + OTI)

Write your residual as a `template<class Scalar>`; instantiate with `double`
(real solve) and with an OTILib scalar (arbitrary-order sensitivity). Nothing is
shared with us.

- `residual.hpp` — the templated residual (edit the physics).
- Build against the seam in [../../include/otilib_scalar.hpp](../../include/otilib_scalar.hpp)
  and your OTILib install (define `RESASM_HAVE_OTILIB` and point it at OTILib).

## Two modes from one source

```cpp
#include "residual.hpp"

// real solve
double R_real[N];
eval_residual<double>(U, params, N, R_real);

// OTI sensitivity: seed each parameter a_i + e_i, evaluate, read coefficients
OtiScalar Up[N], pp[M], Rp[N];
for (int i=0;i<M;++i) pp[i] = oti_seed(params[i], i+1);   // all seeded at once
eval_residual<OtiScalar>(Up, pp, N, Rp);
// R^(p)[j] column for direction kappa:
//   oti_coeff(Rp[j], exponents_to_index_list(kappa));
```

## Order-by-order

For `p > 1`, inject solved `U^(k)` (k<p) into `Up` before the order-`p`
evaluation, then extract order-`p` coefficients and solve `T U^(p) = -R^(p)`.

## Install OTILib

Genuine OTILib builds from source: <https://github.com/mauriaristi/otilib.git>
(branch `master`, **GPLv3**, Windows = WSL only). **Not** the squatted PyPI
`pyoti`. See [../../docs/otilib_partner_quickstart.md](../../docs/otilib_partner_quickstart.md).

## Build

```bash
c++ -std=c++17 -I../../include -I<otilib_include> -DRESASM_HAVE_OTILIB \
    your_driver.cpp -o your_solver -L<otilib_lib> -loti
```
