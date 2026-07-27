# OTILib Integration

OTILib (**O**rder **T**runcated **I**maginary numbers) is the production
hypercomplex backend for arbitrary-order, multi-parameter residual sensitivity.
It is the real HYPAD algebra behind

$$\mathbf{T}\,\mathbf{U}^{(p)} = -\mathbf{R}^{(p)}$$

where $\mathbf{R}^{(p)}$ are the $p$-th order imaginary coefficients of an OTI
residual evaluation. **Dual1 was only a first-order smoke test; OTILib is the
production path.**

## Where OTILib comes from

- **Library**: OTILib / `pyoti` — Order Truncated Imaginary numbers, by
  Mauricio Aristizabal (UTSA HYPAD group). It implements the algebra used in
  *Non-Intrusive Arbitrary-order Sensitivity Analysis of Non-Linear
  Residual-Based Problems*.
- **Repository**: <https://github.com/mauriaristi/otilib.git>
- **Branch**: `master`
- **License**: **GPLv3** (`LICENSE` in the repo).
- **Languages**: C (C99) core routines; Python 3 bindings via **Cython**;
  Fortran (F95+, static-dense). Three implementations: dynamic-dense (outdated),
  static-dense (`pyoti.static`, fastest), dynamic-sparse (`pyoti.sparse`, most
  versatile — the one this adapter targets).
- **Python module**: `import pyoti.sparse as oti` after a source build.

> ⚠️ Do **not** `pip install pyoti` expecting OTILib — that PyPI name is squatted
> by an unrelated email/URL/IP-validation library. Install the genuine library
> from the GitHub repository above. The adapter accepts a module only if it
> exposes the real OTI surface (`e` + `number` + `get_im`), so the squat is never
> mistaken for it.

## License status (GPLv3)

OTILib is **GPLv3**. To avoid imposing GPL obligations on this framework's core,
OTILib is treated strictly as an **optional external dependency**:

- No OTILib source is vendored into `residual_core/` (or anywhere in this repo).
- The library is imported/linked **only** when the user explicitly selects
  `--backend otilib`.
- If a downloader is used, it must place OTILib under a clearly-marked GPL folder
  (`external/gpl/otilib/` or `third_party/gpl/otilib/`) and the docs must state it
  is GPLv3 external code.
- Redistribution of OTILib (or a binary linked against it) triggers GPLv3 terms —
  handle deliberately.

## Build / install (external dependency)

The framework does not bundle OTILib. Build it from source (repo README):

```bash
git clone https://github.com/mauriaristi/otilib.git
cd otilib
conda env create -f environment.yml
conda activate pyoti
mkdir build && cd build
cmake ..
make
make gendata          # precomputed direction-helper data
cd ..
conda develop .       # put pyoti on the conda path
```

- **Windows**: works **only under WSL** (per the repo README). Native Windows is
  not supported by the OTILib build.
- **Helper**: `scripts/setup_otilib.sh` automates the clone + conda + cmake + make
  + `make gendata` + `conda develop .` steps and prints the detected import path.
  It is **not** run automatically by the test suite.

### Expected import / library paths

- Python package: `pyoti` (use `pyoti.sparse`), importable after `conda develop .`.
- C headers: `otilib/include/oti/...` and `otilib/include/pyoti/...`.
- Built library / precomputed data: under `otilib/build/`.

### Detection method

`residual_core/algebra/otilib_adapter.py::_probe_backend` tries, in order:

1. **Environment variables** — `OTILIB_ROOT` / `PYOTI_PATH` (python package roots,
   added to `sys.path`; each is probed at its top level and at `src/python` and
   `build`). `OTILIB_INCLUDE_DIR` / `OTILIB_LIBRARY_DIR` are recorded for C linking.
2. **Installed package** — `import pyoti.sparse` (then `pyoti.static`, `pyoti`).
3. **Local external checkouts** — `external/gpl/otilib/` or
   `third_party/gpl/otilib/` at the repo root.

Verify: `python -c "from residual_core.algebra.otilib_adapter import otilib_status; print(otilib_status())"`
should report `available: True`.

If OTILib is absent, `otilib_available()` returns `False`, OTILib tests skip
cleanly, and `--backend otilib` reports the missing dependency with the install
command (there is **no** silent fall back to Dual1).

## Scalar type & order/basis configuration

- Scalar: an OTI number with `m` imaginary bases and truncation order `nt`.
- `m` = number of design parameters selected for sensitivity.
- `nt` = maximum derivative order `q` requested.
- Total coefficients: $N = \binom{m + n_t}{m}$; order-$p$ directions:
  $N^{(p)} = \binom{p + m - 1}{p}$.

The adapter (`residual_core/algebra/otilib_adapter.py`) hides all of this behind
`OtiContext(num_bases=m, order=nt)`.

## How parameters are seeded

All selected parameters are seeded **simultaneously**, each along its own basis:

$$a_i^{*} = a_i + \varepsilon_i \qquad (\text{one OTI evaluation carries every parameter})$$

```python
ctx = OtiContext(num_bases=m, order=q)
a_star = { p: ctx.seed(value_p, basis_index_p) for p in parameters }
```

## How coefficients are extracted

`OtiContext.coeff(x, exponents)` returns the raw OTI coefficient of `x` along a
direction given by an **exponent multi-index** `kappa` (length `m`). The true
partial derivative is recovered by multiplying by the direction's recovery
factor $\prod_i \kappa_i!$:

$$\frac{\partial^{|\kappa|} f}{\partial a^{\kappa}} = \Big(\textstyle\prod_i \kappa_i!\Big)\,\mathrm{coeff}(f, \kappa).$$

## Direction map ↔ OTI basis exponents

Directions are canonical exponent multi-indices, identical to
`core/sensitivity_package.py`. For `m = 3`:

| order | exponents | label | recovery factor |
|---|---|---|---|
| 1 | [1,0,0] | e1 | 1 |
| 1 | [0,1,0] | e2 | 1 |
| 1 | [0,0,1] | e3 | 1 |
| 2 | [2,0,0] | e1^2 | 2 |
| 2 | [1,1,0] | e1*e2 | 1 |
| 2 | [1,0,1] | e1*e3 | 1 |
| 2 | [0,2,0] | e2^2 | 2 |
| 2 | [0,1,1] | e2*e3 | 1 |
| 2 | [0,0,2] | e3^2 | 2 |

`OtiContext._exponents_to_index_list` maps `kappa` to the flat basis-index list
OTILib addresses a direction with, e.g. `(2,0,1) -> [1,1,3]`. This is exactly the
argument that `oti.e([1,1,3])` and `x.get_im([1,1,3])` expect. See
[otilib_api_inventory.md](otilib_api_inventory.md) for the concrete API.

## The order-by-order residual method (what the provider does)

`core/oti_rhs_provider.py::OtiLibRHSProvider` runs:

```
u* = real converged u                 # imaginary coefficients zero
a* = a_i + e_i                        # all parameters seeded at once
for p in 1..q:
    R* = R(u*, a*)                    # evaluate through the SAME formulation backend
    R^(p) = order-p coefficients of R*
    solve T U^(p) = -R^(p)
    inject U^(p) into u*              # BEFORE evaluating order p+1
return all U^(p)
```

The injection is essential: for `p > 1` the residual's order-`p` coefficients
depend on the already-solved lower-order coefficients carried in `u*`. This is a
genuine OTI residual evaluation — not hand-coded derivative extraction.

## Backend requirements on the residual

The residual must be **OTI-safe**: written in plain scalar arithmetic so it
evaluates with OTI numbers (e.g. `nonlinear_spring1`, `nonlinear_bar1`). NumPy
matrix-kernel backends (C3D8) are not OTI-safe as written and are reported as
such — they need an OTI-templated kernel to participate.

## CLI

```bash
resasm sensitivity model.json --params params.json --order 2 --backend otilib
```

Missing OTILib is reported cleanly with install guidance; `--backend dual1` is
offered only for first-order smoke tests, never as an automatic fallback.
