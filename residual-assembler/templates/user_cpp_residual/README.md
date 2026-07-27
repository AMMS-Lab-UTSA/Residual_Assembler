# C++ residual template (Path B — local compiled code)

A complete, runnable job: `residual.cpp` builds to a single binary that is both a
self-test and the black-box responder `resasm run` calls. The model is a 1-DOF
cubic spring `R = k u^3 - f` with parameters `k`, `f`.

## Edit
`residual.cpp`, section 1 — three functions and two constants:

| symbol | meaning |
|---|---|
| `eval_residual` | `R(u, params)` |
| `eval_tangent` | `dR/du`, used for the sensitivity solve |
| `eval_dR_dparam` | `dR/d(parameter i)`, the order-1 right-hand side |
| `MODEL_NDOF` / `MODEL_NPARAM` | unknowns / parameter count |

`eval_dR_dparam` is analytic, so **edit it whenever you edit `eval_residual`**.
The self-test finite-differences `eval_residual` and fails loudly if the two ever
disagree, so you will not ship a silent mismatch.

Everything below section 1 (JSON scanning, the responder) is plumbing you can
leave alone.

## Build

```
cmake -S . -B build
cmake --build build
```

or, with no cmake:

```
g++ -O2 -std=c++17 -o residual residual.cpp
```

Both put the binary next to `resasm.yml` (`residual`, or `residual.exe` on
Windows — the compiler and CMake add `.exe` themselves), which is where
`resasm.yml` expects it.

## Run

```
./residual                      # self-test: R, dR/du, dR/dk, du/dk + an FD cross-check
resasm check resasm.yml
resasm run   resasm.yml         # -> resasm_output/public/summary.md
```

`./residual` with no arguments prints the self-test. With
`--request <in.json> --response <out.npz>` it is the black-box responder: it
reads the framework's request (`u`, `parameters`, `order`, `direction_map`),
computes the tangent and `dR/da_i`, and writes the response. That is exactly what
`residual.command` in `resasm.yml` invokes — nothing else to wire up.

**Response-file detail.** The framework passes a `--response` path ending in
`.npz`. A compiled program has no business writing numpy archives, so the binary
strips the `.npz` and writes `<that>.json`; the framework's reader falls back to
that JSON automatically (`resasm_user/providers.py::_read_response`). Keep this
behaviour if you rewrite the responder.

### Expected numbers
`solution.npy` holds the converged `u = 2` (`2u^3 = 16`). With `k=2, f=16`:
`dR/du = 3k u^2 = 24`, `dR/dk = u^3 = 8`, `dR/df = -1`, hence
`du/dk = -u/(3k) = -1/3` and `du/df = 1/(3k u^2) = 1/24 = 0.0416667`.
`resasm run` reproduces these exactly.

## Scope
This template implements **order 1**. If asked for order >= 2 it exits with a
clear error rather than returning something wrong: higher orders also need the
request's `u_star_coefficients` to rebuild `u*`. For arbitrary order, write the
residual generic in the scalar type (it already is — `template <class Scalar>`)
and instantiate it with an OTILib scalar instead of `double`; see
`../../partner_kit/include/otilib_scalar.hpp` and the commented block in
`CMakeLists.txt`. OTILib is GPLv3 —
https://github.com/mauriaristi/otilib.git (do **not** `pip install pyoti`).

## Note for Windows
`residual.command` is `./residual ...`, resolved relative to the directory you
run `resasm` from. Run `resasm run` from **this** folder (as shown above) and it
works on Windows and Linux alike. If you need to invoke it from elsewhere, put an
absolute path in `resasm.yml`.
