# Residual provider contract

The residual provider is the one thing you must write. It is the object that
makes the whole run possible: something that can evaluate

```
R(u, params, state, time)
```

locally, on your machine, without handing anyone your model. Everything else the
toolkit does — seeding `a_i* = a_i + e_i`, extracting `R^(p)`, solving
`T U^(p) = -R^(p)` — hangs off this contract.

There are three ways to satisfy it.

| Path | `residual.type` | You provide | Who runs the hypercomplex algebra |
|---|---|---|---|
| **A — Python** | `python` | a Python callable | the framework (OTILib), in-process |
| **B — Executable** | `executable` | a command + request/response files | **you** (inside your program) |
| **C — Compiled C++/Fortran** | `executable` | a compiled binary speaking the same request/response contract | **you** (inside your binary) |

Paths B and C are the same contract; C is just the compiled instance of it, with
one extra file-format concession documented below.

Not sure which one? See
[`which_path_should_i_use.md`](which_path_should_i_use.md).

---

## Path A — Python provider

### Signature

```python
def residual(u, params, state=None, time=None):
    ...
    return R          # length-ndof list / sequence

def tangent(u, params, state=None, time=None):
    ...
    return T          # (ndof, ndof) nested list / array
```

Both live in the same module (`residual.module`). Names come from
`residual.function` (default `residual`) and `tangent.function` (default
`tangent`).

```yaml
residual:
  type: python
  module: user_residual.py
  function: residual
tangent:
  type: python
  function: tangent
```

The working template (`templates/user_python_residual/user_residual.py`) is the
whole thing:

```python
def residual(u, params, state=None, time=None):
    k = params["k"]
    f = params["f"]
    return [k * u[0] ** 3 - f]        # cubic spring:  R = k u^3 - f


def tangent(u, params, state=None, time=None):
    # optional. dR/du at the real solution (plain floats here).
    k = params["k"]
    return [[3.0 * k * u[0] ** 2]]
```

### CRITICAL: your residual is called with OTI numbers

Your `residual` is called **twice, with different scalar types**:

1. **Real call** (`runner.py`) — `u` is a NumPy float array, `params` is a
   `dict[str, float]`. Used for the real residual and its free-DOF norm.
2. **OTI call** (`oti_global.py::solve_python`), once per order `p = 1..q`:

```python
ctx    = OtiContext(num_bases=m, order=order)
seeded = {n: ctx.seed(float(parameters[n]), i + 1) for i, n in enumerate(names)}
u_star = [ctx.scalar(float(u[j])) for j in range(ndof)]
...
R_oti  = list(residual_fn(u_star, seeded, state, time))
```

So on the OTI call:

- **`u` is a plain Python `list` of OTI scalars** (length `ndof`) — not a NumPy
  array. Index it (`u[0]`); do not call `.shape`, `.dot`, or `.T` on it.
- **`params` is a `dict` whose values are OTI scalars** — `params["k"]` is
  `k + e1`, not `2.0`.
- The return value is consumed by `list(...)`, and each entry `R_oti[j]` must
  still carry its imaginary coefficients, which are read with
  `ctx.coeff(R_oti[j], exponents)`.

**What that means for your code:**

| Do | Don't |
|---|---|
| Use ordinary arithmetic: `+ - * / **` | `float(x)`, `int(x)`, `round(x)` — these collapse an OTI number to its real part and **silently destroy every derivative** (you get zeros) |
| Return a plain Python `list` (or an object-dtype sequence) | `np.array(R, dtype=float)` / `np.zeros(n)` then assigning OTI values into it — a float array cannot hold an OTI number |
| Keep the code generic in the scalar type: it must work when `u[j]` is a `float` **and** when it is an OTI scalar | `math.sqrt`, `math.exp`, `math.sin`, … — the `math` module coerces to `float` and drops the OTI part |
| Compare/branch on real quantities only if you must (branches are fine; they just do not carry derivatives) | `np.linalg.*`, `np.dot`, and other NumPy kernels that require a numeric dtype |
| Read `params` by name | Assume a parameter order inside `params` — it is a dict; the seed order is the `parameters:` key order in `resasm.yml` |

If a derivative comes out as exactly zero when it should not be, the usual cause
is a `float()` cast or a NumPy float array somewhere in the chain.

`tangent(...)` has no such constraint: it is only ever called with plain floats
(`runner.py::_acquire_tangent`), and the result is immediately
`np.asarray(..., float)`. It must be `(ndof, ndof)`.

### Why the residual is called once per order

`solve_python` runs the HYPAD residual loop:

```
for p = 1..q:
    R*      = residual(u*, a*)        # OTI evaluation
    R^(p)   = order-p coefficients of R*
    solve   T U^(p) = -R^(p)          # free DOFs only
    inject  U^(p) into u*             # BEFORE order p+1
```

Between orders, the framework writes the just-solved `U^(p)` back into the
imaginary part of `u*`. Your function does nothing special — it just gets called
`q` times with an increasingly enriched `u*`.

### The `time` argument

`runner.py::_time_pair` builds `((t, t), dtime)` from `time.time`. The two call
sites pass different things:

| Call | `time` your function receives |
|---|---|
| real residual, tangent, FD check | `time_pair[0]` — a **float** |
| OTI residual (inside `solve_python`) | `time_pair` — a **2-tuple** `(t, t)` |

`dtime` is never passed to a Python residual. Do not depend on the type or
arity of `time`; if you need it, handle both, or carry the value in `state`.

### `residual.type: element`

Accepted by the config, but in the current code `runner.py` builds the same
`PythonResidual` for it — module + function, exactly like `python`. There is no
mesh/element assembly in the run path.

---

## Path B — Executable (black-box) provider

Your program stays entirely yours. The framework only exchanges two files.

```yaml
residual:
  type: executable
  command: ./my_solver --request {request} --response {response}
tangent:
  type: response          # or: type: file / file: tangent.npz
```

`{request}` and `{response}` are substituted with absolute paths
(`providers.py::ExecutableResidual`). Mechanics:

- The command runs with `cwd` = the directory containing `resasm.yml`.
- Request and response live in a temporary directory created **inside** that
  directory and deleted afterwards.
- `subprocess.run(..., timeout=3600)`.
- On Windows a first token ending in `.bat` is routed through `cmd /c`.
- The framework does not inspect the exit code. What it checks is that a
  response file exists; if none does, it raises and prints the tail of your
  stdout and stderr.
- **The executable is invoked once per order** `p = 1..q`, not once for all
  orders.

### Request (framework → you), `request.json`

Fields, verbatim from `providers.py::ExecutableResidual.eval_rhs`:

| Field | Type | Meaning |
|---|---|---|
| `schema` | string | always `"resasm-user-request/1"` |
| `u` | list of `ndof` floats | the converged real solution |
| `parameters` | `{name: float}` | real parameter values, **in seed order** |
| `seed_directions` | `{name: int}` | parameter → 1-based imaginary basis index (`{"k": 1, "f": 2}`) |
| `basis_count` | int | `m` = number of parameters |
| `truncation_order` | int | **the order currently being requested** (same value as `order`) |
| `order` | int | the order `p` currently being requested |
| `direction_map` | `{"<p>": [[exponents], ...]}` | the exponent vectors for **this** order only |
| `u_star_coefficients` | `{"<k>": [[...]]}` | the already-solved `U^(k)` for `k < p` (empty at `p = 1`) |
| `state` | anything JSON-serializable, or `null` | your `state` block, verbatim |
| `time` | `[t, t]` | from `time.time` |
| `dtime` | float | from `time.dtime` |

Note both `order` and `truncation_order` carry the current loop order `p` — the
total requested order is not sent as a separate field.

```json
{ "schema": "resasm-user-request/1",
  "u": [2.0],
  "parameters": { "k": 2.0, "f": 16.0 },
  "seed_directions": { "k": 1, "f": 2 },
  "basis_count": 2,
  "truncation_order": 1,
  "order": 1,
  "direction_map": { "1": [[1,0],[0,1]] },
  "u_star_coefficients": {},
  "state": null, "time": [0.0, 0.0], "dtime": 0.0 }
```

#### `direction_map` — what a column *is*

`direction_map["<p>"]` is the ordered list of **exponent vectors** for order `p`.
Each vector has length `m` (one slot per parameter, in seed order), and its
entries sum to `p`. Column `c` of the `R^(p)` you return must correspond to
`direction_map["<p>"][c]`.

At order 1 the vectors are one-hot, so **the index of the `1` selects the
parameter**:

| direction | exponents (m = 2) | column means |
|---|---|---|
| `e1` | `[1, 0]` | `dR/dk` |
| `e2` | `[0, 1]` | `dR/df` |

Both templates do exactly that lookup:

```python
i = int(np.argmax(exps))        # my_solver.py — which parameter this direction differentiates
```

```cpp
std::size_t ip = 0;                                  // residual.cpp
for (std::size_t e = 0; e < dirs[c].size(); ++e)
  if (dirs[c][e] != 0) { ip = e; break; }
```

At order 2 with `m = 2` you get `[[2,0],[1,1],[0,2]]` — `e1^2`, `e1*e2`, `e2^2`.
The values you return are OTI **coefficients**, i.e. the Taylor coefficients:
the true partial derivative is the coefficient times the product of the
factorials of the exponents (`e1^2` → `×2`, `e1*e2` → `×1`).

#### `u_star_coefficients` — why order ≥ 2 needs it

The residual must be differentiated at `u*`, and `u*` is not just the real `u`:
after each order is solved, the solution sensitivity `U^(p)` is injected into
its imaginary part. The framework does that injection internally on the Python
path; on the black-box path it cannot reach inside your program, so it **feeds
the solved lower orders back to you**:

```json
"u_star_coefficients": { "1": [[-0.33333333, 0.04166667]] }
```

Key `"<k>"` → the `U^(k)` matrix as nested lists, shape `(ndof, N^(k))`, using
the same column order as `direction_map["<k>"]`. (Above: `ndof = 1`, two
parameters, so one row and two order-1 columns — `du/dk` and `du/df`.)

- At `order: 1` this is `{}` — you can ignore it.
- At `order ≥ 2` you **must** consume it to rebuild `u*` before differentiating,
  otherwise your `R^(2)` is wrong. The shipped templates implement order 1 only
  and say so explicitly (the C++/Fortran ones exit with an error for
  `order != 1`).

### Response (you → framework)

> ## ⚠ Black-box response arrays are **COEFFICIENTS**
>
> **Do not multiply by recovery factors yourself.** The framework does that when
> it exports derivatives and writes the public report.
>
> The arrays in `residual_coefficients_by_order` / `R_order_<p>` must be OTI
> **Taylor coefficients**:
>
> ```
>     coefficient = derivative / Π_i (κ_i !)
> ```
>
> For direction `[2,0]` (`d²/dk²`) return `(1/2!)·d²R/dk²`, **not** `d²R/dk²`.
>
> **At order 1 every factor is 1, so this mistake is invisible.** At order ≥ 2 it
> silently doubles your second derivatives (×6 at order 3). If you are writing
> `math.factorial(...)` while building the response, you are almost certainly about
> to get this wrong.
>
> Full rules + a worked, tested reference:
> **[blackbox_order2_contract.md](blackbox_order2_contract.md)** and
> `resasm init --template blackbox-order2`.

Two accepted forms. Both are read by `providers.py::_read_response`.

**`.npz`** (what `--response` literally asks for):

| Key | Shape | Required |
|---|---|---|
| `R_order_<p>` | `(ndof, N^(p))` | yes — the order `p` that was requested |
| `tangent` | `(ndof, ndof)` | only if `tangent.type: response` (i.e. you are not supplying `T` from a file) |
| `diagnostics` | a JSON **string** stored in the archive | optional |

**`.json`**:

```json
{ "residual_coefficients_by_order": { "1": [[8.0, -1.0]] },
  "tangent": [[24.0]],
  "diagnostics": {"solver": "template_cpp", "method": "analytic"} }
```

`residual_coefficients_by_order["<p>"]` is the `R^(p)` matrix as nested lists,
rows = DOFs, columns = directions.

If the requested order is absent from the response, the run fails with
`executable did not return R^(<p>) (response had orders [...])`. If no tangent
is available from anywhere, it fails with `no tangent available: provide
tangent.file or return \`tangent\` in the black-box response`.

Your `diagnostics` are parsed but not written to any output file; only the
engine name reaches `private/metadata.json`.

The working reference is `templates/user_blackbox_residual/my_solver.py`, which
computes both the tangent and `R^(1)` with its **own** finite differences — it
needs nothing from the framework, not even OTILib.

---

## Path C — Compiled C++ / Fortran provider

Same `residual.type: executable` contract, same `request.json`. You write your
model once, generic in the scalar type, and (optionally) instantiate it with an
OTI scalar type for arbitrary-order sensitivities:

```cpp
// templates/user_cpp_residual/residual.cpp
template <class Scalar>
void eval_residual(const Scalar* u, const Scalar* params,
                   std::size_t /*ndof*/, Scalar* R) {
  Scalar k = params[0];
  Scalar f = params[1];
  R[0] = k * (u[0] * u[0] * u[0]) - f;
}

template <class Scalar>
void eval_tangent(const Scalar* u, const Scalar* params,
                  std::size_t /*ndof*/, Scalar* T) {
  Scalar k = params[0];
  T[0] = Scalar(3) * k * u[0] * u[0];
}
```

Instantiate with `double` for the real evaluation; instantiate the same template
with an OTILib scalar (see `partner_kit/include/otilib_scalar.hpp`) for the
hypercomplex one. The Fortran template
(`templates/user_fortran_residual/residual.f90`) has the same three subroutines —
`eval_residual`, `eval_tangent`, `eval_dR_dparam` — with `params(1)`, `params(2)`
in `resasm.yml` key order.

Both templates additionally ship `eval_dR_dparam` (analytic `dR/da_i`) and a
**self-test** mode: run the binary with no arguments and it finite-differences
`eval_residual` and fails loudly if the analytic parameter derivatives have
drifted.

### A compiled program may write `response.json` instead of `response.npz`

This is the one concession the contract makes for compiled code, and it is
implemented in the reader.

`providers.py::_read_response(resp)` does:

1. if the `.npz` path exists → read it as a NumPy archive;
2. **otherwise strip the trailing `.npz` and open `<that>.json`** —
   `resp_json = resp_npz[:-4] + ".json"`.

So a binary handed `--response /tmp/xyz/response.npz` may simply write
`/tmp/xyz/response.json`, and the framework will find it. No NumPy archive
writer is needed in C++ or Fortran.

Both templates do exactly this:

```cpp
// The framework gives us "<dir>/response.npz"; we write "<dir>/response.json".
static std::string response_json_path(const std::string& resp) {
  std::string out = resp;
  if (out.size() > 4 && out.compare(out.size() - 4, 4, ".npz") == 0)
    out.erase(out.size() - 4);
  if (!(out.size() > 5 && out.compare(out.size() - 5, 5, ".json") == 0))
    out += ".json";
  return out;
}
```

```fortran
! templates/user_fortran_residual/residual.f90
outp = resp
n = len(outp)
if (n > 4) then
  if (outp(n-3:n) == '.npz') outp = outp(1:n-4)
end if
```

and then write

```
{"residual_coefficients_by_order": {"1": [[...]]}, "tangent": [[...]],
 "diagnostics": {"solver": "template_cpp", "method": "analytic", ...}}
```

Config is unchanged:

```yaml
residual:
  type: executable
  command: ./residual --request {request} --response {response}
tangent:
  type: response
```

Build them with the shipped `CMakeLists.txt` / `Makefile`. Copy either template
with `resasm init --template cpp --out my_job` or `--template fortran`.

---

## Checklist before you run

1. `resasm check resasm.yml` — it evaluates your residual, loads the tangent,
   generates the order-1 RHS and solves it, stopping at the first blocker.
2. Python path: confirm the derivatives are not silently zero — turn on
   `validation.rhs_finite_difference_check: true` and look at
   `public/summary.md`.
3. Compiled path: run the binary with no arguments (self-test) after every edit
   to `eval_residual`.
4. Black-box path: remember nothing verifies your equilibrium — the public
   summary will say so.
