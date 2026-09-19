# Example 8: second derivatives from your own solver (black box, order 2)

Your solver stays a black box. It is an executable that reads a request file
and writes back the Taylor coefficients of its residual. Residual_Assembler
solves for the first and second derivatives of the solution with respect to
two parameters. No OTILib is needed on this side: the solver carries the
derivatives itself.

| | |
| --- | --- |
| Folder | `templates/user_blackbox_order2_residual/` (what `resasm init --template blackbox-order2` copies) |
| Needs Abaqus | No |
| Needs OTILib | No |
| Run time | under 1 s |

## What it demonstrates

- **The black-box route.** Use it when you can run your solver but cannot, or
  may not, hand over its source. Residual_Assembler never sees the model. It
  sends the converged solution and the parameters and receives residual
  coefficients and the tangent. Your solver's own arithmetic produces the
  derivatives.
- **Coefficients, not derivatives.** At order 2 the solver must return Taylor
  coefficients, `coefficient = derivative / prod_i (kappa_i!)`. The template
  shows how to get this right without deriving anything by hand.
- **The order-by-order solve.** Before asking for order 2, the framework sends
  the solver the order-1 solution it has already solved, so the solver can
  include it in its evaluation.

## The mathematics

The residual inside `my_solver.py`, its unknown `u` and its parameters `k` and
`f`:

    R(u; k, f) = k^2 u^3 - f = 0,     converged at u = 2 for k = 2, f = 32,
    u(k, f) = (f / k^2)^(1/3),        T = dR/du = 3 k^2 u^2 = 48.

The model is `k^2 u^3` rather than the `k u^3` of Example 1 because the
second parameter derivatives of `k u^3 - f` are all zero. With `k^2`, the
order-2 residual coefficients are nonzero.

Differentiating `R(u(p), p) = 0` once and twice gives the equations solved
order by order:

    T U^(1) = -R^(1),     T U^(2) = -R^(2)(U^(1)).

Here `R^(p)` holds the order-`p` Taylor coefficients of `R(u*, p*)`. The
parameters are seeded as `k* = k + e1` and `f* = f + e2`, and
`u* = u + U^(1)` carries the lower orders already solved. The derivative of a
direction is its coefficient times its recovery factor: 1 for `e1`, `e2` and
`e1*e2`, and 2 for `e1^2` and `e2^2`. From the closed form of `u(k, f)`:

    du/dk = -2/3,      du/df = 1/48,
    d2u/dk2 = 5/9,     d2u/dk df = -1/144,     d2u/df2 = -1/2304.

The solver evaluates `R` once per order in a small truncated Taylor algebra
(`class T2` in `my_solver.py`) and reads the coefficients straight off. For
`e1^2` it returns `-40/3`, which is `(1/2!) d2R/dk2` and not the derivative
`d2R/dk2`.

## Inputs

| File | What it is |
| --- | --- |
| `resasm.yml` | the job: `residual: type: executable` with the command `python my_solver.py --request {request} --response {response}`, `tangent: type: response` (the solver returns it), parameters `k: 2.0` and `f: 32.0`, order 2 |
| `my_solver.py` | the black box. You edit `residual(u, k, f)`, which returns `[k*k*u[0]**3 - f]`, and `dR_du(u, k, f)`. The rest (the Taylor algebra and the file exchange) stays as it is |
| `solution.npy` | the converged solution, `[2.0]` |

The file exchange is specified in
[docs/blackbox_order2_contract.md](../../docs/blackbox_order2_contract.md):
`request.json` in, and `response.npz` out with `R_order_<p>` and `tangent`.

## Run it from the command line

From the root of the Residual_Assembler checkout, with the virtual environment
active ([INSTALL.md](../../docs/INSTALL.md)). OTILib is not needed:

```bash
export WORK="$HOME/resasm_work"; mkdir -p "$WORK"

resasm init   --template blackbox-order2 --out "$WORK/bb2"
resasm check  "$WORK/bb2/resasm.yml"
resasm run    "$WORK/bb2/resasm.yml"
resasm report "$WORK/bb2/resasm.yml"
```

Measured output of `resasm check` (exit 0):

```text
[ok] config parsed: blackbox_order2_demo (order 2, backend otilib)
[ok] loaded solution vector: shape (1,)
[ok] loaded parameter map: 2 parameters (k, f)
[ok] black-box executable configured: python my_solver.py --request {request} --response {response}
[warn] tangent will come from the black-box response
[ok] RHS order 1 generated: shape (1, 2)
[ok] sensitivity solve completed
```

The `[warn]` line is expected: this template takes the tangent from the
solver's response. `resasm run` printed (exit 0):

```text
run complete: {'parameters': ['k', 'f'], 'order': 2, 'tangent_source': 'response', 'residual_free_norm': None, 'orders_solved': [1, 2]}
```

`residual_free_norm` is `None` because a black box returns coefficients, not
the real residual. The public summary says so: equilibrium was not verified
by this run. `resasm report` then printed the same run with
`derivative orders : [1, 2]` and `validation status : ok` (exit 0).

To print the signed derivatives against the references:

```bash
python - "$WORK/bb2/resasm_output/private" <<'EOF'
import json, sys, numpy as np
private = sys.argv[1]
reference = {(1, 0): -2/3, (0, 1): 1/48, (2, 0): 5/9, (1, 1): -1/144, (0, 2): -1/2304}
for order in (1, 2):
    columns = json.load(open("%s/direction_map_order%d.json" % (private, order)))["columns"]
    arrays = np.load("%s/solution_sensitivities_order%d.npz" % (private, order))
    for column in columns:
        value = arrays["U_derivatives"][0, column["index"]]
        exact = reference[tuple(column["exponents"])]
        print("%-9s %+.17e  reference %+.17e  relative error %.1e"
              % (column["label"], value, exact, abs(value - exact) / abs(exact)))
EOF
```

## Run it in the GUI

Start `streamlit run scripts/app.py` ([GUI guide](../../docs/GUI_GUIDE.md)):

1. In the sidebar, set **Working directory** to a folder outside the
   repository and press **Create it** if needed.
2. Open **5. Job**. Type a **Job folder** name, for example `bb2`.
3. Choose `blackbox-order2` under **--template** and press **Copy template**.
4. Press **Run check**, **Run job** and **Read report**. The report opens
   `bb2/resasm_output/public/summary.md` below it, with a download button.

Measured on 2026-09-19, without OTILib: all four buttons exited 0 (0.18 s,
0.20 s, 0.38 s and 0.01 s).

## What it writes

Everything goes to `resasm_output/` inside the job folder, or to the folder
named by `output: dir:` in `resasm.yml`.

| Path | Contents | Share? |
| --- | --- | --- |
| `public/summary.md` | the run table, the note that equilibrium was not verified, the parameter ranking and the derivative convention | yes |
| `public/sensitivity_norms.csv` | one row per direction: label, OTI direction, recovery factor, norm of the solution sensitivity and of the right-hand side | yes |
| `public/parameter_ranking.csv` | parameters ranked by their first-order sensitivity | yes |
| `public/validation_summary.json` | status, and `residual_free_norm: null` with its reason | yes |
| `public/timing.json` | run time | yes |
| `private/solution_sensitivities_order{1,2}.npz` | `U_coefficients`, `U_derivatives`, `recovery_factors`, `direction_exponents` | no |
| `private/rhs_order{1,2}.npz` | the right-hand sides as coefficients and as derivatives | no |
| `private/direction_map_order{1,2}.json` | which column is which derivative, with its recovery factor | no |
| `private/tangent.npz`, `private/metadata.json`, `private/parameter_map.json`, `private/dof_map.json`, `private/validation_full.json` | the tangent from the response and the bookkeeping | no |
| `private/residual_real_UNAVAILABLE.json` | why no real residual is stored (a zero vector would falsely suggest equilibrium was checked) | no |

## Expected output (measured on 2026-09-19)

`public/sensitivity_norms.csv`:

```text
order,direction,oti_direction,recovery_factor,solution_sensitivity_norm,rhs_norm
1,d/dk,e1,1,6.666667e-01,3.200000e+01
1,d/df,e2,1,2.083333e-02,1.000000e+00
2,d2/dk2,e1^2,2,5.555556e-01,2.666667e+01
2,d2/dk_df,e1*e2,1,6.944444e-03,3.333333e-01
2,d2/df2,e2^2,2,4.340278e-04,2.083333e-02
```

Signed values from `private/` (the snippet above), against the closed form.
Acceptance: a relative error of at most 1e-12 for each derivative. The
largest measured error was 1.2e-16:

| Derivative | Coefficient returned | Recovery factor | Derivative | Reference | Relative error |
| --- | --- | --- | --- | --- | --- |
| du/dk | -6.6666666666666663e-01 | 1 | -6.6666666666666663e-01 | -2/3 | 0 |
| du/df | 2.0833333333333332e-02 | 1 | 2.0833333333333332e-02 | 1/48 | 0 |
| d2u/dk2 | 2.7777777777777779e-01 | 2 | 5.5555555555555558e-01 | 5/9 | 0 |
| d2u/dk df | -6.9444444444444449e-03 | 1 | -6.9444444444444449e-03 | -1/144 | 1.2e-16 |
| d2u/df2 | -2.1701388888888888e-04 | 2 | -4.3402777777777775e-04 | -1/2304 | 0 |

The order-2 right-hand-side coefficients are `[-13.333333, 0.333333,
0.010416667]`, that is `-40/3`, `1/3` and `1/96`, as the solver's docstring
states.

## How the result is checked independently

- **Closed form.** Each derivative above is compared with the derivatives of
  `u(k, f) = (f/k^2)^(1/3)`, which neither the solver nor the framework uses.
  The largest relative error is 1.2e-16.
- **Finite differences of re-solved equilibria.** This script re-solves
  `k^2 u^3 = f` by Newton at perturbed parameters (step `h = 1e-3`) and
  differences the solutions. No Taylor algebra is involved:

  ```bash
  python - "$WORK/bb2/resasm_output/private" <<'EOF'
  import json, sys, numpy as np
  def solve(k, f, u=2.0):
      for _ in range(50):
          u -= (k * k * u ** 3 - f) / (3 * k * k * u ** 2)
      return u
  k, f, h = 2.0, 32.0, 1e-3
  fd = {(1, 0): (solve(k + h, f) - solve(k - h, f)) / (2 * h),
        (0, 1): (solve(k, f + h) - solve(k, f - h)) / (2 * h),
        (2, 0): (solve(k + h, f) - 2 * solve(k, f) + solve(k - h, f)) / h ** 2,
        (1, 1): (solve(k + h, f + h) - solve(k + h, f - h) - solve(k - h, f + h)
                 + solve(k - h, f - h)) / (4 * h * h),
        (0, 2): (solve(k, f + h) - 2 * solve(k, f) + solve(k, f - h)) / h ** 2}
  for order in (1, 2):
      columns = json.load(open("%s/direction_map_order%d.json" % (sys.argv[1], order)))["columns"]
      values = np.load("%s/solution_sensitivities_order%d.npz" % (sys.argv[1], order))["U_derivatives"][0]
      for c in columns:
          e = tuple(c["exponents"])
          print("%-9s relative difference %.1e" % (c["label"], abs(values[c["index"]] - fd[e]) / abs(fd[e])))
  EOF
  ```

  Measured: 1.9e-7 (`d/dk`), 1.8e-10 (`d/df`), 2.0e-7 (`d2/dk2`), 1.8e-7
  (`d2/dk_df`) and 3.4e-7 (`d2/df2`). This is the size of the `h^2`
  truncation error of the differences.
- **The convention is tested.** `tests/framework/test_blackbox_order2_coefficients.py`
  runs the template end to end. It also checks that a solver returning
  derivatives instead of coefficients gives a wrong answer
  (`python -m pytest -q tests/framework/test_blackbox_order2_coefficients.py`:
  7 passed in 0.8 s).

## Run time

Measured on 2026-09-19 (Linux, Python 3.11, without OTILib): `resasm init`
0.41 s, `resasm check` 0.57 s, `resasm run` 0.73 s wall time, of which the job
itself took 0.28 s (`public/timing.json`).

## Common problems

| Symptom | Cause and fix |
| --- | --- |
| `d2u/dk2` comes out as `1.1111111111111112` and `d2u/df2` as `-8.6805555555555551e-04`, twice the references, while `d2u/dk df` is right | The solver returns derivatives instead of Taylor coefficients (measured by multiplying the returned `e1^2` and `e2^2` columns by 2!). Return `(1/2!) d2R/dk2`. The template's `T2` algebra does this for you, so do not rescale its output. |
| `ERROR: black-box executable produced no response file.` with a traceback ending in `TypeError: float() argument must be a string or a real number, not 'T2'`, exit 2 | The edited `residual()` casts a parameter with `float()`, or calls a NumPy function on it, which destroys the Taylor number. Use only `+ - * **` on `u`, `k` and `f`. |
| `[warn] tangent will come from the black-box response` | Expected: `tangent: type: response` takes the tangent from `response.npz`. |
| `residual norm (free) : None` in `resasm report` | Expected for a black box: the solver returns coefficients, not the real residual, so equilibrium is not checked here. Check convergence in your own solver. |
| `this reference template implements order <= 2` | `sensitivity: order:` in `resasm.yml` is above 2. Extend the `T2` algebra, or use the direct Python route with OTILib ([Example 1](../../examples/user_config_minimal/WALKTHROUGH.md)). |

Previous: [Example 7](../../examples/finite_strain_c3d8/WALKTHROUGH.md).
All examples: [examples/README.md](../../examples/README.md).
