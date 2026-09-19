# Example 1: the smallest sensitivity calculation

A cubic spring with one unknown. You write its residual in ordinary Python,
and Residual_Assembler returns the first and second derivatives of the
solution with respect to both parameters. Everything can be checked by hand.

| | |
| --- | --- |
| Folder | `examples/user_config_minimal/` (the same job that `resasm init --template python` creates) |
| Needs Abaqus | No |
| Needs OTILib | Yes ([installation](../../docs/INSTALL.md#4-install-otilib-for-the-direct-residual-paths)) |
| Run time | about 1 s |

## What it demonstrates

Many models end in the same question: *if I change a parameter slightly, how
much does the converged solution move?* This example answers that question for
the simplest nonlinear problem there is, so that you can see every step of the
machinery on numbers you can check with a pencil.

It also shows the working pattern of the direct-residual route. You supply a
function `residual(u, params)` and a converged solution. Residual_Assembler
evaluates the function once with ordinary numbers and once with hypercomplex
(OTI) numbers, and reads the derivatives from the result. You never write a
derivative yourself.

Use this route when you can evaluate the global residual of your own model.
For an Abaqus analysis, where the residual is not available, use the request
route instead ([Example 3](../presentation_request/WALKTHROUGH.md)).

## The mathematics

The residual of the spring, its unknown `u` and its parameters `k` and `f`:

    R(u; k, f) = k u^3 - f = 0,        converged at u = 2 for k = 2, f = 16.

Differentiating `R(u(p), p) = 0` with respect to a parameter `p` gives the
sensitivity equation that the program solves:

    (dR/du) du/dp = - dR/dp,           dR/du = 3 k u^2 = 24,
    du/dk = - u^3 / 24 = -1/3,         du/df = 1/24.

Because `u = (f/k)^(1/3)` in closed form, the second derivatives are known too:

    d2u/dk2 = 4u/(9k^2) = 2/9,   d2u/dk df = -u/(9kf) = -1/144,   d2u/df2 = -2u/(9f^2) = -1/576.

The program does not use these formulas. It evaluates `R` with
`k + e1` and `f + e2`, where `e1` and `e2` are OTI imaginary directions
truncated at order 2, and solves `T U^(p) = -R^(p)` order by order (`T` is the
tangent `dR/du`). OTI returns Taylor coefficients. The derivative is the
coefficient times a recovery factor, the product of the factorials of the
exponents: 1 for first order and for `e1*e2`, and 2 for `e1^2` and `e2^2`.
Both forms are saved, so you never have to apply the factor yourself.

## Inputs

| File | What it is |
| --- | --- |
| `resasm.yml` | the job: one unknown, parameters `k: 2.0` and `f: 16.0`, order 2, backend `otilib`, and a finite-difference check of the right-hand side |
| `user_residual.py` | `residual(u, params)` returns `[k*u[0]**3 - f]`; the optional `tangent` returns `[[3*k*u[0]**2]]` |
| `solution.npy` | the converged solution, `[2.0]` |

## Run it from the command line

From the root of the Residual_Assembler checkout, with the virtual environment
active (see [INSTALL.md](../../docs/INSTALL.md)):

```bash
export WORK="$HOME/resasm_work"            # all outputs go here, outside the repository
mkdir -p "$WORK"
export PYOTI_PATH=/path/to/otilib/build    # your OTILib build directory
export OTILIB_ROOT="$PYOTI_PATH"

resasm init --template python --out "$WORK/spring"
resasm check  "$WORK/spring/resasm.yml"
resasm run    "$WORK/spring/resasm.yml"
resasm report "$WORK/spring/resasm_output"
```

`resasm init` copies the template into a new folder. Paths inside
`resasm.yml` are resolved against the folder that holds it, so you can run the
three commands from anywhere. To run the committed copy in this folder instead,
copy it first, so that no output is written into the repository:

```bash
cp -r examples/user_config_minimal "$WORK/spring_copy"
resasm run "$WORK/spring_copy/resasm.yml"
```

`resasm check` prints one line per requirement and stops at the first `[fail]`.
On the measured run every line was `[ok]`, including
`[ok] OTILib available: order 2, basis 2` and
`[ok] residual norm on free DOFs is 0.000e+00`.

`resasm run` then printed:

```text
run complete: {'parameters': ['k', 'f'], 'order': 2, 'tangent_source': 'python', 'residual_free_norm': 0.0, 'orders_solved': [1, 2]}
```

and `resasm report` summarised the run:

```text
  residual norm (free) : 0.0
  tangent source       : python
  parameters           : k, f
  derivative orders    : [1, 2]
  validation status    : ok
```

To see the signed derivatives, read the private arrays:

```bash
python - "$WORK/spring/resasm_output/private" <<'EOF'
import sys, numpy as np
p = sys.argv[1]
for order in (1, 2):
    d = np.load(f"{p}/solution_sensitivities_order{order}.npz")
    print(order, d["U_derivatives"][0], "coefficients:", d["U_coefficients"][0])
EOF
```

## Run it in the GUI

Start the GUI with `streamlit run scripts/app.py` ([GUI guide](../../docs/GUI_GUIDE.md)).

1. In the sidebar, set **Working directory** to a folder outside the
   repository (for example `$HOME/resasm_work`) and press **Create it** if
   needed.
2. Open the **5. Job** tab. Type a **Job folder** name, for example `spring`.
3. Choose `python` under **--template** and press **Copy template**.
4. Press **Run check**, then **Run job**. Each shows the exact `resasm` command
   and its exit code.
5. Press **Read report**. It runs `resasm report spring/resasm.yml`, which
   reads the folder the job wrote (`spring/resasm_output`, or the folder named
   by `output: dir:` in `resasm.yml`). The `public/summary.md` panel opens
   below it, with a download button. Measured on 2026-09-18: all four buttons
   exited 0.

The same spring is also available as a built-in element: on **1. Model** pick
`examples/minimal_nonlinear_spring_sensitivity`, then on **4. Sensitivity**
press **Run sensitivity**. The measured output was
`d^1/e1 = -3.333333e-01 [analytic -3.333333e-01, rel 0.00e+00] [FD -3.333330e-01, rel 1.00e-06]`.

## What it writes

Everything goes to `resasm_output/` inside the job folder.

| Path | Contents | Share? |
| --- | --- | --- |
| `public/summary.md` | run table, parameter ranking, the convention for derivatives, the finite-difference check | yes |
| `public/sensitivity_norms.csv` | one row per derivative direction: label, OTI direction, recovery factor, norm of the solution sensitivity, norm of the right-hand side | yes |
| `public/parameter_ranking.csv` | parameters ranked by the size of their first-order sensitivity | yes |
| `public/validation_summary.json` | machine-readable status and the finite-difference check | yes |
| `public/timing.json` | run time | yes |
| `private/solution_sensitivities_order{1,2}.npz` | `U_coefficients`, `U_derivatives`, `recovery_factors`, `direction_exponents` | no |
| `private/rhs_order{1,2}.npz`, `private/residual_real.npz`, `private/tangent.npz` | the right-hand sides, the residual and the tangent | no |
| `private/direction_map_order{1,2}.json` | which column is which derivative, with its recovery factor | no |
| `private/metadata.json`, `private/parameter_map.json`, `private/dof_map.json`, `private/validation_full.json` | bookkeeping and the full validation record | no |

The `public/` folder holds norms, rankings and status only. It is the folder
you may send to someone else. The `private/` folder holds the full arrays.

## Expected output (measured on 2026-09-18)

`public/sensitivity_norms.csv`:

```text
order,direction,oti_direction,recovery_factor,solution_sensitivity_norm,rhs_norm
1,d/dk,e1,1,3.333333e-01,8.000000e+00
1,d/df,e2,1,4.166667e-02,1.000000e+00
2,d2/dk2,e1^2,2,2.222222e-01,5.333333e+00
2,d2/dk_df,e1*e2,1,6.944444e-03,1.666667e-01
2,d2/df2,e2^2,2,1.736111e-03,4.166667e-02
```

Signed values from `private/`, against the closed form:

| Derivative | Program | Closed form | Difference |
| --- | --- | --- | --- |
| du/dk | -3.3333333333333331e-01 | -1/3 | 0 |
| du/df | 4.1666666666666664e-02 | 1/24 | 0 |
| d2u/dk2 | 2.2222222222222224e-01 | 2/9 | 2.8e-17 |
| d2u/dk df | -6.9444444444444449e-03 | -1/144 | 8.7e-19 |
| d2u/df2 | -1.7361111111111110e-03 | -1/576 | 0 |

The raw coefficient of `e1^2` is `1.1111111111111112e-01` (1/9). Multiplied by
its recovery factor 2 it gives the derivative 2/9. This is why the public
report always quotes recovered derivatives.

## How the result is checked independently

- **Closed form.** The table above compares every derivative with the
  analytic value. The largest difference is 2.8e-17.
- **Finite differences of the right-hand side.** With
  `rhs_finite_difference_check: true`, the run differentiates the residual
  with respect to each parameter by finite differences at fixed `u` and solves
  with the same tangent. Measured maximum relative difference: 2.676e-11
  (`status: pass`). This checks the OTI right-hand side and the solve. It does
  not re-solve the nonlinear problem, so it cannot detect an error in the
  tangent; the closed form above covers that.

## Run time

Measured on 2026-09-18 (Linux, Python 3.11): `resasm init` 0.5 s, `resasm run`
0.9 s wall time, of which the sensitivity job itself took 0.48 s
(`public/timing.json`).

## Common problems

| Symptom | Cause and fix |
| --- | --- |
| `resasm check` prints `[fail] OTILib backend requested but genuine OTILib was not found` and exits 1 | OTILib is not visible. Set `PYOTI_PATH` and `OTILIB_ROOT` to your OTILib build directory ([INSTALL.md](../../docs/INSTALL.md#4-install-otilib-for-the-direct-residual-paths)). Never `pip install pyoti`: that PyPI name is an unrelated package. |
| `resasm run` prints `ERROR: OTILib backend requested but genuine OTILib was not found.` and exits 3 | The same cause as above. |
| `ERROR: ... already exists (use --force to overwrite)`, exit 2 | `resasm init` never overwrites a folder. Choose a new `--out`, or add `--force` if you really want to replace it. |
| Derivatives are zero or the evaluation fails after you edit the residual | Keep to ordinary arithmetic (`+ - * / **`) on `u` and `params`, as the comments in `user_residual.py` say: the same function is called with OTI numbers. |
| You need derivatives without OTILib | Use the black-box templates (`resasm init --template blackbox-order2`), where your own executable returns the Taylor coefficients ([CLI guide](../../docs/CLI_GUIDE.md#resasm-init)). |

Next: [Example 2](../../residual_core/examples/minimal_c3d8_stress_driven/WALKTHROUGH.md),
a residual assembled from a finite-element mesh. All examples:
[examples/README.md](../README.md).
