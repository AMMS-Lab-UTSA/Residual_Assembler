# Example 7: finite-strain neo-Hookean assembly and sensitivities

Two distorted C3D8 elements that share a face, stretched and rotated far
beyond small strain. Residual_Assembler assembles the finite-strain residual
and its exact tangent, and computes the derivatives of the solution with
respect to the two material constants with OTILib. Every step is checked
against an independent calculation.

| | |
| --- | --- |
| Folder | `examples/finite_strain_c3d8/` |
| Needs Abaqus | No (an optional Abaqus check is recorded in `verified/abaqus_report.json`) |
| Needs OTILib | Yes ([installation](../../docs/INSTALL.md#4-install-otilib-for-the-direct-residual-paths)) |
| Run time | about 3 s for the benchmark, under 1.5 s per command |

## What it demonstrates

- **Finite-strain assembly.** The kinematics use the deformation gradient, the
  stress is the Cauchy stress, and the tangent includes the geometric
  stiffness. The two elements are distorted, share nodes, and carry a large
  superposed rotation (0.7 rad), so errors in the element geometry or the
  tangent cannot hide.
- **Material-parameter sensitivities through the same assembly.** The
  residual is evaluated with OTI numbers in place of `mu` and `lambda`, so the
  derivatives come from the same code path as the residual itself.
- **The ordinary commands.** It uses `resasm assemble` and
  `resasm sensitivity`, not the request interface of Examples 3 to 5.

## The mathematics

The compressible neo-Hookean strain energy (dimensionless benchmark units):

    W(F) = mu/2 (tr(F F^T) - 3) - mu ln J + lambda/2 (ln J)^2,     J = det F > 0,
    mu = 2.3,  lambda = 4.1.

The residual `R(u; mu, lambda) = F_int(u; mu, lambda) - F_ext` is assembled
with full integration. The external nodal loads are manufactured: they are
computed once, by an independent reference-volume first-Piola quadrature, at a
prescribed target solution, and they stay fixed when the parameters change.
So the target solution is an exact equilibrium, and differentiating
`R(u(p), p) = 0` gives

    K du/dp = - dR/dp,     p in {mu, lambda},

with `K` the exact linearisation (material plus geometric stiffness).
`dR/dp` is read from the OTI evaluation of `R`.

## Inputs

`benchmark.py` writes all inputs into its output folder, and a copy of them is
committed in `verified/`:

| File | What it is |
| --- | --- |
| `model.json` | 12 nodes, 2 C3D8 elements bound to `solid_c3d8_finite_strain`, material `compressible_neo_hookean` with constants `[2.3, 4.1]`, fixed supports and 12 concentrated loads; 36 DOF, 12 free |
| `config.json` | `mode: material-replay` and the target solution (36 values) under `options.solution` |
| `params.json` | the parameters `solid.mu` and `solid.lambda`, order 1, backend `otilib` |
| `report.json` | the benchmark's own verification results |

## Run it from the command line

From the Residual_Assembler root, with OTILib visible
([INSTALL.md](../../docs/INSTALL.md)):

```bash
export WORK="$HOME/resasm_work"; mkdir -p "$WORK"
export PYOTI_PATH=/path/to/otilib/build OTILIB_ROOT=/path/to/otilib/build

python examples/finite_strain_c3d8/benchmark.py --out "$WORK/finite"

F="$WORK/finite"
resasm --config "$F/config.json" assemble    "$F/model.json" --mode material-replay --tangent
resasm --config "$F/config.json" sensitivity "$F/model.json" --params "$F/params.json" --out "$F/sens"
```

`--config` is a global option: it goes before the subcommand. You may also use
the committed copies in `examples/finite_strain_c3d8/verified/` in the same
way; only `--out` writes files.

Measured output. The benchmark prints `report.json`, ending in
`"passed": true` (exit 0). The model written is identical to
`verified/model.json`. Then:

```text
assembled mode 'material-replay': ndof=36  ||R||=1.744643e+00  max|R|=5.938249e-01
```

```text
hypercomplex backend     : otilib
basis count (m)          : 2
truncation order (nt)    : 1
total coefficients (N)   : 3
sensitivity analysis (mode=material-replay, order=1)
  residual norm ||R_free|| = 1.056172e-15
  tangent source           = backend-assembled
  parameters               = ['solid.mu', 'solid.lambda']
  R^(1) shape              = (36, 2)  (2 directions)
  hypercomplex ready       = True
  solved derivative orders : [1]
  order 1:
    d^1/e1         = +1.872427e-01   [FD +1.872427e-01, rel 3.07e-10]
    d^1/e2         = +3.398245e-02   [FD +3.398245e-02, rel 7.12e-10]
```

`||R|| = 1.744643` is not an equilibrium error: the full residual includes the
support reactions. On the free degrees of freedom it is `1.06e-15`. The lines
`d^1/e1` and `d^1/e2` are the Euclidean norms of `du/dmu` and `du/dlambda`
over all 36 DOF; the arrays are in the saved package. `rel` is the relative
error of the whole 36-entry column against the re-solved finite differences;
if it reached `1e-4` for any parameter, the command would say
`finite-difference check FAILED` and exit 1.

## Run it in the GUI

Start `streamlit run scripts/app.py` ([GUI guide](../../docs/GUI_GUIDE.md)):

1. Set the sidebar's **Working directory** to a folder outside the repository.
2. On **1. Model**, open **Your own model**, type the absolute path of
   `examples/finite_strain_c3d8/verified/model.json` and press **Use this path**.
3. On **4. Sensitivity**, the `--params` field is filled with the
   `params.json` beside the model. Choose `--mode material-replay` (or leave it
   to the params file) and press **Run sensitivity**. Measured: the same two
   derivatives, `+1.872427e-01` and `+3.398245e-02`, with finite-difference
   agreement 2.49e-10 and 6.79e-10.

The GUI tabs do not pass `--config`. The sensitivity command then finds the
equilibrium itself (measured free residual 4.5e-15) and gives the same
derivatives, but **3. Assemble** evaluates R at zero displacement
(measured `||R||=1.553269e+00`) instead of at the target solution. Use the
command line for the assembly step of this example.

## What it writes

| Path | Contents |
| --- | --- |
| `$WORK/finite/model.json`, `config.json`, `params.json` | the inputs, as described above |
| `$WORK/finite/report.json` | the benchmark's checks and their numbers (below) |
| `$WORK/finite/sens.json`, `sens.md`, `sens.npz` | the sensitivity package: what was run, whether it is runnable, the tangent source, and the arrays |
| `sens_residual.*`, `sens_sensitivity.*`, `sens_state.*`, `sens_validation.*` | the parts of the package: the residual and reactions, the right-hand sides `R^(1)` and their direction map, the state, and the validation report |

These files contain the full model and its arrays. There is no separate
public folder for this command; share only what the model owner allows.

## Expected output (measured on 2026-09-18)

`report.json`:

| Check | Steps | Relative error | Acceptance |
| --- | --- | --- | --- |
| residual against the independent first-Piola quadrature | fixed state | 5.70e-16 | 1e-12 |
| free residual at the target solution | | 1.06e-15 (absolute) | |
| complete 36-column nodal tangent against finite differences | 1e-4, 1e-5, 1e-6 | 4.44e-9, 4.55e-11, 1.09e-10 | 2e-8 at every step |
| OTILib `du/d(mu, lambda)` against nonlinear central re-solves | 1e-3, 1e-4, 1e-5 (relative) | 9.48e-7, 9.48e-9, 1.01e-10 | 2e-5 at every step |

`sens_validation.json`: `residual-vs-solver` 1.06e-15 (pass),
`sensitivity-vs-finite-difference` 3.27e-10 (pass, tolerance 1e-4); the
reaction check is marked not applicable because no exported reactions were
given.

## How the result is checked independently

- **The residual** is compared with an independent first-Piola quadrature
  over the reference volume, a different formulation of the same integral
  (5.7e-16).
- **The tangent** is compared column by column with finite differences of the
  residual at three step sizes; the errors fall with the step and then level
  off at the roundoff floor, which is what a correct tangent does.
- **The sensitivities** are compared with the solutions of the full nonlinear
  problem at perturbed `mu` and `lambda` (central differences at three steps),
  not merely with finite differences of the right-hand side.
- **Abaqus.** An independent UHYPER (`neo_hookean_uhyper.for`) and one Abaqus
  job were used to check reactions and one tangent action; the result is
  recorded in `verified/abaqus_report.json` and described in the
  [README](README.md). It was not re-run for this guide.

## Run time

Measured on 2026-09-18: benchmark 2.7 s, `assemble` 0.5 s, `sensitivity`
1.2 s (wall time, including Python start-up).

## Common problems

| Symptom | Cause and fix |
| --- | --- |
| The benchmark ends with `RuntimeError: OTILib backend requested but genuine OTILib was not found` (exit 1) | Set `PYOTI_PATH` and `OTILIB_ROOT` to your OTILib build. Never `pip install pyoti`. |
| `resasm sensitivity` prints the same message and exits 3 | The same cause. `--backend dual1` is not an alternative for this model: it cannot seed the material constants of the finite-strain formulation and refuses with `ERROR: backend='dual1' cannot differentiate the finite-strain C3D8 formulation (solid_c3d8_finite_strain): only the OTILib backend seeds its material constants. Use backend='otilib'.` (exit 2). |
| `resasm: error: unrecognized arguments: --config` | `--config` must come before the subcommand: `resasm --config config.json assemble ...`. |
| `||R||` is 1.553269 instead of 1.744643 | The command ran without `--config`, so R was evaluated at zero displacement. |
| A different material or a history-dependent law is refused | This path supports only isotropic, stateless total hyperelasticity, first-order material parameters and fixed (non-follower) loads. The refusal names the missing capability. |

Previous: [Example 6](../bounded_j2_c3d8/WALKTHROUGH.md).
All examples: [examples/README.md](../README.md).
