# Example 6: the provider-to-sensitivity pipeline in one command

One command builds the compiled J2 material from the companion repository,
solves a small cyclic finite-element history with it, replays the history
with the OTI material, and verifies every derivative against finite
differences of the unmodified material. It needs neither Abaqus nor OTILib.

| | |
| --- | --- |
| Folder | `examples/bounded_j2_c3d8/` |
| Needs Abaqus | No |
| Needs OTILib | No |
| Run time | about 10 s |

## What it demonstrates

The two repositories work as one system. UMAT-OTI turns a UMAT into a
provider object that computes derivatives exactly; Residual_Assembler replays
that object over a history and solves for the sensitivities. This example
exercises both halves from a clean start and checks the result with an
independent method, so it is the quickest way to confirm that an installation
of both packages is complete.

It also shows why the *total-history* derivative matters. Under cyclic
loading, the plastic state at the end depends on the whole path. A derivative
that ignores how earlier increments change with the parameter (a "fixed-path"
local solve) is wrong as soon as plasticity is involved, and the plot the
example writes shows by how much.

## The mathematics

At every increment `n` the equilibrium residual on the free degrees of freedom
is zero, `R_n(u_n(p), p) = 0`, with parameter-independent loads. Its total
derivative gives

    K_n du_n/dp = - dR_n/dp |_(u_n fixed),

where the right-hand side contains the derivative of the stress with respect
to `p` *including* its dependence, through the stored state, on `du_k/dp` for
all earlier increments `k < n`. The compiled provider carries these history
derivatives exactly (first-order OTI).

The reference re-solves the whole model with the ORIGINAL routine at
`p (1 + h)` and `p (1 - h)` and forms `(u(p + h p) - u(p - h p)) / (2 h p)` for
`h` = 1e-4, 3e-5 and 1e-5. The derivative passes when the scaled error is
below 2e-6 at all three steps and the last two steps agree (a plateau).

## Inputs

| File | What it is |
| --- | --- |
| `model.json` | one unit-cube C3D8 (eight integration points, selective-reduced integration), bottom nodes 1 to 4 clamped, top nodes 5 to 8 loaded by 100 in z times the load factors 0.1, 0.4, 1.0, 1.2, 0.8, -1.0, -1.1 (seven increments: loading, unloading and reversed loading). Material constants `[210000, 0.3, 250, 2000]` for `E, nu, SIGY0, H` |
| the m3_j2 contract of UMAT-OTI | `parameter_sensitivity/models/m3_j2/contract_v2.json` in the companion repository, found automatically when the two repositories sit side by side |

The history is a converged finite-element solution computed by the script
itself with the ORIGINAL material; it is labelled `synthetic_converged_fe`, not
an Abaqus result. A genuine archived elastic Abaqus export is checked in the
same run.

## Run it from the command line

From the Residual_Assembler root, with both packages installed
([INSTALL.md](../../docs/INSTALL.md)):

```bash
export WORK="$HOME/resasm_work"; mkdir -p "$WORK"
python scripts/reproduce_connected_pipeline.py --skip-abaqus --out "$WORK/pipeline"
```

Measured output (exit 0):

```text
verified bounded J2 pipeline: .../pipeline/manifest.json
```

and `manifest.json`:

```json
{
  "passed": true,
  "private_manifest": "private/manifest.json",
  "public_results": "public",
  "error": null
}
```

`--out` must be a new folder. `--skip-abaqus` uses the archived elastic
Abaqus export instead of starting an Abaqus job. By default
(`--imports installed`) the script runs its steps with the installed packages
in an isolated Python. If you work from source checkouts through `PYTHONPATH`
instead of an installation, add `--imports environment`.

The script runs, in order: `umat-oti-provider build` on the m3_j2 contract,
then `resasm replay model.json --object ... --contract ... --solve --verify`,
then a second replay of the solved record, then the elastic Abaqus fixture
check. You can run the central step yourself:

```bash
umat-oti-provider build ../UMAT_source_transformation/parameter_sensitivity/models/m3_j2/contract_v2.json \
    --out "$WORK/provider_j2"
resasm replay examples/bounded_j2_c3d8/model.json \
    --object "$WORK/provider_j2/umat_m3_j2_oti.obj" --contract "$WORK/provider_j2/umat_m3_j2_oti.json" \
    --solve --verify --out "$WORK/replay_j2"
```

Measured: `replay: 7 increments, 8 IPs; total-history du/dp; verified=True`
(2.1 s), with the summary in `public/summary.json`.

## Run it in the GUI

The pipeline script has no GUI of its own, but its central step does. Start
`streamlit run scripts/app.py` ([GUI guide](../../docs/GUI_GUIDE.md)) and open
**Advanced Replay**:

1. **Record or model JSON** is already `examples/bounded_j2_c3d8/model.json`.
2. **Compiled provider object**: `$WORK/provider_j2/umat_m3_j2_oti.obj`.
   **Provider contract JSON**: `$WORK/provider_j2/umat_m3_j2_oti.json`.
3. **Replay output directory**: a new folder outside the repository (the
   default is inside the checkout, in the ignored `resasm_gui_workspace/`).
4. Keep **Solve model** and **Verify with ORIGINAL finite differences**
   ticked, and press **Run replay**.

Measured on 2026-09-18: exit code 0 in 1.60 s, output
`replay: 7 increments, 8 IPs; total-history du/dp; verified=True`.

## What it writes

| Path in `$WORK/pipeline/` | Contents | Share? |
| --- | --- | --- |
| `manifest.json` | `passed`, `error`, and where the rest is | yes |
| `public/verification.json` | status, `verified`, parameters, increments, points, provenance, and one row per parameter and step with the scaled errors of displacement, stress and state (total history and fixed path) | yes |
| `public/fd_table.csv` | the same rows as a table | yes |
| `public/j2_history.png` | left: the vertical displacement of the last node over the seven increments; right: the total-history `du/dnu` against the fixed-path local solve | yes |
| `public/abaqus_fixture.json` | the elastic Abaqus export check (errors and limits) | yes |
| `private/manifest.json` | environment, Python and compiler versions, repository commits, source hashes, every command with its log | no: it records local paths |
| `private/command-*.log` | the output of each command | no |
| `private/provider/` | the freshly built object, its mapping and its build folder | no |
| `private/solve/`, `private/replay/` | the solved record and the replay results, each with `public/` and `private/` | no |

## Expected output (measured on 2026-09-18)

`public/fd_table.csv`, largest scaled error of the total-history derivatives
over all parameters and the three steps:

| Quantity | Largest scaled error | Tolerance |
| --- | --- | --- |
| displacement `du/dp` | 1.02e-8 (E, h = 1e-4) | 2e-6 |
| stress `dsigma/dp` | 9.72e-7 (H, h = 1e-5) | 2e-6 |
| state `d eqplas/dp` | 2.84e-8 (E, h = 1e-5) | 2e-6 |

`public/abaqus_fixture.json` (genuine Abaqus 2021 elastic export,
`passed: true`):

| Check | Error | Tolerance |
| --- | --- | --- |
| `du/dE` against Abaqus finite differences | 2.08e-5 | 3e-5 |
| `du/dnu` against Abaqus finite differences | 5.45e-6 | 1e-5 |
| integration-point stress reconstructed from the Abaqus displacements | 7.6e-8 | 1e-4 |
| equilibrium with the Abaqus stresses | 1.6e-8 | 1e-5 |
| the same with the integration points deliberately permuted (must fail) | 0.475 | must exceed 0.1 |

The plot shows the point of the example: up to increment 4 the total-history
and fixed-path `du/dnu` coincide; at increment 5, where unloading begins, the
fixed-path value jumps to about -0.0018 while the total-history value is about
-0.0005, and the two never meet again. Open `public/j2_history.png` in your
output folder to see it.

## How the result is checked independently

- **Whole-model finite differences of the ORIGINAL routine**, described above:
  the finite differences never touch the OTI code, and the model is
  re-equilibrated at every perturbed parameter value.
- **A genuine Abaqus result** (the elastic fixture): displacement derivatives
  against Abaqus finite differences, the stress reconstruction and the
  equilibrium check, with a negative control that must fail.
- **Negative result shown on purpose**: the fixed-path local solve is written
  next to the total-history derivative so that the difference is visible.

## Run time

Measured on 2026-09-18: 10.7 s with `--imports environment` in a development
setup, and 9.9 s with the default `--imports installed` in a fresh virtual
environment with both packages installed and no OTILib. The central
`resasm replay` step alone takes 1.6 to 2.1 s.

## Common problems

| Symptom | Cause and fix |
| --- | --- |
| `reproduction failed: command failed (1): [... '-I', '-m', 'umat_oti.provider', 'build', ...]` with `No module named umat_oti.provider` | The isolated Python does not see the UMAT-OTI package you intend to use (for example, you work through `PYTHONPATH`, or an older installation is active). Install both packages ([INSTALL.md](../../docs/INSTALL.md)) or add `--imports environment`. |
| `reproduction failed: Command '['git', '-C', ..., 'rev-parse', 'HEAD']' returned non-zero exit status 128` | The script records the commit of each repository, so both must be git clones, not unpacked archives. |
| `--out must be a fresh directory` | Choose a new `--out`; a failed run leaves its folder behind. |
| The provider build fails | `gfortran` is missing ([INSTALL.md](../../docs/INSTALL.md)). |
| `manifest.json` says `"passed": false` (exit 2) | Its `error` field and `private/command-*.log` name the failing step. |

This example is a verification of the algebra on one element. It makes no
claim about large deformations (the cyclic loads accumulate about 10 %
displacement on the unit cube) and it is not an Abaqus J2 analysis; for those, see
[Example 3](../presentation_request/WALKTHROUGH.md) and
[Example 5](../cantilevers/WALKTHROUGH.md). All examples:
[examples/README.md](../README.md).
