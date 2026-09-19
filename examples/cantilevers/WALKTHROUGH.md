# Example 5: full-size cantilevers, J2 and FCC crystal plasticity

Two nonlinear cantilevers of realistic size, bent into the plastic range by a
prescribed tip displacement. From each Abaqus analysis the program returns,
at every one of 40 (J2) or 25 (FCC) increments, the sensitivities of the tip
reaction, of stresses and plastic strain, and the full von Mises field with
its derivative with respect to every material parameter, together with the
weighted share of each parameter.

| | J2 cantilever | FCC crystal cantilever |
| --- | --- | --- |
| Mesh | 48 x 16 x 2 C3D8: 1,536 elements, 2,499 nodes, 7,497 DOF | 24 x 8 x 2 C3D8: 384 elements, 675 nodes, 2,025 DOF |
| Integration points | 12,288 | 3,072 |
| Increments | 40 (tip pushed 0.7 mm) | 25 (tip pushed 0.25 mm) |
| Parameters | E, nu, SIGY0, H = 200000, 0.3, 250, 2000 | C11, C12, C44, g0, gsat, h0, a, q, gd0, m = 168000, 121000, 75000, 13, 55, 800, 2, 1.4, 0.001, 0.05 |
| Material | ORIGINAL UMAT `parameter_sensitivity/models/m3_j2/umat.for` of UMAT-OTI | `parameter_sensitivity/models/m6_fcc/umat.for`, 12 slip systems {111}<110>, 12 state variables |

| | |
| --- | --- |
| Folder | `examples/cantilevers/` |
| Needs Abaqus | Yes, once: to run the two analyses and read their ODBs. The sensitivity computation itself does not call Abaqus. |
| Needs OTILib | No |
| Run time | J2 25 to 40 s, FCC 25 to 80 s (see below) |

## What it demonstrates

- **Scale.** The same four-file interface as [Example 3](../presentation_request/WALKTHROUGH.md)
  on models with thousands of degrees of freedom, tens of increments and up to
  ten parameters. Assembly is sparse, and each increment's tangent is factorised
  once for all parameters.
- **Any compiled material.** The J2 and the crystal-plasticity UMATs are
  handled by the same engine; only the provider object differs.
- **Full-field sensitivities.** `dMISES/dp` at every integration point and
  increment, stored in `fields.npz`, together with `U`, `RF`, `S`, `SDV` and
  their derivatives.
- **Weighted per-increment shares**, which show which parameter governs the
  response as the loading moves from elastic to plastic.

## The mathematics

The engine is the one of [Example 4](../replay_history/WALKTHROUGH.md): at each
increment `K_ff du_f/dp = -dR_f/dp` with the history carried from increment to
increment, then

    dsigma_n/dp = dsigma^A/dp + DDSDDE B du_n/dp,     dRF/dp = dR_c/dp + K_cf du_f/dp,
    dMISES/dp   = (d MISES / d sigma) : dsigma_n/dp.

Weighted share of parameter `j` at increment `n` (von Mises field, volumes `V_i`):

    W_j(n) = sum_i V_i |p_j dMISES_i/dp_j| / sum_i V_i,       share_j = 100 W_j / sum_k W_k.

`sensitivity_shares.csv` gives this *field share* and, next to it, the *scalar
share*, the same formula applied to the volume-mean von Mises stress.

Independent check without finite differences: both materials are homogeneous
of degree one in their stress-dimensioned parameters, J2 in `(E, SIGY0, H)` at
fixed `nu`, and the FCC crystal in `(C11, C12, C44, g0, gsat, h0)`, because its
slip rate depends only on the ratio of resolved shear stress to slip
resistance. Under prescribed displacements, at every increment,

    sum_p p dQ/dp = Q    for reactions, stresses and von Mises,
    sum_p p dQ/dp = 0    for displacements and plastic strain.

## Inputs

| File | What it is |
| --- | --- |
| `gen_cantilever.py` | writes either deck: clamped root (`ROOT, ENCASTRE`), every tip node (`TIP`) pushed in -y, one static step with one fixed increment per load step (`*Static, direct`), NLGEOM=NO, every increment written to the ODB. The FCC deck also raises the iteration limits (`*Controls, parameters=time incrementation`), because that UMAT returns its elastic stiffness as the tangent |
| `j2_request.json` | six outputs at every increment: `tip_RF2` (sum over `TIP`), `tiptop_U1` (node 833), `mises_mean` (volume mean), `mises_root_max` (max over `ROOTEL`, the elements at the root), `S11_e1_ip1`, `eqplas_max` (max plastic strain); all parameters; weighted shares of the von Mises field; `"full_field": true` |
| `fcc_request.json` | four outputs: `tip_RF2`, `mises_mean`, `mises_root_mean` (volume mean over `ROOTEL`), `e1_ip1_S11`; all ten parameters; weighted shares; full field |
| `export_odb.py`, `run_fd.sh`, `fd_reference.py` | the ODB exporter for Abaqus Python, and the perturbed Abaqus reruns with their central differences (optional reference, see [README.md](README.md)) |

## Run it from the command line

From the Residual_Assembler root ([INSTALL.md](../../docs/INSTALL.md)):

```bash
export RA="$PWD"
export UMAT="$RA/../UMAT_source_transformation"
export WORK="$HOME/resasm_work"
C="$WORK/cantilevers"; mkdir -p "$C/j2" "$C/fcc"
```

**Step 1. Write the decks.**

```bash
python examples/cantilevers/gen_cantilever.py j2  --out "$C/j2/cantilever_j2_nominal.inp"
python examples/cantilevers/gen_cantilever.py fcc --out "$C/fcc/cantilever_fcc_nominal.inp"
```

Measured output:

```text
wrote .../j2/cantilever_j2_nominal.inp: 1536 C3D8, 2499 nodes, 7497 DOF, 12288 IPs, PROPS=[200000.0, 0.3, 250.0, 2000.0]
wrote .../fcc/cantilever_fcc_nominal.inp: 384 C3D8, 675 nodes, 2025 DOF, 3072 IPs, PROPS=[168000.0, 121000.0, 75000.0, 13.0, 55.0, 800.0, 2.0, 1.4, 0.001, 0.05]
```

**Step 2. Run the two analyses in Abaqus** (once; this is the ordinary
analysis, with the ORIGINAL UMATs):

```bash
cd "$C/j2" && abaqus job=cantilever_j2_nominal input=cantilever_j2_nominal.inp \
    user="$UMAT/parameter_sensitivity/models/m3_j2/umat.for" double=both interactive
cd "$C/fcc" && abaqus job=cantilever_fcc_nominal input=cantilever_fcc_nominal.inp \
    user="$UMAT/parameter_sensitivity/models/m6_fcc/umat.for" double=both interactive
cd "$RA"
grep "COMPLETED SUCCESSFULLY" "$C/j2/cantilever_j2_nominal.sta" "$C/fcc/cantilever_fcc_nominal.sta"
```

The Abaqus analyses were not re-run for this guide. The measured numbers below
come from Abaqus 2021.HF5 results of the decks that `gen_cantilever.py`
writes. Judge each job by the `.sta` file (see Common problems).

**Step 3. Build the two material providers.**

```bash
umat-oti-provider build "$UMAT/parameter_sensitivity/models/m3_j2/contract_v2.json"  --out "$WORK/provider_j2"
umat-oti-provider build "$UMAT/parameter_sensitivity/models/m6_fcc/contract_v2.json" --out "$WORK/provider_fcc"
```

The FCC mapping lists its parameters in the provider's order:
`g0, h0, q, gd0, m, gsat, C11, C12, C44, a`.

**Step 4. Compute the sensitivities from the ODBs.**

```bash
resasm request --model "$C/j2/cantilever_j2_nominal.inp" --odb "$C/j2/cantilever_j2_nominal.odb" \
    --material "$WORK/provider_j2/umat_m3_j2_oti.obj" \
    --request examples/cantilevers/j2_request.json --out "$C/j2_results"

resasm request --model "$C/fcc/cantilever_fcc_nominal.inp" --odb "$C/fcc/cantilever_fcc_nominal.odb" \
    --material "$WORK/provider_fcc/umat_m6_fcc_oti.obj" \
    --request examples/cantilevers/fcc_request.json --out "$C/fcc_results"
```

Each command first says that the model is outside the scope of the bounded
single-material engine and names the reason (here
`unsupported *Static options: ['direct']`), then uses the history engine and
ends with `request executed: history engine, 40 increments; verified=False`
(25 increments for FCC). The ODB export is kept in `private/fields.npz`.

**Step 5. Re-equilibrate every increment** (recommended for final numbers).
This reuses the export of step 4, so it needs no Abaqus:

```bash
resasm history --model "$C/j2/cantilever_j2_nominal.inp" --fields "$C/j2_results/private/fields.npz" \
    --material "$WORK/provider_j2/umat_m3_j2_oti.obj" \
    --request examples/cantilevers/j2_request.json --out "$C/j2_polished" --reequilibrate

resasm history --model "$C/fcc/cantilever_fcc_nominal.inp" --fields "$C/fcc_results/private/fields.npz" \
    --material "$WORK/provider_fcc/umat_m6_fcc_oti.obj" \
    --request examples/cantilevers/fcc_request.json --out "$C/fcc_polished" --reequilibrate
```

`--reequilibrate` Newton-polishes each recorded increment to double-precision
equilibrium before the sensitivities are taken. It matters when the analysis
was converged loosely: the recorded FCC state has a free residual of up to
1.5e-2 N, close to its limit of 3.4e-2 N (see the checks below).

The folder's own exporter gives the same arrays, and its file can be passed
to `resasm history --fields` directly:

```bash
cd "$C/j2" && abaqus python "$RA/examples/cantilevers/export_odb.py" -- \
    cantilever_j2_nominal.odb cantilever_j2_nominal_fields.npz && cd "$RA"
```

Measured: 15.3 s, and every array identical to the export kept by step 4;
the replay from it gave identical results.

## Run it in the GUI

Start `streamlit run scripts/app.py` ([GUI guide](../../docs/GUI_GUIDE.md)) and
stay on **Sensitivity Request**:

1. **OTI_UMAT.obj path**: `$WORK/provider_j2/umat_m3_j2_oti.obj` (its
   `umat_m3_j2_oti.json` is found beside it). **Analysis.inp path** and
   **Analysis.odb path**: the J2 deck and ODB.
2. **sensitivity_request.json path**: `examples/cantilevers/j2_request.json`.
   The screen's own output list offers U, RF, S and SDV1 only, so give the file
   to obtain the von Mises outputs, the shares and the full field. A caption
   confirms that the file is used instead of the controls.
3. **Output directory**: a new folder. Press **Solve**.

Measured on 2026-09-18: Solve took 24.5 s and reported
`Executed: 240 scalar results. Independent validation: not run.` with the three
downloads. The full-field table under the downloads is shown only for runs of
the bounded engine; for these models the full field is the `fields.npz` in the
output directory, read as shown below. Re-equilibration is a command-line
option (step 5).

## What it writes

| Path in the output directory | Contents | Share? |
| --- | --- | --- |
| `sensitivity_results.json` | request, scope, metadata (ODB parity, tolerances, timings, parameter values, input hashes), every scalar result with `derivatives` and `weighted` (`p dQ/dp`), and `weighted_shares` | yes |
| `sensitivity_tables.csv` | one row per output, increment and parameter: 960 rows for J2 (6 x 40 x 4), 1,000 for FCC (4 x 25 x 10) | yes |
| `sensitivity_shares.csv` | per increment and parameter: field share %, scalar share %, field weight, volume-mean von Mises: 160 rows (J2), 250 (FCC) | yes |
| `run_report.txt` | the verdicts, timings and tolerance formulas | yes |
| `fields.npz` | the full field: `U, dU, RF, dRF` (increment, node, component, parameter), `S, dS, SDV, dSDV, MISES, dMISES` (increment, element, point, component, parameter), `ip_volume`, `time`, `node_labels`, `elem_labels`, `parameters`, `parameter_values`. 134 MB (J2), 111 MB (FCC) | only if the model may be shared: it contains the whole solution |
| `private/fields.npz` | the ODB export (14 MB for J2) | no |
| `private/run_details.json`, `private/odb_export.log`, `private/link/` | per-increment residuals, Newton corrections and parity; the export log; the link library | no |

## Expected output (measured on 2026-09-18)

**Replay against the ODB** (`run_report.txt`, recorded state):

| | J2 | FCC |
| --- | --- | --- |
| `Equilibrium passed` | `max|R_free| = 1.386e-02 N at increment 10 (limit 2.572e-01 N there)` | `max|R_free| = 1.546e-02 N at increment 3 (limit 3.375e-02 N there)` |
| stress, reaction, state error / limit | 0.010, 0.006, 0.002 | 0.006, 0.006, 0.002 |
| after `--reequilibrate` | `max|R_free| = 5.935e-10 N at increment 8 (limit 2.202e-09 N there)` | `max|R_free| = 3.477e-11 N at increment 3 (limit 2.088e-10 N there)` |

**J2, last increment** (re-equilibrated):

| Output | Value | d/dE | d/dnu | d/dSIGY0 | d/dH |
| --- | --- | --- | --- | --- | --- |
| tip_RF2 (N) | -776.516 | -1.2998e-04 | -4.9836e-01 | -2.8416e+00 | -2.0057e-02 |
| tiptop_U1 (mm) | 0.137653 | -6.0714e-08 | -2.0956e-03 | 4.2067e-05 | 8.1302e-07 |
| mises_mean (MPa) | 137.074 | 3.3472e-05 | -5.2501e-01 | 4.9312e-01 | 3.5500e-03 |
| mises_root_max (MPa) | 283.316 | 6.8416e-05 | 1.5129e+00 | 9.6635e-01 | 1.4022e-02 |
| eqplas_max | 0.0189376 | 3.9200e-08 | 2.5172e-03 | -1.9704e-05 | -1.4570e-06 |

**J2 weighted shares of the von Mises field** (%):

| Increment | E | nu | SIGY0 | H |
| --- | --- | --- | --- | --- |
| 1 to 6 (elastic) | 98.21 | 1.79 | 0 | 0 |
| 7 (first yield) | 97.63 | 1.78 | 0.59 | 0.00 |
| 10 | 67.81 | 1.35 | 30.74 | 0.10 |
| 16 | 25.09 | 1.04 | 72.90 | 0.96 |
| 40 | 4.83 | 0.69 | 89.05 | 5.43 |

**FCC weighted shares of the von Mises field** (%, largest first):

| Increment | Shares |
| --- | --- |
| 1 (elastic) | C11 57.27, C12 40.06, C44 2.67, all others 0 |
| 3 (hardening begins) | C11 53.73, C12 37.57, g0 5.38, C44 2.82, gd0 0.27, ... |
| 10 | g0 39.95, C11 21.59, C12 15.15, C44 4.64, h0 4.61, q 3.99, gsat 3.17, a 2.72, gd0 2.37, m 1.81 |
| 25 | g0 34.56, C11 15.04, C12 10.64, h0 9.65, q 8.34, gsat 7.65, a 6.42, C44 3.29, gd0 2.54, m 1.86 |

**Full fields.** Read `fields.npz` at the last increment:

```bash
python - "$C/j2_polished/fields.npz" SIGY0 <<'EOF'
import sys, numpy as np
d = np.load(sys.argv[1])
names, p = [str(x) for x in d["parameters"]], d["parameter_values"]
n = len(d["time"]) - 1                                   # last increment
mises, dmises, V = d["MISES"][n], d["dMISES"][n], d["ip_volume"]
print("volume-mean von Mises:", (V * mises).sum() / V.sum())
weighted = np.abs(dmises * p)                            # |p dq/dp| at every point
governs = weighted.argmax(axis=-1)
print("governing parameter:", {names[j]: int((governs == j).sum()) for j in range(len(names))})
j = names.index(sys.argv[2])
e, q = np.unravel_index(np.abs(dmises[..., j]).argmax(), mises.shape)
print(f"largest |dMISES/d{names[j]}| = {abs(dmises[e, q, j]):.4e} at element {d['elem_labels'][e]}, point {q + 1}")
EOF
```

```text
volume-mean von Mises: 137.07355123165271
governing parameter: {'E': 92, 'nu': 4, 'SIGY0': 12192, 'H': 0}
largest |dMISES/dSIGY0| = 1.0002e+00 at element 783, point 1
```

For FCC (`"$C/fcc_polished/fields.npz" g0`): volume mean 24.675457368553978,
governing parameter `g0` at 2,456 of 3,072 points, `C11` at 592, `C44` at 16,
`h0` at 8; largest `|dMISES/dg0| = 3.6134e+00` at element 290, point 1. The
volume means agree with the scalar output `mises_mean` exactly, and their
derivatives to 4.7e-16 (relative).

## How the result is checked independently

1. **Replay against the ODB** at every integration point and increment (table
   above). Any excess over the single-precision limits would stop the run.

2. **Homogeneity identity**, for every output and increment:

   ```bash
   python - "$C/j2_polished/sensitivity_results.json" E,SIGY0,H <<'EOF'
   import json, sys
   rows = json.load(open(sys.argv[1]))["results"]
   stress_like = sys.argv[2].split(",")                     # the stress-dimensioned parameters
   degree = {"RF": 1, "S": 1, "MISES": 1, "U": 0, "SDV": 0}
   worst = 0.0
   for r in rows:
       w = [r["weighted"][p] for p in stress_like]          # weighted = p * dQ/dp
       scale = max(max(abs(x) for x in w), abs(r["value"]))
       if scale:
           worst = max(worst, abs(sum(w) - degree[r["field"]] * r["value"]) / scale)
   print(f"{len(rows)} outputs x increments, largest homogeneity residual {worst:.1e}")
   EOF
   ```

   For FCC, pass `"$C/fcc_polished/sensitivity_results.json" C11,C12,C44,g0,gsat,h0`.

   | Run | J2 (240 values) | FCC (100 values) |
   | --- | --- | --- |
   | recorded state | 8.4e-5 | 5.4e-4 |
   | re-equilibrated | 1.4e-12 | 1.2e-13 |

   On the recorded state the residual is the size of the ODB's own
   equilibrium error; after re-equilibration it is rounding. The weighted
   derivatives of the two runs differ by at most 8.9e-6 (J2) and 3.6e-3 (FCC)
   of the largest weighted derivative of the same output and increment, which
   is why step 5 is recommended for the loosely converged FCC analysis.

3. **Whole-model finite differences of the ORIGINAL UMAT at full size.** The
   script below re-equilibrates the whole J2 cantilever in Python with the
   unmodified UMAT at `p (1 +/- h)`, `h` = 1e-3, 1e-4, 1e-5, and compares every
   field at every increment. Measured for `SIGY0` (two runs, 217 s and 240 s, identical results):

   ```bash
   python scripts/replay_history_fd_fullsize.py --deck "$C/j2/cantilever_j2_nominal.inp" \
       --fields "$C/j2_results/private/fields.npz" --object "$WORK/provider_j2/umat_m3_j2_oti.obj" \
       --parameters SIGY0 --steps 1e-3,1e-4,1e-5 --rtol 1e-12 --out "$C/fd_SIGY0.json"
   ```

   Where the finite difference has a plateau (two adjacent steps agree to
   1e-5), the largest relative difference over each whole field was
   U 2.6e-9, RF 3.0e-9, S 3.1e-9, MISES 7.9e-9 and SDV 9.1e-11. At the other
   increments the finite difference itself changes with the step, and the OTI
   value agrees with its finest step to 1.7e-8 or better, except at increment
   11 (reaction, 1.4e-4 against a finite-difference spread of 8.8e-5) and
   increment 33 (S, SDV and MISES, up to 6.0e-2 against spreads up to 2.6e-2).
   At those two increments integration points sit on the yield surface, and
   `p (1 + h)` and `p (1 - h)` fall on different branches for every `h`: the
   finite difference, not the derivative, is unresolved there.

4. **Abaqus finite differences** (optional; 24 more Abaqus runs for J2 and
   60 for FCC): `run_fd.sh` and `fd_reference.py` in this folder
   ([README.md](README.md)). ODB-based differences are limited by single
   precision and by the Abaqus convergence tolerance. On the J2 tip reaction,
   at the increments where the estimates for `h` = 0.01 and 0.005 agree to
   1e-3, the OTI values matched the `h` = 0.005 estimate to within 1.3e-3 (H),
   1.7e-3 (SIGY0), 3.8e-3 (nu) and 1.5e-2 (E), relative. The whole-model check
   in point 3 is the precise reference.

## Run time

Measured on 2026-09-18 on a 24-core Linux workstation (Python 3.11), while
other jobs were running on it:

| Run | J2 | FCC |
| --- | --- | --- |
| provider build | 5.5 s | 8.6 s |
| `resasm request` from the ODB (includes the export) | 27.3 s wall, 11.7 s engine | 25.3 s wall, 18.0 s engine |
| `resasm history --fields`, recorded state | 23.9 s wall, 12.6 s engine | 25.3 s wall, 19.5 s engine |
| `resasm history --fields --reequilibrate` (two runs) | 36.8 and 39.8 s wall, 27.0 and 28.5 s engine | 71.9 and 79.4 s wall, 66.1 and 74.0 s engine |
| GUI Solve (request file) | 24.5 s | not run |
| whole-model FD, one parameter | 217 to 240 s | not run |

The difference between wall and engine time is mostly writing `fields.npz`
(over 100 MB). The J2 time is dominated by the sparse factorisations, the FCC
time by the material (a crystal-plasticity update at 3,072 points and ten
parameter directions per increment).

## Common problems

| Symptom | Cause and fix |
| --- | --- |
| The Abaqus job returns an error code (signal 6) after the analysis finished | Abaqus 2021.HF5 can abort during teardown on machines whose process IDs exceed 999,999, after writing a complete ODB. Judge the job by `THE ANALYSIS HAS COMPLETED SUCCESSFULLY` in the `.sta` file. |
| You write your own FCC deck | Keep the `*Controls, parameters=time incrementation` line that `gen_cantilever.py` writes: this UMAT returns its elastic stiffness as the tangent, so Abaqus needs more iterations per increment than its defaults allow. |
| `request failed` or `history replay failed` mentioning Abaqus Python | The ODB export needs `abaqus` on `PATH` (or `--abaqus`). Once exported, rerun with `resasm history --fields`. |
| Disk space | Each output directory with `"full_field": true` takes 110 to 135 MB; set `"full_field": false` if you need only the scalar outputs. |
| `Unsupported feature detected` | Only the deck subset in [docs/REPLAY_HISTORY.md](../../docs/REPLAY_HISTORY.md) is replayed; everything else is refused by name. |
| Results differ slightly between the recorded and the re-equilibrated runs | Expected when the ODB was converged loosely (see check 2). Quote the re-equilibrated numbers. |

Previous: [Example 4](../replay_history/WALKTHROUGH.md).
All examples: [examples/README.md](../README.md).
