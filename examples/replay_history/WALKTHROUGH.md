# Example 4: history replay of a J2 beam, offline

A small Abaqus analysis whose output is already exported and committed, so it
runs anywhere in about a second, without Abaqus. A C3D8 beam made of a
J2-plastic material is bent into the plastic range over ten increments.
`resasm history` replays the whole history, returns the sensitivities of
reactions, displacements, stresses and von Mises stress at every increment,
and ranks the parameters by their weighted influence. Three independent checks
come with it.

| | |
| --- | --- |
| Folder | `examples/replay_history/j2_beam/` |
| Needs Abaqus | No (the ODB export and the Abaqus finite differences are committed) |
| Needs OTILib | No |
| Run time | about 1 s (12 s with the whole-model finite-difference check) |

## What it demonstrates

- The **history engine** behind `resasm history`. It accepts any material
  provider built by UMAT-OTI, prescribed displacements, many increments,
  node and element sets and von Mises outputs. `resasm request` hands every
  model of this kind to it automatically.
- **Total-history sensitivities.** In a plastic material the stress at the
  end depends on everything that happened before. The derivatives are carried
  from increment to increment, not recomputed from the last state alone.
- **Weighted shares**, which make parameters with different units comparable
  and show how the governing parameter changes from elastic to plastic
  loading.
- **Re-equilibration** (`--reequilibrate`) of the recorded single-precision
  states to double precision.

## The mathematics

At increment `n`, with the recorded displacements `u_n` and parameters `p`:

    R_n(u_n, p) = sum_e sum_q w_q B_q^T sigma_{n,q} - F(t_n),
    sigma_{n,q} = UMAT(sigma_{n-1,q}, state_{n-1,q}, B_q u_{n-1}, B_q (u_n - u_{n-1}), p).

On the free degrees of freedom `R_n = 0`. The prescribed tip displacement does
not depend on `p`. Differentiating gives one linear solve per increment,

    K_ff du_f/dp = - dR_f/dp,       K = sum w B^T DDSDDE B,
    dsigma_n/dp  = dsigma^A/dp + DDSDDE B du_n/dp,
    dRF/dp       = dR_c/dp + K_cf du_f/dp,

where `dsigma^A/dp` comes from one call of the compiled material per point,
seeded with the derivatives of the previous increment. For the von Mises stress
`q = sqrt(3/2 s:s)` the chain rule gives `dq/dp = (dq/dsigma) : dsigma/dp`.

**Weighted share** of parameter `j` at increment `n`, for the von Mises field
over the element volume `V`:

    W_j(n) = sum_i V_i |p_j dq_i/dp_j| / sum_i V_i,        share_j = 100 W_j / sum_k W_k.

**Homogeneity identity.** J2 plasticity with linear hardening is homogeneous
of degree one in `(E, SIGY0, H)` at fixed `nu`: scaling the three by the same
factor scales every stress and leaves every strain unchanged. Under prescribed
displacements, Euler's theorem then gives, at every increment,

    E dQ/dE + SIGY0 dQ/dSIGY0 + H dQ/dH = Q     for reactions, stresses, von Mises,
    E dQ/dE + SIGY0 dQ/dSIGY0 + H dQ/dH = 0     for displacements and plastic strain.

The engine does not use this identity, so it is an independent check.

## Inputs

| File in `j2_beam/` | What it is |
| --- | --- |
| `Analysis.inp` | 12 x 4 x 2 unit C3D8 (96 elements, 195 nodes, 585 DOF, 768 integration points). Root face `x = 0` clamped (`ROOT, ENCASTRE`). Every node of the tip face (`TIP`) pushed to `U2 = -0.08` mm in one static step of 10 fixed increments. User material with `E = 200000`, `nu = 0.3`, `SIGY0 = 250`, `H = 2000` (MPa), one state variable. Tight Abaqus convergence controls |
| `fields.npz` | the Abaqus 2021.HF5 ODB exported frame by frame (U, RF, S, SDV in single precision, frame 0 included) by `residual_core/replay/odb_export_npz.py` |
| `sensitivity_request.json` | six outputs at every increment: `tip_RF2` (sum over `TIP`), `midtop_U2` (node 59), `tiptop_U1` (node 65), `e1_ip1_S11` (element 1, point 1), `mises_mean` (volume mean over all elements), `mises_root_max` (max over the element set `ROOTEL`); all parameters; `weighted_shares` on the von Mises field |
| `abaqus_fd.json` | central differences of perturbed Abaqus reruns, `(q(p(1+h)) - q(p(1-h))) / (2 h p)` for `h` = 0.01, 0.005, 0.002, 0.001, with the nominal values |

The material provider is built from the companion repository, as in step 1.
`make_beam_deck.py` and `regenerate_example.py` in this folder rebuild the deck
and all committed data with Abaqus; you do not need them to run the example.

## Run it from the command line

From the Residual_Assembler root ([INSTALL.md](../../docs/INSTALL.md)):

```bash
export RA="$PWD"
export UMAT="$RA/../UMAT_source_transformation"
export WORK="$HOME/resasm_work"; mkdir -p "$WORK"
B=examples/replay_history/j2_beam

# 1. the compiled J2 material (once)
umat-oti-provider build "$UMAT/parameter_sensitivity/models/m3_j2/contract_v2.json" --out "$WORK/provider_j2"

# 2. replay the recorded history
resasm history --model "$B/Analysis.inp" --fields "$B/fields.npz" \
    --material "$WORK/provider_j2/umat_m3_j2_oti.obj" \
    --request "$B/sensitivity_request.json" --out "$WORK/beam"

# 3. the same, with every increment re-equilibrated to double precision first
resasm history --model "$B/Analysis.inp" --fields "$B/fields.npz" \
    --material "$WORK/provider_j2/umat_m3_j2_oti.obj" \
    --request "$B/sensitivity_request.json" --out "$WORK/beam_polished" --reequilibrate
```

Measured output of step 2 (exit 0):

```text
history replay: 10 increments, 768 integration points, 4 parameters; yes: max|R_free| = 1.288e-03 N at increment 10 (limit 1.609e-01 N there)
```

and of step 3:

```text
history replay: 10 increments, 768 integration points, 4 parameters; yes: max|R_free| = 3.503e-12 N at increment 9 (limit 2.049e-09 N there)
```

The mapping is found automatically beside the object
(`umat_m3_j2_oti.json`). `--fields` replaces `--odb`: when you have an ODB
instead, give `--odb Analysis.odb` and the command exports it through Abaqus
Python first.

## Run it in the GUI

The **Sensitivity Request** screen reads an `.odb`; it has no field for an
existing `fields.npz` export, so this offline example is a command-line
example. With an ODB, the same request runs from the GUI: give the object, the
deck, the ODB and this `sensitivity_request.json` (the request file is used as
written), then press **Solve**. `resasm request` passes the model to the
history engine, and the three public files appear as downloads. See
[Example 5](../cantilevers/WALKTHROUGH.md) for a measured GUI run of this kind.

## What it writes

| Path in `$WORK/beam/` | Contents | Share? |
| --- | --- | --- |
| `sensitivity_results.json` | request, scope (elements, points, DOF, increments, parameters, provider), metadata (ODB parity numbers, tolerances, timings, parameter values, SHA-256 of every input) and 60 results, each with `value`, `derivatives` and `weighted` (`p dQ/dp`); plus `weighted_shares` | yes |
| `sensitivity_tables.csv` | 240 rows: output, field, component, reduction, domain, increment, time, value, parameter, derivative, weighted derivative | yes |
| `sensitivity_shares.csv` | per increment and parameter: field share %, scalar share %, field weight, volume mean of von Mises | yes |
| `run_report.txt` | the verdicts listed below, timings and the tolerance formulas | yes |
| `private/run_details.json` | per-increment residuals, parity, Newton corrections and sensitivity-solve residuals | no |
| `private/link/` | the generated link library for the object | no |

Add `"full_field": true` to the request to also write `fields.npz` with every
field and its derivatives at every increment ([Example 5](../cantilevers/WALKTHROUGH.md)).

## Expected output (measured on 2026-09-18)

`run_report.txt` of step 2, first lines:

```text
Status: executed successfully
Residual assembled: yes: C3D8 selective-reduced (B-bar), 96 elements, 768 integration points, 585 DOF, 10 increments, sparse assembly
Equilibrium checked: yes: free-DOF residual of the recorded state at all 10 increments
Equilibrium passed: yes: max|R_free| = 1.288e-03 N at increment 10 (limit 1.609e-01 N there)
Tangent available: yes: DDSDDE = dSTRESS/dDSTRAN from the provider's OTI strain directions
Tangent verified: not run
Derivative calculated: yes: total-history du/dp, dRF/dp, dS/dp, dSDV/dp, dMISES/dp for 4 parameters ['E', 'nu', 'SIGY0', 'H'] at 10 increments
Derivative verified: not run
Abaqus comparison available: yes (primal): replayed S, SDV and RF reproduce the ODB at every integration point and increment; largest error/limit ratios: stress 0.005 (max |dS| 7.752e-04 MPa), reaction 0.005 (max |dRF| 5.666e-04 N), state 0.001 (max |dSDV| 2.366e-10); no Abaqus derivative reference is part of this request
Unsupported feature detected: none
```

Tip reaction at the last increment (`tip_RF2 = -191.72617449210222` N):

| Parameter | dRF2/dp | p dRF2/dp (weighted) |
| --- | --- | --- |
| E | -1.7653752262283606e-04 | -35.31 |
| nu | -1.2378979584158118 | -0.37 |
| SIGY0 | -0.6189102241416395 | -154.73 |
| H | -8.459844826609506e-04 | -1.69 |

Weighted shares of the von Mises field (`sensitivity_shares.csv`, field share, %):

| Increment | E | nu | SIGY0 | H |
| --- | --- | --- | --- | --- |
| 1 to 4 (elastic) | 93.99 | 6.01 | 0 | 0 |
| 5 (first yield) | 85.86 | 5.56 | 8.58 | 0.002 |
| 7 | 43.81 | 3.70 | 52.27 | 0.21 |
| 10 | 24.72 | 3.03 | 71.44 | 0.81 |

The elastic increments depend on `E` and `nu` only. From increment 5 the
yield stress takes over, as plasticity spreads from the root.

## How the result is checked independently

1. **Replay against the ODB.** At every integration point and increment the
   replayed stress, state and reactions are compared with the recorded ones,
   with limits set by single-precision storage (see
   [docs/REPLAY_HISTORY.md](../../docs/REPLAY_HISTORY.md)). Any excess stops
   the run. Measured: all within 0.5 % of their limits (report above).

2. **Abaqus finite differences** (`abaqus_fd.json`). For every one of the
   first four outputs, parameters and increments (160 comparisons), the
   difference between the OTI derivative and the Abaqus central difference is
   within the finite difference's own uncertainty: the spread between its two
   closest step sizes plus the single-precision floor `eps32 |q| / (2 h p)`.
   Measured: 115 nonzero references and 45 zero references, no failure. At
   the last increment, for the tip reaction:

   ```bash
   python - "$WORK/beam_polished/sensitivity_results.json" "$B/abaqus_fd.json" <<'EOF'
   import json, sys
   ours = {(r["output"], r["increment"]): r["derivatives"] for r in json.load(open(sys.argv[1]))["results"]}
   fd = json.load(open(sys.argv[2]))
   for name in fd["parameters"]:
       estimates = [fd["fd"][name][str(h)]["tip_RF2"][10] for h in fd["steps"]]
       print(f"{name:6s} OTI {ours[('tip_RF2', 10)][name]: .6e}   Abaqus FD, h = {fd['steps']}: "
             + ", ".join(f"{e: .6e}" for e in estimates))
   EOF
   ```

   ```text
   E      OTI -1.765375e-04   Abaqus FD, h = [0.01, 0.005, 0.002, 0.001]: -1.765490e-04, -1.765370e-04, -1.765454e-04, -1.765347e-04
   nu     OTI -1.237895e+00   Abaqus FD, h = [0.01, 0.005, 0.002, 0.001]: -1.237392e+00, -1.239141e+00, -1.239777e+00, -1.230240e+00
   SIGY0  OTI -6.189102e-01   Abaqus FD, h = [0.01, 0.005, 0.002, 0.001]: -6.189026e-01, -6.189083e-01, -6.189098e-01, -6.189022e-01
   H      OTI -8.459845e-04   Abaqus FD, h = [0.01, 0.005, 0.002, 0.001]: -8.458138e-04, -8.460045e-04, -8.459091e-04, -8.459091e-04
   ```

   The scatter of the Abaqus values across step sizes is the resolution of an
   ODB-based finite difference; the OTI value sits inside it.
   `tests/replay_history/test_history_example.py` asserts all 160 comparisons.

3. **Homogeneity identity.** This snippet evaluates it for every output and
   increment:

   ```bash
   python - "$WORK/beam_polished/sensitivity_results.json" <<'EOF'
   import json, sys
   rows = json.load(open(sys.argv[1]))["results"]
   degree = {"RF": 1, "S": 1, "MISES": 1, "U": 0, "SDV": 0}
   worst = {}
   for r in rows:
       w = [r["weighted"][p] for p in ("E", "SIGY0", "H")]      # weighted = p * dQ/dp
       scale = max(max(abs(x) for x in w), abs(r["value"]))
       if scale:
           err = abs(sum(w) - degree[r["field"]] * r["value"]) / scale
           worst[r["field"]] = max(worst.get(r["field"], 0.0), err)
   print({k: f"{v:.1e}" for k, v in worst.items()})
   EOF
   ```

   Largest residual, relative to the largest term, over all 10 increments:

   | Run | RF | U | S | MISES |
   | --- | --- | --- | --- | --- |
   | recorded state (`$WORK/beam`) | 8.0e-6 | 6.8e-8 | 4.7e-8 | 2.1e-7 |
   | re-equilibrated (`$WORK/beam_polished`) | 9.4e-15 | 1.1e-15 | 1.6e-15 | 1.5e-15 |

   On the recorded state the residual has the size of the ODB's own
   equilibrium error. After re-equilibration it is rounding. Re-equilibration
   changed the weighted derivatives by at most 4.1e-5 of the largest weighted
   derivative of the same output (`midtop_U2`, increment 3, `nu`): this ODB was
   converged tightly.

4. **Whole-model finite differences of the ORIGINAL routine.** Add
   `--verify fd` to step 3 (measured 11.6 s). The material object also contains
   the unmodified UMAT; the command re-solves the whole model with it at
   `p (1 +/- h)` for a ladder of steps and compares. Measured report lines:
   `Tangent verified: yes: max relative error 1.88e-10 vs central FD of the
   ORIGINAL UMAT at 36 points` and `Derivative verified: yes: ... worst
   nonzero-derivative error 1.46e-07 (plateau spread 3.05e-07); zero
   references: |OTI - FD| <= 1.6e-09 on the field scale`.

## Run time

Measured on 2026-09-18 (24-core Linux workstation, Python 3.11): provider
build 5.5 s; step 2 0.96 s wall time (0.12 s in the engine); step 3 1.0 s
(0.20 s); with `--verify fd` 11.6 s. The same replay from a fresh installation
of both packages gave identical tables.

## Common problems

| Symptom | Cause and fix |
| --- | --- |
| `history replay failed: output directory must be empty or new` | Choose a new `--out`. |
| `Abaqus launcher 'abaqus' not found: exporting Analysis.odb needs Abaqus Python (or pass --fields with an existing export)` | You gave `--odb` on a machine without Abaqus. Export once where Abaqus is available (`abaqus python residual_core/replay/odb_export_npz.py -- Analysis.odb fields.npz`) and use `--fields`. |
| `Unsupported feature detected: ...` in `run_report.txt` | The deck uses something the history engine refuses (several steps, distributed loads, NLGEOM, other element types, amplitudes, ...). The report names it. See the supported subset in [docs/REPLAY_HISTORY.md](../../docs/REPLAY_HISTORY.md). |
| An object built by an older UMAT-OTI is refused with a rebuild message | Rebuild the provider with the current `umat-oti-provider build`. |
| `max is attained at N locations whose derivatives differ` | A `max`/`min` reduction whose extreme sits at several points with different derivatives is not differentiable. Use `volume_mean`, a smaller domain, or a `component` output. |

Previous: [Example 3](../presentation_request/WALKTHROUGH.md).
Next: [Example 5](../cantilevers/WALKTHROUGH.md), the same engine at full size.
All examples: [examples/README.md](../README.md).
