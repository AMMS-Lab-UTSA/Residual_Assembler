# residual-assembler — offline material-parameter sensitivities (Program 2)

Get the sensitivity of any FE response to material parameters — `d(response)/d(p)`
— from a **saved** analysis, **without re-running** the production job. It
consumes the OTI-enabled object produced by
[`umat-oti`](https://…/umat-oti) (Program 1), replays the material at each integration
point, assembles the residual and its derivative, and solves the sensitivity
system:

$$ K\,\frac{\partial u}{\partial p} = -\frac{\partial R}{\partial p},\qquad
   \frac{\partial R}{\partial p} = \int_V B^{\mathsf T}\,\underbrace{\frac{\partial\sigma}{\partial p}}_{\texttt{DSIGMA\_DP}}\,dV - \frac{\partial F_\text{ext}}{\partial p} $$

The production analysis is **replayed once**, never re-solved — so all parameter
sensitivities cost roughly one extra analysis instead of `2N+1` finite-difference
re-runs, and the material source stays private (only the compiled object crosses
the boundary).

```
   OTI object + contract (from umat-oti)  +  converged FE record  ──►  residual-assembler  ──►  d(response)/d(parameter)
```

## Install

```bash
pip install -e .            # requires numpy (PyYAML optional). Console script: `resasm`
python scripts/run_tests.py --core   # 3-test offline smoke check (should pass)
```

Linking and loading the distributed material object needs `gfortran` on `PATH`
(or `FC`). Program 2 loads that object through `ctypes` the same way Program 1
does, and works natively on **Windows and Linux** via
[`residual_core.runtime.libload`](residual_core/runtime/libload.py): the linked
library gets the platform's suffix, and on Windows the gfortran runtime directory
is discovered and registered with `os.add_dll_directory`. Both platforms run the
test suite in [CI](.github/workflows/ci.yml).

## Quickstart — a residual request

One high-level request pairs each **output** with its **own** material parameters;
a single `scope` selects where and when; the material is the OTI object (its
completed contract is the sibling `.json`):

```json
{
  "material": "umat_j2_oti.obj",
  "model":    "j2_tension.inp",
  "record":   "j2_tension.resrec.h5",
  "scope":    {"domain": "ALL", "increments": "ALL_CONVERGED"},
  "requests": [
    {"output": "U",      "with_respect_to": ["E", "nu"]},
    {"output": "S",      "with_respect_to": ["E", "SIGY0", "H"]},
    {"output": "EQPLAS", "with_respect_to": ["SIGY0", "H"]}
  ]
}
```

```bash
python -m residual_core.replay.request_contract init     request.json   # starter
python -m residual_core.replay.request_contract validate request.json   # actionable checks
python -m residual_core.replay.request_contract run      request.json   # replay + assemble + solve
```

Validation happens up front — schema, referenced-file existence, output names,
and each request's parameters checked **against the material contract before the
object is linked** — so you get one aggregated, plain-language error instead of a
raw traceback or a Fortran dump.

For fine-grained per-point outputs (a specific node / component / reduction) the
lower-level `sensitivity_job.json` form is also available
(`python -m residual_core.replay.job init/validate/run`). See
[docs/QUICKSTART.md](docs/QUICKSTART.md).

## Validation

The whole framework is validated across **20 material models** (elasticity,
plasticity, viscoplasticity, viscoelasticity, crystal plasticity) — the OTI
`DSIGMA_DP` against finite differences of the original UMAT (Program 1) and the
residual-method `d(σ_vM)/dp` against full-analysis FD (Program 2). See
[`results/`](results): `per_umat_report.html` (per-material contracts, plots and
error tables), `report.html`, and the sweep tables. The residual assembler is
also validated across six element formulations (C3D8/R, C3D20/R, C3D4, C3D10).

## What's here

```
residual_core/           the framework: replay, assembly, formulations, engine, interface
  replay/                request_contract, job runner, engine, path_material, contracts
  formulations/          C3D8 + solid3d element kernels
  interface/versions.py  single source of truth for every schema tag + the ABI hash
resasm_user/  partner_kit/   user-facing layers
schemas/                 JSON schemas (job, request, material package, record, …)
tests/                   framework test suite    scripts/run_tests.py   the runner
results/                 validation scripts, figures and HTML reports
docs/                    QUICKSTART and contract references
```

See [INTERFACE.md](INTERFACE.md) for the exact hand-off from `umat-oti`.

## License

MIT (matches OTILib / pyoti) — see [LICENSE](LICENSE).
