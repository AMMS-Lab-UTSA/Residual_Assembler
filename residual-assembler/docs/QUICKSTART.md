# Residual Assembler — quickstart

Get material-parameter sensitivities from a **saved** analysis, without re-running it.

## Install

```bash
pip install -e .            # console script: `resasm`
python scripts/run_tests.py --core   # 3-test offline smoke check (should pass)
```

## The residual request contract (recommended)

One high-level request pairs **each output with its own material parameters**, a
single `scope` says where and when, and the material is just the OTI object (its
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
python -m residual_core.replay.request_contract init     request.json
python -m residual_core.replay.request_contract validate request.json
python -m residual_core.replay.request_contract run      request.json
```

- Output names: `U` (displacement), `RF` (reaction), `S` (stress), or a
  state-variable name such as `EQPLAS`.
- `with_respect_to` is checked **against the material contract** up front; only
  the union of all requested parameters is seeded, and each output is reported
  against just its own parameters.
- `scope.domain` = `"ALL"` or a named node/element set in the record;
  `scope.increments` = `"ALL_CONVERGED"` or `"LAST"`.

## The lower-level sensitivity job (per-point outputs)

For fine-grained control (a specific node/component/reduction), the
`sensitivity_job.json` form drives per-point outputs. Start from the template,
then let the validator catch every mistake before anything expensive happens.

```bash
# 1. write a ready-to-edit starter job
python -m residual_core.replay.job init job.json

# 2. edit job.json: point material.oti_umat / material.contract at the package
#    you received, analysis.record at the converged run, and list the parameters
#    and outputs you care about.

# 3. validate — reports EVERY problem at once, with the fix
python -m residual_core.replay.job validate job.json
#   3 problem(s) in job.json:
#     - material.oti_umat not found: material/umat_oti.obj
#     - parameter 'YOung' is not in the material contract; available: E, nu
#     - outputs[1].component 'RF9' invalid for reaction_force; use RF1, RF2, RF3

# 4. run — replays the record once, assembles R / dR/dp / K, solves, writes results
python -m residual_core.replay.job run job.json
```

The job is validated up front (schema, required blocks, referenced-file
existence, parameter names checked **against the material contract before the
object is linked**, and each output's type / component / region / reduction).
You get one aggregated, actionable error — never a raw `KeyError` or a Fortran
compile dump.

### Job structure (`resasm_sensitivity_job_v1`)

```json
{
  "schema": "resasm_sensitivity_job_v1",
  "material": { "oti_umat": "material/umat_oti.obj", "contract": "material/umat_oti.json" },
  "analysis": { "model": "analysis/model.inp", "record": "analysis/record.odb" },
  "parameters": [ {"name": "E"}, {"name": "nu"} ],
  "outputs": [
    { "output_id": "tip_displacement", "type": "nodal_displacement",
      "region": {"node_set": "LOADED_FACE"}, "component": "U1", "reduction": "average" },
    { "output_id": "support_reaction", "type": "reaction_force",
      "region": {"node_set": "FIXED_FACE"}, "component": "RF1", "reduction": "sum" }
  ]
}
```

- `parameters` also accepts the shorthand `["E", "nu"]`.
- `reduction` defaults per output type — `average` for displacements, `sum` for
  reaction forces, `volume_average` for integration-point stress — so a face
  reaction with no explicit reduction is summed over the whole face (not just the
  first node).

## Interface versions

Every schema tag and the C-ABI/material-package contract hash live in one place:

```bash
python -m residual_core.interface.versions
```

Pin these when you share a package so a stale binary or a drifted contract is
caught immediately.

## Tests

```bash
python scripts/run_tests.py            # full suite; OTILib/Abaqus tests -> SKIP offline
python scripts/run_tests.py --core     # fast offline smoke set
pytest -q                              # pytest view (otilib/abaqus markers opt-in
                                       #   via RESASM_RUN_OTILIB=1 / RESASM_RUN_ABAQUS=1)
```
