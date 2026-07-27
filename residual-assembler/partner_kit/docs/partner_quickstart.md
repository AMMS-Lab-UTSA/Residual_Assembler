# Partner Quickstart

Run the residual-sensitivity kit locally in three steps. Requires Python 3 +
NumPy. Nothing you provide leaves your machine unless you choose to share the
public report.

## 0. Get the kit on your path

```bash
export PYTHONPATH=/path/to/partner_kit/..:$PYTHONPATH   # so `import partner_kit` works
# (Windows PowerShell)  $env:PYTHONPATH = "C:\path\to\partner_kit\.."
```

## 1. Implement a residual provider (pick one level)

A Python class in your own file (`provider.py`). It must be **scalar-generic**:
parameter values may arrive as floats or as a hypercomplex `Dual1`.

```python
import numpy as np
from partner_kit.python.residual_provider import GlobalResidualProvider

class MyModel(GlobalResidualProvider):
    name = "my_model"
    parameters = ("k",)          # design parameters you want sensitivities for
    ndof = 1

    def parameter_values(self):
        return {"k": 2.0}

    def eval_global_residual(self, u, state, parameters, time, dtime):
        k = parameters["k"]                     # float OR Dual1
        return [k * u[0] ** 3 - 16.0]           # ordinary arithmetic

    def get_tangent(self, u, state, parameters, time, dtime):
        return np.array([[3.0 * parameters["k"] * u[0] ** 2]])
```

Element-level and black-box options: see
[residual_provider_contract.md](residual_provider_contract.md).

## 2. Write a config

`config.json` next to your provider:

```json
{
  "provider_file": "provider.py",
  "provider_class": "MyModel",
  "parameters": ["k"],
  "order": 1,
  "solve": true,
  "u0": [1.0],
  "output_dir": "out"
}
```

Supply `"solution": [...]` (or a `.npy` path) instead of `"solve": true` if you
already have a converged `u`.

## 3. Run

```bash
python -m partner_kit.python.partner_cli --config config.json
```

Output:

```
provider        : my_model (MyModel)
parameters      : ['k']
tangent source  : provider-get_tangent
R^(1) shape     : (1, 1)
hypercomplex ok : True
  dU/d(k         ) = -3.333333e-01
validation      : True
outputs         : out/{sensitivity_package,public_report}
```

- `out/sensitivity_package/` — **private**, full arrays (keep local).
- `out/public_report/` — **shareable**, norms + validation summary only.

## What to send back (optional)

Only if you want to: the contents of `public_report/`
(`validation_summary.md`, `parameter_ranking.csv`, `sensitivity_norms.csv`,
`timing_summary.json`, `errors.json`). See [privacy_model.md](privacy_model.md).

## Black-box (no linking) in one line

```json
{ "blackbox": { "command": ["./your_solver"], "ndof": 1, "parameters": ["k"],
                "parameter_values": {"k": 2.0} },
  "parameters": ["k"], "order": 1, "solution": [2.0], "output_dir": "out" }
```

Your executable reads `request.json` and writes `response.json` — see
[residual_provider_contract.md](residual_provider_contract.md) and
[../templates/blackbox_executable](../templates/blackbox_executable).
