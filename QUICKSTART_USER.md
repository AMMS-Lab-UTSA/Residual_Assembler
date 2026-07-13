# Quickstart (users)

Compute residual-method sensitivities from **one residual + one config + one
command**. This uses the Python path.

## 1. Install OTILib (once)

OTILib is the sensitivity engine (GPLv3, source build):
<https://github.com/mauriaristi/otilib.git>. See
[residual_core/docs/otilib_integration.md](residual_core/docs/otilib_integration.md)
or run `scripts/setup_otilib.sh`. **Do not** `pip install pyoti` (unrelated).
The black-box path does not require OTILib on our side.

## 2. Start from a template

```
pip install -e .                          # exposes the `resasm` command
resasm init --template python --out my_job && cd my_job
```

`--template` is one of `python`, `blackbox`, `cpp`, `fortran`. It just copies the
matching folder from `templates/` — no prompts. (`resasm init` with no `--template`
still runs the interactive wizard.)

## 3. Edit the residual — `user_residual.py`

```python
def residual(u, params, state=None, time=None):
    k = params["k"]; f = params["f"]
    return [k * u[0] ** 3 - f]        # your model here

def tangent(u, params, state=None, time=None):
    return [[3.0 * params["k"] * u[0] ** 2]]
```

## 4. Edit `resasm.yml`

```yaml
problem:   { name: my_job }              # `unknowns` is optional -> inferred from solution.npy
residual:  { type: python, module: user_residual.py, function: residual }
tangent:   { type: python, function: tangent }
parameters: { k: 2.0, f: 16.0 }
solution:  { file: solution.npy }
sensitivity: { order: 2, backend: otilib }
validation: { rhs_finite_difference_check: true }   # checks d(residual)/d(param), not a re-solve
```

## 5. Check it

```
resasm check resasm.yml
```

You should see `[ok]` lines ending in `sensitivity solve completed`. If something
is missing, the first `[fail]` tells you exactly what to add.

## 6. Run it

```
resasm run resasm.yml
```

## 7. Read the result

```
resasm report resasm_output/
open resasm_output/public/summary.md
```

`public/` is safe to share (norms + rankings only). `private/` holds the full
`R^(p)` / `U^(p)` arrays and stays on your machine.

Can't use Python? See [docs/which_path_should_i_use.md](docs/which_path_should_i_use.md)
for the C++/Fortran and black-box paths.
