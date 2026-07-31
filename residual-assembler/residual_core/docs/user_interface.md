# User Interface

Two front ends, both thin orchestration layers over the agnostic core and the
backend registries. Neither contains physics.

- **Python API** — `residual_core.ResidualProblem` (progressive disclosure).
- **Command line** — `resasm` (`python -m residual_core.ui.cli`).

Optional **config** files (`ui/config.py`) exist only for overrides that
auto-detection cannot infer.

---

## Python API — `ResidualProblem`

The common path needs three lines:

```python
from residual_core import ResidualProblem

p = ResidualProblem.from_abaqus("model.inp")
p.inspect()                              # elements, materials, possible modes, missing data
R = p.assemble(mode="stress-driven")     # after a field export is attached
```

### Construction
| Call | Purpose |
|---|---|
| `ResidualProblem.from_abaqus(path)` | parse an Abaqus `.inp` into a neutral model |
| `ResidualProblem.from_neutral(path)` | load a solver-neutral JSON model |
| `p.save_neutral(path)` | write the neutral JSON model |

### Attaching data (only what a mode needs)
| Call | For mode | Purpose |
|---|---|---|
| `p.attach_results(dict_or_json)` | stress-driven | exported element field `{eid: (n_ip,6)}` or `{'stress_ip': {...}}` |
| `p.attach_subroutine(path)` | material-replay / direct | UMAT/UEL source file |
| `p.set_solution(U)` | any | supply the DOF solution vector |

### Inspect / requirements / assemble
| Call | Returns |
|---|---|
| `p.inspect()` | `InspectionReport` (elements, materials, modes, missing data) |
| `p.requirements(mode)` | `RequirementsReport` (have / missing / minimum next) |
| `p.assemble(mode=..., U=None, compute_tangent=False)` | global residual `R` (and `K` if requested) |

`assemble` chooses the formulation policy appropriate to the mode by querying the
registry (a backend is selected only if its `spec.supported_modes` includes the
requested mode). If a required input is missing it raises a clear error naming the
**minimum** missing item — it never silently produces a wrong result.

---

## Command line — `resasm`

```
resasm inspect       MODEL [--detail]
resasm requirements  MODEL --mode M [--fields F] [--subroutine S]
resasm assemble      MODEL --mode M [--fields F] [--subroutine S] [--tangent] [--out R.npy]
resasm verify        MODEL --fields F
resasm doctor        MODEL [--write-config-template PATH]
resasm template      (--formulation NAME | --material NAME)
resasm backends
resasm modes
```

- `inspect` — the model inspection summary (auto-detected element/material
  support and which modes are reachable). `--detail` adds a per-element-type
  backend-selection card (selected backend, verification status, limitations,
  available modes, and the minimum next step for anything unsupported).
- `requirements` — what a mode needs vs. what is present, printing only the
  **minimum** missing item (never a generic checklist).
- `assemble` — assemble the global residual for a mode. `--fields` (alias
  `--odb`) attaches an exported field export (JSON) for stress-driven mode. If a
  required input is missing it prints the clean requirements report and exits
  non-zero rather than crashing.
- `verify` — stress-driven residual + reaction check (equilibrium sanity).
- `doctor` — per-mode readiness for the whole model; optionally emits a
  pre-filled config template.
- `template` — print one backend's full declared contract (element/DOF types,
  supported modes, per-mode required/optional inputs, material-interface need,
  state requirements, tangent support, verification status, limitations).
- `backends` — a user-readable registry audit of **every** registered backend.
- `modes` — list assembly modes and their minimum inputs.

Run with `python -m residual_core.ui.cli ...` if no console script is installed.

### End-to-end walkthroughs

Four runnable, self-contained examples ship under `residual_core/examples/`, each
with a model file, the minimal command, the expected output, and a note of what
was auto-detected vs. supplied manually:

| Example | Command | Mode |
|---|---|---|
| `minimal_truss/` | `resasm assemble minimal_truss/model.json --mode formulation` | formulation |
| `minimal_beam/` | `resasm assemble minimal_beam/model.json --mode formulation` | formulation |
| `minimal_mixed/` | `resasm inspect minimal_mixed/model.json --detail` | formulation |
| `minimal_c3d8_stress_driven/` | `resasm assemble … --mode stress-driven --fields fields.json` | stress-driven |

Regenerate their model/field files with
`python -m residual_core.examples.generate_minimal`.

---

## Config (advanced, optional)

Auto-detection covers the common case. A config (`.yml` or `.json`) is only
needed to override: `mode`, `odb`, `subroutine`, `formulation_policy`
(`{ELEMENT_TYPE: backend}`), `material_backend`, `material_parameters`. Generate a
pre-filled starter with:

```
resasm doctor model.inp --write-config-template config.yml
```

YAML is used if PyYAML is installed; otherwise JSON (identical schema). See
`ui/config.py`.

---

## Runnable examples

`python -m residual_core.ui.examples` builds small models directly (no external
solver) and checks known results — a truss axial force `EA/L`, a cantilever tip
deflection `PL³/3EI`, and a **mixed truss + beam** model proving the physics-blind
core dispatches heterogeneous-DOF elements. These double as living documentation.
