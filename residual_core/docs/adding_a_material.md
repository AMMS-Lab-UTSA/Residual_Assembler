# Adding a Material Backend

A **material** answers the constitutive question at a point only: given a
kinematic input and the previous state, return stress, tangent, and updated
state. It knows nothing about elements, shape functions, or assembly — that is the
formulation's job. Adding a material means writing a backend that satisfies the
contract below and registering it.

> **A material is NOT always a UMAT.** It may be a continuum stress–strain law, a
> beam moment–curvature law, a cohesive traction–separation law, a thermal/
> diffusion flux law, or a Python return map. The contract is deliberately general
> so any of these present the same face to a formulation.

Worked references: `materials/elastic_adapter.py` (`IsotropicElastic`, runnable
offline), `materials/umat_adapter.py` (generic Abaqus-UMAT bridge),
`materials/crystal_plasticity_adapter.py` (**one example backend**, not the core).

---

## 1. Subclass `Material` and declare capabilities

```python
from .base import Material

class IsotropicElastic(Material):
    name              = "isotropic_elastic"
    n_state_vars      = 0
    stress_measure    = "cauchy"          # 'cauchy' | 'pk2' | 'pk1'
    tangent_measure   = "ddsdde"          # 'ddsdde' | 'material' | None
    kinematic_input   = "small_strain"    # 'small_strain' | 'deformation_gradient'

    # capability declaration (reported by the registry / inspector)
    constitutive_kind      = "stress_strain"
    input_variables        = ("strain", "dstrain")
    output_variables       = ("stress", "tangent")
    parameters             = ("E", "nu")               # PROPS names, if known
    supported_formulations = ("solid_c3d8_small_strain",)  # empty = any compatible
    limitations            = ()
    notes                  = ""
```

Declare honestly:
- **`constitutive_kind`** — what law family this is (not necessarily mechanics).
- **`input_variables` / `output_variables`** — what `evaluate` reads and returns.
- **`kinematic_input`** — `small_strain` needs `strain`, `dstrain`;
  `deformation_gradient` needs `F0`, `F1`.
- **`stress_measure` / `tangent_measure`** — so the formulation bridges correctly;
  `tangent_measure = None` means no tangent (the `tangent_available` property
  reflects this).
- **`n_state_vars`** — mirrors Abaqus NSTATV; the formulation supplies/receives
  the per-IP state slab.
- **`parameters`** — PROPS names.
- **`supported_formulations`** — leave empty for "any compatible".

---

## 2. Implement `evaluate`

```python
def evaluate(self, kinematics, state_prev, binding, time, dtime, fields, options):
    ...
    return (stress_voigt,   # (6,) Abaqus Voigt order (11,22,33,12,13,23), in stress_measure
            tangent,        # (6,6) or None
            state_new,      # (n_state_vars,)
            diagnostics)    # dict
```

- Read parameters (PROPS) from the `MaterialBinding` passed in `binding`.
- History-dependent laws must use `state_prev` and return `state_new`; the
  `StateManager` enforces committed/trial discipline across increments (never jump
  to the final step — see `limitations.md`).
- Override `init_state(binding, coords, n_ip)` to seed a non-zero initial state.

---

## 3. Register it

Add the backend in `materials/registry.py::build_material_registry()`. The
inspector will then recognize it, and `resasm template --material NAME` will print
its declared contract.

---

## 4. Materials vs. modes (which data is needed)

The requirements engine (see `minimal_input_contract.md`) reflects material
capability:
- A **stateless built-in law** (e.g. `IsotropicElastic`) needs no history and no
  state — material-replay is runnable offline with just PROPS.
- A **history-dependent UMAT** (crystal plasticity) needs the full increment
  sequence and, for genuine stress, the compiled Fortran (Intel `ifort` +
  Abaqus). Offline it runs only via a Python mock backend, or use the
  stress-driven mode with an exported field instead. The framework **reports**
  this rather than pretending.
