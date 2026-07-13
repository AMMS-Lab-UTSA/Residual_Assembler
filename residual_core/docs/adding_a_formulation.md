# Adding a Formulation Backend

A **formulation** is an element weak form. The core assembler dispatches each
element to the formulation bound to it and scatters the returned element residual
(and tangent) into the global system. Adding support for a new element/physics
means writing a backend that satisfies the contract below and registering it — no
core change.

> The framework is formulation-agnostic *by architecture*. A formulation becomes
> **supported** once a backend satisfying this contract is registered **and
> verified**.

Worked reference backends: `formulations/truss2.py` (3 DOF/node, translational),
`formulations/beam2.py` (6 DOF/node, rotations), `formulations/shell_placeholder.py`
(contract only, not yet runnable).

---

## 1. Subclass `Formulation` and declare metadata

```python
from .base import Formulation

class Truss2(Formulation):
    name            = "truss2"                       # registry key (unique)
    element_types   = ("T3D2", "T2D2", "TRUSS2", "BAR2")
    dof_types       = ("UX", "UY", "UZ")             # per-node DOFs, in order
    sign_convention = "residual"                     # returns +F_internal
    verification_levels = (0, 1, 4)

    # backend spec metadata (used by the registry, inspector and CLI `template`)
    supported_modes    = ("formulation",)            # which assembly modes apply
    required_inputs    = ("coords", "connectivity", "dofs", "section: E, A")
    optional_inputs    = ()
    limitations        = ("small-strain linear-elastic axial bar",
                          "2 nodes, axial stiffness only, no bending")
    verification_tests = ("zero-field", "axial-force", "FD-tangent")
    notes              = "Proof backend: shows the core is not C3D8/CP-specific."
```

The metadata becomes a `BackendSpec` (via the `spec` property on the base class)
that the registry, the model inspector, and `resasm template --formulation NAME`
report. `dof_types` is authoritative: the DOF manager builds each node's DOF set
as the ordered union of the `dof_types` of the formulations touching it, so a
truss node gets 3 DOFs and a beam node 6 — automatically, in the same model.

---

## 2. Implement `eval_element`

```python
def eval_element(self, element_id, element_type, coords, dofs,
                 solution_state, material_state, properties,
                 time, dtime, fields, options):
    ...
    return (element_residual,   # (ndof,)      +F_internal contribution to R
            element_tangent,    # (ndof, ndof) or None
            updated_state,      # (n_ip, n_state) or None
            diagnostics)        # dict
```

Rules:
- Return the **`+F_internal`** contribution (the core subtracts `F_external`).
- `element_tangent` must be `d(element_residual)/d(dofs)` — that is what the
  Level-4 finite-difference check verifies. Return `None` if unavailable.
- Read section/material data from `properties` (a `MaterialBinding` with a
  `.section` dict, a plain dict, or `.constants`). Raise a clear error if a
  required property is missing — do not guess.
- Order element DOFs node-major in your declared `dof_types` order; the DOF
  manager gathers exactly those DOFs per node.

---

## 3. Register it

Add the backend in `formulations/registry.py::build_formulation_registry()` (and,
if it should be part of the default set, in
`formulations/__init__.py::default_formulations()`). After that the model
inspector auto-selects it for its `element_types`, and the CLI/API can dispatch to
it.

---

## 4. Verify it (required before claiming support)

Add a test under `tests/framework/` proving, at minimum:
- **zero-field**: residual vanishes for the trivial state (e.g. `U = 0`);
- **known result**: a closed-form check (truss → axial force `EA/L`; beam →
  cantilever tip deflection `PL³/3EI`);
- **FD tangent**: the returned tangent matches a finite-difference of the
  residual.

See `tests/framework/test_truss2_backend.py`, `test_beam2_backend.py`, and
`test_mixed_model_dispatch.py` for the pattern. Only after verification should the
backend be listed as a supported example.

---

## 5. Shell contract (placeholder pattern)

For a formulation you want to declare but not yet implement, follow
`formulations/shell_base.py` + `formulations/shell_placeholder.py`: declare the
full contract (element types, 6 DOF/node, required inputs) with
`supported_modes = ()` so the inspector lists it as *recognized but not runnable*,
and have `eval_element` raise `NotImplementedError` spelling out the contract. This
documents the extension point honestly without pretending it works.
