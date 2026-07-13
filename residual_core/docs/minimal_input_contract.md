# Minimal Input Contract

The framework asks for **the least data needed** to assemble a residual in a
chosen *mode*, then reports the single most important missing item rather than a
generic checklist. This contract defines, per mode, what is required and what is
optional. The requirements engine (`core/requirements.py`) enforces it; the model
inspector (`core/diagnostics.py`) reports which modes a given model can reach.

The **element field that enters the weak form is mode- and physics-dependent**
and intentionally generic:

- solids → integration-point **stress** (Voigt `S`),
- beams/shells → **section resultants** (`N`, `M`, `Q`),
- thermal → **heat flux**.

The engine reasons over declared input *keys* only; it never assumes mechanics.

---

## Mode 1 — Stress-driven (field supplied externally)

Assemble `f_int = Σ Bᵀ σ dV` from a field exported by another solver. Verifies
the finite-element assembly **independently of any material model**.

| Requirement | Key | Required? |
|---|---|---|
| Mesh: connectivity + coordinates | `mesh` | required |
| DOF solution field (U, or U+θ for beams/shells, T for thermal) | `dof_field` | required |
| Element field entering the weak form (IP stress / resultants / flux) | `element_field` | required |
| Config (`small` vs `finite`) | `config` | optional (default per backend) |

No material model, no material parameters, no state, no tangent needed.

**Attach with:** `problem.attach_results(fields_or_path)` where the export is a
dict `{eid: (n_ip, 6)}` (or `{'stress_ip': {...}}`), or a JSON file. Binary ODB
reading requires Abaqus (`io/abaqus_odb_export.py`).

---

## Mode 2 — Material-replay (material update drives the field)

Recompute the field by marching the constitutive law over the loading history,
then assemble. This is the full constitutive + assembly pipeline offline.

| Requirement | Key | Required? |
|---|---|---|
| Mesh: connectivity + coordinates | `mesh` | required |
| Solution history (increment sequence, **not just the final step**) | `solution_history` | required for history-dependent materials |
| Material model or adapter (UMAT source, built-in law, or Python law) | `material_model` | required |
| Material parameters (PROPS) | `material_parameters` | required |
| Previous state variables (STATEV) or an initial state | `state_prev` | required for stateful materials |
| Time increments (`dtime` per step) | `time_increments` | required |

**Honesty note.** A material is **not** always a UMAT. A stateless built-in
elastic law needs no history and no state. A crystal-plasticity **UMAT** is
history-dependent and, for genuine stress, needs the compiled Fortran
(Intel `ifort` + Abaqus) — offline the framework runs it only via a Python mock
backend, or you use Mode 1 instead. The framework reports this rather than
pretending.

---

## Mode 3 — Direct residual / UEL (element returns its own residual)

The element math lives entirely inside a user-supplied routine (an Abaqus-UEL-like
callable). No material update, no exported IP field.

| Requirement | Key | Required? |
|---|---|---|
| Mesh: connectivity + coordinates | `mesh` | required |
| Element DOF layout | `element_dofs` | required |
| Callable residual routine returning `(RHS, AMATRX, SVARS_new)` | `uel_routine` | required |

Sign convention is normalized by `formulations/uel_adapter.py` (Abaqus
`RHS = −R`; the adapter returns `+R`).

---

## Mode 4 — Formulation / manufactured (backend computes everything)

A pure formulation backend (e.g. truss, beam) computes the residual from geometry
and section properties alone — used for verification and for demonstrating new
backends. Aliased as `manufactured` / `manual`.

| Requirement | Key | Required? |
|---|---|---|
| A formulation backend implementing the base interface | `formulation_backend` | required |
| Section / constitutive properties (e.g. `E`, `A`, `I` for a beam) | `section_properties` | required |

---

## Reading a requirements report

```
resasm requirements model.inp --mode stress-driven
```

prints what is present, what is missing, and the **minimum next item** to
provide. `resasm doctor model.inp` prints per-mode readiness for every mode at
once and can emit a pre-filled config template with
`--write-config-template config.yml`.
