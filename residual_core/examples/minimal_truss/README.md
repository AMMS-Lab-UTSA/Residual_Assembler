# Minimal example — 2-node truss (bar)

A two-bar truss, authored as a solver-neutral JSON model. Shows the whole
inspect → assemble path with **no Abaqus and no low-level Python**.

- Model: [model.json](model.json) — 3 nodes, 2 `T3D2` bar elements, one section
  material (`E`, `A`).

## Minimal commands

```bash
resasm inspect  residual_core/examples/minimal_truss/model.json
resasm assemble residual_core/examples/minimal_truss/model.json --mode formulation
```

(Use `python -m residual_core.ui.cli ...` if the `resasm` console script is not
installed.)

## Expected output

```
Model inspection summary
------------------------
Elements:
  - T3D2 elements (x2): supported by truss2 backend

Materials:
  - steel_bar: built-in material adapter available

Required user inputs:
  (none — model is ready to assemble)

Possible modes:
  - formulation (self-contained): available for T3D2 (needs section/constitutive properties)
  - stress-driven residual: not applicable
  - material replay: not applicable
  - UEL-direct: not applicable
```

```
assembled mode 'formulation': ndof=9  ||R||=0.000000e+00  max|R|=0.000000e+00
```

(The residual is zero at the undeformed state `U = 0`, as it must be.)

## What was auto-detected

- Element type `T3D2` → `truss2` backend (3 translational DOFs/node).
- The assembly **mode** (`formulation`) is self-contained: no material update,
  no exported field, no UMAT needed.
- Global DOF numbering (9 DOFs = 3 nodes × 3).

## What had to be supplied manually

- The **section properties** `E` and `A`, carried on the material binding in the
  model file. That is the only required input for this mode (see
  `resasm requirements ... --mode formulation`).
