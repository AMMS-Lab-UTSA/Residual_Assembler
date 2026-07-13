# Minimal example — 2-node 3D beam (frame)

A single cantilever beam, authored as a solver-neutral JSON model. Exercises the
**rotational DOFs** (6 per node) that prove the core is not displacement-only.

- Model: [model.json](model.json) — 2 nodes, 1 `B31` beam element, one section
  material (`E`, `A`, `Iz`, `Iy`, `G`, `J`).

## Minimal commands

```bash
resasm inspect  residual_core/examples/minimal_beam/model.json
resasm assemble residual_core/examples/minimal_beam/model.json --mode formulation
```

## Expected output

```
Model inspection summary
------------------------
Elements:
  - B31 elements (x1): supported by beam2 backend

Materials:
  - beam_section: built-in material adapter available

Required user inputs:
  (none — model is ready to assemble)

Possible modes:
  - formulation (self-contained): available for B31 (needs section/constitutive properties)
  - stress-driven residual: not applicable
  - material replay: not applicable
  - UEL-direct: not applicable
```

```
assembled mode 'formulation': ndof=12  ||R||=0.000000e+00  max|R|=0.000000e+00
```

(12 DOFs = 2 nodes × 6; residual is zero at `U = 0`.)

## What was auto-detected

- Element type `B31` → `beam2` backend (`UX,UY,UZ,RX,RY,RZ` per node).
- Heterogeneous DOF set: each node carries 6 DOFs, unlike the truss's 3.
- The `formulation` mode (self-contained, no material update / field / UMAT).

## What had to be supplied manually

- The **section properties** (`E`, `A`, `Iz` required; `Iy`, `G`, `J` optional)
  on the material binding. Run
  `resasm template --formulation beam2` to see the declared required/optional
  inputs.
