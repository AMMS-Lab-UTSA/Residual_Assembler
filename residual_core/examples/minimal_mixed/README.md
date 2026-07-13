# Minimal example — mixed truss + beam (heterogeneous DOFs)

One model containing **both** a `T3D2` bar (3 DOFs/node) and a `B31` beam
(6 DOFs/node) on disjoint nodes. Proves the physics-blind core dispatches each
element to its own backend and builds a per-node DOF set from the union of the
formulations touching that node.

- Model: [model.json](model.json) — 4 nodes, 1 bar (nodes 1–2), 1 beam (3–4).

## Minimal commands

```bash
resasm inspect  residual_core/examples/minimal_mixed/model.json --detail
resasm assemble residual_core/examples/minimal_mixed/model.json --mode formulation
```

## Expected output

```
Model inspection summary
------------------------
Elements:
  - B31 elements (x1): supported by beam2 backend
  - T3D2 elements (x1): supported by truss2 backend

Materials:
  - bar: built-in material adapter available
  - frame: built-in material adapter available

Required user inputs:
  (none — model is ready to assemble)

Possible modes:
  - formulation (self-contained): available for B31, T3D2 (needs section/constitutive properties)
  ...
```

```
assembled mode 'formulation': ndof=18  ||R||=0.000000e+00  max|R|=0.000000e+00
```

(18 DOFs = 2 truss nodes × 3 + 2 beam nodes × 6. See `--detail` for the
per-element backend selection.)

## What was auto-detected

- `T3D2` → `truss2` (3 DOFs/node); `B31` → `beam2` (6 DOFs/node), in **one** model.
- The DOF manager assigns 3 DOFs to the truss nodes and 6 to the beam nodes
  automatically — no assumption that all nodes have the same DOFs.
- Physics-blind dispatch + scatter into the correct global indices.

## What had to be supplied manually

- Section properties for each material (bar: `E`, `A`; frame: `E`, `A`, `Iz`, …).
