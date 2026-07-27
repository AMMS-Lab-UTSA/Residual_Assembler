# Minimal example — C3D8 solid, stress-driven

A single unit-cube `C3D8` solid. The residual is assembled from an **externally
supplied** integration-point stress field (uniaxial `S11 = 100`), with no
material model and no UMAT — the stress-driven verification path.

- Model: [model.json](model.json) — 8 nodes, 1 `C3D8` element.
- Field: [fields.json](fields.json) — `stress_ip` for element 1, `(8, 6)` Voigt.
  The model references it via `field_refs` but the data lives in the field file.

## Minimal commands

```bash
resasm inspect      residual_core/examples/minimal_c3d8_stress_driven/model.json --detail
resasm requirements residual_core/examples/minimal_c3d8_stress_driven/model.json --mode stress-driven
resasm assemble     residual_core/examples/minimal_c3d8_stress_driven/model.json \
                    --mode stress-driven --fields residual_core/examples/minimal_c3d8_stress_driven/fields.json
```

## Expected output

Without a field, the tool reports the single minimum missing input (it does not
crash and does not dump a generic checklist):

```
Cannot assemble in stress-driven mode.
Available:
  mesh: yes
  solution field (U / U+rotation / T): yes
  stress / resultant field: no

Minimum missing input:
  provide integration-point stress field S, or an ODB/CSV export (section
  resultants N/M/Q for beams/shells; heat flux for thermal).
```

With the field attached:

```
assembled mode 'stress-driven': ndof=24  ||R||=7.071068e+01  max|R|=2.500000e+01
```

(24 DOFs = 8 nodes × 3. The nonzero residual is the internal nodal force from the
uniaxial stress — the surface tractions on the cube.)

## What was auto-detected

- Element type `C3D8` → `stress_driven_c3d8` backend for this mode.
- The stress-driven **mode** requires only mesh + the exported field; no material,
  no PROPS, no state, no tangent.
- The `field_refs` entry in the model points at the export; the CLI reads it via
  `--fields`.

## What had to be supplied manually

- The **exported stress field** (`fields.json`). This is the one required input
  the framework cannot infer — it comes from a prior solver run/export.
