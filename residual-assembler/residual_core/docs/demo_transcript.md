# Demo Transcript

Real, captured output of the framework's command line — no Abaqus, no external
code. Reproduce with the commands shown. See
[`scripts/demo_framework.sh`](../../scripts/demo_framework.sh) for the full
scripted run (`python -m residual_core.ui.cli ...` works if the `resasm` console
script is not installed).

## 1. Inspect the registered backends (agnostic by construction)

```text
$ resasm backends
Registered formulation backends
===============================

beam2
  element types : B31, B33, B31_LIKE, BEAM2, FRAME3D
  dof types     : UX, UY, UZ, RX, RY, RZ
  status        : implemented/simple
  available modes: formulation
  material iface : False
  state         : none
  tangent       : analytic
  limitations   : small-strain Euler-Bernoulli; 2 nodes, 6 DOF/node; no shear (no Timoshenko)
... (truss2, solid_c3d8_small_strain, solid_c3d8_finite_strain, stress_driven_c3d8,
     uel_direct, shell_placeholder) ...
Materials: crystal_plasticity, isotropic_elastic, umat
```

Crystal plasticity is **one entry** among seven formulation backends.

## 2. Inspect a truss model (auto-detection)

```text
$ resasm inspect residual_core/examples/minimal_truss/model.json
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

## 3. Assemble the truss residual (no Abaqus, no low-level objects)

```text
$ resasm assemble residual_core/examples/minimal_truss/model.json --mode formulation
assembled mode 'formulation': ndof=9  ||R||=0.000000e+00  max|R|=0.000000e+00
```

Residual is zero at the undeformed state, as required.

## 4. Requirements engine — the minimum missing input, not a checklist

```text
$ resasm requirements residual_core/examples/minimal_c3d8_stress_driven/model.json --mode stress-driven
Cannot assemble in stress-driven mode.
Available:
  mesh: yes
  solution field (U / U+rotation / T): yes
  stress / resultant field: no

Minimum missing input:
  provide integration-point stress field S, or an ODB/CSV export (section
  resultants N/M/Q for beams/shells; heat flux for thermal).
```

## 5. Assemble a solid residual once the field is attached

```text
$ resasm assemble residual_core/examples/minimal_c3d8_stress_driven/model.json \
        --mode stress-driven --fields residual_core/examples/minimal_c3d8_stress_driven/fields.json
assembled mode 'stress-driven': ndof=24  ||R||=7.071068e+01  max|R|=2.500000e+01
```

The nonzero residual is the internal nodal force from the supplied uniaxial
stress — assembled with no knowledge of any material model.

## 6. External residual comparison (offline, once a field export exists)

```text
$ python scripts/compare_residuals.py --model .../model.json --fields .../fields.json
residual comparison (mode=stress-driven)
  ndof                 = 24
  ||R_free||           = 7.071068e+01
  ...
```

## 7. Abaqus scripts skip cleanly when Abaqus is absent

```text
$ python scripts/extract_odb_fields.py --odb job.odb --out fields.json
Abaqus not available: validation pending
$ python scripts/compare_umat_replay.py --reference ref.json
Abaqus not available: validation pending (no reference export: ref.json)
```

Nothing here fails the offline suite.
