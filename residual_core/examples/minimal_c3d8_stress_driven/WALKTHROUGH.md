# Example 2: a residual assembled from a supplied stress field

One unit-cube C3D8 element. Its eight integration points carry a given
uniaxial stress. Residual_Assembler assembles the nodal residual from that
stress and the mesh alone, and the result is compared with the face tractions
you can compute by hand.

| | |
| --- | --- |
| Folder | `residual_core/examples/minimal_c3d8_stress_driven/` |
| Needs Abaqus | No |
| Needs OTILib | No |
| Run time | under 1 s per command |

## What it demonstrates

A commercial solver reports displacements and stresses, but it keeps its
global residual to itself. Residual_Assembler rebuilds that residual from its
ingredients. The simplest ingredient set is a mesh plus the stress at every
integration point, which is exactly what an output database contains. This is
the *stress-driven* mode: no material model, no parameters, no tangent.

The example matters for two reasons. First, it isolates the assembly: if the
residual is right here, the element integration, the node numbering and the
sign convention are right. Second, it shows how the tool behaves when an
ingredient is missing: it names the single missing item instead of failing.

## The mathematics

The residual is internal minus external force, assembled element by element:

    R = F_int - F_ext,     F_int = sum_e sum_q  w_q det(J_q) B_q^T sigma_q.

Here the cube `[0,1]^3` has `sigma_11 = 100` and every other component zero at
all eight points, and there are no loads, so `F_ext = 0`. For a uniform
stress, the divergence theorem turns the volume integral into face
tractions: `integral B^T sigma dV = integral N^T (sigma n) dA`. Only the faces
`x = 0` and `x = 1` carry a traction, `-100` and `+100` in the x direction.
Each face has area 1 and four nodes, so each node receives a quarter:

    R_x = -25 at the four nodes with x = 0,   R_x = +25 at the four nodes with x = 1,
    R_y = R_z = 0,   ||R|| = sqrt(8 * 25^2) = 70.71067811865476,   max|R| = 25.

## Inputs

| File | What it is |
| --- | --- |
| `model.json` | neutral model format: 8 nodes, one `C3D8` element, no boundary conditions, no loads |
| `fields.json` | `stress_ip`: for element 1, eight rows (one per integration point) of Voigt stress `[S11, S22, S33, S12, S13, S23]` = `[100, 0, 0, 0, 0, 0]` |

## Run it from the command line

From the repository root, with the environment of [INSTALL.md](../../../docs/INSTALL.md):

```bash
export WORK="$HOME/resasm_work"; mkdir -p "$WORK"
M=residual_core/examples/minimal_c3d8_stress_driven

resasm inspect      "$M/model.json" --detail
resasm requirements "$M/model.json" --mode stress-driven
resasm requirements "$M/model.json" --mode stress-driven --fields "$M/fields.json"
resasm assemble     "$M/model.json" --mode stress-driven --fields "$M/fields.json" --out "$WORK/R_cube.npy"
resasm verify       "$M/model.json" --fields "$M/fields.json"
```

What each step tells you (measured output, trimmed):

1. `inspect` finds the element and the backend and lists the possible modes:

   ```text
   Elements:
     - C3D8 elements (x1): supported by solid_c3d8_finite_strain backend
   ...
   Possible modes:
     - stress-driven residual: available if exported stress/resultant fields are provided (C3D8)
   ```

2. `requirements` without the field names the one missing ingredient (exit 0):

   ```text
   Cannot assemble in stress-driven mode.
   Available:
     mesh: yes
     formulation backend: yes
     solution field (U / U+rotation / T): yes
     stress / resultant field: no
     one-step deck in scope: yes

   Minimum missing input:
     provide integration-point stress field S, or an ODB/CSV export (section resultants N/M/Q for beams/shells; heat flux for thermal).
   ```

3. With `--fields` it answers `Ready to assemble in stress-driven mode.` and
   `All minimum inputs are present.`

4. `assemble` builds R and saves it (exit 0):

   ```text
   assembled mode 'stress-driven': ndof=24  ||R||=7.071068e+01  max|R|=2.500000e+01
   ```

5. `verify` splits R into free and constrained parts. This cube has no
   supports and no loads, so the stress is not in equilibrium with anything,
   and the check correctly fails (exit 1):

   ```text
   stress-driven verification:
     ndof                 = 24
     ||R_free||           = 7.071068e+01   (should be ~0 at equilibrium)
     ||reaction (at BC)|| = 0.000000e+00
     absolute tolerance   = 1.000000e-06 (model force units)
     equilibrium          = FAIL
     reaction reference   = NOT CHECKED (assembled reactions only)
   ```

   On a converged analysis with its supports and loads, `||R_free||` would be
   small and the rest of R would sit on the supports as reactions. Choose
   `--atol` for the units and precision of your analysis.

## Run it in the GUI

Start the GUI with `streamlit run scripts/app.py`
([GUI guide](../../../docs/GUI_GUIDE.md)).

- **Fastest:** open the **Start here** tab and press **Run the demo**. It runs
  `inspect`, `requirements`, `assemble` and `verify` on this example, shows
  each command with its output, and offers `demo_R.npy` for download. On the
  measured run the first three were green (exit 0) and `verify` was yellow
  (exit 1), for the reason explained above.
- **Step by step:** on **1. Model**, choose
  `examples/minimal_c3d8_stress_driven` under **Example model** (the field
  file is filled in for you) and press **Run inspect**. On
  **2. Requirements**, keep **Mode** = `stress-driven` and press
  **Run requirements**. On **3. Assemble**, keep `--mode stress-driven`, set
  `--out` (default `R.npy`, written in the sidebar's working directory) and
  press **Run assemble**, then **Download R.npy**, then **Run verify**.

Set the sidebar's **Working directory** to a folder outside the repository
before you assemble, so that `R.npy` is not written into the checkout.

## What it writes

| File | Contents |
| --- | --- |
| `$WORK/R_cube.npy` | the assembled residual, a NumPy array of 24 values ordered node by node (`x, y, z` for node 1, then node 2, ...) |

The other commands print to the terminal only. Nothing in this example is
private: the model and the field are both in the repository.

## Expected output (measured on 2026-09-18)

`R_cube.npy` reshaped to one row per node:

```text
node 1 (0,0,0)  [-25.  0.  0.]
node 2 (1,0,0)  [ 25.  0.  0.]
node 3 (1,1,0)  [ 25.  0.  0.]
node 4 (0,1,0)  [-25.  0.  0.]
node 5 (0,0,1)  [-25.  0.  0.]
node 6 (1,0,1)  [ 25.  0.  0.]
node 7 (1,1,1)  [ 25.  0.  0.]
node 8 (0,1,1)  [-25.  0.  0.]
```

`||R|| = 70.71067811865474`, `max|R| = 25`, as predicted.

## How the result is checked independently

The reference is the hand calculation above, which does not use the element
code at all. This snippet builds it from the node coordinates and compares:

```bash
python - "$WORK/R_cube.npy" residual_core/examples/minimal_c3d8_stress_driven/model.json <<'EOF'
import json, sys, numpy as np
R = np.load(sys.argv[1])
nodes = json.load(open(sys.argv[2]))["nodes"]
ref = np.concatenate([[25.0 if nodes[n][0] == 1.0 else -25.0, 0.0, 0.0]
                      for n in sorted(nodes, key=int)])
print("max |R - face tractions| =", np.abs(R - ref).max())
EOF
```

Measured: `max |R - face tractions| = 3.552713678800501e-15`, that is, agreement
to rounding. The offline test suite checks the same element further with
patch tests and a frame-objectivity test ([STATUS.md](../../../STATUS.md)).

## Run time

Measured on 2026-09-18: `resasm assemble` 0.46 s wall time; the other commands
are similar. The Start-here demo in the GUI reported 0.01 to 0.03 s per
command, because the GUI calls the command line in-process.

## Common problems

| Symptom | Cause and fix |
| --- | --- |
| `Cannot assemble in stress-driven mode ... Minimum missing input: provide integration-point stress field S`, exit 2 | `--fields` is missing or points to the wrong file. Give the JSON field export. |
| `equilibrium = FAIL`, exit 1 | Expected for this unloaded cube. For your own model, check that the field and the mesh describe the same converged state, and choose an `--atol` suited to its units. |
| `ERROR: cannot read the field export ...: No such file or directory`, exit 2 | The file given to `--fields` does not exist. Check the path. |
| `ERROR: field export ...: key 'schema' is not an element id; expected {"stress_ip": ...}`, exit 2 (or `is not valid JSON`, `holds no integration-point stress`) | The file given to `--fields` is not in the `stress_ip` layout (for example, it is a model file or the export written by `resasm request`). The layout is `{"stress_ip": {"<element id>": [[S11, S22, S33, S12, S13, S23], ... one row per point]}}`. |
| You want to use an `.odb` directly | `--odb` on `assemble`, `requirements` and `verify` is only another name for `--fields`; given a binary ODB it stops with `ERROR: field export 'my_job.odb' is not a .json file: binary ODB reading needs Abaqus` (exit 2). Export the ODB once with `abaqus python scripts/extract_odb_fields.py --odb my_job.odb --out fields.json` (Abaqus 2021 Python 2.7 or Python 3) and give that file. The ODB-reading sensitivity workflow is `resasm request` ([Example 3](../../../examples/presentation_request/WALKTHROUGH.md)). |

Previous: [Example 1](../../../examples/user_config_minimal/WALKTHROUGH.md).
Next: [Example 3](../../../examples/presentation_request/WALKTHROUGH.md).
All examples: [examples/README.md](../../../examples/README.md).
