# Worked examples

Seven complete examples, ordered from the smallest calculation to full-size
Abaqus models. Each has a walkthrough that explains what it shows, gives the
exact commands and the equivalent GUI steps, lists the files it writes, quotes
the output measured on 2026-09-18, and says how the result is checked
independently.

| # | Example | What it shows | Needs Abaqus? | One-line command (from the repository root) |
| --- | --- | --- | --- | --- |
| 1 | [The smallest sensitivity](user_config_minimal/WALKTHROUGH.md) | a residual you write yourself (a cubic spring); first and second derivatives against the closed form | No (needs OTILib) | `resasm init --template python --out "$WORK/spring" && resasm run "$WORK/spring/resasm.yml"` |
| 2 | [Stress-driven C3D8 assembly](../residual_core/examples/minimal_c3d8_stress_driven/WALKTHROUGH.md) | the residual of one element from a supplied stress field, against the analytic face tractions | No | `resasm assemble residual_core/examples/minimal_c3d8_stress_driven/model.json --mode stress-driven --fields residual_core/examples/minimal_c3d8_stress_driven/fields.json --out "$WORK/R_cube.npy"` |
| 3 | [The four-file request on one element](presentation_request/WALKTHROUGH.md) | `Analysis.inp` + `Analysis.odb` + `OTI_UMAT.obj` + `sensitivity_request.json` to sensitivities; J2 plasticity against the uniaxial closed form | Yes (the ODB is not shipped; Abaqus Python reads it) | in the folder holding the five input files: `resasm request --model Analysis.inp --odb Analysis.odb --material OTI_UMAT.obj --request sensitivity_request.json --out results` |
| 4 | [History replay of a J2 beam, offline](replay_history/WALKTHROUGH.md) | 10 plastic increments; ODB parity, Abaqus finite differences, the homogeneity identity, weighted shares, `--reequilibrate` | No (the Abaqus export is committed) | `resasm history --model examples/replay_history/j2_beam/Analysis.inp --fields examples/replay_history/j2_beam/fields.npz --material "$WORK/provider_j2/umat_m3_j2_oti.obj" --request examples/replay_history/j2_beam/sensitivity_request.json --out "$WORK/beam"` |
| 5 | [Full-size cantilevers](cantilevers/WALKTHROUGH.md) | J2 (1,536 C3D8, 40 increments) and FCC crystal plasticity (384 C3D8, 25 increments, 10 parameters); full-field von Mises sensitivities and per-increment weighted shares | Yes, once (to run the two analyses) | with `C` the work folder of the walkthrough: `resasm request --model "$C/j2/cantilever_j2_nominal.inp" --odb "$C/j2/cantilever_j2_nominal.odb" --material "$WORK/provider_j2/umat_m3_j2_oti.obj" --request examples/cantilevers/j2_request.json --out "$C/j2_results"` |
| 6 | [Provider-to-sensitivity pipeline](bounded_j2_c3d8/WALKTHROUGH.md) | builds the compiled J2 material, solves a cyclic one-element history, verifies every derivative against whole-model finite differences | No | `python scripts/reproduce_connected_pipeline.py --skip-abaqus --out "$WORK/pipeline"` |
| 7 | [Finite-strain neo-Hookean C3D8](finite_strain_c3d8/WALKTHROUGH.md) | finite-strain assembly with the exact tangent and OTILib parameter sensitivities, against nonlinear re-solves | No (needs OTILib) | `python examples/finite_strain_c3d8/benchmark.py --out "$WORK/finite"` |

## Before you run an example

1. Install both packages as described in [docs/INSTALL.md](../docs/INSTALL.md):
   Residual_Assembler and UMAT-OTI side by side, in one virtual environment.
2. Run the commands from the root of the Residual_Assembler checkout, and send
   every output to a folder outside the repository:

   ```bash
   export RA="$PWD"
   export UMAT="$RA/../UMAT_source_transformation"
   export WORK="$HOME/resasm_work"
   mkdir -p "$WORK"
   ```

3. Examples 3, 4 and 5 use a compiled material provider. Build the J2 one once
   (about 6 s):

   ```bash
   umat-oti-provider build "$UMAT/parameter_sensitivity/models/m3_j2/contract_v2.json" --out "$WORK/provider_j2"
   ```

4. Examples 1 and 7 use OTILib. Point `PYOTI_PATH` and `OTILIB_ROOT` to its
   build directory ([INSTALL.md](../docs/INSTALL.md#4-install-otilib-for-the-direct-residual-paths)).
5. Examples 3 and 5 read Abaqus output databases. Reading an `.odb` needs
   Abaqus Python (`abaqus python`, tested with Abaqus 2021.HF5). Once a model
   is exported, `resasm history --fields` works without Abaqus.

## Which example first?

- To see the idea on one line of algebra, start with **Example 1**.
- To see how a residual is built from a mesh, run **Example 2**.
- To check that both packages work together, run **Example 6** (no Abaqus,
  no OTILib, 10 s).
- For the workflow on Abaqus results, read **Example 3**, then run
  **Example 4**, which needs no Abaqus, and finally **Example 5**.

## How each result is checked

| # | Independent reference | Measured agreement (2026-09-18) |
| --- | --- | --- |
| 1 | closed form of `u = (f/k)^(1/3)` | 2.8e-17 or better, all five derivatives |
| 2 | face tractions of a uniform stress | 3.6e-15 |
| 3 | uniaxial J2 closed form | 1.7e-7 relative (`dU1/dE`); to about 1e-14 after `resasm history --reequilibrate` |
| 4 | ODB parity, Abaqus finite differences, homogeneity identity, whole-model finite differences | all 160 Abaqus comparisons within their uncertainty; identity 9.4e-15; finite differences 1.5e-7 |
| 5 | ODB parity, homogeneity identity, whole-model finite differences at full size | identity 1.4e-12 (J2), 1.2e-13 (FCC); `SIGY0` finite differences 7.9e-9 or better where resolved |
| 6 | whole-model finite differences of the ORIGINAL routine; an archived elastic Abaqus export | 9.7e-7 or better (tolerance 2e-6); Abaqus fixture passed |
| 7 | first-Piola quadrature, tangent finite differences, nonlinear re-solves | 5.7e-16; 4.6e-11; 1.0e-10 |

## Run times

Measured on 2026-09-18 on a 24-core Linux workstation with Python 3.11:
Example 1, 1 s; Example 2, 0.5 s per command; Example 3, 1.1 s (plus 5.5 s
to build the provider); Example 4, 1 s; Example 5, 25 to 80 s per model and
run; Example 6, 10 s; Example 7, 3 s.

## More small models

The folders under `residual_core/examples/` hold further neutral-format models
that the GUI's **1. Model** tab offers directly:
[minimal_truss](../residual_core/examples/minimal_truss/README.md),
[minimal_beam](../residual_core/examples/minimal_beam/README.md),
[minimal_mixed](../residual_core/examples/minimal_mixed/README.md) (all three
assemble with `resasm assemble <model.json> --mode formulation`) and
[minimal_nonlinear_spring_sensitivity](../residual_core/examples/minimal_nonlinear_spring_sensitivity/README.md)
(`resasm sensitivity <model.json> --params <params.json>`; measured
`du/dk = -3.333333e-01`, equal to the analytic value). The job templates for
black-box, C++ and Fortran residuals are in [templates/](../templates/).
`examples/residual_sensitivity_c3d8/sensitivity_engine.py` is a short,
self-checking script of the sensitivity equation on one elastic C3D8
(`python examples/residual_sensitivity_c3d8/sensitivity_engine.py`).

The commands are explained in [docs/CLI_GUIDE.md](../docs/CLI_GUIDE.md) and the
graphical interface in [docs/GUI_GUIDE.md](../docs/GUI_GUIDE.md).
