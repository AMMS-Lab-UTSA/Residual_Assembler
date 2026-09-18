# Bounded J2 C3D8 Example

`model.json` is a self-contained one-element, eight-IP, seven-increment cyclic
model, not the presentation cantilever or an Abaqus export. Coordinates are a
unit cube; the bottom is fixed and the top has parameter-independent nodal
forces. Material constants `[210000, 0.3, 250, 2000]` come from the W3 m3_j2
contract. Consistent units must be used throughout.

From the repository root, with the environment in
[the workflow](../../docs/IMQCAM_WORKFLOW.md):

```sh
"$PY" scripts/reproduce_imqcam_pipeline.py --skip-abaqus \
  --out .pytest_cache/j2_example
```

Expected: exit 0, `manifest.json` reports `passed=true`, 7 increments and 8 IPs,
four verified parameter columns, and FD scaled errors below `2e-6` at all three
steps with a last-two-step plateau check. The output includes a plot, CSV table,
fresh compiled sources/object, converged record and K/R/Rp/state/derivative arrays.
Failure exits 2 and records the diagnostic in the manifest.

The independent references are ORIGINAL whole-material paths, ORIGINAL
whole-model re-equilibration, and the separately shipped genuine elastic Abaqus
fixture. J2 results are explicitly labelled `synthetic_converged_fe`. The
fixed-path local solve is retained only to demonstrate that it is not the
total-history displacement derivative.

This is a dense small-strain algebra verification. Strong cyclic loads can
accumulate about 10% displacement on the unit cube; no finite-strain or physical
large-deformation accuracy is claimed. Exact measured results and remaining
dataset limitations are in [recovery evidence](../../docs/evidence/recovery_integration.md).