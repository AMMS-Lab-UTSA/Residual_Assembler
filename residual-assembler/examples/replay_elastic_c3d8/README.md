# Collaborator replay sensitivity — runnable example (Program 2)

The complete **collaborator** pipeline: given a JHU package (opaque OTI binary +
manifest) and a saved elastic analysis, compute sensitivities of chosen outputs
w.r.t. chosen material parameters — **without any JHU source**.

```bash
python examples/replay_elastic_c3d8/run_replay.py            # displacement-driven
python examples/replay_elastic_c3d8/run_replay.py --drive force
```

For a self-contained demo the script first *simulates* what a collaborator would
already have — the shipped package and the saved base/perturbed analyses — using
the isolated test fixture (`tests/fixtures/replay_elastic`). The real Program-2
work is the single `run_request(record, manifest, request)` call; results are
finite-difference-validated against the perturbed analyses (worst rel err ~1e-9).

The provider that BUILDS the binary is a different program in a different repo
(`UMAT_source_transformation/oti_provider`); this example never imports it.

Contract documents written next to this file: `replay_record.json`,
`material_manifest.json` (via the fixture), `sensitivity_request.json`,
`sensitivity_result.json`.
