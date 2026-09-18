# Presentation Interface

The primary collaborator interface implements slide 12 of the original
`IMQCAM_Annual_Meeting_20260812_V4.pptx`. XML for slides 9-12 and 14-18 was read.
PPTX SHA-256: `55688498a351ea403788b49faf7a7ac26758c84483588fa54dd3e95f0e8ae81d`.
Exact interface names do **not** imply support for FCC on slides 15-18.
This implementation is bounded to the fingerprint-pinned m3_j2 model.

## Collaborator Command

For a wheel installation, use the installed console command from a directory
containing the shared artifacts (no source checkout or PYTHONPATH):

```sh
resasm request --model Analysis.inp --odb Analysis.odb \
  --material OTI_UMAT.obj --request sensitivity_request.json --out new_results
```

Keep the generated `Mapping.json` beside the object. The four explicit inputs
and three public output names remain unchanged. Installation, installed GUI
launch commands, and the executable joint wheel gate are documented in
[evidence/recovery_install.md](evidence/recovery_install.md).

For development from the recovery checkout only:

```sh
export PY=/home/ammslab3/softwarex_work/.venv/bin/python
export PYTHONPATH="$PWD:../imq-umat-recovery/src"
"$PY" -m residual_core.ui.cli request \
  --model /path/to/Analysis.inp --odb /path/to/Analysis.odb \
  --material /path/to/OTI_UMAT.obj \
  --request /path/to/sensitivity_request.json --out /path/to/new_results
```

The output directory must be empty or new. No production solve, synthetic
solution, source transformation, or manual replay record is involved.
Licensed `abaqus python` is invoked for extraction only. Missing binary error:
`Abaqus executable '<name>' not found on PATH; licensed Abaqus Python with
odbAccess is required to read Analysis.odb`. `--abaqus PATH` overrides the
executable. Failed extraction includes the command, exit code and diagnostic.
The binary-compatible linker/compiler is required (gfortran on verified Linux),
but **material source is not needed**.

GUI: `"$PY" -m streamlit run scripts/app.py` with the same environment.
**Sensitivity Request** is the primary screen, with the four exact file labels,
uploads and local paths. Advanced contains the optional mapping selector.
Downloads have the three exact slide names. It calls the existing CLI bridge:
`cli.main -> cmd_request -> presentation.run_request`, the same application
service as the CLI. Legacy functionality remains under **Advanced Replay**.

## Slide-To-Implementation Map

| Input/output | Actual field/tool | Test in `tests/integration/test_presentation_request.py` |
| --- | --- | --- |
| `OTI_UMAT.obj` | `--material`; `mapping_for`, `PathMaterial`, `require_j2`; compiled EVAL/MARCH | `test_mapping_discovery_hash_and_alias`, `test_real_archived_j2_export` |
| `Analysis.inp` | `--model`; `presentation_inputs.read_model` wraps existing `parse_inp`; mesh, properties, BCs, CLOAD, time | `test_no_dropped_inp_physics`, `test_malformed_boundaries_not_coerced` |
| `Analysis.odb` | `--odb`; existing `abaqus_odb_export.py --frames all --strict yes`; U/RF/CF/S/SDV1 | `test_conversion_rejects_mismatches`; separate genuine ODB proof |
| `sensitivity_request.json` | `--request`; `validate_request`, `scalar_results`, `reduce_scalar` | `test_selection_after_total_history`, `test_invalid_request_scope`, `test_scalar_reductions` |
| `sensitivity_results.json` | Root output: original request, resolved scope, metadata, scalar values/derivatives | `test_compiled_offline_conversion_and_all_outputs` |
| `sensitivity_tables.csv` | Root output: output, field, component, reduction, domain, increment, time, value, parameter, derivative | Same test checks 16 rows and correct derivative |
| `run_report.txt` | Root output: executed versus verified, residual/solver checks, tolerances, limitations; failure diagnostic | CLI malformed-input test and genuine source-denied run |

Application code is in `residual_core/replay/presentation.py`, input conversion
in `presentation_inputs.py`, scalar math in `request_reductions.py`.
Full fields, K/R/derivatives, automatically generated record, link shim/library,
export command and log stay in `private/`. Public JSON contains no full arrays.

## Provider And Mapping

Slide 11 has the developer share `REAL_UMAT.obj`, `OTI_UMAT.obj` and completed
`Mapping.json`; source stays private. The request consumes only the OTI object.
The current recovery producer bundles ORIGINAL binary and OTI symbols. ORIGINAL
is used for the virgin elastic tangent, not source access or production FE rerun.

Generator: `../imq-umat-recovery/src/umat_oti/provider/build.py`.
Contract: `resasm_umat_oti_contract_v1`, with dimensions, symbols, PROPS indices,
OTI directions, derivative/Voigt layouts, source fingerprint and object SHA-256.
Use the generated completed sidecar, not the input transformation contract.

Discovery checks `<object-stem>.json` and adjacent `Mapping.json`. `Mapping.json`
is an alias of the **unchanged generated sidecar**, not a new schema. Renaming
the object requires matching `object.sha256_full`; source fingerprint, layouts,
dimensions and parameter order/directions are checked too. Conflicting sidecars
require `--mapping PATH`. Provenance hashes are not certification of untrusted
binaries; run only trusted providers. This task did not modify the generator.

## Request Format

See `examples/presentation_request/sensitivity_request.json`. Top-level keys
are exactly `outputs`, `parameters`, `domain`, and `increments`.

```json
{
  "outputs": [
    {"name": "loaded_U1", "field": "U", "component": 1, "reduction": "mean"}
  ],
  "parameters": ["E", "SIGY0", "H"],
  "domain": {"nodes": [2, 3, 6, 7]},
  "increments": "LAST"
}
```

- Parameters: nonempty unique subset of `E`, `nu`, `SIGY0`, `H`.
- Fields: `U`, `RF`, `S`, `SDV`. One-based component or `ALL`.
  Stress order: 11,22,33,12,13,23. SDV1: equivalent plastic strain.
- Domain: `nodes` for U/RF, `elements` for S/SDV; `ALL` or explicit existing
  integer ids. Each output may override `domain`. All eight IPs are included.
  No implicit node/element conversion; no IP-subset selector yet.
- Increments: `ALL`, `LAST`, or a unique list such as `[2,4]`, one-based after
  virgin increment 0. **Replay the entire preceding history before selection.**
- Reductions: `component` requires one scalar location; `sum` and `mean` are
  unweighted; `L2` is Euclidean norm, not von Mises; `max` is signed maximum.
  L2 derivative is `values @ derivatives / norm`; zero norm fails. Max uses
  the unique active derivative; ties/near ties fail (`rtol=1e-12, atol=1e-14`).
  Mean is not volume weighted; sum is not a spatial integral.

Unknown fields/options, invalid ids, duplicate names/parameters/increments,
unsupported reductions and nondifferentiable cases fail, never disappear.
Von Mises, mixed-output arithmetic, shape/load/BC and higher derivatives are
not supported. Ordinary execution never runs FD and says `verified=false`.
Explicit `--validate` requests independent ORIGINAL whole-history FD; its
strict double-precision gate may reject float32 ODB data. It launches no Abaqus
production job. The analytic real-ODB proof is separate from ordinary execution.

## Physics And Precision

Supported: one named `*Static`, NLGEOM=NO step; one untransformed instance;
C3D8/B-bar; homogeneous pinned m3_j2 (four constants, one SDV); virgin state;
zero fixed displacements; parameter-independent ramped nodal concentrated loads.
ODB must save all nodes/elements and every increment through the step end,
including frame 0, with U/RF/CF/S/SDV1. CF verifies actual loads against the INP.
Mesh ids/coordinates/connectivity, time, loads, BCs, IP stress/state, reactions
and free equilibrium are checked. Missing frames are not interpolated.

Pressure/body loads, amplitudes, includes, equations, contact, dynamics, finite
strain, multiple steps/instances/materials, nonzero BCs and initial state fail.
This is a dense bounded engine, not a generic UMAT or scalability claim.

The connected engine carries total solution-history derivatives:
`du_dp (ndof,4)`, `dsigma_dp (nelem,8,6,4)`, `dstate_dp (nelem,8,1,4)`.
Reactions use `Rp_history + K @ du_dp`. The old single-step
`request_contract.replay_sensitivities` is never called.

ODB float32 displacements are replayed, not re-equilibrated. Acceptance:
scaled free residual below `1e-5`; stress/reaction `rtol=2e-5` and absolute
tolerance `2e-5 * max(max_abs_recorded_field,1)`; state `rtol=2e-5, atol=1e-8`.
Zero-BC residues below `1e-12 * max_mesh_extent` are checked then set to zero.
All tolerances are public. In the genuine reference, displacement rounding
produces transverse stress error `1.23e-4`; the expected amplification
`eps32 * E * max(strain)` is approximately `6e-4`. Synthetic default replay
retains `1e-14` equilibrium tolerance. Unit-dependent absolute floors and ODB
precision bound accuracy; inputs outside gates fail without convergence repair.

## Reproduction

`scripts/reproduce_presentation_request.py prepare --work .../imqrp_reference`
builds the provider and launches one **developer-only** licensed job. Use an
Intel compiler environment, a new `imqrp_` directory under
`../imq_abaqus/recovery_presentation`, and never prepare for ordinary use.
`consume --work ... --out NEW` uses the five shared artifacts, denying material
source reads in the collaborator process. Only its generated ABI shim is
exempt; this audit guard is not an OS sandbox. No source is deleted.

Small genuine exported fields are retained in `tests/fixtures/presentation_j2`,
not an ODB binary. `scripts/check_presentation_browser.py` tests actual browser
execution/downloads and stops its exact server PID. Exact commands, counts and
measured limitations: `docs/evidence/recovery_presentation.md`.