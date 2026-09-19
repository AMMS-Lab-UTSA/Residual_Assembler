# The `resasm request` interface

This page specifies `resasm request`, the four-input command that computes
parameter sensitivities of a finished Abaqus analysis: its inputs and outputs,
the request format, the scope and precision of its bounded engine, and how a
model outside that scope is handed to the history engine. It is for users
preparing a request and for reviewers checking what the command accepts. For a
first run see [QUICKSTART_USER.md](../QUICKSTART_USER.md); for every option,
[CLI_GUIDE.md](CLI_GUIDE.md).

## The command

With the packages installed, from a directory holding the shared files (no
source checkout or `PYTHONPATH` needed):

```sh
resasm request --model Analysis.inp --odb Analysis.odb \
  --material OTI_UMAT.obj --request sensitivity_request.json --out new_results
```

Keep the generated `Mapping.json` beside the object. Installation, including
the clean-install gate that checks this command from built wheels, is in
[INSTALL.md](INSTALL.md).

From source checkouts (the Residual_Assembler root, with
UMAT_source_transformation beside it):

```sh
export PYTHONPATH="$PWD:../UMAT_source_transformation/src"
python -m residual_core.ui.cli request \
  --model /path/to/Analysis.inp --odb /path/to/Analysis.odb \
  --material /path/to/OTI_UMAT.obj \
  --request /path/to/sensitivity_request.json --out /path/to/new_results
```

The output directory must be empty or new. No production solve, synthetic
solution, source transformation or hand-made replay record is involved.
Licensed `abaqus python` is invoked only to extract the ODB. If it is
missing, the command stops with `request failed: Category: odb_export`, and
`private/error_report.txt` says `Abaqus executable '<name>' not found on PATH;
licensed Abaqus Python with odbAccess is required to read Analysis.odb`.
`--abaqus PATH` selects another executable, and a failed extraction reports
the command, exit code and diagnostic. A binary-compatible linker and compiler are required
(gfortran on the verified Linux platform), but **the material source is not**.

## Which engine runs

`resasm request` decides which engine a model needs
(`residual_core/ui/cmd_history.py::route_request`) and says which one it used:

- **Bounded engine** (`residual_core/replay/presentation.py`): one homogeneous
  pinned `m3_j2` J2 model with zero-valued boundaries, concentrated loads and
  the four-key request below (fields U/RF/S/SDV; reductions component, sum,
  mean, L2, max). Its scope and precision are described on this page.
- **History engine** (`resasm history`): models outside that scope, for
  example nonzero prescribed displacements, many increments, `*Controls`,
  a provider other than the pinned `m3_j2`, node and element sets, von Mises
  outputs and the other request extensions, and the full-size cantilevers of
  [examples/cantilevers](../examples/cantilevers/README.md). The first line of
  output names the reason for the hand-over, and `--validate` becomes
  `--verify fd`. The exact rule is in
  [REPLAY_HISTORY.md](REPLAY_HISTORY.md#how-resasm-request-chooses-this-engine),
  with the engine's scope, mathematics and tolerances.

Features neither engine supports (pressure and body loads, contact,
amplitudes, several steps, materials or instances, initial state, finite
strain) are refused and named in `run_report.txt`.

## Inputs and outputs

| Input or output | Field and implementation | Test in `tests/integration/test_presentation_request.py` |
| --- | --- | --- |
| `OTI_UMAT.obj` | `--material`; `mapping_for`, `PathMaterial`, `require_j2`; compiled EVAL/MARCH | `test_mapping_discovery_hash_and_alias`, `test_real_archived_j2_export` |
| `Analysis.inp` | `--model`; `presentation_inputs.read_model` wraps `parse_inp`: mesh, properties, BCs, CLOAD, time | `test_no_dropped_inp_physics`, `test_malformed_boundaries_not_coerced` |
| `Analysis.odb` | `--odb`; `residual_core/io/abaqus_odb_export.py --frames all --strict yes`: U/RF/CF/S/SDV1 | `test_conversion_rejects_mismatches`; separate genuine ODB proof |
| `sensitivity_request.json` | `--request`; `validate_request`, `scalar_results`, `reduce_scalar` | `test_selection_after_total_history`, `test_invalid_request_scope`, `test_scalar_reductions` |
| `sensitivity_results.json` | the original request, resolved scope, metadata, scalar values and derivatives | `test_compiled_offline_conversion_and_all_outputs` |
| `sensitivity_tables.csv` | output, field, component, reduction, domain, increment, time, value, parameter, derivative | same test; checks 16 rows and the derivative |
| `run_report.txt` | executed versus verified, residual and solver checks, tolerances, limitations; failure diagnostic | malformed-input CLI test and the genuine source-denied run |

The application code is in `residual_core/replay/presentation.py`, the input
conversion in `presentation_inputs.py` and the scalar mathematics in
`request_reductions.py`. Full fields, `K`, `R` and derivatives, the generated
replay record, the link shim and library, and the export command and log stay
in `private/`. The public JSON contains no full arrays.

On failure, `run_report.txt` gives a fixed category and action and the
relative path of `private/error_report.txt`, which keeps the original exception
and traceback. Missing input paths are identified by their role without
echoing the path text.

## Provider and mapping

The material developer builds the provider with UMAT-OTI
(`umat-oti-provider build CONTRACT --out DIR`, implemented in
`umat_oti/provider/build.py`) and shares the compiled object and its completed
`Mapping.json`; the source stays private. The request consumes only the OTI
object, which bundles the ORIGINAL routine and the OTI symbols in binary form.
The ORIGINAL routine supplies the virgin elastic tangent; it is never used for
source access or a production rerun.

The mapping follows the contract `resasm_umat_oti_contract_v1`: dimensions,
symbols, PROPS indices, OTI directions, derivative and Voigt layouts, the
source fingerprint and the object SHA-256. Use the generated completed
sidecar, not the input transformation contract.

Discovery checks `<object-stem>.json` and an adjacent `Mapping.json`.
`Mapping.json` is an alias of the **unchanged generated sidecar**, not a new
schema. Renaming the object requires a matching `object.sha256_full`; the
source fingerprint, layouts, dimensions and parameter order and directions are
checked too. Conflicting sidecars require `--mapping PATH`. Provenance hashes
do not certify an untrusted binary: run only providers you trust.

## Request format

See `examples/presentation_request/sensitivity_request.json`. The top-level
keys are exactly `outputs`, `parameters`, `domain` and `increments`.

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

- **Parameters:** a nonempty unique subset of `E`, `nu`, `SIGY0`, `H`.
- **Fields:** `U`, `RF`, `S`, `SDV`, with a one-based component or `ALL`. Stress
  order is 11, 22, 33, 12, 13, 23; SDV1 is the equivalent plastic strain.
- **Domain:** `nodes` for U and RF, `elements` for S and SDV; `ALL` or explicit
  existing integer ids. Each output may override `domain`. All eight
  integration points are included; there is no implicit node/element
  conversion and no integration-point subset selector.
- **Increments:** `ALL`, `LAST`, or a unique list such as `[2,4]`, one-based
  after the virgin increment 0. **The entire preceding history is replayed
  before selection.**
- **Reductions:** `component` requires one scalar location; `sum` and `mean`
  are unweighted; `L2` is the Euclidean norm (not von Mises); `max` is the
  signed maximum. The L2 derivative is `values @ derivatives / norm`, and a
  zero norm fails. Max uses the unique active derivative; ties and near ties
  fail (`rtol=1e-12, atol=1e-14`). Mean is not volume-weighted and sum is not a
  spatial integral.

Unknown fields or options, invalid ids, duplicate names, parameters or
increments, unsupported reductions and nondifferentiable cases fail; they never
disappear. Von Mises outputs, volume-weighted means, element sets, weighted
shares and full fields are available through the history engine
([REPLAY_HISTORY.md](REPLAY_HISTORY.md)). Mixed-output arithmetic,
shape, load and boundary sensitivities, and higher derivatives are not
supported.

Ordinary execution never runs finite differences and reports
`verified=false`. `--validate` requests an independent whole-history finite
difference of the ORIGINAL routine and launches no Abaqus production job. How
it behaves on single-precision ODB data, and the alternative check with the
history engine, are described in
[Example 3](../examples/presentation_request/WALKTHROUGH.md#common-problems).

## Bounded engine: physics and precision

Supported: one named `*Static` step with NLGEOM=NO; one untransformed
instance; C3D8 with B-bar; the homogeneous pinned `m3_j2` material (four
constants, one SDV); virgin state; zero fixed displacements;
parameter-independent ramped nodal concentrated loads. The ODB must save all
nodes and elements and every increment through the step end, including frame
0, with U/RF/CF/S/SDV1; CF verifies the actual loads against the input deck.
Mesh ids, coordinates and connectivity, time, loads, boundaries, integration
point stress and state, reactions and free equilibrium are checked. Missing
frames are not interpolated.

The engine carries total solution-history derivatives: `du_dp (ndof,4)`,
`dsigma_dp (nelem,8,6,4)` and `dstate_dp (nelem,8,1,4)`. Reactions use
`Rp_history + K @ du_dp`. It is a dense bounded engine, not a generic UMAT or
scalability claim; the history engine is the scalable one.

The float32 displacements of the ODB are replayed, not re-equilibrated.
Acceptance: scaled free residual below `1e-5`; stress and reaction
`rtol=2e-5` with absolute tolerance `2e-5 * max(max_abs_recorded_field, 1)`;
state `rtol=2e-5, atol=1e-8`. Zero-boundary residues below
`1e-12 * max_mesh_extent` are checked, then set to zero. All tolerances are
public. In the genuine reference, displacement rounding produces a transverse
stress error of `1.23e-4`, against an expected amplification
`eps32 * E * max(strain)` of about `6e-4`. The synthetic default replay keeps a
`1e-14` equilibrium tolerance. Unit-dependent absolute floors and ODB precision
bound the accuracy; inputs outside the gates fail without convergence repair.

## GUI

`streamlit run scripts/app.py` in the same environment opens the GUI.
**Sensitivity Request** is the first screen: the four file inputs (upload or
local path), parameter ticks read from `Mapping.json`, an output and a region
(written to `sensitivity_request.json` when no request file is given), an
optional mapping selector under Advanced, a **Solve** button, a full-field
table after the run, and downloads of the three public files. It calls the
same CLI function as the command line (`cli.main`, then `route_request`). The
older replay functions remain under **Advanced Replay**. Screen by screen:
[GUI_GUIDE.md](GUI_GUIDE.md).

## Reproduction

`scripts/reproduce_presentation_request.py` reproduces the genuine reference
case of `examples/presentation_request/Analysis.inp` in two steps:

- `prepare --work DIR` is **developer-only**: it builds the provider and runs
  one licensed Abaqus job. It needs `abaqus` and `ifort` on `PATH` (source the
  Intel compiler environment first) and a new work directory in the scratch
  location the script checks (a directory elsewhere is refused, and the error
  names the expected location). Never run it for ordinary use.
- `consume --work DIR --out NEW` runs the analysis owner's side on the five
  shared files and denies reads of the material source in that process. Only the
  generated ABI shim is exempt. This audit guard is not an operating-system
  sandbox, and no source is deleted.

Small genuine exported fields, not an ODB, are kept in
`tests/fixtures/presentation_j2`. `scripts/check_presentation_browser.py`
tests execution and downloads in a real browser and stops the exact server
process it started. The exact commands, counts and measured limits are in the
[reference-case evidence](evidence/recovery_presentation.md).
