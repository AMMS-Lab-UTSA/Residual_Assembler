# Residual_Assembler: usage report

This page describes what Residual_Assembler does, how to install and use it,
and how its results are checked. It is the page to read first as a new user
or as a reviewer. The other guides give the details:
[INSTALL.md](INSTALL.md), [CLI_GUIDE.md](CLI_GUIDE.md),
[GUI_GUIDE.md](GUI_GUIDE.md) and the seven [worked examples](../examples/README.md).

Every number below was measured on 2026-09-18 on Linux (Ubuntu 20.04), Python
3.11.7, GNU Fortran 9.4.0, with Abaqus 2021.HF5 used only to read output
databases. Every command runs from the root of the Residual_Assembler
checkout, with the companion repository beside it and the environment of
[INSTALL.md](INSTALL.md) active. Outputs go to a folder outside the
repository:

```bash
export WORK="$HOME/resasm_work"; mkdir -p "$WORK"
```

Contents:

1. [What it does](#1-what-it-does)
2. [What works now, and the known limits](#2-what-works-now-and-the-known-limits)
3. [Installation](#3-installation)
4. [The command line](#4-the-command-line)
5. [Worked examples](#5-worked-examples)
6. [The GUI](#6-the-gui)
7. [The connected workflow with UMAT-OTI](#7-the-connected-workflow-with-umat-oti)
8. [Full-size models: history replay](#8-full-size-models-history-replay)
9. [How results are verified](#9-how-results-are-verified)
10. [The current transform generation](#10-the-current-transform-generation)
11. [The clean-install gate](#11-the-clean-install-gate)
12. [Status and evidence](#12-status-and-evidence)

## 1. What it does

Residual_Assembler computes **parameter sensitivities of a converged
finite-element analysis without re-running it**. A commercial solver such as
Abaqus reports displacements and stresses, but not its global residual `R`.
Residual_Assembler rebuilds `R` from the saved analysis, the material and the
loads, and solves the sensitivity equation

    K du/dp = -dR/dp

at every increment, where `K` is the tangent and `p` the material
parameters. The material enters as a compiled *provider*: the ORIGINAL UMAT
and its OTI (order-truncated imaginary) version in one object, built by the
companion package [UMAT-OTI](https://github.com/AMMS-Lab-UTSA/UMAT_source_transformation).
OTI arithmetic carries exact derivatives through the unchanged material
logic, so no derivative is written by hand. The material's source code is
not needed by the person computing the sensitivities.

The same package also assembles residuals from their ingredients (mesh,
stress or material, loads and constraints), and computes sensitivities of
residuals you supply yourself, in Python or through your own executable.

## 2. What works now, and the known limits

| Route | Commands | What it does | Verified by |
| --- | --- | --- | --- |
| Sensitivities of a finished Abaqus analysis | `resasm request`, `resasm history` | total-history first derivatives of displacements, reactions, stresses, state variables and von Mises stress with respect to the material parameters, for small-strain C3D8 analyses and any provider built by UMAT-OTI; prescribed displacements or concentrated loads; many increments; full fields and weighted shares | closed forms, replay against the ODB, whole-model finite differences of the ORIGINAL UMAT, Abaqus finite differences, a homogeneity identity (Examples 3 to 6) |
| Residual assembly from a stress field | `resasm assemble --mode stress-driven`, `resasm verify` | the global residual of a C3D8 model from supplied integration-point stresses | analytic face tractions (Example 2); patch tests in the offline suite |
| Finite-strain assembly and sensitivities | `resasm assemble --mode material-replay`, `resasm sensitivity` | compressible neo-Hookean C3D8 with the exact tangent; first-order material sensitivities with OTILib | independent quadrature, tangent and solution finite differences (Example 7) |
| Your own residual | `resasm init`, `check`, `run`, `report` | derivatives of any order of a residual written in Python (OTILib) or returned by your own executable | closed forms (Example 1) |

**Known limits.** Features outside these limits are refused with a named
reason rather than approximated.

- **Material replay.** Small strain (NLGEOM=NO), C3D8 with Abaqus's default
  selective-reduced integration, one `*Static` step, one user material on
  every element, a virgin initial state.
- **Loads and boundaries.** Concentrated loads and prescribed displacements,
  ramped over the step. Distributed and body loads, amplitudes, contact,
  several steps, materials or instances, and nonzero initial state are
  refused.
- **Derivatives.** First derivatives with respect to material parameters in
  the replay. There are no load, boundary-condition or shape sensitivities.
- **Precision.** An ODB stores single-precision fields. The replay checks
  itself against those limits, and `resasm history --reequilibrate` removes
  the equilibrium part of the rounding
  ([REPLAY_HISTORY.md](REPLAY_HISTORY.md)).
- **Assembly recipe.** The `resasm.yml` assembly recipe for C3D8 reports
  `OTI-differentiate R: NO`. Sensitivities of C3D8 models with a UMAT go
  through the compiled provider (`request`, `history`).
- **Other elements.** Only C3D8 is supported among Abaqus solids. The truss,
  beam and spring backends are small reference elements.
- **Platform.** Verified on Linux only; on Windows use WSL
  ([INSTALL.md](INSTALL.md#9-windows-use-wsl)).
- **External tools.** Abaqus is needed only to run analyses and to read
  `.odb` files. OTILib is needed only for the direct-residual routes.

## 3. Installation

On Linux with Python 3.10 or newer and `gfortran`, clone both repositories
side by side and install them into one virtual environment:

```bash
git clone https://github.com/AMMS-Lab-UTSA/Residual_Assembler.git
git clone https://github.com/AMMS-Lab-UTSA/UMAT_source_transformation.git
python3 -m venv .venv && . .venv/bin/activate
pip install -e "./Residual_Assembler[gui,yaml,test]" -e "./UMAT_source_transformation[test]"
pip check
```

Measured in a new environment: the installation took 40 s and `pip check`
found no broken requirements. `resasm` and `umat-oti-provider` are then on
`PATH`. OTILib, needed only for the direct-residual routes, is built from
source; its build took about 9 minutes. [INSTALL.md](INSTALL.md) covers:

- requirements;
- the OTILib build (never `pip install pyoti`, an unrelated package);
- the permissive test sources (`./scripts/init_permissive_sources.sh`);
- how to verify an installation;
- the clean-install gate;
- troubleshooting and WSL.

## 4. The command line

All workflows are subcommands of `resasm`:

| Group | Subcommands |
| --- | --- |
| Sensitivities of a finished analysis | `request` (the four-file interface), `history` (the history engine), `replay` (the bounded J2 replay) |
| Inspection | `inspect`, `inspect-model`, `requirements`, `doctor`, `modes`, `backends`, `template` |
| Assembly | `assemble`, `verify`, `sensitivity` |
| Jobs (`resasm.yml`) | `init`, `init-assembly`, `check`, `run`, `report` |

The companion command `umat-oti-provider build` compiles a material provider.
The exit codes follow one pattern:

- `0`: success;
- `1`: the command ran and the answer is negative;
- `2`: the command could not run with these inputs;
- `3`: OTILib was requested and is missing.

Every subcommand, with options, a worked invocation, its real output and the
`sensitivity_request.json` format, is in [CLI_GUIDE.md](CLI_GUIDE.md).

## 5. Worked examples

Seven examples, each with a walkthrough of the commands, the GUI steps, the
files written, the measured output and an independent check
([examples/README.md](../examples/README.md)):

| # | Example | Needs Abaqus? | Measured headline result |
| --- | --- | --- | --- |
| 1 | [The smallest sensitivity](../examples/user_config_minimal/WALKTHROUGH.md): a cubic spring written in Python | No (OTILib) | `du/dk = -1/3`, `du/df = 1/24` and all three second derivatives equal the closed form to 2.8e-17 or better |
| 2 | [Stress-driven C3D8 assembly](../residual_core/examples/minimal_c3d8_stress_driven/WALKTHROUGH.md) | No | `||R|| = 70.71067811865474`; the difference from the analytic face tractions is 3.6e-15 |
| 3 | [The four-file request on one element](../examples/presentation_request/WALKTHROUGH.md) | Yes (the ODB) | `dU1/dE`, `dU1/dSIGY0` and `dU1/dH` match the uniaxial J2 closed form to 1.7e-7, 8.7e-16 and 8.8e-9, relative; to about 1e-14 after re-equilibration |
| 4 | [History replay of a J2 beam, offline](../examples/replay_history/WALKTHROUGH.md) | No | all 160 Abaqus finite-difference comparisons within their uncertainty; homogeneity identity 9.4e-15; whole-model finite differences 1.5e-7 |
| 5 | [Full-size cantilevers, J2 and FCC](../examples/cantilevers/WALKTHROUGH.md) | Yes, once | homogeneity identity 1.4e-12 (J2) and 1.2e-13 (FCC) at every increment; full-size finite differences for `SIGY0` 7.9e-9 or better where resolved |
| 6 | [Provider-to-sensitivity pipeline](../examples/bounded_j2_c3d8/WALKTHROUGH.md) | No | every derivative within 9.7e-7 of whole-model finite differences (tolerance 2e-6), in about 10 s |
| 7 | [Finite-strain neo-Hookean C3D8](../examples/finite_strain_c3d8/WALKTHROUGH.md) | No (OTILib) | residual 5.7e-16 from an independent quadrature; sensitivities 1.0e-10 from nonlinear re-solves |

Examples 2, 4 and 6 need neither Abaqus nor OTILib. Example 6 checks both
packages together:

```bash
python scripts/reproduce_connected_pipeline.py --skip-abaqus --out "$WORK/pipeline"
```

Expected: `verified bounded J2 pipeline: .../manifest.json`, with
`"passed": true` in that file (measured 9.9 s in a new environment).

## 6. The GUI

```bash
streamlit run scripts/app.py
```

A Streamlit web application opens on **Sensitivity Request**, the GUI form of
`resasm request`. On that screen you:

1. give the object, the deck and the ODB (path or upload);
2. tick the parameters, which are read from `Mapping.json`;
3. choose an output, a region, a summary and the increments, or give a request
   file;
4. press **Solve**;
5. download the three public files and look at the results in the chosen
   region ([what appears after Solve](GUI_GUIDE.md#5-what-appears-after-solve)).

The other tabs (**Start here**, **1. Model** to **6. Backends**, **Advanced
Replay**) are the assembly console. Every button runs the real `resasm`
command in the same process and shows its exit code and output, so the GUI
cannot show a number the command line would not produce.

Each Solve needs a new or empty output folder. A second Solve into the same
folder is refused (`Category: output_directory`) so that results of two runs
are never mixed. Measured Solve times: 1.0 s on the one-element analysis and
24.5 s on the full-size J2 cantilever. Screen by screen, with screenshots:
[GUI_GUIDE.md](GUI_GUIDE.md).

## 7. The connected workflow with UMAT-OTI

The two packages divide the work between the owner of a material and the
person who ran the analysis:

1. **Material owner.** Builds the provider from the UMAT and its contract with
   `umat-oti-provider build`, which writes the object and its completed
   mapping. The owner hands over only those two files (renaming them
   `OTI_UMAT.obj` and `Mapping.json` is allowed). The mapping records the parameter names,
   their PROPS positions, the layouts and the object's SHA-256, which the
   request checks.
2. **Analysis owner.** Runs `resasm request` with the deck, the ODB, the
   object and a request. The material source is never read. The three public
   files `sensitivity_results.json`, `sensitivity_tables.csv` and
   `run_report.txt` are written at the top of the output folder. Full fields,
   the ODB export and the link library stay in `private/`.

From the repository root, with the one-element example (Example 3):

```bash
umat-oti-provider build ../UMAT_source_transformation/parameter_sensitivity/models/m3_j2/contract_v2.json \
    --out "$WORK/provider_j2"
resasm request --model examples/presentation_request/Analysis.inp --odb /path/to/Analysis.odb \
    --material "$WORK/provider_j2/umat_m3_j2_oti.obj" \
    --request examples/presentation_request/sensitivity_request.json --out "$WORK/one_element_results"
```

The ODB is not shipped. It is the Abaqus result of
`examples/presentation_request/Analysis.inp`, and reading it needs Abaqus
Python (`abaqus` on `PATH`, or `--abaqus`). Measured: the build took 5.5 s;
the request, on an Abaqus 2021.HF5 ODB of this deck, took 1.1 s including
the export and printed `request executed: 4 scalar results; verified=False`.
`verified=False` means that no independent check was requested; the checks
are in
[section 9](#9-how-results-are-verified).

`resasm request` solves models like this one with its bounded
single-material engine. It hands models outside that engine's scope to the
history engine (`resasm history`), and its first line of output names the
reason: prescribed displacements, deck options the bounded reader does not
accept, a provider other than the pinned m3_j2, sets, von Mises outputs and
the other request extensions. The exact rule is in
[REPLAY_HISTORY.md](REPLAY_HISTORY.md#how-resasm-request-chooses-this-engine).
The request format is specified in
[CLI_GUIDE.md](CLI_GUIDE.md#the-request-file). The interface is described in
detail in [REQUEST_INTERFACE.md](REQUEST_INTERFACE.md).

## 8. Full-size models: history replay

`resasm history` replays the whole recorded history of an analysis with any
provider built by UMAT-OTI. It assembles the residual sparsely and factorises
the tangent once per increment for all parameters. It accepts the ODB, or an
export of it, which then needs no Abaqus. The committed example runs anywhere:

```bash
resasm history --model examples/replay_history/j2_beam/Analysis.inp \
    --fields examples/replay_history/j2_beam/fields.npz \
    --material "$WORK/provider_j2/umat_m3_j2_oti.obj" \
    --request examples/replay_history/j2_beam/sensitivity_request.json --out "$WORK/beam"
```

Measured (1.0 s):

```text
history replay: 10 increments, 768 integration points, 4 parameters; yes: max|R_free| = 1.288e-03 N at increment 10 (limit 1.609e-01 N there)
```

Two full-size cantilevers, from the ODBs of the decks written by
`examples/cantilevers/gen_cantilever.py`
([Example 5](../examples/cantilevers/WALKTHROUGH.md)):

| | J2 cantilever | FCC crystal-plasticity cantilever |
| --- | --- | --- |
| mesh, DOF, integration points | 1,536 C3D8, 7,497, 12,288 | 384 C3D8, 2,025, 3,072 |
| increments, parameters | 40, 4 | 25, 10 |
| `resasm request` from the ODB (includes the export) | 22 to 27 s | 24 to 25 s |
| `resasm history --reequilibrate` from the export | 33 to 40 s | 71 to 79 s |
| replayed stress, state and reactions against the ODB, worst error/limit | 0.010, 0.002, 0.006 | 0.006, 0.002, 0.006 |
| homogeneity identity at every increment, re-equilibrated | 1.4e-12 | 1.2e-13 |

The times were measured on a 24-core workstation while other jobs were
running on it. They include writing the full-field `fields.npz`, which is
110 to 135 MB. The J2 weighted shares of the von Mises field change as
plasticity spreads: `E` holds 98.21 % while the beam is elastic (increments 1
to 6), and `SIGY0` holds 89.05 % at increment 40. The mathematics, the
supported deck subset and the tolerances are in
[REPLAY_HISTORY.md](REPLAY_HISTORY.md).

## 9. How results are verified

An ordinary run computes derivatives but does not claim to have verified
them. `run_report.txt` states each verdict separately:

- command executed;
- residual assembled;
- equilibrium checked and passed;
- tangent available and verified;
- derivative calculated and verified;
- reference resolved;
- Abaqus comparison available;
- unsupported feature detected;
- public and private outputs separated.

Only an independent reference turns *calculated* into *verified*. The
references used across the examples:

| Reference | What it checks | Where |
| --- | --- | --- |
| Closed forms | the whole chain on problems with a known answer | Examples 1, 2, 3 |
| Replay against the ODB | at every integration point and increment, the replayed stress, state and reactions are compared with the recorded values; any excess over the single-precision limits stops the run | every `request`/`history` run |
| Whole-model finite differences of the ORIGINAL UMAT | `resasm history --verify fd` re-solves the model with the unmodified routine at `p (1 +/- h)` over a ladder of steps | Examples 3, 4, 5, 6 |
| Abaqus finite differences | perturbed Abaqus reruns | Examples 4, 5 |
| Homogeneity identity | J2 with linear hardening, and the FCC crystal, are homogeneous of degree one in their stress-dimensioned parameters, so `sum_p p dQ/dp` equals `Q` for stresses and reactions and zero for displacements, at every increment; the engine does not use this | Examples 4, 5 |

For example, adding the finite-difference check to the committed beam
(measured 11.1 s):

```bash
resasm history --model examples/replay_history/j2_beam/Analysis.inp \
    --fields examples/replay_history/j2_beam/fields.npz \
    --material "$WORK/provider_j2/umat_m3_j2_oti.obj" \
    --request examples/replay_history/j2_beam/sensitivity_request.json \
    --out "$WORK/beam_checked" --reequilibrate --verify fd
```

Its `run_report.txt` then reads `Tangent verified: yes: max relative error
1.88e-10 ...` and `Derivative verified: yes: whole-model central FD of the
ORIGINAL UMAT re-equilibrated in Python; worst nonzero-derivative error
1.46e-07 ...`.

[VERIFICATION_RECORD.md](VERIFICATION_RECORD.md) lists every quantitative
claim with the command that reproduces it, its independent reference and the
measured value.

## 10. The current transform generation

Evidence about a transformed material belongs to the transform code that
produced it. The current transform generation is **`dbe9f928191e1d43`**.
It is recorded once, in `schemas/transform_generation.json`, a file that is
identical in both repositories. The two regression fixtures in
`tests/fixtures/verified/` (isotropic elasticity and J2) were regenerated in
Abaqus at this generation. Their tangents agree over a step-size plateau to
1.1e-14 and 8.3e-11. Older fixtures, under `tests/fixtures/historical/`, are
kept as history and refused as regression baselines. How the re-freeze was
done: [evidence/final_refreeze.md](evidence/final_refreeze.md).

To check that the recorded generation is the one the installed UMAT-OTI
computes, and that the fixtures carry it:

```bash
python -c "import json; print(json.load(open('schemas/transform_generation.json'))['transform_fingerprint'])"
python -c "from umat_oti.store import transform_fingerprint; print(transform_fingerprint())"
python -m pytest -q tests/contract/test_the_two_repositories_speak_one_contract.py \
    tests/framework/test_a_fixture_is_held_to_the_rule_that_froze_it.py
```

Measured: both commands print `dbe9f928191e1d43`, and the tests report
`26 passed`.

## 11. The clean-install gate

`scripts/clean_install_gate.py` is the acceptance test of an installation
([INSTALL.md](INSTALL.md#7-the-clean-install-gate)). It builds wheels of both
repositories from clean trees and installs them into a new environment. From
the installed commands only, it then runs:

- the provider build;
- the four-file request on a genuine ODB, against the uniaxial closed form;
- the same request with every read of a Fortran source denied;
- both GUIs up to HTTP readiness;
- with `--cantilever`, the full-size J2 cantilever with the homogeneity
  identity.

The recorded run is in
[evidence/final_clean_clone.md](evidence/final_clean_clone.md). It ran on
fresh clones of the published branches, at Residual_Assembler `3504a02` and
UMAT_source_transformation `1352114`:

- all 22 commands exited 0 and the gate passed;
- the request matched the closed form to 1.7e-7 (`E`), 8.7e-16 (`SIGY0`) and
  8.8e-9 (`H`);
- the source-denied outputs were byte-identical;
- the cantilever request took 22.7 s and the re-equilibrated replay 33.0 s,
  with the identity at 1.4e-12.

From the same clones, the UMAT-OTI suite passed (3,370 passed, 158 skipped).
The offline Residual_Assembler suite had 494 passed, 20 skipped and 6 failed,
all six from one check in the verification scripts, which has since been
corrected. That record describes an **earlier commit pair**. A new gate run
on the published commits follows the push of this version.

## 12. Status and evidence

The requirement-by-requirement status is in
[COMPLETION_LEDGER.md](COMPLETION_LEDGER.md); the evidence files are in
[evidence/](evidence/).
