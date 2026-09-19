# Installation guide

This page takes you from an empty machine to a working installation of
Residual_Assembler and its companion UMAT-OTI, and shows how to confirm that
everything works. The two packages install in about a minute; the optional
OTILib build takes about ten minutes more.

Contents:

1. [Requirements](#1-requirements)
2. [Get the two repositories](#2-get-the-two-repositories)
3. [Create a virtual environment and install](#3-create-a-virtual-environment-and-install)
4. [Install OTILib (for the direct-residual paths)](#4-install-otilib-for-the-direct-residual-paths)
5. [Fetch the permissive test sources](#5-fetch-the-permissive-test-sources)
6. [Verify the installation](#6-verify-the-installation)
7. [The clean-install gate](#7-the-clean-install-gate)
8. [Troubleshooting](#8-troubleshooting)
9. [Windows: use WSL](#9-windows-use-wsl)

## 1. Requirements

| Item | Needed for | Notes |
| --- | --- | --- |
| Linux (x86-64) | everything | verified on Ubuntu 20.04 (glibc 2.31). macOS and native Windows are not established; on Windows use WSL ([section 9](#9-windows-use-wsl)) |
| Python 3.10 or newer, with `venv`, `ctypes` and `ssl` | everything | Residual_Assembler alone accepts 3.9, but UMAT-OTI needs 3.10. Verified with Python 3.11.7 |
| `gfortran` | building material providers (`umat-oti-provider build`) and linking them at run time | verified with GNU Fortran 9.4.0 |
| `git` | cloning, and the pipeline script of Example 6 (it records commits) | |
| NumPy, SciPy | everything | installed automatically by `pip` |
| Streamlit | the GUI | installed by the `gui` extra |
| Abaqus (tested: 2021.HF5) | only to run analyses and to read `.odb` files (`abaqus python`) | the sensitivity computation itself never calls Abaqus; exported fields can be replayed without it |
| OTILib (GPLv3, external) | only the direct-residual paths: `resasm run` with `backend: otilib`, `resasm sensitivity`, the finite-strain example | built from source, [section 4](#4-install-otilib-for-the-direct-residual-paths). **Never** `pip install pyoti`: that PyPI name is an unrelated package |

What needs what, in short: the collaborator workflow (`resasm request`,
`resasm history`) needs Python, `gfortran` and the two packages, plus Abaqus
Python only when it has to read an `.odb`. It does not need OTILib.

## 2. Get the two repositories

Clone both repositories into the same parent folder, side by side. Several
scripts find the companion there by default.

```bash
mkdir -p ~/resasm && cd ~/resasm
git clone https://github.com/AMMS-Lab-UTSA/Residual_Assembler.git
git clone https://github.com/AMMS-Lab-UTSA/UMAT_source_transformation.git
```

The result:

```text
~/resasm/
    Residual_Assembler/            this repository
    UMAT_source_transformation/    the companion, UMAT-OTI
```

The pair of commits known to work together is recorded in
[COMPATIBILITY.md](COMPATIBILITY.md).

## 3. Create a virtual environment and install

From the parent folder:

```bash
python3 -m venv .venv
. .venv/bin/activate
pip install -e "./Residual_Assembler[gui,yaml,test]" -e "./UMAT_source_transformation[test]"
pip check
```

`-e` installs both packages in editable mode, so the commands use the code in
your checkouts and the examples and templates are found where they are.
The extras add the GUI (`gui`: Streamlit), YAML support for `resasm.yml`
(`yaml`: PyYAML; a minimal reader is bundled if you leave it out) and the test
tools (`test`: pytest and friends).

Measured on 2026-09-18 in a new virtual environment made from Python 3.11.7:
the `pip install` took 40 s and `pip check` reported
`No broken requirements found.` The two commands below are then on your `PATH`:

| Command | Package | What it does |
| --- | --- | --- |
| `resasm` | Residual_Assembler | every workflow of this repository ([CLI guide](CLI_GUIDE.md)) |
| `umat-oti-provider` | UMAT-OTI | builds a compiled material provider (`OTI_UMAT.obj` + its mapping) from a UMAT and its contract |

Activate the environment (`. .venv/bin/activate`) in every new shell.

## 4. Install OTILib (for the direct-residual paths)

Skip this section if you only need `resasm request` and `resasm history`.

OTILib is the hypercomplex (OTI) arithmetic library used by the Python
residual paths. It is GPLv3, developed separately
(<https://github.com/mauriaristi/otilib>), and never copied into this
repository. It must be built for the same Python version as your virtual
environment.

**Without Conda (recommended with a virtual environment).** With your
virtual environment active (the build uses the `python` on `PATH`), in a
folder outside both checkouts:

```bash
python -m pip install cmake==3.31.10 Cython==3.3.0
git clone https://github.com/mauriaristi/otilib.git "$HOME/otilib"
git -C "$HOME/otilib" checkout --detach a4b7a05ca275e8d441b0b717b7b272e96728ffcc
cmake -S "$HOME/otilib" -B "$HOME/otilib/build" -DBUILD_TESTING=OFF
cmake --build "$HOME/otilib/build" --target oticython --parallel 2
cmake --build "$HOME/otilib/build" --target gendata --parallel 2
```

Measured on 2026-09-18 in the fresh environment of section 3: every command
exited 0 (the compiler prints many warnings); the `oticython` build took
8 min 45 s and `gendata` 8 s. The pinned commit and tool versions are those
of the validated procedure in [OTILIB_VENV.md](OTILIB_VENV.md), which also
shows how to build in a fully isolated shell.

**With Conda.** `bash Residual_Assembler/scripts/setup_otilib.sh [TARGET_DIR]`
clones OTILib, creates a Conda environment `pyoti`, builds it with CMake and
generates its data tables. It needs `git`, `cmake` and `conda` (this script
was not re-run for this guide).

**Tell Residual_Assembler where the build is.** Point both variables to the
**build directory** (not the source folder), in every shell where you use the
OTILib paths:

```bash
export PYOTI_PATH="$HOME/otilib/build"
export OTILIB_ROOT="$HOME/otilib/build"
python -c "from residual_core.algebra.otilib_adapter import otilib_status; print(otilib_status())"
```

A working build prints a dictionary with `'available': True` and
`'api_module': 'pyoti.sparse'`. Without OTILib, the commands that need it stop
with `OTILib backend requested but genuine OTILib was not found` instead of
computing anything else; they never fall back silently.

## 5. Fetch the permissive test sources

The offline test suite and a few inspection examples read Abaqus input files
from external, pinned git submodules. Only the permissively licensed ones
(MIT, BSD-3) are fetched by this script:

```bash
cd Residual_Assembler
./scripts/init_permissive_sources.sh           # the sources the tests need
./scripts/init_permissive_sources.sh --check   # show the state, fetch nothing
```

Copyleft and licence-unknown sources are marked `update = none` and are never
fetched by setup ([sources/LICENSING.md](../sources/LICENSING.md)). None of the
eight worked examples needs these sources. Without them, the tests that read
them are skipped, each naming the file it needs and this command.

## 6. Verify the installation

Run these from the `Residual_Assembler` folder with the environment active.

**1. The commands answer.**

```bash
resasm --help
umat-oti-provider --help
```

`resasm --help` lists the 18 subcommands (`replay`, `request`, `history`,
`inspect`, `requirements`, `assemble`, `verify`, `doctor`, `template`,
`backends`, `sensitivity`, `modes`, `init`, `inspect-model`, `init-assembly`,
`check`, `run`, `report`).

**2. A first sensitivity, with no Abaqus and no OTILib.** Build the J2 material
provider and replay the committed Abaqus beam
([Example 4](../examples/replay_history/WALKTHROUGH.md)):

```bash
export WORK="$HOME/resasm_work"; mkdir -p "$WORK"
umat-oti-provider build ../UMAT_source_transformation/parameter_sensitivity/models/m3_j2/contract_v2.json \
    --out "$WORK/provider_j2"
resasm history --model examples/replay_history/j2_beam/Analysis.inp \
    --fields examples/replay_history/j2_beam/fields.npz \
    --material "$WORK/provider_j2/umat_m3_j2_oti.obj" \
    --request examples/replay_history/j2_beam/sensitivity_request.json --out "$WORK/beam"
```

Expected (measured in the fresh environment of section 3):

```text
history replay: 10 increments, 768 integration points, 4 parameters; yes: max|R_free| = 1.288e-03 N at increment 10 (limit 1.609e-01 N there)
```

and `$WORK/beam/run_report.txt` begins with `Status: executed successfully`.

**3. Both packages together, with an independent check.** This builds the
provider again, solves a cyclic history and verifies it against finite
differences ([Example 6](../examples/bounded_j2_c3d8/WALKTHROUGH.md)):

```bash
python scripts/reproduce_connected_pipeline.py --skip-abaqus --out "$WORK/pipeline"
```

Expected: `verified bounded J2 pipeline: .../manifest.json` with
`"passed": true` (measured 9.9 s in the fresh environment, without OTILib).

**4. OTILib, if you installed it.**

```bash
resasm init --template python --out "$WORK/spring"
resasm check "$WORK/spring/resasm.yml"
```

Every line must read `[ok]`, including `[ok] OTILib available: order 2, basis 2`
([Example 1](../examples/user_config_minimal/WALKTHROUGH.md)).

**5. The GUI starts.**

```bash
streamlit run scripts/app.py
```

Open the address it prints (by default <http://localhost:8501>). The first tab
is **Sensitivity Request** ([GUI guide](GUI_GUIDE.md)). Measured: in the fresh
environment the server answered its health check (`/_stcore/health` returned
`ok`).

**6. The offline test suite (optional, several minutes).**

```bash
python -m pytest -q -m "not abaqus and not arc and not network"
```

The tests that need OTILib run whenever `PYOTI_PATH` and `OTILIB_ROOT`
(section 4) point to a working build, and are skipped by name when it is
missing. Set `RUN_OTILIB_TESTS=1` whenever OTILib is meant to be present: a
missing or broken build then fails those tests instead of skipping them.

If the two packages are installed from wheels rather than with `-e` (the
clean-install gate of section 7 installs them that way), first run
`export UMAT_OTI_REPO="$PWD/../UMAT_source_transformation"`. The history-replay
tests (`tests/replay_history/`) and the verification tests (`tests/verification/`)
read the UMAT models from that checkout, and a wheel does not contain them.
With the `-e` install of section 3 they find the checkout themselves.

The full suite was not re-run for this guide. The tests that exercise the
examples were: `tests/integration/test_presentation_request.py -k "real_archived
or compiled_offline"` (2 passed) and `tests/replay_history/test_history_example.py`
with `tests/gui/test_solve_screen.py` (16 passed).

**7. Every example that needs no Abaqus, in one command (optional, 4 min).**
With OTILib (section 4):

```bash
python scripts/audit_recovery_usage.py --umat ../UMAT_source_transformation --phase examples \
    --work "$WORK/examples" --evidence-dir "$WORK/examples_record"
```

It first runs the two documented checks, `otilib_status()` (section 4; it must
report `'available': True`) and `transform_fingerprint()` (it must equal the
generation in `schemas/transform_generation.json`). Then, in this order, it
runs the commands each example's walkthrough gives, checks each result against
the example's reference, and stops at the first failure:

| Record | Example | Commands | Check |
| --- | --- | --- | --- |
| R-X1 | 1 | `resasm init --template python`, `resasm check`, `resasm run`, `resasm report` | the five derivatives against the closed form, to 1e-8 |
| R-X2 | 8 | the same with `--template blackbox-order2` | the same |
| R-X3 | 2 | `resasm assemble` of `minimal_c3d8_stress_driven` | the norm and the largest entry against the face tractions, to 1e-10 |
| Example 4 | 4 | `umat-oti-provider build` of the J2 contract, then `resasm history` on the committed beam | `run_report.txt` begins with `Status: executed successfully` |
| R-X4 | 6 | `python scripts/reproduce_connected_pipeline.py --skip-abaqus` | `"passed": true` in its manifest |
| R-X5 | 7 | `python examples/finite_strain_c3d8/benchmark.py`, `resasm --config ... assemble` and `resasm --config ... sensitivity` | `"passed": true` in the benchmark's report |
| U-X1 to U-X5 | UMAT-OTI 1 to 5 | the command blocks of their READMEs | `RESULT: PASS` from each example's script; the verifier's `"passed": true` for 3, 4 and Part A of 5 |
| Example 6 | UMAT-OTI 6 | the twenty-model sweep | every model reproduces the original's stress; no comparison row disagrees |

R-X1 to R-X5 and U-X1 to U-X5 are the ten reference examples that
[COMPLETION_LEDGER.md](COMPLETION_LEDGER.md) names. The Residual_Assembler
commands run from this folder, the UMAT-OTI commands from the `--umat` checkout,
and every output goes to `--work` (the UMAT-OTI examples write there instead of
to `umat_oti_workspace/examples/`). Each command, exit code, output and check is
recorded in `RA/usage_examples.json` and `UMAT/usage_examples.json` under
`--evidence-dir` (without it, in `docs/evidence/` of each checkout).

It sets `UMAT_OTI_REPO` to the `--umat` checkout, from which the pipeline of
Example 6 builds its provider. It puts its own Python's `bin` folder first on
`PATH`, so that `resasm`, `umat-oti` and `umat-oti-provider` are those of the
same environment, and runs every Python command as `python -I`, so that the
installed packages are imported rather than the checkouts (`--imports
environment` runs plain `python` with your `PYTHONPATH` instead). Examples 3
and 5 need Abaqus and are not part of it; the clean-install gate (section 7)
runs Example 3's request and, with `--cantilever`, Example 5's J2 model.

## 7. The clean-install gate

`scripts/clean_install_gate.py` is the acceptance test of an installation. It
does from scratch what sections 3 and 6 describe, in an isolated environment,
and writes a machine-readable report.

```bash
python Residual_Assembler/scripts/clean_install_gate.py \
    --umat-repo UMAT_source_transformation --python python3.11 \
    --odb /path/to/Analysis.odb --work /path/to/new/folder \
    [--branch main] [--abaqus abaqus] [--cantilever /path/to/cantilever/work]
```

What it checks, in order:

1. Both repositories are git checkouts with no uncommitted or untracked files
   (so the result describes a commit), on the named branch if `--branch` is
   given; it records whether each commit is the published head of that branch.
2. The given Python is 3.10 or newer with `ctypes`, `ssl` and `venv`;
   `gfortran` and the Abaqus launcher are on `PATH`.
3. It creates a new virtual environment with a scratch `HOME` and no inherited
   Python path, builds a wheel of each repository, installs both with their
   extras and runs `pip check`. It then checks that the packages
   (`residual_core`, `resasm_user`, `umat_oti` and both GUIs) are imported
   from the new environment's site-packages and not from the checkouts, that
   no installed distribution is editable, and that the packaged data files
   (among them the Fortran driver and stubs and the contract schemas) are
   present.
4. With the installed `umat-oti-provider`, it builds the J2 provider from a
   copy of the public contract and UMAT.
5. With the installed `resasm`, it runs the four-file request of
   [Example 3](../examples/presentation_request/WALKTHROUGH.md) on the ODB you
   give (`--odb`, the Abaqus result of `examples/presentation_request/Analysis.inp`),
   checks the three public files, and compares the derivatives with the
   uniaxial closed form (relative error below 2e-5, `|dU1/dnu|` below 1e-8).
6. It repeats the request with every read of a Fortran source denied and
   requires identical public files: the material source is not needed.
7. It renders both installed GUIs headlessly with Streamlit's `AppTest` (no
   exception, a title or header present), then starts each with
   `python -m streamlit run` and waits until it answers over HTTP
   (`/_stcore/health`).
8. With `--cantilever` (a folder holding `j2/cantilever_j2_nominal.inp` and
   `.odb`, see [Example 5](../examples/cantilevers/WALKTHROUGH.md)), it runs
   the request on the full-size J2 cantilever, re-equilibrates it, and
   requires the homogeneity identity at every increment to 1e-10.

**The two Abaqus inputs.** `--odb` is the Abaqus result of
`examples/presentation_request/Analysis.inp` with the ORIGINAL J2 UMAT
([Example 3](../examples/presentation_request/WALKTHROUGH.md)), and
`--cantilever` a folder whose `j2/` holds the J2 cantilever deck and its
result ([Example 5](../examples/cantilevers/README.md)). Make both once, from
the folder that holds the two checkouts, into a folder outside them:

```bash
RA="$PWD/Residual_Assembler"; UMAT="$PWD/UMAT_source_transformation"
IN="$HOME/gate_inputs"; mkdir -p "$IN/one_element" "$IN/cantilever/j2"
cp "$RA/examples/presentation_request/Analysis.inp" "$IN/one_element/"
python "$RA/examples/cantilevers/gen_cantilever.py" j2 --out "$IN/cantilever/j2/cantilever_j2_nominal.inp"
cd "$IN/one_element" && abaqus job=Analysis input=Analysis.inp \
    user="$UMAT/parameter_sensitivity/models/m3_j2/umat.for" cpus=1 interactive
cd "$IN/cantilever/j2" && abaqus job=cantilever_j2_nominal input=cantilever_j2_nominal.inp \
    user="$UMAT/parameter_sensitivity/models/m3_j2/umat.for" double=both interactive
grep "COMPLETED SUCCESSFULLY" "$IN/one_element/Analysis.sta" "$IN/cantilever/j2/cantilever_j2_nominal.sta"
```

Then give the gate `--odb "$IN/one_element/Analysis.odb" --cantilever "$IN/cantilever"`.
Judge each job by `THE ANALYSIS HAS COMPLETED SUCCESSFULLY` in its `.sta`
file: Abaqus 2021.HF5 can abort with signal 6 during teardown after writing a
complete ODB. The two Abaqus jobs were not re-run for this guide; the other
lines were (2026-09-19).

`--work` must be a new folder outside both repositories. It receives
`report.json` (every command, exit code, log path, the wheel digests, the
versions and the verdict `passed`) and the logs. The gate was not re-run for
this guide; its recorded results are in
[evidence/final_clean_clone.md](evidence/final_clean_clone.md).

**The environment the gate made.** `<work>/env` holds both packages installed
from their wheels. To run the offline test suites (section 6, item 6) and the
examples check (section 6, item 7) against exactly those wheels, activate it
with `. <work>/env/bin/activate`, set the OTILib variables of section 4,
`RUN_OTILIB_TESTS=1` and `UMAT_OTI_REPO` (section 6, item 6), and run them
from the checkouts.

**From fresh clones, all at once.** `scripts/reproduce_from_clean_clones.sh`
does all of this from fresh clones of the published `main` branches. With the
OTILib build of section 4 and the two Abaqus inputs above:

```bash
export PYOTI_PATH="$HOME/otilib/build" OTILIB_ROOT="$HOME/otilib/build"
bash Residual_Assembler/scripts/reproduce_from_clean_clones.sh --python python3.11 \
    --odb "$IN/one_element/Analysis.odb" --cantilever "$IN/cantilever" /path/to/new/folder
```

It first drops the calling shell's Python settings (`PYTHON*`, `PIP_*`,
`CONDA*`, `VIRTUAL_ENV`), so that no other Python installation reaches the
run. It then clones both repositories into the new folder (section 2), runs
the gate with `--branch main` and `--cantilever`, and, in the gate's
environment as above, runs both offline suites and the examples check. Each
step's exit code goes to `summary.tsv`, next to the suites' JUnit XML
(`ra_suite.xml`, `umat_suite.xml`), the gate's `gate/report.json` and the
examples' records (`usage/`). It stops before the gate if the OTILib build is
not at the commit section 4 pins, and records whether the run left the clones
unmodified. It needs Abaqus (the gate reads both ODBs) and takes about
30 minutes.

## 8. Troubleshooting

| Symptom | Cause and fix |
| --- | --- |
| `resasm: command not found` | The virtual environment is not active: `. .venv/bin/activate`. |
| `pip` refuses `umat-oti` because it requires a different Python | UMAT-OTI needs Python 3.10 or newer. Create the environment with such a Python (`python3.11 -m venv .venv`). |
| `umat-oti-provider build` fails, or a request stops while linking the material | `gfortran` is missing: `sudo apt-get install gfortran` (Debian, Ubuntu). |
| `OTILib backend requested but genuine OTILib was not found` | Set `PYOTI_PATH` and `OTILIB_ROOT` to the OTILib **build** directory ([section 4](#4-install-otilib-for-the-direct-residual-paths)). If it still fails, the build is for another Python version: rebuild it with the environment's Python. |
| `import pyoti` works but the adapter reports it unavailable, or behaves strangely | You installed the PyPI package `pyoti`, which is unrelated. `pip uninstall pyoti`, then build OTILib from source. |
| `request failed: Category: odb_export` (the private `error_report.txt` says `Abaqus executable 'abaqus' not found on PATH`), or `Abaqus launcher 'abaqus' not found` from `resasm history` | Reading an `.odb` needs Abaqus Python. Pass `--abaqus /path/to/abaqus`, or export once where Abaqus is installed and use `resasm history --fields` elsewhere. |
| `No module named umat_oti.provider` from `reproduce_connected_pipeline.py` | The isolated Python sees an older UMAT-OTI installation. Reinstall both packages as in section 3, or pass `--imports environment` when you work through `PYTHONPATH`. |
| Tests skip with a message naming `sources/permissive/...` | Run `./scripts/init_permissive_sources.sh` (section 5). |
| The GUI port is taken | `streamlit run scripts/app.py --server.port 8502`. |

## 9. Windows: use WSL

The compiled material providers and OTILib are built with the GNU toolchain
and are verified on Linux only. On Windows, install WSL 2 with an Ubuntu
distribution and follow this guide inside it:

```powershell
wsl --install -d Ubuntu
```

Then, in the Ubuntu shell, install the tools (`sudo apt-get install python3
python3-venv gfortran git`), check that `python3 --version` reports 3.10 or
newer (install a newer Python if it does not), and continue from
[section 2](#2-get-the-two-repositories).
Keep the repositories and outputs inside the Linux file system (for example
under `~/resasm`), not under `/mnt/c`, for speed. WSL 2 forwards local ports,
so the GUI started inside WSL normally opens in the Windows browser at the
address it prints.
`scripts/run_otilib_tests_wsl.sh` (and `scripts/run_otilib_tests_wsl.ps1`,
started from Windows) run the OTILib tests inside WSL and stop with a
message naming anything that is missing. Abaqus itself runs on Windows; export
the `.odb` there with `abaqus python residual_core/replay/odb_export_npz.py --
Analysis.odb fields.npz` and replay the export inside WSL with
`resasm history --fields`.

Next: the [worked examples](../examples/README.md), the
[CLI guide](CLI_GUIDE.md) and the [GUI guide](GUI_GUIDE.md).
