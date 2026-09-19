# Contributing

This page is for anyone who wants to change Residual_Assembler: how to set up a
development checkout, the testing rules every change must follow, and what to
run before opening a pull request. Users who only want to run the program
should start with [START_HERE.md](START_HERE.md) and
[docs/INSTALL.md](docs/INSTALL.md).

## Setting up

```bash
git clone https://github.com/AMMS-Lab-UTSA/Residual_Assembler.git
cd Residual_Assembler
python -m venv .venv && . .venv/bin/activate
pip install -e ".[gui,yaml,test]"
./scripts/init_permissive_sources.sh
python -m pytest -q -m "not abaqus and not arc and not network"
```

`init_permissive_sources.sh` fetches only `sources/permissive/` (MIT and
BSD-3) and refuses any path outside that tier. Without it, eight tests skip
with a message naming the file, the owning submodule, its tier and the command
that fetches it.

The connected workflow tests (`resasm request`, `resasm history`,
`resasm replay`) also need the companion
[UMAT-OTI](https://github.com/AMMS-Lab-UTSA/UMAT_source_transformation)
installed beside this checkout, `gfortran`, and a genuine OTILib build for the
OTI tests; [docs/INSTALL.md](docs/INSTALL.md) and
[docs/OTILIB_VENV.md](docs/OTILIB_VENV.md) give the exact steps, and
[docs/COMPATIBILITY.md](docs/COMPATIBILITY.md) the pinned companion version.

## The rule that matters most

**A misconfigured checkout must go red, not green-but-smaller.**

Skips are how a test suite quietly stops testing things. If an external source
is missing, the affected tests skip and the suite still passes, which looks
identical to a healthy run. Set `REQUIRE_EXTERNAL_TEST_SOURCES=1` to turn those
skips into failures; CI does exactly that. Likewise, set `RUN_OTILIB_TESTS=1`
when OTILib is meant to be present, so that its tests fail rather than skip if
it is not.

The same principle governs Abaqus. A missing Abaqus installation must skip or
block only `abaqus`-marked tests. It must never shrink the offline suite.

## Test markers

Mark every new test with what it needs:

| Marker | Means |
|---|---|
| `unit` | pure Python, milliseconds |
| `integration` | several components, still offline and redistributable |
| `regression` | pins a behaviour that was once wrong; the docstring says which |
| `slow` | minutes rather than milliseconds |
| `network` | reaches the public internet |
| `fortran` | needs a Fortran compiler on PATH |
| `abaqus` | needs a licensed Abaqus; cannot run in public CI |
| `arc` | needs the documented HPC environment |
| `publication` | regenerates or checks a published artefact |
| `gui` | drives a GUI workflow in a headless browser; runs only when selected with `-m gui` |

The offline suite is `-m "not abaqus and not arc and not network"`.

## Licensing

Do not vendor external source. External Abaqus UMAT/UEL sources are pinned
submodules organised by licence tier; see
[sources/LICENSING.md](sources/LICENSING.md) and
[sources/SUBMODULES.md](sources/SUBMODULES.md). Never include, compile or link
anything from `sources/copyleft/` or `sources/license-unknown/` into the
redistributable core.

This repository is licensed **GPL-3.0-only**. By contributing you agree that
your contribution is licensed under those terms. See [LICENSE](LICENSE).

## Before opening a pull request

```bash
python -m pytest -q -m "not abaqus and not arc and not network"
python tools/audit_repository_standards.py
REQUIRE_EXTERNAL_TEST_SOURCES=1 python -m pytest -q -m "not abaqus and not arc and not network"
```

`audit_repository_standards.py` checks the required project files, that no
build products, caches, secrets or absolute home-directory paths are tracked,
and the packaging metadata. CI runs the suite in strict mode, checks that
every remaining skip names a declared external dependency (OTILib, the private
corpus, the third-party flow-rule sources, a browser for `-m gui`), proves that
strict mode fails when the sources are absent, and runs the audit.

Explain *why* in the commit message, not just what. A commit that fixes a
defect should say what the defect would have caused.

## Reporting problems

Open an issue with the input files (or the recipe), the solution data and the
command you ran. A wrong residual or sensitivity is a real finding and we want
to know about it.
