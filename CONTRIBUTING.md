# Contributing

## Setting up

```bash
git clone https://github.com/AMMS-Lab-UTSA/Residual_Assembler.git
cd Residual_Assembler
python -m venv .venv && . .venv/bin/activate
pip install -e ".[yaml]" pytest
./scripts/init_permissive_sources.sh
python -m pytest -q
```

`init_permissive_sources.sh` fetches only `sources/permissive/` (MIT / BSD-3).
It refuses any path outside that tier. Without it, eight tests skip with a
message naming the file, the owning submodule, its tier and the command that
fetches it.

## The rule that matters most

**A misconfigured checkout must go red, not green-but-smaller.**

Skips are how a test suite quietly stops testing things. If an external source
is missing, the affected tests skip and the suite still passes, which looks
identical to a healthy run. Set `REQUIRE_EXTERNAL_TEST_SOURCES=1` to turn those
skips into failures; CI does exactly that.

The same principle governs Abaqus. A missing Abaqus installation must skip or
block only `abaqus`-marked tests. It must never shrink the offline suite.

## Test markers

Mark every new test with what it needs:

| Marker | Means |
|---|---|
| `unit` | pure Python, milliseconds |
| `integration` | several components, still offline and redistributable |
| `slow` | minutes rather than milliseconds |
| `network` | reaches the public internet |
| `fortran` | needs a Fortran compiler on PATH |
| `abaqus` | needs a licensed Abaqus; cannot run in public CI |
| `arc` | needs the documented HPC environment |
| `publication` | regenerates or checks a published artefact |

The offline suite is `-m "not abaqus and not arc and not network"`.

## Licensing

Do not vendor external source. External Abaqus UMAT/UEL sources are pinned
submodules organised by licence tier; see
[`sources/LICENSING.md`](sources/LICENSING.md). Never `include`, compile or link
anything from `sources/copyleft/` or `sources/license-unknown/` into the
redistributable core.

Note that **this repository does not yet declare a licence for its own code**.
Until it does, contributions cannot be redistributed under any grant. See
[`AUTHORS.md`](AUTHORS.md).

## Before opening a pull request

```bash
python -m pytest -q -m "not abaqus and not arc and not network"
python tools/audit_repository_standards.py
REQUIRE_EXTERNAL_TEST_SOURCES=1 python -m pytest -q
```

Explain *why* in the commit message, not just what. A commit that fixes a defect
should say what the defect would have caused.

## Reporting problems

Open an issue with the recipe, the solution data and the command you ran. A
wrong residual or sensitivity is a real finding and we want to know about it.
