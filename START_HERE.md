# Start here

This page is the entry point for anyone opening Residual_Assembler for the
first time: what the program does, how to get it running, and which document
to read next.

## What Residual_Assembler does

Residual_Assembler computes **parameter sensitivities of a converged
finite-element analysis without re-running it**. You give it four files:

| Input | What it is |
|---|---|
| `Analysis.inp` | the Abaqus input deck of the finished analysis |
| `Analysis.odb` | its output database |
| `OTI_UMAT.obj` + `Mapping.json` | the material, compiled as an OTI provider by the companion [UMAT-OTI](https://github.com/AMMS-Lab-UTSA/UMAT_source_transformation) |
| `sensitivity_request.json` | which outputs, parameters, region and increments you want |

It replays the material at every integration point and increment, assembles
the residual `R`, its tangent `K` and `dR/dp`, solves `K du/dp = -dR/dp`, and
writes the requested sensitivities:

```bash
resasm request --model Analysis.inp --odb Analysis.odb \
    --material OTI_UMAT.obj --request sensitivity_request.json --out results
```

The material's source code is not needed to run this, and nothing is sent over
the network.

The same package also assembles residuals from their ingredients (mesh,
material, solution fields, loads and constraints) and computes sensitivities
of residuals you supply yourself, directly or through a private black-box
executable. The ready-made templates for those routes are in
[`templates/`](templates/).

## Get it running

[docs/INSTALL.md](docs/INSTALL.md) is the installation guide. In short, on
Linux with Python 3.10 or newer and `gfortran`:

```bash
git clone https://github.com/AMMS-Lab-UTSA/Residual_Assembler.git
git clone https://github.com/AMMS-Lab-UTSA/UMAT_source_transformation.git
python3 -m venv .venv && . .venv/bin/activate
pip install -e "./Residual_Assembler[gui,yaml,test]" "./UMAT_source_transformation[test]"
resasm --help
```

Abaqus (tested: 2021.HF5) is needed only to export an ODB; the sensitivity
computation itself does not call it.

To run the offline test suite from a checkout:

```bash
cd Residual_Assembler
./scripts/init_permissive_sources.sh
python -m pytest -q -m "not abaqus and not arc and not network"
```

`init_permissive_sources.sh` fetches only the MIT/BSD-3 tier of the external
Abaqus sources that eight tests read. Copyleft and licence-unknown tiers are
marked `update = none` and are never fetched by ordinary setup. If you skip
that step, those eight tests skip, each naming the file it needs, the
submodule that owns it, its tier and the command that fetches it.

## What is supported today

- **Any UMAT-OTI provider, full-size models.** `resasm history` replays
  small-strain C3D8 analyses with prescribed displacements and many
  increments. `resasm request` hands every model outside its bounded
  single-material scope to it and says so in `run_report.txt`. Both
  full-size cantilevers in [examples/cantilevers](examples/cantilevers/README.md)
  (J2 plasticity, 1,536 C3D8, 40 increments; FCC crystal plasticity, 384 C3D8,
  25 increments) run this way.
- **Residual assembly** of C3D8 models from exported stresses, and bounded
  finite-strain neo-Hookean assembly with first-order parameter sensitivities.
- **Direct and black-box residuals** at arbitrary order through `resasm.yml`.

Features outside that scope (finite-strain plasticity, other element types,
several steps or materials, distributed loads) are refused with a named
reason rather than approximated. [STATUS.md](STATUS.md) lists what is
verified and how.

## Where to go next

| You want to | Read |
|---|---|
| An overview of the whole interface | [README.md](README.md) |
| Install the two packages and OTILib | [docs/INSTALL.md](docs/INSTALL.md) |
| Every command, with real output | [docs/CLI_GUIDE.md](docs/CLI_GUIDE.md) |
| The graphical interface | [docs/GUI_GUIDE.md](docs/GUI_GUIDE.md) |
| Runnable examples | [examples/README.md](examples/README.md) |
| The complete user guide and troubleshooting | [docs/USAGE_REPORT.md](docs/USAGE_REPORT.md) |
| A guided first run | [QUICKSTART_USER.md](QUICKSTART_USER.md) |
| What is verified, and how | [STATUS.md](STATUS.md) |
| To review this code | [REVIEW_GUIDE.md](REVIEW_GUIDE.md) |
| To contribute | [CONTRIBUTING.md](CONTRIBUTING.md) |
| Licence tiers of the external sources | [sources/LICENSING.md](sources/LICENSING.md) |

## The companion product

[UMAT-OTI](https://github.com/AMMS-Lab-UTSA/UMAT_source_transformation)
transforms Abaqus UMAT Fortran so that its derivatives are computed with
order-truncated imaginary (OTI) arithmetic, and compiles it into the provider
object this program consumes. The two are separate products, connected by a
versioned contract ([docs/COMPATIBILITY.md](docs/COMPATIBILITY.md)) and
released independently.

## Licence

GPL-3.0-only; see [LICENSE](LICENSE). External Abaqus sources are pinned
submodules that keep their own licences and are never vendored; see
[THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md).
