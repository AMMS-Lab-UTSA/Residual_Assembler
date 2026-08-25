# Start here

Abaqus owns the global residual internally and never exposes it. You cannot hand
it a residual function. Residual_Assembler builds R from the ingredients instead
— mesh, element kernels, material response, converged solution, loads,
constraints — and differentiates it.

```
R(u, a) = F_internal(u, a, q) - F_external(a, t) + F_constraints(u, t)
```

**Your model never leaves your machine.** Assembly and sensitivity computation
make no network calls.

## Five minutes

```bash
git clone https://github.com/AMMS-Lab-UTSA/Residual_Assembler.git
cd Residual_Assembler
python -m venv .venv && . .venv/bin/activate
pip install -e ".[yaml]" pytest
./scripts/init_permissive_sources.sh
python -m pytest -q -m "not abaqus and not arc and not network"
```

Python 3.9 or newer. No Abaqus and no licence needed for the offline suite.

`init_permissive_sources.sh` fetches only the MIT/BSD-3 tier of external Abaqus
sources. Copyleft and licence-unknown tiers are marked `update = none` and are
never fetched by ordinary setup. Skip that step and eight tests skip, each
naming the file it needs, the submodule that owns it, its tier, and the command
that fetches it.

## Read the status before planning around this

[`STATUS.md`](STATUS.md) says exactly what works. In short: the assembly beats
are verified for C3D8 — you can assemble and verify R from exported Abaqus
ingredients — but the *sensitivity* beat does not yet run through the C3D8
assembly path.

That distinction matters and is not softened anywhere in this repository.

## Where to go next

| You want to | Read |
|---|---|
| The whole interface | [`README.md`](README.md) |
| A guided first run | [`QUICKSTART_USER.md`](QUICKSTART_USER.md) |
| What actually works today | [`STATUS.md`](STATUS.md) |
| To review this code | [`REVIEW_GUIDE.md`](REVIEW_GUIDE.md) |
| Licence tiers of external source | [`sources/LICENSING.md`](sources/LICENSING.md) |
| To contribute | [`CONTRIBUTING.md`](CONTRIBUTING.md) |

## The companion product

[UMAT-OTI](https://github.com/AMMS-Lab-UTSA/UMAT_source_transformation)
transforms Abaqus UMAT Fortran so its derivatives are computed by
order-truncated imaginary arithmetic. It supplies the differentiated material
this framework consumes. The two are separate products connected by a versioned
contract and are released independently.

## Licence

GPL-3.0-only; see [`LICENSE`](LICENSE). External Abaqus sources are pinned
submodules that keep their own licences and are never vendored -- see
[`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md).
