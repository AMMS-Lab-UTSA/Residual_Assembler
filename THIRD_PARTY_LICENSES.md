# Third-party licences

No third-party source is vendored into this repository. External Abaqus
UMAT/UEL sources are **pinned git submodules**, organised by licence tier, and
the restricted tiers are never fetched by ordinary setup.

The authoritative tier policy is [`sources/LICENSING.md`](sources/LICENSING.md)
and [`sources/SUBMODULES.md`](sources/SUBMODULES.md); this file summarises it.

## Tiers

| Tier | Directory | Licences | `update` | Fetched by setup? |
|---|---|---|---|---|
| Permissive | `sources/permissive/` | MIT, BSD-3-Clause | `checkout` | yes |
| Copyleft | `sources/copyleft/` | GPL / AGPL | `none` | no |
| Licence-unknown | `sources/license-unknown/` | none granted | `none` | no |

`update = none` means `git submodule update` skips the entry even with
`--init`. Fetching a restricted source is therefore a deliberate act by whoever
does it, and they own the licence consequences.

## Pinned upstreams

| Path | Upstream | Licence |
|---|---|---|
| `permissive/ngrilli_Oxford_Crystal_Plasticity` | ngrilli/Oxford_Crystal_Plasticity | MIT |
| `permissive/bibekanandadatta_Abaqus-UEL-Elasticity` | bibekananda-datta/Abaqus-UEL-Elasticity | BSD-3-Clause |
| `permissive/bibekanandadatta_Abaqus-UEL-Hyperelasticity` | bibekananda-datta/Abaqus-UEL-Hyperelasticity | BSD-3-Clause |
| `permissive/jgomezc1_ABAQUS-US` | jgomezc1/ABAQUS-US | MIT |
| `copyleft/ICAMS_Crystal_Plasticity_UMAT` | ICAMS/crystal_plasticity_UMAT | AGPL-3.0 |
| `license-unknown/TarletonGroup_CrystalPlasticity` | TarletonGroup/CrystalPlasticity | none |

Exact pinned commit SHAs are recorded in `sources/SUBMODULES.md` and in the
gitlinks themselves. Each upstream `LICENSE` file travels with its submodule
content, so attribution and licence text are retained wherever the source is used.

## External runtime dependencies

**OTILib** is GPL-3.0 and is *not* vendored or linked into the redistributable
core. It is used as an independent, separately distributed component.

**Abaqus** is commercial software from Dassault Systèmes. This repository
contains no Abaqus code; it reads and writes Abaqus file formats and invokes an
installation the user already licenses.
