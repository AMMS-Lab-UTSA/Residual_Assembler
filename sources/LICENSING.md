# Licensing & Isolation Manifest

This workspace mirrors third-party code. The `sources/` tree is organized by
**license tier** so that a future redistributable wrapper can be built against a
clean permissive core without contaminating it with copyleft or unlicensed code.

## Tiers

### `sources/permissive/` — MIT / BSD-3-Clause  ✅ redistributable core
| Repo | License |
|------|---------|
| `bibekanandadatta_Abaqus-UEL-Hyperelasticity` | BSD-3-Clause |
| `bibekanandadatta_Abaqus-UEL-Elasticity` | BSD-3-Clause |
| `jgomezc1_ABAQUS-US` | MIT |
| `ngrilli_Oxford_Crystal_Plasticity` | MIT |

Safe to combine, modify, and redistribute **with attribution** and retention of the
original license texts.

### `sources/copyleft/` — AGPL-3.0  ⛔ isolated
| Repo | License |
|------|---------|
| `ICAMS_Crystal_Plasticity_UMAT` | AGPL-3.0 (docs: CC BY-NC-SA 4.0) |

Do **not** `include`, compile, statically/dynamically link, or otherwise combine with the
permissive core. AGPL copyleft (including the network-use clause) would relicense the
entire combined work. Use only as an independent, separately distributed component or for
study.

### `sources/license-unknown/` — no license  ⛔ reference only
| Repo | Status |
|------|--------|
| `TarletonGroup_CrystalPlasticity` | No LICENSE file → all rights reserved |
| `Huang_Kysar_Single_Crystal_UMAT` | Academic courtesy-ware, no redistribution grant |

No redistribution rights granted. Study/benchmark reference only until explicit
permission is obtained from the authors.

## Rule of thumb for the wrapper
Build the residual assembler against `sources/permissive/` **only**. If any algorithm from
a copyleft or unlicensed source is needed, re-implement it independently from the published
equations/paper rather than copying code.
