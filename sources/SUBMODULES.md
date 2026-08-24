# External source submodules — policy and setup

External Abaqus UMAT/UEL sources are **pinned git submodules**, never vendored
copies. That keeps this framework's licensing clean (it is not GPL) and keeps
provenance exact: every external file is identified by an upstream repository
plus a commit SHA recorded in our tree.

## The three tiers

| Tier | Directory | Licenses | `update` | Fetched by setup? |
|---|---|---|---|---|
| Permissive | `sources/permissive/` | MIT, BSD-3 | `checkout` | Yes |
| Copyleft | `sources/copyleft/` | GPL / AGPL | `none` | **No** |
| License-unknown | `sources/license-unknown/` | none granted | `none` | **No** |

`update = none` in [`.gitmodules`](../.gitmodules) means `git submodule update`
skips the entry **even with `--init`**. Git prints `Skipping submodule ...` and
the directory stays empty. Restricted source therefore cannot arrive as a side
effect of ordinary setup; fetching one has to be a deliberate act by whoever
does it, and they own the license consequences.

The restricted entries are still *mapped* so that `git submodule status` works
and the provenance (URL + pinned SHA) is recorded. Mapping is not vendoring.

## Pinned commits

| Path | Upstream | Pinned commit | License |
|---|---|---|---|
| `permissive/ngrilli_Oxford_Crystal_Plasticity` | ngrilli/Oxford_Crystal_Plasticity | `85102ed35dc4592edc5fd4aaec542fa63e7162fd` | MIT |
| `permissive/bibekanandadatta_Abaqus-UEL-Elasticity` | bibekananda-datta/Abaqus-UEL-Elasticity | `9187e54a4069ed0e6196a29e97e92185d14724a5` | BSD-3 |
| `permissive/bibekanandadatta_Abaqus-UEL-Hyperelasticity` | bibekananda-datta/Abaqus-UEL-Hyperelasticity | `ba5bf018b71e0717d3a3136a351637a27e6ff0ed` | BSD-3 |
| `permissive/jgomezc1_ABAQUS-US` | jgomezc1/ABAQUS-US | `54181407aa7aa23055e33d354d0b2a3abc266365` | MIT |
| `copyleft/ICAMS_Crystal_Plasticity_UMAT` | ICAMS/crystal_plasticity_UMAT | `469464b646a6f250d996828559342b81a1dfbbd7` | AGPL-3.0 |
| `license-unknown/TarletonGroup_CrystalPlasticity` | TarletonGroup/CrystalPlasticity | `2f0909472f4cdf1b0b71da6ba97db900dacc6f05` | none |

Each pinned SHA was confirmed to exist in the named upstream repository. The
upstream `LICENSE` file travels with the submodule content, so attribution and
license text are retained wherever the source is used.

## Setup

```bash
./scripts/init_permissive_sources.sh                  # test dependency only
./scripts/init_permissive_sources.sh --all-permissive # every permissive entry
./scripts/init_permissive_sources.sh --check          # report state, fetch nothing
```

The script refuses any path outside `sources/permissive/` and verifies that each
submodule landed on its pinned commit.

## What the offline test suite needs

Exactly one external source:
`sources/permissive/ngrilli_Oxford_Crystal_Plasticity` — its
`ExampleInputFiles/HCPnoTwin/Compression111.inp` is the real C3D8
crystal-plasticity mesh used by the assembler, recipe and neutral-IO tests.

Without it, `pytest -q` reports **60 passed, 17 skipped**; the eight affected
tests skip with a message naming the file, the owning submodule, its tier and
the command that fetches it. With it, the suite reports **68 passed, 9 skipped**
(the 9 remaining skips are OTILib-dependent and unrelated).

No test depends on a copyleft or license-unknown submodule, so the restricted
tiers can stay permanently empty.

## Verification

```bash
python scripts/verify_source_submodules.py      # offline structural check
./scripts/verify_clean_clone.sh                 # full fresh-clone procedure
```

`verify_source_submodules.py` needs no network: it checks that every gitlink has
a mapping and vice versa, that pinned SHAs are unchanged, that restricted tiers
are `update = none` and unpopulated, and that any populated submodule sits on its
pinned commit with its license file intact.

`verify_clean_clone.sh` clones into a temporary directory, proves a bare
`git submodule update --init` skips the restricted tiers, bootstraps the
permissive dependency, runs the suite, and checks the parent worktree is clean.

## History

Before this was fixed, six gitlinks existed with **no committed `.gitmodules`**.
`git submodule status` failed outright, a fresh clone could not obtain any
external source, and the tests that read one died with `FileNotFoundError` and
assertion failures unrelated to their actual subject. The pinned SHAs were never
changed by the repair — only the missing mapping, the tier policy, the bootstrap
path and the diagnostics were added.
