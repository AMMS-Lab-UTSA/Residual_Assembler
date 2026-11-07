# Licence-unknown sources: removed from the current tree

`sources/license-unknown/` held two third-party crystal-plasticity codes for
which no redistribution grant exists. They were removed on 2026-11-06. **They had
already been distributed: they remain in this repository's git history.**
Removing them from the current tree does not withdraw those copies; rewriting the
history is a separate decision.

## TarletonGroup/CrystalPlasticity (git submodule)

- Upstream: https://github.com/TarletonGroup/CrystalPlasticity
- Pinned commit: `2f0909472f4cdf1b0b71da6ba97db900dacc6f05`
- Licence: no LICENSE file, so all rights reserved.
- What was removed: the gitlink and its `.gitmodules` entry (`update = none`).
  The submodule was never fetched by setup, so its files were not in this
  repository's tree; `scripts/verify_source_submodules.py` no longer expects
  the pin.

## Huang / Kysar single-crystal UMAT (vendored files)

- Upstream: the Kysar Research Laboratory page
  http://www.columbia.edu/~jk2079/Kysar_Research_Laboratory/Single_Crystal_UMAT.html
  (since removed), retrieved through the Internet Archive snapshot
  `20210219050056` of http://www.columbia.edu/~jk2079/fem/.
- Licence: none stated; academic courtesy-ware with no redistribution grant.
- What was removed, from `sources/license-unknown/Huang_Kysar_Single_Crystal_UMAT/`:

| File | sha256 (prefix) | What it was |
|---|---|---|
| `umatcrystal_mod.f` | `e0e5384d9580` | the single-crystal plasticity UMAT source |
| `umatcryspl_mod.inp` | `a9dd516734e5` | the sample Abaqus input deck |
| `umat_documentation.pdf` | `7d34bb979eb8` | Huang (1991), Mech Report 178 |
| `Addendum_to_umat.pdf` | `2b2e34826530` | Kysar (1997), addendum to Mech Report 178 |
| `PROVENANCE.md` | `7d119f07181f` | this project's provenance note for the folder |

Citation: Huang, Y. (1991), "A user-material subroutine incorporating single
crystal plasticity in the ABAQUS finite element program", Mech Report 178,
Harvard University; Kysar, J. W. (1997), addendum to Mech Report 178.

## What changed with them

- The `.inp` manifest (`residual_core/example_manifest.json` and `.md`) drops
  the 27 decks of these two sources (105 -> 78 files); every other row and the
  classification invariants are unchanged.
- No test, script or library code read these files. The tier stays defined in
  `.gitmodules`, `sources/LICENSING.md` and the setup and CI checks, so anything
  placed there in future is still reference-only and never fetched by setup.
