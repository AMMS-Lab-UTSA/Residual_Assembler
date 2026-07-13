# Source Inventory — Crystal Plasticity + Brick UEL Codes

Snapshot date: 2026-07-09. All repos cloned `--depth 1` from their default branches
unless noted. Physical layout uses **license tiers** so a permissive redistributable
core stays isolated from copyleft and license-unknown code (see
[LICENSING.md](LICENSING.md)).

```
sources/
  permissive/         MIT / BSD-3  → safe to build a redistributable wrapper against
  copyleft/           AGPL-3.0     → ISOLATED, do not statically/dynamically fuse into a permissive core
  license-unknown/    no license   → reference only, all-rights-reserved by default
```

## Inventory table

| # | Repo (folder) | Role | License | Element types | Integration scheme | Solution variable | Exposes element residual `r`? | Notes |
|---|---------------|------|---------|---------------|--------------------|-------------------|-------------------------------|-------|
| 1 | `permissive/bibekanandadatta_Abaqus-UEL-Hyperelasticity` | Brick UEL (finite strain) | **BSD-3-Clause** | 2D: Tri3, Tri6, Quad4, Quad8 · 3D: Tet4, Tet10, **Hex8**, Hex20 | Implicit, total Lagrangian, Newton | Displacement `u` (nodal) | **Yes** — `RHS = -r`, `AMATRX = dr/du` in `elem_nlmech` | PK-II and PK-I variants. Clean modular Fortran. **Candidate skeleton.** |
| 2 | `permissive/bibekanandadatta_Abaqus-UEL-Elasticity` | Brick UEL (small strain) | **BSD-3-Clause** | Same Lagrange library (2D/3D, incl. Hex8/Hex20) | Implicit, small strain (linear, Newton-structured) | Displacement `u` | **Yes** — `RHS`/`AMATRX` | Simpler sibling of #1; good for first residual plumbing sanity checks. |
| 3 | `permissive/jgomezc1_ABAQUS-US` | Brick/plane UELs + UMATs | **MIT** | UEL8 / UEL9 = 2D 8- & 9-node quads (classical + Cosserat); assorted UMATs | Implicit | Displacement `u` (+ micro-rotation for Cosserat) | **Yes** (UELs build `RHS`/`AMATRX`) | 2D only; useful reference for UEL/UMAT coupling patterns, not a 3D brick source. |
| 4 | `permissive/ngrilli_Oxford_Crystal_Plasticity` | CP UMAT + cohesive UEL | **MIT** | Continuum UMAT (C3D8 etc.) + cohesive/twin UEL | Implicit | Displacement `u`; UMAT returns stress + DDSDDE | UMAT: **No** (stress+tangent only). Cohesive UEL: yes | Grilli/Tarleton, based on Dunne's UEL. HCP/BCC/carbide/olivine. Permissively licensed CP stress update. |
| 5 | `copyleft/ICAMS_Crystal_Plasticity_UMAT` | CP UMAT | **AGPL-3.0** ⚠️ | Continuum UMAT | Implicit | Displacement `u`; stress + DDSDDE | **No** (UMAT only) | Docs are CC BY-NC-SA 4.0. **Keep isolated** — see flag below. |
| 6 | `license-unknown/TarletonGroup_CrystalPlasticity` | CP UMAT (+ CZM UEL) | **None specified** ⚠️ | Continuum UMAT (C3D8/C3D8R) + cohesive-zone UEL | Implicit (Newton; `cpsolver`/`innerloop`) | Displacement `u`; stress + DDSDDE | UMAT: **No**. CZM UEL: yes | OXFORD-UMAT v3.4, paper: doi 10.1016/j.ijsolstr.2024.113110. No LICENSE file → all-rights-reserved; reference only. |
| 7 | `license-unknown/Huang_Kysar_Single_Crystal_UMAT` | Single-crystal UMAT | **None (courtesy-ware)** ⚠️ | Continuum UMAT (works with any solid element) | Implicit (Peirce–Asaro–Needleman rate-tangent) | Displacement `u`; stress + DDSDDE | **No** (UMAT only) | Huang 1991 / Kysar 1997 Harvard report. Downloaded from Wayback (see folder `PROVENANCE.md`). Canonical FCC single-crystal reference. |

### Legend
- **Exposes element residual `r`?** — "Yes" means the code computes an element-level
  internal-force residual vector and its Jacobian (Abaqus `UEL`: `RHS`, `AMATRX`).
  "No" means it is a material-point routine (Abaqus `UMAT`: `STRESS`, `DDSDDE`) and the
  element residual is assembled by Abaqus internally and never exposed.
- All CP material routines here are **UMATs** → they give the *constitutive* response
  (stress + consistent tangent), not the element residual. The brick **UELs** are where
  `r(u,a)` becomes explicit.

## GPL-3.0 / AGPL-3.0 flag (task step 3)

| Repo | Copyleft class | Action |
|------|----------------|--------|
| `copyleft/ICAMS_Crystal_Plasticity_UMAT` | **AGPL-3.0** | **Flagged.** Isolated under `sources/copyleft/`. Must **not** be compiled, `include`d, or linked into any permissive redistributable core. AGPL network-copyleft would virally relicense the combined work. Treat as study/reference or as a fully separate, independently distributed component only. |

- **No GPL-3.0 repos** were found in this set.
- Two repos carry **no license at all** (`TarletonGroup_CrystalPlasticity`,
  `Huang_Kysar_Single_Crystal_UMAT`). Absence of a license means **all rights reserved**
  by default — for redistribution purposes these are *more* restrictive than AGPL. They
  are isolated under `sources/license-unknown/` and must be treated as reference material
  only until explicit permission is obtained.

## Redistribution-safe permissive core

For a wrapper we intend to redistribute, only these are clean to build against:

- `permissive/bibekanandadatta_Abaqus-UEL-Hyperelasticity` (BSD-3) — UEL skeleton
- `permissive/bibekanandadatta_Abaqus-UEL-Elasticity` (BSD-3)
- `permissive/jgomezc1_ABAQUS-US` (MIT)
- `permissive/ngrilli_Oxford_Crystal_Plasticity` (MIT) — CP stress update

BSD-3 and MIT are compatible with each other and permit redistribution with attribution,
so a `Hex8 UEL (BSD-3) + CP UMAT stress update (MIT)` combination is redistributable
provided both license/attribution notices are retained.
