# PROVENANCE — Mazars_UMAT (concrete damage)

- **Source name**: Mazars_UMAT
- **URL**: https://github.com/marioruiarruda/Mazars_UMAT
- **Date reviewed**: 2026-07-13
- **License**: BSD-3-Clause (LICENSE file) — **with caveat** (see flag)
- **Source category**: permissive-with-caveat
- **Subroutine type**: UMAT
- **Physics**: isotropic concrete damage (modified Mazars model with fracture-energy regularization); softening / history-dependent
- **Element types**: continuum solids
- **Input files available**: partial — `2-Inputs_Outputs.txt` documents I/O; no
  ready `.inp` in the listing
- **Abaqus validation runnable**: yes (needs ifort + Abaqus)
- **Framework backend required**: `solid_c3d8_small_strain` (Mode 2, softening)
- **Residual class**: B
- **Support status**: planned (**hold** pending terms review)
- **Next missing item**: review `1-Permissions_&_Instructions.txt` for
  author-imposed conditions **beyond** BSD-3 before any reuse; then build a
  single-element `.inp`

## Flag
The repository ships a BSD-3 `LICENSE` **and** a separate
`1-Permissions_&_Instructions.txt`. The apparent double-terms must be resolved
(which governs?) before this source is used for anything but reference. Code is
**not vendored** here.
