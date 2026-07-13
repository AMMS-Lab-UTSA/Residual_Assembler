# PROVENANCE — Bilinear-CZM-UMAT (cohesive) — REFERENCE ONLY

- **Source name**: Bilinear-CZM-UMAT
- **URL**: https://github.com/harshaa765/Bilinear-CZM-UMAT
- **Date reviewed**: 2026-07-13
- **License**: **GPL-3.0 (copyleft)**
- **Source category**: copyleft → **reference only**
- **Subroutine type**: UMAT (cohesive / traction-separation material)
- **Physics**: mixed-mode bilinear cohesive zone model (Mode I/II/mixed), damage
- **Element types**: 2D cohesive elements
- **Input files available**: yes (DCB, ENF, MMB benchmarks + single-element patch)
- **Abaqus validation runnable**: yes
- **Framework backend required**: **none yet** — a cohesive/traction-separation
  formulation (separation kinematics + TSL) is not implemented
- **Residual class**: G (unsupported but useful later)
- **Support status**: **unsupported / reference-only** (license-limited)
- **Next missing item**: a cohesive formulation backend; and GPL-3.0 forbids code
  reuse (reference only)

## Isolation note
**Do NOT copy/adapt/link/import.** Useful conceptual reference for a future
cohesive backend and for classic DCB/ENF/MMB verification cases. The
single-element patch `.inp` would be an ideal first traction-separation check
once an independent cohesive backend exists.
