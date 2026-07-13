# CASE — Bilinear cohesive-zone UMAT (unsupported)

- **Source**: Bilinear-CZM-UMAT — https://github.com/harshaa765/Bilinear-CZM-UMAT
- **License**: **GPL-3.0 (copyleft)** → reference only
- **Physics**: mixed-mode bilinear cohesive zone (traction-separation), damage
- **Element type**: 2D cohesive
- **User subroutine type**: UMAT (cohesive material)
- **Required residual backend**: **none** — no cohesive/traction-separation formulation exists
- **Residual class**: G — unsupported but useful later
- **Required inputs**: separation kinematics (not strain), a TSL, damage state
- **Expected outputs**: cohesive traction from separation; element residual from the interface weak form
- **Offline tests possible**: no (no backend; and code is copyleft)
- **Abaqus tests required**: yes (for any future reference comparison)
- **Current status**: **unsupported / reference-only**
- **Next missing item**: implement an independent cohesive formulation backend
  (separation kinematics + TSL) — then the single-element patch and DCB/ENF/MMB
  benchmarks become good verification cases. GPL forbids reusing this code.
