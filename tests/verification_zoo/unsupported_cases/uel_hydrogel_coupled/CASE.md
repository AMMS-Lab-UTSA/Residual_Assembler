# CASE — Coupled chemo-mechanics UEL (hydrogel, unsupported)

- **Source**: Abaqus-UEL-Hydrogel — https://github.com/bibekananda-datta/Abaqus-UEL-Hydrogel
- **License**: **custom / non-SPDX** → reference only until terms confirmed
- **Physics**: coupled chemo-mechanics (species diffusion + finite deformation)
- **Element type**: coupled continuum (displacement `u` + chemical-potential `μ` DOFs)
- **User subroutine type**: UEL (multi-field)
- **Required residual backend**: **coupled-field formulation (Mode 5)** — not implemented
- **Residual class**: E — coupled-field formulation
- **Required inputs**: mesh + heterogeneous DOF layout (u + μ) + callable/exported coupled `RHS`/`AMATRX`; PROPS
- **Expected outputs**: coupled residual over both fields; block-structured tangent
- **Offline tests possible**: no (no coupled-field DOF wiring / weak form yet)
- **Abaqus tests required**: yes (future reference)
- **Current status**: **reference-only** (license custom + backend missing)
- **Next missing item**: confirm the custom license; implement coupled-field DOFs
  (the `DofManager` already supports heterogeneous DOF sets) and the coupled weak
  form. This is the primary motivator for the Mode-5 roadmap item.
