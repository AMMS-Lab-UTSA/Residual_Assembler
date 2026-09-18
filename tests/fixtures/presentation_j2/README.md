# Genuine Abaqus J2 Export

`fields.json` is copied unchanged from the strict exporter applied to the one
licensed Abaqus 2021.HF5 job `imqrp_j2`, 2026-09-18. It is not synthesized.
The developer job used the original m3_j2 UMAT and the example
`examples/presentation_request/Analysis.inp`, with ifort 2021.10.0.
It contains virgin frame 0 and four converged increments, eight nodes, one
C3D8 element and eight integration points. No ODB or compiled binary is stored here.

Provenance retained outside the repository:
`../imq_abaqus/recovery_presentation/imqrp_reference/`.

SHA-256:

- ODB: `54745b97b515de159c3251e8a3ba0ed89ea1c65452e388ef34c8c81794e3b736`
- INP: `1cde25bd934403a4c3b341d7deca8d5a47897ae700f4329ab448d72d995b74d3`
- Original UMAT: `9b779f0c6cadf9c493068d84cd18bbf91f91ae14cb08aa0e7c8ac4ca119a9b40`

The offline test freshly builds the genuine provider, replays these exported
fields, and compares displacement sensitivities against the independent
uniaxial J2 formula `U1=300/E+(300-SIGY0)/H`. This offline test verifies
conversion/replay, not a new licensed ODB extraction.