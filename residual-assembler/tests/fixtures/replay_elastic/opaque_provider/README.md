# FIXTURE ONLY — not part of the collaborator application

This `.f90` is a **stand-in for JHU's shipped OTI binary**, kept here solely so
CI can compile an *opaque* `.so` to exercise the Program-2 consumer through the
real C ABI. It is:

- **excluded** from the installed/distributed collaborator package (see the
  packaging excludes),
- **never imported** by `residual_core/` (the application),
- authored/owned by Program 1 (`UMAT_source_transformation/oti_provider/`).

The consumer tests treat the built `.so` as an opaque binary: they load it only
by path via the public ABI, exactly as they would JHU's real binary.
