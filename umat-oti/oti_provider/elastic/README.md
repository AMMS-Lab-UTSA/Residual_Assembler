# Elastic OTI provider (reference)

Program-1 vertical slice: build a distributable elastic OTI material.

```bash
python build_provider.py        # -> dist/{libmat_elastic_oti.so, material_manifest.json, build_report.json}
```

- `elastic_reference.f90` — the provider material implementing the v1 C ABI
  (`../contract/resasm_mat_abi_v1.h`). Isotropic linear elasticity; derivatives
  wrt E, nu computed **analytically**. This is a REFERENCE, not a real OTI
  toolchain build.
- `build_provider.py` — compiles the `.so`, records hashes/build id, verifies
  the OTI real part matches the regular material (parity), emits the manifest
  and a build report, and packages `dist/`.

Only `dist/` is distributed to collaborators. The `.f90` source stays here.

Requires `gfortran`.
