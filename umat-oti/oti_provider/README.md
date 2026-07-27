# `oti_provider` — Program 1: JHU-side OTI material provider

**Audience:** material-model developers who hold the private source.
**Purpose:** turn private material source into a *matched, OTI-enabled, compiled
material artifact* a collaborator can use **without receiving the source**.

This is one half of a two-program system. The other half — the collaborator's
residual/sensitivity tool — lives in the **`Residual_Assembler`** repo and never
sees any source. The two share only a versioned public **contract** (`contract/`).

```
private source ──▶ OTI/provider build ──▶ dist/ { opaque .so, material_manifest.json, build_report.json }
                                              │
                                              ▼  (only dist/ is distributed)
                                       collaborator tool (Program 2)
```

## What this program does / does not

Does: promote source to OTI, generate the ABI wrapper, compile the OTI (and
identify the regular) material, verify the OTI real part matches the regular
material, produce the manifest, and package only distributable artifacts.

Does **not**: FE mesh handling, ODB extraction, replay-record generation,
residual assembly, `K du/dp` solves, or `dq/dp` — those are the collaborator's
(Program 2). The source never leaves this program.

## Contract (`contract/`)

`resasm_mat_abi_v1.h`, `resasm_material_package_v1.schema.json`, and
`CONTRACT_VERSION.json` (a combined hash). These files are shared byte-for-byte
with Program 2; the manifest carries `contract_version` so a drift is caught.
**Program 1 owns the provider implementation; Program 2 owns the loader.**

## Elastic vertical slice

See `elastic/`. Run:

```bash
python oti_provider/elastic/build_provider.py
```

It compiles `elastic/elastic_reference.f90` (the ABI-conformant provider),
checks real-part parity against the regular material, and writes `elastic/dist/`.
The reference computes derivatives analytically and is **not** a real OTI
toolchain build — it exists to exercise the contract end-to-end.
