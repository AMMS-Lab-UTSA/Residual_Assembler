# Interface — the hand-off between the two programs

The framework is deliberately split into two independently-shippable programs so
that a material developer can share differentiable material objects **without
exposing the source**:

```
 ┌─────────────────────────┐        object + contract         ┌──────────────────────────┐
 │  umat-oti  (Program 1)   │  ───────────────────────────►    │ residual-assembler (Prog 2)│
 │  material developer side │   umat_<name>_oti.obj  +  .json   │   collaborator side        │
 └─────────────────────────┘                                   └──────────────────────────┘
```

Only two kinds of artifact cross the boundary, and **no OTI / dual object ever
does** — the compiled object returns *real* values plus, in separate arrays, the
first-order derivative coefficients for each seeded parameter.

**The object is platform-specific; the JSON contract is not.** A compiled object
is COFF (links into a `.dll`) on Windows and ELF (links into a `.so`) on Linux —
one binary cannot serve both. Ship the object built for the collaborator's
platform (or ship both) alongside the single shared JSON. See
[docs/platform_packaging.md](docs/platform_packaging.md).

## What Program 1 emits

1. **`umat_<name>_oti.obj`** — the compiled OTI-enabled material, **per
   platform** (COFF on Windows → `.dll`, ELF on Linux → `.so`). Beyond the
   standard UMAT outputs (`STRESS`, `STATEV`, `DDSDDE`) it returns:
   - `DSIGMA_DP` = ∂STRESS/∂p  (ntens × nparam)
   - `DSTATEV_DP` = ∂STATEV/∂p  (nstatev × nparam), for path-dependent models
   - and exposes `UMAT_OTI_MARCH` for a fast state-alive whole-path replay.

2. **`umat_<name>_oti.json`** — the *completed interface contract*
   (`resasm_umat_oti_contract_v1`), **platform-independent**: `dimensions`
   (ntens/nprops/nstatev/nparam), `parameters` (name → PROPS index → OTI
   direction), `symbols`, `history` (path-dependent? `dstatev_dp` returned?), the
   `object` file + `sha256`, the shared `contract_version`, and a **`binary`
   metadata block** (OS, arch, compiler + version, binary format, ABI version,
   build id, source/transform hashes) that lets Program 2 reject a foreign-platform
   object *before* linking.

Note: the transient `.dll`/`.so` files Program 1 and Program 2 build under
`build/` to call a material through `ctypes` are **temporary validation
libraries**, not the deliverable — they never cross the boundary.

## The shared C-ABI

The boundary is the versioned material C-ABI `resasm_mat_abi_v1` (header:
`resasm_mat_abi_v1.h`, shipped with Program 2). It is identical for the reference
provider and a real JHU OTI binary — only the shared-library path differs. A
combined hash of the header + the material-package schema
(`CONTRACT_VERSION.json`) lets either side detect drift immediately.

## What Program 2 does with them

`residual-assembler` reads the contract, links the object, replays the material at
each recorded integration point to obtain `DSIGMA_DP` (and `DSTATEV_DP`), then
assembles `R`, `∂R/∂p`, `K`, solves `K ∂u/∂p = -∂R/∂p`, and evaluates the
requested `d(response)/dp`. The production analysis is replayed, never re-solved.

## Versioning

Both sides pin the same schema tags and ABI hash. On Program 2 they live in
`residual_core/interface/versions.py`; the contract hash is in
`residual_core/replay/contract/CONTRACT_VERSION.json`. A mismatched tag or a
stale binary is rejected up front with an actionable message.
