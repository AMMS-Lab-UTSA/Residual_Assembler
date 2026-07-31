# The transformation contract (and the verification case)

The material developer writes a **compact contract** carrying only genuine
decisions. Everything else is inferred or generated. Property *values* are not a
transformation decision — they say *where* to test the transform — so they live
in a separate **verification case**.

> **The rule:** the contract defines *what* to transform; the verification case
> defines *where* to test it.

## What the developer writes

### The transformation contract

```json
{
  "material": "j2_linear_hardening",
  "source": "umat_j2.for",
  "kinematics": "small_strain",
  "dimensions": {
    "stress_components": 6,
    "material_properties": 4,
    "state_variables": 1
  },
  "parameters": {
    "E": 1,
    "nu": 2,
    "SIGY0": 3,
    "H": 4
  },
  "derivatives": {
    "stress": true,
    "statev": "auto"
  }
}
```

| field | meaning | notes |
|---|---|---|
| `material` | package identifier / label | optional; defaults to the material directory name |
| `source` | the UMAT file | `"umat.for"` or `{ "main": "umat.for", "additional_files": [...] }` |
| `kinematics` | `small_strain` or `finite_strain` | defaults to `small_strain` |
| `dimensions` | interface sizes | `stress_components`, `material_properties`, `state_variables` (Fortran names `ntens`/`nprops`/`nstatev` also accepted) |
| `parameters` | name → PROPS index | the listed order fixes each parameter's OTI direction |
| `derivatives` | requested outputs | `stress: true`; `statev`: `"auto"` \| `true` \| `false` |

**`nparam` is never written.** It is inferred as `len(parameters)`. A contract
that also declares `nparam` (anywhere) is rejected, because a manual value can
contradict the parameter list. `material_properties` (`nprops`) is the size of
the whole PROPS array; `nparam` is how many of those you differentiate — the two
are usually equal but need not be:

```json
"dimensions": { "material_properties": 8 },
"parameters": { "E": 1, "SIGY0": 3, "H": 4 }
```
→ `nprops = 8`, `nparam = 3` (the UMAT uses eight properties; sensitivities are
requested for three).

### The verification case (separate file)

Property values and the loading history go beside the contract, in
`verification.json` (or `<name>.verify.json`):

```json
{
  "props": [210000.0, 0.30, 250.0, 2000.0],
  "loading_path": "uniaxial_tension"
}
```

- `props` — the complete baseline PROPS vector the transform is validated at
  (length must equal `material_properties`; all entries finite).
- `loading_path` — optional; recorded in the validation provenance. When omitted
  the built-in verification path is used (elastic → plastic → unload → reload for
  path-dependent models; a single strain state for elastic ones).

Property values are resolved by precedence:
`--props` → `--verify <file>` → sibling `verification.json` → props embedded in an
older contract (back-compat) → an explicit error. There is deliberately **no**
`[1.0]*nprops` placeholder — validating at meaningless parameters certifies the
transform at an operating point the model never sees (for crystal plasticity it
diverges to NaN).

### A single-file elastic material

```json
{
  "material": "isotropic_elastic",
  "source": "umat_elastic.for",
  "kinematics": "small_strain",
  "dimensions": { "stress_components": 6, "material_properties": 2, "state_variables": 0 },
  "parameters": { "E": 1, "nu": 2 },
  "derivatives": { "stress": true, "statev": false }
}
```

### Unusual UMATs — the optional `options` block

When auto-detection needs a hand, an advanced block overrides it (all optional):

```json
{
  "material": "custom",
  "source": "umat_finite.for",
  "dimensions": { "stress_components": 6, "material_properties": 3, "state_variables": 0 },
  "parameters": { "C11": 1, "C12": 2, "C44": 3 },
  "options": {
    "kinematics": "finite_strain",
    "entry_point": "UMAT",
    "stress_update_line": 41,
    "ddsdde_block": "44-46"
  }
}
```

## What the program infers or generates

The developer never writes any of this — Program 1 produces it:

- the UMAT entry point;
- `nparam`, and each parameter's OTI direction and OTI order;
- `DSIGMA_DP` shape `[ntens, nparam]`, and `DSTATEV_DP` shape `[nstatev, nparam]`
  when `statev` resolves to needed;
- whether `DSTATEV_DP` is required (`statev: "auto"` ⇔ `state_variables > 0`);
- output object / contract filenames;
- the regular and OTI build/compile commands;
- platform, ABI, and binary metadata (OS, arch, compiler + version, format,
  build id, hashes);
- validation tolerances and the finite-difference ladder;
- the **completed interface contract** and the diagnostic report.

## The completed contract (generated, not written)

The build emits `umat_<name>_oti.json` — the detailed
`resasm_umat_oti_contract_v1` the collaborator's tool consumes. It records the
inferred/generated interface, for example:

```json
{
  "material": "j2_linear_hardening",
  "dimensions": { "ntens": 6, "nprops": 4, "nstatev": 1, "nparam": 4 },
  "derived": {
    "nparam": 4,
    "DSIGMA_DP_shape": [6, 4],
    "DSTATEV_DP_shape": [1, 4],
    "oti_directions": { "E": 1, "nu": 2, "SIGY0": 3, "H": 4 }
  },
  "binary": { "os": "...", "arch": "...", "compiler": "...", ... },
  "validation": { ... converged FD evidence, provenance ... }
}
```

## Accepted input shapes (backward compatibility)

All three normalize to the same internal representation, so any of them builds:

| shape | how it's recognised | example material |
|---|---|---|
| **new compact** | `parameters` is a dict (documented above) | `m1_elastic`, `m3_j2` |
| **compact v2** | `parameters` is a list of `{name, props_index}` | most `sweep_*` |
| **verbose** | `interface` + `derivative_requests` | `m2_elastic3d` |

`validate-all` discovers and validates all three in one pass. Discovery is by
*structure*, not a fixed filename or schema tag, so a new-format contract that
omits `schema` is still found; completed contracts (`umat_*`) and verification
cases (a `props` file with no `parameters`) are excluded.
