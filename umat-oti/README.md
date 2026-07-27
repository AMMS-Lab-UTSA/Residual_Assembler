# umat-oti — UMAT → OTI transformer (Program 1)

Take a normal Abaqus **UMAT** and automatically produce an **OTI-enabled** material
object that returns, on top of its usual outputs, the exact derivatives of the
stress and state variables with respect to chosen material parameters — computed
by hypercomplex automatic differentiation (OTI / HYPAD), with **no** hand-coded
tangent and **no** finite differences.

This is the **material-developer side** of a two-program framework. It produces a
compiled object plus a JSON *interface contract* that the second program,
[`residual-assembler`](https://…/residual-assembler), consumes to compute offline
material-parameter sensitivities — while the original source stays private.

```
   real UMAT  +  transform contract  ──►  umat-oti  ──►  OTI object (.obj)  +  completed contract (.json)
     (STRESS, STATEV, DDSDDE)                            (… + DSIGMA_DP, DSTATEV_DP)
```

## Install

```bash
pip install -e .            # requires numpy; a Fortran compiler (gfortran) for the build step
```

The transformer emits OTI-promoted Fortran and links it with the bundled
OTILib / MultiZ module generator in [`UMATs/OTI/`](UMATs/OTI); `gfortran` must be
on `PATH` (or pointed at by `FC`).

**Runs natively on Windows and Linux from the same repo and CLI.** The build
compiles two throwaway libraries and loads them through `ctypes`; that load is
made to work identically on both platforms by
[`umat_oti.runtime.libload`](src/umat_oti/runtime/libload.py) — the library gets
the platform's suffix (`.dll` / `.so`), and on Windows the active gfortran
runtime directory is discovered dynamically and registered with
`os.add_dll_directory` (Python 3.8+ ignores `PATH` for *dependent* DLLs, so a
gfortran on `PATH` is necessary but not sufficient). No compiler path is
hard-coded. If a load ever fails you get a diagnostic naming the platform,
compiler, library and missing dependency; `UMAT_OTI_STATIC_FORTRAN=1` links the
runtime in statically as a fallback. Both platforms are exercised in
[CI](.github/workflows/ci.yml).

## Quickstart — transform a material

Each material lives in `oti_provider/materials/<name>/` as a `umat.for` plus a
`contract.json`. Build it:

```bash
python oti_provider/umat_transform.py build oti_provider/materials/m1_elastic/contract.json
```

A successful build prints the parity checks and validates the parameter
derivatives against a **converged** central-finite-difference reference of the
original UMAT, with the per-parameter evidence:

```
 validated at PROPS          : [210000, 0.3]
   source of those values    : contract validation.props_values
 stress real parity  vs orig : 2.10e-16
 DDSDDE parity       vs orig : 5.15e-17
 DSIGMA_DP vs converged FD   : 5.45e-11
 per-parameter FD evidence   :
   param    block        h tested   sel. h    FD self  converged   OTI vs FD  status
   E        stress   1e-03..1e-08    3e-04    3.5e-13   3.5e-14R    3.22e-13  PASS
   nu       stress   1e-03..1e-08    1e-06    1.8e-11   2.3e-12R    5.45e-11  PASS
 -> PASS
```

The finite-difference step is chosen **per parameter, from the FD sequence's own
convergence** (Richardson-checked) — never from agreement with OTI. Validation
runs at real material parameters: `--props` overrides, else the contract's
`validation.props_values`, else an explicit error (there is no `[1.0]*nprops`
placeholder). A `NaN`/`Inf` anywhere is an unconditional failure, never reported
as `0.00e+00`. The methodology is one shared module,
[`umat_oti.validation.fd_reference`](src/umat_oti/validation/fd_reference.py),
used identically by this self-check and by the deck's figure/table scripts; see
[docs/validation_methodology.md](docs/validation_methodology.md).

The build writes the two hand-off artifacts next to the material:
`umat_<name>_oti.obj` (the object) and `umat_<name>_oti.json` (the completed
contract, including the full per-parameter validation record).

To build and validate **every** material at once — all 21 contracts, discovered
by schema (both `contract.json` and the legacy `transform_contract_v2.json`
layouts), with a discovered/executed/passed/failed/skipped accounting:

```bash
python oti_provider/umat_transform.py validate-all
```

## The transform contract (compact)

The developer writes only genuine decisions — the UMAT, the interface sizes, the
parameters to differentiate, and which derivatives are wanted. `nparam`, OTI
directions, array shapes, output names, ABI/binary metadata, tolerances and the
completed contract are all inferred or generated.

```json
{
  "material": "isotropic_elastic",
  "source": "umat.for",
  "kinematics": "small_strain",
  "dimensions": {
    "stress_components": 6,
    "material_properties": 2,
    "state_variables": 0
  },
  "parameters": { "E": 1, "nu": 2 },
  "derivatives": { "stress": true, "statev": "auto" }
}
```

| field | meaning |
|---|---|
| `material` | package label (optional; defaults to the directory name) |
| `source` | the UMAT file — `"umat.for"` or `{ "main": ..., "additional_files": [...] }` |
| `kinematics` | `small_strain` or `finite_strain` |
| `dimensions` | interface sizes: `stress_components`, `material_properties`, `state_variables` |
| `parameters` | name → PROPS index; the order fixes the OTI directions; **`nparam` is inferred** |
| `derivatives` | `stress: true`; `statev`: `"auto"` \| `true` \| `false` |

**Property values live in a separate verification case** (`verification.json`) —
they say *where* to test the transform, not *how* to transform it:

```json
{ "props": [210000.0, 0.30], "loading_path": "uniaxial_tension" }
```

Property values are resolved by precedence `--props` → `--verify <file>` →
sibling `verification.json` → props embedded in an older contract → an explicit
error (no placeholder fallback). The **compact-v2** (`parameters` as a list) and
**verbose** (`interface` + `derivative_requests`) layouts are still accepted and
normalize to the same internal representation. See
[docs/contract_formats.md](docs/contract_formats.md).

## What it produces

- **`umat_<name>_oti.obj`** — the OTI-enabled object. It still returns `STRESS`,
  `STATEV` and `DDSDDE`, and additionally `DSIGMA_DP = ∂STRESS/∂p` and (for
  path-dependent models) `DSTATEV_DP = ∂STATEV/∂p`. For fast whole-path replay it
  exposes `UMAT_OTI_MARCH` (state-alive march).
- **`umat_<name>_oti.json`** (`resasm_umat_oti_contract_v1`) — the completed
  interface contract: dimensions, parameters + OTI directions, symbols, the C-ABI
  version and the object hash. This is the file `residual-assembler` reads.

See [INTERFACE.md](INTERFACE.md) for the exact hand-off between the two
programs.

## What's here

```
src/umat_oti/            the transform engine (parse → promote to OTI → emit → verify)
oti_provider/            build_resmat.py, umat_transform.py, and reference materials/
UMATs/OTI/               the OTILib / MultiZ Fortran module generator
examples/  docs/         examples and documentation
```

Reference materials in `oti_provider/materials/` (isotropic/cubic/orthotropic
elasticity, J2 with several hardening laws, Drucker–Prager, Perzyna, Maxwell,
crystal plasticity, …) double as the transform test suite: each is original,
self-contained, and validated.

## License

MIT (matches OTILib / pyoti) — see [LICENSE](LICENSE). The bundled OTILib / MultiZ
module generator under `UMATs/OTI/` is used under the same terms.
