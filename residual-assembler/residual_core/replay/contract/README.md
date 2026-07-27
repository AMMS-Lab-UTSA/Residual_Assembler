# Material C ABI (`resasm_mat_abi.h`, v1)

The stable, compiler-independent boundary between a collaborator's residual
tool and a material developer's (e.g. JHU's) **closed-source** constitutive
binary. The header `resasm_mat_abi.h` is the normative spec; this file explains
it.

## Why a C ABI and not "just call the UMAT"

A compiled `.obj`/`.o` is a *linkable* object, not a safely-callable model:
symbol names, calling conventions, and Fortran runtime dependencies are
compiler-specific. So the material developer wraps the private model behind a
**single, versioned, C-callable symbol** and the collaborator links only that.

## The two rules that make it private *and* correct

1. **No OTI object ever crosses the boundary.** The binary keeps OTI/dual types
   internal. It returns the real response in one array and, *separately*, the
   first-order derivative coefficient for each seeded parameter in another. A
   `double*` never "contains an OTI number."
2. **Nothing is anonymous.** Every dimension is passed explicitly in
   `resasm_mat_desc_t`; the manifest (see `schemas/resasm_material_package_v1`)
   maps parameter name ↔ PROPS slot ↔ OTI direction ↔ state variable. Without
   that map the tool cannot know what it is differentiating.

## Symbols

| symbol | purpose |
|---|---|
| `mat_describe_v1(desc_out, model_id, cap)` | metadata only: dimensions + model id, so the caller can size buffers and check the twin *before* evaluating. |
| `mat_eval_v1(desc, props, seed_indices, nseed, kin, dkin_dseed, state_in, dstate_dseed_in, time_data, stress, dstress_dseed, state_out, dstate_dseed_out, ddsdde, status)` | evaluate one material point for one increment; return real response **and** the derivative coefficients for every seeded direction. |

## Array layout (row-major, 0-based)

| array | shape | meaning |
|---|---|---|
| `props` | `[nprops]` | parameter values |
| `seed_indices` | `[nseed]` | **1-based** PROPS indices being differentiated |
| `kin` | `[ntens]` (small) or `[18]` (finite: `F0(9),F1(9)`) | kinematic input |
| `stress` | `[ntens]` | real stress |
| `dstress_dseed` | `[ntens*nseed]`, `[i*nseed+s]` | `d stress_i / d props[seed_indices[s]]` |
| `ddsdde` | `[ntens*ntens]`, `[i*ntens+j]` | consistent tangent `d stress_i/d dstrain_j` |
| `state_out` / `dstate_dseed_out` | `[nstatev]` / `[nstatev*nseed]` | updated state + its derivatives |

Voigt order `[11,22,33,12,13,23]`, engineering shear.

## Contract

- **Ownership:** every array is caller-allocated; the callee only writes `*_out`.
- **Optional inputs** (`dkin_dseed`, `state_in`, `dstate_dseed_in`, and the state
  outputs) **may be NULL**; the callee must test before dereferencing. Elastic
  small-strain uses none.
- **Failure is explicit:** a non-zero return means the outputs are undefined.
  The callee must not return partially-valid results. Codes: `ERR_ABI_VERSION`,
  `ERR_DIMS`, `ERR_NULL`, `ERR_SEED`, `ERR_KINEMATICS`, `ERR_STATE`,
  `ERR_CONVERGENCE`, `ERR_INTERNAL`.
- **Versioning:** the symbol is `mat_eval_v1`. An incompatible future contract
  is `mat_eval_v2` — never a silent change to v1.

## Path dependence

The caller replays increment by increment, feeding the previous increment's
`state_out`/`dstate_dseed_out` back in as `state_in`/`dstate_dseed_in`, so the
binary propagates `d(state)/d(p)` forward. Small-strain elasticity is a single
evaluation at the converged strain.

## Reference vs real binary

`../reference/elastic_reference.f90` implements this ABI for isotropic
elasticity with **analytic** derivatives — a REFERENCE to exercise the pipeline,
**not** the real OTI toolchain. The Python side (`../abi.py`) is identical for
both; the real JHU OTI binary drops in by pointing the material package's
`binaries.oti.path` at it.
