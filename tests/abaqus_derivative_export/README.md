# M3 — real Abaqus ODB → derivative fields → residual sensitivities

Milestone M3 closes the loop from a **real Abaqus analysis** to the field-driven
residual method: an Abaqus ODB is exported to the canonical
`resasm_derivative_fields_v1` JSON that the M2 recipe already consumes, and the
resulting `du/da` is validated against full finite-difference Abaqus reruns.

```
abaqus job=model user=elastic_export_umat.for      # solve C3D8 + verification UMAT
abaqus python export_derivative_fields.py -- \      # ODB -> canonical JSON
    --odb model.odb --step Step-1 --frame -1 \
    --layout derivative_layout.json --output derivative_fields.json
resasm run sensitivity.yaml                          # K, R_,a -> du/da
```

## Files

| file | purpose |
|---|---|
| `elastic_export_umat.for` | verification UMAT (M3.4): isotropic elasticity; writes exact `DDSDDE`→SDV1–36, `dσ/dE`→SDV37–42, `dσ/dν`→SDV43–48 |
| `nonuniform_c3d8.inp` | distorted single-C3D8, base fixed, mixed top loads → **distinct strain at all 8 IPs** (a uniform test cannot check IP order) |
| `derivative_layout.json` | `resasm_sdv_layout_v1` sidecar (M3.2): the SDV→semantics map (single source of truth) |
| `run_chain.sh` | runs the full chain into a paren-free temp dir |
| `fixtures/` | outputs of one real **Abaqus 2021** run (bit-reproducible), so the tests validate offline |

The exporter itself is `residual_core/io/export_derivative_fields.py` (reuses the
per-IP extraction in `residual_core/io/abaqus_odb_export.py`, **not** the flat
`scripts/extract_odb_fields.py`).

## What is proven (`tests/framework/test_abaqus_derivative_export.py`)

- **M3.1/M3.2** exporter + sidecar produce a valid schema and reject malformed input.
- **M3.5** `K` from the ODB-exported `DDSDDE` vs the analytic C3D8 tangent: **2.6e-8**.
- **M3.3** per-IP stress recovered from the export matches Abaqus `S` at the same IP
  label (**7.6e-8**); an intentional IP permutation is detected (rel ~1.1); and the
  exported stress reproduces the applied load (equilibrium), which a permuted order
  breaks. Abaqus IP label *k* ↔ kernel Gauss point *k*.
- **M3.6** `run_field_recipe` reproduces the Abaqus displacement sensitivities.
  `du/dν` vs central-difference Abaqus reruns is **5.5e-6** (< 1e-5). `du/dE` is
  validated two ways because its finite-difference oracle is precision-fragile
  (see below): the **accurate elastic oracle** `du/dE = -u/E` gives **6.3e-8**
  (correct to the ODB float32 floor), and the FD comparison reaches **9.4e-6** at
  the optimal step (h/E≈0.002).

### du/dE finite-difference step-size plateau

For linear elasticity `u ∝ 1/E`, so `du/dE = -u/E` is tiny (~1e-8 here) and its
central finite difference subtracts two nearly-equal **single-precision** ODB
displacement vectors. The error is therefore a numerical U-curve (real Abaqus
E-perturbation reruns; `fixtures/fd_sweep_E.json`):

| h/E | rel(du/dE) | regime |
|---|---|---|
| 0.0005 | 3.9e-5 | cancellation floor (h too small) |
| 0.0010 | 3.5e-5 | |
| **0.0020** | **9.4e-6** | **minimum — meets the 1e-5 target** |
| 0.0040 | 2.1e-5 | |
| 0.0080 | 6.4e-5 | truncation growth (h too large) |
| 0.0160 | 2.6e-4 | |

The minimum (9.4e-6) meets the original 1e-5 target; the growth at large h
confirms this is a numerical (storage-precision) plateau, not a residual
formulation offset. The precision-independent `du/dE = -u/E` oracle (6.3e-8) is
the authoritative elastic check; the Abaqus FD is a storage-limited end-to-end
cross-check. `du/dν` is well-resolved by FD directly because it is O(1e-3).

A gated `test_live_abaqus_regenerates_fixture` re-runs the whole chain and confirms
it reproduces the fixtures (`RESASM_RUN_ABAQUS=1`, Abaqus on PATH). Verified here:
regen vs fixture = **0.0** (bit-identical).

## M4 — OTI-parameter-seeded UMAT (automatic differentiation)

`elastic_oti_umat.for` is the transformed counterpart of the analytic UMAT: it
obtains `dσ/dE` and `dσ/dν` by **OTI automatic differentiation** instead of by
hand. Using the `OTIM4N1` library (`~/MultiZ_f/oti`, type `ONUMM4N1`, m=4
independent first-order directions):

```fortran
      USE OTIM4N1
      E_OTI  = PROPS(1) + E1        ! seed E along direction 1
      NU_OTI = PROPS(2) + E2        ! seed nu along direction 2 (simultaneously)
      ... ordinary constitutive equations in ONUMM4N1 arithmetic ...
      SIG_OTI = MATMUL(D_OTI, EPS)  ! eps = STRAN + DSTRAN (total strain)
      STRESS(I)    = SIG_OTI(I)%R   ! real stress from the same OTI result
      STATEV(36+I) = SIG_OTI(I)%E1  ! d sigma_i / dE
      STATEV(42+I) = SIG_OTI(I)%E2  ! d sigma_i / dnu
      DDSDDE(I,J)  = D_OTI(I,J)%R   ! exact real tangent (unchanged from M3)
```

Link recipe (Abaqus 2021 uses ifort 2021.10.0 — the same compiler that built the
`.mod`, so no version mismatch): a local `abaqus_v6.env` adds `-I` (for the `.mod`)
and `-L -l otim4n1`; `run_chain.sh <workdir> elastic_oti_umat.for` does it all.

Acceptance (tests `test_oti_umat_*`): the OTI export reproduces the analytic M3
export **bit-for-bit** — `dσ/dE`, `dσ/dν`, `DDSDDE`, and nodal `U` all match to
0.0, and the OTI fields drive the residual pipeline to the identical `du/da`. The
OTI automatic derivative *is* the analytic derivative.

### FCC / cubic C11,C12,C44 (`cubic_oti_umat.for`, `fcc_cubic_c3d8.inp`)

The same pattern with **three** parameters: a cubic single crystal ([100]) seeds
C11/C12/C44 along OTI directions E1/E2/E3 in one pass and writes `dσ/dC11`,
`dσ/dC12`, `dσ/dC44` to SDV37–54 (`fcc_layout.json`; constants are Cu:
168400/121400/75400 MPa). Validated (`test_fcc_cubic_oti.py`, fixtures `fcc_*`)
on precision-independent identities that need no separate analytic UMAT:

- **Euler-on-stress** `C11·∂σ/∂C11 + C12·∂σ/∂C12 + C44·∂σ/∂C44 = σ`: 4.5e-8.
- **Analytic** `∂σ/∂Cij = (∂D/∂Cij)·ε`: 4e-8 each.
- **Euler-on-displacement** `ΣCᵢ·du/dCᵢ = −u` (the finite-difference-free
  sensitivity oracle, exact because K is degree-1 homogeneous in Cij): 2.6e-8.
- `du/dCij` vs full Abaqus PROPS-perturbation reruns: 1.3e-5–3.1e-5 (float32 floor).

`run_chain.sh <workdir> cubic_oti_umat.for fcc_cubic_c3d8.inp fcc_layout.json`
drives the live chain. This is the configuration the FCC `simpleTension` case uses.

### Oriented cubic crystal (`cubic_oriented_oti_umat.for`)

The aligned-[100] test above is *blind to a bad rotation*. `cubic_oriented_oti_umat.for`
adds a non-trivial Bunge orientation (phi1,PHI,phi2 = 25,40,15 deg) using a
**single OTI-rotated operator** for STRESS, DDSDDE and every dsigma/dCij (rotate
strain into the crystal frame `eps_c = G eps G^T`, cubic law, rotate stress back
`sig = G^T sig_c G`; the real tangent is the same operator's real part). Validated
(`test_fcc_oriented_oti.py`, fixtures `ori_*`) against an **independent 4th-order
tensor rotation** `Cs_ijkl = G_mi G_nj G_ka G_lb Cc_mnkl` (a different method than
the UMAT):

- rotated **DDSDDE** vs independent 4th-order rotation: 2.9e-8 (with genuine
  normal-shear coupling `D[0,3]=-19513` that is zero when unrotated);
- Euler-on-stress 4.2e-8; rotated analytic `dσ/dCij` 7e-8; Euler-on-displacement
  7.6e-8; `du/dCij` vs Abaqus FD 1.9e-5-3.2e-5 (float32 floor).

This closes the crystal-orientation gap -- the last cheap checkpoint before a
history-dependent crystal-plasticity UMAT.

## Two findings that shaped M3

1. **Abaqus C3D8 uses selective reduced (mean-dilatation B-bar) volumetric
   integration.** A full-integration kernel does *not* reproduce it — stresses
   differ by a per-IP hydrostatic term, and `du/da` was off ~28%. Mean-dilatation
   B-bar recovery matches Abaqus `S` to ~1e-7. The engine now has an additive
   `integration="selective_reduced"` path (`c3d8_kernel.element_tangent_bbar` /
   `element_internal_force_bbar`, selected via the field metadata); the M1 default
   full-integration path is unchanged. This is what makes `du/da` match Abaqus.

2. **The Abaqus ODB stores field output (S, U, SDV) in single precision.** Every
   ODB-exported quantity is good to ~1e-7 relative, not double. So M3.5's `<1e-10`
   target is physically unreachable through an ODB (2.6e-8 is the float32 floor),
   and `du/dE` — which is intrinsically ~1e-8 in magnitude — bottoms out near 2e-5
   in the finite-difference comparison. `du/dν` (larger) resolves to 5.5e-6.

## Environment note

Abaqus 2021 on a newer glibc prints a benign `buffer overflow detected` abort
during job **wrap-up**, *after* the ODB is written (the `.sta` reports
`COMPLETED SUCCESSFULLY`). Scripts check the `.sta`/ODB, not the exit code.
Also: **UMAT compilation fails in any path containing parentheses** (this repo is
`.../Residual_Assembler(2)/...`) — always build in a paren-free temp dir.
