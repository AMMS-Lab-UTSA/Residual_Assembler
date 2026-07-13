# Verification Mode 1 — stress-driven residual vs Abaqus reactions

**Goal.** Take the stress field Abaqus computed (integration-point Cauchy stress
`S`) plus the displaced geometry (`U`), assemble the internal nodal force
**externally** with the finite-strain C3D8 residual assembler, and show it
reproduces Abaqus equilibrium:

* on **free** DOFs, the assembled internal force ≈ 0 (nothing external acts
  there), and
* on **prescribed-displacement** DOFs, the assembled internal force equals the
  Abaqus **reaction force** `RF` up to sign.

This checks the *residual assembler + field exchange*, taking the UMAT stress as
ground truth (the UMAT itself is checked separately in Verification Mode 2,
`tests/umat_replay_vs_abaqus/`).

> Abaqus is **not** installed in this environment. The commands below are
> ready-to-run; the numeric pass/fail must be produced by an actual Abaqus run.

Convention (CONTRACT §0): these are `nlgeom=YES` C3D8 elements, so Abaqus
assembles `f_int,e = ∫_v B_spatial(x)^T · σ_cauchy dv` in the **current**
configuration. The assembler must use `mode='finite'`
(`element_internal_force_finite_strain`), `x = X + U`, and the Cauchy stress from
`S`.

---

## 0. One-time input edit — request `RF`

The shipped `Compression111.inp` writes only `U` at nodes. Add `RF` (needed for
the reaction comparison) to the field output, e.g.:

```
*Output, field, frequency=200
*Node Output
U, RF
*Element Output, directions=YES
S, SDV
```

(Requesting all `SDV` also helps Mode 2. `frequency=200` is fine for Mode 1,
which only needs converged frames.)

---

## 1. Run the Abaqus job

Copy the model-specific includes next to the input and run the UMAT job. From
`sources/permissive/ngrilli_Oxford_Crystal_Plasticity/ExampleInputFiles/HCPnoTwin/`
(it needs `mycommon.f` and `kMaterialParam.f`; the repo-root copies work — set
`nElements >= 125` in `mycommon.f`):

```bash
# Windows / Linux, inside the Abaqus environment:
abaqus job=Compression111 user=umat.for double=both interactive
```

Produces `Compression111.odb`.

## 2. Extract the fields

```bash
abaqus python /path/to/residual_core/extract_abaqus_fields.py -- \
    --odb Compression111.odb \
    --instance PART-1-1 \
    --step Step-1 \
    --frames all \
    --out fields.json
```

`fields.json` follows CONTRACT §5: `nodes` (undeformed), `elements` (C3D8
connectivity), and per frame `U`, `RF`, `S` (8 IPs, Abaqus order), `SDV`.

## 3. Assemble the residual externally and compare

The integration driver `residual_core/stress_driven_residual.py` is authored
separately (by the orchestrator); it consumes the parsed model
(`abaqus_inp_parser.parse_inp('Compression111.inp')`) + `fields.json` and uses
`c3d8_residual.assemble_global_internal_force(..., mode='finite')`. Invoke it on
the last converged frame, e.g.:

```bash
python ../../residual_core/stress_driven_residual.py --inp <path-to>/Compression111.inp --fields fields.json --frame last --report mode1_report.json
```

It should, for the selected converged frame:

1. build the global internal force `F_int` from `S` (Cauchy) and the deformed
   coordinates `x = X + U`, node-major DOF order (CONTRACT §1);
2. classify DOFs using the model's boundary conditions
   (`Set-1 YSYMM`, `Set-2 XSYMM`, `Set-3 ZSYMM`, `Set-4` prescribed `u3`);
3. compare against `RF`.

---

## Pass criteria

Let `f = F_int` (assembled), `R = RF` (from ODB), at a **converged** frame.
Let `Fchar = max_I |R_I|` be a characteristic reaction magnitude.

| DOF class | Check | Tolerance |
|---|---|---|
| **Free** (no BC, no `CLOAD`) | `|f_I| ≈ 0` | `max_free |f_I| ≤ 1e-3 · Fchar` (loosen to `1e-2` for coarse increments) |
| **Prescribed displacement** (`Set-4` dof 3; symmetry-fixed dofs) | `f_I = ± R_I` | `|f_I − s·R_I| ≤ 1e-3 · Fchar`, `s ∈ {+1,−1}` fixed for the whole model |
| **Global balance** | `Σ_I f_I ≈ 0` and `Σ_I R_I ≈ 0` per component | `≤ 1e-3 · Fchar` |

Notes
* Determine the single global sign `s` once (Abaqus `RF` is the external reaction
  balancing the internal force; with the assembler's sign convention the two are
  equal up to that one sign). Report it.
* "Free DOF" = a DOF that is neither in a `*Boundary` set nor loaded. On these,
  equilibrium says internal force must vanish; a non-zero value flags a bug in
  `B_spatial`, the IP/stress pairing, or the deformed-geometry integration.
* Use a **converged** frame only. Mid-increment/cutback frames are not in
  equilibrium.

## Independent (Abaqus-free) sub-checks already covered elsewhere
The assembler's kinematics are unit-tested without Abaqus in
`c3d8_residual.py` (divergence-theorem patch test, uniform-uniaxial sanity,
tangent finite-difference — CONTRACT §4 / `tests/c3d8_tangent_fd_check/`). Mode 1
adds the end-to-end check against real Abaqus stress/reaction data.
```
