# Verification Mode 2 — UMAT replay vs Abaqus

**Goal.** Reproduce the Abaqus UMAT result **outside Abaqus**: drive the Grilli
crystal-plasticity UMAT at a chosen integration point through the *same*
deformation-gradient history Abaqus saw, **increment by increment**, and check

1. replayed `STRESS` matches the ODB integration-point stress `S`, and
2. replayed `STATEV` matches the ODB `SDV`,

then (optionally) assemble a residual from the *replayed* stress field and check
it against Abaqus `RF` (closing the loop with Mode 1, but now with a stress field
we regenerated ourselves).

> **Plasticity is history-dependent — replay the whole path.** The UMAT state
> (`STATEV` + the `/UMPS/` common block) evolves increment to increment. You must
> march the entire loading history; you may **never** jump to the final step and
> "solve" for it. The driver marches the full history in one process so state
> persists (see `residual_core/umat_adapter_fortran/README.md`).

> Abaqus is **not** installed here, and the **real** UMAT object cannot be built
> with gfortran (two ifort-only constructs in the untouchable sources — see the
> adapter README "Build status"). So today this test runs end-to-end only with
> the **mock** material (plumbing/known-answer). The crystal-plasticity numbers
> require building `umat_driver` with **Intel ifort + Abaqus (MKL)**.

---

## 1. Build the driver

```bash
cd residual_core/umat_adapter_fortran
./build.sh                     # or build.bat (Windows cmd)
```

* `build/umat_driver_mock` — always builds (gfortran). Isotropic-elastic mock.
* `build/umat_driver` — the **real** crystal-plasticity driver; only appears when
  `umat.for` is compiled with ifort+Abaqus (or a patched tree) and linked with
  `umat_driver.o` (+ `aba_stubs/aba_stubs.f`, and MKL or `lapack_stub.f`).

Sanity-check the plumbing (known answer: identity deformation → zero stress; a
history-dependent state witness `STATEV(35)` accumulates):

```bash
python umat_replay.py --dry-run --mock
```

## 2. Run Abaqus and extract fields — with DENSE output

Exact replay needs `DFGRD1` at **every** increment. Set field output to **every
increment** and dump **all** SDVs, so each ODB frame equals one solver increment:

```
*Output, field, frequency=1
*Node Output
U, RF
*Element Output, directions=YES
S, SDV
```

```bash
abaqus job=Compression111 user=umat.for double=both interactive
abaqus python /path/to/residual_core/extract_abaqus_fields.py -- \
    --odb Compression111.odb --frames all --out fields.json
```

`extract_abaqus_fields.py` records the exported SDV numbers in `fields.json`
under `sdv_labels`, so the replay lines them up with the right `STATEV` indices.

If you can only afford sparse output (e.g. `frequency=200`), the replay still
runs but frame-to-frame marching becomes a **sub-stepping approximation** of the
true per-increment path; expect drift in the plastic regime and treat the
comparison as qualitative.

## 3. Replay one integration point and compare

```bash
cd residual_core/umat_adapter_fortran
python umat_replay.py \
    --fields  /path/to/fields.json \
    --element 63 --ip 1 \
    --props   "0.,0.89931,-0.4373,0.,0.26422,0.54336,-0.79684,0.34846,0.71661,0.6042,1." \
    --nstatv  125
# add --mock to test the wiring before the real umat_driver exists
```

The orchestrator:
* reads `nodes`/`elements`/`U`/`S`/`SDV` from `fields.json`,
* rebuilds `DFGRD1` at (element, IP) from nodal `U` and the C3D8 shape-function
  gradients (CONTRACT §2/§4),
* marches from the undeformed reference (`F=I`, `t=0`) frame-to-frame with
  `KSTEP=1, KINC=1..` (the first increment triggers the UMAT init block),
* prints per-increment `STRESS` vs `S` and `STATEV` vs `SDV` (abs/rel error).

`--props` must match the material constants of that element's grain
(`PROPS(1)=0` HCP, `PROPS(2:10)` rotation matrix, `PROPS(11)` grain index).

## Pass criteria (real `umat_driver`)

At every replayed increment for which the ODB has a frame:

| Quantity | Check | Tolerance |
|---|---|---|
| `STRESS` vs `S` (6 comps) | `max_i |σ_rep,i − S_i|` | `≤ 1e-6 · max_i|S_i| + 1e-3` (elastic); `≤ 1e-3` relative through plasticity with `frequency=1` |
| `STATEV(48:53)` (Cauchy stress slot) | equals `STRESS` | exact to write precision |
| `STATEV` vs `SDV` (`sdv_labels`) | per-component rel error | `≤ 1e-3` with per-increment history |
| `STATEV(35)` cumulative slip | monotone, matches ODB `SDV35` if exported | `≤ 1e-3` rel |

Larger errors that **grow with accumulated plastic slip** almost always mean the
history was not fully marched (sparse frames, or a skipped `KINC=1` init) — fix
the sampling, do not loosen the tolerance.

## 4. (Optional) residual from replayed stress vs `RF`

Collect replayed `STRESS` at **all** IPs of **all** elements at a converged
frame into `sigma_all[eid] -> (8,6)`, then reuse the Mode-1 assembler:

```python
from c3d8_residual import assemble_global_internal_force
F_int, id2idx = assemble_global_internal_force(
    node_ids, coords, connectivity, U_frame, sigma_all_replayed, mode='finite')
```

and apply the **same** pass criteria as
`tests/abaqus_stress_driven_residual/README.md` (free DOFs ≈ 0; prescribed DOFs =
`±RF`). Passing here means the externally-regenerated stress field is
equilibrium-consistent with Abaqus — the UMAT replay and the residual assembler
agree with Abaqus simultaneously.

---

### Current status in this environment
* `umat_driver_mock` + `umat_replay.py --dry-run/--fields --mock`: **working**
  end-to-end (verified).
* Real `umat_driver`: **blocked on toolchain** (needs ifort+Abaqus/MKL). All
  scripts are ready; only the compiled crystal-plasticity object is missing.
```
