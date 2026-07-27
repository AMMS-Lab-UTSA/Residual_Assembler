# Validation methodology — what each check proves

The framework makes several different correctness claims, and they are **not
interchangeable**. This note names each one, says exactly what it compares, and
points at the code that computes it, so a reader can tell "the transformed UMAT
reproduces the original stress" from "the assembled FE sensitivity matches an
independent full-analysis finite difference" — two very different guarantees.

The checks form a ladder. Each rung assumes the ones below it and adds one new
thing:

| # | Check | Question it answers | Where |
|---|-------|---------------------|-------|
| 1 | **Primal equivalence** | Does the OTI object return the *same stress* as the original UMAT? | `_validate_path` / `_validate`, `stress_parity_max_rel` |
| 2 | **History equivalence** | Does it return the *same STATEV* (internal state) along the path? | `_validate_path`, `statev_parity_max_rel` |
| 3a | **Local `DSIGMA_DP`** | Is `∂STRESS/∂p` correct at the material point? | `fd_reference.verify_derivatives`, block `stress` |
| 3b | **Local `DSTATEV_DP`** | Is `∂STATEV/∂p` correct at the material point? | `fd_reference.verify_derivatives`, block `state` |
| 3c | **Consistent tangent** | Is `DDSDDE = ∂STRESS/∂DSTRAN` correct? | `_validate_path`, `ddsdde_parity_max_rel` |
| 4 | **Residual assembly** | Does `∂R/∂p = ∫ Bᵀ (∂σ/∂p) dV` assemble without adding error? | Program 2, `results/build_error_tables.py` (P2 rows) |
| 5 | **Full FE sensitivity** | Does `d(response)/dp` from the residual solve match a full perturbed re-analysis? | Program 2 examples, `residual-assembler` |

Rungs 1–3 live in **Program 1** (this repo). Rungs 4–5 live in **Program 2**
(`residual-assembler`). Keeping them distinct matters: a transform can reproduce
the stress perfectly (rung 1) and still have a wrong parameter derivative (rung
3), and a correct material derivative (rung 3) can still be assembled wrong into
the residual (rung 4).

## 1–2. Primal and history equivalence

The OTI object still contains the original UMAT. These checks march *both* the
original UMAT and the OTI object over the same strain path
(elastic → plastic → unload → reload for path-dependent models) and compare the
real outputs pointwise:

- **stress parity** — `max‖σ_OTI − σ_orig‖ / ‖σ_orig‖`, must be ≤ 1e-8
  (it is machine-precision in practice: the *same* stress-update code runs).
- **STATEV parity** — the same, for the internal state vector.

These are exact-reproduction checks, not derivative checks. A non-finite value
anywhere is an unconditional failure (see below).

## 3. Local derivative verification (`DSIGMA_DP`, `DSTATEV_DP`)

This is the claim the whole framework rests on, and it has **one canonical
implementation**: [`umat_oti.validation.fd_reference`](../src/umat_oti/validation/fd_reference.py).
The transformer's build-time self-check and the deck's figure/table scripts in
`residual-assembler/results/` both call it, so a number in a slide and a number
in a build log are the same quantity computed the same way. (They used to each
roll their own fixed finite-difference step, which is how the FCC crystal model
came to "pass" in one place and "fail" in another.)

### The finite-difference reference

For parameter `p` at operating point `p₀`, the reference is a **parameter-scaled
centered difference**

```
F(h) = ( f(p₀ + h·|p₀|) − f(p₀ − h·|p₀|) ) / ( 2 h |p₀| )
```

evaluated over a ladder of *relative* steps `h ∈ {1e-3, 3e-4, …, 3e-8, 1e-8}`.
Scaling by `|p₀|` is what lets one ladder serve parameters that differ by orders
of magnitude (`C11 = 168000` beside `gd0 = 0.001`).

### Choosing the step — from the FD sequence alone

Centered differences carry truncation error `O(h²)` and roundoff error
`O(ε/h)`, so accuracy improves as `h` falls, bottoms out, then degrades. The
step is chosen by **two-sided self-consistency of the FD sequence**: a rung is a
plateau only if it agrees with its coarser neighbour *and* the next finer rung
agrees with it (score = `max(δₖ, δₖ₊₁)`, minimised). A one-sided "smallest
successive change" rule can lock onto an isolated cancellation dip — two
fine-step estimates of a near-zero derivative that happen to differ by ~1e-16 —
and certify a step where roundoff, not truncation, dominates. Requiring the
agreement to persist across two refinements rejects that. **Richardson
extrapolation** of the selected pair is then the reported convergence metric.

**The selection never looks at the OTI value.** Choosing the step that best
agrees with the quantity under test would make the check circular; it would also
have hidden the real FCC story, whose FD reference simply has not converged at
`1e-4`. Only after a converged reference is established is OTI compared to it.

### What is reported, per parameter and per response block

- perturbation sizes tested (the ladder span),
- the FD self-consistency at each rung,
- the selected step and how convergence was judged (Richardson gap `R`, or
  successive change `d`),
- the OTI-vs-FD relative error,
- PASS / FAIL with the reason.

A parameter **PASSes** only when its FD reference converged (Richardson gap ≤ a
tenth of the comparison tolerance) *and* OTI agrees with that converged reference
to within tolerance (1e-4 for path-dependent stress/state, 1e-5 for the elastic
single-point stress). A converged reference that OTI *disagrees* with is a
failure of the transform; a reference that never converges is reported as such
and cannot certify anything.

## Non-finite values are always failures

Any `NaN` or `Inf` — in the original response, the OTI response, an FD estimate,
or an OTI derivative — is an **unconditional failure**. It is never divided away,
clamped, or reported as `0.00e+00`. `fd_reference.require_finite` raises
`NonFiniteResult`, and the build reports `NON-FINITE` and fails. (Validating a
crystal-plasticity model at a placeholder `[1.0]*nprops` produced exactly such
NaNs; that is now impossible — see *operating point* below.)

## Operating point — no placeholder fallback

Validation runs at real material parameters, resolved in this order:

1. `--props v1,v2,…` on the command line, else
2. the contract's `validation.props_values` (verbose layout:
   `resasm_provider.props_values`), else
3. an **explicit error**.

There is deliberately no `[1.0]*nprops` fallback: validating at all-ones
certifies the transform at an operating point the model never sees (for FCC
crystal plasticity that means a zero-width hardening range and a unit rate
exponent, which diverges to NaN). The build prints the exact vector used and
where it came from.

## 4–5. Residual assembly and full FE sensitivity (Program 2)

These are **not** in this repo; they are what `residual-assembler` adds on top of
a correct `DSIGMA_DP`:

- **Residual assembly (rung 4)** — with the material derivative held fixed,
  assemble `∂R/∂p = ∫ Bᵀ (∂σ/∂p) dV` and confirm it introduces no error of its
  own beyond the material derivative. Under controlled homogeneous deformation
  `‖∂u/∂p‖ → 0`, so the volume-averaged stress sensitivity equals the
  material-point value and the assembly is checked directly against it.
- **Full FE sensitivity (rung 5)** — solve `K ∂u/∂p = −∂R/∂p` and compare the
  requested `d(response)/dp` against centered finite differences of an
  *independent full re-analysis* at perturbed parameters. This is the end-to-end
  claim: the residual method reproduces what re-running the whole analysis would
  give, at a fraction of the cost.

Both use the *same* canonical FD reference as rung 3 (imported from this repo via
`residual-assembler/results/_p1_root.py`), so "the residual assembly adds no
error" is measured against a converged reference, not a coincidentally-agreeing
fixed step.

## Cross-platform note

Every check above compiles two throwaway libraries with `gfortran` and loads
them through `ctypes`. That load is made to work identically on native Windows
and Linux by [`umat_oti.runtime.libload`](../src/umat_oti/runtime/libload.py):
the library gets the platform's suffix (`.dll` / `.so`), and on Windows the
active gfortran runtime directory is discovered dynamically and registered with
`os.add_dll_directory` (Python 3.8+ ignores `PATH` for dependent DLLs). No
compiler path is hard-coded. A failed load raises a diagnostic naming the
platform, compiler, library and the missing dependency — not a bare
`FileNotFoundError`.
