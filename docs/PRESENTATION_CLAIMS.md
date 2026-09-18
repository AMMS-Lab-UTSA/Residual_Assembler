# IMQCAM Annual Meeting presentation: claims, commands, references, measured values

Deck: `IMQCAM_Annual_Meeting_20260812_V4.pptx` (11 August 2026), 42 slides.
Every quantitative claim of the deck is listed below with the command that
reproduces it from the two repositories, the independent reference it is
checked against, the value measured on 2026-09-18 and its status. A slide number
that did not reproduce is reported as measured, with the reason; the slide's own
number is never copied into a result.

Status words: **reproduced** (the measured value supports the slide's statement),
**reproduced, differs** (the capability works and is verified, but a number on
the slide is not what the current code measures), **not reproduced** (with the
reason), **owned by ...** (a claim another lane of the 2026-09-18 programme is
responsible for; see `COORDINATION.md` in the work directory).

## How to run

```sh
PY=/path/to/python3.11
export UMAT_OTI_REPO=/path/to/UMAT_source_transformation
export PYTHONPATH=$PWD:$UMAT_OTI_REPO/src:/path/to/otilib/build_py311
export PRESENTATION_WORK=/scratch/presentation_work   # large intermediates (never in the repo)
$PY presentation/run_all.py            # offline parts of every claim (about 15 min)
$PY presentation/run_all.py --abaqus   # plus the Abaqus jobs of claims 3 and 4 (one job at a time)
$PY presentation/run_all.py --quick    # reduced cases (a few minutes)
$PY -m pytest -q tests/presentation -m "not abaqus and not slow"
```

Each claim also runs on its own (`presentation/claimN_*.py --help`). Results are
written to `presentation/results/` (git-ignored; small JSON/CSV) and summarised in
`presentation/results/presentation_summary.json`; the committed copy of the run
reported here is `presentation/expected/presentation_summary.json`. Requirements:
gfortran (all claims), ifort (claim 2 timing with the original benchmark's
compiler), Abaqus 2021 with ifort (claims 3 and 4 with `--abaqus`), and, for claim
2 only, the flow-rule sources that are not part of either repository (below).

## Summary table

| Slide | Claim | Command | Independent reference | Measured (2026-09-18) | Status |
|---|---|---|---|---|---|
| 13 | 18/20 models, 76 parameter directions, DSIGMA_DP and dσ_vM/dp within 1.6e-7 of FD, both < 1e-5 | `presentation/claim1_sensitivity_sweep.py` | centred FD of the ORIGINAL UMAT re-marched over the whole path (step ladder, plateau) | **20/20 models, 84 directions**, worst 1.56e-7 (m6_fcc DSIGMA_DP; inside its FD plateau uncertainty 3.2e-7); strict max-over-path metric also < 1e-5 for all 20 | reproduced, differs (more models pass than the slide shows) |
| 26–27 | OTI vs hand-coded flow-rule Jacobian: max rel. error ~1e-15; FD error ~1e7× larger; FD ~9× slower | `presentation/claim2_flowrule_jacobian.py` | the hand-coded analytical derivatives; centred/forward FD of the analytical routine | OTI vs analytical **7.1e-15**; FD/OTI error ratio **4.5e5** (centred, best step) … **2.9e6** (centred, plateau step) … **1.7e8** (forward, best step); FD time **6.0×** OTI (ifort -O2, known step), **41.9×** with the step search | first statement reproduced; FD ratios depend on the unrecorded FD protocol: reproduced, differs |
| 28–32 | single C3D8 shear, six σ_vM sensitivities in one OTI run; regimes; OTI vs analytic chain rule NRMSE < 1e-8; same on a larger mesh; plain vs HYPAD vs FD cost | `presentation/claim3_cp_residual_c3d8.py [--abaqus]` | hand-derived chain rule (implicit-function theorem) of the same update; centred FD of the ORIGINAL with the C3D8 re-solved; Abaqus | NRMSE **2.9e-16**; 4×4×4 vs 1 element **1.0e-14**; Abaqus σ_vM path vs this solver **1.1e-15**, FD of Abaqus vs OTI **≤ 6.2e-11**; HYPAD **9.2–10.0×** a plain UMAT pass (compiled kernel; 2.5× in the whole Python analysis) vs **13×** for central FD; regime statements: 3 of 5 hold, "H small" and "τ0, ΔG, q dominate" do not | reproduced for the surrogate model m5_cpflow; the slide's own 12-slip model is not in the repositories |
| 8 | 18/18 benchmark DDSDDE verified in Abaqus, 12 exact, 6 differ (notching/rounding/error in original, spin_elas_def 740 / 2.6e-3) | `presentation/claim4_benchmark_ddsdde.py --abaqus --variants` | the original UMAT's own DDSDDE in a paired Abaqus run (claudeP_ jobs) | committed contracts: **18 of 18 run, 17 pass (12 exact, 5 within tolerance)**; NKH fails because its source reads DTHTA unset; with NKH's labelled variant row (PROPS(1) = 0 and an initial temperature) **18/18 verified, 12 exact, 6 within tolerance**; absolute differences 0.7234 (NKH), 0.2370 (VPDCL), 0.3716 (VPDCO), 0.03125, 0.015625 as on the slide; spin_elas_def 868.3 (slide 740) | reproduced (NKH through a labelled source-defect variant); four transformer/contract gaps fixed on the way |
| 25 | 19 internal Jacobian entries, 14 exact, 5 within 1e-5 | `presentation/claim5_constitutive_jacobians.py` | centred FD of the ORIGINAL at the converged local solve; the hand-coded value | all **21** (UMAT, symbol) pairs present in the sources measured; OTI agrees with FD in every one (≤ 2e-9); hand-coded vs OTI: **5 exact, 10 within 1e-5 (8 of them ≤ 1e-12), 6 differ** — NKH ANP1P/BNP1P 1.9e-4, VPDCO and VPDCL_R FJAC 2.6e-3 and GDIA 6.7e-5, the hand-coded value being the wrong one (FD) | reproduced, differs (hand-coded Jacobians of NKH, VPDCO, VPDCL_R are wrong where the slide says Exact/Pass) |

The other slides are listed in [All slides](#all-slides).

## Claim 1 – slide 13: parameter sensitivities of the model collection

**What produced the slide.** RA `origin/cross-platform-hardening`,
`results/build_sweep_tables.py` and `results/figures/sweep_error_tables.json`:
20 materials were run and the 18 on the slide are those below 1e-5 on both
programs (sweep_drucker_prager Program-1 6.7e-3 and sweep_perzyna_linear
Program-2 1.5e-4 were dropped; 84 − 8 = 76 directions). Program 1 was the OTI
DSIGMA_DP of the compiled provider marched over a 150-increment uniaxial strain
ramp to 3 % (dt = 1/150), scored as RMSE(OTI − FD)/max|FD| per parameter.
Program 2 was dσ_vM/dp **at one material point** (the script's docstring said
"single C3D8", but the code marched one integration point: under prescribed
homogeneous strain du/dp = 0) at the last increment of a 120-increment ramp.

**What is measured now.** For every model of the UMAT repository's own list
(`tools/run_parameter_sensitivity_sweep.py --list`, 20 models):

* Program 1, source transform: the UMAT repository's sweep tool (transform + OTI
  material-point driver vs centred FD of the separately compiled ORIGINAL), all 20
  models: funnel 20 transformed, 20 compiled, 20 primal parity, **19 verified**,
  **83/84 directions verified**, 14 539 of 14 540 comparison rows agree; the one
  unresolved row (m6_fcc) is a noise-floor row the tool withholds, not a
  disagreement.
* Program 1, compiled provider: `umat-oti-provider build` (the module
  `umat_oti.provider`) for every model, then DSIGMA_DP through `UMAT_OTI_EVAL`
  (incoming derivatives carried for path-dependent models; stateless models at
  the total strain from the virgin state) over the slide's 150-increment path.
* Program 2, the Residual Assembler replay: `residual_core.replay.path_material.
  PathMaterial` links the compiled object; dσ_vM/dp = (∂σ_vM/∂σ):dσ/dp at every
  increment.
* Reference: centred FD of the ORIGINAL UMAT (the untransformed source bundled in
  the provider object, symbol `umat`), the **whole path re-marched** for every
  perturbation, steps 1e-3 … 1e-7 relative, the step chosen per parameter by the
  plateau of consecutive steps (never by agreement with OTI). Where a perturbed
  run switches between the elastic and inelastic branch at an increment boundary
  (a kink: the two one-sided derivatives differ), the second-order one-sided
  stencil on the nominal side is used; this happens only for
  sweep_drucker_prager, which first yields exactly at the end of increment 5.
* Cross-checks: the UMAT sweep tool's transform driver vs the provider's
  `UMAT_OTI_MARCH` on the identical path: ≤ 4.2e-16 for all 20 models (separate
  builds and drivers of the same transformer). `UMAT_OTI_MARCH` vs the EVAL replay:
  identical (≤ 4e-15) except sweep_real_ECL_TEMP, whose stress depends on TEMP:
  MARCH (like the Program-1 driver) runs at TEMP = 293.15, the replay shims at
  TEMP = 0 — a convention difference between the two entry points, reported, not a
  derivative error.

Slide metric (Program 1: RMSE over path and components / max|FD|; Program 2:
|OTI − FD| / |FD| at the final increment), and the stricter max-over-path metric:

| Model | Physics | Params | DSIGMA_DP (P1) | dσ_vM/dp (P2) | P1 max over path | P2 max over path | UMAT sweep tool | Slide P1 / P2 |
|---|---|---:|---:|---:|---:|---:|---|---|
| m1_elastic | Isotropic elasticity | 2 | 7.1e-11 | 5.1e-10 | 2.0e-10 | 5.1e-10 | verified | 1.0e-10 / 1.6e-10 |
| m2_cubic | Cubic anisotropic elasticity | 3 | 3.6e-13 | 1.4e-12 | 1.8e-12 | 1.8e-12 | verified | 1.1e-10 / 6.3e-11 |
| m3_j2 | J2 plasticity, linear hardening | 4 | 7.7e-11 | 1.3e-08 | 1.6e-10 | 5.3e-10 | verified | 7.7e-11 / 4.5e-08 |
| m5_cpflow | Viscoplastic flow | 6 | 1.8e-09 | 3.1e-10 | 5.2e-09 | 7.6e-10 | verified | 1.4e-10 / 2.8e-10 |
| m6_fcc | FCC crystal, 12 slip systems | 10 | 1.6e-07 | 3.5e-09 | 3.6e-07 | 5.2e-09 | unresolved (1 row) | 1.6e-07 / 2.2e-08 |
| sweep_aniso_ortho | Orthotropic elasticity | 9 | 3.9e-13 | 2.2e-10 | 2.3e-12 | 2.2e-10 | verified | 1.4e-10 / 1.4e-10 |
| sweep_damage_elastic | Elasticity + scalar damage | 3 | 1.0e-10 | 5.5e-10 | 4.5e-10 | 5.5e-10 | verified | 2.8e-10 / 4.9e-10 |
| sweep_drucker_prager | Drucker-Prager plasticity | 4 | 1.1e-10 | 2.5e-10 | 2.9e-10 | 2.6e-10 | verified | not on slide (old P1 6.7e-3) |
| sweep_eco | Cosserat elasticity (ECO) | 2 | 7.1e-11 | 5.1e-10 | 2.0e-10 | 5.1e-10 | verified | 8.6e-11 / 1.6e-10 |
| sweep_j2_bilinear | J2 plasticity, bilinear | 4 | 1.1e-10 | 3.8e-10 | 2.6e-10 | 5.3e-10 | verified | 1.1e-10 / 8.5e-09 |
| sweep_j2_combined | J2 plasticity, combined | 5 | 1.1e-10 | 1.3e-09 | 2.7e-10 | 5.3e-10 | verified | 1.1e-10 / 1.1e-08 |
| sweep_j2_kinematic | J2 plasticity, kinematic | 4 | 1.1e-10 | 2.7e-09 | 2.8e-10 | 5.3e-10 | verified | 1.1e-10 / 6.9e-08 |
| sweep_lame_elastic | Isotropic elasticity (Lamé) | 2 | 1.3e-12 | 4.3e-09 | 5.8e-12 | 4.5e-09 | verified | 3.2e-10 / 2.8e-10 |
| sweep_maxwell_ve | Maxwell viscoelasticity | 3 | 3.0e-10 | 1.2e-09 | 9.7e-10 | 1.2e-09 | verified | 3.0e-10 / 2.2e-10 |
| sweep_mooney_small | Neo-Hookean (small strain) | 2 | 1.2e-12 | 1.2e-09 | 6.1e-12 | 1.2e-09 | verified | 3.7e-10 / 6.8e-08 |
| sweep_perzyna_linear | Perzyna viscoplasticity | 4 | 1.0e-10 | 8.7e-08 | 2.6e-10 | 5.3e-10 | verified | not on slide (old P2 1.5e-4) |
| sweep_real_ECL_TEMP | Thermo-elasticity (ECL_TEMP) | 5 | 1.7e-09 | 2.0e-08 | 4.2e-09 | 2.0e-08 | verified | 1.1e-10 / 5.5e-08 |
| sweep_real_PCO | Couple-stress plasticity (PCO) | 4 | 8.0e-11 | 1.3e-08 | 1.7e-10 | 5.3e-10 | verified | 8.0e-11 / 4.5e-08 |
| sweep_thermoelastic | Thermo-elasticity | 3 | 7.1e-11 | 5.1e-10 | 2.0e-10 | 5.1e-10 | verified | 8.6e-11 / 1.6e-10 |
| sweep_transiso | Transversely isotropic | 5 | 4.9e-10 | 5.5e-09 | 1.6e-09 | 6.2e-09 | verified | 5.4e-09 / 1.5e-08 |

**Measured vs slide.** 20 models and 84 directions below 1e-5 on both programs
(slide: 18 and 76); worst 1.56e-7 (slide: 1.6e-7, the same m6_fcc entry). The two
models the slide dropped now pass: Drucker-Prager's old 6.7e-3 was the centred FD
straddling the yield kink at increment 5 (the one-sided stencil resolves it to
1.1e-10), and Perzyna's old 1.5e-4 came from the old fixed FD step. The m6_fcc
1.56e-7 is inside the reference's own plateau uncertainty (3.2e-7 for gd0), i.e.
it measures the FD reference, not OTI. Numbers are not identical to the slide's
because the FD reference is different (plateau-selected ladder vs the old fixed
step) and Program 2 is scored along the whole path, not only at a 120-increment
ramp's last increment. The residual method "on a single C3D8" of the old table
was one integration point (du/dp ≡ 0); the genuine C3D8 residual solves are in
claim 3 and in `tests/integration/test_connected_j2.py`.

## Claim 2 – slides 26–27: OTI vs a hand-coded crystal-plasticity Jacobian

**Which routine.** The "crystal plasticity UMAT" of slide 26 is the flow-rule
subroutine `computeFlowRule` of the UTSA constitutive driver (12 FCC slip
systems, thermally activated slip, γ̇0 = 1e7, ΔG = 9.5e-19 J, p = 0.78,
q = 1.15, T = 1123 K) with hand-coded derivatives of Δγ^α with respect to the
resolved shear stress, back stress, thermal SSD, cross-slip and athermal SSD
resistances; its OTI version `computeFlowRule_otis.f90` (module OTIM6N1) was
written by hand. They live in the git repository
`https://github.com/santiagarcia/OTI_computeflowrule` (local clone
`~/Downloads/driver_utsa`, HEAD `fde918b`) and are **not part of either
repository** (no licence stated); the script builds them from `--source-dir` /
`CLAIM2_FLOWRULE_DIR` and records their SHA-256. Its benchmark compared OTI with
analytical only (`results_*.csv`, `timing_*.csv`: OTI/analytical 1.42× per call).
**No finite-difference script, data or timing record for the slide-27 FD
statements exists anywhere on this machine**, so FD is measured here from scratch.

**Protocol.** One driver (`presentation/claim2_flowrule_driver.f90`) compiled with
the flow-rule sources by one compiler with one flag set for every method: ifort
`-O2` (the original benchmark's flags) and gfortran `-O2`. Accuracy: every active
slip system (12) and variable (5); FD of the analytical routine's primal output,
one slip system perturbed at a time, h = s·max(|z|, |τ^α|), s = 1e-2 … 1e-8.
Timing: `system_clock` around N = 200 000 calls after one untimed call, R = 5
repetitions, median; single thread; the machine was shared (load average ≈ 8).

| Quantity | ifort -O2 | gfortran -O2 | Slide |
|---|---:|---:|---|
| OTI vs hand-coded, max rel. error | 7.1e-15 | 7.2e-15 | "10E−15" |
| centred FD, best step (smooth variables) | 3.2e-9 (ratio 4.5e5) | 3.2e-9 (4.4e5) | "about 10E7 times larger" |
| centred FD, plateau-chosen step | 2.0e-8 (ratio 2.9e6) | 2.0e-8 (2.8e6) | |
| forward FD, best step (all variables) | 1.2e-6 (ratio 1.7e8) | 1.2e-6 (1.6e8) | |
| time: OTI / analytical (with its Jacobian) | 1.49 | 2.25 | |
| time: centred FD, known step (10 primal calls) / OTI | 5.97 | 3.49 | "about 9.00×" |
| time: forward FD, known step (6 calls) / OTI | 3.33 | 2.10 | |
| time: centred FD incl. the 7-step search / OTI | 41.9 | 24.6 | "(getting the correct perturbation)" |

The athermal SSD resistance is 0 on every slip system and enters as
√(ssd² + gnd²) with gnd = 0, i.e. |ssd|: the hand-coded and OTI derivatives are the
one-sided (ssd > 0) slope, which the centred difference cannot reproduce (it
averages the two slopes, error 100 %); the forward difference does. The 1e7 and
9.00× of the slide cannot be tied to a recorded FD protocol; the measured values
bracket 1e7 (centred 4.5e5 … forward 1.7e8), and the FD time is 6.0× OTI with a
known step and 41.9× when the step has to be found.

## Claim 3 – slides 28–32: residual-method sensitivities of a crystal-plasticity flow model

**Model.** The six-parameter thermally activated flow model
`parameter_sensitivity/models/m5_cpflow` (τ0, ΔG, p, q, γ̇0, H;
props 200000, 0.3, 1500, 25, 0.4, 1.6, 0.1, 60000). The study on the slides was
done with a different, 12-slip-system single-crystal UMAT ("the original
constitutive law from the driver", ΔG in joules) whose source, driver and plotting
scripts are **not on this machine**; its figures existed by 2026-06-19, before
m5_cpflow was written (2026-07-20) as a von Mises surrogate with the same six
parameter names. Slide-specific numbers of that model (σ_vM ≈ 2.1 GPa at the end,
its share pattern, ×1.1 cost) therefore cannot be reproduced; the capabilities are
reproduced with m5_cpflow. The slide-31 "4×4×4 mesh" figure is pixel-identical to
the slide-29 single-element figure apart from its title and legend box, so it is
not independent evidence; the mesh claim is re-measured here.

**Set-up.** One C3D8 (unit cube, Abaqus Gauss points, B-bar), bottom face y = 0
fixed, top face u1 = γ·H, u2 = 0, u3 free; γ12 0 → 0.1 (tensor E12 0 → 0.05, the
old study's "5 %") in 50 increments of dt = 1. Nonlinear solve with the ORIGINAL
UMAT; the residual method (`presentation/c3d8_residual.py`: K du/dp = −dR/dp with
the provider's carried DSIGMA_DP) gives all six sensitivities in one enriched run.

| Check | Measured | Slide |
|---|---|---|
| σ_vM at the end | 1582.677 MPa (Abaqus: 1582.677, rel. diff 1.1e-15) | ≈ 2.1 GPa (other model) |
| six sensitivities in one run | yes; max\|du/dp\| 5.8e-19 (homogeneous shear) | yes |
| OTI vs hand-derived chain rule, NRMSE (RMSE / max\|ref\|) | 2.9e-16 worst of six | < 1e-8 |
| chain rule vs FD of its own primal (reference self-check) | NRMSE ≤ 6.2e-11 | |
| central FD of the ORIGINAL (C3D8 re-solved), best step | NRMSE 1.3e-11 … 5.1e-11 | FD "best case, h optimised" |
| FD of Abaqus (26 jobs, double-precision .fil) vs OTI | NRMSE ≤ 6.2e-11 | |
| 4×4×4 mesh (affine field on the boundary, 81 free dof) vs 1 element | 1.0e-14 (Abaqus 4×4×4 vs 1 element: 2.9e-16) | "remain the same" (old 1.0e-14) |
| cost, compiled kernel (8 IPs × 50 incs, same object and flags; 5 samples × 2000 repetitions, median) | plain 1, HYPAD 9.2 (MARCH) / 10.0 (EVAL with carry), central FD 13 | HYPAD "small overhead" (old study figure ×1.1), FD 13 runs |
| cost, whole Python analysis driver | plain 1, HYPAD 2.5, central FD 13 | |

Regimes (plastic fraction dε̄p/(dγ/√3): elastic < 0.05, plastic flow ≥ 0.95×steady,
transition between), measured on a 10× finer path at the same strain rate:

| Slide 29 statement | Measured | Holds? |
|---|---|---|
| elastic regime: sensitivities near zero | Σ\|p dσ_vM/dp\| at increment 1 is 3.5e-7 of its final value | yes |
| yield transition: τ0 becomes more important | τ0 share 15.9 % → 21.3 % | yes |
| plastic flow: ΔG decreases, q grows | ΔG 18.75 → 18.43 %, q 31.27 → 31.41 % (over the whole path 45 → 18 % and 16 → 31 %) | yes |
| H remains small | H share rises to 14.9 % (4th of six at the end) | **no** |
| τ0, ΔG and q dominate | q 31.4, p 23.7, ΔG 18.4, τ0 12.2, H 10.2, γ̇0 4.1 % (mean in plastic flow) | **no** (p, not τ0) |

The two failures are statements about the slide's model; m5_cpflow with the
repository's parameters behaves differently.

## Claim 4 – slide 8: benchmark DDSDDE in Abaqus

**How the slide table was produced.** `~/Documents/UMAT_source_transformation/
validate_all_local.py` transformed every completed benchmark contract, built the
paired Abaqus workspace with DDSDDE **forced** into the compared outputs (the
committed spin_elas_def contract does not request it), ran the original and the
transformed UMAT in Abaqus locally and wrote `umat_oti_workspace/validate_all/
summary.json`; RA commit `3521485` (`results/make_transform_table.py`) rendered the
slide table from it (19 cases; UMAT_VPDCL_R failed to run and was left off). The
ARC run of the same batch (`paper_results/arc_791506`) has the same absolute
differences for the rows it compared; it did not request DDSDDE for spin_elas_def.

**What is run now.** For each `benchmarks/*.json` of the UMAT repository:
`run_config_transform` and `build_validation_workspace` (what
`tools/run_completed_json_batch.py --validate` calls) with STRESS, STATEV,
DDSDDE, CONVERGENCE compared; both Abaqus 2021.HF5 jobs one at a time, named
`claudeP_c4_<case>_{orig,oti}`, `double=both`; a job counts only with the
completion mark (the documented teardown abort is recognised by the repository's
`completed_despite_teardown_abort`); `extract_results`,
`compare_validation_results`. Material: the paired-validation probe vector (unit
constants, 0.3 for Poisson-like names) the job builder writes — as on the slide.
The DDSDDE "max error" is the largest entry-wise |original − transformed| over the
compared increments; the relative column of the current comparator divides by the
largest entry of the increment's matrices (the slide's relative column used the
older entry-wise definition, which is why the absolute values agree and the
relative ones do not).

Two rows that failed on the current code were product gaps and are fixed in the
UMAT repository, so they run from their committed contracts:

* UMAT_HIN: the contract promotes ONE, TWO, ZERO, which the source's helpers set
  by DATA and never assign; the transformer refused DATA-initialised promoted
  names. A DATA-initialised name that nothing assigns is now kept real (it is a
  compile-time constant nothing seeded can reach); an assigned one is still
  refused (`src/umat_oti/transform/source_transform.py`, test
  `tests/test_a_data_constant_in_the_promote_list_stays_real.py`).
* UMAT_PCO: `UMAT_PCO.for` calls KCLEAR, KMMULT, KSMULT, KMTRAN, KMAVEC, KUPDVEC,
  KCLEARV and KMATSUB and defines none of them. The committed contract
  `benchmarks/UMAT_PCO.json` now declares `"dependency_roots": ["../UMATs/UMATs/ICP"]`;
  the transformation service resolves the routine closure (here from the sibling
  `UMAT_ECL_TEMP.for`), writes it entry file first so every line anchor holds,
  transforms it and records the closure; the paired validation compiles the
  original from the same resolved file (`src/umat_oti/services/transformation.py`,
  one line in `tools/run_completed_json_batch.py`, test
  `tests/test_a_contract_resolves_its_helper_closure.py`). The helper sources are
  in the repository: `UMATs/UMATs/ICP/*.for` are, after line-ending
  normalisation, byte-identical to `UMATS/*.for` of
  `https://github.com/jgomezc1/ABAQUS-US` (Juan Gómez, Universidad EAFIT; MIT
  licence, "Copyright (c) 2015 Juan Gomez"), which the upstream itself pairs with
  its UEL files (`UELS/UEL8_PCOR.for` carries the same helpers). MIT permits
  inclusion provided the copyright and permission notice are kept; note that
  `THIRD_PARTY_NOTICES.md` currently describes these files as the authors' own
  under GPL-3.0-only, which does not match their MIT upstream (left for the
  lead: licence text is not this lane's to change).

One row stays a labelled variant, because the defect is in the source:

* UMAT_NKH_1.02: with the probe's PROPS(1) = 1 the source takes THTA = PROPS(1)
  and never sets DTHTA, then uses it in the thermal strain (line 111; a
  `-finit-real=snan` build traps there), so each build computes with whatever
  memory holds (the committed-contract run differs by 2657 in DDSDDE and 0.9 % in
  stress between the two builds). The variant sets PROPS(1) = 0 and gives every
  node the initial temperature 1.0, i.e. the probe's material with DTHTA = 0
  defined.

Four transformer or contract defects were found by this run and fixed in the
UMAT repository (see the end of this file): lost statement labels (PCL, PCLI,
PCLI_R, PCLK did not compile in Abaqus), skipped predictor-stiffness inputs
(VPDCL, NKH wrong stress), DATA constants listed under promote (HIN refused) and
the missing helper closure (PCO refused). Before the first two fixes the same run
gave 10 passes of 16 compared slide cases
(`imq_abaqus/claude_P/claim4_before_transformer_fix.json`).

| UMAT | Committed contract: DDSDDE max abs / rel, verdict | Labelled variant (source defect) | Slide (abs / rel, explanation) |
|---|---|---|---|
| UMAT_ECL_TEMP | 0 / 0 pass | – | 0 / 0, Exact |
| UMAT_ECO | 0 / 0 pass | – | 0 / 0, Exact |
| UMAT_HIN | 0 / 0 pass | – | 0 / 0, Exact |
| UMAT_NKH_1.02 | 2657 / 8.9e-01 FAIL (uninitialised DTHTA) | **0.7234** / 7.5e-04 pass | 0.72 / 2.1e-3, notching |
| UMAT_PCL | 0 / 0 pass | – | 0 / 0, Exact |
| UMAT_PCLI | 0 / 0 pass | – | 0 / 0, Exact |
| UMAT_PCLI_R | 0 / 0 pass | – | 0 / 0, Exact |
| UMAT_PCLK | 0 / 0 pass | – | 0 / 0, Exact |
| UMAT_PCO | 0 / 0 pass (resolved helper closure) | – | 0 / 0, Exact |
| UMAT_VPDCL | **0.2370** / 2.5e-04 pass | – | 0.24 / 1.2e-3, notching |
| UMAT_VPDCO | **0.3716** / 3.6e-04 pass | – | 0.37 / 1.1e-3, notching |
| code_exp | **0.03125** / 1.1e-07 pass | – | 0.031 / 1.6e-7, rounding |
| code_imp | **0.015625** / 6.2e-08 pass | – | 0.016 / 8.2e-8, rounding |
| elastic | 0 / 0 pass | – | 0 / 0, Exact |
| spin_elas_def | **868.3** / 3.1e-03 pass | – | 740 / 2.6e-3, "error in original UMAT" |
| spin_elastic | 0 / 0 pass | – | 0 / 0, Exact |
| visco_beam | 0 / 0 pass | – | 0 / 0, Exact |
| visco_imp | 0 / 0 pass | – | 0 / 0, Exact |
| UMAT_VPDCL_R (not on the slide) | both jobs fail (3-D deck; source in bounds only for NTENS = 4) | 0.3673 / 3.2e-04 pass (plane-strain deck) | – |

**Measured vs slide.** Committed contracts: all 18 slide cases run, **17 pass
(12 exact, 5 differ within tolerance)**; NKH fails on its uninitialised DTHTA.
With NKH's labelled variant row: **18/18 compared and verified, 12 exact, 6
differing within tolerance** — the slide's counts; the absolute differences of
NKH, VPDCL, VPDCO, code_exp and code_imp are the slide's to the digits it shows.
spin_elas_def differs by 868.3 (slide 740): its stress is integrated from DFGRD0/DFGRD1 with a spin correction,
so the OTI tangent carries stress-dependent terms (0.1–0.3 % of the diagonal) that
the hand-coded elastic DDSDDE omits; the shear diagonal differs by 2.6e-4
relative, **not by a factor of two** (the slide's "stray factor of two in a shear
term" is not what this comparison shows for spin_elas_def, and no other benchmark
row shows one). Whether the OTI tangent is the one Abaqus's co-rotational update
needs was not checked against an independent finite-strain reference here.

## Claim 5 – slide 25: internal constitutive Jacobians

**What existed.** No executed record of the 19-entry table exists (the compliance
matrix calls the old Table-3 artefact header-only); `tools/run_internal_jacobian_round.py`
verifies FJAC only, and only for sources with a declared property vector.

**What is measured.** For every (UMAT, symbol) pair present in the source, at the
converged state of the model's own local Newton solve (the repository's probe:
the solve's iterate GAM_PAR is overridden on one increment by a seeded PROPS
slot, `umat_oti.transform.local_jacobian_probe`), the hand-coded value, the OTI
value of d(response)/d(GAM_PAR) from the OTI build of the same injected source, and
centred FD of the ORIGINAL (plateau of a 5-step ladder). The derivative each symbol
stands for is declared: FJAC = dFGAM/dGAM_PAR (discovered from the Newton
update), DETDG = dPHIINV/dGAM_PAR, ANP1P = dANP1/dGAM_PAR, BNP1P = dBNP1/dGAM_PAR,
GDIA = dDIAG/dGAM_PAR (each diagonal entry). CEVPI (HIN) is the inverse of HIN's
local operator and enters DDSDDE = (1 − D)·CEVPI, so it is scored by claim 4's HIN
DDSDDE comparison. Material: the paired-validation probe vector (unit constants,
0.3 for Poisson-like names) with the repository's own material for PCL and PCLK;
declared corrections, each with its reason in the JSON: NTENS = 4 for VPDCL_R,
PCLI, PCLI_R (plane-strain sources), PROPS(1) = 0 for NKH (otherwise DTHTA is read
unset), unit-constant PROPS headroom for the KUHARD table read. Every source first
passes a gate that compiles the ORIGINAL with `-finit-real=snan -ffpe-trap=invalid`
and runs the chosen path (a trap means it computes with unwritten memory).

Cell: slide classification of hand-coded vs OTI (Exact = identical; Pass = relative
difference ≤ 1e-5; otherwise the measured difference), then the relative
differences of OTI and of the hand-coded value against FD, then the slide's cell.

| UMAT | FJAC | DETDG | GDIA | ANP1P | BNP1P | CEVPI |
|---|---|---|---|---|---|---|
| HIN | – | – | – | – | – | Exact (claim 4 DDSDDE, committed contract); slide Exact |
| NKH_1.02 | Pass 7e-14 (FD: OTI 6e-14, hand 1e-13); slide Pass | Exact (FD 4e-14); slide Exact | – | **Differs 1.9e-4** (FD: OTI 4e-11, hand 1.9e-4); slide Exact | **Differs 1.9e-4** (FD: OTI 7e-11, hand 1.9e-4); slide Exact | – |
| PCL | Pass 3.5e-16 (FD 4e-13 both); slide Exact | – | – | – | – | – |
| PCLI | Pass 2.2e-13 (FD 2e-9 both); slide Exact | – | Pass 1.5e-16 (FD 2e-11 both); slide Exact | – | – | – |
| PCLI_R | Pass 2.2e-13 (FD 2e-9 both); slide Exact | – | Pass 1.5e-16 (FD 2e-11 both); not on slide | – | – | – |
| PCLK | Pass 3.5e-16 (FD 7e-13 both); slide Exact | – | – | – | – | – |
| PCO | Pass 2.2e-13 (FD 2e-9 both); slide Exact | – | Pass 1.5e-16 (FD 2e-11 both); slide Exact | – | – | – |
| VPDCL | Pass 1.1e-16 (FD 3e-14 both); slide Pass | Exact (FD 2e-14); slide Exact | – | – | – | – |
| VPDCL_R | **Differs 2.6e-3** (FD: OTI 4e-11, hand 2.6e-3); slide Pass | Exact (FD 8e-14); slide Exact | **Differs 6.7e-5** (FD: OTI 6e-11, hand 6.7e-5); not on slide | – | – | – |
| VPDCO | **Differs 2.6e-3** (FD: OTI 4e-11, hand 2.6e-3); slide Pass | Exact (FD 8e-14); slide Exact | **Differs 6.7e-5** (FD: OTI 6e-11, hand 6.7e-5); slide Pass | – | – | – |

**Measured vs slide.** The sources contain 21 (UMAT, symbol) pairs (the slide
shows 19; PCLI_R and VPDCL_R also carry GDIA). All 21 were measured; the OTI
value agrees with FD in every one (worst 2e-9, the FD reference's own limit).
Against the hand-coded value: 5 exact, 10 within 1e-5 (eight of them at 1e-13 or
below, i.e. one to a few units in the last place, which the slide would have
called Exact), and 6 that differ — every one of the six because the hand-coded
Jacobian is wrong (it disagrees with FD by the same amount): NKH's ANP1P and BNP1P
drop a factor (1 − D), VPDCO's and VPDCL_R's GDIA(3,3) drop (1 − D), and their
FJAC, which uses GDIA, is 2.6e-3 off. With damage D = 0 these would be exact; the
probe path activates damage. The slide's "14 exact, 5 within tolerance" is
therefore not what the current code measures for these hand-coded values; the
statement "all needed derivatives can be computed with our framework" is
supported (OTI = FD everywhere).

## All slides

| Slide | Content | Status |
|---|---|---|
| 1–5 | title, objectives, motivation | no quantitative claim |
| 6 | UMAT inputs/outputs; DDSDDE must be derived by hand | background |
| 7 | transformer: real UMAT + minimal contract → HYPAD UMAT, completed mapping, transform report | exercised by claim 4 (19 benchmark contracts). Actual CLI outputs (`run_config_transform`): `<name>_oti.for` and `<name>_oti_combined.f90` (≙ HYPAD_UMAT.f), `derivative_manifest.json` (≙ completed_mapping.json), `transform_report.json`/`.txt`, `compile_hint.sh`, `compile_order.txt`, support modules `otim<n>n<k>.f90`, `master_parameters.f90`, `real_utils.f90`, `oti_intrinsics.f90`, `umat_oti_helpers.f90`. Owned by the transformer lanes. |
| 8 | 18/18 benchmark DDSDDE | claim 4 |
| 9–12 | developer/collaborator split (REAL_UMAT.obj, OTI_UMAT.obj, Mapping.json, transform_report.txt; Analysis.inp/.odb + sensitivity_request.json → sensitivity_results.json, sensitivity_tables.csv, run_report.txt) | owned by Agent G (GUI), Agent C (cantilever-scale replay engine) and Copilot (`resasm request`). Measured here: `umat-oti-provider build` produces `umat_<model>_oti.obj` + `umat_<model>_oti.json` (≙ OTI_UMAT.obj + Mapping.json) for 20 of the 21 model directories (m2_elastic3d's contract is a different schema); the ORIGINAL UMAT is bundled inside the object (symbol `umat`), not shipped as a separate REAL_UMAT.obj. |
| 13 | sensitivity verification table | claim 1 |
| 14 | section title | – |
| 15, 33, 39 | FCC and J2 cantilevers (384 / 1 536 C3D8, 25 / 40 steps, weighted sensitivities E 96 %, σ_y0 71 %) | owned by Agent C; Abaqus data by the lead in `imq_abaqus/claude_cantilevers/` |
| 16–18, 40–42 | GUI walk-through (transform, provider build, residual solve) | owned by Agent G |
| 19 | 207 UMATs collected by web scraping, full set re-run after every fix | owned by the lead: pass 13 has 391 acquired sources, 244 transformed, 43 of 260 adequately specified genuine UMATs clear every Abaqus acceptance gate; the slide's 207 is an older count |
| 20, 22 | section titles | – |
| 21 | take-home statements | Jacobians generated automatically: claims 2, 4, 5; parameter sensitivities by the residual method: claims 1, 3; source never leaves the developer: slides 9–12 lane (the provider object bundles the ORIGINAL and the OTI lift; `resasm request` consumes object + mapping only); "verified against independent references": every row of this file |
| 23 | OTI: arbitrary order, machine accuracy, step independence | background; demonstrated by RA `tests/framework/test_otilib_adapter.py` (exact polynomial coefficients to order 2), `tests/framework/test_otilib_fe_sensitivity.py` (third-order C3D8 force derivatives against analytic values, 1e-13), UMAT `tests/test_actual_umat_higher_order.py` (orders 2–4 of transformed UMATs against independent references), and claims 2–3 here (OTI at 1e-15 vs FD whose error depends on the step) |
| 24 | Jacobians in FEM (residual, Newton, B-matrix) | background |
| 25 | internal constitutive Jacobians | claim 5 |
| 26–27 | CP flow rule: OTI vs analytical vs FD | claim 2 |
| 28–32 | CP residual-method sensitivities | claim 3 |
| 34–38 | motivation (repeated) | background |

## Changes to existing files

See `docs/evidence/claude_P.md`. In the UMAT repository four defects found by
claim 4 were fixed, each with a regression test that fails without the fix:

| Commit | File | Defect |
|---|---|---|
| `884d39e` | `src/umat_oti/transform/source_transform.py` | labelled statements lost their label (PCL, PCLI, PCLI_R, PCLK did not compile) |
| `884d39e` | `src/umat_oti/fortran/regions.py` | inputs of a DDSDDE used as the predictor stiffness were skipped (VPDCL, NKH wrong stress) |
| `1cd2e58` | `src/umat_oti/transform/source_transform.py` | a DATA constant listed under promote refused the file (HIN) |
| `0b075b4` | `src/umat_oti/services/transformation.py`, `benchmarks/UMAT_PCO.json`, `tools/run_completed_json_batch.py` | a contract could not name the published helper sources its UMAT calls (PCO) |

Tests: `tests/test_benchmark_transforms_keep_labels_and_predictor_inputs.py`,
`tests/test_a_data_constant_in_the_promote_list_stays_real.py`,
`tests/test_a_contract_resolves_its_helper_closure.py`. The transformer edits move
the transform fingerprint, so
`tests/test_contract_fixtures.py::test_the_recorded_generation_is_this_worktrees_actual_transform`
fails until the lead re-freezes `transform_generation.json` (reserved to the lead).
