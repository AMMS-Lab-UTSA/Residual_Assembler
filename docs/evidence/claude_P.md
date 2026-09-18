# Evidence: Claude agent P, presentation claims (2026-09-18)

Branch `claude/P-2026-09-18` in both repositories (worktrees `claude-ra-P`,
`claude-umat-P`), from the snapshot commits RA `d8ea804`, UMAT `f6fcc42`. Mission:
make the quantitative claims of the IMQCAM Annual Meeting deck (V4, 11 August
2026) reproducible from the repositories, each with one command, an independent
reference, a machine-readable result and documentation. The claim-by-claim record
is `docs/PRESENTATION_CLAIMS.md`; this file records how the work was done, what
changed in existing files and what remains.

## Environment of the recorded run

Python 3.11.7 (`/home/ammslab3/softwarex_work/.venv`), NumPy 2.4.6, gfortran 9.4.0,
ifort 2021.10.0 (oneAPI 2023.2.1), Abaqus 2021.HF5. Imports verified to come from
the two worktrees (every claim script refuses to run otherwise and records the
import paths in its JSON). Large intermediates were written to
`/home/ammslab3/softwarex_work/imq_abaqus/claude_P/` (never inside a worktree).
The machine was shared with other agents during the run (load average 8–21 on 24
cores); timing ratios are measured with every method under the same load, and
the load average is recorded beside them.

Abaqus: one job at a time, every job named `claudeP_*`, judged by
"THE ANALYSIS HAS COMPLETED SUCCESSFULLY" in the .sta. Jobs run: claim 3, 26 per
run (nominal, 24 perturbed, 4×4×4), two runs; claim 4, 2 per case (19 cases + the
labelled variants) over three full runs and the NKH variant reruns.

## Recovered provenance of each slide

* Slide 13: RA `origin/cross-platform-hardening:residual-assembler/results/build_sweep_tables.py`,
  `sweep_extra.py`, `figures/sweep_error_tables.json` (20 materials; the 18 on the
  slide are those below 1e-5; the old "Program 2 on a single C3D8" was one
  integration point).
* Slide 8: `~/Documents/UMAT_source_transformation/validate_all_local.py` →
  `umat_oti_workspace/validate_all/summary.json` → RA commit `3521485`
  `results/make_transform_table.py`.
* Slides 26–27: `computeFlowRule.f90` / `computeFlowRule_otis.f90` of
  `https://github.com/santiagarcia/OTI_computeflowrule` (local `~/Downloads/driver_utsa`,
  HEAD `fde918b`); the benchmark compared OTI with analytical only; no FD script,
  data or timing record exists on this machine.
* Slides 28–32: the original 12-slip "Constitutive_UMAT.f" study (figures dated by
  2026-06-19 in `IMQCAM_BZM_project_presenation_2026024_V9.pptx`) is not on this
  machine; RA `results/cp_residual_sensitivities.py`, `cp_mesh_residual.py`,
  `cp_analysis_timing.py` (July 2026) re-created it with the surrogate m5_cpflow;
  the slide-31 "4×4×4" figure is pixel-identical to the slide-29 single-element one.
* Slide 25: no executed record of the 19-entry table (compliance matrix: the old
  Table-3 artefact is header-only).

## New files

RA (`claude-ra-P`):

| File | Purpose |
|---|---|
| `presentation/__init__.py`, `common.py` | provenance, provider build (`python -m umat_oti.provider build`), FD ladder and plateau |
| `presentation/claim1_sensitivity_sweep.py` | slide 13 |
| `presentation/claim2_flowrule_jacobian.py`, `claim2_flowrule_driver.f90` | slides 26–27 |
| `presentation/claim3_cp_residual_c3d8.py`, `c3d8_residual.py`, `cpflow_analytic.py`, `cp_timing_driver.f90`, `claim3_abaqus.py` | slides 28–32 (generic small-strain C3D8 residual-method driver with nonzero prescribed displacements; hand-derived chain rule of m5_cpflow; compiled timing driver; Abaqus deck + ASCII .fil reader) |
| `presentation/claim4_benchmark_ddsdde.py` | slide 8 (paired Abaqus, claudeP job names, documented variants) |
| `presentation/claim5_constitutive_jacobians.py` | slide 25 (all six internal-Jacobian symbols; uninitialised-read gate) |
| `presentation/run_all.py` | runs all claims, writes `presentation_summary.json` |
| `presentation/expected/presentation_summary.json` | committed summary of the recorded run |
| `tests/presentation/*` | one end-to-end test module per claim (conftest + 5 modules) |
| `docs/PRESENTATION_CLAIMS.md`, `docs/evidence/claude_P.md` | documentation |

UMAT (`claude-umat-P`): `docs/PRESENTATION_CLAIMS.md`,
`tests/test_benchmark_transforms_keep_labels_and_predictor_inputs.py`,
`tests/test_a_data_constant_in_the_promote_list_stays_real.py`,
`tests/test_a_contract_resolves_its_helper_closure.py`.

## Edits to existing files (future merge points)

| Repository | File | Change | Why |
|---|---|---|---|
| RA | `.gitignore` | `presentation/results/` ignored | run-time results |
| UMAT | `src/umat_oti/transform/source_transform.py` | new `_logical_line_prefix`, `_statement_label`, `_restore_statement_label`; the three logical-line helpers use the prefix; the branch and assignment rewrite sites restore the label | a labelled IF reading a promoted variable (`  802 IF (IFLAG.EQ.1) THEN`) lost its label and started in the label field: UMAT_PCL, PCLI, PCLI_R, PCLK did not compile in Abaqus |
| UMAT | `src/umat_oti/fortran/regions.py` | in `_dependency_summary`, inputs of DDSDDE writes that precede the first stress use of DDSDDE count as stress inputs | ELAM/EBULK3 were classed tangent-only and skipped while the kept predictor block read them: UMAT_VPDCL and UMAT_NKH_1.02 returned a stress off by a hydrostatic offset |

| UMAT | `src/umat_oti/transform/source_transform.py` (`1cd2e58`) | in `_roles_with_stress_path_promotions`, a DATA-initialised name that nothing assigns is moved from `promote` to `constant` (10 lines) | UMAT_HIN's committed contract promotes ONE, TWO, ZERO (DATA constants of its helpers) and the DATA blocker refused the file |
| UMAT | `src/umat_oti/services/transformation.py` (`0b075b4`) | new `_payload_with_resolved_closure`: a compact contract may declare `dependency_roots`; the closure is resolved, written entry-first to `<out>/<stem>_resolved.<ext>` and transformed; recorded as `dependency_closure` in the summary | UMAT_PCO.for calls helpers it does not define |
| UMAT | `benchmarks/UMAT_PCO.json` (`0b075b4`) | `"dependency_roots": ["../UMATs/UMATs/ICP"]` | as above |
| UMAT | `tools/run_completed_json_batch.py` (`0b075b4`) | the original UMAT of the paired validation is the resolved closure when the summary has one (one statement) | the original needs the helpers to link in Abaqus |

All UMAT transformer edits move the transform fingerprint. The lead owns
`transform_generation.json`; until it is re-frozen,
`tests/test_contract_fixtures.py::test_the_recorded_generation_is_this_worktrees_actual_transform`
fails (it passes on the snapshot without these edits). Copilot edits the same two
files in `imq-umat-recovery`; the edits are small and self-contained.

## Findings about the sources (not changed; they are inputs)

* UMAT_NKH_1.02 uses DTHTA unset when PROPS(1) ≠ 0 (line 111; trapped by a
  `-finit-real=snan -ffpe-trap=invalid` build).
* UMAT_VPDCL_R: `DO K1=10,2*NTENS+5` into `XBACK(NTENS)` is in bounds only for
  NTENS = 4 (`-fcheck=bounds` stops at XBACK(7) with NTENS = 6); the benchmark
  contract declares 6.
* UMAT_PCLI, PCLI_R: plane-strain sources (shear stiffness for component 4 only);
  at NTENS = 6 the probe returns NaN state.
* UMAT_PCO calls helpers it does not define; its family defines them in sibling
  files and the committed contract now names that root (fixed, above).
* Licence: `UMATs/UMATs/ICP/*.for` are, after line-ending normalisation,
  byte-identical to `UMATS/*.for` of the MIT-licensed
  `https://github.com/jgomezc1/ABAQUS-US` ("Copyright (c) 2015 Juan Gomez",
  Universidad EAFIT); `THIRD_PARTY_NOTICES.md` describes them as the authors' own
  under GPL-3.0-only. MIT allows their inclusion with its notice kept; the notice
  text is left for the lead.
* Hand-coded internal Jacobians wrong when damage is active (confirmed by FD, OTI
  agrees with FD): NKH ANP1P/BNP1P (1.9e-4), VPDCO and VPDCL_R GDIA(3,3) (6.7e-5)
  and FJAC (2.6e-3).
* The probe vector's Poisson-name rule misses `ENU=MIN(PROPS(2),ENUMAX)` (PCLI,
  PCLI_R, PCO), so PROPS(2) = 1.0 is capped near 0.5.

## Tests (counts from JUnit/pytest output of the recorded runs)

| Suite | Command | Result |
|---|---|---|
| RA `tests/presentation`, offline incl. `slow` | `pytest -q tests/presentation -m "not abaqus"` | 7 passed (4 min 53 s) |
| RA `tests/presentation`, Abaqus | `pytest -q tests/presentation -m abaqus` | 2 passed (claim-3 nominal job, claim-4 `elastic` pair) |
| UMAT new regression file | `pytest -q tests/test_benchmark_transforms_keep_labels_and_predictor_inputs.py` | 4 passed; the same 4 fail without the two fixes |
| UMAT new regression files (HIN, PCO) | `pytest -q tests/test_a_data_constant_in_the_promote_list_stays_real.py tests/test_a_contract_resolves_its_helper_closure.py` | 8 passed; 3 of the 4 DATA tests fail without the change (the fourth checks the assigned case stays refused) |
| UMAT offline suite, with all four fixes | same command | 3327 passed, 3 failed, 125 skipped (same three as below) |
| UMAT offline suite, with the first two fixes | `pytest -q -m "not abaqus and not corpus_pass and not network and not browser and not arc"` | 3315 passed, 3 failed, 125 skipped: `test_the_recorded_generation_is_this_worktrees_actual_transform` (the fixes move the fingerprint; re-freeze reserved to the lead) and two `test_repository_standards.py` checks that fail on files this agent did not touch (absolute home paths and stale references in `docs/PROVIDER.md`, `docs/evidence/recovery_*.md`) |
| RA offline suite | `pytest -q -m "not abaqus and not arc and not network" --continue-on-collection-errors` | 299 passed, 57 failed, 19 skipped, 5 errors; 56 failures and the 5 errors are `FixtureError` (the verified fixtures are frozen at transform fingerprint `94a92c01814f107a`, `schemas/transform_generation.json` of the snapshot says `6aa20d22e37f14c9`), present in the snapshot and unrelated to this work; the remaining failure was the claim-5 test's FD threshold, fixed (noise-aware plateau) and passing in the presentation run above. Without `--continue-on-collection-errors` the suite stops at the collection error of `tests/framework/test_a_verified_deck_drives_the_global_assembly.py` (same cause) |

## Measured values of the recorded run

`presentation/expected/presentation_summary.json` (and the per-claim JSON in
`presentation/results/`). Headline numbers: claim 1, 20/20 models, 84 directions,
worst 1.56e-7; claim 2, OTI vs hand-coded 7.1e-15, FD/OTI error 4.5e5–1.7e8, FD
time 6.0× OTI with a known step (41.9× with the step search); claim 3, OTI vs
chain rule 2.9e-16, mesh 1.0e-14, Abaqus 1.1e-15, FD of Abaqus 6.2e-11, HYPAD
9.2–10.0× plain vs FD 13×; claim 4, 17 of 18 slide cases pass from the committed contracts (12
exact, 5 within tolerance), 18/18 with NKH's labelled source-defect variant (12
exact, 6 within tolerance); claim 5, 21 pairs,
OTI = FD everywhere, hand-coded: 5 exact, 10 within 1e-5, 6 wrong.
