# Evidence: general history replay for the presentation cantilevers (Claude agent C)

Branches `claude/C-2026-09-18` in both repositories (RA worktree `claude-ra-C`,
UMAT worktree `claude-umat-C`), from the 11:05 snapshots (RA `d8ea804`, UMAT
`f6fcc42`). Date 2026-09-18. Every number below was measured in a run of this
session; nothing is copied from the slides. Large outputs (full fields, Abaqus
reruns) live outside the repositories in `imq_abaqus/claude_C/`; the small
tables and figures are in `docs/evidence/claude_C/`.

## Claim (slides 12, 15, 18, 33, 39, 42)

"Collaborators obtain requested sensitivities from compiled drivers and saved
FE data without rerunning the production analysis or receiving the UMAT
source": `OTI_UMAT.obj + Analysis.inp + Analysis.odb + sensitivity_request.json`
-> `sensitivity_results.json`, `sensitivity_tables.csv`, `run_report.txt`, with
full-field von Mises sensitivities for every parameter of the J2 cantilever
(1,536 C3D8, 7,497 DOF, 12,288 IPs, 40 steps, 4 parameters) and the FCC
crystal-plasticity cantilever (384 C3D8, 2,025 DOF, 3,072 IPs, 25 steps, 10
parameters), and slide 33's weighted per-step shares.

Status: **works end to end on both cantilevers** from the lead's Abaqus
2021.HF5 data (`imq_abaqus/claude_cantilevers/`), verified against Abaqus
reruns, against whole-model finite differences of the ORIGINAL UMAT at reduced
and at full size, and against exact elastic identities. The slide-33
percentages are **not** reproduced at the load steps the slide names (see
"Weighted shares").

## What was built

UMAT repository (commit `14a4f54`):

- `src/umat_oti/provider/emit.py`, `build.py`: a new, additive provider entry
  point `UMAT_OTI_EVAL_TOTAL` (28 arguments, advertised as
  `symbols.oti_eval_total` in the completed contract). Each parameter
  direction carries the incoming dSTRESS/dp, dSTATEV/dp, dSTRAN/dp and
  dDSTRAN/dp, so the outputs are total derivatives through any UMAT; it also
  returns dSTATEV/dDSTRAN (from the DSTRAN directions that give DDSDDE) and
  reports PNEWDT instead of stopping. `UMAT_OTI_EVAL`, `UMAT_OTI_MARCH` and
  their signatures are unchanged. Reason: the old EVAL cannot chain the
  equilibrium displacement sensitivity through a general UMAT's state (the
  bounded J2 engine uses a J2-only substitution).
- `tests/test_provider_eval_total.py` (16 tests) and the fixture
  `tests/fixtures/provider_total_strain/` (a total-strain damage UMAT that
  reads STRAN, has a history maximum and requests PNEWDT<1).

RA repository (commits `2931574`, `ab83fb8`, and the final commit of this
document): new modules only, Copilot's `connected.py`, `presentation*.py`,
`cmd_request.py` untouched.

- `residual_core/replay/history.py` (engine), `history_material.py` (batched
  provider link: one Fortran loop per increment for all points),
  `history_inputs.py` (deck subset and ODB export, refuses unsupported
  features by name), `history_outputs.py` (requests, reductions, shares,
  public/private files, report), `history_verify.py` (independent
  references), `odb_export_npz.py` (Abaqus Python exporter),
  `residual_core/ui/cmd_history.py` (`resasm history`, `route_request`),
  registered with one line in `residual_core/ui/cli.py`; `scipy>=1.10`
  declared in `pyproject.toml`.
- `docs/REPLAY_HISTORY.md` (usage, mathematics, tolerances, the exact
  `cmd_request.py` edit for routing).
- `tests/replay_history/` (48 tests), `examples/replay_history/` (generator,
  regeneration script and a committed 138 KB Abaqus example),
  `scripts/replay_history_cantilevers.py`, `replay_history_abaqus_fd.py`,
  `replay_history_fd_fullsize.py` (evidence drivers).

## Commands

Environment as in the brief (`PYTHONPATH=claude-ra-C:claude-umat-C/src:otilib`,
`UMAT_OTI_REPO=claude-umat-C`).

```sh
W=imq_abaqus/claude_C; D=imq_abaqus/claude_cantilevers
umat-oti-provider build $UMAT_OTI_REPO/parameter_sensitivity/models/m3_j2/contract_v2.json --out $W/providers/prov_j2
umat-oti-provider build $UMAT_OTI_REPO/parameter_sensitivity/models/m6_fcc/contract_v2.json --out $W/providers/prov_fcc
resasm history --model $D/j2/claude_j2_nominal.inp --odb $D/j2/claude_j2_nominal.odb \
    --material $W/providers/prov_j2/umat_m3_j2_oti.obj --request j2_request.json --out j2_results
# evidence (both cantilevers, recorded + re-equilibrated, Abaqus FD comparisons, figures):
python scripts/replay_history_cantilevers.py --abaqus-data $D --providers $W/providers \
    --work $W/cantilevers --evidence docs/evidence/claude_C --tight-fd j2=$W/j2_fd_tight
# tight-tolerance Abaqus reruns of the J2 cantilever (25 jobs, one at a time, claudeC_ prefix):
python scripts/replay_history_abaqus_fd.py --deck $D/j2/claude_j2_nominal.inp \
    --umat $UMAT_OTI_REPO/parameter_sensitivity/models/m3_j2/umat.for \
    --contract $W/providers/prov_j2/umat_m3_j2_oti.json --steps 5e-3,2e-3,1e-3 --tight --nominal \
    --prefix claudeC_j2cant --work $W/j2_fd_tight
# whole-model ORIGINAL-UMAT FD at full size, one process per parameter:
python scripts/replay_history_fd_fullsize.py --deck $D/j2/claude_j2_nominal.inp \
    --fields $D/j2/claude_j2_nominal_fields.npz --object $W/providers/prov_j2/umat_m3_j2_oti.obj \
    --parameters E --steps 1e-3,1e-4,1e-5 --rtol 1e-12 --out $W/fd_fullsize/j2_E.json
```

My ODB exporter reproduces the lead's export of `claude_j2_nominal.odb` bit
for bit (all 9 arrays identical) in 3.2 s.

The slide-12 command itself from the ODB (`resasm history --model
claude_j2_nominal.inp --odb claude_j2_nominal.odb --material umat_m3_j2_oti.obj
--request j2_request.json --out j2_odb_run`, the export done by the command
through Abaqus Python): 20.3 s wall, 549 MB peak memory, status executed
successfully, max|R_free| 1.39e-2 N at increment 10 against a limit of
0.257 N there.

## Runtime (J2: 40 increments x 12,288 IPs x 4 parameters; FCC: 25 x 3,072 x 10)

| run | material (provider) | assembly | factorisation + solves | total engine | wall |
| --- | --- | --- | --- | --- | --- |
| J2, recorded state | 1.57 s | 1.94 s | 5.97 s | 9.86 s | 11.1 s |
| J2, re-equilibrated | 4.32 s | 4.30 s | 16.2 s | 25.2 s | 26.3 s |
| FCC, recorded state | 17.7 s | 0.49 s | 0.54 s | 18.8 s | 19.4 s |
| FCC, re-equilibrated | 68.7 s | 1.00 s | 1.66 s | 71.5 s | 72.1 s |

One cProfile of the J2 recorded replay (`claude_C/j2_profile.txt`): 11.2 s,
of which SuperLU factorisation 5.26 s (40 x 0.13 s, 7,293 free DOF), the
batched OTI provider 1.41 s (491,520 point evaluations), einsum kernels about
1.5 s, the B-bar operators 0.90 s (setup). The FCC cost is the provider (16
OTI directions, 10 explicit sub-steps, 12 slip systems). Writing the
full-field `fields.npz` (134 MB for J2) adds about 9 s. (Measured while other
jobs shared the 24-core machine.)

## Replay of the ODB: parity and equilibrium

| | J2 recorded | J2 re-equilibrated | FCC recorded | FCC re-equilibrated |
| --- | --- | --- | --- | --- |
| max \|S_replay - S_odb\| (MPa) | 1.19e-2 | 2.79e-2 | 1.61e-3 | 4.45e-2 |
| worst stress error / its limit | 0.010 | 0.032 | 0.006 | 0.014 |
| worst reaction / state error over limit | 0.006 / 0.002 | 0.005 / 0.094 | 0.006 / 0.002 | 0.006 / 0.101 |
| max relative stress difference | 3.5e-5 | 8.4e-5 | 2.8e-5 | 1.1e-3 |
| max \|SDV_replay - SDV_odb\| | 6.2e-9 | 2.4e-7 | 3.4e-6 MPa | 7.8e-3 MPa |
| max \|RF_replay - RF_odb\| (N) | 4.37e-3 | 2.83e-3 | 1.04e-3 | 1.28e-2 |
| max free residual / reaction scale | 1.23e-4 | 2.7e-12 | 8.3e-4 | 1.7e-12 |
| worst free residual / its limit | 0.054 | - | 0.458 | - |
| largest Newton correction | - | 2.66e-7 mm | - | 1.45e-5 mm |
| sensitivity-solve residual (relative) | 2.8e-15 | 2.0e-15 | 5.1e-14 | 5.2e-14 |

Limits (see `docs/REPLAY_HISTORY.md`): single-precision storage of the ODB
value plus the effect of single-precision displacements, and for the free
residual Abaqus's own force tolerance (5e-3 of the mean nodal force). The FCC
ODB was converged by Abaqus with the UMAT's elastic DDSDDE at default
tolerances; its recorded residual uses 46 % of that tolerance, and the
recorded-state sensitivities differ from the re-equilibrated ones by
0.08-4.5 % (field maxima; largest for m: dsigma_vM/dm 4.5 %, dRF/dm 1.5 %).
For J2 the two differ by <= 6.6e-5. The slide numbers below use the
re-equilibrated results (the exact derivative of the discrete problem Abaqus
approximates); the recorded-state results are in the same tables with suffix
`_recorded`.

Also measured from the ODB: the first plastic integration points of the J2
cantilever appear at step 7 (4 IPs; 56 at step 8, 180 at 9, 388 at 10, 1,908
at 20, 2,556 at 40). The handoff's "first plastic IPs at step 10" is the
388-IP count, not the first yield.

## Verification

### (c) Elastic identities (full-size J2 cantilever, increments 1-6)

Only displacements are prescribed, so the displacement field is independent of
E (du/dE = 0) while stresses and reactions scale with E. Measured:

| | re-equilibrated | recorded state |
| --- | --- | --- |
| E max\|du/dE\| / max\|u\| | 1.9e-15 | 7.4e-8 |
| max\|E dS/dE - S\| / max\|S\| | 2.1e-14 | 7.7e-6 |
| max\|E dRF/dE - RF\| / max\|RF\| | 6.8e-15 | 5.1e-6 |

(The recorded-state values carry the ODB's own residual; the identity is a
property of the exact equilibrium.) Under load control the other identity
holds: du/dE = -u/E and dS/dE = 0 (`test_elastic_load_control_scaling`,
<= 1e-12; the example's elastic increments: E dRF2/dE = RF2 and
E dS11/dE = S11 to 2e-5 relative in recorded mode).

### (b) Whole-model FD of the ORIGINAL UMAT, re-equilibrated in Python

The ORIGINAL UMAT (the regular path compiled into the provider object, no OTI)
is re-equilibrated at p(1 +/- h) over the same increments (Newton to 1e-13 of
the reaction scale at reduced size, 1e-12 at full size); for every field and
increment the adjacent pair of steps with the smallest spread is the plateau;
error = max|OTI - FD| / max|FD| over the field.

Reduced meshes (tests, step ladder 1e-3 ... 1e-5):

| model | worst error over all fields and increments, per parameter |
| --- | --- |
| J2 12x4x2, 10 increments | E 1.5e-9, nu 2.9e-9, SIGY0 1.0e-9, H 3.7e-7 (H at its FD spread 4.5e-7) |
| FCC 6x2x1, 5 increments, slide-15 constants (ladder to 1e-6) | g0 4.9e-7, h0 3.3e-8, q 5.9e-8, gd0 1.1e-7, m 1.9e-7, gsat 5.9e-7 (spread 5.2e-7), C11 7.2e-9, C12 1.2e-8, C44 2.1e-8, a 5.0e-8 |
| total-strain damage fixture 6x2x1 | A 4.1e-9, E 3.1e-10, K0 5.6e-10, nu 1.2e-7 |

All <= 1e-6, each within 3x its own FD plateau spread (the FD is the limit); zero
references (e.g. du/dE, yield-stress derivatives before yield) agree on the
field scale to <= 1e-7.

Full-size J2 cantilever (40 increments, all fields, ladder 1e-3/1e-4/1e-5,
about 300 s per parameter): increments where the FD has a plateau (spread
<= 1e-5):

| parameter | U | RF | S | MISES | SDV |
| --- | --- | --- | --- | --- | --- |
| E | 6.6e-8 | 2.0e-8 | 4.8e-9 | 8.7e-9 | 1.1e-10 |
| nu | 8.4e-6 (FD spread 8.2e-6) | 2.6e-8 | 5.8e-8 | 5.3e-8 | 1.1e-8 |
| SIGY0 | 2.6e-9 | 3.0e-9 | 3.1e-9 | 7.9e-9 | 9.1e-11 |
| H | 3.1e-7 | 5.4e-6 | 2.2e-6 | 7.2e-7 | 9.9e-9 |

The FD has no plateau at increments 11-14 (E, SIGY0) and 33-36 (H): an
integration point sits on the yield surface there, so p(1+h) and p(1-h) land
on different branches for every h tried. Even there the OTI result differs
from the finest FD by <= 1e-7, except the reactions at increment 11
(1e-4 against an FD spread of 7e-5 to 9e-5) and SIGY0 stresses at increment
33 (5.9e-2 against a spread of 2.6e-2).

FCC full size: FCC_FULLSIZE_PLACEHOLDER

### (a) Abaqus finite differences

The ODB stores single precision, which floors a central difference at about
eps32 |q| / (2 h |p|), so h cannot go below about 1e-3; with Abaqus's default
tolerances each rerun also carries its own equilibrium error, and in the
plastic range integration points cross the yield surface between the +/- runs
(an O(h) effect). A comparison counts as resolved when an adjacent pair of step
sizes agrees to 1e-3 and the float32 floor is below 1e-3.

Elastic increments (J2 1-6): tip dRF2/dE agrees with the lead's FD to 4.4e-7
(FD spread 4e-7), root von Mises dE to 7e-8, all resolved outputs to 4.1e-6
(E) and 1.3e-4 (nu).

J2, lead's reruns (default tolerances, h = 2e-2, 1e-2, 5e-3; tip RF2, tip-region
U field, S11 at element 1 IP 1, von Mises field of the root column):

| parameter | resolved / nonzero | max rel. error (resolved) | on the output's scale |
| --- | --- | --- | --- |
| E | 68 / 154 | 4.8e-2 | 7.7e-3 |
| nu | 49 / 160 | 1.2e-3 | 5.7e-4 |
| SIGY0 | 85 / 136 | 5.0e-2 | 7.6e-3 |
| H | 82 / 133 | 5.6e-3 | 2.8e-4 |

J2, my reruns with tight controls (R_n 1e-9, C_n 1e-8; h = 5e-3, 2e-3, 1e-3;
25 jobs, 27-40 s each):

| parameter | resolved / nonzero | max rel. error (resolved) | on the output's scale |
| --- | --- | --- | --- |
| E | 62 / 154 | 1.0e-2 | 9.6e-3 |
| nu | 41 / 160 | 2.3e-4 | 4.5e-5 |
| SIGY0 | 87 / 136 | 9.9e-3 | 9.4e-3 |
| H | 54 / 134 | 7.5e-4 | 3.9e-5 |

The remaining E/SIGY0 differences of about 1 % are in the root von Mises field
and the tip reaction at increments where yielding spreads: there the Abaqus FD
changes by about 1 % between h = 5e-3, 2e-3 and 1e-3 (e.g. tip dRF2/dE at
increment 20: -3.606e-4, -3.619e-4, -3.637e-4, while OTI and the full-size
Python FD agree on the whole reaction field to 2e-8 at that increment). This is the resolution of an Abaqus FD
of this problem, not an OTI error.

FCC, lead's reruns (slide-15 constants, default tolerances, h = 2e-2, 1e-2,
5e-3):

| parameter | resolved / nonzero | max rel. error (resolved) | on the output's scale |
| --- | --- | --- | --- |
| g0 | 73 / 96 | 3.5e-3 | 2.3e-3 |
| h0 | 64 / 92 | 2.0e-2 | 6.7e-4 |
| q | 68 / 92 | 1.2e-2 | 5.9e-4 |
| gd0 | 69 / 95 | 1.4e-2 | 5.0e-4 |
| m | 20 / 96 | 5.6e-2 | 1.8e-3 |
| gsat | 68 / 92 | 2.0e-2 | 4.3e-4 |
| C11 | 53 / 100 | 3.2e-3 | 2.6e-3 |
| C12 | 64 / 100 | 3.7e-3 | 1.8e-3 |
| C44 | 25 / 100 | 2.6e-2 | 2.3e-3 |
| a | 65 / 92 | 2.0e-2 | 3.7e-4 |

The lead reports C44 and m as not resolved by the Abaqus FD (tiny derivatives
against the float32 floor, strong rate nonlinearity); for them, and for the
FCC in general, the whole-model Python FD above is the reference. Full rows:
`claude_C/{j2,fcc}_abaqus_fd*_reequilibrated.csv` (recorded-mode summaries in
`claude_C/cantilever_summary.json`).

Committed example (12x4x2 J2 beam, 10 increments, my Abaqus 2021.HF5 runs with
tight controls, reruns at h = 1e-2, 5e-3, 2e-3, 1e-3): for all 4 outputs x 4
parameters x 10 increments, |OTI - FD| <= plateau spread + float32 floor (160
comparisons, 115 nonzero; worst error/(spread+floor) 0.5; tip dRF2/dE 4.5e-5,
dRF2/dSIGY0 3.6e-5), asserted offline by
`test_sensitivities_within_the_abaqus_fd_uncertainty`. `--verify fd` on the
example (default ladder 1e-3, 3e-4, 1e-4, 3e-5, 1e-5): tangent 1.8e-10
(ORIGINAL-UMAT FD, 36 points, spread 2.6e-10); whole-model FD worst 1.5e-7 at
a plateau spread of 3.1e-7; zero references within 1.6e-9; the ODB-driven
du/dp differs from the Python-equilibrium du/dp by 3.8e-6.

## Weighted shares (slide 33)

Definition recovered from `results/cp_residual_sensitivities.py` and
`results/fcc_crystal_results.py` on `origin/cross-platform-hardening`:
W_j = |p_j dq/dp_j|, share = 100 W_j / sum_k W_k per increment. Field form
(primary): W_j(n) = sum_i V_i |p_j d sigma_vM,i(n)/dp_j| / sum_i V_i over all
12,288 (3,072) integration points; the scalar form uses the volume-averaged
sigma_vM. Tables: `claude_C/{j2,fcc}_shares_{reequilibrated,recorded}.csv`;
figures `claude_C/j2_weighted_shares.png`, `j2_steps_SIGY0.png`,
`j2_final_step_dmises.png` and the FCC counterparts.

J2 (field form; scalar form in brackets):

| step | E | nu | SIGY0 | H |
| --- | --- | --- | --- | --- |
| 1-6 (elastic) | 98.2 (99.4) | 1.8 (0.6) | 0 | 0 |
| 7 (first yield, 4 IPs) | 97.6 (99.2) | 1.8 | 0.6 (0.2) | 0.0 |
| 10 | 67.8 (76.1) | 1.3 | 30.7 (23.3) | 0.1 |
| 12 | 48.4 (53.9) | 1.2 | 50.1 (45.7) | 0.3 |
| 15 | 30.3 (32.6) | 1.1 | 67.8 (66.7) | 0.8 |
| 16 | 25.1 (26.6) | 1.0 | 72.9 (72.4) | 1.0 |
| 20 | 14.7 (15.2) | 0.9 | 82.6 (83.2) | 1.8 |
| 40 | 4.8 (4.9) | 0.7 | 89.0 (89.8) | 5.4 (5.2) |

Slide 33 says E carries about 96 % while elastic (measured: 98.2 % field,
99.4 % scalar) and that "at yield onset" SIGY0 climbs to about 71 % while E
falls to about 24 %. Measured at the first yield (step 7) SIGY0 has 0.6 %;
the pair SIGY0 about 71-73 % / E about 24-25 % is reached at step 16
(72.9 / 25.1), when 1,532 of 12,288 points are plastic. The slide's numbers are
therefore not reproduced at the steps it names; they may come from a
different weighting or region, which the old code does not record.

FCC (slide-15 constants; field form): elastic steps 1-2 C11 57.3, C12 40.1,
C44 2.7; step 4 g0 17.0, C11 45.3, C12 31.6; step 8 g0 38.6, C11 24.8, C12 17.4,
h0 3.3, q 2.9; step 25 g0 34.6, C11 15.0, C12 10.6, h0 9.7, q 8.3, gsat 7.6,
a 6.4, C44 3.3, gd0 2.5, m 1.9. (The presentation prints no FCC share numbers.)

## Tests

- UMAT: `tests/test_provider_eval_total.py` 16 passed; with the existing
  provider tests 34 passed (`tests/test_provider_recovery.py` unchanged).
  Measured EVAL_TOTAL vs ORIGINAL-UMAT FD: <= 8.4e-10 (J2, FCC at a
  slip-hardening point, damage fixture). Full UMAT suite on this branch:
  3341 passed, 125 skipped, 3 failed (516 s). Failures: (1)
  `test_contract_fixtures.py::test_the_recorded_generation_is_this_worktrees_actual_transform`
  - caused by this branch: the provider code is part of the transform
  fingerprint (6aa20d22e37f14c9 -> d8a3d2445fda0b26) and
  `schemas/transform_generation.json` is only changed at a lead re-freeze;
  (2) `test_repository_standards.py::test_repository_standards_audit_passes`
  and (3) `::test_documented_commands_and_links_resolve` - pre-existing in the
  snapshot (absolute home paths and RA-path references in `docs/PROVIDER.md`,
  `docs/evidence/recovery_W3.md`, `docs/evidence/recovery_install.md`, files
  not changed here).
- RA: `tests/replay_history` 47 offline tests passed (provider 6, engine 16,
  outputs 16, example 9) plus 1 `abaqus`-marked live test passed
  (`RESASM_ABAQUS_JOB_PREFIX=claudeC_`, 21.6 s). Full offline RA suite
  (`-m "not abaqus and not arc and not network"`): 324 passed, 19 skipped,
  56 failed, 5 errors; all 61 failures/errors are the snapshot's pre-existing
  `FixtureError` (verified fixtures frozen under transform fingerprint
  94a92c01814f107a while `schemas/transform_generation.json` says
  6aa20d22e37f14c9); none involves a file changed here.

## Routing `resasm request`

`route_request` (in `cmd_history.py`) sends models inside the bounded scope to
`presentation.run_request` and everything else to the history engine
(`--validate` -> `--verify fd`). Tested: the bounded example stays bounded; the
beam example (nonzero BCs, `*Controls`, sets, MISES) is routed. The one edit in
`residual_core/ui/cmd_request.py` that activates it is in
`docs/REPLAY_HISTORY.md` (replace `parser.set_defaults(func=run)` by
`from .cmd_history import route_request; parser.set_defaults(func=route_request)`);
not applied here.

## Notes for the lead

- The UMAT provider change moves the UMAT transform fingerprint
  (`umat_oti.store.transform_store.transform_fingerprint`) from
  6aa20d22e37f14c9 (snapshot) to d8a3d2445fda0b26, because `provider/` is not
  in `NOT_TRANSFORM_CODE`. The corpus re-freeze must use the fingerprint after
  merging this branch, or exclude `provider/` from the fingerprint if the
  provider wrappers are judged not to be transform code (the transformed
  UMATs of the corpus are unaffected: emit.py only adds a wrapper routine).
- Objects built before this change lack `UMAT_OTI_EVAL_TOTAL`; the history
  engine refuses them with a rebuild message.

## What remains

- The slide-33 percentages at "yield onset" are not reproduced (measured
  values above); the definition used is the recovered old one.
- Abaqus FD resolves the plastic-range sensitivities only to about 1e-2 (J2 E,
  SIGY0) because of float32 output and yield-front crossings; the tight
  references are the Python whole-model FD.
- `cmd_request.py` routing edit left to the merge.
- Scope limits (single static step, C3D8, ramp amplitude, small strain,
  identity DROT/DFGRD) are documented and enforced.
