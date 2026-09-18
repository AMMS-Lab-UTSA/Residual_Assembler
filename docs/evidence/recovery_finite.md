# W6 Recovery: Bounded Finite Strain

Date: 2026-09-18. Worktree: `/home/ammslab3/softwarex_work/imq-ra-recovery`.
No commits, pushes, branch changes, subagents or concurrent solver jobs.
Existing worktree changes were preserved, including W1 scalar-preserving
`c3d8_kernel` edits. No fixture, fingerprint or contract-lock edits were made.

## Initial Evidence And Implementation

Read the shared `imq_BRIEF.md`, `docs/BRANCH_IMPLEMENTATION_AUDIT.md`, working-tree
status, local formulation/kernel/material contract and the old files using:

```bash
git show origin/cross-platform-hardening:residual-assembler/residual_core/formulations/solid_c3d8_finite_strain.py
git show origin/cross-platform-hardening:residual-assembler/residual_core/formulations/c3d8_kernel.py
```

The historical and current implementations used raw DDSDDE with geometric
stiffness without a stress correction. The current formulation explicitly
called this approximate but declared `verification_status = verified`, accepted
incompatible materials and stored live material results in float arrays.

Implemented the bounded direct `compressible_neo_hookean` material, exact
formulation conversion, explicit compatibility/geometry/history guards, registry
entry, scalar-preserving global scatter, and a first-order genuine OTILib
material-parameter provider. Generic finite UMAT/plasticity support is not claimed.
The status now says `bounded-hyperelastic` rather than unqualified `verified`.

Minimal public workflow repairs preserve explicit neutral material/formulation
bindings, honor backend config selection and `options.solution`, and return a
CLI diagnostic for unsupported measures. The bounded Newton path applies real
prescribed values, starts from the available solution, and rejects nonconvergence.
The sensitivity path uses the existing package and public CLI, with a tested
GUI-to-CLI assembly bridge. Presentation requests remain unchanged and J2-pinned.

## Exact Weak Linearization

Let $g_a=\nabla_x N_a$, $L=\delta F F^{-1}=d+w$, $d=d^T$, $w=-w^T$,
$\tau=J\sigma$. Abaqus's specified tangent is

$$D:d=\frac{\delta\tau-w\tau+\tau w}{J}.$$

Therefore

$$\delta\sigma=D:d+w\sigma-\sigma w-\sigma\operatorname{tr}(d).$$

For $r_a=\int_v\sigma g_a\,dv$, use
$\delta g_a=-L^Tg_a$ and $\delta dv=\operatorname{tr}(L)dv$:

$$\delta r_a=\int_v[(D:d-d\sigma-\sigma d)+L\sigma]g_a\,dv.$$

Thus $c:d=D:d-d\sigma-\sigma d$ and

$$K_{ab}=\int_v B_a^T c B_b\,dv+
\int_v(g_a^T\sigma g_b)I\,dv.$$

Engineering shear is used in Abaqus order `(11,22,33,12,13,23)`, so shear
basis tensors have two entries of `1/2`. The new kernel conversion does not
cast material scalars to float. The low-level `element_tangent` still accepts
the spatial modulus and does not guess an unspecified material measure.

For $W=\mu(I_1-3)/2-\mu\log J+\lambda(\log J)^2/2$ and $b=FF^T$:

$$\sigma=[\mu(b-I)+\lambda\log(J)I]/J,$$
$$D:d=[\mu(db+bd)+\lambda\operatorname{tr}(d)I]/J,$$
$$c:d=[2(\mu-\lambda\log J)d+\lambda\operatorname{tr}(d)I]/J.$$

The independent reference force uses
$P=\mu(F-F^{-T})+\lambda\log J F^{-T}$ integrated on reference volume.
FD checks recompute material stress and residual under perturbed nodal DOFs;
they do not compare the tangent against itself or hold stress artificially fixed.

F0/F1 use reference gradients and explicit previous/current displacements.
This total isotropic law validates F0 but depends only on F1; both stress and
tangent are in global axes. DROT is unused, which is valid only for this
stateless total response. Nonempty state or undeclared/rotating history models
are refused. No generic plastic history transport is implemented.

## Reproduction Environment

Every Python/test command used the shared healthy Python 3.11.7 venv explicitly,
from the RA recovery root. The import probe printed this checkout's
`residual_core/__init__.py`.

```bash
cd /home/ammslab3/softwarex_work/imq-ra-recovery
export PYTHONPATH="$PWD:/home/ammslab3/softwarex_work/imq-umat-recovery/src:/home/ammslab3/otilib/build_py311"
export PYOTI_PATH=/home/ammslab3/otilib/build_py311
export OTILIB_ROOT=/home/ammslab3/otilib/build_py311
export UMAT_OTI_REPO=/home/ammslab3/softwarex_work/imq-umat-recovery
export RUN_OTILIB_TESTS=1
PY=/home/ammslab3/softwarex_work/.venv/bin/python
"$PY" examples/finite_strain_c3d8/benchmark.py
```

Both OTI environment hints must target `build_py311`. An initial run with
`OTILIB_ROOT=/home/ammslab3/otilib` put an older source package ahead of the
build and produced 11 loading failures. Correcting only the environment gave
26/26 passes on the same checks; no substitute algebra or skip was introduced.

## Numerical Evidence

Machine-readable offline report:
[../../examples/finite_strain_c3d8/verified/report.json](../../examples/finite_strain_c3d8/verified/report.json).
Inputs, public CLI commands and acceptance thresholds:
[../../examples/finite_strain_c3d8/README.md](../../examples/finite_strain_c3d8/README.md).

| Independent check | Measured relative error |
|---|---:|
| Two distorted elements, spatial vs first-Piola force | 5.7034478275e-16 |
| Complete global K vs DOF FD, h=1e-4 | 4.4391866729e-9 |
| Complete global K vs DOF FD, h=1e-5 | 4.5520094802e-11 |
| Complete global K vs DOF FD, h=1e-6 | 1.0909865102e-10 |
| Genuine OTI du/dp vs central nonlinear re-solves, relative h=1e-3 | 9.4777210456e-7 |
| Same, relative h=1e-4 | 9.4779262310e-9 |
| Same, relative h=1e-5 | 1.0127558267e-10 |

Actual parameter perturbations in re-solves are `h * max(1, abs(parameter))`.
The free equilibrium residual norm is `1.0561716829e-15` and recovered target
solution error is `8.2321817675e-17`. Tests also cover distorted single elements,
large superposed rotations (residual covariance, invariant norm and tangent
covariance), rigid motion, the small-strain force/tangent limit, direct OTI
stress and stiffness parameter derivatives, invalid geometry/parameters,
incompatible measures, missing tangent and unsupported history/order requests.

## Real Abaqus Evidence

Exactly one sequential licensed analysis, `imqr6_hyper`, completed under Abaqus
2021.HF5 with ifort. No whole corpus, concurrent analysis or AMATRX export.
Commands, in order (the prepare step refuses to overwrite the existing job):

```bash
"$PY" examples/finite_strain_c3d8/abaqus_check.py prepare --work /home/ammslab3/softwarex_work/imq_abaqus/recovery_finite
cd /home/ammslab3/softwarex_work/imq_abaqus/recovery_finite
/usr/bin/abaqus job=imqr6_hyper input=imqr6_hyper.inp user=/home/ammslab3/softwarex_work/imq-ra-recovery/examples/finite_strain_c3d8/neo_hookean_uhyper.for cpus=1 interactive output_precision=full
/usr/bin/abaqus python /home/ammslab3/softwarex_work/imq-ra-recovery/examples/finite_strain_c3d8/extract_odb.py imqr6_hyper.odb fields.json
cd /home/ammslab3/softwarex_work/imq-ra-recovery
"$PY" examples/finite_strain_c3d8/abaqus_check.py compare --work /home/ammslab3/softwarex_work/imq_abaqus/recovery_finite
```

The ODB extractor necessarily uses Abaqus Python 2.7, not the venv.
UHYPER implements the same energy in modified invariants and supplies its
analytic energy derivatives, leaving Abaqus to form the constitutive response.
Its derivative ordering was checked against the published UHYPER interface:
https://abaqus-docs.mit.edu/2017/English/SIMACAESUBRefMap/simasub-c-uhyper.htm

The single homogeneous affine C3D8 has finite strain and a 0.8-radian rotation.
All 24 DOFs are prescribed, so RF is the internal force. Six further steps
perturb one affine direction by +/- h. All seven steps reached time 1.
Prescribed-U max error: `4.4408920985e-16`.
Reaction relative error: `3.9685038232e-16`.
Directional reaction derivative errors for `h=1e-3,3e-4,1e-4`:
`2.2686070112e-7, 2.0416638597e-8, 2.2721394357e-9`.

Retained numerical report with input/ODB/status/field/reference/UHYPER hashes:
[../../examples/finite_strain_c3d8/verified/abaqus_report.json](../../examples/finite_strain_c3d8/verified/abaqus_report.json).
Original solver artifacts remain in the designated workspace directory.
The homogeneous comparison avoids Abaqus C3D8 mean-dilatation differences;
arbitrary distorted/nonhomogeneous equivalence to Abaqus's element is NOT
established. The offline implementation uses ordinary full quadrature.

## Exact Test Results

Focused final regression command:

```bash
"$PY" -m pytest -q tests/framework/test_finite_strain_hyperelastic.py tests/framework/test_c3d8_sensitivity.py tests/framework/test_otilib_adapter.py tests/framework/test_otilib_fe_sensitivity.py tests/framework/test_otilib_spring_sensitivity.py tests/framework/test_gui_is_a_thin_cli_front_end.py tests/framework/test_neutral_io.py tests/framework/test_assembler.py
```

Result: **61 passed**, no skips/failures, 5.94 s. This includes **24 new finite
tests**. Genuine `pyoti.sparse` ran; GUI bridge tests did not skip.

Required full offline command:

```bash
"$PY" -m pytest -q -m 'not abaqus and not arc and not network' --tb=short --junitxml=docs/evidence/recovery_finite_offline.xml
```

Result: **1 collection error**, stopping before test execution. The exact
collector is `tests/framework/test_a_verified_deck_drives_the_global_assembly.py`.
Its `hex_fixtures()` rejects `compresibleneohookean--d78e61195f.json`: fixture
fingerprint `94a92c01814f107a` versus required `6aa20d22e37f14c9`.

To exercise the remaining tests without weakening the refusal:

```bash
"$PY" -m pytest -q -m 'not abaqus and not arc and not network' --continue-on-collection-errors --tb=line --junitxml=docs/evidence/recovery_finite_offline_final.xml
```

Final result: **325 passed, 56 failed, 11 skipped, 5 errors**, 43.555 s.
All 61 failures/errors contain the stale-transform-fingerprint refusal.
No finite-strain tests failed or skipped. Exact test names, outcomes and
messages are in [recovery_finite_failures.json](recovery_finite_failures.json);
complete tracebacks and skips are in
[recovery_finite_offline_final.xml](recovery_finite_offline_final.xml).

An earlier continued run without `UMAT_OTI_REPO` compared against the main
producing checkout: 323 passed, 58 failed, 11 skipped, 5 errors. Its two extra
failures were `test_the_shared_files_are_byte_identical_to_the_producing_repository`
and `test_both_repositories_record_the_same_lock`. Correcting the existing
environment override resolved both; no schema or lock was changed.

The offline suite is therefore **not green**. Fixture regeneration belongs to
the evidence workflow; stale fixtures were neither relabelled nor bypassed.
No completion-ledger claim is made for generic finite UMAT replay, plastic
rotation/history, higher-order geometric OTI, general Newton robustness,
follower loads, arbitrary Abaqus element equivalence or the full legacy suite.