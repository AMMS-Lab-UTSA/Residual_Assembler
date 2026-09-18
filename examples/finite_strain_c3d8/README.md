# Bounded Finite-Strain C3D8

This example uses the existing public `assemble` and `sensitivity` workflows,
not the presentation-request workflow (which remains pinned to J2).
The two distorted, shared-node C3D8 elements have 36 DOFs, 12 free DOFs, large
strain and a 0.7-radian superposed rotation. Manufactured nodal loads come from
independent reference-volume first-Piola integration at the prescribed target
solution. They remain fixed when material parameters are perturbed.

## Constitutive Law

The direct Python material `compressible_neo_hookean` uses explicit constants
`[mu, lambda] = [2.3, 4.1]` in consistent dimensionless benchmark units:

$$W(F)=\frac{\mu}{2}(\operatorname{tr}(FF^T)-3)-\mu\log J
       +\frac{\lambda}{2}(\log J)^2,\qquad J=\det F>0.$$

It returns Cauchy stress and the analytic Kirchhoff-Jaumann/J DDSDDE. The
formulation converts that tangent to the spatial elasticity modulus before
adding geometric stiffness. This is an exact linearization, not an approximate
`B^T DDSDDE B + Kgeo` shortcut.

## Reproduce Offline

Run commands sequentially from the recovery checkout. Genuine OTILib is
required; the unrelated PyPI `pyoti` is not a substitute.

```bash
cd /home/ammslab3/softwarex_work/imq-ra-recovery
export PYTHONPATH="$PWD:/home/ammslab3/softwarex_work/imq-umat-recovery/src:/home/ammslab3/otilib/build_py311"
export PYOTI_PATH=/home/ammslab3/otilib/build_py311
export OTILIB_ROOT=/home/ammslab3/otilib/build_py311
export RUN_OTILIB_TESTS=1
export UMAT_OTI_REPO=/home/ammslab3/softwarex_work/imq-umat-recovery
PY=/home/ammslab3/softwarex_work/.venv/bin/python
"$PY" examples/finite_strain_c3d8/benchmark.py
"$PY" -m residual_core.ui.cli --config examples/finite_strain_c3d8/verified/config.json assemble examples/finite_strain_c3d8/verified/model.json --mode material-replay --tangent
"$PY" -m residual_core.ui.cli --config examples/finite_strain_c3d8/verified/config.json sensitivity examples/finite_strain_c3d8/verified/model.json --params examples/finite_strain_c3d8/verified/params.json
```

The neutral [verified/model.json](verified/model.json) explicitly binds the finite
formulation and material. [verified/config.json](verified/config.json) carries
the real target solution under `options.solution`.
The config is a global CLI option and precedes the subcommand.
In the existing GUI, select this model/config, then `material-replay` in
Assemble and Sensitivity. The GUI invokes the same CLI backend; the finite
assembly bridge is tested. No new GUI or presentation interface is introduced.

Expected assembly: `ndof=36`, `||R||=1.744643e+00`. The nonzero full norm includes
reactions; the free residual norm is approximately `1.06e-15`.
[verified/report.json](verified/report.json) contains measured results:

| Check | Steps | Relative errors |
|---|---|---|
| Spatial force vs independent first-Piola integration | fixed state | 5.70e-16 |
| Complete 36-column nodal residual FD vs analytic K | 1e-4, 1e-5, 1e-6 | 4.44e-9, 4.55e-11, 1.09e-10 |
| OTILib du/d(mu,lambda) vs nonlinear central re-solves | relative 1e-3, 1e-4, 1e-5 | 9.48e-7, 9.48e-9, 1.01e-10 |

Acceptance thresholds are respectively `1e-12`, `2e-8`, and `2e-5`; all
three steps must meet the sweep thresholds. The nodal FD sweep shows the
roundoff floor, rather than claiming agreement from one lucky step.
The script exits nonzero if the report fails. Parameter sensitivity is generated
using actual `pyoti.sparse` scalars through the same material, formulation and
global assembler, then solved using the existing sensitivity package.

## Abaqus Check

One real Abaqus 2021.HF5 job, `imqr6_hyper`, completed in
`/home/ammslab3/softwarex_work/imq_abaqus/recovery_finite`.
The independent [neo_hookean_uhyper.for](neo_hookean_uhyper.for) provides energy
derivatives in modified invariants, not DDSDDE from the Python implementation.
Seven sequential steps prescribe a homogeneous finite deformation, then three
central perturbation pairs along one affine nodal direction. This avoids
confounding the comparison with Abaqus C3D8's mean-dilatation treatment of
nonhomogeneous volume changes. Offline distorted-element tests do perturb
every nodal DOF independently.

[verified/abaqus_report.json](verified/abaqus_report.json) records actual artifact
hashes, reaction error `3.97e-16`, and reaction-derivative errors `2.27e-7`,
`2.04e-8`, `2.27e-9` for steps `1e-3`, `3e-4`, `1e-4`. Tolerances are `2e-6`
for reactions and `2e-5` at every derivative step. AMATRX was not exported;
this is one independently observed tangent action, not a whole-matrix Abaqus
comparison. The exact commands are in the recovery evidence document.

## Scope And Diagnostics

- Full-integration C3D8; finite, real geometry with positive reference, previous
  and current Jacobians at all integration points. No geometric OTI directions.
- Isotropic, stateless total hyperelasticity. F0 and F1 are explicitly computed;
  F0 is validated but unused by this law. Stress is in global axes, so DROT is
  unused. Generic rotating plasticity/crystal-plasticity history is refused.
- Required material metadata: `deformation_gradient`, `cauchy`,
  `kirchhoff_jaumann_over_j`, `isotropic_total_hyperelastic`. Merely declaring
  `ddsdde` is insufficient. Incompatibility names the offending attribute.
- `mu > 0`, `lambda >= 0`; invalid or missing parameters are errors. The law is
  the logarithmic compressible neo-Hookean variant, not Abaqus's built-in
  isochoric/volumetric neo-Hookean potential.
- First-order material-parameter sensitivities only in the public finite path;
  higher-order equilibrium derivatives require live geometry differentiation
  and raise an explicit diagnostic. Parameters are `solid.mu`, `solid.lambda`.
- No follower loads, general history replay, mixed finite sensitivity model,
  arbitrary anisotropy, general load stepping, or global Newton convergence
  guarantee is established. Nonconvergence is reported, not returned as success.

See [../../docs/evidence/recovery_finite.md](../../docs/evidence/recovery_finite.md)
for the derivation, exact regression counts and unchanged stale-fixture failures.