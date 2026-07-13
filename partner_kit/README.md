# Partner Kit — local residual-sensitivity SDK

A self-contained kit a collaborator runs **on their own machine** to generate
residual-sensitivity right-hand-sides and solve/export

$$\mathbf{T}\,\mathbf{U}^{(p)} = -\mathbf{R}^{(p)}$$

**without sharing** source code, mesh, input deck, material model, or the full
residual/tangent. The partner keeps everything private and sends back only a
small, non-sensitive report if they choose to.

## The honest claim

> We do **not** need to see your residual. But the residual must be **evaluable
> locally** by you, through one of the provider interfaces below. We cannot
> compute sensitivities without a locally evaluable residual.

## How it works

1. You implement a **residual provider** locally (element, global, or a black-box
   executable) — see [docs/residual_provider_contract.md](docs/residual_provider_contract.md).
2. The kit obtains a converged solution `u` (you supply it or the kit Newton-solves
   using your residual + tangent).
3. It **seeds a design parameter** with a hypercomplex/dual number and evaluates
   your residual, reading $\partial R/\partial a_i$ from the imaginary part → this
   is $\mathbf{R}^{(1)}$. (Black-box providers do the seeding inside the
   executable and just return the coefficients.)
4. It acquires the tangent $\mathbf{T}$ (assembled, `get_tangent`, matrix-free, or
   an external file), solves $\mathbf{T}\,\mathbf{U}^{(1)} = -\mathbf{R}^{(1)}$,
   and validates locally.
5. It writes a **private** `sensitivity_package/` (full arrays, stays with you)
   and a **public** `public_report/` (norms, rankings, validation summary — no
   proprietary data).

## Quick start

```bash
python -m partner_kit.python.partner_cli --config examples/nonlinear_spring_private/config.json
```

More: [docs/partner_quickstart.md](docs/partner_quickstart.md).

## Layout

```
partner_kit/
  README.md
  docs/        quickstart, privacy model, provider contract, output contract,
               validation checklist, history-dependent models
  include/     C++ contract headers (templated on Scalar); otilib_scalar.hpp seam
  python/      the runnable kit: hypercomplex, provider contract, blackbox runner,
               validators, sensitivity core, partner_cli
  templates/   starter skeletons: C++/Fortran OTILib residual, C++ element/global
               residual, Fortran UEL wrapper, Abaqus UMAT replay, black-box exe
  examples/    nonlinear_spring_private, nonlinear_bar_private, blackbox_residual_demo
```

## Privacy in one line

Private (never leaves you): source, mesh, material model, state, full residual,
full tangent. Optionally shareable: validation report, timings, derivative norms,
parameter ranking, selected sensitivities. See
[docs/privacy_model.md](docs/privacy_model.md).

## Scope / honesty

- **OTILib is the production backend** for arbitrary-order, multi-parameter
  $\mathbf{R}^{(p)}$ — install it locally from
  <https://github.com/mauriaristi/otilib.git> (branch `master`, **GPLv3**,
  Windows = WSL only) and compile your residual once as `template<class Scalar>`;
  see [docs/otilib_partner_quickstart.md](docs/otilib_partner_quickstart.md) and
  `templates/cpp_otilib_residual` / `templates/fortran_otilib_residual`.
- The tiny dual number (`python/hypercomplex.py`) is a **legacy first-order smoke
  test** so the kit runs standalone; it is not the production path and must not be
  grown toward OTI.
- **History-dependent models** (plasticity, damage, viscoelasticity, crystal
  plasticity) require replaying the load history — final-step overloading is not
  valid. See [docs/history_dependent_models.md](docs/history_dependent_models.md).
- OTILib is **not bundled** (external GPLv3 dependency;
  <https://github.com/mauriaristi/otilib.git>; do not `pip install pyoti` — that
  name is squatted). This kit defines the contracts and the templated seam.
