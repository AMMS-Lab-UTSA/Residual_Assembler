# History-Dependent Models

**Essential for plasticity, damage, viscoelasticity, and crystal plasticity.**

For a history-dependent material the stress at an increment depends on the
accumulated state, not just the current strain. Therefore:

> Overloading the residual only at the **final step** is **not** valid. The
> parameter perturbation must be propagated through the **entire load/time
> history**, because the state at step *s* depends on the parameter through all
> previous steps.

Formally, the transient residual is $r(u^s, u^{s-1}, a) = 0$ and the state
$q^s(a)$ evolves increment by increment. A single final-step evaluation ignores
$\partial q^{s-1}/\partial a$, giving a wrong $\partial u^s/\partial a$.

## Supported approaches (choose one)

### A. Replay the full history locally with hypercomplex seeding

Re-run the load/time history on the partner machine with the design parameter
seeded as a hypercomplex/dual number **from the first increment**. The state
carries its imaginary part forward automatically, so at the final step the
residual's imaginary part already includes the full history contribution.

- Provider shape: a `GlobalResidualProvider` (or element) whose `state` is a
  hypercomplex-capable object, driven in a loop over increments by your harness.
- Kit role: seed once, evaluate per increment, accumulate — the residual stays
  scalar-generic.

### B. Propagate state sensitivities increment by increment

Keep state real, but carry $\partial q^s/\partial a$ alongside $q^s$:

$$\frac{\partial q^{s}}{\partial a} = \frac{\partial q^{s}}{\partial a}\Big|_{\text{explicit}}
  + \frac{\partial q^{s}}{\partial q^{s-1}}\,\frac{\partial q^{s-1}}{\partial a}
  + \frac{\partial q^{s}}{\partial u^{s}}\,\frac{\partial u^{s}}{\partial a}.$$

At each step solve $T^s\,\partial u^s/\partial a = -\partial r^s/\partial a$ using
the *current* state sensitivity, then update $\partial q^s/\partial a$. This is
the classic direct-differentiation-of-history scheme.

### C. Expose a coupled residual including state/history variables

Treat the state variables as additional unknowns and expose one coupled residual
$R([u, q]) = 0$ with a block tangent. The kit then seeds parameters and solves the
coupled system — no special history bookkeeping, at the cost of a larger system.

## What the kit does today

- The order-1 dual pathway and the provider contract are **ready** for approaches
  A and C: your residual/state just need to be scalar-generic (A) or include the
  state DOFs (C).
- The kit does **not** yet ship an increment-marching driver or a state-sensitivity
  propagator (B). The `abaqus_umat_replay` template
  ([../templates/abaqus_umat_replay](../templates/abaqus_umat_replay)) sketches
  the single-element, increment-by-increment replay you would implement locally.

## Honest limitation

If you provide only a final-step residual for a history-dependent model, the
computed sensitivity is **not** correct. The kit will still run, but you must use
one of A/B/C above for a valid result. This is called out in the validation notes
and should be recorded in your `diagnostics.json`.

## Relationship to our core framework

The core framework's `StateManager` enforces committed/trial state across
increments and its `docs/umat_validation_plan.md` states the same rule: replay the
increment sequence in order. The partner kit inherits that discipline for
partner-side execution.
