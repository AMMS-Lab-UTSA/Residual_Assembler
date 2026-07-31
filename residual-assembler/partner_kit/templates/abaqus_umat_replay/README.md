# Template — Abaqus UMAT single-element replay (history-aware)

For a UMAT-driven material the residual is history-dependent, so parameter
sensitivities require replaying the **whole increment sequence** with the
parameter seeded — not just the final step. This template sketches the local,
single-element replay you implement on your machine.

## Idea

1. Drive one element (or one material point) through the exact strain/time history
   your job saw (export it once from your solver).
2. Replay increment by increment. Seed the design parameter with a dual/OTI value
   from the **first** increment so the state carries its derivative forward.
3. At the final increment, assemble the element residual; its imaginary part is
   `dR/dp` for that element — scatter into `R^(1)`.
4. Provide the element tangent (from your UMAT's DDSDDE) so the kit can solve
   `T U^(1) = -R^(1)`.

## Why final-step-only is wrong

The stress at step *s* depends on `STATEV` accumulated over all prior steps. A
single differentiation at the last step omits `d(STATEV)/dp` and gives the wrong
sensitivity. See
[../../docs/history_dependent_models.md](../../docs/history_dependent_models.md).

## Skeleton

`replay.py` shows the loop shape as a Python `GlobalResidualProvider` whose state
is marched over increments. Wire your compiled UMAT (via f2py/ctypes) into
`umat_stress_update(...)`. Keep `strain`, `STATEV`, and the material model private
— only the element residual coefficients are ever surfaced.

## Validation

Cross-check the replayed STRESS/STATEV against your ODB export per increment
before trusting the sensitivity (see
[../../docs/validation_checklist.md](../../docs/validation_checklist.md)).
