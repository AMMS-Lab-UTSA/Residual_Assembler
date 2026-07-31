"""replay.py — history-aware single-element UMAT replay skeleton (edit me).

Sketches the increment-by-increment replay required for history-dependent
materials. Marches a per-increment strain history, seeding the design parameter
from the first increment so the state carries its derivative forward. Wire your
compiled UMAT into ``umat_stress_update`` (via f2py / ctypes) — keep it private.

This is a skeleton: it shows the loop shape and the provider surface, not a
working UMAT bridge.
"""

import numpy as np

from partner_kit.python.residual_provider import GlobalResidualProvider
from partner_kit.python.hypercomplex import Dual1, seed


def umat_stress_update(strain, dstrain, statev, props, dtime):
    """Replace with a call into your compiled UMAT. Must be scalar-generic if you
    want dual seeding (or return numeric STRESS/STATEV and use approach B/C from
    docs/history_dependent_models.md). Returns (stress, statev_new, ddsdde)."""
    raise NotImplementedError("bridge your UMAT here (f2py/ctypes)")


class HistoryReplayProvider(GlobalResidualProvider):
    name = "umat_replay"
    parameters = ("p1",)
    ndof = 1                      # single element / material point demo

    def __init__(self, strain_history, dtimes, props, statev0):
        self._hist = list(strain_history)   # per-increment total strain
        self._dt = list(dtimes)
        self._props = dict(props)
        self._statev0 = np.asarray(statev0, float)

    def parameter_values(self):
        return {p: float(self._props.get(p, 0.0)) for p in self.parameters}

    def eval_global_residual(self, U, state, parameters, time, dtime):
        # Replay the WHOLE history with the (possibly seeded) parameter.
        props = dict(self._props)
        props.update(parameters)            # parameters may carry a Dual1 seed
        statev = self._statev0.copy()
        prev = 0.0
        stress = 0.0
        for eps, dt in zip(self._hist, self._dt):
            deps = eps - prev
            stress, statev, _ddsdde = umat_stress_update(eps, deps, statev, props, dt)
            prev = eps
        # element residual at the final increment (1-DOF illustrative):
        # R = internal_force(stress) - external. Fill in your weak form.
        return [stress - props.get("external", 0.0)]


# NOTE: because the residual replays the full history internally, seeding the
# parameter (Dual1) at the top propagates its derivative through every increment —
# this is approach A from docs/history_dependent_models.md.
