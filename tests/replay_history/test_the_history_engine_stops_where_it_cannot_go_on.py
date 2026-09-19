"""The history engine stops, and says why, where going on would mean nothing.

Two corrections of the history engine had no test of their own:

* a UMAT that returns non-finite stresses at the state it is handed is refused
  by name (a7271cd); before, a replay whose last recorded increment overflowed
  the UMAT completed without an error;
* Newton stops when no step along its direction lowers the residual (15a5240),
  instead of taking a step that does not help and iterating on roundoff up to
  the iteration limit.

Both run the ORIGINAL J2 UMAT through the provider built from the UMAT-OTI
repository; nothing constitutive is replaced.
"""
import numpy as np
import pytest

from residual_core.replay.history import HistoryEngine, ReplayMismatch, run_history
from residual_core.replay.history_inputs import RecordedFields

from history_support import beam_model

pytestmark = pytest.mark.regression


def _engine(provider_factory, tmp_path):
    _, _, material = provider_factory("m3_j2")
    _, model = beam_model(tmp_path, n=(4, 2, 2), push=0.03, steps=2)
    return model, HistoryEngine(model, material)


def test_a_umat_that_returns_non_finite_stresses_is_refused_by_name(provider_factory, tmp_path, times):
    model, engine = _engine(provider_factory, tmp_path)
    solved = run_history(engine, times=times(2), rtol=1e-13)
    count = len(solved.increments) + 1
    U = np.zeros((count, engine.ndof))
    RF = np.zeros((count, engine.ndof))
    S = np.zeros((count, engine.ne, 8, 6))
    SDV = np.zeros((count, engine.ne, 8, engine.material.nstatev))
    for n, record in enumerate(solved.increments, 1):
        U[n], S[n], SDV[n] = record.u, record.stress, record.state
        RF[n, engine.constrained] = record.reaction[engine.constrained]
    # finite, so the provider takes it; the elastic stress E * strain overflows.
    # The last increment: no later increment hands the state back to the provider.
    U[2, engine.free[0]] = 1e306
    nodes = len(model.node_ids)
    recorded = RecordedFields(time=solved.times, U=U.reshape(count, nodes, 3), RF=RF.reshape(count, nodes, 3),
                              S=S, SDV=SDV, source="python", precision="float64")
    with pytest.raises(ReplayMismatch,
                       match="increment 2: the OTI UMAT returned non-finite stresses at the recorded state"):
        run_history(engine, fields=recorded)


def test_newton_stops_when_no_step_lowers_the_residual(provider_factory, tmp_path, times):
    _, engine = _engine(provider_factory, tmp_path)
    # 1e-30 of the residual scale is far below the roundoff floor of a double-precision solve
    with pytest.raises(ReplayMismatch, match=r"increment 1: Newton stagnated at max\|R_free\| = .*roundoff floor"):
        run_history(engine, times=times(2), rtol=1e-30)
