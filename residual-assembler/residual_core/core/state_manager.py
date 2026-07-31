"""Per-element / per-integration-point material state across increments.

History-dependent materials (plasticity) require the state at the START of an
increment and produce a TRIAL state; only when the increment is accepted is the
trial committed. This manager keeps that discipline so a replay marches the full
loading history and never jumps to a final step (see docs/limitations.md).
"""

from __future__ import annotations

from typing import Dict, Optional

import numpy as np


class StateManager:
    def __init__(self, n_ip: int, n_state: int):
        self.n_ip = n_ip
        self.n_state = n_state
        self._committed: Dict[int, np.ndarray] = {}
        self._trial: Dict[int, np.ndarray] = {}

    def zeros(self) -> np.ndarray:
        return np.zeros((self.n_ip, self.n_state), dtype=float)

    def init_element(self, eid: int, state: Optional[np.ndarray] = None):
        self._committed[eid] = self.zeros() if state is None else np.asarray(state, float)

    def get(self, eid: int) -> np.ndarray:
        """State at increment start (committed)."""
        if eid not in self._committed:
            self._committed[eid] = self.zeros()
        return self._committed[eid]

    def set_trial(self, eid: int, state):
        if state is not None:
            self._trial[eid] = np.asarray(state, float)

    def commit(self):
        """Accept the increment: trial -> committed."""
        self._committed.update(self._trial)
        self._trial = {}

    def rollback(self):
        """Reject the increment: discard trial."""
        self._trial = {}
