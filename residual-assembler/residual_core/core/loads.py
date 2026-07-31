"""External loads: F_external in R = F_int - F_external + F_constraints.

Currently: concentrated nodal loads (*Cload). Distributed loads (*Dsload) and
body forces are placeholders (parsed elsewhere; not yet assembled) -- a documented
limitation. An amplitude/load factor can scale the applied loads for stepping.
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from .constraints import resolve_target_nodes


def external_force(model, dof_manager, load_factor: float = 1.0,
                   amplitudes: Optional[dict] = None) -> np.ndarray:
    """Assemble the concentrated-load vector F_external (ndof,)."""
    F = np.zeros(dof_manager.ndof)
    for cl in model.cloads:
        scale = load_factor
        if amplitudes and getattr(cl, "amplitude", None) in amplitudes:
            scale *= amplitudes[cl.amplitude]
        for nid in resolve_target_nodes(model, cl.target):
            if nid in dof_manager.node_index and dof_manager.has_dof(nid, cl.dof):
                F[dof_manager.dof_index(nid, cl.dof)] += scale * float(cl.value)
    return F
