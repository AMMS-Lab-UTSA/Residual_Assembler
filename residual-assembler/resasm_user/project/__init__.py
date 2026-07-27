"""Project action engine (M5-A): resumable, idempotent sensitivity projects.

A thin, reusable layer the CLI, a desktop UI, or a future API all call. It owns a
persistent project (manifest + artifact directories), an explicit pipeline state
machine, and a set of logged, idempotent actions that wrap the proven M1-M4
backend. A failed action keeps the last good state so the pipeline can resume.

    from resasm_user.project import Project, run_pipeline, State
    p = Project.create("MyProj", name="demo", inp="model.inp", umat="umat.for")
    run_pipeline(p, State.INSPECTED)
"""

from __future__ import annotations

from .manifest import Project, ProjectError
from .state import State
from .engine import (Action, ActionResult, run_action, run_pipeline,
                     ordered_actions, get_action, register)

# importing these registers their actions in the engine registry
from . import inspect as _inspect          # noqa: F401,E402
from . import pipeline as _pipeline        # noqa: F401,E402
from .pipeline import set_parameter_selection, set_umat_oti_config   # noqa: E402

__all__ = [
    "Project", "ProjectError", "State", "Action", "ActionResult",
    "run_action", "run_pipeline", "ordered_actions", "get_action", "register",
    "set_parameter_selection", "set_umat_oti_config",
]
