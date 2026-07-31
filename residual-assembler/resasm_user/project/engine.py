"""The action engine (M5-A): idempotent, resumable, logged pipeline steps.

Every pipeline step is an :class:`Action` with

  * a minimum required state and the state it produces,
  * a deterministic set of output artifacts under the project directory,
  * a ``skip because already complete`` check,
  * structured logs, and
  * recoverable errors (a failure leaves the project at its last good state).

The UI, CLI and any future API all call :func:`run_action` / :func:`run_pipeline`
against the same actions -- no Abaqus or transformer logic lives in a UI callback.
"""

from __future__ import annotations

import time
import traceback
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional

from .manifest import Project, ProjectError
from .state import State


# --------------------------------------------------------------------------- #
class ActionLog:
    """Append-only structured log for one action run (also mirrored to a file)."""

    def __init__(self, path: str):
        self.path = path
        self.lines: List[str] = []
        self._fh = open(path, "w", encoding="utf-8")

    def __call__(self, message: str, level: str = "info") -> None:
        line = "%s [%s] %s" % (time.strftime("%H:%M:%S"), level.upper(), message)
        self.lines.append(line)
        self._fh.write(line + "\n")
        self._fh.flush()

    def close(self) -> None:
        try:
            self._fh.close()
        except Exception:              # noqa: BLE001
            pass


@dataclass
class ActionResult:
    action: str
    status: str                        # "complete" | "skipped" | "failed"
    message: str = ""
    outputs: Dict[str, Any] = field(default_factory=dict)
    artifacts: List[str] = field(default_factory=list)
    error: Optional[str] = None

    @property
    def ok(self) -> bool:
        return self.status in ("complete", "skipped")


# --------------------------------------------------------------------------- #
class Action:
    """Base class. Subclasses set ``name``/``requires``/``produces`` and
    implement :meth:`perform`."""

    name: str = "action"
    requires: State = State.NEW
    produces: State = State.NEW

    def artifacts(self, project: Project) -> List[str]:
        """Absolute paths this action produces (used by the completeness check)."""
        return []

    def is_complete(self, project: Project) -> bool:
        import os
        rec = project.action_record(self.name)
        if rec.get("status") != "complete":
            return False
        arts = self.artifacts(project)
        return all(os.path.exists(a) for a in arts)

    def perform(self, project: Project, log: ActionLog) -> Dict[str, Any]:
        raise NotImplementedError


# --------------------------------------------------------------------------- #
_REGISTRY: Dict[str, Action] = {}
_ORDER: List[str] = []


def register(action: Action) -> Action:
    _REGISTRY[action.name] = action
    _ORDER.append(action.name)
    return action


def get_action(name: str) -> Action:
    if name not in _REGISTRY:
        raise ProjectError("unknown action %r (known: %s)"
                           % (name, ", ".join(_ORDER)))
    return _REGISTRY[name]


def ordered_actions() -> List[Action]:
    return [_REGISTRY[n] for n in _ORDER]


# --------------------------------------------------------------------------- #
def run_action(project: Project, name: str, force: bool = False,
               **kwargs: Any) -> ActionResult:
    """Run one action with idempotency, structured logging and error recovery."""
    action = get_action(name)

    if project.state < action.requires:
        return ActionResult(name, "failed",
                            error="requires state >= %s but project is at %s"
                            % (action.requires, project.state))

    if not force and action.is_complete(project):
        return ActionResult(name, "skipped", message="already complete")

    log = ActionLog(project.logfile(name))
    log("action '%s' started (state=%s)" % (name, project.state))
    t0 = time.time()
    try:
        outputs = action.perform(project, log) or {}
        arts = action.artifacts(project)
        project.set_action(name, status="complete", finished=time.strftime(
            "%Y-%m-%dT%H:%M:%S"), seconds=round(time.time() - t0, 3),
            outputs=outputs, artifacts=[_rel(project, a) for a in arts], error=None)
        project.advance_to(action.produces)
        project.save()
        log("action '%s' complete -> state=%s" % (name, project.state))
        log.close()
        return ActionResult(name, "complete", outputs=outputs, artifacts=arts)
    except Exception as exc:                                   # noqa: BLE001
        tb = traceback.format_exc()
        log("action '%s' FAILED: %s" % (name, exc), level="error")
        log(tb, level="error")
        log.close()
        project.set_action(name, status="failed", finished=time.strftime(
            "%Y-%m-%dT%H:%M:%S"), error=str(exc))
        project.save()                        # state NOT advanced: resume point kept
        return ActionResult(name, "failed", error=str(exc),
                            message="see %s" % project.logfile(name))


def run_pipeline(project: Project, upto: State, force: bool = False,
                 on_result: Optional[Callable[[ActionResult], None]] = None,
                 ) -> List[ActionResult]:
    """Run actions in order until the project reaches ``upto`` (or an action
    fails). Completed actions are skipped unless ``force``."""
    results: List[ActionResult] = []
    for action in ordered_actions():
        if action.produces > upto:
            break
        res = run_action(project, action.name, force=force)
        results.append(res)
        if on_result:
            on_result(res)
        if res.status == "failed":
            break
    return results


def _rel(project: Project, path: str) -> str:
    import os
    try:
        return os.path.relpath(path, project.root)
    except ValueError:
        return path
