"""Persistent project manifest + artifact directory layout (M5-A).

A project is a self-contained directory the tool owns::

    MyProject/
      project.json          <- the manifest (this module)
      inputs/               <- copies of the user's .inp and UMAT (never modified)
      generated/            <- transformed UMAT, patched .inp, layout, env, ...
      runs/                 <- Abaqus working dirs (one per run)
      results/              <- npz / manifest / report
      logs/                 <- one structured log per action

The manifest records the inputs, detected toolchain, inspection results, the
parameter configuration, the current pipeline state, and per-action status
(so actions can be skipped when already complete and the pipeline can resume
after a failure). It is plain JSON — machine-managed, not hand-edited.
"""

from __future__ import annotations

import json
import os
import time
from typing import Any, Dict, List, Optional

from .state import State

SCHEMA = "resasm_project_v1"
_SUBDIRS = ("inputs", "generated", "runs", "results", "logs")


class ProjectError(Exception):
    """A user-facing problem with a project or one of its actions."""


class Project:
    """Owns a project directory and its manifest."""

    def __init__(self, root: str, data: Dict[str, Any]):
        self.root = os.path.abspath(root)
        self.data = data

    # ---- lifecycle ------------------------------------------------------ #
    @classmethod
    def create(cls, root: str, name: str, inp: str, umat: Optional[str] = None) -> "Project":
        root = os.path.abspath(root)
        if os.path.exists(os.path.join(root, "project.json")):
            raise ProjectError("a project already exists at %s" % root)
        for d in _SUBDIRS:
            os.makedirs(os.path.join(root, d), exist_ok=True)
        data: Dict[str, Any] = {
            "schema": SCHEMA, "name": name, "created": _now(),
            "state": str(State.NEW), "inputs": {}, "toolchain": {},
            "model": {}, "umat": {}, "parameters": [], "actions": {},
            "assumptions": {"parameter_independent_geometry": True,
                            "parameter_independent_loads": True,
                            "parameter_independent_boundaries": True},
        }
        proj = cls(root, data)
        proj.data["inputs"]["inp"] = proj._ingest(inp, "inputs")
        if umat:
            proj.data["inputs"]["umat"] = proj._ingest(umat, "inputs")
        proj.save()
        return proj

    @classmethod
    def load(cls, root: str) -> "Project":
        path = os.path.join(os.path.abspath(root), "project.json")
        if not os.path.exists(path):
            raise ProjectError("no project at %s (expected project.json)" % root)
        with open(path, "r", encoding="utf-8") as fh:
            data = json.load(fh)
        if data.get("schema") != SCHEMA:
            raise ProjectError("unsupported project schema %r" % data.get("schema"))
        return cls(os.path.abspath(root), data)

    def save(self) -> None:
        with open(self.path("project.json"), "w", encoding="utf-8") as fh:
            json.dump(self.data, fh, indent=2)

    def _ingest(self, src: str, subdir: str) -> str:
        """Copy an input into the project (originals are never modified)."""
        if not os.path.exists(src):
            raise ProjectError("input file not found: %s" % src)
        import shutil
        dest = os.path.join(self.root, subdir, os.path.basename(src))
        shutil.copy(src, dest)
        return os.path.relpath(dest, self.root)

    # ---- paths ---------------------------------------------------------- #
    def path(self, *parts: str) -> str:
        return os.path.join(self.root, *parts)

    def artifact(self, subdir: str, name: str) -> str:
        return os.path.join(self.root, subdir, name)

    def input_path(self, key: str) -> Optional[str]:
        rel = self.data.get("inputs", {}).get(key)
        return self.path(rel) if rel else None

    # ---- state ---------------------------------------------------------- #
    @property
    def state(self) -> State:
        return State.parse(self.data.get("state", "NEW"))

    def advance_to(self, state: State) -> None:
        if state > self.state:
            self.data["state"] = str(state)

    # ---- per-action status ---------------------------------------------- #
    def action_record(self, name: str) -> Dict[str, Any]:
        return self.data.setdefault("actions", {}).get(name, {})

    def set_action(self, name: str, **fields: Any) -> None:
        rec = self.data.setdefault("actions", {}).setdefault(name, {})
        rec.update(fields)

    def logfile(self, name: str) -> str:
        return self.artifact("logs", "%s.log" % name)

    def summary(self) -> Dict[str, Any]:
        acts = self.data.get("actions", {})
        return {
            "name": self.data.get("name"), "root": self.root,
            "state": str(self.state),
            "inputs": self.data.get("inputs", {}),
            "parameters": [p.get("name") for p in self.data.get("parameters", [])],
            "actions": {k: v.get("status") for k, v in acts.items()},
        }


def _now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S")
