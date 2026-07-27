"""resasm_user — the minimal user-facing layer.

Four steps: implement one residual, write one small `resasm.yml`, run one command,
get a private/ + public/ sensitivity package.

Public API (this is all a user imports)::

    from resasm_user import run_from_config, check_config
    report = check_config("resasm.yml")   # readiness, actionable errors
    result = run_from_config("resasm.yml") # full sensitivity job

Internal framework classes (ResidualProblem, OtiContext, SensitivityPackage, ...)
are NOT part of this API.
"""

from __future__ import annotations

import json
import os
from typing import Any, Dict

from .config import ConfigError, load_config           # noqa: F401
from .checks import check_config, CheckReport          # noqa: F401
from .runner import run_from_config, RunResult         # noqa: F401

__all__ = ["run_from_config", "check_config", "read_report", "ConfigError",
           "RunResult", "CheckReport"]


def read_report(output_dir: str) -> Dict[str, Any]:
    """Load the machine-readable summary from a completed run's output dir."""
    pub = os.path.join(output_dir, "public")
    priv = os.path.join(output_dir, "private")
    out: Dict[str, Any] = {"output_dir": output_dir,
                           "private_dir": priv, "public_dir": pub}
    vs = os.path.join(pub, "validation_summary.json")
    meta = os.path.join(priv, "metadata.json")
    if os.path.exists(vs):
        with open(vs, "r", encoding="utf-8") as fh:
            out["validation_summary"] = json.load(fh)
    if os.path.exists(meta):
        with open(meta, "r", encoding="utf-8") as fh:
            out["metadata"] = json.load(fh)
    return out
