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

__all__ = ["run_from_config", "check_config", "read_report", "job_output_dir",
           "ConfigError", "RunResult", "CheckReport"]


def job_output_dir(config_path: str) -> str:
    """The folder the job of ``config_path`` writes to, resolved exactly as
    :func:`run_from_config` resolves it: ``output: dir:`` (default
    ``resasm_output``) relative to the folder of the resasm.yml."""
    from .checks import _is_recipe
    if _is_recipe(config_path):
        from .recipe import load_recipe
        recipe = load_recipe(config_path)
        return recipe.path(recipe.output_dir)
    cfg = load_config(config_path)
    return cfg.path(cfg.output_dir)


def read_report(output_dir: str) -> Dict[str, Any]:
    """Load the machine-readable summary from a completed run's output dir."""
    pub = os.path.join(output_dir, "public")
    priv = os.path.join(output_dir, "private")
    out: Dict[str, Any] = {"output_dir": output_dir,
                           "private_dir": priv, "public_dir": pub}
    vs = os.path.join(pub, "validation_summary.json")
    meta = os.path.join(priv, "metadata.json")
    for key, path in (("validation_summary", vs), ("metadata", meta)):
        try:
            with open(path, "r", encoding="utf-8") as fh:
                payload = json.load(fh)
        except (OSError, ValueError) as exc:
            raise ConfigError("Cannot read report file %s: %s" % (path, exc)) from exc
        if not isinstance(payload, dict) or not payload:
            raise ConfigError("Report file %s must contain a non-empty JSON object" % path)
        out[key] = payload
    return out
