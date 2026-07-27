"""abaqus_job.py -- single-file residual verification against an Abaqus run.

The user writes **one** small JSON file that names the ingredients, and runs
**one** command::

    resasm verify-job job.json

A job file (JSON; YAML also accepted if PyYAML is installed)::

    {
      "name": "elastic_c3d8",
      "model": "elastic_c3d8.inp",      # Abaqus .inp OR neutral .json (mesh + BCs)
      "odb":   "elastic_c3d8.odb",      # OR  "fields": "fields.json" (already exported)
      "mode":  "stress-driven",         # optional (default stress-driven)
      "compare_reactions": true,          # optional (default true)
      "tol": 1e-6,                        # optional relative tolerance
      "output": "out"                     # optional output dir (default resasm_output)
    }

What it does, in order:

1. If ``odb`` is given and ``fields`` is not, export the ODB fields by calling
   ``abaqus python scripts/extract_odb_fields.py`` (needs Abaqus on this machine).
   If ``fields`` is given, use it directly -- no Abaqus needed (offline path).
2. Load the model (``.inp`` -> ``from_abaqus``; ``.json`` -> ``from_neutral``).
3. Attach the exported integration-point stress and assemble the residual.
4. Split free vs. prescribed DOFs; on free DOFs the residual must be ~0, on
   prescribed DOFs it must equal the Abaqus reaction ``RF`` up to one global sign.
5. Write ``<output>/verification_report.json`` and print a PASS/FAIL summary.

Relative paths inside the job file are resolved against the job file's folder.
"""

from __future__ import annotations

import json
import os
import shlex
import subprocess
import sys
from dataclasses import dataclass, field
from typing import Any, Dict, Optional

import numpy as np


class JobError(Exception):
    """A user-facing problem with the job file or its inputs."""


# --------------------------------------------------------------------------- #
# job file loading
# --------------------------------------------------------------------------- #
def load_job(path: str) -> Dict[str, Any]:
    """Read and lightly validate a job file. Returns the raw dict plus a
    resolved ``base_dir`` (the folder the job file lives in)."""
    if not os.path.exists(path):
        raise JobError("job file not found: %s" % path)
    with open(path, "r", encoding="utf-8") as fh:
        text = fh.read()
    data: Optional[Dict[str, Any]] = None
    try:
        data = json.loads(text)
    except json.JSONDecodeError:
        try:
            import yaml  # optional
            data = yaml.safe_load(text)
        except Exception:
            raise JobError(
                "cannot parse %s as JSON (and PyYAML is not installed for YAML). "
                "Write a JSON object with at least 'model' and 'odb' or 'fields'." % path)
    if not isinstance(data, dict):
        raise JobError("job file must be a JSON object, got %s" % type(data).__name__)
    if "model" not in data:
        raise JobError("job file '%s' is missing required key 'model' "
                       "(path to the .inp or neutral .json)." % path)
    if "odb" not in data and "fields" not in data:
        raise JobError("job file '%s' needs either 'odb' (an Abaqus .odb, which is "
                       "auto-exported) or 'fields' (a pre-exported fields.json)." % path)
    data["base_dir"] = os.path.dirname(os.path.abspath(path)) or "."
    return data


def _resolve(base_dir: str, p: Optional[str]) -> Optional[str]:
    if not p:
        return None
    return p if os.path.isabs(p) else os.path.normpath(os.path.join(base_dir, p))


# --------------------------------------------------------------------------- #
# ODB export (only when 'odb' is given and 'fields' is not)
# --------------------------------------------------------------------------- #
def _extractor_path() -> str:
    here = os.path.dirname(os.path.abspath(__file__))     # resasm_user/
    root = os.path.dirname(here)                           # repo root
    return os.path.join(root, "scripts", "extract_odb_fields.py")


def _abaqus_cmd() -> str:
    return os.environ.get("ABAQUS_CMD") or "abaqus"


def export_fields_from_odb(odb_path: str, out_path: str) -> str:
    """Run ``abaqus python extract_odb_fields.py`` to produce ``out_path``.
    Raises JobError with an actionable message if Abaqus is missing or the
    export produced no stress."""
    if not os.path.exists(odb_path):
        raise JobError("odb not found: %s" % odb_path)
    extractor = _extractor_path()
    if not os.path.exists(extractor):
        raise JobError("cannot locate the ODB extractor at %s "
                       "(run from a checkout or `pip install -e .`)." % extractor)
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    q = (lambda s: '"%s"' % s) if os.name == "nt" else shlex.quote
    cmd = '%s python %s --odb %s --out %s' % (
        _abaqus_cmd(), q(extractor), q(odb_path), q(out_path))
    try:
        # shell=True so a Windows abaqus.bat resolves via PATH like it does in a shell
        rc = subprocess.call(cmd, shell=True)
    except OSError as exc:
        raise JobError("could not launch Abaqus ('%s'): %s. Set ABAQUS_CMD or "
                       "provide 'fields' instead of 'odb'." % (_abaqus_cmd(), exc))
    if rc != 0:
        raise JobError("ODB export failed (exit %d). Check that '%s' runs and the "
                       "ODB is readable." % (rc, _abaqus_cmd()))
    if not os.path.exists(out_path):
        raise JobError("ODB export produced no file. Abaqus may be unavailable; "
                       "provide 'fields' (a pre-exported fields.json) instead.")
    with open(out_path, "r", encoding="utf-8") as fh:
        d = json.load(fh)
    if not d.get("stress_ip"):
        raise JobError("exported %s has no integration-point stress ('stress_ip'). "
                       "Request S in the Abaqus field output." % out_path)
    return out_path


# --------------------------------------------------------------------------- #
# assembly + reaction comparison
# --------------------------------------------------------------------------- #
@dataclass
class JobResult:
    name: str
    mode: str
    ndof: int
    free_residual_norm: float
    reaction_norm: float
    compared_reactions: bool = False
    sign: Optional[str] = None
    rel_reaction_error: Optional[float] = None
    reaction_pass: Optional[bool] = None
    free_pass: bool = False
    tol: float = 1e-6
    report_path: Optional[str] = None
    extras: Dict[str, Any] = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        rp = True if self.reaction_pass is None else self.reaction_pass
        return bool(self.free_pass and rp)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name, "mode": self.mode, "ndof": self.ndof,
            "free_residual_norm": self.free_residual_norm,
            "reaction_norm": self.reaction_norm,
            "compared_reactions": self.compared_reactions,
            "sign_convention": self.sign,
            "rel_reaction_error": self.rel_reaction_error,
            "reaction_match": self.reaction_pass,
            "free_residual_ok": self.free_pass,
            "tolerance": self.tol, "pass": self.ok, "notes": self.extras,
        }


def _load_problem(model_path: str):
    from residual_core import ResidualProblem
    if model_path.lower().endswith(".json"):
        return ResidualProblem.from_neutral(model_path)
    return ResidualProblem.from_abaqus(model_path)


def _rf_vector(problem, reactions, ndof):
    dm = problem.dof_manager
    rf = np.zeros(ndof)
    for node_s, comps in reactions.items():
        try:
            gdofs = dm.node_dofs(int(node_s))
        except Exception:
            continue
        for k, val in enumerate(comps):
            if k < len(gdofs):
                rf[gdofs[k]] = float(val)
    return rf


def run_job(path: str) -> JobResult:
    """Execute a verification job file and return a :class:`JobResult`."""
    job = load_job(path)
    base = job["base_dir"]
    name = job.get("name") or os.path.splitext(os.path.basename(path))[0]
    mode = job.get("mode", "stress-driven")
    tol = float(job.get("tol", 1e-6))
    compare_reactions = bool(job.get("compare_reactions", True))
    out_dir = _resolve(base, job.get("output", "resasm_output"))
    os.makedirs(out_dir, exist_ok=True)

    model_path = _resolve(base, job["model"])
    if not os.path.exists(model_path):
        raise JobError("model not found: %s" % model_path)

    # 1) obtain the exported fields (extract from ODB, or use the given file)
    fields_path = _resolve(base, job.get("fields"))
    if fields_path is None:
        odb_path = _resolve(base, job["odb"])
        fields_path = export_fields_from_odb(odb_path, os.path.join(out_dir, "fields.json"))
    elif not os.path.exists(fields_path):
        raise JobError("fields file not found: %s" % fields_path)

    with open(fields_path, "r", encoding="utf-8") as fh:
        fields = json.load(fh)
    stress_raw = fields.get("stress_ip", fields)
    stress = {int(k): np.asarray(v, float) for k, v in stress_raw.items()}

    # 2-3) load model, attach stress, assemble
    problem = _load_problem(model_path)
    problem.attach_results({"stress_ip": stress})
    R = problem.assemble(mode=mode)

    # 4) partition and compare
    from residual_core.core import constraints as _constraints
    free_mask, pres_idx, _ = _constraints.partition(problem.model, problem.dof_manager)
    free_R = R[free_mask]
    reac = R[pres_idx] if pres_idx.size else np.zeros(0)

    res = JobResult(name=name, mode=mode, ndof=int(R.size),
                    free_residual_norm=float(np.linalg.norm(free_R)),
                    reaction_norm=float(np.linalg.norm(reac)), tol=tol)
    res.extras["fields"] = fields_path
    res.extras["model"] = model_path

    Rscale = max(1.0, float(np.linalg.norm(R)))
    res.free_pass = res.free_residual_norm <= max(tol, 1e-6) * Rscale

    reactions = fields.get("reactions")
    if compare_reactions and reactions and pres_idx.size:
        rf = _rf_vector(problem, reactions, R.size)[pres_idx]
        denom = max(float(np.linalg.norm(rf)), 1e-30)
        rel_same = float(np.linalg.norm(reac - rf) / denom)
        rel_opp = float(np.linalg.norm(reac + rf) / denom)
        res.compared_reactions = True
        res.sign = "R = +RF" if rel_same <= rel_opp else "R = -RF"
        res.rel_reaction_error = min(rel_same, rel_opp)
        res.reaction_pass = res.rel_reaction_error < tol

    # 5) report
    report_path = os.path.join(out_dir, "verification_report.json")
    with open(report_path, "w", encoding="utf-8") as fh:
        json.dump(res.to_dict(), fh, indent=2)
    res.report_path = report_path
    return res


def render(res: JobResult) -> str:
    lines = ["residual verification: %s  (mode=%s)" % (res.name, res.mode),
             "  ndof                 = %d" % res.ndof,
             "  ||R_free||           = %.6e   (expect ~0 at equilibrium)"
             % res.free_residual_norm,
             "  free-residual check  = %s" % ("PASS" if res.free_pass else "CHECK")]
    if res.compared_reactions:
        lines += ["  ||reaction (at BC)|| = %.6e" % res.reaction_norm,
                  "  sign convention      = %s" % res.sign,
                  "  rel |R - RF|         = %.6e" % res.rel_reaction_error,
                  "  reaction match       = %s (tol=%g)"
                  % ("PASS" if res.reaction_pass else "CHECK", res.tol)]
    else:
        lines += ["  reaction comparison  = skipped (no RF, no prescribed DOFs, "
                  "or compare_reactions=false)"]
    lines += ["  overall              = %s" % ("PASS" if res.ok else "CHECK"),
              "  report               = %s" % res.report_path]
    return os.linesep.join(lines)


def main(argv=None) -> int:
    import argparse
    p = argparse.ArgumentParser(
        prog="resasm verify-job",
        description="Assemble the residual and verify it against an Abaqus run, "
                    "driven by a single job file.")
    p.add_argument("job", help="path to the job file (.json/.yml)")
    args = p.parse_args(argv)
    try:
        res = run_job(args.job)
    except JobError as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        return 2
    print(render(res))
    return 0 if res.ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
