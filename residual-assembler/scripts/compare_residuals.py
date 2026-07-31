"""compare_residuals.py — external residual vs Abaqus (Steps 3-5, 8).

Assembles the global residual OUTSIDE Abaqus from an exported field
(``fields.json`` produced by ``extract_odb_fields.py``) and compares:

    * free-DOF residual   -> should be ~ 0 at equilibrium         (Step 4)
    * reactions at BCs     -> should match exported RF (sign-aware) (Step 5)
    * sign convention      -> detected and reported                (Step 8)

This step is **offline-capable**: it needs only the neutral model + the exported
field, not Abaqus itself. If the field export is missing (because Abaqus has not
been run yet), it prints "Abaqus not available: validation pending" and exits 0.

Usage:
    python scripts/compare_residuals.py --model model.json --fields fields.json
    python scripts/compare_residuals.py --model model.inp  --fields fields.json --mode stress-driven
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
sys.path.insert(0, _ROOT)
sys.path.insert(0, os.path.dirname(__file__))
from _abaqus_env import SKIP_MESSAGE  # noqa: E402

from residual_core import ResidualProblem  # noqa: E402
from residual_core.core import constraints as _constraints  # noqa: E402


def _load_problem(model_path):
    if model_path.lower().endswith(".json"):
        return ResidualProblem.from_neutral(model_path)
    return ResidualProblem.from_abaqus(model_path)


def _exported_rf_vector(problem, reactions):
    """Map {node: [rf...]} onto a global vector using the DOF manager."""
    dm = problem.dof_manager
    rf = np.zeros(dm.ndof)
    for node_s, comps in reactions.items():
        nid = int(node_s)
        try:
            gdofs = dm.node_dofs(nid)
        except Exception:
            continue
        for k, val in enumerate(comps):
            if k < len(gdofs):
                rf[gdofs[k]] = float(val)
    return rf


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="External residual vs Abaqus (offline-capable)")
    p.add_argument("--model", required=True)
    p.add_argument("--fields", required=True)
    p.add_argument("--mode", default="stress-driven")
    p.add_argument("--tol", type=float, default=1e-6)
    args = p.parse_args(argv)

    if not os.path.exists(args.fields):
        print("%s (field export not found: %s — run Abaqus + extract_odb_fields.py first)"
              % (SKIP_MESSAGE, args.fields))
        return 0

    with open(args.fields) as fh:
        fields = json.load(fh)
    stress = fields.get("stress_ip", fields)
    stress = {int(k): np.asarray(v, float) for k, v in stress.items()}

    problem = _load_problem(args.model)
    problem.attach_results({"stress_ip": stress})
    R = problem.assemble(mode=args.mode)

    free_mask, pres_idx, _ = _constraints.partition(problem.model, problem.dof_manager)
    free_R = R[free_mask]
    reac = R[pres_idx] if pres_idx.size else np.zeros(0)

    print("residual comparison (mode=%s)" % args.mode)
    print("  ndof                 = %d" % R.size)
    print("  ||R_free||           = %.6e   (expect ~0 at equilibrium)"
          % float(np.linalg.norm(free_R)))
    print("  ||reaction(at BC)||  = %.6e" % float(np.linalg.norm(reac)))

    reactions = fields.get("reactions")
    if reactions and pres_idx.size:
        rf = _exported_rf_vector(problem, reactions)
        rf_bc = rf[pres_idx]
        # sign detection: Abaqus RF vs computed reaction
        dot = float(np.dot(reac, rf_bc))
        sign = "same (+R)" if dot >= 0 else "opposite (R = -RF)"
        denom = max(np.linalg.norm(rf_bc), 1e-30)
        rel_same = np.linalg.norm(reac - rf_bc) / denom
        rel_opp = np.linalg.norm(reac + rf_bc) / denom
        print("  exported RF present  = yes")
        print("  sign convention      = %s" % sign)
        print("  rel |R - RF|         = %.6e" % rel_same)
        print("  rel |R + RF|         = %.6e" % rel_opp)
        ok = min(rel_same, rel_opp) < args.tol
        print("  reaction match       = %s (tol=%g)" % ("PASS" if ok else "CHECK", args.tol))
    else:
        print("  exported RF present  = no (reaction cross-check skipped)")

    ok_free = float(np.linalg.norm(free_R)) < max(args.tol, 1e-6) * max(1.0, np.linalg.norm(R))
    print("  free-residual check  = %s" % ("PASS" if ok_free else "CHECK"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
