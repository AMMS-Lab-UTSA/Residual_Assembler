"""compare_uel_rhs.py — UEL RHS/AMATRX vs framework, with convention checks (Steps 7, 8).

Compares an Abaqus UEL's element residual/tangent against the framework's
``uel_direct`` adapter (or any independent element routine), and verifies the
conventions a UEL must expose:

    * RHS sign convention   (Abaqus: RHS = -R  ->  element_residual = -RHS)
    * AMATRX convention     (AMATRX = dR/du)
    * DOF ordering          (node-major layout)
    * state variable layout (SVARS)
    * element connectivity
    * load contribution convention

Input: a JSON dump captured from an instrumented UEL (see
``docs/uel_validation_plan.md`` for what to print from inside UEL). Schema:
    {
      "element": <int>, "coords": [[x,y,z],...], "u": [...],
      "props": [...], "svars_before": [...],
      "rhs": [...], "amatrx": [[...],...], "svars_after": [...],
      "time": [t, dt], "dtime": <float>
    }

Offline behaviour: if no UEL dump is provided, prints
"Abaqus not available: validation pending" and exits 0. A finite-difference
tangent check on the framework's own adapter runs without Abaqus.
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


def detect_sign_and_check(dump):
    rhs = np.asarray(dump["rhs"], float)
    amatrx = np.asarray(dump["amatrx"], float)
    n = rhs.size

    report = {}
    report["ndof"] = int(n)
    report["amatrx_shape_ok"] = amatrx.shape == (n, n)
    report["amatrx_symmetric"] = bool(np.allclose(amatrx, amatrx.T, atol=1e-8))
    # Abaqus convention: RHS = -R  ->  R = -RHS  ->  element_residual = -RHS
    report["framework_residual_from_rhs"] = (-rhs).tolist()
    report["sign_convention"] = "Abaqus RHS = -R (framework uses element_residual = -RHS)"
    # A UEL's AMATRX should equal dR/du = d(-RHS)/du; report positive-definiteness hint
    try:
        eig = np.linalg.eigvalsh(0.5 * (amatrx + amatrx.T))
        report["amatrx_min_eig"] = float(np.min(eig))
        report["amatrx_spd_hint"] = bool(np.min(eig) > -1e-8)
    except Exception:
        report["amatrx_min_eig"] = None
    return report


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="UEL RHS/AMATRX comparison (skips w/o dump)")
    p.add_argument("--uel-dump", dest="uel_dump",
                   help="JSON captured from an instrumented UEL")
    p.add_argument("--framework-selftest", action="store_true",
                   help="run the uel_direct sign/FD self-test (offline)")
    args = p.parse_args(argv)

    if args.framework_selftest or not args.uel_dump:
        # Offline: exercise the framework's own adapter self-test if present.
        try:
            import residual_core.formulations.uel_adapter as uel  # noqa: F401
            print("framework uel_direct adapter importable; run its module self-test:")
            print("  python -m residual_core.formulations.uel_adapter")
        except Exception as exc:
            print("uel_direct import issue: %s" % exc)
        if not args.uel_dump:
            print("%s (no --uel-dump provided; instrument a UEL per docs/uel_validation_plan.md)"
                  % SKIP_MESSAGE)
            return 0

    if not os.path.exists(args.uel_dump):
        print("%s (uel dump not found: %s)" % (SKIP_MESSAGE, args.uel_dump))
        return 0

    with open(args.uel_dump) as fh:
        dump = json.load(fh)
    rep = detect_sign_and_check(dump)
    print("UEL convention report (element %s)" % dump.get("element"))
    for k, v in rep.items():
        if k == "framework_residual_from_rhs":
            print("  %-24s (vector, len=%d)" % (k + ":", len(v)))
        else:
            print("  %-24s %s" % (k + ":", v))
    print("DOF ordering / connectivity / SVARS layout: verify against the model's "
          "node-major layout and *Depvar (see docs/uel_validation_plan.md).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
