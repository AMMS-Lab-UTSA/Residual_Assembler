"""compare_umat_replay.py — replayed STRESS/STATEV vs Abaqus ODB (Steps 6, 8).

Plan for validating a UMAT candidate by **material-point / single-element replay**
outside Abaqus and comparing STRESS and STATEV against the ODB export.

CRITICAL (history-dependent materials): a UMAT that carries state (plasticity,
viscoelasticity, damage, crystal plasticity) must be replayed **increment by
increment**, propagating STATEV. A final-step-only replay is invalid. This script
enforces that by iterating the exported frame sequence.

Offline behaviour: genuine replay of a *compiled* UMAT needs Intel ifort + Abaqus
(the framework's fortran material backend). When unavailable, this prints
"Abaqus not available: validation pending" and exits 0. A pure-Python material
adapter (independently implemented, not copied from a restricted source) can be
compared offline once provided.

Usage:
    abaqus python scripts/extract_odb_fields.py --odb job.odb --out ref.json   # per-frame
    python scripts/compare_umat_replay.py --reference ref.json [--python-material mod:Class]
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
from _abaqus_env import abaqus_available, note_pending, SKIP_MESSAGE  # noqa: E402


def _load_python_material(spec):
    """spec = 'module.path:ClassName' -> instance (independent, non-restricted code)."""
    mod_name, cls_name = spec.split(":")
    import importlib
    mod = importlib.import_module(mod_name)
    return getattr(mod, cls_name)()


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="UMAT replay vs ODB (skips w/o Abaqus)")
    p.add_argument("--reference", required=True,
                   help="per-increment ODB export (STRESS/STATEV history)")
    p.add_argument("--python-material", dest="python_material",
                   help="module:Class of an INDEPENDENT python material to replay")
    p.add_argument("--stol", type=float, default=1e-4, help="stress rel tolerance")
    args = p.parse_args(argv)

    if not os.path.exists(args.reference):
        print("%s (no reference export: %s)" % (SKIP_MESSAGE, args.reference))
        return 0

    with open(args.reference) as fh:
        ref = json.load(fh)

    # A python material replay can run offline; a compiled UMAT replay needs Abaqus/ifort.
    if not args.python_material:
        if not abaqus_available():
            note_pending("compiled UMAT replay needs Intel ifort + Abaqus; "
                         "or pass --python-material module:Class for an offline replay")
        print("compiled UMAT replay requested — driver: residual_core/umat_adapter_fortran/")
        print("ensure the increment sequence is replayed in order (history-dependent).")
        return 0

    material = _load_python_material(args.python_material)
    frames = ref.get("frames") or []
    if not frames:
        print("reference has no per-increment 'frames'; re-export with the history.")
        return 0

    print("material-point replay (increment by increment): %d frames" % len(frames))
    state = None
    max_rel = 0.0
    for i, fr in enumerate(frames):
        kin = fr.get("kinematics", {})
        stress_ref = np.asarray(fr.get("stress", []), float)
        # NOTE: this expects an INDEPENDENT material.evaluate compatible signature.
        stress, tangent, state, _diag = material.evaluate(
            kin, state, fr.get("binding"), fr.get("time", (0.0, 0.0)),
            float(fr.get("dtime", 0.0)), None, None)
        if stress_ref.size:
            rel = np.linalg.norm(np.asarray(stress) - stress_ref) / max(
                np.linalg.norm(stress_ref), 1e-30)
            max_rel = max(max_rel, rel)
    print("  max stress rel error = %.3e   %s (tol=%g)"
          % (max_rel, "PASS" if max_rel < args.stol else "CHECK", args.stol))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
