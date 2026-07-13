"""User-facing layer (resasm_user) — config validation, actionable errors, the
black-box path end-to-end, and the private/public output split.

The black-box path runs fully offline (no OTILib on the framework side). The
Python OTI path is exercised when OTILib is available; otherwise the actionable
'install OTILib' message is asserted.

Run: python tests/framework/test_user_layer.py
"""

import os
import shutil
import sys
import tempfile

import numpy as np
import pytest

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from resasm_user import run_from_config, check_config, read_report, ConfigError
from resasm_user.config import load_config


def _script_run(fn, *args):
    """Script-mode runner: a test passes unless it raises AssertionError."""
    try:
        fn(*args)
    except AssertionError as exc:
        print("  FAIL: %s" % exc)
        return False
    except pytest.skip.Exception as exc:
        print("  SKIP: %s" % exc)
    return True


def _report(checks):
    """Print every check, then assert that none failed (same condition as all())."""
    for label, v in checks:
        print("  %-38s -> %s" % (label, "PASS" if v else "FAIL"))
    failed = [label for label, v in checks if not v]
    assert not failed, "failed checks: " + "; ".join(failed)


_BLACKBOX_SOLVER = '''\
import argparse, json
import numpy as np

def R(u, p):
    return np.array([p["k"] * u[0]**3 - p["f"]], float)

ap = argparse.ArgumentParser()
ap.add_argument("--request"); ap.add_argument("--response")
a = ap.parse_args()
req = json.load(open(a.request))
u = np.array(req["u"], float); p = dict(req["parameters"]); names = list(p.keys())
dm = req["direction_map"]; ndof = u.size
T = np.array([[3.0*p["k"]*u[0]**2]])
dirs1 = dm.get("1", [])
R1 = np.zeros((ndof, len(dirs1)))
for c, exps in enumerate(dirs1):
    i = int(np.argmax(exps)); v0 = p[names[i]]; h = 1e-6*max(1.0, abs(v0))
    pp = dict(p); pp[names[i]] = v0+h; pm = dict(p); pm[names[i]] = v0-h
    R1[:, c] = (R(u, pp) - R(u, pm)) / (2*h)
np.savez(a.response, R_order_1=R1, tangent=T)
'''

_PY_RESIDUAL = '''\
def residual(u, params, state=None, time=None):
    return [params["k"] * u[0]**3 - params["f"]]

def tangent(u, params, state=None, time=None):
    return [[3.0 * params["k"] * u[0]**2]]
'''


def _write(d, name, text):
    with open(os.path.join(d, name), "w", encoding="utf-8") as fh:
        fh.write(text)


def _cfg(kind):
    common = ("parameters:\n  k: 2.0\n  f: 16.0\nsolution:\n  file: solution.npy\n")
    if kind == "blackbox":
        return ("problem:\n  name: bb\n  unknowns: 1\n"
                "residual:\n  type: executable\n"
                "  command: python my_solver.py --request {request} --response {response}\n"
                "tangent:\n  type: response\n"
                + common + "sensitivity:\n  order: 1\n  backend: otilib\n")
    return ("problem:\n  name: py\n  unknowns: 1\n"
            "residual:\n  type: python\n  module: user_residual.py\n  function: residual\n"
            "tangent:\n  type: python\n  function: tangent\n"
            + common + "sensitivity:\n  order: 1\n  backend: otilib\n")


def test_actionable_errors():
    checks = []
    d = tempfile.mkdtemp()
    try:
        # missing sensitivity.order
        _write(d, "resasm.yml",
               "problem:\n  name: x\n  unknowns: 1\nresidual:\n  type: python\n"
               "  module: user_residual.py\nparameters:\n  k: 2.0\n"
               "solution:\n  file: solution.npy\nsensitivity:\n  backend: otilib\n")
        try:
            load_config(os.path.join(d, "resasm.yml"))
            checks.append(("missing order raises", False))
        except ConfigError as exc:
            msg = str(exc)
            checks.append(("missing order actionable",
                           "sensitivity.order" in msg and "order: 1" in msg))
        # missing residual block
        _write(d, "r2.yml", "problem:\n  name: x\n  unknowns: 1\nparameters:\n"
               "  k: 2.0\nsolution:\n  file: s.npy\nsensitivity:\n  order: 1\n")
        try:
            load_config(os.path.join(d, "r2.yml"))
            checks.append(("missing residual raises", False))
        except ConfigError as exc:
            checks.append(("missing residual actionable",
                           "residual" in str(exc) and "type: python" in str(exc)))
        # bad backend order combo
        _write(d, "solution.npy", "")  # placeholder
    finally:
        shutil.rmtree(d, ignore_errors=True)
    _report(checks)


def test_blackbox_end_to_end():
    d = tempfile.mkdtemp()
    try:
        _write(d, "my_solver.py", _BLACKBOX_SOLVER)
        _write(d, "resasm.yml", _cfg("blackbox"))
        np.save(os.path.join(d, "solution.npy"), np.array([2.0]))
        res = run_from_config(os.path.join(d, "resasm.yml"))
        U1 = np.load(os.path.join(res.private_dir,
                                  "solution_sensitivities_order1.npz"))["U"].ravel()
        pub = os.listdir(res.public_dir)
        priv = os.listdir(res.private_dir)
        rep = read_report(res.output_dir)
        checks = [
            ("du/dk = -1/3", abs(U1[0] + 1.0 / 3.0) < 1e-5),
            ("du/df = 1/24", abs(U1[1] - 1.0 / 24.0) < 1e-5),
            ("public has no npz arrays", not any(f.endswith(".npz") for f in pub)),
            ("public has summary.md", "summary.md" in pub),
            ("public has rankings", "parameter_ranking.csv" in pub),
            ("private has U^(1)", "solution_sensitivities_order1.npz" in priv),
            ("private has tangent", "tangent.npz" in priv),
            ("report readable", rep["validation_summary"]["orders_solved"] == [1]),
        ]
    finally:
        shutil.rmtree(d, ignore_errors=True)
    _report(checks)


def test_python_path_check():
    """Python path: runs with OTILib; otherwise reports an actionable install
    message (never a raw KeyError/AttributeError)."""
    from residual_core.algebra.otilib_adapter import otilib_available
    d = tempfile.mkdtemp()
    try:
        _write(d, "user_residual.py", _PY_RESIDUAL)
        _write(d, "resasm.yml", _cfg("python"))
        np.save(os.path.join(d, "solution.npy"), np.array([2.0]))
        rep = check_config(os.path.join(d, "resasm.yml"))
        text = rep.render()
        if otilib_available():
            ok = rep.ok and "sensitivity solve completed" in text
            print("  python check (OTILib present) -> %s" % ("PASS" if ok else "FAIL"))
            detail = ("python path with OTILib: report.ok=%r, "
                      "'sensitivity solve completed' in report=%r"
                      % (rep.ok, "sensitivity solve completed" in text))
        else:
            actionable = (not rep.ok) and "mauriaristi/otilib" in text \
                and "pyoti" in text
            print("  python check OTILib-missing is actionable -> %s"
                  % ("PASS" if actionable else "FAIL"))
            ok = actionable
            detail = ("python path without OTILib: report.ok=%r (want False), "
                      "install hint present=%r, pyoti warning present=%r"
                      % (rep.ok, "mauriaristi/otilib" in text, "pyoti" in text))
    finally:
        shutil.rmtree(d, ignore_errors=True)
    assert ok, detail


def main():
    print("resasm_user user-facing layer")
    res = [_script_run(fn) for fn in (test_actionable_errors,
                                      test_blackbox_end_to_end,
                                      test_python_path_check)]
    ok = all(res)
    print("OVERALL: %s" % ("ALL PASS" if ok else "FAILURE"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
