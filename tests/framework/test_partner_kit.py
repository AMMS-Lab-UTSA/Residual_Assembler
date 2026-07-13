"""Partner kit — end-to-end local sensitivity generation (privacy-preserving).

Runs the three shipped examples through the kit and checks the sensitivities,
validation, and the private/public output split. Proves a partner can generate
R^(1) and solve T U^(1) = -R^(1) locally without sharing code — element-level,
global-level, and black-box executable paths.

Run: python tests/framework/test_partner_kit.py
"""

import json
import os
import sys
import tempfile

import numpy as np
import pytest

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from partner_kit.python import partner_cli
from partner_kit.python.hypercomplex import Dual1, seed, imag_part
from partner_kit.python import validators as _val

_EX = os.path.join(_ROOT, "partner_kit", "examples")


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
        print("  %-40s -> %s" % (label, "PASS" if v else "FAIL"))
    failed = [label for label, v in checks if not v]
    assert not failed, "failed checks: " + "; ".join(failed)


def _run(example_dir, overrides=None):
    cfg = partner_cli._load_config(os.path.join(example_dir, "config.json"))
    out = tempfile.mkdtemp()
    cfg["output_dir"] = out
    if overrides:
        cfg.update(overrides)
    ctx = partner_cli.run_kit(cfg, base_dir=example_dir)
    return ctx, out


def test_dual_algebra():
    x = seed(2.0)
    r = x ** 3 - 16.0            # d/dx = 3x^2 = 12 at x=2
    ok = abs(imag_part(r) - 12.0) < 1e-12 and abs(r.real - (8.0 - 16.0)) < 1e-12
    print("  %-40s -> %s" % ("dual1 imag = derivative", "PASS" if ok else "FAIL"))
    assert ok, ("dual1: imag=%.12g (need 12), real=%.12g (need -8)"
                % (imag_part(r), r.real))


def test_spring_global():
    ctx, out = _run(os.path.join(_EX, "nonlinear_spring_private"))
    dudk = float(ctx["U1"][0, 0])
    checks = [
        ("du/dk = -1/3", abs(dudk + 1.0 / 3.0) < 1e-9),
        ("hypercomplex ready", ctx["rhs_diag"].get("hypercomplex_ready") is True),
        ("validation overall pass", _val.overall(ctx["checks"]) is True),
        ("tangent available", ctx["T"] is not None),
        ("private package written",
         os.path.exists(os.path.join(out, "sensitivity_package", "rhs_order_1.npz"))),
        ("public report written",
         os.path.exists(os.path.join(out, "public_report", "validation_summary.md"))),
    ]
    # privacy: public report has NO npz arrays by default
    pub = os.listdir(os.path.join(out, "public_report"))
    checks.append(("public report has no arrays", not any(f.endswith(".npz") for f in pub)))
    with open(os.path.join(out, "sensitivity_package", "metadata.json"),
              encoding="utf-8") as fh:
        meta = json.load(fh)
    checks.append(("target system recorded", meta["target_system"] == "T U^(p) = -R^(p)"))
    _report(checks)


def test_bar_element():
    ctx, _out = _run(os.path.join(_EX, "nonlinear_bar_private"))
    U1 = ctx["U1"]
    checks = [
        ("R^(1) shape (3,2)", np.asarray(ctx["R1"]).shape == (3, 2)),
        ("tangent backend-assembled", ctx["T_source"] == "backend-assembled"),
        ("dU/dk1 ~ -1/6", abs(float(U1[1, 0]) + 1.0 / 6.0) < 1e-6),
        ("dU/dk2 ~ +1/6", abs(float(U1[1, 1]) - 1.0 / 6.0) < 1e-6),
        ("prescribed rows zero",
         np.allclose(U1[~ctx["free"]], 0.0)),
        ("validation overall pass", _val.overall(ctx["checks"]) is True),
    ]
    _report(checks)


def test_blackbox_demo():
    ctx, _out = _run(os.path.join(_EX, "blackbox_residual_demo"))
    dudk = float(ctx["U1"][0, 0])
    checks = [
        ("du/dk = -1/3 (blackbox)", abs(dudk + 1.0 / 3.0) < 1e-9),
        ("tangent from response", ctx["T_source"] == "blackbox-response"),
        ("FD checks pending (no local residual)",
         any(c["name"] == "rhs-vs-FD" and c["passed"] is None for c in ctx["checks"])),
        ("sensitivity-solve passes",
         any(c["name"] == "sensitivity-solve" and c["passed"] for c in ctx["checks"])),
    ]
    _report(checks)


def test_blackbox_reference_template():
    tmpl = os.path.join(_ROOT, "partner_kit", "templates", "blackbox_executable")
    cfg = {
        "blackbox": {"command": ["python", "reference_solver.py"], "ndof": 1,
                     "parameters": ["k"],
                     "parameter_values": {"k": 2.0, "f": 16.0}, "free_mask": [True]},
        "parameters": ["k"], "order": 1, "solution": [2.0],
        "output_dir": tempfile.mkdtemp(),
    }
    ctx = partner_cli.run_kit(cfg, base_dir=tmpl)
    dudk = float(ctx["U1"][0, 0])
    ok = abs(dudk + 1.0 / 3.0) < 1e-5      # FD-based executable -> looser tol
    print("  %-40s -> %s" % ("reference_solver template du/dk~-1/3",
                             "PASS" if ok else "FAIL"))
    assert ok, ("reference_solver template: du/dk=%.9g (need -1/3, tol 1e-5)" % dudk)


def main():
    print("Partner kit end-to-end (local sensitivity, no code sharing)")
    results = [_script_run(fn) for fn in (
        test_dual_algebra,
        test_spring_global,
        test_bar_element,
        test_blackbox_demo,
        test_blackbox_reference_template,
    )]
    ok = all(results)
    print("OVERALL: %s" % ("ALL PASS" if ok else "FAILURE"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
