"""OTILib adapter tests.

The direction-map / coefficient-count checks run ALWAYS (backend-independent
combinatorics — acceptance item 6). The scalar/arithmetic/extraction checks
require a genuine OTILib install and SKIP cleanly when it is absent (the PyPI
`pyoti` squat is NOT accepted).

Run: python tests/framework/test_otilib_adapter.py
"""

import math
import os
import sys

import pytest

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from residual_core.algebra import otilib_adapter as A


def _require_otilib():  # noqa: D401
    """Skip when OTILib is absent — unless RUN_OTILIB_TESTS=1 demands it."""
    if A.otilib_available():
        return
    assert os.environ.get("RUN_OTILIB_TESTS") != "1", (
        "RUN_OTILIB_TESTS=1 but OTILib is not installed; install: "
        "scripts/setup_otilib.sh (https://github.com/mauriaristi/otilib.git)")
    pytest.skip("OTILib not installed: %s"
                % A.otilib_status()["error"].splitlines()[0])


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


def test_direction_counts():
    checks = []
    for m in (1, 2, 3):
        for nt in (1, 2, 3):
            n_total = A.num_coefficients_total(m, nt)
            checks.append(("N_total m=%d nt=%d = C(m+nt,m)" % (m, nt),
                           n_total == math.comb(m + nt, m)))
        for p in (1, 2, 3):
            n_p = A.num_coefficients_order(m, p)
            checks.append(("N_order m=%d p=%d = C(p+m-1,p)" % (m, p),
                           n_p == math.comb(p + m - 1, p)))
    _report(checks)


def test_direction_map_and_recovery():
    # m=3, order 2 canonical directions + recovery factors
    dirs = A.order_directions(3, 2)
    exps = [d["exponents"] for d in dirs]
    labels = [d["label"] for d in dirs]
    factors = [d["recovery_factor"] for d in dirs]
    checks = [
        ("order-2 count = 6", len(dirs) == 6),
        ("exponents", exps == [(2, 0, 0), (1, 1, 0), (1, 0, 1),
                               (0, 2, 0), (0, 1, 1), (0, 0, 2)]),
        ("labels", labels == ["e1^2", "e1*e2", "e1*e3", "e2^2", "e2*e3", "e3^2"]),
        ("recovery factors (prod kappa!)", factors == [2, 1, 1, 2, 1, 2]),
    ]
    # order 3 single variable: [3] recovery 3! = 6
    d3 = A.order_directions(1, 3)
    checks.append(("m=1 order-3 recovery = 6", d3[0]["recovery_factor"] == 6))
    _report(checks)


def test_otilib_scalars():
    _require_otilib()
    ctx = A.OtiContext(num_bases=2, order=2)
    checks = []
    x = ctx.scalar(3.0)
    checks.append(("scalar real", abs(ctx.real(x) - 3.0) < 1e-12))
    s = ctx.seed(2.0, 1)                      # 2 + e1
    checks.append(("seed real", abs(ctx.real(s) - 2.0) < 1e-12))
    checks.append(("seed e1 coeff = 1", abs(ctx.coeff(s, (1, 0)) - 1.0) < 1e-12))
    # arithmetic: (2+e1)^3 -> real 8, d/de1 = 3*2^2 = 12 along [1,0]
    y = s ** 3
    checks.append(("(2+e1)^3 real = 8", abs(ctx.real(y) - 8.0) < 1e-9))
    checks.append(("(2+e1)^3 e1 coeff = 12", abs(ctx.coeff(y, (1, 0)) - 12.0) < 1e-9))
    # second order: e1^2 coeff = (1/2!) d2/de1^2 (x^3) = (1/2)*6*2 = 6
    checks.append(("(2+e1)^3 e1^2 coeff = 6", abs(ctx.coeff(y, (2, 0)) - 6.0) < 1e-9))
    # multiplication / division sanity
    checks.append(("mul", abs(ctx.real(ctx.seed(2.0, 1) * ctx.seed(3.0, 2)) - 6.0) < 1e-9))
    _report(checks)


def main():
    print("OTILib adapter tests (counts always; scalars gated on OTILib)")
    results = [_script_run(fn) for fn in (test_direction_counts,
                                          test_direction_map_and_recovery,
                                          test_otilib_scalars)]
    ok = all(results)
    print("OVERALL: %s" % ("ALL PASS" if ok else "FAILURE"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
