"""OTILib nonlinear-spring sensitivity — order 2, single parameter.

Proves the order-by-order residual method with genuine OTILib scalars:
R(u,k,f) = k u^3 - f, u=(f/k)^(1/3), du/dk=-u/(3k), d2u/dk2=4u/(9k^2).

SKIPS cleanly (returns 0) when OTILib is not installed. When present, the residual
is genuinely evaluated with OTI numbers — no hand-coded derivative extraction.

Run: python tests/framework/test_otilib_spring_sensitivity.py
"""

import os
import sys
import tempfile

import numpy as np
import pytest

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from residual_core.algebra import otilib_adapter as A
from residual_core import ResidualProblem
from residual_core.core.model import Model, Element
from residual_core.materials.base import MaterialBinding

K, F = 2.0, 16.0
U_EXACT = (F / K) ** (1.0 / 3.0)
DUDK = -U_EXACT / (3.0 * K)               # -1/3
D2UDK2 = 4.0 * U_EXACT / (9.0 * K ** 2)   # 2/9


def _require_otilib():
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


def _spring():
    m = Model(nodes={1: (0.0, 0.0, 0.0)})
    m.elements[1] = Element(1, "SPRING1", [1])
    m.element_material = {1: "spring"}
    b = MaterialBinding(material=None, name="spring")
    b.section = {"k": K, "f": F}
    m.materials = {"spring": b}
    return ResidualProblem(m)


def test_otilib_spring_order2():
    """Order-by-order OTI sensitivities of k u^3 - f — requires OTILib."""
    _require_otilib()
    p = _spring()
    u = p.solve_newton("formulation")
    pkg = p.sensitivity_package(mode="formulation", parameters=["spring.k"],
                                max_order=2, generate_rhs=True, backend="otilib",
                                fd_check=True)
    s = pkg.sensitivity
    checks = [
        ("real solve R~0", pkg.residual.free_residual_norm < 1e-9),
        ("u = (f/k)^(1/3)", abs(u[0] - U_EXACT) < 1e-9),
        ("backend otilib", s.algebra.algebra == "otilib"),
        ("R^(1) present", s.R(1) is not None),
        ("R^(2) present", s.R(2) is not None),
    ]
    # order 1: U^(1) recovery factor 1
    U1 = pkg.solve(1)
    d1 = s.direction_map(1)[0]["recovery_factor"] * float(U1[0, 0])
    checks.append(("du/dk = -u/(3k)", abs(d1 - DUDK) < 1e-8))
    # order 2: recovery factor 2 for direction [2]
    U2 = pkg.solve(2)
    d2 = s.direction_map(2)[0]["recovery_factor"] * float(U2[0, 0])
    checks.append(("d2u/dk2 = 4u/(9k^2)", abs(d2 - D2UDK2) < 1e-7))
    # FD cross-check on first derivative present + passed
    fd = [c for c in pkg.validation.checks if c.name == "sensitivity-vs-finite-difference"]
    checks.append(("FD cross-check passed", bool(fd) and fd[0].passed))
    # serialization still works
    d = tempfile.mkdtemp()
    pkg.save(os.path.join(d, "spring_oti"))
    checks.append(("serialization", os.path.exists(os.path.join(d, "spring_oti.json"))
                   and os.path.exists(os.path.join(d, "spring_oti_sensitivity.npz"))))

    for label, v in checks:
        print("  %-40s -> %s" % (label, "PASS" if v else "FAIL"))
    failed = [label for label, v in checks if not v]
    assert not failed, "failed checks: " + "; ".join(failed)


def main():
    print("OTILib nonlinear-spring sensitivity (order 2)")
    ok = _script_run(test_otilib_spring_order2)
    if ok and not A.otilib_available():
        print("OVERALL: SKIPPED (OTILib unavailable)")
    else:
        print("OVERALL: %s" % ("ALL PASS" if ok else "FAILURE"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
