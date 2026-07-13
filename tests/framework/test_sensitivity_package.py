"""Output contract / sensitivity-package tests (the HYPAD deliverable).

Verifies the framework emits the ingredients of  T U^(p) = -R^(p)  and that the
serializable result objects round-trip. No OTI/HYPAD algebra is exercised — the
R^(p) columns are supplied synthetically to test the contract and the solve.

Run: python tests/framework/test_sensitivity_package.py
"""

import json
import os
import sys
import tempfile
from types import SimpleNamespace

import numpy as np
import pytest

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from residual_core import ResidualProblem, SensitivityPackage
from residual_core.core.model import Model, Element
from residual_core.materials.base import MaterialBinding

_EX = os.path.join(_ROOT, "residual_core", "examples")


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
        print("  %-42s -> %s" % (label, "PASS" if v else "FAIL"))
    failed = [label for label, v in checks if not v]
    assert not failed, "failed checks: " + "; ".join(failed)


def _fixed_beam_problem():
    m = Model(nodes={1: (0, 0, 0), 2: (1, 0, 0)})
    m.elements[1] = Element(1, "B31", [1, 2])
    m.element_material = {1: "sec"}
    b = MaterialBinding(material=None, name="sec")
    b.section = {"E": 210000.0, "A": 100.0, "Iz": 833.0, "Iy": 833.0,
                 "G": 80000.0, "J": 1666.0}
    m.materials = {"sec": b}
    m.boundaries = [SimpleNamespace(target=1, dof_start=1, dof_end=6,
                                    value=0.0, kind="fixed")]
    return ResidualProblem(m)


def test_real_side_and_tangent():
    p = _fixed_beam_problem()
    res = p.result(mode="formulation")
    checks = [
        ("R has ndof=12", res.ndof == 12),
        ("tangent available", res.tangent.available),
        ("tangent source backend-assembled",
         res.tangent.source == "backend-assembled"),
        ("tangent symmetric", res.tangent.symmetric is True),
        ("dof ordering labelled", len(res.tangent.dof_labels) == 12
         and res.tangent.dof_labels[0].endswith(":UX")),
        ("element contributions captured", len(res.element_ids) == 1),
        ("free/prescribed partition", int(res.free_mask.sum()) == 6),
        ("||R_free|| ~ 0 at U=0", res.free_residual_norm < 1e-9),
    ]
    _report(checks)


def test_algebra_and_direction_maps():
    p = _fixed_beam_problem()
    pkg = p.sensitivity_package(mode="formulation",
                                parameters=["E", "A"], max_order=3)
    s = pkg.sensitivity
    # N^(p) = C(p+m-1, p) for m=2: 2, 3, 4
    checks = [
        ("parameter map E->1,A->2", s.parameter_map == {"E": 1, "A": 2}),
        ("N^(1)=2", s.algebra.n_directions(1) == 2),
        ("N^(2)=3", s.algebra.n_directions(2) == 3),
        ("N^(3)=4", s.algebra.n_directions(3) == 4),
        ("N coeffs C(2+3,2)=10", s.algebra.n_coefficients() == 10),
        ("expected R^(2) shape (12,3)", s.expected_shape(2) == (12, 3)),
    ]
    dmap2 = s.direction_map(2)
    # order-2 directions: e1^2 (factor 2), e1*e2 (factor 1), e2^2 (factor 2)
    labels = [d["label"] for d in dmap2]
    factors = [d["recovery_factor"] for d in dmap2]
    checks.append(("order-2 labels", labels == ["e1^2", "e1*e2", "e2^2"]))
    checks.append(("order-2 recovery factors", factors == [2, 1, 2]))
    _report(checks)


def test_set_rhs_and_solve():
    p = _fixed_beam_problem()
    pkg = p.sensitivity_package(mode="formulation",
                                parameters=["E", "A"], max_order=2)
    s = pkg.sensitivity
    checks = [("system not ready before R^(1)", not pkg.system_ready(1))]
    np.random.seed(1)
    R1 = np.random.rand(*s.expected_shape(1))
    s.set_residual_order(1, R1)
    checks.append(("rhs^(1) = -R^(1)", np.allclose(s.rhs(1), -R1)))
    checks.append(("system ready after R^(1)", pkg.system_ready(1)))
    U1 = pkg.solve(1)
    checks.append(("U^(1) shape (12,2)", U1.shape == (12, 2)))
    # prescribed DOFs carry zero sensitivity
    checks.append(("prescribed rows zero",
                   np.allclose(U1[~pkg.residual.free_mask], 0.0)))
    # residual of the solved system on free part: T_ff U_f = -R_f
    T = pkg.residual.tangent.T
    free = pkg.residual.free_mask
    resid = T[np.ix_(free, free)] @ U1[free] + R1[free]
    checks.append(("T_ff U_f + R_f ~ 0", np.linalg.norm(resid) < 1e-6))
    # shape guard
    try:
        s.set_residual_order(1, np.zeros((3, 3)))
        checks.append(("shape guard", False))
    except ValueError:
        checks.append(("shape guard", True))
    _report(checks)


def test_stress_driven_tangent_unavailable():
    model = os.path.join(_EX, "minimal_c3d8_stress_driven", "model.json")
    fields = os.path.join(_EX, "minimal_c3d8_stress_driven", "fields.json")
    if not os.path.exists(model):
        from residual_core.examples import generate_minimal
        generate_minimal.main()
    p = ResidualProblem.from_neutral(model)
    p.attach_results(fields)
    pkg = p.sensitivity_package(mode="stress-driven")
    checks = [
        ("runnable", pkg.runnable),
        ("tangent unavailable", not pkg.tangent_available()),
        ("system not ready", not pkg.system_ready(1)),
        ("||R_free|| > 0 (nonzero traction)", pkg.residual.free_residual_norm > 1.0),
    ]
    try:
        pkg.solve(1)
        checks.append(("solve blocked", False))
    except RuntimeError:
        checks.append(("solve blocked with message", True))
    _report(checks)


def test_not_runnable_reports_missing():
    model = os.path.join(_EX, "minimal_c3d8_stress_driven", "model.json")
    p = ResidualProblem.from_neutral(model)   # no field attached
    pkg = p.sensitivity_package(mode="stress-driven")
    ok = (not pkg.runnable) and bool(pkg.minimum_missing) and pkg.residual is None
    print("  %-42s -> %s" % ("not-runnable + minimum_missing", "PASS" if ok else "FAIL"))
    assert ok, ("runnable=%r minimum_missing=%r residual=%r (need not-runnable, "
                "a minimum_missing item, and no residual)"
                % (pkg.runnable, pkg.minimum_missing, pkg.residual))


def test_serialization_roundtrip():
    p = _fixed_beam_problem()
    pkg = p.sensitivity_package(mode="formulation", parameters=["E", "A"], max_order=1)
    pkg.sensitivity.set_residual_order(1, np.ones(pkg.sensitivity.expected_shape(1)))
    d = tempfile.mkdtemp()
    written = pkg.save(os.path.join(d, "pkg"))
    files = set(os.listdir(d))
    checks = [
        ("json written", "pkg.json" in files),
        ("npz written", "pkg.npz" in files),
        ("md written", "pkg.md" in files),
        ("sub-results written", "pkg_residual.json" in files
         and "pkg_sensitivity.npz" in files and "pkg_validation.md" in files),
    ]
    with open(os.path.join(d, "pkg.json"), encoding="utf-8") as fh:
        meta = json.load(fh)
    checks.append(("target system recorded", meta["target_system"] == "T U^(p) = -R^(p)"))
    npz = np.load(os.path.join(d, "pkg.npz"), allow_pickle=True)
    checks.append(("R array present", "residual_R" in npz))
    checks.append(("rhs_order_1 present", "sensitivity_rhs_order_1" in npz))
    checks.append(("R matches", np.allclose(npz["residual_R"], pkg.residual.R)))
    _report(checks)


def main():
    print("Output contract / sensitivity-package tests")
    results = [_script_run(fn) for fn in (
        test_real_side_and_tangent,
        test_algebra_and_direction_maps,
        test_set_rhs_and_solve,
        test_stress_driven_tangent_unavailable,
        test_not_runnable_reports_missing,
        test_serialization_roundtrip,
    )]
    ok = all(results)
    print("OVERALL: %s" % ("ALL PASS" if ok else "FAILURE"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
