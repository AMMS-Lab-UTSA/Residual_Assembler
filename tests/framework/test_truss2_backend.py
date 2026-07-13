"""Truss2 backend verification (non-solid, non-CP proof of agnosticism).

Checks:
  1. Axial pull of a single bar gives the analytical internal force EA/L * dL.
  2. Rigid translation -> zero residual (Level 1).
  3. Consistent tangent == finite-difference tangent (Level 4).
  4. The GENERIC assembler drives the truss with no truss knowledge.

Run: python tests/framework/test_truss2_backend.py
"""

import os
import sys

import numpy as np
import pytest

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from residual_core.formulations.truss2 import Truss2
from residual_core.core.model import Model, Element
from residual_core.core.dof_manager import DofManager
from residual_core.core.assembler import Assembler
from residual_core.materials.base import MaterialBinding


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


def _section(E, A):
    b = MaterialBinding(material=None, name="bar")
    b.section = {"E": E, "A": A}
    return b


def test_axial_force():
    E, A, L = 210e3, 25.0, 2.0
    coords = np.array([[0, 0, 0], [L, 0, 0]], float)
    t = Truss2()
    d = 1e-3                                   # stretch node 2 along x by d
    u = np.array([0, 0, 0, d, 0, 0], float)
    r, K, _, diag = t.eval_element(1, "T3D2", coords, u, {}, None, _section(E, A),
                                   (0, 0), 0.0, None, {"compute_tangent": True})
    N_exact = E * A / L * d
    ok = (abs(diag["axial_force"] - N_exact) < 1e-9 and
          abs(r[3] - N_exact) < 1e-9 and abs(r[0] + N_exact) < 1e-9)
    print("  axial force: N=%.6g expected=%.6g -> %s"
          % (diag["axial_force"], N_exact, "PASS" if ok else "FAIL"))
    assert ok, ("axial force: N=%.9g expected=%.9g, r[3]=%.9g, r[0]=%.9g (tol 1e-9)"
                % (diag["axial_force"], N_exact, r[3], r[0]))


def test_rigid_translation():
    coords = np.array([[0, 0, 0], [1.5, 0.5, 0.2]], float)
    t = Truss2()
    u = np.array([0.3, -0.2, 0.1] * 2, float)   # same translation both nodes
    r, _, _, _ = t.eval_element(1, "T3D2", coords, u, {}, None,
                                _section(70e3, 10.0), (0, 0), 0.0, None, {})
    ok = np.linalg.norm(r) < 1e-9
    print("  rigid translation: |r|=%.3e -> %s" % (np.linalg.norm(r), "PASS" if ok else "FAIL"))
    assert ok, "rigid translation: |r|=%.3e (need < 1e-9)" % np.linalg.norm(r)


def test_fd_tangent():
    coords = np.array([[0, 0, 0], [1.0, 1.0, 1.0]], float)
    t = Truss2()
    b = _section(200e3, 12.0)
    rng = np.random.default_rng(0)
    u = rng.uniform(-1e-3, 1e-3, 6)
    r0, K, _, _ = t.eval_element(1, "T3D2", coords, u, {}, None, b, (0, 0), 0.0,
                                 None, {"compute_tangent": True})
    Kfd = np.zeros((6, 6))
    h = 1e-7
    for j in range(6):
        up = u.copy(); up[j] += h
        rp, _, _, _ = t.eval_element(1, "T3D2", coords, up, {}, None, b, (0, 0),
                                     0.0, None, {"compute_tangent": False})
        Kfd[:, j] = (rp - r0) / h
    rel = np.linalg.norm(K - Kfd) / max(np.linalg.norm(K), 1e-30)
    ok = rel < 1e-6
    print("  FD tangent: rel=%.3e -> %s" % (rel, "PASS" if ok else "FAIL"))
    assert ok, "FD tangent: rel=%.3e (need < 1e-6)" % rel


def test_through_assembler():
    """Two colinear bars pulled at the end, driven by the generic assembler."""
    E, A = 100e3, 5.0
    m = Model()
    m.nodes = {1: (0, 0, 0), 2: (1, 0, 0), 3: (2, 0, 0)}
    m.elements = {1: Element(1, "T3D2", [1, 2]), 2: Element(2, "T3D2", [2, 3])}
    m.element_formulation = {1: "truss2", 2: "truss2"}
    m.element_material = {1: "bar", 2: "bar"}
    m.materials = {"bar": _section(E, A)}
    forms = {"truss2": Truss2()}
    dm = DofManager.for_model(m, forms)
    asm = Assembler(m, dm, forms)
    U = np.zeros(dm.ndof)
    U[dm.dof_index(3, 1)] = 1e-3               # pull node 3 in x
    R, K, diag = asm.assemble(U, compute_tangent=True)
    # internal force at the pulled node = EA/L * u3 (series of 2 unit bars, but
    # node 2 free -> here we only check the assembler ran and dof map is right)
    ok = (diag["elements"] == 2 and dm.ndof == 9 and K is not None)
    print("  assembler dispatch: elems=%d ndof=%d -> %s"
          % (diag["elements"], dm.ndof, "PASS" if ok else "FAIL"))
    assert ok, ("assembler dispatch: elems=%d (need 2), ndof=%d (need 9), "
                "tangent present=%r" % (diag["elements"], dm.ndof, K is not None))


def main():
    print("Truss2 backend verification")
    res = [_script_run(fn) for fn in (test_axial_force, test_rigid_translation,
                                      test_fd_tangent, test_through_assembler)]
    ok = all(res)
    print("OVERALL: %s" % ("ALL PASS" if ok else "FAILURE"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
