"""Beam2 backend verification (rotational DOFs -> heterogeneous DOF proof).

Checks:
  1. Cantilever with a transverse tip load: tip deflection == P L^3 / (3 E I)
     (Euler-Bernoulli closed form) by solving K_ff u_f = f on the free DOFs.
  2. Rigid translation -> zero residual (Level 1).
  3. Consistent tangent == finite-difference tangent (Level 4).

Run: python tests/framework/test_beam2_backend.py
"""

import os
import sys

import numpy as np
import pytest

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from residual_core.formulations.beam2 import Beam2
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


def _section(**kw):
    b = MaterialBinding(material=None, name="sec")
    b.section = dict(kw)
    return b


def test_cantilever_tip_deflection():
    """Single beam element, fixed at node 1, transverse point load P at node 2.
    Tip deflection in local/global y = P L^3 / (3 E Iz)."""
    E, A, Iz, L, P = 210e3, 100.0, 833.33, 3.0, 10.0
    coords = np.array([[0, 0, 0], [L, 0, 0]], float)
    beam = Beam2()
    b = _section(E=E, A=A, Iz=Iz, Iy=Iz)
    K, _ = beam.global_stiffness(coords, b)

    # DOFs 0..5 node1 (fixed), 6..11 node2 (free). Apply Fy at dof 7.
    free = [6, 7, 8, 9, 10, 11]
    f = np.zeros(12); f[7] = P
    Kff = K[np.ix_(free, free)]
    uf = np.linalg.solve(Kff, f[free])
    tip_v = uf[free.index(7)]
    exact = P * L**3 / (3.0 * E * Iz)
    rel = abs(tip_v - exact) / abs(exact)
    ok = rel < 1e-9
    print("  cantilever tip: v=%.6g exact=%.6g rel=%.2e -> %s"
          % (tip_v, exact, rel, "PASS" if ok else "FAIL"))
    assert ok, ("cantilever tip: v=%.6g exact=%.6g rel=%.2e (need < 1e-9)"
                % (tip_v, exact, rel))


def test_rigid_translation():
    coords = np.array([[0, 0, 0], [2.0, 1.0, 0.5]], float)
    beam = Beam2()
    b = _section(E=70e3, A=10.0, Iz=100.0)
    u = np.zeros(12)
    trans = np.array([0.2, -0.1, 0.05])
    u[0:3] = trans; u[6:9] = trans          # pure translation, no rotation
    r, _, _, _ = beam.eval_element(1, "B31", coords, u, {}, None, b, (0, 0), 0.0,
                                   None, {})
    ok = np.linalg.norm(r) < 1e-7
    print("  rigid translation: |r|=%.3e -> %s" % (np.linalg.norm(r), "PASS" if ok else "FAIL"))
    assert ok, "rigid translation: |r|=%.3e (need < 1e-7)" % np.linalg.norm(r)


def test_fd_tangent():
    coords = np.array([[0, 0, 0], [1.0, 0.7, 0.3]], float)
    beam = Beam2()
    b = _section(E=200e3, A=12.0, Iz=90.0, Iy=60.0, J=120.0)
    rng = np.random.default_rng(2)
    u = rng.uniform(-1e-4, 1e-4, 12)
    r0, K, _, _ = beam.eval_element(1, "B31", coords, u, {}, None, b, (0, 0), 0.0,
                                    None, {"compute_tangent": True})
    Kfd = np.zeros((12, 12))
    h = 1e-7
    for j in range(12):
        up = u.copy(); up[j] += h
        rp, _, _, _ = beam.eval_element(1, "B31", coords, up, {}, None, b, (0, 0),
                                        0.0, None, {"compute_tangent": False})
        Kfd[:, j] = (rp - r0) / h
    rel = np.linalg.norm(K - Kfd) / max(np.linalg.norm(K), 1e-30)
    ok = rel < 1e-5
    print("  FD tangent: rel=%.3e -> %s" % (rel, "PASS" if ok else "FAIL"))
    assert ok, "FD tangent: rel=%.3e (need < 1e-5)" % rel


def main():
    print("Beam2 backend verification")
    res = [_script_run(fn) for fn in (test_cantilever_tip_deflection,
                                      test_rigid_translation, test_fd_tangent)]
    ok = all(res)
    print("OVERALL: %s" % ("ALL PASS" if ok else "FAILURE"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
