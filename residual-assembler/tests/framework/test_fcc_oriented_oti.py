"""M4 oriented-FCC checkpoint: OTI-seeded cubic elasticity with crystal ORIENTATION.

`cubic_oriented_oti_umat.for` rotates the cubic single crystal by a non-trivial
Bunge orientation (phi1,PHI,phi2 = 25,40,15 deg) using a single OTI-rotated
operator for STRESS, DDSDDE and all three dsigma/dCij. This proves crystal-
orientation handling -- the last cheap checkpoint before crystal plasticity. An
aligned-[100] test (test_fcc_cubic_oti.py) is blind to a bad rotation; these
tests validate the rotation against an INDEPENDENT 4th-order tensor rotation in
Python (Cs_ijkl = G_mi G_nj G_ka G_lb Cc_mnkl), which uses a different method than
the UMAT's rotate-strain-in / rotate-stress-out operator.

Fixtures `ori_*` are the export + PROPS-perturbation reruns of a real Abaqus 2021
analysis of the oriented single-C3D8 cubic crystal.

Checks (float32-floor exact for the identities):
* DDSDDE (rotated) vs the independent Python 4th-order rotation.
* Euler-on-stress and per-component analytic dsig/dCij (rotated structure tensors).
* Euler-on-displacement oracle sum Ci du/dCi = -u.
* du/dCij vs full Abaqus PROPS-perturbation reruns (storage-precision-limited).
"""

import json
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, os.pardir, os.pardir))
_ASSET = os.path.join(_ROOT, "tests", "abaqus_derivative_export")
_FIX = os.path.join(_ASSET, "fixtures")
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from residual_core.io import abaqus_inp_parser as P
from residual_core.core.model import from_abaqus
from residual_core.core.dof_manager import DofManager
from residual_core.io.derivative_fields import load_derivative_fields
from residual_core.core.field_sensitivity import (
    fields_from_statev, solve_field_sensitivities)

C0 = [168400.0, 121400.0, 75400.0]
PARAMS = ["C11", "C12", "C44"]
EULER_DEG = (25.0, 40.0, 15.0)
INP = os.path.join(_ASSET, "fcc_cubic_c3d8.inp")
_VP = {0: (0, 0), 1: (1, 1), 2: (2, 2), 3: (0, 1), 4: (0, 2), 5: (1, 2)}


# --------------------------------------------------------------------------- #
# independent orientation math (different method than the UMAT)
# --------------------------------------------------------------------------- #
def _G():
    p1, ph, p2 = np.radians(EULER_DEG)

    def Rz(a):
        c, s = np.cos(a), np.sin(a)
        return np.array([[c, s, 0], [-s, c, 0], [0, 0, 1.0]])

    def Rx(a):
        c, s = np.cos(a), np.sin(a)
        return np.array([[1.0, 0, 0], [0, c, s], [0, -s, c]])

    return Rz(p2) @ Rx(ph) @ Rz(p1)          # Bunge g = Rz(phi2) Rx(PHI) Rz(phi1)


def _rotate_to_voigt(Cc):
    """4th-order rotation Cs_ijab = G_mi G_nj G_ka G_lb Cc_mnkl -> 6x6 Voigt."""
    G = _G()
    Cs = np.einsum("mi,nj,ka,lb,mnkl->ijab", G, G, G, G, Cc)
    return np.array([[Cs[_VP[I][0], _VP[I][1], _VP[J][0], _VP[J][1]]
                      for J in range(6)] for I in range(6)])


def _cubic_tensor(c11, c12, c44):
    Cc = np.zeros((3, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    if i == j == k == l:
                        Cc[i, j, k, l] = c11
                    elif i == j and k == l:
                        Cc[i, j, k, l] = c12
                    elif (i == k and j == l) or (i == l and j == k):
                        Cc[i, j, k, l] = c44
    return Cc


def _struct_voigt(param):
    """Rotated dD/dCij structure matrix (6x6) via the independent 4th-order path."""
    c = {"C11": (1, 0, 0), "C12": (0, 1, 0), "C44": (0, 0, 1)}[param]
    return _rotate_to_voigt(_cubic_tensor(*c))


def _uvec(d, dm):
    u = np.zeros(dm.ndof)
    for n, v in d.items():
        gd = dm.node_dofs(int(n))
        for i in range(3):
            u[gd[i]] = v[i]
    return u


def _fields():
    df = load_derivative_fields(os.path.join(_FIX, "ori_fields.json"))
    return df, df.statev[1]


# --------------------------------------------------------------------------- #
def test_oriented_ddsdde_vs_independent_rotation():
    _df, sv = _fields()
    D_odb = sv[:, 0:36].reshape(8, 6, 6)
    Dpy = _rotate_to_voigt(_cubic_tensor(*C0))
    rel = np.max(np.abs(D_odb[0] - Dpy)) / np.max(np.abs(Dpy))
    # rotation must be genuinely non-trivial (normal-shear coupling present)
    assert abs(Dpy[0, 3]) > 1e3 and abs(Dpy[3, 4]) > 1e3
    print("  ORI DDSDDE vs independent 4th-order rotation rel(inf)=%.3e" % rel)
    assert rel < 1e-6


def test_oriented_euler_stress_identity():
    _df, sv = _fields()
    dC = {"C11": sv[:, 36:42], "C12": sv[:, 42:48], "C44": sv[:, 48:54]}
    S = json.load(open(os.path.join(_FIX, "ori_S_ip.json")))
    Sabq = np.array([S[str(i + 1)] for i in range(8)])
    euler = C0[0] * dC["C11"] + C0[1] * dC["C12"] + C0[2] * dC["C44"]
    rel = np.max(np.abs(euler - Sabq)) / np.max(np.abs(Sabq))
    print("  ORI Euler stress rel(inf)=%.3e" % rel)
    assert rel < 1e-6


def test_oriented_analytic_stress_derivatives():
    _df, sv = _fields()
    D_ip = sv[:, 0:36].reshape(8, 6, 6)
    dC = {"C11": sv[:, 36:42], "C12": sv[:, 42:48], "C44": sv[:, 48:54]}
    S = json.load(open(os.path.join(_FIX, "ori_S_ip.json")))
    Sabq = np.array([S[str(i + 1)] for i in range(8)])
    scale = np.max(np.abs(np.concatenate([dC[p] for p in PARAMS])))
    for p in PARAMS:
        A = _struct_voigt(p)
        err = max(np.max(np.abs(A @ np.linalg.solve(D_ip[q], Sabq[q]) - dC[p][q]))
                  for q in range(8))
        print("  ORI analytic dsig/d%-3s (rotated) rel(inf)=%.3e" % (p, err / scale))
        assert err / scale < 1e-6


def _solve():
    model = from_abaqus(P.parse_inp(INP))
    dm = DofManager(model.nodes.keys())
    df = load_derivative_fields(os.path.join(_FIX, "ori_fields.json"))
    tf, sf = fields_from_statev(df.statev, df.sdv_layout, parameters=PARAMS)
    res = solve_field_sensitivities(
        model=model, tangent_fields=tf, stress_derivative_fields=sf,
        parameters=PARAMS, dof_manager=dm, integration="selective_reduced")
    return dm, res


def test_oriented_euler_displacement_oracle():
    dm, res = _solve()
    ub = _uvec(json.load(open(os.path.join(_FIX, "ori_U_base.json"))), dm)
    euler = sum(C0[i] * res.du_da(PARAMS[i]) for i in range(3))
    rel = np.max(np.abs(euler + ub)) / max(np.max(np.abs(ub)), 1e-30)
    print("  ORI Euler displacement (sum Ci du/dCi == -u) rel(inf)=%.3e" % rel)
    assert rel < 1e-6


def test_oriented_du_dCij_vs_abaqus_fd():
    dm, res = _solve()
    fd = json.load(open(os.path.join(_FIX, "ori_fd_U.json")))
    h = dict(zip(fd["meta"]["params"], fd["meta"]["h"]))
    for p in PARAMS:
        du_fd = (_uvec(fd["steps"][p]["p"], dm) - _uvec(fd["steps"][p]["m"], dm)) / (2 * h[p])
        rel = np.max(np.abs(res.du_da(p) - du_fd)) / max(np.max(np.abs(du_fd)), 1e-30)
        print("  ORI du/d%-3s vs Abaqus-FD rel(inf)=%.3e" % (p, rel))
        assert rel < 5e-5


_IN_SCRIPT = False


def test_live_abaqus_oriented_regenerates_fixture():
    import shutil
    import subprocess
    import tempfile
    oti_lib = os.path.expanduser(
        os.environ.get("OTI_DIR", "~/MultiZ_f/oti")) + "/libotim4n1.a"
    if not (shutil.which(os.environ.get("ABAQUS_CMD", "abaqus"))
            and os.environ.get("RESASM_RUN_ABAQUS") and os.path.exists(oti_lib)):
        msg = "needs RESASM_RUN_ABAQUS=1, Abaqus on PATH, and the OTI lib"
        if _IN_SCRIPT:
            print("  SKIP: %s" % msg)
            return
        import pytest
        pytest.skip(msg)
    work = tempfile.mkdtemp(prefix="resasm_ori_live_")
    try:
        rc = subprocess.call(["bash", os.path.join(_ASSET, "run_chain.sh"), work,
                              "cubic_oriented_oti_umat.for", "fcc_cubic_c3d8.inp",
                              "fcc_layout.json"])
        assert rc == 0
        got = load_derivative_fields(os.path.join(work, "derivative_fields.json"))
        ref = load_derivative_fields(os.path.join(_FIX, "ori_fields.json"))
        rel = np.max(np.abs(got.statev[1] - ref.statev[1])) / \
            max(np.max(np.abs(ref.statev[1])), 1e-30)
        print("  live oriented regen vs fixture rel(inf)=%.3e" % rel)
        assert rel < 1e-5
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main():
    global _IN_SCRIPT
    _IN_SCRIPT = True
    tests = [test_oriented_ddsdde_vs_independent_rotation,
             test_oriented_euler_stress_identity,
             test_oriented_analytic_stress_derivatives,
             test_oriented_euler_displacement_oracle,
             test_oriented_du_dCij_vs_abaqus_fd,
             test_live_abaqus_oriented_regenerates_fixture]
    ok = True
    for t in tests:
        try:
            t()
            print("  PASS: %s" % t.__name__)
        except Exception as exc:            # noqa: BLE001
            ok = False
            print("  FAIL: %s -> %s" % (t.__name__, exc))
    print("OVERALL: %s" % ("ALL PASS" if ok else "FAILURE"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
