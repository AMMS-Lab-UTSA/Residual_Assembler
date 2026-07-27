"""M4 FCC checkpoint: OTI-seeded cubic (C11,C12,C44) elasticity, from real Abaqus.

`cubic_oti_umat.for` seeds the three independent cubic elastic constants
simultaneously with OTIM4N1 directions E1/E2/E3 and writes dsigma/dC11,
dsigma/dC12, dsigma/dC44 to SDV37..54. Committed fixtures under
`tests/abaqus_derivative_export/fixtures/fcc_*` are the export + FD reruns of a
real Abaqus 2021 analysis of the nonuniform single-C3D8 cubic-crystal job, so
these tests validate the C11/C12/C44 seeding OFFLINE.

Checks (all precision-independent identities are float32-floor exact):
* Euler-on-stress:  C11 dsig/dC11 + C12 dsig/dC12 + C44 dsig/dC44 = sigma.
* Analytic per-component:  dsig/dCij = (dD/dCij) . eps  (eps = Dinv . sigma).
* Euler-on-displacement:  C11 du/dC11 + C12 du/dC12 + C44 du/dC44 = -u  (exact
  because K is degree-1 homogeneous in (C11,C12,C44); this is the accurate,
  finite-difference-free sensitivity oracle).
* du/dCij vs full Abaqus PROPS-perturbation reruns (storage-precision-limited).

Run:  pytest tests/framework/test_fcc_cubic_oti.py
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

C0 = [168400.0, 121400.0, 75400.0]         # Cu cubic constants (MPa): C11, C12, C44
PARAMS = ["C11", "C12", "C44"]
INP = os.path.join(_ASSET, "fcc_cubic_c3d8.inp")


def _model_dm():
    model = from_abaqus(P.parse_inp(INP))
    return model, DofManager(model.nodes.keys())


def _uvec(d, dm):
    u = np.zeros(dm.ndof)
    for n, v in d.items():
        gd = dm.node_dofs(int(n))
        for i in range(3):
            u[gd[i]] = v[i]
    return u


def _fields():
    df = load_derivative_fields(os.path.join(_FIX, "fcc_fields.json"))
    sv = df.statev[1]                      # (8, 54)
    return df, sv


def _dD(param):
    """Structure matrix dD/dCij for cubic elasticity (Voigt, engineering shear)."""
    A = np.zeros((6, 6))
    if param == "C11":
        for i in range(3):
            A[i, i] = 1.0
    elif param == "C12":
        for i in range(3):
            for j in range(3):
                if i != j:
                    A[i, j] = 1.0
    elif param == "C44":
        for i in range(3, 6):
            A[i, i] = 1.0
    return A


# --------------------------------------------------------------------------- #
def test_fcc_euler_stress_identity():
    _df, sv = _fields()
    dC = {"C11": sv[:, 36:42], "C12": sv[:, 42:48], "C44": sv[:, 48:54]}
    S = json.load(open(os.path.join(_FIX, "fcc_S_ip.json")))
    Sabq = np.array([S[str(i + 1)] for i in range(8)])
    euler = C0[0] * dC["C11"] + C0[1] * dC["C12"] + C0[2] * dC["C44"]
    rel = np.max(np.abs(euler - Sabq)) / np.max(np.abs(Sabq))
    print("  FCC Euler stress (C.dsig/dC == sigma) rel(inf)=%.3e" % rel)
    assert rel < 1e-6                      # OTI cubic derivatives satisfy homogeneity


def test_fcc_analytic_stress_derivatives():
    _df, sv = _fields()
    D_ip = sv[:, 0:36].reshape(8, 6, 6)
    dC = {"C11": sv[:, 36:42], "C12": sv[:, 42:48], "C44": sv[:, 48:54]}
    S = json.load(open(os.path.join(_FIX, "fcc_S_ip.json")))
    Sabq = np.array([S[str(i + 1)] for i in range(8)])
    scale = np.max(np.abs(np.concatenate([dC[p] for p in PARAMS])))
    for p in PARAMS:
        A = _dD(p)
        err = 0.0
        for q in range(8):
            eps = np.linalg.solve(D_ip[q], Sabq[q])
            err = max(err, np.max(np.abs(A @ eps - dC[p][q])))
        rel = err / scale
        print("  FCC analytic dsig/d%-3s vs OTI rel(inf)=%.3e" % (p, rel))
        assert rel < 1e-6


def _solve():
    model, dm = _model_dm()
    df = load_derivative_fields(os.path.join(_FIX, "fcc_fields.json"))
    tf, sf = fields_from_statev(df.statev, df.sdv_layout, parameters=PARAMS)
    res = solve_field_sensitivities(
        model=model, tangent_fields=tf, stress_derivative_fields=sf,
        parameters=PARAMS, dof_manager=dm, integration="selective_reduced")
    return model, dm, res


def test_fcc_euler_displacement_oracle():
    _model, dm, res = _solve()
    ub = _uvec(json.load(open(os.path.join(_FIX, "fcc_U_base.json"))), dm)
    euler = sum(C0[i] * res.du_da(PARAMS[i]) for i in range(3))
    rel = np.max(np.abs(euler + ub)) / max(np.max(np.abs(ub)), 1e-30)
    print("  FCC Euler displacement (sum Ci du/dCi == -u) rel(inf)=%.3e" % rel)
    assert rel < 1e-6                      # precision-independent global oracle


def test_fcc_du_dCij_vs_abaqus_fd():
    _model, dm, res = _solve()
    fd = json.load(open(os.path.join(_FIX, "fcc_fd_U.json")))
    h = dict(zip(fd["meta"]["params"], fd["meta"]["h"]))
    for p in PARAMS:
        up = _uvec(fd["steps"][p]["p"], dm)
        um = _uvec(fd["steps"][p]["m"], dm)
        du_fd = (up - um) / (2 * h[p])
        rel = np.max(np.abs(res.du_da(p) - du_fd)) / max(np.max(np.abs(du_fd)), 1e-30)
        print("  FCC du/d%-3s vs Abaqus-FD rel(inf)=%.3e  ||fd||inf=%.3e"
              % (p, rel, np.max(np.abs(du_fd))))
        assert rel < 5e-5                  # single-precision-ODB finite-difference floor


_IN_SCRIPT = False


def test_live_abaqus_fcc_regenerates_fixture():
    """Run the live chain with the cubic OTI UMAT and confirm it reproduces the
    committed FCC fixture. Skipped unless RESASM_RUN_ABAQUS=1, Abaqus on PATH, and
    the OTI Fortran lib is present."""
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
    work = tempfile.mkdtemp(prefix="resasm_fcc_live_")
    try:
        rc = subprocess.call(["bash", os.path.join(_ASSET, "run_chain.sh"), work,
                              "cubic_oti_umat.for", "fcc_cubic_c3d8.inp",
                              "fcc_layout.json"])
        assert rc == 0, "run_chain.sh (FCC) failed (rc=%d)" % rc
        got = load_derivative_fields(os.path.join(work, "derivative_fields.json"))
        ref = load_derivative_fields(os.path.join(_FIX, "fcc_fields.json"))
        rel = np.max(np.abs(got.statev[1] - ref.statev[1])) / \
            max(np.max(np.abs(ref.statev[1])), 1e-30)
        print("  live FCC regen vs fixture statev rel(inf)=%.3e" % rel)
        assert rel < 1e-5
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main():
    global _IN_SCRIPT
    _IN_SCRIPT = True
    tests = [test_fcc_euler_stress_identity, test_fcc_analytic_stress_derivatives,
             test_fcc_euler_displacement_oracle, test_fcc_du_dCij_vs_abaqus_fd,
             test_live_abaqus_fcc_regenerates_fixture]
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
