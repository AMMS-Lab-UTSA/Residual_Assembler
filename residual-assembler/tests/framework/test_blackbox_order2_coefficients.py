"""Black-box order-2: the response must carry OTI Taylor COEFFICIENTS.

    coefficient = derivative / prod_i factorial(kappa_i)

At order 1 the factor is always 1, so a provider that returns *derivatives*
looks perfectly correct. At order >= 2 it is silently wrong by a factorial.
This file pins the contract down on a closed-form problem AND proves the test
actually catches the mistake (a deliberately-wrong provider must FAIL).

Model (deliberately nonlinear in the parameter, unlike the plain spring whose
second parameter-derivatives all vanish):

    R(u, k, f) = k^2 u^3 - f        =>   u(k, f) = (f / k^2)^(1/3)

At k = 2, f = 32  ->  u = 2, T = dR/du = 3 k^2 u^2 = 48, and exactly:

    du/dk    = -2/3        du/df    = 1/48
    d2u/dk2  =  5/9        d2u/dkdf = -1/144      d2u/df2 = -1/2304

The order-2 objects (verified with exact rational arithmetic):

    direction   recovery   R^(2) coeff   U^(2) coeff   recovered derivative
    k^2            2!         -40/3          5/18            5/9
    k*f            1            1/3         -1/144          -1/144
    f^2            2!          1/96        -1/4608         -1/2304

Note U^(2) coefficient for k^2 is 5/18 while d2u/dk2 is 5/9 -- a factor of 2.

The black-box path does NOT require OTILib (the framework only does the linear
solve), so this test runs everywhere.
"""

import json
import os
import shutil
import subprocess
import sys

import numpy as np
import pytest

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from resasm_user import run_from_config  # noqa: E402

TEMPLATE = os.path.join(_ROOT, "templates", "user_blackbox_order2_residual")

# exact values (see module docstring)
DU_DK, DU_DF = -2.0 / 3.0, 1.0 / 48.0
D2U_DK2, D2U_DKDF, D2U_DF2 = 5.0 / 9.0, -1.0 / 144.0, -1.0 / 2304.0
R2_KK, R2_KF, R2_FF = -40.0 / 3.0, 1.0 / 3.0, 1.0 / 96.0
U2_KK_COEFF = 5.0 / 18.0                     # NOT the derivative (that is 5/9)
TOL = 1e-8


def _job(tmpdir, solver_src=None):
    """Copy the shipped order-2 template into tmpdir; optionally override the solver."""
    dst = os.path.join(str(tmpdir), "job")
    shutil.copytree(TEMPLATE, dst,
                    ignore=shutil.ignore_patterns("__pycache__", "resasm_output"))
    if solver_src is not None:
        with open(os.path.join(dst, "my_solver.py"), "w") as fh:
            fh.write(solver_src)
    return run_from_config(os.path.join(dst, "resasm.yml"))


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    return _job(tmp_path_factory.mktemp("bb2"))


def _cols(res, p):
    with open(os.path.join(res.private_dir, "direction_map_order%d.json" % p)) as fh:
        return json.load(fh)["columns"]


def _idx(res, p, exps):
    for c in _cols(res, p):
        if c["exponents"] == list(exps):
            return c["index"]
    raise AssertionError("no column %s at order %d" % (exps, p))


# --------------------------------------------------------------------------- #
def test_template_runs_end_to_end_at_order_2(run):
    assert os.path.exists(os.path.join(run.private_dir,
                                       "solution_sensitivities_order2.npz"))
    assert run.summary["orders_solved"] == [1, 2]


def test_blackbox_returned_coefficients_not_derivatives(run):
    """The R^(2) the provider handed back must be Taylor COEFFICIENTS."""
    r2 = np.load(os.path.join(run.private_dir, "rhs_order2.npz"))
    Rc = r2["residual_coefficients"]
    kk, kf, ff = (_idx(run, 2, e) for e in ((2, 0), (1, 1), (0, 2)))

    assert abs(Rc[0, kk] - R2_KK) < TOL, Rc[0, kk]
    assert abs(Rc[0, kf] - R2_KF) < TOL, Rc[0, kf]
    assert abs(Rc[0, ff] - R2_FF) < TOL, Rc[0, ff]

    # if the provider had (wrongly) returned derivatives, R^(2)[k^2] would be
    # 2! * (-40/3) = -80/3. Assert it did NOT.
    assert abs(Rc[0, kk] - 2.0 * R2_KK) > 1.0


def test_direction_map_has_recovery_factor_2_for_k_squared(run):
    got = {tuple(c["exponents"]): c["recovery_factor"] for c in _cols(run, 2)}
    assert got == {(2, 0): 2.0, (1, 1): 1.0, (0, 2): 2.0}, got
    lab = {tuple(c["exponents"]): c["label"] for c in _cols(run, 2)}
    assert lab[(2, 0)] == "d2/dk2"


def test_rhs_order2_npz_has_both_coefficient_and_derivative(run):
    z = np.load(os.path.join(run.private_dir, "rhs_order2.npz"))
    for key in ("residual_coefficients", "rhs_coefficients",
                "residual_derivatives", "rhs_derivatives",
                "recovery_factors", "direction_exponents"):
        assert key in z, key
    kk = _idx(run, 2, (2, 0))
    # derivative == recovery_factor * coefficient, and they differ for k^2
    assert abs(z["residual_derivatives"][0, kk]
               - 2.0 * z["residual_coefficients"][0, kk]) < TOL
    assert abs(z["residual_derivatives"][0, kk]
               - z["residual_coefficients"][0, kk]) > 1.0
    # rhs == -residual, in both conventions
    np.testing.assert_allclose(z["rhs_coefficients"], -z["residual_coefficients"])
    np.testing.assert_allclose(z["rhs_derivatives"], -z["residual_derivatives"])


def test_solution_sensitivities_match_analytic(run):
    z1 = np.load(os.path.join(run.private_dir, "solution_sensitivities_order1.npz"))
    z2 = np.load(os.path.join(run.private_dir, "solution_sensitivities_order2.npz"))

    d1 = z1["U_derivatives"]
    assert abs(d1[0, _idx(run, 1, (1, 0))] - DU_DK) < TOL
    assert abs(d1[0, _idx(run, 1, (0, 1))] - DU_DF) < TOL

    kk, kf, ff = (_idx(run, 2, e) for e in ((2, 0), (1, 1), (0, 2)))
    d2, c2 = z2["U_derivatives"], z2["U_coefficients"]

    # RECOVERED derivatives are the true second derivatives
    assert abs(d2[0, kk] - D2U_DK2) < TOL, d2[0, kk]
    assert abs(d2[0, kf] - D2U_DKDF) < TOL, d2[0, kf]
    assert abs(d2[0, ff] - D2U_DF2) < TOL, d2[0, ff]

    # the raw COEFFICIENT for k^2 is half of it -- reading it as a derivative
    # would silently under-report by 2x
    assert abs(c2[0, kk] - U2_KK_COEFF) < TOL, c2[0, kk]
    assert abs(d2[0, kk] - 2.0 * c2[0, kk]) < TOL


def test_public_report_uses_recovered_derivatives(run):
    path = os.path.join(run.public_dir, "sensitivity_norms.csv")
    with open(path) as fh:
        rows = [l.strip().split(",") for l in fh if l.strip()]
    h, data = rows[0], rows[1:]
    i_o, i_d = h.index("order"), h.index("direction")
    i_f, i_n = h.index("recovery_factor"), h.index("solution_sensitivity_norm")

    row = next(r for r in data if r[i_o] == "2" and r[i_d] == "d2/dk2")
    assert float(row[i_f]) == 2.0
    # 1 DOF -> the norm IS |d2u/dk2| = 5/9, NOT the raw coefficient 5/18
    assert abs(float(row[i_n]) - abs(D2U_DK2)) < 1e-6, row
    assert abs(float(row[i_n]) - abs(U2_KK_COEFF)) > 0.1, "public quoted the COEFFICIENT!"


# --------------------------------------------------------------------------- #
# The guard on the guard: a provider that returns DERIVATIVES must produce a
# wrong answer. If this ever passes, the checks above have stopped meaning
# anything.
# --------------------------------------------------------------------------- #
_BUGGY_SOLVER = r'''
"""WRONG ON PURPOSE: returns recovered DERIVATIVES instead of Taylor coefficients."""
import argparse, itertools, json, math, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from my_solver_correct import T2, directions, residual, dR_du   # reuse the algebra

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--request", required=True)
    ap.add_argument("--response", required=True)
    a = ap.parse_args()
    req = json.load(open(a.request))
    u = [float(x) for x in req["u"]]
    params = req["parameters"]
    seeds = req["seed_directions"]
    names = sorted(seeds, key=lambda n: seeds[n]); m = len(names)
    p = int(req["order"]); dmap = req["direction_map"][str(p)]
    usc = req.get("u_star_coefficients") or {}
    ustar = [T2(m, u[j]) for j in range(len(u))]
    if p >= 2 and "1" in usc:
        for j in range(len(u)):
            ustar[j].o1 = [float(usc["1"][j][i]) for i in range(m)]
    pstar = {n: T2.seed(m, float(params[n]), seeds[n]-1) for n in names}
    R = residual(ustar, *[pstar[n] for n in names])
    R_p = np.zeros((len(u), len(dmap)))
    for col, exps in enumerate(dmap):
        idx = [i for i, e in enumerate(exps) for _ in range(e)]
        fac = 1.0
        for e in exps:
            fac *= math.factorial(e)
        for j in range(len(u)):
            c = R[j].o1[idx[0]] if p == 1 else R[j].get2(idx[0], idx[1])
            R_p[j, col] = fac * c          # <-- THE BUG: pre-multiplied by kappa!
    np.savez(a.response, **{"R_order_%d" % p: R_p,
                            "tangent": np.asarray(dR_du(u, *[float(params[n]) for n in names]), float)})

if __name__ == "__main__":
    main()
'''


def test_a_provider_returning_derivatives_is_WRONG(tmp_path):
    """Proves the contract is load-bearing: pre-multiplying by kappa! breaks it."""
    dst = os.path.join(str(tmp_path), "job")
    shutil.copytree(TEMPLATE, dst,
                    ignore=shutil.ignore_patterns("__pycache__", "resasm_output"))
    # keep the correct algebra importable, then overwrite the solver with the buggy one
    shutil.copy(os.path.join(dst, "my_solver.py"),
                os.path.join(dst, "my_solver_correct.py"))
    with open(os.path.join(dst, "my_solver.py"), "w") as fh:
        fh.write(_BUGGY_SOLVER)

    res = run_from_config(os.path.join(dst, "resasm.yml"))
    z2 = np.load(os.path.join(res.private_dir, "solution_sensitivities_order2.npz"))
    with open(os.path.join(res.private_dir, "direction_map_order2.json")) as fh:
        kk = next(c["index"] for c in json.load(fh)["columns"]
                  if c["exponents"] == [2, 0])

    wrong = z2["U_derivatives"][0, kk]
    # order 1 still looks fine (factor 1) -- which is exactly why this is a footgun
    z1 = np.load(os.path.join(res.private_dir, "solution_sensitivities_order1.npz"))
    assert abs(z1["U_derivatives"][0, 0] - DU_DK) < TOL, "order 1 should still be right"

    # ...but order 2 is off by exactly the factorial (2x)
    assert abs(wrong - 2.0 * D2U_DK2) < TOL, wrong
    assert abs(wrong - D2U_DK2) > 0.5, "the buggy provider was NOT caught!"
