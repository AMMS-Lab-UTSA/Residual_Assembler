"""Order-2 spring: raw OTI coefficients vs recovered partial derivatives.

The OTI evaluation yields Taylor COEFFICIENTS. The true partial derivative is

    derivative = recovery_factor * coefficient ,   recovery_factor = prod_i (kappa_i!)

At order 1 every factor is 1, so nothing changes. At order >= 2 a REPEATED
direction is scaled by a factorial, and forgetting it silently under-reports the
derivative. This test pins that down on a closed-form problem.

Spring:   R(u, k, f) = k u^3 - f      =>   u(k, f) = (f/k)^(1/3)

At k = 2, f = 16  ->  u = 2, and analytically:

    du/dk    = -u / (3k)            = -1/3
    du/df    =  1 / (3 k u^2)       =  1/24
    d2u/dk2  =  (4/9) * u / k^2     =  2/9      <-- repeated direction, factor 2! = 2
    d2u/df2  = -2 / (9 k^2 u^5)     = -1/576    <-- repeated direction, factor 2! = 2
    d2u/dkdf                        -> checked against FD of u(k,f), factor 1

Therefore the RAW coefficient for the k^2 direction must be (1/2!) * 2/9 = 1/9,
and the RECOVERED derivative must be 2 * 1/9 = 2/9.

Requires OTILib (the Python path). Skips cleanly when it is unavailable; see
scripts/run_otilib_tests_wsl.sh.
"""

import json
import os
import shutil
import sys
import tempfile

import numpy as np
import pytest

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from residual_core.algebra.otilib_adapter import otilib_available  # noqa: E402
from resasm_user import run_from_config  # noqa: E402

pytestmark = pytest.mark.skipif(
    not otilib_available(),
    reason="OTILib not installed (Python path needs it; see scripts/run_otilib_tests_wsl.sh)")

_RESIDUAL = '''
def residual(u, params, state=None, time=None):
    k = params["k"]; f = params["f"]
    return [k * u[0] ** 3 - f]

def tangent(u, params, state=None, time=None):
    k = params["k"]
    return [[3.0 * k * u[0] ** 2]]
'''

_CFG = """
problem:
  name: recovery_factor_spring
residual:
  type: python
  module: user_residual.py
  function: residual
tangent:
  type: python
  function: tangent
parameters:
  k: 2.0
  f: 16.0
solution:
  file: solution.npy
sensitivity:
  order: 2
  backend: otilib
"""

# analytic, at k=2 f=16 u=2
DU_DK = -1.0 / 3.0
DU_DF = 1.0 / 24.0
D2U_DK2 = 2.0 / 9.0
D2U_DF2 = -1.0 / 576.0


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    d = str(tmp_path_factory.mktemp("recov"))
    with open(os.path.join(d, "user_residual.py"), "w") as fh:
        fh.write(_RESIDUAL)
    with open(os.path.join(d, "resasm.yml"), "w") as fh:
        fh.write(_CFG)
    np.save(os.path.join(d, "solution.npy"), np.array([2.0]))
    res = run_from_config(os.path.join(d, "resasm.yml"))
    return res


def _cols(res, p):
    with open(os.path.join(res.private_dir, "direction_map_order%d.json" % p)) as fh:
        return json.load(fh)["columns"]


def _col_of(res, p, exponents):
    for c in _cols(res, p):
        if c["exponents"] == list(exponents):
            return c
    raise AssertionError("no column with exponents %s at order %d" % (exponents, p))


def test_direction_map_is_exported_with_recovery_factors(run):
    cols = _cols(run, 2)
    # parameters are (k, f) -> order-2 directions: k^2, k*f, f^2
    got = {tuple(c["exponents"]): c["recovery_factor"] for c in cols}
    assert got == {(2, 0): 2.0, (1, 1): 1.0, (0, 2): 2.0}, got
    labels = {tuple(c["exponents"]): c["label"] for c in cols}
    assert labels[(2, 0)] == "d2/dk2"
    assert labels[(1, 1)] == "d2/dk_df"
    assert labels[(0, 2)] == "d2/df2"


def test_order1_unchanged_coefficients_equal_derivatives(run):
    z = np.load(os.path.join(run.private_dir, "solution_sensitivities_order1.npz"))
    # legacy key preserved
    assert "U" in z
    coeff = z["U_coefficients"]
    deriv = z["U_derivatives"]
    np.testing.assert_allclose(z["recovery_factors"], [1.0, 1.0])
    # at order 1 the two conventions coincide
    np.testing.assert_allclose(coeff, deriv, rtol=0, atol=1e-14)
    np.testing.assert_allclose(z["U"], coeff, rtol=0, atol=1e-14)
    ck = _col_of(run, 1, (1, 0))["index"]
    cf = _col_of(run, 1, (0, 1))["index"]
    assert abs(deriv[0, ck] - DU_DK) < 1e-8, deriv[0, ck]
    assert abs(deriv[0, cf] - DU_DF) < 1e-8, deriv[0, cf]


def test_order2_repeated_direction_needs_the_factorial(run):
    """THE point of this file: coefficient(k^2) = 1/9, derivative(k^2) = 2/9."""
    z = np.load(os.path.join(run.private_dir, "solution_sensitivities_order2.npz"))
    ck2 = _col_of(run, 2, (2, 0))["index"]
    cf2 = _col_of(run, 2, (0, 2))["index"]
    ckf = _col_of(run, 2, (1, 1))["index"]

    coeff = z["U_coefficients"]
    deriv = z["U_derivatives"]

    # raw coefficient is the derivative divided by 2!
    assert abs(coeff[0, ck2] - D2U_DK2 / 2.0) < 1e-8, coeff[0, ck2]   # 1/9
    # recovered derivative is the real second derivative
    assert abs(deriv[0, ck2] - D2U_DK2) < 1e-8, deriv[0, ck2]         # 2/9
    # and they genuinely differ (this is the bug the export prevents)
    assert abs(deriv[0, ck2] - coeff[0, ck2]) > 0.1

    assert abs(coeff[0, cf2] - D2U_DF2 / 2.0) < 1e-8, coeff[0, cf2]
    assert abs(deriv[0, cf2] - D2U_DF2) < 1e-8, deriv[0, cf2]

    # mixed direction: factor 1, so both conventions agree
    assert abs(deriv[0, ckf] - coeff[0, ckf]) < 1e-14


def test_derivatives_match_finite_difference_of_the_exact_solution(run):
    """Independent check: FD the closed-form u(k,f) and compare to U_derivatives."""
    def u_exact(k, f):
        return (f / k) ** (1.0 / 3.0)

    h = 1e-4
    fd_kk = (u_exact(2 + h, 16) - 2 * u_exact(2, 16) + u_exact(2 - h, 16)) / h ** 2
    fd_ff = (u_exact(2, 16 + h) - 2 * u_exact(2, 16) + u_exact(2, 16 - h)) / h ** 2
    fd_kf = (u_exact(2 + h, 16 + h) - u_exact(2 + h, 16 - h)
             - u_exact(2 - h, 16 + h) + u_exact(2 - h, 16 - h)) / (4 * h ** 2)

    z = np.load(os.path.join(run.private_dir, "solution_sensitivities_order2.npz"))
    deriv = z["U_derivatives"]
    assert abs(deriv[0, _col_of(run, 2, (2, 0))["index"]] - fd_kk) < 1e-5
    assert abs(deriv[0, _col_of(run, 2, (0, 2))["index"]] - fd_ff) < 1e-5
    assert abs(deriv[0, _col_of(run, 2, (1, 1))["index"]] - fd_kf) < 1e-5


def test_public_report_quotes_recovered_derivatives(run):
    """public/sensitivity_norms.csv must carry derivative norms + the factor used."""
    path = os.path.join(run.public_dir, "sensitivity_norms.csv")
    with open(path) as fh:
        rows = [l.strip().split(",") for l in fh if l.strip()]
    header, data = rows[0], rows[1:]
    assert "recovery_factor" in header
    i_ord, i_dir = header.index("order"), header.index("direction")
    i_fac, i_norm = header.index("recovery_factor"), header.index("solution_sensitivity_norm")

    row = next(r for r in data if r[i_ord] == "2" and r[i_dir] == "d2/dk2")
    assert float(row[i_fac]) == 2.0
    # 1-DOF problem -> the norm IS |d2u/dk2| = 2/9, NOT the raw coefficient 1/9
    assert abs(float(row[i_norm]) - abs(D2U_DK2)) < 1e-6, row   # %.6e formatting
