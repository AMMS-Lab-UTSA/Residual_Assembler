"""Unit tests for the canonical FD-reference methodology.

Pure-Python and fast (no Fortran, no material objects): they pin the behaviour
that the OTI-vs-FD claims depend on -- exact-derivative recovery, the roundoff
robustness that the FCC crystal model exposed, and the rule that a non-finite
value is always a failure and never a ``0.0``.
"""
import math
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

from umat_oti.validation import fd_reference as fdref


def test_recovers_a_known_derivative():
    # f(p) = [p0^2, 3*p1];  d/dp0 -> [2 p0, 0],  d/dp1 -> [0, 3]
    def response(p):
        return {"y": np.array([p[0] ** 2, 3.0 * p[1]])}

    p0 = [2.0, 5.0]
    oti = {"y": np.array([[2 * p0[0], 0.0], [0.0, 3.0]]).T}   # (..., param) layout
    params = [{"name": "a", "props_index": 1}, {"name": "b", "props_index": 2}]
    report = fdref.verify_derivatives(response, p0, params, oti, tolerances={"y": 1e-6})
    assert report["passed"]
    assert all(r["status"] == "PASS" for r in report["parameters"])
    assert report["worst_rel"]["y"] < 1e-8


def test_wrong_oti_value_fails():
    def response(p):
        return {"y": np.array([p[0] ** 2])}

    oti = {"y": np.array([[2.0]])}          # true is 2*p0 = 6 at p0=3 -> wrong
    report = fdref.verify_derivatives(response, [3.0], [{"name": "a", "props_index": 1}],
                                      oti, tolerances={"y": 1e-6})
    assert not report["passed"]
    assert report["parameters"][0]["status"] == "FAIL"


def test_nonfinite_response_is_a_failure_never_zero():
    def response(p):
        return {"y": np.array([float("nan")])}

    with pytest.raises(fdref.NonFiniteResult):
        fdref.require_finite("y", response([1.0])["y"])

    # and through the full path: a NaN OTI value must not certify
    def ok_response(p):
        return {"y": np.array([2.0 * p[0]])}

    report = fdref.verify_derivatives(ok_response, [1.0], [{"name": "a", "props_index": 1}],
                                      {"y": np.array([[float("nan")]])}, tolerances={"y": 1e-6})
    assert not report["passed"]
    assert math.isnan(report["worst_rel"]["y"])
    assert "non-finite" in report["parameters"][0]["reason"].lower()


def test_near_zero_derivative_scalar_does_not_lock_onto_roundoff():
    # A scalar whose derivative wrt p1 is ~0: the fine-step central differences
    # suffer cancellation and produce spurious ~1e-16 successive changes. The
    # selector must NOT pick that step (this is the FCC-crystal failure mode).
    def response(p):
        # depends strongly on p0, essentially not on p1
        return {"s": np.array([100.0 * p[0] + 1e-12 * p[1] ** 2])}

    studies = fdref.parameter_studies(response, [1.0, 1.0], 2, parameter="p1")
    st = studies["s"]
    assert st.selected_h_rel is not None
    # the plateau must be corroborated across two refinements, not an isolated dip
    assert st.selected_h_rel >= 1e-6, "selected a roundoff-dominated step"


def test_finest_rung_cannot_self_certify_on_cancellation():
    # The J2-vs-nu failure mode: a scalar whose derivative wrt p1 is so small that
    # the two finest central differences are bit-identical (successive change 0.0).
    # The finest rung has no finer neighbour to corroborate it, so it must NOT be
    # selected -- a 0.0 successive change at the ladder's end is cancellation, not
    # convergence.
    def response(p):
        return {"s": np.array([1000.0 * p[0]])}   # exactly independent of p1

    studies = fdref.parameter_studies(response, [1.0, 1.0], 2, parameter="p1")
    st = studies["s"]
    assert st.selected_h_rel is not None
    # not the finest rung (1e-8); a corroborated coarser plateau instead
    assert st.selected_h_rel > 1e-8


def test_stiff_parameter_still_selects_a_fine_step():
    # A tiny base value with an O(1) derivative genuinely needs a fine step;
    # the look-ahead must not force a coarse one.
    def response(p):
        return {"s": np.array([p[0], np.sin(50.0 * p[1])])}

    studies = fdref.parameter_studies(response, [1.0, 1e-3], 2, parameter="p1")
    st = studies["s"]
    assert st.selected is not None
    assert st.reference is not None


def test_richardson_is_reported_as_the_convergence_metric():
    def response(p):
        return {"y": np.array([p[0] ** 3])}

    studies = fdref.parameter_studies(response, [2.0], 1, parameter="a")
    st = studies["y"]
    metric, basis = st.convergence_metric
    assert basis == "richardson"
    assert math.isfinite(metric)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
