"""Reductions, request validation and the weighted-share definition."""
import numpy as np
import pytest

from residual_core.replay.history import HistoryEngine, run_history
from residual_core.replay.history_outputs import (Fields, RequestError, reduce, select,
                                                  validate_request, weighted_shares)

from conftest import beam_model


def test_reductions_and_their_gradients():
    values = np.array([1.0, 3.0, 2.0])
    derivatives = np.array([[1.0, 0.0], [0.0, 2.0], [4.0, 4.0]])
    weights = np.array([1.0, 1.0, 2.0])
    assert reduce(values, derivatives, weights, "sum")[0] == 6.0
    value, gradient, _ = reduce(values, derivatives, weights, "volume_mean")
    assert value == pytest.approx((1 + 3 + 4) / 4)
    np.testing.assert_allclose(gradient, [(1 + 0 + 8) / 4, (0 + 2 + 8) / 4])
    value, gradient, info = reduce(values, derivatives, weights, "max")
    assert value == 3.0 and gradient.tolist() == [0.0, 2.0] and info["arg_locations"] == 1
    value, gradient, _ = reduce(values, derivatives, weights, "min")
    assert value == 1.0 and gradient.tolist() == [1.0, 0.0]
    value, gradient, _ = reduce(values, derivatives, weights, "L2")
    assert value == pytest.approx(np.sqrt(14.0))
    np.testing.assert_allclose(gradient, values @ derivatives / np.sqrt(14.0))


def test_tied_extremes_need_equal_derivatives():
    values = np.array([2.0, 2.0, 1.0])
    same = np.array([[1.0, 5.0], [1.0, 5.0], [0.0, 0.0]])
    value, gradient, info = reduce(values, same, None, "max")
    assert info["arg_locations"] == 2 and gradient.tolist() == [1.0, 5.0]
    different = np.array([[1.0, 5.0], [2.0, 5.0], [0.0, 0.0]])
    with pytest.raises(RequestError, match="not differentiable"):
        reduce(values, different, None, "max")
    with pytest.raises(RequestError, match="exactly one"):
        reduce(values, same, None, "component")


@pytest.mark.parametrize("change, message", [
    (lambda r: r.pop("increments"), "lacks"),
    (lambda r: r.update(extra=1), "unsupported request keys"),
    (lambda r: r.update(parameters=["E", "E"]), "unique"),
    (lambda r: r.update(parameters=["Young"]), "drawn from"),
    (lambda r: r["outputs"][0].update(field="LE"), "unsupported field"),
    (lambda r: r["outputs"][0].update(component=7), "component"),
    (lambda r: r["outputs"][0].update(reduction="median"), "reduction"),
    (lambda r: r["outputs"][0].update(field="U", reduction="volume_mean"), "integration-point"),
    (lambda r: r["outputs"][0].update(domain={"points": [9]}), "points"),
    (lambda r: r.update(increments=[0]), "increments"),
    (lambda r: r.update(weighted_shares={"field": "U"}), "weighted_shares"),
])
def test_request_validation(change, message):
    request = {"outputs": [{"name": "a", "field": "S", "component": 1, "reduction": "max"}],
               "parameters": ["E", "nu"], "domain": {"elements": "ALL"}, "increments": "ALL"}
    change(request)
    with pytest.raises(RequestError, match=message):
        validate_request(request, ["E", "nu", "SIGY0", "H"], 1)


def test_parameters_all_follow_the_provider_order():
    request = {"outputs": [], "parameters": "ALL", "domain": {}, "increments": "LAST"}
    assert validate_request(request, ["H", "E"], 1) == ["H", "E"]


def test_field_share_reduces_to_the_single_point_definition(provider_factory, tmp_path):
    """One integration point in the domain: share_j = |p_j dq/dp_j| / sum_k |p_k dq/dp_k|."""
    _, _, material = provider_factory("m3_j2")
    _, model = beam_model(tmp_path, n=(4, 2, 1), push=0.05, steps=3)
    engine = HistoryEngine(model, material)
    fields = Fields(engine, run_history(engine, times=np.linspace(0, 1, 4), rtol=1e-12))
    p = fields.result.parameter_values
    domain = {"elements": [1], "points": [1]}
    rows = weighted_shares(fields, p, "MISES", domain)
    for row in rows:
        _, derivative, _, _, _ = select(fields, "MISES", 1, domain, row["increment"])
        weights = np.abs(derivative[0] * p)
        np.testing.assert_allclose(row["field_share"], 100 * weights / weights.sum())
        np.testing.assert_allclose(row["scalar_share"], row["field_share"])
    assert rows[-1]["field_share"][2] > 0                 # element 1 (root) has yielded


def test_volume_weights_come_from_the_geometry(provider_factory, tmp_path):
    """A stretched mesh: volume_mean weighs by det(J) w_q, not by point count."""
    _, _, material = provider_factory("m3_j2")

    def stretch(text):
        lines = text.splitlines()
        start, end = lines.index("*Node") + 1, next(i for i, l in enumerate(lines) if l.startswith("*Element"))
        for i in range(start, end):
            node, x, y, z = [s.strip() for s in lines[i].split(",")]
            x = float(x)
            lines[i] = "%s, %r, %s, %s" % (node, x + 0.5 * x * x / 6.0, y, z)   # graded along x
        return "\n".join(lines) + "\n"

    _, model = beam_model(tmp_path, n=(6, 2, 1), push=0.03, steps=2, extra=stretch)
    engine = HistoryEngine(model, material)
    assert engine.w.std() > 0.1 * engine.w.mean()
    fields = Fields(engine, run_history(engine, times=np.linspace(0, 1, 3), rtol=1e-12))
    values, derivatives, weights, _, _ = select(fields, "MISES", 1, {"elements": "ALL"}, 2)
    value, gradient, _ = reduce(values, derivatives, weights, "volume_mean")
    assert value == pytest.approx(np.sum(values * weights) / np.sum(weights))
    assert value != pytest.approx(values.mean(), rel=1e-3)
    assert weights.sum() == pytest.approx(np.sum(engine.volume))
