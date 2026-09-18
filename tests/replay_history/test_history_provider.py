"""The batched provider link: equality with the single-point path, refusals."""
import json

import numpy as np
import pytest

from residual_core.replay.history_material import HistoryMaterial, ProviderError
from residual_core.replay.path_material import PathMaterial


def _points(material, count, seed=3):
    rng = np.random.default_rng(seed)
    stress = rng.normal(size=(count, 6)) * 80.0
    state = np.abs(rng.normal(size=(count, material.nstatev))) * 1e-3
    stran = rng.normal(size=(count, 6)) * 1e-3
    dstran = rng.normal(size=(count, 6)) * 1.5e-3
    return stress, state, stran, dstran


def _common(count, dtime=1.0):
    return dict(time=[0.0, 0.0], dtime=dtime, coords=np.zeros((count, 3)), celent=np.ones(count),
                noel=np.arange(1, count + 1), npt=np.ones(count))


def test_batched_total_equals_single_point_eval_with_history_carry(provider_factory, tmp_path):
    obj, contract, material = provider_factory("m3_j2")
    single = PathMaterial(str(obj), contract, str(tmp_path / "path"))
    props = np.array([200000.0, 0.3, 250.0, 2000.0])
    count = 40
    stress, state, stran, dstran = _points(material, count)
    rng = np.random.default_rng(5)
    dsig = rng.normal(size=(count, 6, 4))
    dstv = rng.normal(size=(count, 1, 4)) * 1e-6
    zero = np.zeros((count, 6, 4))
    out = material.total(props, stress, state, stran, dstran, dstress_in=dsig, dstate_in=dstv,
                         stran_dp=zero, dstran_dp=zero, **_common(count))
    regular = material.regular(props, stress, state, stran, dstran, **_common(count))
    plastic = 0
    for i in range(count):
        s, st, dd, ds, dv = single._oti_step(props, stress[i], state[i], dsig[i], dstv[i], dstran[i], 1.0)
        np.testing.assert_array_equal(out["stress"][i], s)
        np.testing.assert_array_equal(out["state"][i], st)
        np.testing.assert_array_equal(out["ddsdde"][i], dd)
        np.testing.assert_array_equal(out["dstress_dp"][i], ds)
        np.testing.assert_array_equal(out["dstate_dp"][i], dv)
        rs, rst, _ = single._reg_step(props, stress[i], state[i], dstran[i], 1.0)
        np.testing.assert_array_equal(regular["stress"][i], rs)
        plastic += int(st[0] > state[i, 0])
    assert 0 < plastic < count            # both branches of the return map are exercised


def test_strain_seed_is_the_directional_derivative_of_the_original_umat(provider_factory):
    """dSTRAN/dp and dDSTRAN/dp seeds vs central FD of the ORIGINAL UMAT (damage fixture reads STRAN)."""
    _, _, material = provider_factory("damage")
    props = np.array([70000.0, 0.25, 5e-4, 400.0, 0.05])
    count = 8
    rng = np.random.default_rng(11)
    stran = rng.normal(size=(count, 6)) * 6e-4
    dstran = rng.normal(size=(count, 6)) * 2e-4
    state = np.full((count, 1), 4e-4)
    stress = np.zeros((count, 6))
    values = material.parameter_values(props)
    # seeds per unit parameter: a relative strain change equal to the relative parameter change
    e0 = rng.normal(size=(count, 6, 4)) * 6e-4 / np.abs(values)
    de = rng.normal(size=(count, 6, 4)) * 2e-4 / np.abs(values)
    zero = np.zeros((count, 6, 4))
    out = material.total(props, stress, state, stran, dstran, dstress_in=zero,
                         dstate_in=np.zeros((count, 1, 4)), stran_dp=e0, dstran_dp=de, **_common(count))
    for j, name in enumerate(material.params):
        estimates = []
        for h in (1e-6, 5e-7):
            step = h * abs(values[j])
            plus = material.regular(material.props_with(props, {name: values[j] + step}), stress, state,
                                    stran + step * e0[..., j], dstran + step * de[..., j], **_common(count))
            minus = material.regular(material.props_with(props, {name: values[j] - step}), stress, state,
                                     stran - step * e0[..., j], dstran - step * de[..., j], **_common(count))
            estimates.append(((plus["stress"] - minus["stress"]) / (2 * step),
                              (plus["state"] - minus["state"]) / (2 * step)))
        for k, key in enumerate(("dstress_dp", "dstate_dp")):
            reference = estimates[1][k]
            scale = max(np.abs(reference).max(), 1e-30)
            assert np.abs(estimates[0][k] - reference).max() / scale < 1e-6
            assert np.abs(out[key][..., j] - reference).max() / scale < 1e-6, (name, key)


def test_cutback_request_is_refused_with_a_message(provider_factory):
    _, _, material = provider_factory("damage")
    props = np.array([70000.0, 0.25, 5e-4, 400.0, 0.05])
    dstran = np.array([[0.08, 0, 0, 0, 0, 0]])
    with pytest.raises(ProviderError, match="cut-back"):
        material.regular(props, np.zeros((1, 6)), np.zeros((1, 1)), np.zeros((1, 6)), dstran, **_common(1))


def test_contract_without_total_entry_point_is_refused(provider_factory, tmp_path):
    obj, contract, _ = provider_factory("m3_j2")
    legacy = json.loads(json.dumps(contract))
    del legacy["symbols"]["oti_eval_total"]
    with pytest.raises(ProviderError, match="UMAT_OTI_EVAL_TOTAL"):
        HistoryMaterial(obj, legacy, str(tmp_path / "legacy"))


def test_object_hash_is_bound_to_the_contract(provider_factory, tmp_path):
    obj, contract, _ = provider_factory("m3_j2")
    changed = json.loads(json.dumps(contract))
    changed["object"]["sha256_full"] = "0" * 64
    with pytest.raises(ProviderError, match="sha256"):
        HistoryMaterial(obj, changed, str(tmp_path / "hash"))


def test_parameter_order_and_slots_come_from_the_contract(provider_factory):
    _, contract, material = provider_factory("m6_fcc")
    assert material.params == [p["name"] for p in contract["parameters"]]
    assert material.params[:3] == ["g0", "h0", "q"]             # provider order, not PROPS order
    props = np.arange(1.0, 11.0)
    assert material.parameter_values(props).tolist() == [props[p["props_index"] - 1]
                                                         for p in contract["parameters"]]
