import json
from pathlib import Path

import numpy as np
import pytest

from residual_core.replay.path_material import PathMaterial
from repository_paths import umat_repo_root


ROOT = Path(__file__).resolve().parents[2]
PROVIDER = umat_repo_root()
CONTRACT = PROVIDER / "parameter_sensitivity/models/m3_j2/contract_v2.json"
PROPS = np.array([210000., 0.3, 250., 2000.])
PATH = np.array([[0.0002, 0, 0, 0, 0, 0],
                 [0.0004, 0, 0, 0.0001, 0, 0],
                 [0.002, 0, 0, 0.0004, 0, 0],
                 [0.001, 0, 0, -0.0002, 0, 0],
                 [-0.0003, 0, 0, 0, 0, 0],
                 [-0.004, 0, 0, -0.0005, 0, 0]])


@pytest.fixture(scope="module")
def material(tmp_path_factory):
    from umat_oti.provider.build import build_provider
    output = tmp_path_factory.mktemp("connected-j2")
    built = build_provider(CONTRACT, output)
    metadata = json.loads(Path(built["contract"]).read_text())
    return PathMaterial(built["object"], metadata, str(output / "consumer"))


def test_compiled_eval_march_original_whole_path(material):
    times = np.ones(len(PATH))
    actual = material.march_oti(PROPS, PATH, times)
    fast = material.march_fast(PROPS, PATH, times)
    original = material.march_regular(PROPS, PATH, times)
    assert material.has_march
    assert material.nstatev == 1
    for step, reference, march in zip(actual, original, fast):
        np.testing.assert_allclose(step["stress"], reference["stress"], atol=1e-10)
        np.testing.assert_allclose(step["statev"], reference["statev"], atol=1e-14)
        np.testing.assert_allclose(step["dsigma_dp"], march["dsigma_dp"], atol=1e-10)
    for column in range(4):
        for relative_step in (1e-4, 3e-5, 1e-5):
            step = relative_step * max(abs(PROPS[column]), 1)
            plus, minus = PROPS.copy(), PROPS.copy()
            plus[column] += step
            minus[column] -= step
            upper = material.march_regular(plus, PATH, times)
            lower = material.march_regular(minus, PATH, times)
            for value, derivative in (("stress", "dsigma_dp"), ("statev", "dstatev_dp")):
                finite_difference = np.array([(hi[value] - lo[value]) / (2 * step)
                                              for hi, lo in zip(upper, lower)])
                lifted = np.array([row[derivative][:, column] for row in actual])
                scale = max(np.max(np.abs(finite_difference)), 1e-12)
                assert np.max(np.abs(lifted - finite_difference)) / scale < 2e-6


def test_strain_history_chain_includes_state(material):
    from residual_core.replay.j2_history import step_j2
    path_derivative = np.zeros((len(PATH), 6, 4))
    path_derivative[:, :, 1] = PATH * 0.37
    stress, state = np.zeros(6), np.zeros(1)
    dstress, dstate = np.zeros((6, 4)), np.zeros((1, 4))
    for increment, derivative in zip(PATH, path_derivative):
        stress, state, tangent, dstress, dstate = step_j2(
            material, PROPS, stress, state, dstress, dstate, increment, derivative, 1.)
    for step in (1e-4, 3e-5, 1e-5):
        upper, lower = PROPS.copy(), PROPS.copy()
        upper[1] += step
        lower[1] -= step
        plus = material.march_regular(upper, PATH + step * path_derivative[:, :, 1], np.ones(len(PATH)))[-1]
        minus = material.march_regular(lower, PATH - step * path_derivative[:, :, 1], np.ones(len(PATH)))[-1]
        for value, derivative in (("stress", dstress), ("statev", dstate)):
            reference = (plus[value] - minus[value]) / (2 * step)
            assert np.max(np.abs(derivative[:, 1] - reference)) / np.max(np.abs(reference)) < 2e-6


@pytest.mark.parametrize("method", ["march_oti", "march_fast", "march_regular"])
def test_path_validation(material, method):
    with pytest.raises(ValueError, match="times"):
        getattr(material, method)(PROPS, PATH, [1.])
    with pytest.raises(ValueError, match="positive"):
        getattr(material, method)(PROPS, PATH, np.zeros(len(PATH)))


def test_no_live_oti_cast(material):
    from residual_core.algebra.otilib_adapter import OtiContext
    from residual_core.replay.j2_history import real_array
    seeded = OtiContext(1, 1).seed(2., 1)
    with pytest.raises(ValueError, match="live OTI"):
        real_array([seeded], (1,), "input")


def test_equilibrated_history_against_original_whole_model(material):
    import copy
    from residual_core.replay.connected import solve_history
    from residual_core.replay.record import ReplayRecord
    record = ReplayRecord.load(str(ROOT / "examples/bounded_j2_c3d8/model.json"))
    result = solve_history(record, material)
    repeated = solve_history(ReplayRecord(result["record"]), material, replay=True)
    actual = np.array([row["du_dp"] for row in result["increments"]])
    np.testing.assert_allclose(actual, [row["du_dp"] for row in repeated["increments"]], atol=1e-13)
    for column in range(4):
        previous_fd = None
        for relative in (1e-4, 3e-5, 1e-5):
            step = relative * max(abs(PROPS[column]), 1.)
            plus, minus = copy.deepcopy(record.raw), copy.deepcopy(record.raw)
            plus["material"]["props"][column] += step
            minus["material"]["props"][column] -= step
            upper = solve_history(ReplayRecord(plus), material, original=True)
            lower = solve_history(ReplayRecord(minus), material, original=True)
            reference = (np.array([row["u"] for row in upper["increments"]]) -
                         np.array([row["u"] for row in lower["increments"]])) / (2 * step)
            scale = max(np.max(np.abs(reference)), 1e-12)
            assert np.max(np.abs(actual[:, :, column] - reference)) / scale < 2e-6
            if relative == 1e-5:
                assert np.max(np.abs(reference - previous_fd)) / scale < 2e-6
            previous_fd = reference
    naive = np.array([row["fixed_path_local_solve"] for row in result["increments"]])
    assert np.max(np.abs(naive - actual)) > 1e-5
    assert max(row["equilibrium_error"] for row in result["increments"]) < 1e-11
    assert np.max(result["increments"][-1]["state"]) > 0


def test_verification_including_state_and_real_abaqus(material):
    from residual_core.replay.connected import solve_history
    from residual_core.replay.record import ReplayRecord
    from residual_core.replay.verification import verify_connected, verify_abaqus_fixture
    record = ReplayRecord.load(str(ROOT / "examples/bounded_j2_c3d8/model.json"))
    report = verify_connected(record, material, solve_history(record, material))
    assert report["passed"] and len(report["rows"]) == 12
    assert verify_abaqus_fixture(ROOT / "tests/abaqus_derivative_export")["passed"]


def test_public_cli_solve_replay_and_failure(material, tmp_path):
    from residual_core.ui.cli import main
    object_path = Path(material.work).parent / material.contract["object"]["file"]
    contract_path = object_path.with_suffix(".json")
    common = ["--object", str(object_path), "--contract", str(contract_path)]
    output = tmp_path / "solved"
    assert main(["replay", str(ROOT / "examples/bounded_j2_c3d8/model.json"), *common,
                 "--out", str(output), "--solve", "--verify"]) == 0
    assert json.loads((output / "public/summary.json").read_text())["verified"]
    assert main(["replay", str(output / "private/record.json"), *common,
                 "--out", str(tmp_path / "replayed")]) == 0
    assert main(["replay", str(tmp_path / "missing.json"), *common,
                 "--out", str(tmp_path / "failure")]) == 2
    assert json.loads((tmp_path / "failure/public/summary.json").read_text())["status"] == "failed"


def test_gui_replay_uses_real_cli(material, tmp_path):
    from streamlit.testing.v1 import AppTest
    object_path = Path(material.work).parent / material.contract["object"]["file"]
    app = AppTest.from_string("from residual_core.app.streamlit_app import _init_state, _tab_replay\n_init_state()\n_tab_replay()")
    app.run()
    app.text_input(key="replay_object").set_value(str(object_path))
    app.text_input(key="replay_contract").set_value(str(object_path.with_suffix(".json")))
    app.text_input(key="replay_output").set_value(str(tmp_path / "gui"))
    app.run()
    app.button(key="btn_replay_run").click().run(timeout=60)
    assert not app.exception
    summary = json.loads((tmp_path / "gui/public/summary.json").read_text())
    assert summary["status"] == "completed" and summary["verified"]


@pytest.mark.parametrize("mutation, message", [
    ("fingerprint", "regular_source_hash"), ("pressure", "concentrated"),
    ("initial_state", "unsupported record"), ("boundary", "homogeneous"),
    ("element", "C3D8"), ("time", "positive"), ("material", "one homogeneous")])
def test_unsupported_physics_is_not_silently_ignored(material, mutation, message):
    from residual_core.replay.connected import solve_history
    from residual_core.replay.record import ReplayRecord
    record = ReplayRecord.load(str(ROOT / "examples/bounded_j2_c3d8/model.json"))
    if mutation == "fingerprint":
        record.raw["provenance"]["regular_source_hash"] = "wrong"
    elif mutation == "pressure":
        record.raw["loads"]["pressure"] = 10
    elif mutation == "initial_state":
        record.raw["initial_state"] = [1]
    elif mutation == "boundary":
        record.raw["boundaries"][0]["value"] = .1
    elif mutation == "element":
        record.raw["mesh"]["elements"]["1"]["type"] = "C3D8R"
    elif mutation == "time":
        record.raw["increments"][0]["dt"] = 0
    elif mutation == "material":
        record.raw["material"]["temperature"] = 300
    with pytest.raises(ValueError, match=message):
        solve_history(record, material)


def test_recorded_physical_state_is_checked(material):
    from residual_core.replay.connected import solve_history
    from residual_core.replay.record import ReplayRecord
    record = ReplayRecord.load(str(ROOT / "examples/bounded_j2_c3d8/model.json"))
    result = solve_history(record, material)
    result["record"]["increments"][0]["state_ip"]["1"][0][0] += .01
    with pytest.raises(AssertionError):
        solve_history(ReplayRecord(result["record"]), material, replay=True)


def test_field_inputs_reject_live_oti():
    from residual_core.algebra.otilib_adapter import OtiContext
    from residual_core.core.field_sensitivity import _as_tangent_ip, _as_dsigma_ip, fields_from_statev
    seeded = OtiContext(1, 1).seed(2., 1)
    for shape, validate in (((8, 6, 6), lambda values: _as_tangent_ip(values, 1)),
                            ((8, 6), lambda values: _as_dsigma_ip(values, 1, "E")),
                            ((8, 42), lambda values: fields_from_statev({1: values}, {"ddsdde": [1, 36], "parameters": {"E": [37, 42]}}))):
        values = np.full(shape, seeded, dtype=object)
        with pytest.raises(ValueError, match="live OTI"):
            validate(values)


def test_reproducer_fresh_build_and_nonzero_failure(tmp_path):
    import subprocess
    import sys
    script = ROOT / "scripts/reproduce_connected_pipeline.py"
    output = tmp_path / "reproduction"
    command = [sys.executable, str(script), "--skip-abaqus", "--provider-repo", str(PROVIDER), "--out", str(output)]
    completed = subprocess.run(command, cwd=tmp_path, capture_output=True, text=True)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    manifest = json.loads((output / "private/manifest.json").read_text())
    assert manifest["passed"]
    assert manifest["imports"] == "installed"
    assert all("site-packages" in path for path in manifest["backend_modules"].values())
    assert manifest["environment"]["PYTHONPATH"] is None
    assert all(len(repo["commit"]) == 40 for repo in manifest["repositories"].values())
    assert len(manifest["provider"]["object_sha256"]) == 64
    assert manifest["abaqus_fixture"]["passed"]
    assert (output / "public/j2_history.png").stat().st_size > 1000
    assert (output / "public/fd_table.csv").is_file()
    failed = subprocess.run([sys.executable, str(script), "--skip-abaqus", "--model", str(tmp_path / "absent.json"),
                             "--out", str(tmp_path / "failure")], cwd=ROOT, capture_output=True, text=True)
    assert failed.returncode == 2
    assert not json.loads((tmp_path / "failure/manifest.json").read_text())["passed"]