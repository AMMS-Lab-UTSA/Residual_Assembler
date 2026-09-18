import numpy as np
import pytest
import copy
import json
from pathlib import Path
import shutil

from residual_core.replay.request_reductions import reduce_scalar


@pytest.mark.parametrize("reduction,value,derivative", [
    ("sum", 7., [4., 6.]), ("mean", 3.5, [2., 3.]),
    ("L2", 5., [3., 4.4]), ("max", 4., [3., 4.]),
])
def test_scalar_reductions(reduction, value, derivative):
    actual, gradient = reduce_scalar([3., 4.], [[1., 2.], [3., 4.]], reduction)
    assert actual == pytest.approx(value)
    np.testing.assert_allclose(gradient, derivative)


@pytest.mark.parametrize("values,reduction,message", [
    ([0., 0.], "L2", "nondifferentiable"),
    ([2., 2.], "max", "tied"),
    ([1., 2.], "component", "exactly one"),
    ([1., 2.], "median", "unsupported"),
])
def test_invalid_reductions(values, reduction, message):
    with pytest.raises(ValueError, match=message):
        reduce_scalar(values, np.ones((2, 2)), reduction)


ROOT = Path(__file__).resolve().parents[2]
EXAMPLE = ROOT / "examples/presentation_request"


@pytest.mark.unit
@pytest.mark.parametrize("stage,category,action", [
    ("export_odb", "odb_export", "Check --odb and --abaqus; licensed Abaqus Python with odbAccess and a complete field history are required."),
    ("solve_history", "history_replay", "Check the saved history, material mapping and equilibrium consistency in the private diagnostics."),
])
def test_unit_failure_diagnostics_are_private(tmp_path, monkeypatch, capsys, stage, category, action):
    """Injected failures test disclosure boundaries, not numerical correctness."""
    from residual_core.replay import verification
    from residual_core.ui.cli import main
    from streamlit.testing.v1 import AppTest
    from subprocess import CompletedProcess

    marker = "SENSITIVE_FIELD_ARRAY_918 [123456.789, 987654.321] material source / exporter stdout"
    export_odb = presentation.export_odb
    inputs = {"model": tmp_path / "Analysis.inp", "odb": tmp_path / "Analysis.odb",
              "material": tmp_path / "OTI_UMAT.obj", "request": tmp_path / "sensitivity_request.json"}
    for path in inputs.values():
        path.write_text("{}")
    monkeypatch.setattr(presentation, "mapping_for", lambda *args: (None, {"parameters": []}))
    monkeypatch.setattr(presentation, "validate_request", lambda *args: None)
    monkeypatch.setattr(presentation, "read_model", lambda *args: None)
    monkeypatch.setattr(presentation, "PathMaterial", lambda *args: None)
    monkeypatch.setattr(presentation, "require_j2", lambda *args: None)
    monkeypatch.setattr(presentation, "export_odb", lambda *args: ({}, []))
    monkeypatch.setattr(presentation, "convert_fields", lambda *args: None)

    def fail(*args, **kwargs):
        raise AssertionError(marker)

    if stage == "export_odb":
        monkeypatch.setattr(presentation, "export_odb", export_odb)
        monkeypatch.setattr(presentation.shutil, "which", lambda *args: "/test/abaqus")
        monkeypatch.setattr(presentation.subprocess, "run", lambda *args, **kwargs:
                            CompletedProcess(args[0], 7, marker, "stderr: " + marker))
    else:
        monkeypatch.setattr(presentation, stage, fail)
    expected = "Category: %s\nAction: %s\nPrivate diagnostics: private/error_report.txt" % (category, action)
    direct_out = tmp_path / "direct"
    with pytest.raises(presentation.RequestFailure) as caught:
        presentation.run_request(**inputs, out=direct_out)
    assert str(caught.value) == expected
    assert marker in str(caught.value.__cause__)

    cli_out = tmp_path / "cli"
    command = ["request"] + [part for key, value in {**inputs, "out": cli_out}.items()
                             for part in ("--" + key, str(value))]
    assert main(command) == 2
    captured = capsys.readouterr()
    assert captured.err == "request failed: " + expected + "\n"
    assert marker not in captured.out + captured.err

    gui_out = tmp_path / "gui"
    app = AppTest.from_string("from residual_core.app.streamlit_app import _init_state, _tab_request\n_init_state()\n_tab_request()")
    app.run()
    for key, path in inputs.items():
        app.text_input(key="request_path_" + key).set_value(str(path))
    app.text_input(key="request_output").set_value(str(gui_out))
    app.run()
    app.button(key="btn_request_run").click().run(timeout=60)
    assert not app.exception
    assert marker not in str(app)
    assert app.error
    visible = "\n".join(str(element.value) for kind in ("error", "code", "text", "markdown")
                        for element in app.get(kind))
    assert marker not in visible
    assert category in visible and "private/error_report.txt" in visible
    for output in (direct_out, cli_out, gui_out):
        report = (output / "run_report.txt").read_text()
        assert expected in report and "Status: failed" in report
        for name in presentation.PUBLIC_FILES:
            path = output / name
            assert not path.exists() or marker not in path.read_text()
        private = (output / "private/error_report.txt").read_text()
        assert marker in private and "Traceback (most recent call last)" in private
        if stage == "export_odb":
            assert marker in (output / "private/odb_export.log").read_text()


@pytest.mark.unit
@pytest.mark.parametrize("role,suffix", [("model", ".inp"), ("odb", ".odb"),
                                         ("material", ".obj"), ("request", ".json")])
def test_unit_missing_input_diagnostic_names_role_not_private_path(tmp_path, role, suffix):
    from residual_core.replay.presentation import RequestFailure, run_request
    marker = "PRIVATE_PATH_918"
    inputs = {"model": tmp_path / "Analysis.inp", "odb": tmp_path / "Analysis.odb",
              "material": tmp_path / "OTI_UMAT.obj", "request": tmp_path / "request.json"}
    for path in inputs.values():
        path.touch()
    inputs[role] = tmp_path / (marker + suffix)
    output = tmp_path / "output"
    with pytest.raises(RequestFailure) as caught:
        run_request(**inputs, out=output)
    assert caught.value.category == "input_path"
    assert caught.value.action == "Check --%s: provide an existing readable %s file." % (role, suffix)
    assert marker not in str(caught.value)
    assert marker in (output / "private/error_report.txt").read_text()


@pytest.fixture(scope="module")
def compiled(tmp_path_factory):
    from umat_oti.provider.build import build_provider
    from residual_core.replay.connected import solve_history
    from residual_core.replay.path_material import PathMaterial
    from residual_core.replay.record import ReplayRecord
    from residual_core.replay.presentation_inputs import read_model, target_nodes
    directory = tmp_path_factory.mktemp("presentation_provider")
    from repository_paths import umat_repo_root
    built = build_provider(umat_repo_root() / "parameter_sensitivity/models/m3_j2/contract_v2.json", directory)
    contract = json.loads(Path(built["contract"]).read_text())
    material = PathMaterial(built["object"], contract, str(directory / "link"))
    model = read_model(EXAMPLE / "Analysis.inp")
    raw = {"schema": "resasm_replay_record_v1", "kinematics": "small_strain", "integration": "selective_reduced",
           "mesh": {"nodes": {str(node): list(coords) for node, coords in model.nodes.items()},
                    "elements": {"1": {"type": "C3D8", "connectivity": list(range(1, 9))}}},
           "material": {"props": [210000., .3, 250., 2000.]},
           "boundaries": [{"target": node, "dof": boundary.dof_start, "value": 0.}
                          for boundary in model.boundaries for node in target_nodes(model, boundary.target)],
           "loads": {"cload": [{"node": node, "dof": 1, "value": 75.} for node in (2, 3, 6, 7)]},
           "increments": [{"dt": .25, "load_factor": number / 4.} for number in range(1, 5)],
           "provenance": {"regular_source_hash": contract["regular_source_hash"]}}
    result = solve_history(ReplayRecord(raw), material)
    fields = {"export_mode": "strict_request", "step_names": ["Loading"], "instance_names": ["PART-1-1"],
              "element_type": "C3D8", "nodes": raw["mesh"]["nodes"], "elements": {"1": list(range(1, 9))},
              "sdv_labels": [1], "frames": []}
    for number in range(5):
        row = result["increments"][number - 1] if number else {"u": np.zeros(24), "R": np.zeros(24),
                                                               "stress": np.zeros((1, 8, 6)), "state": np.zeros((1, 8, 1))}
        fields["frames"].append({"increment": number, "step": "Loading", "time": number / 4.,
            "U": {str(node): np.asarray(row["u"])[3 * (node - 1):3 * node].tolist() for node in range(1, 9)},
            "RF": {str(node): np.asarray(row["R"])[3 * (node - 1):3 * node].tolist() for node in range(1, 9)},
            "CF": {str(node): [75. * number / 4., 0., 0.] for node in (2, 3, 6, 7)},
            "S": {"1": np.asarray(row["stress"])[0].tolist()}, "SDV": {"1": np.asarray(row["state"])[0].tolist()}})
    return built, contract, result, fields


@pytest.fixture
def offline_inputs(compiled, tmp_path, monkeypatch):
    from residual_core.replay import verification
    built, contract, result, fields = compiled
    shutil.copyfile(built["object"], tmp_path / "OTI_UMAT.obj")
    (tmp_path / "Mapping.json").write_text(json.dumps(contract))
    shutil.copyfile(EXAMPLE / "Analysis.inp", tmp_path / "Analysis.inp")
    shutil.copyfile(EXAMPLE / "sensitivity_request.json", tmp_path / "sensitivity_request.json")
    (tmp_path / "Analysis.odb").write_bytes(b"offline transport placeholder, NOT an ODB")
    def exported(odb, destination, abaqus):
        destination.write_text(json.dumps(fields))
        return copy.deepcopy(fields), ["offline-test-transport"]
    monkeypatch.setattr(presentation, "export_odb", exported)
    return {"model": tmp_path / "Analysis.inp", "odb": tmp_path / "Analysis.odb",
            "material": tmp_path / "OTI_UMAT.obj", "request": tmp_path / "sensitivity_request.json", "out": tmp_path / "output"}


def test_compiled_offline_conversion_and_all_outputs(offline_inputs):
    from residual_core.replay.presentation import PUBLIC_FILES, run_request
    report = run_request(**offline_inputs)
    assert report["scope"]["history_increments_replayed"] == 4
    assert report["scope"]["output_increments"] == [4]
    assert not report["metadata"]["verified"]
    assert report["metadata"]["verification"]["status"] == "not_run"
    assert {path.name for path in offline_inputs["out"].iterdir()} == {*PUBLIC_FILES, "private"}
    displacement, reaction, stress, state = report["results"]
    assert displacement["value"] == pytest.approx(300 / 210000 + (300 - 250) / 2000)
    assert displacement["derivatives"]["E"] == pytest.approx(-300 / 210000 ** 2)
    assert displacement["derivatives"]["SIGY0"] == pytest.approx(-1 / 2000)
    assert state["derivatives"]["H"] == pytest.approx(-50 / 2000 ** 2)
    assert reaction["value"] == pytest.approx(-300)
    assert max(abs(value) for value in reaction["derivatives"].values()) < 1e-9
    assert stress["value"] == pytest.approx(300)
    import csv
    rows = list(csv.DictReader((offline_inputs["out"] / "sensitivity_tables.csv").open()))
    assert len(rows) == 16
    assert float(rows[0]["derivative"]) == displacement["derivatives"]["E"]
    assert "NOT RUN" in (offline_inputs["out"] / "run_report.txt").read_text()
    with pytest.raises(ValueError, match="stale"):
        run_request(**offline_inputs)


def test_request_cli_and_malformed_input(offline_inputs):
    from residual_core.ui.cli import main
    command = ["request"] + [part for key, value in offline_inputs.items() for part in ("--" + key, str(value))]
    assert main(command) == 0
    offline_inputs["request"].write_text('{"wrong": true}')
    command[-1] = str(offline_inputs["out"].parent / "failed")
    assert main(command) == 2
    assert "Status: failed" in (Path(command[-1]) / "run_report.txt").read_text()


@pytest.mark.parametrize("field", ["U", "RF", "CF", "S", "SDV"])
def test_missing_odb_field_has_safe_action(offline_inputs, compiled, monkeypatch, field):
    from residual_core.replay import verification
    fields = copy.deepcopy(compiled[3])
    fields["frames"][1].pop(field)
    monkeypatch.setattr(presentation, "export_odb", lambda *args: (fields, ["offline-test-transport"]))
    with pytest.raises(presentation.RequestFailure) as caught:
        presentation.run_request(**offline_inputs)
    assert caught.value.category == "odb_missing_field"
    assert caught.value.action == "Request complete ODB %s output for every required entity and saved increment." % ("SDV1" if field == "SDV" else field)
    assert "Traceback" in (offline_inputs["out"] / "private/error_report.txt").read_text()


def test_mapping_discovery_hash_and_alias(offline_inputs):
    from residual_core.replay.presentation import mapping_for
    material = offline_inputs["material"]
    path, contract = mapping_for(material)
    assert path.name == "Mapping.json"
    (material.with_suffix(".json")).write_text(json.dumps(contract))
    assert mapping_for(material)[1] == contract
    changed = copy.deepcopy(contract)
    changed["regular_source_hash"] = "wrong"
    material.with_suffix(".json").write_text(json.dumps(changed))
    with pytest.raises(ValueError, match="ambiguous"):
        mapping_for(material)
    assert mapping_for(material, path)[1] == contract
    material.write_bytes(material.read_bytes() + b"corruption")
    with pytest.raises(ValueError, match="sha256_full"):
        mapping_for(material, path)


@pytest.mark.parametrize("selection,expected", [("ALL", [1, 2, 3, 4]), ("LAST", [4]), ([2, 4], [2, 4])])
def test_selection_after_total_history(compiled, selection, expected):
    from residual_core.replay.presentation import scalar_results
    request = json.loads((EXAMPLE / "sensitivity_request.json").read_text())
    request["increments"] = selection
    request["parameters"] = ["SIGY0", "H"]
    rows = scalar_results(compiled[2], request)
    assert sorted(set(row["increment"] for row in rows)) == expected
    assert set(rows[0]["derivatives"]) == {"SIGY0", "H"}
    assert rows[-1]["derivatives"]["SIGY0"] == pytest.approx(-1 / 2000)


@pytest.mark.parametrize("mutation,message", [
    ("node", "node ids"), ("connectivity", "connectivity"), ("load", "loads"),
    ("history", "every increment"), ("time", "time period"), ("initial", "initial"), ("state", "SDV"),
])
def test_conversion_rejects_mismatches(compiled, mutation, message):
    from residual_core.replay.presentation_inputs import convert_fields, read_model
    fields = copy.deepcopy(compiled[3])
    if mutation == "node":
        fields["nodes"].pop("1")
    elif mutation == "connectivity":
        fields["elements"]["1"].reverse()
    elif mutation == "load":
        fields["frames"][1]["CF"]["2"][0] += 1
    elif mutation == "history":
        fields["frames"].pop(2)
    elif mutation == "time":
        fields["frames"].pop()
    elif mutation == "initial":
        fields["frames"][0]["SDV"]["1"][0][0] = .01
    elif mutation == "state":
        fields["frames"][1].pop("SDV")
    with pytest.raises((ValueError, AssertionError), match=message):
        convert_fields(read_model(EXAMPLE / "Analysis.inp"), fields, compiled[1])


@pytest.mark.parametrize("keyword", ["*Dload\n1, P1, 3.", "*Amplitude, name=A\n0.,0.,1.,1.", "*Boundary, op=NEW\n1,1,1,0."])
def test_no_dropped_inp_physics(tmp_path, keyword):
    from residual_core.replay.presentation_inputs import read_model
    path = tmp_path / "bad.inp"
    path.write_text((EXAMPLE / "Analysis.inp").read_text() + keyword + "\n")
    with pytest.raises(ValueError, match="unsupported"):
        read_model(path)


def test_missing_abaqus_reports_exact_binary(tmp_path):
    from residual_core.replay.presentation import export_odb
    with pytest.raises(ValueError, match="imqrp_missing_abaqus.*not found"):
        export_odb(tmp_path / "Analysis.odb", tmp_path / "fields.json", "imqrp_missing_abaqus")


def test_gui_shared_service_with_compiled_core(offline_inputs):
    from streamlit.testing.v1 import AppTest
    app = AppTest.from_string("from residual_core.app.streamlit_app import _init_state, _tab_request\n_init_state()\n_tab_request()")
    app.run()
    for key in ("model", "odb", "material", "request"):
        app.text_input(key="request_path_" + key).set_value(str(offline_inputs[key]))
    app.text_input(key="request_output").set_value(str(offline_inputs["out"]))
    assert not app.checkbox(key="request_validate").value
    app.run()
    app.button(key="btn_request_run").click().run(timeout=60)
    assert not app.exception and not app.error
    assert len(app.get("download_button")) == 3
    report = json.loads((offline_inputs["out"] / "sensitivity_results.json").read_text())
    assert report["scope"]["history_increments_replayed"] == 4


def test_odb_boundary_roundoff_is_bounded(compiled, tmp_path):
    from residual_core.replay.connected import solve_history
    from residual_core.replay.path_material import PathMaterial
    from residual_core.replay.record import ReplayRecord
    material = PathMaterial(compiled[0]["object"], compiled[1], str(tmp_path / "link"))
    raw = copy.deepcopy(compiled[2]["record"])
    raw["increments"][0]["u"][0] = 2e-35
    result = solve_history(ReplayRecord(raw), material, replay=True, odb_tolerances=True)
    assert result["increments"][0]["u"][0] == 0.
    with pytest.raises(ValueError, match="boundaries"):
        solve_history(ReplayRecord(raw), material, replay=True)
    raw["increments"][0]["u"][0] = 1e-8
    with pytest.raises(ValueError, match="boundaries"):
        solve_history(ReplayRecord(raw), material, replay=True, odb_tolerances=True)


def test_odb_float32_displacement_stress_roundoff(compiled, tmp_path):
    """eps32*E*max(strain) is ~6e-4 for this plastic path, including zero stress components."""
    from residual_core.replay.connected import solve_history
    from residual_core.replay.path_material import PathMaterial
    from residual_core.replay.record import ReplayRecord
    material = PathMaterial(compiled[0]["object"], compiled[1], str(tmp_path / "link"))
    raw = copy.deepcopy(compiled[2]["record"])
    for increment in raw["increments"]:
        increment["u"] = np.asarray(increment["u"], dtype=np.float32).astype(float).tolist()
    solve_history(ReplayRecord(raw), material, replay=True, odb_tolerances=True)
    raw["increments"][-1]["stress_ip"]["1"][0][1] += 1.
    with pytest.raises(AssertionError):
        solve_history(ReplayRecord(raw), material, replay=True, odb_tolerances=True)


@pytest.mark.parametrize("key,value", [("increments", [5]), ("increments", "FINAL"),
                                      ("parameters", ["temperature"]), ("domain", {"nodes": [999]})])
def test_invalid_request_scope(compiled, key, value):
    from residual_core.replay.presentation import scalar_results
    request = json.loads((EXAMPLE / "sensitivity_request.json").read_text())
    request[key] = value
    with pytest.raises(ValueError):
        scalar_results(compiled[2], request)


@pytest.mark.parametrize("field,component,reduction", [("MISES", 1, "mean"), ("S", 7, "mean"),
                                                       ("U", 1, "median"), ("U", 1, "component")])
def test_invalid_output_selection(compiled, field, component, reduction):
    from residual_core.replay.presentation import scalar_results
    request = json.loads((EXAMPLE / "sensitivity_request.json").read_text())
    request["outputs"] = [{"name": "bad", "field": field, "component": component, "reduction": reduction}]
    with pytest.raises(ValueError):
        scalar_results(compiled[2], request)


def test_mapping_missing_and_derivative_layout(offline_inputs):
    from residual_core.replay.presentation import mapping_for
    path = offline_inputs["material"].parent / "Mapping.json"
    contract = json.loads(path.read_text())
    path.rename(path.with_suffix(".saved"))
    with pytest.raises(ValueError, match="missing completed mapping"):
        mapping_for(offline_inputs["material"])
    contract["layouts"]["voigt"].reverse()
    path.write_text(json.dumps(contract))
    with pytest.raises(ValueError, match="Voigt"):
        mapping_for(offline_inputs["material"])


@pytest.mark.parametrize("line", ["1, 1.5, 1, 0.", "1, 1, 1, 0., 3", "1, XASYMM"])
def test_malformed_boundaries_not_coerced(tmp_path, line):
    from residual_core.replay.presentation_inputs import read_model
    path = tmp_path / "bad.inp"
    path.write_text((EXAMPLE / "Analysis.inp").read_text().replace("XZERO, 1, 1, 0.", line))
    with pytest.raises(ValueError):
        read_model(path)


def test_real_archived_j2_export(compiled, tmp_path):
    from residual_core.replay.connected import solve_history
    from residual_core.replay.path_material import PathMaterial
    from residual_core.replay.presentation_inputs import read_model, convert_fields
    from residual_core.replay.presentation import scalar_results
    fixture = ROOT / "tests/fixtures/presentation_j2/fields.json"
    fields = json.loads(fixture.read_text())
    material = PathMaterial(compiled[0]["object"], compiled[1], str(tmp_path / "link"))
    record = convert_fields(read_model(EXAMPLE / "Analysis.inp"), fields, compiled[1])
    result = solve_history(record, material, replay=True, odb_tolerances=True)
    assert max(row["equilibrium_error"] for row in result["increments"]) < 1e-5
    assert np.asarray(result["increments"][-1]["dsigma_dp"]).shape == (1, 8, 6, 4)
    assert np.asarray(result["increments"][-1]["dstate_dp"]).shape == (1, 8, 1, 4)
    request = json.loads((EXAMPLE / "sensitivity_request.json").read_text())
    rows = scalar_results(result, request)
    assert rows[0]["derivatives"]["E"] == pytest.approx(-300 / 210000 ** 2, rel=2e-5)
    assert abs(rows[0]["derivatives"]["nu"]) < 1e-8
    assert rows[0]["derivatives"]["SIGY0"] == pytest.approx(-1 / 2000, rel=2e-5)
    assert rows[0]["derivatives"]["H"] == pytest.approx(-50 / 2000 ** 2, rel=2e-5)
    assert rows[1]["value"] == pytest.approx(-300, abs=1e-3)
    assert max(abs(value) for value in rows[1]["derivatives"].values()) < 1e-9