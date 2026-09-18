"""Slide 18/42's Solve screen, driven through AppTest (offline suite).

"Point at the compiled OTI object from the previous build. Point at the saved
analysis .odb and its setup .inp. Tick the parameters, then choose the output
and the region. Click Solve."

The object is built fresh from the UMAT repository's m3_j2 contract, and the
mapping beside it is the build's own generated contract, exactly as the
provider screen hands them over (OTI_UMAT.obj + Mapping.json). The .inp is the
presentation example. Reading an ODB needs licensed Abaqus Python, so here the
ODB export is replaced by the genuine export of that ODB kept in
tests/fixtures/presentation_j2/fields.json (Abaqus 2021.HF5, job imqrp_j2);
everything after the export -- conversion, replay, solve, reductions -- is the
real ``resasm request``. The browser variant reads the real ODB
(``test_imqcam_solve_screen_browser.py``, ``pytest -m gui``).

References: the one-element model is uniaxial, load-controlled and plastic at
the last increment, so U1 = 300/E + (300 - SIGY0)/H and EQPLAS = (300 - SIGY0)/H
with E = 210000, SIGY0 = 250, H = 2000 from the .inp, and the stress does not
depend on the material at all.
"""

from __future__ import annotations

import copy
import json
import os
import shutil
from pathlib import Path

import numpy as np
import pytest

from streamlit.testing.v1 import AppTest

ROOT = Path(__file__).resolve().parents[2]
EXAMPLE = ROOT / "examples" / "presentation_request"
FIELDS = ROOT / "tests" / "fixtures" / "presentation_j2" / "fields.json"
E, SIGY0, H = 210000.0, 250.0, 2000.0
SCREEN = ("from residual_core.app.streamlit_app import _init_state, _tab_request\n"
          "_init_state()\n_tab_request()")

pytestmark = [pytest.mark.integration, pytest.mark.fortran, pytest.mark.slow]


def _umat_repository() -> Path:
    configured = os.environ.get("UMAT_OTI_REPO")
    for candidate in ([Path(configured)] if configured else []) + [ROOT.parent / "UMAT_source_transformation"]:
        if (candidate / "parameter_sensitivity" / "models" / "m3_j2" / "contract_v2.json").is_file():
            return candidate
    pytest.fail("the UMAT repository (UMAT_OTI_REPO) with parameter_sensitivity/models/m3_j2 is needed "
                "to build the OTI object this screen consumes")


@pytest.fixture(scope="module")
def provider(tmp_path_factory):
    from umat_oti.provider.build import build_provider

    directory = tmp_path_factory.mktemp("provider")
    built = build_provider(_umat_repository() / "parameter_sensitivity/models/m3_j2/contract_v2.json",
                           directory / "build")
    shared = directory / "shared"
    shared.mkdir()
    shutil.copyfile(built["object"], shared / "OTI_UMAT.obj")
    shutil.copyfile(built["contract"], shared / "Mapping.json")
    return shared


@pytest.fixture
def inputs(provider, tmp_path, monkeypatch):
    from residual_core.replay import presentation

    for name in ("OTI_UMAT.obj", "Mapping.json"):
        shutil.copyfile(provider / name, tmp_path / name)
    shutil.copyfile(EXAMPLE / "Analysis.inp", tmp_path / "Analysis.inp")
    (tmp_path / "Analysis.odb").write_bytes(b"offline transport: the export below stands in for the ODB")
    fields = json.loads(FIELDS.read_text())

    def exported(odb, destination, abaqus):
        destination.write_text(json.dumps(fields))
        return copy.deepcopy(fields), ["offline-genuine-export"]

    monkeypatch.setattr(presentation, "export_odb", exported)
    return tmp_path


def _solve(inputs: Path, out: Path, *, untick=(), output="displacement U1", region="whole mesh",
           reduction="mean") -> AppTest:
    app = AppTest.from_string(SCREEN, default_timeout=180).run()
    for key, name in (("material", "OTI_UMAT.obj"), ("model", "Analysis.inp"), ("odb", "Analysis.odb")):
        app.text_input(key="request_path_" + key).set_value(str(inputs / name))
    app.text_input(key="request_output").set_value(str(out))
    app.run()
    assert not app.exception and not app.error, [e.value for e in app.error]
    # The tick boxes are the parameters the Mapping.json names, all ticked.
    ticks = {box.label: box.value for box in app.checkbox if box.key.startswith("request_param_")}
    assert ticks == {"E": True, "nu": True, "SIGY0": True, "H": True}
    for name in untick:
        app.checkbox(key=f"request_param_{name}").uncheck()
    app.selectbox(key="request_output_choice").set_value(output).run()
    region_box = next(box for box in app.selectbox if box.key.startswith("request_region_"))
    assert region in region_box.options, region_box.options
    region_box.set_value(region)
    app.selectbox(key="request_reduction").set_value(reduction).run()
    assert not app.button(key="btn_request_run").disabled
    app.button(key="btn_request_run").click().run(timeout=180)
    assert not app.exception and not app.error, [e.value for e in app.error]
    return app


def _cli(inputs: Path, request: dict, out: Path) -> dict:
    """``resasm request`` on the same four files."""
    from residual_core.ui.cli import main

    path = out.parent / (out.name + "_request.json")
    path.write_text(json.dumps(request))
    code = main(["request", "--model", str(inputs / "Analysis.inp"), "--odb", str(inputs / "Analysis.odb"),
                 "--material", str(inputs / "OTI_UMAT.obj"), "--request", str(path), "--out", str(out)])
    assert code == 0
    return json.loads((out / "sensitivity_results.json").read_text())


def test_tick_choose_solve_on_the_loaded_nodes(inputs, tmp_path):
    out = tmp_path / "gui"
    app = _solve(inputs, out, untick=("nu",), output="displacement U1", region="node set LOADED")
    summary = json.loads((out / "sensitivity_results.json").read_text())
    request = summary["request"]
    assert request["parameters"] == ["E", "SIGY0", "H"]
    assert request["outputs"][0]["domain"] == {"nodes": [2, 3, 6, 7]}
    (row,) = summary["results"]
    assert row["value"] == pytest.approx(300 / E + (300 - SIGY0) / H)
    assert row["derivatives"]["E"] == pytest.approx(-300 / E ** 2)
    assert row["derivatives"]["SIGY0"] == pytest.approx(-1 / H)
    assert row["derivatives"]["H"] == pytest.approx(-(300 - SIGY0) / H ** 2)
    assert "nu" not in row["derivatives"]
    # The same request through the command line gives the same numbers.
    assert _cli(inputs, request, tmp_path / "cli")["results"] == summary["results"]
    # "... see the results of the sensitivity in the defined domain": one row
    # per node of the region, whose mean is the reduced result.
    shown = app.dataframe[0].value
    assert sorted(shown["location"]) == ["node 2", "node 3", "node 6", "node 7"]
    for name in ("E", "SIGY0", "H"):
        assert np.mean(shown[f"d/d{name}"]) == pytest.approx(row["derivatives"][name], rel=1e-12)
    assert set(shown["governing parameter"]) == {"SIGY0"}
    successes = [s.value for s in app.success]
    assert "Executed: 1 scalar results. Independent validation: not run." in successes
    assert any(s.startswith("Solved — d(mean_U1_LOADED)/dp at all 4 locations") for s in successes)
    metrics = {m.label: m.value for m in app.metric}
    assert float(metrics["stress reproduced vs. the .odb"]) < 2e-5
    labels = [button.label for button in app.get("download_button")]
    assert labels == ["sensitivity_results.json", "sensitivity_tables.csv", "run_report.txt"]


@pytest.mark.parametrize("output,reference", [
    ("stress S11", {"value": 300.0, "E": 0.0, "SIGY0": 0.0, "H": 0.0}),
    ("state variable SDV1", {"value": (300 - SIGY0) / H, "E": 0.0, "SIGY0": -1 / H,
                             "H": -(300 - SIGY0) / H ** 2}),
])
def test_every_integration_point_of_the_whole_mesh(inputs, tmp_path, output, reference):
    app = _solve(inputs, tmp_path / "gui", untick=("nu",), output=output, region="whole mesh")
    shown = app.dataframe[0].value
    assert len(shown) == 8 and all(location.startswith("element 1 IP ") for location in shown["location"])
    np.testing.assert_allclose(shown["value"], reference["value"], rtol=2e-5)
    for name in ("E", "SIGY0", "H"):
        scale = abs(reference[name]) or 1e-9 * abs(reference["value"])
        np.testing.assert_allclose(shown[f"d/d{name}"], reference[name], rtol=1e-6, atol=scale)
    summary = json.loads((tmp_path / "gui" / "sensitivity_results.json").read_text())
    assert _cli(inputs, summary["request"], tmp_path / "cli")["results"] == summary["results"]


def test_a_supplied_request_still_wins(inputs, tmp_path):
    """Copilot's request-file route is unchanged: the file is used as given."""
    shutil.copyfile(EXAMPLE / "sensitivity_request.json", inputs / "sensitivity_request.json")
    app = AppTest.from_string(SCREEN, default_timeout=180).run()
    for key, name in (("material", "OTI_UMAT.obj"), ("model", "Analysis.inp"), ("odb", "Analysis.odb"),
                      ("request", "sensitivity_request.json")):
        app.text_input(key="request_path_" + key).set_value(str(inputs / name))
    app.text_input(key="request_output").set_value(str(tmp_path / "gui"))
    app.run()
    app.button(key="btn_request_run").click().run(timeout=180)
    assert not app.exception and not app.error
    summary = json.loads((tmp_path / "gui" / "sensitivity_results.json").read_text())
    assert summary["request"] == json.loads((EXAMPLE / "sensitivity_request.json").read_text())
    assert len(summary["results"]) == 4


def test_without_a_mapping_there_is_nothing_to_tick(inputs, tmp_path):
    (inputs / "Mapping.json").unlink()
    app = AppTest.from_string(SCREEN, default_timeout=180).run()
    app.text_input(key="request_path_material").set_value(str(inputs / "OTI_UMAT.obj")).run()
    assert not [box for box in app.checkbox if box.key.startswith("request_param_")]
    assert app.button(key="btn_request_run").disabled


def test_a_reaction_output_is_solved_and_says_what_is_not_stored(inputs, tmp_path):
    """Equilibrium: the support reaction is -300 whatever the material."""
    app = _solve(inputs, tmp_path / "gui", output="reaction force RF1", region="node set XZERO",
                 reduction="sum")
    summary = json.loads((tmp_path / "gui" / "sensitivity_results.json").read_text())
    (row,) = summary["results"]
    assert row["domain"] == {"nodes": [1, 4, 5, 8]}
    assert row["value"] == pytest.approx(-300.0, rel=2e-5)
    assert max(abs(value) for value in row["derivatives"].values()) < 1e-9
    assert any(s.value.startswith("Solved (increment 4). Reaction derivatives are not stored")
               for s in app.success)
    assert not app.dataframe
    assert _cli(inputs, summary["request"], tmp_path / "cli")["results"] == summary["results"]
