"""Slide 18/42's Solve screen in a real browser, on the saved Abaqus analysis.

``pytest -m gui tests/gui`` (also marked ``abaqus``: licensed Abaqus Python
reads the ODB). The OTI object and its Mapping.json come from the developer
side's provider build -- ``python -m umat_oti.provider.collaborator``, the same
function the UMAT GUI's Build button calls -- so this is the hand-over of
slide 17 consumed by slide 18. The saved analysis is the presentation's J2
example: the Abaqus 2021.HF5 job ``imqrp_j2`` (one C3D8, four increments),
its ``Analysis.odb`` and the ``Analysis.inp`` that produced it.

The browser points at the three files, ticks E, nu, SIGY0 and H, chooses the
output and the region, clicks Solve, reads the full-field table and downloads
the results; ``resasm request`` is then run on the same files and the same
request and must agree number for number.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from gui_helpers import REPO_ROOT, save_screenshot, settle, streamlit_server

pytest.importorskip("playwright.sync_api", reason="playwright is not installed")
from playwright.sync_api import expect, sync_playwright  # noqa: E402

pytestmark = [pytest.mark.gui, pytest.mark.abaqus, pytest.mark.slow, pytest.mark.fortran]

APP = REPO_ROOT / "scripts" / "app.py"
#: The saved analysis. Kept outside the repository (an ODB is a binary run
#: product); RESASM_PRESENTATION_ANALYSIS points elsewhere if it moved.
ANALYSIS = Path(os.environ.get(
    "RESASM_PRESENTATION_ANALYSIS",
    REPO_ROOT.parent / "imq_abaqus" / "recovery_presentation" / "imqrp_reference" / "collaborator"))
E, SIGY0, H = 210000.0, 250.0, 2000.0


def _umat_repository() -> Path:
    configured = os.environ.get("UMAT_OTI_REPO")
    for candidate in ([Path(configured)] if configured else []) + [REPO_ROOT.parent / "UMAT_source_transformation"]:
        if (candidate / "parameter_sensitivity" / "models" / "m3_j2" / "umat.for").is_file():
            return candidate
    pytest.fail("the UMAT repository (UMAT_OTI_REPO) is needed to build the OTI object")


@pytest.fixture(scope="module")
def handed_over(tmp_path_factory):
    """Slide 17's hand-over: OTI_UMAT.obj and Mapping.json from the provider build."""
    for name in ("Analysis.odb", "Analysis.inp"):
        if not (ANALYSIS / name).is_file():
            pytest.skip(f"the saved analysis {ANALYSIS / name} is not on this machine")
    if shutil.which("abaqus") is None:
        pytest.skip("licensed Abaqus Python is needed to read the ODB")
    from umat_oti.provider import collaborator

    work = tmp_path_factory.mktemp("handover")
    parameters = collaborator.parse_parameters([
        {"parameter": "E", "PROPS index": 1, "value": E}, {"parameter": "nu", "PROPS index": 2, "value": 0.3},
        {"parameter": "SIGY0", "PROPS index": 3, "value": SIGY0}, {"parameter": "H", "PROPS index": 4, "value": H}])
    contract = collaborator.stage_model(
        _umat_repository() / "parameter_sensitivity/models/m3_j2/umat.for", work / "m3_j2",
        collaborator.provider_contract("umat.for", name="m3_j2", nstatev=1, parameters=parameters))
    completed = subprocess.run([sys.executable, "-m", "umat_oti.provider.collaborator", str(contract),
                                "--out", str(work / "out"), "--j2-branches"], capture_output=True, text=True)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    inputs = work / "collaborator_inputs"
    inputs.mkdir()
    for name in ("OTI_UMAT.obj", "Mapping.json"):
        shutil.copyfile(work / "out" / "collaborator" / name, inputs / name)
    for name in ("Analysis.odb", "Analysis.inp"):
        shutil.copyfile(ANALYSIS / name, inputs / name)
    return inputs


@pytest.fixture(scope="module")
def server(tmp_path_factory):
    work = tmp_path_factory.mktemp("resasm_gui")
    env = dict(os.environ, PYTHONPATH=os.pathsep.join(filter(None, [str(REPO_ROOT),
                                                                   os.environ.get("PYTHONPATH")])))
    with streamlit_server(APP, env=env, log=work / "server.log") as (url, pid):
        yield {"url": url, "pid": pid}


def test_slide18_point_tick_choose_solve(server, handed_over, tmp_path):
    out = tmp_path / "gui_results"
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        page = browser.new_context(viewport={"width": 1440, "height": 2200}, accept_downloads=True).new_page()
        page.goto(server["url"], wait_until="networkidle", timeout=120000)
        expect(page.get_by_role("tab", name="Sensitivity Request", exact=True)).to_have_attribute(
            "aria-selected", "true")
        panel = page.get_by_role("tabpanel", name="Sensitivity Request")
        # Point at the compiled object, the saved .odb and its setup .inp.
        for name in ("OTI_UMAT.obj", "Analysis.inp", "Analysis.odb"):
            panel.get_by_label(name + " path", exact=True).fill(str(handed_over / name))
            panel.get_by_label(name + " path", exact=True).press("Enter")
            settle(page)
        panel.get_by_label("Output directory", exact=True).fill(str(out))
        panel.get_by_label("Output directory", exact=True).press("Enter")
        settle(page)
        # Tick E, nu, SIGY0, H (read from Mapping.json beside the object).
        for name in ("E", "nu", "SIGY0", "H"):
            expect(panel.get_by_role("checkbox", name=name, exact=True)).to_be_checked()
        # Choose the output and the region.
        panel.get_by_role("combobox", name="output", exact=True).click()
        page.get_by_role("option", name="displacement U1", exact=True).click()
        settle(page)
        panel.get_by_role("combobox", name="region", exact=True).click()
        page.get_by_role("option", name="node set LOADED", exact=True).click()
        settle(page)
        panel.get_by_role("button", name="Solve", exact=True).click()
        expect(panel.get_by_text(re.compile(r"^Solved — d\(mean_U1_LOADED\)/dp at all 4 locations"))).to_be_visible(
            timeout=600000)
        settle(page)
        expect(panel.get_by_text("Executed: 1 scalar results. Independent validation: not run.",
                                 exact=True)).to_be_visible()
        downloads = {}
        for name in ("sensitivity_results.json", "run_report.txt"):
            with page.expect_download() as event:
                panel.get_by_role("button", name=name, exact=True).click()
            downloads[name] = tmp_path / name
            event.value.save_as(downloads[name])
        save_screenshot(panel, "resasm_solve.png")

        # The same analysis, the equivalent plastic strain at every integration
        # point of the whole mesh: a new output directory, one more Solve.
        panel.get_by_label("Output directory", exact=True).fill(str(tmp_path / "gui_sdv"))
        panel.get_by_label("Output directory", exact=True).press("Enter")
        settle(page)
        output = panel.get_by_role("combobox", name="output", exact=True)
        output.click()
        output.fill("state variable SDV1")  # the list is long; type to filter it
        page.get_by_role("option", name="state variable SDV1", exact=True).click()
        settle(page)
        expect(panel.get_by_role("combobox", name="region", exact=True)).to_match_aria_snapshot(
            '- combobox "region": whole mesh')
        panel.get_by_role("button", name="Solve", exact=True).click()
        expect(panel.get_by_text(re.compile(r"^Solved — d\(mean_SDV1_all\)/dp at all 8 locations"))).to_be_visible(
            timeout=600000)
        settle(page)
        save_screenshot(panel, "resasm_solve_sdv_all_points.png")
        browser.close()
    sdv = json.loads((tmp_path / "gui_sdv" / "sensitivity_results.json").read_text())
    (sdv_row,) = sdv["results"]
    assert sdv_row["value"] == pytest.approx((300 - SIGY0) / H, rel=1e-5)
    assert sdv_row["derivatives"]["SIGY0"] == pytest.approx(-1 / H, rel=1e-5)
    assert sdv_row["derivatives"]["H"] == pytest.approx(-(300 - SIGY0) / H ** 2, rel=1e-5)

    summary = json.loads(downloads["sensitivity_results.json"].read_text())
    assert "Status: executed successfully" in downloads["run_report.txt"].read_text()
    (row,) = summary["results"]
    # The uniaxial reference (see test_solve_screen.py); ODB float32
    # displacements bound the agreement.
    assert row["value"] == pytest.approx(300 / E + (300 - SIGY0) / H, rel=1e-5)
    assert row["derivatives"]["SIGY0"] == pytest.approx(-1 / H, rel=1e-5)
    assert row["derivatives"]["H"] == pytest.approx(-(300 - SIGY0) / H ** 2, rel=1e-5)
    assert row["derivatives"]["E"] == pytest.approx(-300 / E ** 2, rel=1e-4)
    # The command line on the same four inputs and the same request.
    request = tmp_path / "sensitivity_request.json"
    request.write_text(json.dumps(summary["request"]))
    completed = subprocess.run(
        [sys.executable, "-m", "residual_core.ui.cli", "request", "--model", str(handed_over / "Analysis.inp"),
         "--odb", str(handed_over / "Analysis.odb"), "--material", str(handed_over / "OTI_UMAT.obj"),
         "--request", str(request), "--out", str(tmp_path / "cli_results")],
        capture_output=True, text=True,
        env=dict(os.environ, PYTHONPATH=os.pathsep.join(filter(None, [str(REPO_ROOT),
                                                                     os.environ.get("PYTHONPATH")]))))
    assert completed.returncode == 0, completed.stdout + completed.stderr
    cli = json.loads((tmp_path / "cli_results" / "sensitivity_results.json").read_text())
    assert cli["results"] == summary["results"]
