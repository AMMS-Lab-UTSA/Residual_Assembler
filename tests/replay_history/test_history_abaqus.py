"""Live Abaqus: run a small analysis, export its ODB, replay it (marker ``abaqus``).

Deselected by the offline suite (``-m "not abaqus"``). Needs a licensed
``abaqus`` launcher with a Fortran compiler for user subroutines. Job names
start with ``$RESASM_ABAQUS_JOB_PREFIX`` (default ``resasm_``).
"""
import json
import os
import shutil
import subprocess

import numpy as np
import pytest

from residual_core.ui.cli import main as resasm

from conftest import EXAMPLE, RA_ROOT, umat_repository

pytestmark = pytest.mark.abaqus


def _run(work, job, deck, umat):
    (work / (job + ".inp")).write_text(deck)
    subprocess.run(["abaqus", "job=" + job, "input=" + job + ".inp", "user=" + str(umat), "double=both",
                    "cpus=1", "interactive"], cwd=work, capture_output=True, text=True)
    status = (work / (job + ".sta")).read_text()
    assert "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" in status


def test_odb_export_and_command_line_replay(provider_factory, tmp_path):
    if shutil.which("abaqus") is None:
        pytest.fail("the abaqus marker was selected but no abaqus launcher is on PATH")
    import sys
    sys.path.insert(0, str(EXAMPLE))
    from make_beam_deck import deck
    prefix = os.environ.get("RESASM_ABAQUS_JOB_PREFIX", "resasm_")
    work = tmp_path / "abaqus"
    work.mkdir()
    umat = umat_repository() / "parameter_sensitivity" / "models" / "m3_j2" / "umat.for"
    job = prefix + "history_replay_test"
    text = deck((6, 2, 2), 0.04, 4, [200000.0, 0.3, 250.0, 2000.0], 1)
    _run(work, job, text, umat)
    export = work / "fields.npz"
    completed = subprocess.run(["abaqus", "python", str(RA_ROOT / "residual_core/replay/odb_export_npz.py"),
                                "--", str(work / (job + ".odb")), str(export)], cwd=work,
                               capture_output=True, text=True)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    data = np.load(export)
    assert data["U"].shape == (5, 63, 3) and data["S"].shape == (5, 24, 8, 6)
    label = lambda key: np.asarray(data[key]).astype(str).item()      # Abaqus Python 2 writes bytes
    assert label("precision") == "float32" and label("step_name").upper() == "PUSH"
    obj, contract, _ = provider_factory("m3_j2")
    request = tmp_path / "request.json"
    request.write_text(json.dumps({
        "outputs": [{"name": "tip_RF2", "field": "RF", "component": 2, "reduction": "sum",
                     "domain": {"nset": "TIP"}}],
        "parameters": "ALL", "domain": {"nodes": "ALL", "elements": "ALL"}, "increments": "ALL"}))
    (work / "Analysis.inp").write_text(text)
    out = tmp_path / "results"
    code = resasm(["history", "--model", str(work / "Analysis.inp"), "--odb", str(work / (job + ".odb")),
                   "--material", str(obj), "--request", str(request), "--out", str(out)])
    assert code == 0
    report = (out / "run_report.txt").read_text()
    assert "Equilibrium passed: yes" in report and "Abaqus comparison available: yes" in report
    assert (out / "private" / "fields.npz").is_file()
    exported = np.load(out / "private" / "fields.npz")
    np.testing.assert_array_equal(exported["S"], data["S"])
