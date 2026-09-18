"""Slide 8: paired Abaqus DDSDDE of the benchmark UMATs (claim 4)."""

import json
import shutil

import pytest

from presentation import claim4_benchmark_ddsdde as claim4

pytestmark = [pytest.mark.integration]


def test_every_benchmark_is_prepared_or_carries_its_reason(out_dirs):
    out, work = out_dirs
    assert claim4.main(["--out", str(out), "--work", str(work)]) == 0
    payload = json.loads((out / "claim4_benchmark_ddsdde.json").read_text())
    assert payload["summary"]["cases"] == 19
    for row in payload["rows"]:
        if row["status"] == "prepared":
            assert row["compare_outputs_used"] == ["STRESS", "STATEV", "DDSDDE", "CONVERGENCE"]
        else:
            assert row["status"] == "transform_failed"
            assert row["transform_success"] is False


@pytest.mark.abaqus
def test_elastic_benchmark_end_to_end_in_abaqus(out_dirs):
    if not shutil.which("abaqus"):
        pytest.fail("abaqus-marked test run without Abaqus on PATH")
    out, work = out_dirs
    assert claim4.main(["--abaqus", "--case", "elastic", "--out", str(out), "--work", str(work)]) == 0
    row = json.loads((out / "claim4_benchmark_ddsdde.json").read_text())["rows"][0]
    assert row["runs"]["orig"]["sta_completed"] and row["runs"]["oti"]["sta_completed"]
    assert row["compared"]["overall_pass"] is True
    assert row["compared"]["ddsdde_max_rel"] == 0.0
