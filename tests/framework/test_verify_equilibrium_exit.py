"""Public verification must fail on an unloaded cube with nonzero stress."""
import json
import subprocess
import sys
from pathlib import Path

import pytest


EXAMPLE = (Path(__file__).resolve().parents[2] / "residual_core" / "examples"
           / "minimal_c3d8_stress_driven")


def verify(fields, *extra):
    return subprocess.run(
        [sys.executable, "-m", "residual_core.ui.cli", "verify",
         str(EXAMPLE / "model.json"), "--fields", str(fields), *extra],
        capture_output=True, text=True, check=False)


def test_unbalanced_stressed_cube_has_failing_exit_code():
    result = verify(EXAMPLE / "fields.json")
    assert result.returncode == 1, result.stdout + result.stderr
    assert "equilibrium          = FAIL" in result.stdout
    assert "7.071068e+01" in result.stdout
    assert "reaction reference   = NOT CHECKED" in result.stdout


def test_zero_stress_unloaded_cube_is_equilibrated(tmp_path):
    fields = tmp_path / "fields.json"
    fields.write_text(json.dumps({"stress_ip": {"1": [[0.] * 6] * 8}}))
    result = verify(fields)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "equilibrium          = PASS" in result.stdout


@pytest.mark.parametrize("tolerance", ["-1", "nan", "inf"])
def test_invalid_tolerance_cannot_turn_failure_into_success(tolerance):
    result = verify(EXAMPLE / "fields.json", "--atol", tolerance)
    assert result.returncode == 2
    assert "finite and nonnegative" in result.stderr


def test_nonfinite_stress_never_passes(tmp_path):
    fields = tmp_path / "fields.json"
    fields.write_text(json.dumps({"stress_ip": {"1": [[float("nan")] * 6] * 8}}))
    result = verify(fields)
    assert result.returncode != 0
    assert "equilibrium          = PASS" not in result.stdout


@pytest.mark.parametrize("command", ["assemble", "verify"])
@pytest.mark.parametrize("content, message", [
    (None, "cannot read the field export"),
    ("not json", "is not valid JSON"),
    ("[1, 2, 3]", "holds no integration-point stress"),
    ('{"schema": "resasm-neutral-model/1"}', "key 'schema' is not an element id"),
    ('{"stress_ip": {"1": "S11"}}', "is not a numeric table"),
])
def test_an_unusable_field_export_is_a_one_line_error(tmp_path, command, content, message):
    """Regression: a missing --fields file, or one in another layout, ended in
    a Python traceback (exit 1)."""
    fields = tmp_path / "fields.json"
    if content is not None:
        fields.write_text(content)
    result = subprocess.run(
        [sys.executable, "-m", "residual_core.ui.cli", command,
         str(EXAMPLE / "model.json"), "--fields", str(fields)],
        capture_output=True, text=True, check=False)
    assert result.returncode == 2, result.stdout + result.stderr
    assert "Traceback" not in result.stderr
    lines = result.stderr.strip().splitlines()
    assert len(lines) == 1 and lines[0].startswith("ERROR:") and message in lines[0], lines


def test_gui_propagates_failed_equilibrium():
    pytest.importorskip("streamlit")
    from residual_core.app.streamlit_app import _run
    result = _run(["verify", str(EXAMPLE / "model.json"),
                   "--fields", str(EXAMPLE / "fields.json")])
    assert result.code == 1
    assert "equilibrium          = FAIL" in result.stdout
