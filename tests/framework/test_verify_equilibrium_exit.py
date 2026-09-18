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


def test_gui_propagates_failed_equilibrium():
    pytest.importorskip("streamlit")
    from residual_core.app.streamlit_app import _run
    result = _run(["verify", str(EXAMPLE / "model.json"),
                   "--fields", str(EXAMPLE / "fields.json")])
    assert result.code == 1
    assert "equilibrium          = FAIL" in result.stdout
