"""Guard isolation failures that editable installs otherwise hide."""

import importlib.util
from pathlib import Path

import pytest


SPEC = importlib.util.spec_from_file_location(
    "clean_install_gate", Path(__file__).resolve().parents[2] / "scripts/clean_install_gate.py")
gate = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(gate)


def test_environment_removes_development_imports(monkeypatch, tmp_path):
    for name in ("PYTHONPATH", "PYTHONHOME", "VIRTUAL_ENV", "PYOTI_PATH", "OTILIB_ROOT", "UMAT_OTI_REPO", "PIP_TARGET"):
        monkeypatch.setenv(name, "/private/development")
    environment = gate.clean_environment(tmp_path)
    assert environment["HOME"] == str(tmp_path)
    assert environment["PYTHONNOUSERSITE"] == "1"
    assert "/private/development" not in environment.values()


def test_work_requires_new_external_directory(tmp_path):
    repository = tmp_path / "repo"
    repository.mkdir()
    with pytest.raises(ValueError, match="outside"):
        gate.external_work(repository / "new", [repository])
    with pytest.raises(ValueError, match="NEW"):
        gate.external_work(tmp_path, [repository])
    assert gate.external_work(tmp_path / "new", [repository]) == tmp_path / "new"


def test_command_failure_is_recorded_not_passed(tmp_path):
    import sys
    runner = gate.Gate(tmp_path)
    with pytest.raises(RuntimeError, match="command failed"):
        runner.run([sys.executable, "-I", "-c", "raise SystemExit(7)"])
    assert runner.report["commands"][0]["returncode"] == 7
    assert runner.report["passed"] is False
    assert runner.report["final_branch_clean_clone"] is False


def test_json_stdout_is_not_corrupted_by_diagnostics(tmp_path):
    import json
    import sys
    runner = gate.Gate(tmp_path)
    output = runner.run([sys.executable, "-I", "-c",
                         "import sys; print('{}'); print('diagnostic', file=sys.stderr)"])
    assert json.loads(output) == {}
    assert "diagnostic" in (tmp_path / "logs/001.log").read_text()


def test_source_denial_rejects_private_material(tmp_path):
    import sys
    runner = gate.Gate(tmp_path)
    console = tmp_path / "console.py"
    console.write_text("open('private_umat.f90')")
    with pytest.raises(RuntimeError, match="material source denied"):
        runner.run([sys.executable, "-I", "-c", gate.SOURCE_DENIED_PROBE, console])


def test_source_denial_allows_only_generated_link_shim(tmp_path):
    import sys
    runner = gate.Gate(tmp_path)
    link = tmp_path / "source_denied_results/private/link"
    link.mkdir(parents=True)
    console = tmp_path / "console.py"
    console.write_text("from pathlib import Path\nPath('source_denied_results/private/link/path_shim.for').write_text('shim')")
    runner.run([sys.executable, "-I", "-c", gate.SOURCE_DENIED_PROBE, console])
    assert (link / "path_shim.for").read_text() == "shim"