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


def test_provider_checkout_override_is_independent_of_names(tmp_path, monkeypatch):
    from repository_paths import umat_repo_root

    checkout = tmp_path / "arbitrary-producer-name"
    (checkout / "src/umat_oti").mkdir(parents=True)
    monkeypatch.setenv("UMAT_OTI_REPO", str(checkout))
    monkeypatch.chdir(tmp_path)
    assert umat_repo_root() == checkout
    monkeypatch.setenv("UMAT_OTI_REPO", str(tmp_path / "missing"))
    with pytest.raises(FileNotFoundError, match="set UMAT_OTI_REPO"):
        umat_repo_root()


def test_provider_checkout_default_is_canonical(tmp_path, monkeypatch):
    import repository_paths

    checkout = tmp_path / "UMAT_source_transformation"
    (checkout / "src/umat_oti").mkdir(parents=True)
    monkeypatch.delenv("UMAT_OTI_REPO", raising=False)
    monkeypatch.setattr(repository_paths, "__file__", str(tmp_path / "consumer/tests/integration/repository_paths.py"))
    assert repository_paths.umat_repo_root() == checkout


@pytest.mark.parametrize("template", ["python", "blackbox", "blackbox-order2", "cpp", "fortran"])
def test_installed_template_data_is_declared_and_discoverable(tmp_path, monkeypatch, template):
    import shutil
    import sysconfig
    try:
        import tomllib
    except ModuleNotFoundError:
        import tomli as tomllib
    from residual_core.ui import cli

    repository = Path(__file__).resolve().parents[2]
    settings = tomllib.loads((repository / "pyproject.toml").read_text())
    folder = cli._TEMPLATES[template]
    relative = "share/residual-assembler/templates/" + folder
    declared = settings["tool"]["setuptools"]["data-files"][relative]
    assert all((repository / filename).is_file() for filename in declared)
    assert "templates/" + folder + "/resasm.yml" in declared
    assert "templates/" + folder + "/solution.npy" in declared
    destination = tmp_path / relative
    destination.mkdir(parents=True)
    for filename in declared:
        shutil.copy2(repository / filename, destination)
    original_get_path = sysconfig.get_path
    monkeypatch.setattr(sysconfig, "get_path", lambda name: str(tmp_path) if name == "data" else original_get_path(name))
    monkeypatch.chdir(tmp_path)
    assert cli._templates_root() == str(destination.parent)
    output = tmp_path / "user-job"
    assert cli.main(["init", "--template", template, "--out", str(output)]) == 0
    assert (output / "resasm.yml").read_bytes() == (destination / "resasm.yml").read_bytes()
    assert (output / "solution.npy").read_bytes() == (destination / "solution.npy").read_bytes()