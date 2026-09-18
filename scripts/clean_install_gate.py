"""Verify working-tree wheels together, never claim a committed clean-clone gate."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
import socket
import urllib.error
import urllib.request


BRANCH = "integration/imqcam-recovery-2026-09-18"
PUBLIC_OUTPUTS = ("sensitivity_results.json", "sensitivity_tables.csv", "run_report.txt")


def clean_environment(home: Path) -> dict[str, str]:
    environment = {key: value for key, value in os.environ.items()
                   if not key.startswith(("PYTHON", "PIP_"))
                   and key not in {"VIRTUAL_ENV", "CONDA_PREFIX", "PYOTI_PATH", "OTILIB_ROOT", "UMAT_OTI_REPO"}}
    environment.update(HOME=str(home), XDG_CONFIG_HOME=str(home / "config"),
                       XDG_CACHE_HOME=str(home / "cache"), PYTHONNOUSERSITE="1",
                       PIP_CONFIG_FILE=os.devnull, PIP_DISABLE_PIP_VERSION_CHECK="1")
    return environment


def external_work(path: Path, repositories: list[Path]) -> Path:
    path = path.resolve()
    if any(path == root.resolve() or root.resolve() in path.parents for root in repositories):
        raise ValueError("gate work directory must be outside both repositories")
    if path.exists():
        raise ValueError("gate requires a NEW directory, not an existing environment")
    return path


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Gate:
    def __init__(self, work: Path):
        self.work = work
        self.report = {"kind": "wheel-from-working-tree", "final_branch_clean_clone": False,
                       "passed": False, "commands": [], "work": str(work)}
        self.environment = clean_environment(work / "home")
        (work / "home").mkdir()
        (work / "logs").mkdir()

    def save(self):
        (self.work / "report.json").write_text(json.dumps(self.report, indent=2) + "\n")

    def run(self, command, *, cwd=None):
        command = [str(argument) for argument in command]
        entry = {"argv": command, "cwd": str(cwd or self.work)}
        self.report["commands"].append(entry)
        self.save()
        completed = subprocess.run(command, cwd=cwd or self.work, env=self.environment,
                       text=True, capture_output=True)
        log = self.work / "logs" / f"{len(self.report['commands']):03d}.log"
        log.write_text(completed.stdout + completed.stderr)
        entry.update(returncode=completed.returncode, log=str(log))
        self.save()
        displayed = command[:command.index("-c") + 1] + ["<probe; see log>"] if "-c" in command else command
        print(f"[{completed.returncode}] {' '.join(displayed)}", flush=True)
        if completed.returncode:
            raise RuntimeError(f"command failed; see {log}\n{(completed.stdout + completed.stderr)[-6000:]}")
        return completed.stdout


INSTALLED_PROBE = r'''
import importlib, importlib.metadata, importlib.resources, json, pathlib, site, sys
roots = [pathlib.Path(path).resolve() for path in site.getsitepackages()]
origins = {}
for name in ('residual_core', 'resasm_user', 'umat_oti', 'umat_oti.provider.build',
             'residual_core.replay.presentation', 'residual_core.app.streamlit_app',
             'umat_oti.app.streamlit_app'):
    path = pathlib.Path(importlib.import_module(name).__file__).resolve()
    assert any(root in path.parents for root in roots), (name, path, roots)
    origins[name] = str(path)
for distribution in importlib.metadata.distributions():
    direct = distribution.read_text('direct_url.json')
    assert not direct or not json.loads(direct).get('dir_info', {}).get('editable'), direct
resources = {}
requirements = {
    'residual_core': ['example_manifest.json', 'umat_adapter_fortran/umat_driver.f90',
                      'umat_adapter_fortran/aba_stubs/*.f*'],
    'residual_core.replay': ['contract/CONTRACT_VERSION.json', 'contract/resasm_mat_abi_v1.h'],
    'umat_oti.oti': ['support/*.f90', 'support/pyoti_templates/*.f90'],
    'umat_oti.contract': ['schemas/*.json', 'contract_lock.json'],
}
for package, patterns in requirements.items():
    root = pathlib.Path(str(importlib.resources.files(package)))
    resources[package] = []
    for pattern in patterns:
        matches = list(root.glob(pattern))
        assert matches, (package, pattern)
        assert all(path.stat().st_size for path in matches)
        resources[package].extend(str(path) for path in matches)
assert sys.prefix != sys.base_prefix
assert not site.ENABLE_USER_SITE
for path in sys.path:
    assert not any(pathlib.Path(root).resolve() == pathlib.Path(path).resolve()
                   or pathlib.Path(root).resolve() in pathlib.Path(path).resolve().parents
                   for root in sys.argv[1:]), (path, sys.argv[1:])
print(json.dumps(dict(origins=origins, resources=resources, sys_path=sys.path,
                     executable=sys.executable, version=sys.version), indent=2))
'''


GUI_PROBE = r'''
import importlib.util, json
from streamlit.testing.v1 import AppTest
results = {}
for name in ('residual_core.app.streamlit_app', 'umat_oti.app.streamlit_app'):
    path = importlib.util.find_spec(name).origin
    app = AppTest.from_file(path, default_timeout=60).run()
    assert not app.exception, [str(error) for error in app.exception]
    assert len(app.title) or len(app.header), name
    results[name] = dict(path=path, titles=[item.value for item in app.title],
                         tabs=[item.label for item in app.tabs])
print(json.dumps(results, indent=2))
'''


SOURCE_DENIED_PROBE = r'''
import os, pathlib, runpy, sys
console = sys.argv.pop(1)
allowed_shim = (pathlib.Path.cwd() / 'source_denied_results/private/link/path_shim.for').resolve()
def deny_source(event, arguments):
    if event == 'open' and isinstance(arguments[0], (str, bytes)):
        path = pathlib.Path(os.fsdecode(arguments[0])).resolve()
        if path.suffix.lower() in ('.for', '.f', '.f90') and path != allowed_shim:
            raise PermissionError('material source denied: ' + str(path))
sys.addaudithook(deny_source)
try:
    open('unavailable_material.for')
except PermissionError:
    print('PRIVATE SOURCE READ DENIAL ACTIVE', flush=True)
else:
    raise AssertionError('source denial is inactive')
runpy.run_path(console, run_name='__main__')
'''


def check_servers(gate, python, modules):
    reports = {}
    for name, module in modules.items():
        with socket.socket() as listener:
            listener.bind(("127.0.0.1", 0))
            port = listener.getsockname()[1]
        command = [str(python), "-I", "-m", "streamlit", "run", module["path"],
                   "--server.address=127.0.0.1", f"--server.port={port}",
                   "--server.headless=true", "--browser.gatherUsageStats=false"]
        report = {"argv": command, "url": f"http://127.0.0.1:{port}", "ready": False}
        reports[name] = report
        gate.report["servers"] = reports
        with (gate.work / "logs" / (name + ".log")).open("w") as log:
            server = subprocess.Popen(command, cwd=gate.work, env=gate.environment,
                                      stdout=log, stderr=subprocess.STDOUT)
            report["pid"] = server.pid
            try:
                deadline = time.monotonic() + 30
                while time.monotonic() < deadline:
                    if server.poll() is not None:
                        raise RuntimeError(f"installed GUI exited before readiness: {name}")
                    try:
                        with urllib.request.urlopen(report["url"] + "/_stcore/health", timeout=1) as response:
                            if response.status == 200:
                                report["ready"] = True
                                break
                    except (urllib.error.URLError, TimeoutError):
                        time.sleep(0.1)
                if not report["ready"]:
                    raise RuntimeError(f"installed GUI readiness timeout: {name}")
            finally:
                server.terminate()
                try:
                    server.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    server.kill()
                    server.wait(timeout=5)
                report["stopped"] = server.poll() is not None
                gate.save()


def execute(args, gate: Gate):
    repositories = [args.ra_repo.resolve(), args.umat_repo.resolve()]
    gate.report["repositories"] = {}
    for root in repositories:
        branch = gate.run(["git", "-C", root, "branch", "--show-current"]).strip()
        if branch != BRANCH:
            raise ValueError(f"expected {BRANCH}, found {branch} at {root}")
        gate.report["repositories"][str(root)] = {
            "head": gate.run(["git", "-C", root, "rev-parse", "HEAD"]).strip(),
            "status": gate.run(["git", "-C", root, "status", "--short"])}
    gate.run([args.python, "-I", "-c",
              "import ctypes, ssl, venv, sys; assert sys.version_info >= (3, 10); print(sys.version)"])
    for executable in ("gfortran", args.abaqus):
        if not shutil.which(executable, path=gate.environment.get("PATH")):
            raise ValueError(f"required external executable not found: {executable}")
    gate.run(["gfortran", "--version"])
    gate.run([args.python, "-I", "-m", "venv", gate.work / "env"])
    python = gate.work / "env/bin/python"
    gate.run([python, "-I", "-m", "pip", "install", "--upgrade", "pip"])
    wheels = gate.work / "wheels"
    wheels.mkdir()
    for root in repositories:
        gate.run([python, "-I", "-m", "pip", "wheel", "--no-deps", "--wheel-dir", wheels, root])
    wheel_paths = sorted(wheels.glob("*.whl"))
    if len(wheel_paths) != 2:
        raise ValueError(f"expected exactly two project wheels, found {wheel_paths}")
    gate.report["wheels"] = {path.name: digest(path) for path in wheel_paths}
    extras = {"residual_assembler": "gui,test,yaml", "umat_oti": "test"}
    requirements = [str(path) + "[" + next(value for prefix, value in extras.items()
                                          if path.name.startswith(prefix + "-")) + "]"
                    for path in wheel_paths]
    gate.run([python, "-I", "-m", "pip", "install", *requirements])
    gate.run([python, "-I", "-m", "pip", "check"])
    gate.report["installed"] = json.loads(gate.run([python, "-I", "-c", INSTALLED_PROBE, *repositories]))
    gate.report["versions"] = json.loads(gate.run([python, "-I", "-m", "pip", "list", "--format=json"]))
    developer = gate.work / "public_model"
    developer.mkdir()
    source = repositories[1] / "parameter_sensitivity/models/m3_j2"
    for name in ("contract_v2.json", "umat.for"):
        shutil.copy2(source / name, developer / name)
    built = json.loads(gate.run([gate.work / "env/bin/umat-oti-provider", "build",
                                developer / "contract_v2.json", "--out", gate.work / "provider"]))
    collaborator = gate.work / "collaborator"
    collaborator.mkdir()
    for source_path, name in ((repositories[0] / "examples/presentation_request/Analysis.inp", "Analysis.inp"),
                              (repositories[0] / "examples/presentation_request/sensitivity_request.json", "sensitivity_request.json"),
                              (args.odb.resolve(), "Analysis.odb"),
                              (Path(built["object"]), "OTI_UMAT.obj"),
                              (Path(built["contract"]), "Mapping.json")):
        shutil.copy2(source_path, collaborator / name)
    gate.report["inputs"] = {path.name: digest(path) for path in collaborator.iterdir()}
    gate.run([gate.work / "env/bin/resasm", "request", "--model", "Analysis.inp", "--odb", "Analysis.odb",
              "--material", "OTI_UMAT.obj", "--request", "sensitivity_request.json", "--out", "results",
              "--abaqus", args.abaqus], cwd=collaborator)
    output = collaborator / "results"
    for name in PUBLIC_OUTPUTS:
        if not (output / name).is_file() or not (output / name).stat().st_size:
            raise ValueError(f"missing presentation output: {name}")
    result = json.loads((output / PUBLIC_OUTPUTS[0]).read_text())
    if result["status"] != "completed" or result["metadata"]["verified"] is not False:
        raise ValueError("ordinary execution must be completed, not labelled independently verified")
    if result["scope"]["history_increments_replayed"] != 4 or len(result["results"]) != 4:
        raise ValueError("presentation fixture must replay all four increments and return four outputs")
    displacement = result["results"][0]
    expected = {"E": -300 / 210000 ** 2, "SIGY0": -1 / 2000, "H": -50 / 2000 ** 2}
    errors = {name: abs(displacement["derivatives"][name] - value) / abs(value)
              for name, value in expected.items()}
    if max(errors.values()) >= 2e-5 or abs(displacement["derivatives"]["nu"]) >= 1e-8:
        raise ValueError(f"independent uniaxial J2 derivative check failed: {errors}")
    gate.report["presentation"] = {"outputs": list(PUBLIC_OUTPUTS), "analytic_errors": errors,
                                   "production_jobs_launched": 0, "odb_source": str(args.odb.resolve())}
    gate.run([python, "-I", "-c", SOURCE_DENIED_PROBE, gate.work / "env/bin/resasm", "request",
              "--model", "Analysis.inp", "--odb", "Analysis.odb", "--material", "OTI_UMAT.obj",
              "--request", "sensitivity_request.json", "--out", "source_denied_results",
              "--abaqus", args.abaqus], cwd=collaborator)
    for name in PUBLIC_OUTPUTS:
        if (output / name).read_bytes() != (collaborator / "source_denied_results" / name).read_bytes():
            raise ValueError(f"source-denied console output differs: {name}")
    gate.report["presentation"]["source_denied_outputs_identical"] = True
    gate.report["gui"] = json.loads(gate.run([python, "-I", "-c", GUI_PROBE]))
    check_servers(gate, python, gate.report["gui"])
    gate.report["passed"] = True


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ra-repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--umat-repo", type=Path, required=True)
    parser.add_argument("--python", required=True, help="healthy standalone Python >=3.10 with venv, ctypes, ssl")
    parser.add_argument("--odb", type=Path, required=True, help="genuine ODB matching the public presentation deck")
    parser.add_argument("--abaqus", default="abaqus")
    parser.add_argument("--work", type=Path)
    args = parser.parse_args(argv)
    repositories = [args.ra_repo, args.umat_repo]
    if args.work:
        work = external_work(args.work, repositories)
        work.mkdir(parents=True)
    else:
        work = Path(tempfile.mkdtemp(prefix="imqr_install_"))
    gate = Gate(work)
    try:
        execute(args, gate)
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        gate.report["error"] = str(error)
        print(str(error), file=sys.stderr)
        return 1
    finally:
        gate.save()
        print(f"Evidence: {work / 'report.json'}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())