"""Build both wheels from clean trees, install them in a new venv, and run the workflows from them."""

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
    hasher = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            hasher.update(block)
    return hasher.hexdigest()


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
    published = []
    for root in repositories:
        branch = gate.run(["git", "-C", root, "branch", "--show-current"]).strip()
        if args.branch and branch != args.branch:
            raise ValueError(f"expected branch {args.branch}, found {branch or 'a detached HEAD'} at {root}")
        status = gate.run(["git", "-C", root, "status", "--porcelain", "--untracked-files=all"])
        if status.strip():
            # wheels are built from the tree, so an uncommitted or untracked
            # file would be installed and the result would describe no commit
            raise ValueError(f"working tree is not clean at {root}:\n{status}")
        head = gate.run(["git", "-C", root, "rev-parse", "HEAD"]).strip()
        remote_head = subprocess.run(["git", "-C", str(root), "rev-parse", "--verify", "--quiet",
                                      f"refs/remotes/origin/{branch}"],
                                     text=True, capture_output=True).stdout.strip() if branch else ""
        remote_url = subprocess.run(["git", "-C", str(root), "remote", "get-url", "origin"],
                                    text=True, capture_output=True).stdout.strip()
        gate.report["repositories"][str(root)] = {
            "branch": branch, "head": head, "status": status, "origin": remote_url,
            "origin_branch_head": remote_head}
        published.append(bool(args.branch) and head == remote_head)
    # True only when both wheels come from clean trees whose commit is the
    # remote-tracking head of the named branch (a fresh clone of it, or a
    # checkout identical to it as of its last fetch)
    gate.report["final_branch_clean_clone"] = all(published)
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
    if args.cantilever:
        check_cantilever(args, gate, built, repositories)
    gate.report["passed"] = True


#: Largest homogeneity residual accepted on the re-equilibrated J2 cantilever,
#: relative to the largest term. Measured 1.0e-12 on 2026-09-18 (roundoff of a
#: 40-increment chain through a 7,293-DOF solve); a single derivative wrong by
#: 1e-6 leaves about 1e-9.
HOMOGENEITY_BOUND = 1e-10


def check_cantilever(args, gate: Gate, built, repositories):
    """The full-size J2 cantilever of slide 39, from the installed wheels.

    ``args.cantilever`` holds ``j2/claude_j2_nominal.inp`` and its ODB, as
    ``examples/presentation_cantilevers/README.md`` produces them in Abaqus.
    The collaborator command runs on them (routed to the history engine),
    then the history engine re-equilibrates the recorded increments and every
    output must obey Euler's identity for the J2 model -- homogeneous of
    degree one in (E, SIGY0, H) at fixed nu, so under prescribed displacements
    sum p dQ/dp is Q for reactions, stresses and von Mises and 0 for
    displacements and plastic strain, at every increment. Nothing in the
    engine uses that identity.
    """
    source = args.cantilever.resolve() / "j2"
    deck, odb = source / "claude_j2_nominal.inp", source / "claude_j2_nominal.odb"
    for path in (deck, odb):
        if not path.is_file():
            raise ValueError(f"cantilever input missing: {path}")
    request = repositories[0] / "examples/presentation_cantilevers/j2_request.json"
    folder = gate.work / "cantilever"
    folder.mkdir()
    shutil.copy2(deck, folder / "Analysis.inp")
    (folder / "Analysis.odb").symlink_to(odb)          # 0.6 GB: linked, digested
    shutil.copy2(request, folder / "sensitivity_request.json")
    shutil.copy2(built["object"], folder / "OTI_UMAT.obj")
    shutil.copy2(built["contract"], folder / "Mapping.json")
    record = {"inputs": {name: digest(folder / name) for name in
                         ("Analysis.inp", "Analysis.odb", "sensitivity_request.json", "OTI_UMAT.obj",
                          "Mapping.json")}, "odb_source": str(odb)}
    gate.report["cantilever"] = record
    started = time.monotonic()
    gate.run([gate.work / "env/bin/resasm", "request", "--model", "Analysis.inp", "--odb", "Analysis.odb",
              "--material", "OTI_UMAT.obj", "--request", "sensitivity_request.json", "--out", "results",
              "--abaqus", args.abaqus], cwd=folder)
    record["request_seconds"] = round(time.monotonic() - started, 1)
    result = json.loads((folder / "results/sensitivity_results.json").read_text())
    scope = result["scope"]
    if result["status"] != "completed" or scope.get("increments_replayed") != 40 \
            or scope.get("integration_points") != 12288 or scope.get("dof") != 7497:
        raise ValueError(f"cantilever request did not replay the slide-39 model: {result['status']} {scope}")
    started = time.monotonic()
    gate.run([gate.work / "env/bin/resasm", "history", "--model", "Analysis.inp",
              "--fields", "results/private/fields.npz", "--material", "OTI_UMAT.obj",
              "--request", "sensitivity_request.json", "--out", "reequilibrated", "--reequilibrate"],
             cwd=folder)
    record["reequilibrated_seconds"] = round(time.monotonic() - started, 1)
    polished = json.loads((folder / "reequilibrated/sensitivity_results.json").read_text())
    degree = {"RF": 1, "S": 1, "MISES": 1, "U": 0, "SDV": 0}
    worst, plastic = 0.0, 0
    for row in polished["results"]:
        weighted = [row["weighted"][name] for name in ("E", "SIGY0", "H")]
        largest = max(max(abs(value) for value in weighted), abs(row["value"]))
        if largest:
            worst = max(worst, abs(sum(weighted) - degree[row["field"]] * row["value"]) / largest)
        plastic += row["output"] == "eqplas_max" and row["value"] > 0
    record.update(homogeneity_residual=worst, homogeneity_bound=HOMOGENEITY_BOUND,
                  plastic_increments=plastic, outputs_checked=len(polished["results"]))
    if not plastic or worst > HOMOGENEITY_BOUND:
        raise ValueError(f"cantilever homogeneity check failed: residual {worst:.3e}, "
                         f"{plastic} plastic increments")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ra-repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--umat-repo", type=Path, required=True)
    parser.add_argument("--python", required=True, help="healthy standalone Python >=3.10 with venv, ctypes, ssl")
    parser.add_argument("--odb", type=Path, required=True, help="genuine ODB matching the public presentation deck")
    parser.add_argument("--abaqus", default="abaqus")
    parser.add_argument("--work", type=Path)
    parser.add_argument("--branch", help="branch both repositories must be on; with it, the report says "
                                         "whether both commits are that branch's published head")
    parser.add_argument("--cantilever", type=Path,
                        help="directory holding j2/claude_j2_nominal.inp and .odb "
                             "(examples/presentation_cantilevers); adds the full-size J2 check")
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