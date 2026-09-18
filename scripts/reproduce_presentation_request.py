"""One licensed developer reference job, then source-denied collaborator replay."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT.parent / "imq_abaqus/recovery_presentation"
JOB = "imqrp_j2"


def prepare(work):
    from umat_oti.provider.build import build_provider
    if work.exists():
        raise ValueError("developer preparation requires a new work directory; no second job is launched")
    if not work.resolve().is_relative_to(BASE.resolve()) or not work.name.startswith("imqrp_"):
        raise ValueError("work must be under imq_abaqus/recovery_presentation with an imqrp_ prefix")
    executable = shutil.which("abaqus")
    compiler = shutil.which("ifort")
    if not executable or not compiler:
        raise ValueError("reference needs abaqus and ifort on PATH; source Intel compiler env/vars.sh first")
    work.mkdir(parents=True)
    developer = work / "developer"
    developer.mkdir()
    source = ROOT.parent / "imq-umat-recovery/parameter_sensitivity/models/m3_j2/umat.for"
    built = build_provider(source.with_name("contract_v2.json"), developer / "provider")
    job_directory = developer / JOB
    job_directory.mkdir()
    deck = ROOT / "examples/presentation_request/Analysis.inp"
    command = [executable, "job=" + JOB, "input=" + str(deck), "user=" + str(source),
               "cpus=1", "interactive"]
    completed = subprocess.run(command, cwd=job_directory, capture_output=True, text=True)
    (developer / "abaqus.log").write_text(completed.stdout + completed.stderr)
    manifest = {"command": command, "returncode": completed.returncode, "ifort": compiler,
                "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(), "provider": built}
    (developer / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(completed.stdout + completed.stderr, flush=True)
    status = job_directory / (JOB + ".sta")
    if completed.returncode or not status.is_file() or "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" not in status.read_text():
        raise ValueError("reference job did not complete; see developer/abaqus.log and job .msg/.dat/.sta")
    collaborator = work / "collaborator"
    collaborator.mkdir()
    for source_path, filename in ((deck, "Analysis.inp"), (job_directory / (JOB + ".odb"), "Analysis.odb"),
                                  (Path(built["object"]), "OTI_UMAT.obj"), (Path(built["contract"]), "Mapping.json"),
                                  (deck.with_name("sensitivity_request.json"), "sensitivity_request.json")):
        shutil.copyfile(source_path, collaborator / filename)
    print("Prepared artifact-only collaborator directory:", collaborator)


def consume(work, output):
    collaborator = work / "collaborator"
    code = '''import json, os, sys
from pathlib import Path
def deny_private_source(event, args):
    if event == "open" and isinstance(args[0], (str, bytes)):
        path = os.fsdecode(args[0])
        generated_shim = path.endswith("/private/link/path_shim.for")
        if (path.lower().endswith(".for") and not generated_shim) or "/parameter_sensitivity/models/" in path:
            raise PermissionError("private material source unavailable in collaborator invocation")
sys.addaudithook(deny_private_source)
try:
    open("unavailable_private_source.for")
except PermissionError:
    print("PRIVATE SOURCE READ DENIAL ACTIVE", flush=True)
else:
    raise AssertionError("private-source denial not active")
print("COLLABORATOR CWD", os.getcwd(), flush=True)
print("COLLABORATOR INPUTS", sorted(path.name for path in Path.cwd().iterdir() if path.is_file()), flush=True)
from residual_core.ui.cli import main
sys.exit(main(sys.argv[1:]))
'''
    command = [sys.executable, "-c", code, "request", "--model", "Analysis.inp", "--odb", "Analysis.odb",
               "--material", "OTI_UMAT.obj", "--request", "sensitivity_request.json", "--out", output]
    environment = dict(os.environ, PYTHONPATH=str(ROOT) + ":" + str(ROOT.parent / "imq-umat-recovery/src"))
    completed = subprocess.run(command, cwd=collaborator, env=environment, capture_output=True, text=True)
    (work / "collaborator_invocation.log").write_text(completed.stdout + completed.stderr)
    (work / "collaborator_command.json").write_text(json.dumps(command, indent=2) + "\n")
    print(completed.stdout + completed.stderr, flush=True)
    if completed.returncode:
        raise ValueError("collaborator request failed with exit %s" % completed.returncode)
    report = json.loads((collaborator / output / "sensitivity_results.json").read_text())
    displacement, reaction, stress, state = report["results"]
    expected = {"E": -300 / 210000 ** 2, "nu": 0., "SIGY0": -1 / 2000, "H": -50 / 2000 ** 2}
    errors = {name: abs(displacement["derivatives"][name] - value) / abs(value)
              for name, value in expected.items() if value != 0}
    zero_errors = {name: abs(displacement["derivatives"][name])
                   for name, value in expected.items() if value == 0}
    checks = {"source_read_denial": True, "production_rerun_during_consume": False,
              "independent_reference": "uniaxial J2: U1=F/E+(F-SIGY0)/H", "scaled_derivative_errors": errors,
              "zero_reference_absolute_errors": zero_errors,
              "zero_reference_absolute_tolerance": 1e-8,
              "max_equilibrium_error": report["metadata"]["max_equilibrium_error"],
              "reaction_value_error": abs(reaction["value"] + 300),
              "stress_value_error": abs(stress["value"] - 300),
              "state_value_error": abs(state["value"] - .025)}
    checks["passed"] = (max(errors.values()) < 2e-5 and max(zero_errors.values()) < 1e-8
                        and checks["reaction_value_error"] < 1e-3
                        and checks["stress_value_error"] < 1e-3 and checks["state_value_error"] < 1e-6)
    (work / "analytic_checks.json").write_text(json.dumps(checks, indent=2) + "\n")
    print(json.dumps(checks, indent=2))
    if not checks["passed"]:
        raise ValueError("independent analytic proof failed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=["prepare", "consume"])
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--out", default="results")
    args = parser.parse_args()
    try:
        if args.phase == "prepare":
            prepare(args.work.resolve())
        else:
            consume(args.work.resolve(), args.out)
        return 0
    except (OSError, ValueError, RuntimeError) as error:
        print(str(error), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())