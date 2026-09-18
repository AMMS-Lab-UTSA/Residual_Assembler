"""Rebuild Program 1, run Program 2, and independently verify a bounded J2 slice."""

from __future__ import annotations

import argparse
import csv
import hashlib
from importlib import metadata
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def repository_info(path, patterns):
    commit = subprocess.check_output(["git", "-C", str(path), "rev-parse", "HEAD"], text=True).strip()
    status = subprocess.check_output(["git", "-C", str(path), "status", "--short"], text=True)
    sources = {}
    for pattern in patterns:
        for source in sorted(path.glob(pattern)):
            if source.is_file():
                sources[str(source.relative_to(path))] = digest(source)
    fingerprint = hashlib.sha256(json.dumps(sources, sort_keys=True).encode()).hexdigest()
    return {"path": str(path), "commit": commit, "status": status,
            "source_tree_sha256": fingerprint, "source_sha256": sources}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True, help="fresh output directory; existing paths are refused")
    parser.add_argument("--provider-repo", type=Path, default=ROOT.parent / "imq-umat-recovery")
    parser.add_argument("--model", type=Path, default=ROOT / "examples/imqcam_j2_cantilever/model.json")
    parser.add_argument("--skip-abaqus", action="store_true", help="use genuine archived elastic Abaqus exports; J2 is a synthetic converged FE solve")
    args = parser.parse_args(argv)
    output = args.out.resolve()
    if output.exists():
        parser.error("--out must be a fresh directory")
    output.mkdir(parents=True)
    private = output / "private"
    public = output / "public"
    private.mkdir()
    public.mkdir()
    manifest = {"passed": False, "scope": "one C3D8, eight IPs, seven cyclic J2 increments; first-order material parameters",
                "commands": [], "abaqus_available": shutil.which("abaqus"),
                "abaqus_execution": "not_run", "j2_primal_origin": "synthetic_converged_fe"}

    def command(arguments, environment):
        entry = {"argv": [str(value) for value in arguments], "cwd": str(ROOT)}
        manifest["commands"].append(entry)
        completed = subprocess.run(entry["argv"], cwd=ROOT, env=environment,
                                   capture_output=True, text=True)
        entry["returncode"] = completed.returncode
        log = private / f"command-{len(manifest['commands'])}.log"
        log.write_text(completed.stdout + completed.stderr)
        entry["log"] = str(log.relative_to(output))
        if completed.returncode:
            raise RuntimeError(f"command failed ({completed.returncode}): {entry['argv']}\n{completed.stdout}\n{completed.stderr}")
        return completed.stdout

    try:
        if not args.skip_abaqus:
            raise ValueError("This bounded reproducer does not launch Abaqus. Pass --skip-abaqus to verify the genuine archived elastic export; it does not claim J2 Abaqus validation.")
        provider = args.provider_repo.resolve(strict=True)
        model_path = args.model.resolve(strict=True)
        if ROOT != Path.cwd().resolve():
            raise ValueError(f"run with explicit recovery cwd: {ROOT}")
        if shutil.which("gfortran") is None:
            raise ValueError("gfortran is required; no prebuilt-object fallback")
        environment = os.environ.copy()
        environment["PYTHONPATH"] = os.pathsep.join([str(ROOT), str(provider / "src"),
                                                      environment.get("PYTHONPATH", "")])
        environment["UMAT_OTI_REPO"] = str(provider)
        manifest["environment"] = {key: environment.get(key) for key in
                       ("PYTHONPATH", "UMAT_OTI_REPO", "PYOTI_PATH", "OTILIB_ROOT", "RUN_OTILIB_TESTS", "FC")}
        manifest["python"] = {"executable": sys.executable, "version": sys.version, "platform": platform.platform()}
        manifest["dependencies"] = {name: metadata.version(name) for name in ("numpy", "matplotlib", "pytest", "streamlit")}
        manifest["compiler"] = command(["gfortran", "--version"], environment).splitlines()[0]
        manifest["repositories"] = {
            "residual_assembler": repository_info(ROOT, ["residual_core/**/*.py", "scripts/reproduce_imqcam_pipeline.py",
                                                         "examples/imqcam_j2_cantilever/*", "tests/abaqus_derivative_export/**/*"]),
            "umat_oti": repository_info(provider, ["src/umat_oti/**/*.py", "parameter_sensitivity/models/m3_j2/*"])}
        manifest["model"] = {"path": str(model_path), "sha256": digest(model_path)}
        provider_contract = provider / "parameter_sensitivity/models/m3_j2/contract_v2.json"
        built = json.loads(command([sys.executable, "-m", "umat_oti.provider", "build", provider_contract,
                                    "--out", private / "provider"], environment))
        contract = json.loads(Path(built["contract"]).read_text())
        manifest["provider"] = {"object_sha256": digest(built["object"]),
                                "contract_sha256": digest(built["contract"]),
                                "original_source_sha256": digest(provider_contract.parent / "umat.for"),
                                "abi": contract["symbols"], "march": contract["march"]}
        job = private / "solve"
        common = ["--object", built["object"], "--contract", built["contract"]]
        command([sys.executable, "-m", "residual_core.ui.cli", "replay", model_path, *common,
                 "--solve", "--verify", "--out", job], environment)
        command([sys.executable, "-m", "residual_core.ui.cli", "replay", job / "private/record.json",
                 *common, "--out", private / "replay"], environment)
        sys.path.insert(0, str(ROOT))
        from residual_core.replay.verification import verify_abaqus_fixture
        fixture_report = verify_abaqus_fixture(ROOT / "tests/abaqus_derivative_export")
        manifest["abaqus_fixture"] = fixture_report
        summary = json.loads((job / "public/summary.json").read_text())
        result = json.loads((job / "private/result.json").read_text())
        (public / "verification.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
        (public / "abaqus_fixture.json").write_text(json.dumps(fixture_report, indent=2, allow_nan=False) + "\n")
        rows = summary["verification"]["rows"]
        columns = list(dict.fromkeys(key for row in rows for key in row))
        with (public / "fd_table.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=columns)
            writer.writeheader()
            writer.writerows(rows)
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        displacement = np.array([row["u"] for row in result["increments"]])
        sensitivity = np.array([row["du_dp"] for row in result["increments"]])
        naive = np.array([row["fixed_path_local_solve"] for row in result["increments"]])
        figure, axes = plt.subplots(1, 2, figsize=(10, 4))
        increments = np.arange(1, len(displacement) + 1)
        axes[0].plot(increments, displacement[:, -1], "o-")
        axes[0].set(xlabel="Increment", ylabel="Last node UZ", title="Synthetic converged J2 FE")
        axes[1].plot(increments, sensitivity[:, -1, 1], "o-", label="Total-history du/dnu")
        axes[1].plot(increments, naive[:, -1, 1], "x--", label="Fixed-path local solve (not du/dnu)")
        axes[1].set(xlabel="Increment", ylabel="Last node derivative")
        axes[1].legend(fontsize=8)
        figure.tight_layout()
        figure.savefig(public / "j2_history.png", dpi=160)
        plt.close(figure)
        manifest["verification"] = summary["verification"]
        manifest["max_equilibrium_error"] = summary["max_equilibrium_error"]
        manifest["passed"] = True
        manifest["outputs"] = {str(path.relative_to(output)): digest(path)
                               for path in public.iterdir() if path.is_file()}
    except (OSError, ValueError, RuntimeError, ImportError, AssertionError, subprocess.CalledProcessError) as error:
        manifest["error"] = str(error)
        print(f"reproduction failed: {error}", file=sys.stderr)
    finally:
        (private / "manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
        (output / "manifest.json").write_text(json.dumps({
            "passed": manifest["passed"], "private_manifest": "private/manifest.json",
            "public_results": "public", "error": manifest.get("error")}, indent=2) + "\n")
    if manifest["passed"]:
        print(f"verified bounded J2 pipeline: {output / 'manifest.json'}")
    return 0 if manifest["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())