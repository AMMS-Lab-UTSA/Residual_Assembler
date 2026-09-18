"""Capture bounded recovery help/examples/GUI evidence; never launch an Abaqus job."""

import argparse
from collections import Counter
import fnmatch
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
import xml.etree.ElementTree as ET


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--umat", type=Path, required=True)
    parser.add_argument("--phase", choices=("help", "examples", "gui", "presentation", "ledger", "docs"), required=True)
    parser.add_argument("--work", type=Path)
    parser.add_argument("--imports", choices=("installed", "environment"), default="installed")
    parser.add_argument("--evidence-dir", type=Path, help="write separate evidence instead of replacing historical usage reports")
    args = parser.parse_args(argv)
    ra = Path(__file__).resolve().parents[1]
    umat = args.umat.resolve()
    work = args.work.resolve() if args.work else Path(tempfile.mkdtemp(prefix="recovery_usage_"))
    work.mkdir(parents=True, exist_ok=True)
    environment = dict(os.environ)
    environment["PATH"] = str(Path(sys.executable).parent) + os.pathsep + environment["PATH"]
    if args.imports == "installed":
        environment.pop("PYTHONPATH", None)
        environment.pop("PYTHONHOME", None)
    environment["UMAT_OTI_REPO"] = str(umat)
    environment["RUN_OTILIB_TESTS"] = "1"
    records = {"RA": [], "UMAT": []}
    roots = {"RA": ra, "UMAT": umat}
    previous_failures = {}
    for name, root in roots.items():
        previous = root / "docs/evidence" / ("usage_" + args.phase + ".json")
        data = json.loads(previous.read_text()) if previous.exists() else {}
        previous_failures[name] = data.get("previous_failures", []) + [
            record for record in data.get("commands", []) if record.get("returncode", 0)]

    def save():
        for name, root in roots.items():
            destination = root / "docs/evidence" / ("usage_" + args.phase + ".json")
            if args.evidence_dir:
                destination = args.evidence_dir.resolve() / name / destination.name
                destination.parent.mkdir(parents=True, exist_ok=True)
            payload = {"clean_install_verified": False, "phase": args.phase,
                       "python": sys.executable, "work": str(work),
                       "imports": args.imports,
                       "environment": {key: environment.get(key) for key in
                                       ("PYTHONPATH", "UMAT_OTI_REPO", "PYOTI_PATH", "OTILIB_ROOT", "RUN_OTILIB_TESTS")},
                       "commands": records[name], "previous_failures": previous_failures[name]}
            destination.write_text(json.dumps(payload, indent=2) + "\n")

    def run(repo, name, arguments, cwd=None):
        command = [str(value) for value in arguments]
        if args.imports == "installed" and command[0] == sys.executable:
            command.insert(1, "-I")
        started = time.monotonic()
        completed = subprocess.run(command, cwd=cwd or roots[repo], env=environment,
                                   capture_output=True, text=True)
        record = {"name": name, "argv": command, "cwd": str(cwd or roots[repo]),
                  "returncode": completed.returncode, "stdout": completed.stdout,
                  "stderr": completed.stderr, "seconds": time.monotonic() - started}
        records[repo].append(record)
        save()
        print(repo, name, "exit", completed.returncode, flush=True)
        if completed.returncode:
            raise RuntimeError(completed.stdout + completed.stderr)
        return record

    def proof(repo, record, path):
        content = path.read_bytes()
        record["proof"] = {"path": str(path), "sha256": hashlib.sha256(content).hexdigest(),
                           "data": json.loads(content)}
        save()

    python = sys.executable
    cli = [python, "-m", "residual_core.ui.cli"]
    run("RA", "environment", [python, "-c", "from residual_core.algebra.otilib_adapter import OtiContext; OtiContext(1,1); "
        "import sys,residual_core,umat_oti,pyoti.sparse; "
        "from umat_oti.store import transform_fingerprint; "
        "print(sys.version); print(residual_core.__file__); print(umat_oti.__file__); "
        "print(pyoti.sparse.__file__); print(transform_fingerprint())"])
    run("UMAT", "compiler", ["gfortran", "--version"])

    if args.phase == "help":
        for command in ([], *[[name] for name in
                ("backends", "modes", "inspect", "inspect-model", "requirements", "assemble",
                 "verify", "sensitivity", "init", "init-assembly", "check", "run", "report",
                 "doctor", "template", "replay", "request")]):
            run("RA", "resasm " + " ".join(command) + " --help", cli + command + ["--help"])
        for module in ("umat_oti.cli", "umat_oti.cli_json", "umat_oti.cli_batch",
                       "umat_oti.pipeline.cli", "umat_oti.provider", "umat_oti.reproduce",
                       "umat_oti.validation.parameter_sensitivity_provider"):
            run("UMAT", module + " --help", [python, "-m", module, "--help"])
        for subcommand in ("transform", "config"):
            run("UMAT", "umat-oti " + subcommand + " --help", [python, "-m", "umat_oti.cli", subcommand, "--help"])
        run("UMAT", "provider build --help", [python, "-m", "umat_oti.provider", "build", "--help"])
        run("UMAT", "internal Jacobian --help", [python, "examples/verify_internal_jacobian.py", "--help"])
        for script in ("scripts/reproduce_imqcam_pipeline.py", "scripts/reproduce_presentation_request.py",
                       "scripts/check_presentation_browser.py", "scripts/clean_install_gate.py",
                       "examples/finite_strain_c3d8/benchmark.py"):
            run("RA", script + " --help", [python, script, "--help"])
        run("UMAT", "clean install --help", [python, "scripts/clean_install_gate.py", "--ra-repo", ra, "--help"])
        run("RA", "streamlit --help", [python, "-m", "streamlit", "run", "--help"])
        run("UMAT", "pipeline stages", [python, "-m", "umat_oti.pipeline.cli", "--list-stages"])
        run("RA", "backends", cli + ["backends"])
        run("RA", "modes", cli + ["modes"])

    elif args.phase == "examples":
        import numpy as np

        for template, label, expected in (
            ("user_python_residual", "R-X1", {(1, 0): -1/3, (0, 1): 1/24,
               (2, 0): 2/9, (1, 1): -1/144, (0, 2): -1/576}),
            ("user_blackbox_order2_residual", "R-X2", {(1, 0): -2/3, (0, 1): 1/48,
               (2, 0): 5/9, (1, 1): -1/144, (0, 2): -1/2304}),
        ):
            destination = work / label
            template_name = "python" if template == "user_python_residual" else "blackbox-order2"
            run("RA", label + " init", cli + ["init", "--template", template_name, "--out", destination])
            run("RA", label + " check", cli + ["check", destination / "resasm.yml"])
            record = run("RA", label, cli + ["run", destination / "resasm.yml"])
            private = destination / "resasm_output/private"
            checks = []
            for order in (1, 2):
                mapping = json.loads((private / f"direction_map_order{order}.json").read_text())
                values = np.load(private / f"solution_sensitivities_order{order}.npz")
                for column in mapping["columns"]:
                    reference = expected[tuple(column["exponents"])]
                    measured = float(values["U_derivatives"][0, column["index"]])
                    checks.append({"direction": column["exponents"], "measured": measured,
                                   "reference": reference, "absolute_error": abs(measured-reference)})
            record["verification"] = {"reference": "analytic cubic solution", "tolerance": 1e-8,
                                      "checks": checks, "passed": all(item["absolute_error"] < 1e-8 for item in checks)}
            save()
            assert record["verification"]["passed"]
            run("RA", label + " report", cli + ["report", destination / "resasm_output"])
        model = ra / "residual_core/examples/minimal_c3d8_stress_driven"
        record = run("RA", "R-X3", cli + ["assemble", model / "model.json", "--mode", "stress-driven",
                     "--fields", model / "fields.json", "--out", work / "stress.npy"])
        residual = np.load(work / "stress.npy")
        record["verification"] = {"norm": float(np.linalg.norm(residual)),
            "max_abs": float(np.max(np.abs(residual))), "reference_norm": 100 / np.sqrt(2),
            "reference_max_abs": 25.0, "tolerance": 1e-10,
            "passed": bool(np.allclose(np.linalg.norm(residual), 100 / np.sqrt(2), atol=1e-10, rtol=0)
                           and np.allclose(np.max(np.abs(residual)), 25, atol=1e-10, rtol=0))}
        save()
        assert record["verification"]["passed"]
        record = run("RA", "R-X4", [python, "scripts/reproduce_imqcam_pipeline.py", "--skip-abaqus",
                 "--provider-repo", umat, "--imports", args.imports, "--out", work / "R-X4"])
        proof("RA", record, work / "R-X4/private/manifest.json")
        assert record["proof"]["data"]["passed"]
        record = run("RA", "R-X5", [python, "examples/finite_strain_c3d8/benchmark.py", "--out", work / "R-X5"])
        proof("RA", record, work / "R-X5/report.json")
        assert record["proof"]["data"]["passed"]
        run("RA", "R-X5 public assembly", cli + ["--config", work / "R-X5/config.json", "assemble",
            work / "R-X5/model.json", "--mode", "material-replay", "--tangent"])
        run("RA", "R-X5 public sensitivity", cli + ["--config", work / "R-X5/config.json", "sensitivity",
            work / "R-X5/model.json", "--params", work / "R-X5/params.json", "--out", work / "R-X5/sensitivity"])
        for label, model_name, elastic in (("U-X1", "m1_elastic", True), ("U-X2", "m3_j2", False),
                                           ("U-X3", "m3_j2", False), ("U-X4", "m3_j2", False)):
            command = [python, "-m", "umat_oti.validation.parameter_sensitivity_provider",
                       f"parameter_sensitivity/models/{model_name}/contract_v2.json", "--out", work / label]
            if elastic:
                command.append("--elastic")
            record = run("UMAT", label, command)
            proof("UMAT", record, work / label / "verification.json")
            assert record["proof"]["data"]["passed"]
        record = run("UMAT", "U-X5", [python, "examples/verify_internal_jacobian.py", "--out", work / "U-X5"])
        proof("UMAT", record, work / "U-X5/verification.json")
        assert record["proof"]["data"]["passed"]

    elif args.phase == "gui":
        code = "from streamlit.testing.v1 import AppTest; import sys; app=AppTest.from_file(sys.argv[1]).run(timeout=90); assert not app.exception, list(app.exception); print('rendered', len(app.title), 'titles', len(app.tabs), 'tabs')"
        for repo, module in (("RA", "residual_core/app/streamlit_app.py"),
                             ("UMAT", "src/umat_oti/app/streamlit_app.py")):
            run(repo, "primary GUI render", [python, "-c", code, roots[repo] / module])
        for script in ("src/umat_oti/app/unified_app.py", "src/umat_oti/app/workbench_app.py"):
            if (umat / script).exists():
                run("UMAT", script + " render", [python, "-c", code, umat / script])
        server_code = (
            "import json,os,pathlib,sys; sys.path.insert(0,sys.argv[1]); "
            "from clean_install_gate import Gate,check_servers; "
            "work=pathlib.Path(sys.argv[2]); work.mkdir(); gate=Gate(work); "
            "gate.environment.update(os.environ); gate.report['kind']='source-tree-gui-health'; "
            "check_servers(gate,sys.executable,{'primary':{'path':sys.argv[3]}}); "
            "print(json.dumps(gate.report['servers']))"
        )
        for repo, root in roots.items():
            run(repo, "primary GUI HTTP launch", [python, "-c", server_code, ra / "scripts",
                work / (repo + "_server"), root / "scripts/app.py"])

    elif args.phase == "presentation":
        source = ra.parent / "imq_abaqus/recovery_presentation/imqrp_reference/collaborator"
        destination = work / "presentation/collaborator"
        destination.mkdir(parents=True)
        hashes = {}
        for name in ("Analysis.inp", "Analysis.odb", "OTI_UMAT.obj", "Mapping.json", "sensitivity_request.json"):
            shutil.copyfile(source / name, destination / name)
            hashes[name] = hashlib.sha256((destination / name).read_bytes()).hexdigest()
        record = run("RA", "source-denied presentation", [python, "scripts/reproduce_presentation_request.py",
                     "consume", "--work", destination.parent, "--out", "results"])
        proof("RA", record, destination.parent / "analytic_checks.json")
        record["inputs_sha256"] = hashes
        record["public_files"] = sorted(path.name for path in (destination / "results").iterdir() if path.is_file())
        record["public_result"] = json.loads((destination / "results/sensitivity_results.json").read_text())
        assert record["public_files"] == ["run_report.txt", "sensitivity_results.json", "sensitivity_tables.csv"]
        assert record["proof"]["data"]["passed"]
        records["UMAT"].append(record)
        save()

    elif args.phase == "ledger":
        groups = json.loads((ra / "docs/evidence/usage_coverage.json").read_text())["groups"]
        statuses = {"implemented": "IMPLEMENTED - not yet reproduced from clean install",
                    "partial": "IN PROGRESS (recovery audit: bounded evidence)",
                    "unestablished": "NOT STARTED"}
        for repo, root in roots.items():
            path = root / "docs/COMPLETION_LEDGER.md"
            lines = path.read_text().splitlines()
            rows = []
            for line in lines:
                cells = [cell.strip() for cell in line.split("|")]
                if len(cells) == 9 and cells[1].startswith(("U-", "R-", "X-", "T-", "CI-")):
                    rows.append(cells)
            assert len(rows) == 274
            normalized = []
            for number, group in enumerate(groups):
                item = dict(group, id=number, clean_install_verified=False)
                for field in ("implementation", "tests", "proof"):
                    resolved = [roots[group["repo"]] / value for value in group[field]]
                    assert all(value.exists() for value in resolved), resolved
                    item[field] = [os.path.relpath(value, root) for value in resolved]
                item["cwd"] = os.path.relpath(roots[group["repo"]], root)
                junit = roots[group["repo"]] / "docs/evidence/usage_focused_final.xml"
                if junit.exists():
                    modules = {Path(value).stem for value in group["tests"]}
                    cases = [case for case in ET.parse(junit).findall(".//testcase")
                             if any(module in case.get("classname", "").split(".") for module in modules)]
                    if cases:
                        item["fresh_test_evidence"] = {
                            "junit": os.path.relpath(junit, root), "executed": len(cases),
                            "failed": sum(case.find("failure") is not None for case in cases),
                            "errors": sum(case.find("error") is not None for case in cases),
                            "skipped": sum(case.find("skipped") is not None for case in cases),
                        }
                normalized.append(item)
            entries = []
            rendered = []
            for cells in rows:
                matches = [group for group in normalized if any(fnmatch.fnmatchcase(cells[1], pattern)
                                                                for pattern in group["ids"])]
                group = matches[-1] if matches else None
                status = statuses[group["status"]] if group else "NOT STARTED"
                entry = {"id": cells[1], "requirement": cells[2], "repository": cells[3],
                         "evidence_group": group["id"] if group else None,
                         "status": status, "clean_install_verified": False}
                entries.append(entry)
                if group:
                    def links(field):
                        return "; ".join(f"[{value}](../{value})" for value in group[field])
                    cells[4] = links("implementation")
                    cells[5] = links("tests") + "; [bounded scope/proof](evidence/usage_requirements.json)"
                    cells[6] = f"cwd `{group['cwd']}`: `{group['command']}`"
                else:
                    cells[4:7] = ["Not established by this audit", "No row-specific proof", "No verified reproduction"]
                cells[7] = status
                rendered.append("| " + " | ".join(cells[1:8]) + " |")
            counts = dict(Counter(entry["status"] for entry in entries))
            payload = {"clean_install_verified": False, "requirements": 274,
                       "completed": 0, "outstanding": 274, "status_counts": counts,
                       "groups": normalized, "rows": entries,
                       "unmapped": [entry["id"] for entry in entries if entry["evidence_group"] is None]}
            (root / "docs/evidence/usage_requirements.json").write_text(json.dumps(payload, indent=2) + "\n")
            header_end = next(index for index, line in enumerate(lines) if line.startswith("| ID |"))
            header = lines[:header_end]
            marker = "## Current Recovery Audit"
            if marker in header:
                header = header[:header.index(marker)]
            header.extend([marker, "", "The rows below describe recovery working-tree evidence, not final-branch completion.",
                "Every linked proof has a bounded scope in [the executable index](evidence/usage_requirements.json).",
                "NOT STARTED means not established against this requirement, not proof that code is absent.",
                "Earlier wheel isolation is not a final clean clone. No row is clean-install PASS.", "",
                f"Counts: 274 requirements; 0 completed; 274 outstanding. {counts}", "",
                "| ID | Requirement | Repository | Implementation | Test | Reproduction command | Status |",
                "| -- | ----------- | ---------- | -------------- | ---- | -------------------- | ------ |"])
            path.write_text("\n".join(header + rendered + ["", "Rows: 274", ""]))
            records[repo].append({"name": "ledger rendering", "requirements": 274,
                                  "completed": 0, "outstanding": 274, "status_counts": counts})
            print(repo, counts, "unmapped", payload["unmapped"], flush=True)
            save()

    else:
        spec = importlib.util.spec_from_file_location("command_audit", umat / "tools/audit_documentation_commands.py")
        audit = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(audit)
        for repo, root in roots.items():
            audit.REPO_ROOT = root
            audit.doc_files = lambda root=root: [root / "docs/USAGE_REPORT.md"]
            problems = audit.audit()
            records[repo].append({"name": "usage report links and commands", "problems": problems})
            save()
            assert not problems, problems
    save()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())