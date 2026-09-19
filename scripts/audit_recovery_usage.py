"""Run the documented commands of one phase and record each one; never launch an Abaqus job.

``--phase examples`` is the one-command check of the worked examples that need
no Abaqus: Residual_Assembler Examples 1, 2, 4, 6, 7 and 8 and UMAT-OTI
Examples 1 to 6, each with the commands its walkthrough gives, after the
documented OTILib and transform-generation checks. It runs from the
Residual_Assembler checkout, with the UMAT-OTI checkout given by ``--umat``
(Residual_Assembler docs/INSTALL.md, section 6, item 7).
"""

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
    # the checkout the connected pipeline (Example 6) builds its provider from
    environment["UMAT_OTI_REPO"] = str(umat)
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

    def expect(record, text):
        # a documented check prints its verdict: the record passes only when it does
        record["verification"] = {"expected": text, "passed": text in record["stdout"]}
        save()
        assert record["verification"]["passed"], (record["name"], text)

    python = sys.executable
    # the console commands of the environment running this script (its bin folder is first on PATH)
    cli = ["resasm"]
    # the documented checks: OTILib (Residual_Assembler docs/INSTALL.md section 4) and the
    # transform generation (both usage reports)
    expect(run("RA", "environment", [python, "-c",
           "from residual_core.algebra.otilib_adapter import otilib_status; print(otilib_status())"]),
           "'available': True")
    expect(run("UMAT", "fingerprint", [python, "-c",
           "from umat_oti.store import transform_fingerprint; print(transform_fingerprint())"]),
           json.loads((ra / "schemas/transform_generation.json").read_text())["transform_fingerprint"])
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
        for script in ("scripts/reproduce_connected_pipeline.py", "scripts/reproduce_presentation_request.py",
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

        # Residual_Assembler (cwd: its root). R-X1..R-X5 are the ten reference examples'
        # labels; the documented example each one runs is named in the comment.
        for template, label, expected in (
            # Example 1 (examples/user_config_minimal/WALKTHROUGH.md)
            ("python", "R-X1", {(1, 0): -1/3, (0, 1): 1/24,
               (2, 0): 2/9, (1, 1): -1/144, (0, 2): -1/576}),
            # Example 8 (templates/user_blackbox_order2_residual/WALKTHROUGH.md)
            ("blackbox-order2", "R-X2", {(1, 0): -2/3, (0, 1): 1/48,
               (2, 0): 5/9, (1, 1): -1/144, (0, 2): -1/2304}),
        ):
            destination = work / label
            run("RA", label + " init", cli + ["init", "--template", template, "--out", destination])
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
        # Example 2 (residual_core/examples/minimal_c3d8_stress_driven/WALKTHROUGH.md)
        model = "residual_core/examples/minimal_c3d8_stress_driven"
        record = run("RA", "R-X3", cli + ["assemble", model + "/model.json", "--mode", "stress-driven",
                     "--fields", model + "/fields.json", "--out", work / "R_cube.npy"])
        residual = np.load(work / "R_cube.npy")
        record["verification"] = {"norm": float(np.linalg.norm(residual)),
            "max_abs": float(np.max(np.abs(residual))), "reference_norm": 100 / np.sqrt(2),
            "reference_max_abs": 25.0, "tolerance": 1e-10,
            "passed": bool(np.allclose(np.linalg.norm(residual), 100 / np.sqrt(2), atol=1e-10, rtol=0)
                           and np.allclose(np.max(np.abs(residual)), 25, atol=1e-10, rtol=0))}
        save()
        assert record["verification"]["passed"]
        # Example 4 (examples/replay_history/WALKTHROUGH.md; docs/INSTALL.md section 6, item 2)
        run("RA", "Example 4 provider", ["umat-oti-provider", "build",
            umat / "parameter_sensitivity/models/m3_j2/contract_v2.json", "--out", work / "provider_j2"])
        beam = "examples/replay_history/j2_beam"
        record = run("RA", "Example 4", cli + ["history", "--model", beam + "/Analysis.inp",
                     "--fields", beam + "/fields.npz", "--material", work / "provider_j2/umat_m3_j2_oti.obj",
                     "--request", beam + "/sensitivity_request.json", "--out", work / "beam"])
        record["verification"] = {"expected": "run_report.txt begins with 'Status: executed successfully'",
                                  "passed": (work / "beam/run_report.txt").read_text().startswith(
                                      "Status: executed successfully")}
        save()
        assert record["verification"]["passed"]
        # Example 6 (examples/bounded_j2_c3d8/WALKTHROUGH.md)
        record = run("RA", "R-X4", [python, "scripts/reproduce_connected_pipeline.py", "--skip-abaqus",
                     "--out", work / "R-X4"])
        proof("RA", record, work / "R-X4/private/manifest.json")
        assert record["proof"]["data"]["passed"]
        # Example 7 (examples/finite_strain_c3d8/WALKTHROUGH.md)
        record = run("RA", "R-X5", [python, "examples/finite_strain_c3d8/benchmark.py", "--out", work / "R-X5"])
        proof("RA", record, work / "R-X5/report.json")
        assert record["proof"]["data"]["passed"]
        run("RA", "R-X5 public assembly", cli + ["--config", work / "R-X5/config.json", "assemble",
            work / "R-X5/model.json", "--mode", "material-replay", "--tangent"])
        run("RA", "R-X5 public sensitivity", cli + ["--config", work / "R-X5/config.json", "sensitivity",
            work / "R-X5/model.json", "--params", work / "R-X5/params.json", "--out", work / "R-X5/sensitivity"])

        # UMAT-OTI (cwd: its root), the command blocks of examples/0N_*/README.md with the
        # outputs in the work folder instead of umat_oti_workspace/examples/.
        # Examples 1 and 2: the tangent from four fields, then the comparison script
        for label, model_name, script in (
                ("U-X1", "m1_elastic", "examples/01_elastic_tangent/run.py"),
                ("U-X2", "m3_j2", "examples/02_j2_plasticity_tangent/run.py")):
            run("UMAT", label + " jacobian", ["umat-oti", "jacobian",
                f"parameter_sensitivity/models/{model_name}/umat.for", "--ntens", "6",
                "--out", work / label, "--compile"])
            expect(run("UMAT", label, [python, script, "--jacobian-dir", work / label]), "RESULT: PASS")
        # Examples 3 and 4: the provider, the verified hand-off package, the package reader
        for label, contract, options, script in (
                ("U-X3", "parameter_sensitivity/models/m3_j2/contract_v2.json", ["--j2-branches"],
                 "examples/03_j2_parameter_sensitivities/run.py"),
                ("U-X4", "examples/04_fcc_crystal_plasticity_provider/contract_tension_shear.json", [],
                 "examples/04_fcc_crystal_plasticity_provider/run.py")):
            run("UMAT", label + " build", ["umat-oti-provider", "build", contract,
                "--out", work / label / "build", "--regular-object", "REAL_UMAT.obj"])
            record = run("UMAT", label + " package", [python, "-m", "umat_oti.provider.collaborator", contract,
                         "--out", work / label / "package", *options])
            proof("UMAT", record, work / label / "package/verification/verification.json")
            assert record["proof"]["data"]["passed"]
            expect(run("UMAT", label, [python, script, "--package", work / label / "package"]), "RESULT: PASS")
        # Example 5: Part A (the bundled flow model), then Part B (the damage UMAT)
        record = run("UMAT", "U-X5 part A", [python, "examples/verify_internal_jacobian.py",
                     "--out", work / "U-X5/05_cpflow"])
        proof("UMAT", record, work / "U-X5/05_cpflow/verification.json")
        assert record["proof"]["data"]["passed"]
        expect(run("UMAT", "U-X5", [python, "examples/05_internal_newton_jacobian/run.py",
                                    "--out", work / "U-X5/05_vpdco"]), "RESULT: PASS")
        # Example 6: the twenty-model sweep (absolute directories, as its README says)
        record = run("UMAT", "Example 6", [python, "tools/run_parameter_sensitivity_sweep.py",
                     "--work-dir", work / "06_sweep/work", "--results-dir", work / "06_sweep/results"])
        sweep = json.loads((work / "06_sweep/results/parameter_sensitivity_round.json").read_text())
        disagreeing = sum((model["stages"].get("derivatives_verified") or {}).get("rows_disagreeing", 0)
                          for model in sweep["models"])
        # every model reproduces the original's stress and no comparison row disagrees
        # (a row the reference cannot resolve is reported as unresolved, not as agreement)
        record["verification"] = {"funnel": sweep["funnel"], "rows_disagreeing": disagreeing,
                                  "passed": sweep["funnel"]["primal_parity"] == sweep["funnel"]["attempted"]
                                  and disagreeing == 0}
        save()
        assert record["verification"]["passed"]

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