#!/usr/bin/env python3
"""Slide 8: the benchmark UMATs' automatic DDSDDE against the original tangent in Abaqus.

What the slide says
    "Abaqus tangent (DDSDDE) 18/18 verified. 12 reproduce the original tangent
    exactly. The other 6 differ only from notching, rounding, or an error in the
    original UMAT" with a per-UMAT "DDSDDE max error (abs / rel)" column, e.g.
    UMAT_NKH 0.72 / 2.1e-3, spin_elas_def 740 / 2.6e-3 ("error in original UMAT").

How that table was produced (recovered)
    ~/Documents/UMAT_source_transformation/validate_all_local.py transformed every
    completed benchmark contract, built the paired Abaqus validation workspace with
    DDSDDE FORCED into the compared outputs (the benchmark JSON of spin_elas_def
    does not request it), ran the original and the transformed UMAT in Abaqus
    locally, and wrote umat_oti_workspace/validate_all/summary.json; RA commit
    3521485 (results/make_transform_table.py) rendered the slide table from it
    (19 cases, UMAT_VPDCL_R failed to run and was left off the slide). The ARC
    run of the same batch (paper_results/arc_791506) matches the slide for the
    five nonzero rows it has; it did not request DDSDDE for spin_elas_def.

What this script does
    The same procedure with the current repository: for every
    ``benchmarks/*.json`` of the UMAT repository, ``run_config_transform`` then
    ``build_validation_workspace`` (the functions ``tools/run_completed_json_batch.py
    --validate`` uses) with compare_outputs STRESS, STATEV, DDSDDE, CONVERGENCE,
    then both Abaqus jobs, one at a time, named claudeP_c4_<case>_{orig,oti}
    (the outputs are renamed to the names the extractor expects afterwards),
    ``extract_results`` and ``compare_validation_results``. A job counts only if
    Abaqus wrote its completion mark; the documented 2021.HF5 teardown abort is
    recognised by the repository's own ``completed_despite_teardown_abort``.

    python presentation/claim4_benchmark_ddsdde.py --abaqus            # all cases
    python presentation/claim4_benchmark_ddsdde.py --abaqus --case UMAT_HIN
    python presentation/claim4_benchmark_ddsdde.py                     # transform + workspace only
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from presentation.common import (  # noqa: E402
    default_out, default_work, format_e, fresh_dir, provenance, require_worktree_imports, umat_repo,
    write_json,
)

SLIDE = {
    "UMAT_ECL_TEMP": ("0 / 0", "Exact"), "UMAT_ECO": ("0 / 0", "Exact"), "UMAT_HIN": ("0 / 0", "Exact"),
    "UMAT_NKH_1.02": ("0.72 / 2.1e-3", "Notching of derivative terms"), "UMAT_PCL": ("0 / 0", "Exact"),
    "UMAT_PCLI": ("0 / 0", "Exact"), "UMAT_PCLI_R": ("0 / 0", "Exact"), "UMAT_PCLK": ("0 / 0", "Exact"),
    "UMAT_PCO": ("0 / 0", "Exact"), "UMAT_VPDCL": ("0.24 / 1.2e-3", "Notching of derivative terms"),
    "UMAT_VPDCO": ("0.37 / 1.1e-3", "Notching of derivative terms"),
    "code_exp": ("0.031 / 1.6e-7", "Floating-point rounding"),
    "code_imp": ("0.016 / 8.2e-8", "Floating-point rounding"), "elastic": ("0 / 0", "Exact"),
    "spin_elas_def": ("740 / 2.6e-3", "Error in original UMAT"), "spin_elastic": ("0 / 0", "Exact"),
    "visco_beam": ("0 / 0", "Exact"), "visco_imp": ("0 / 0", "Exact"),
}
FORCED_OUTPUTS = ["STRESS", "STATEV", "DDSDDE", "CONVERGENCE"]

#: Documented input corrections, run only with --variants and reported as
#: separate rows; the benchmark rows above are always the contracts as committed.
VARIANTS = {
    "UMAT_HIN": {
        "demote": ["ONE", "TWO", "ZERO"],
        "why": ("the contract promotes ONE, TWO, ZERO, which the source initialises by DATA "
                "(KFORMC, KDLT2) and the current transformer refuses DATA-initialised promoted "
                "variables; they are the numeric constants 1, 2, 0 and carry no derivative, so the "
                "variant lists them as constants instead"),
    },
    "UMAT_PCO": {
        "closure": True,
        "why": ("UMAT_PCO.for calls helpers it does not define (KCLEAR, KMMULT, ...); the variant "
                "transforms the source's resolved routine closure (umat_oti.transform."
                "dependency_resolution, entry file first so the contract's line anchors hold)"),
    },
    "UMAT_NKH_1.02": {
        "props": {1: 0.0},
        "why": ("with the probe's PROPS(1)=1.0 the source takes THTA=PROPS(1) and never sets DTHTA, "
                "then reads it for the thermal strain (line 111; a -finit-real=snan build traps "
                "there), so each build computes with whatever memory holds; the variant sets "
                "PROPS(1)=0, which takes the temperature and its increment from Abaqus"),
    },
    "UMAT_VPDCL_R": {
        "ntens": 4,
        "why": ("the source reads its back stress with DO K1=10,2*NTENS+5 into XBACK(NTENS), in "
                "bounds only for NTENS=4 (both archived 3-D Abaqus runs of it failed); the variant "
                "runs the plane-strain deck"),
    },
}
JOB_PREFIX = os.environ.get("PRESENTATION_ABAQUS_PREFIX4", "claudeP_c4")


def _safe(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]", "_", name)


def variant_config(config_path: Path, work: Path) -> Path:
    """A copy of the benchmark contract with the documented correction applied."""
    spec = VARIANTS[config_path.stem]
    raw = json.loads(config_path.read_text())
    source = (config_path.parent / raw["source"]).resolve()
    directory = fresh_dir(work / f"{config_path.stem}__variant_input")
    if spec.get("closure"):
        from umat_oti.transform.dependency_resolution import combined_source, resolve_closure
        graph = resolve_closure(source, entry="UMAT", roots=[source.parent])
        if graph.missing:
            raise RuntimeError(f"closure unresolved: {[m.symbol for m in graph.missing]}")
        resolved = directory / source.name
        resolved.write_text(combined_source(graph), encoding="utf-8")
        source = resolved
    raw["source"] = str(source)
    for name in spec.get("demote", []):
        raw["promote"] = [v for v in raw.get("promote", []) if v != name]
        raw["constant"] = sorted(set(raw.get("constant", [])) | {name})
    if "ntens" in spec:
        raw["ntens"] = spec["ntens"]
    target = directory / f"{config_path.stem}__variant.json"
    target.write_text(json.dumps(raw, indent=1))
    return target


def probe_props(source_path: Path, ntens: int, overrides: Dict[int, float]) -> List[float]:
    """The paired-validation probe vector (unit constants, 0.3 for Poisson-like
    names), exactly as the job builder writes it, with declared overrides."""
    from umat_oti.validation import job_builder as jb
    text = source_path.read_text(errors="replace")
    _nstatv, nprops = jb.infer_validation_dimensions_from_source(text, ntens=ntens)
    values = [float(v) for v in jb._props_lines(nprops, text).replace("\n", "").split(",") if v.strip()]
    for index, value in overrides.items():
        values[int(index) - 1] = float(value)
    return values


def prepare(config_path: Path, work: Path, props_override: Dict[int, float] | None = None) -> Dict[str, object]:
    """Transform one benchmark contract and build its paired validation workspace."""
    from umat_oti.cli_json import run_config_transform
    from umat_oti.core.config_loader import load_project_config_json
    from umat_oti.core.transformation_anchors import merge_completed_anchors_into_config
    from umat_oti.validation.job_builder import build_validation_workspace

    sys.path.insert(0, str(umat_repo() / "tools"))
    import run_completed_json_batch as batch  # the tool's own mode/tolerance helpers

    name = config_path.stem
    config = load_project_config_json(config_path.read_bytes(), origin_path=config_path)
    source_path = Path(str(config.get("source", {}).get("selected_umat_file", "")))
    if source_path.is_file():
        config = merge_completed_anchors_into_config(
            config, source_path.read_text(encoding="utf-8", errors="replace"))
    ntens = batch._as_int((config.get("transformation_settings", {}) or {}).get("ntens"))
    try:
        shown = str(config_path.relative_to(umat_repo()))
    except ValueError:
        shown = str(config_path)
    record: Dict[str, object] = {"case": name, "config": shown, "source": str(source_path)}
    transform_dir = fresh_dir(work / name / "oti_transform")
    summary, exit_code = run_config_transform(config_path, transform_dir)
    record["transform_exit"] = exit_code
    record["transform_success"] = bool(summary.get("transform_success"))
    record["blockers"] = [str(b) for b in summary.get("blockers", [])]
    if not record["transform_success"] or record["blockers"]:
        record["status"] = "transform_failed"
        return record
    validation_dir = fresh_dir(work / name / "validation")
    material_props = probe_props(source_path, ntens, props_override) if props_override else None
    record["material_props_override"] = props_override
    requested = batch._validation_compare_outputs(config)
    record["compare_outputs_in_contract"] = requested
    record["compare_outputs_used"] = FORCED_OUTPUTS
    build_validation_workspace(
        validation_dir=validation_dir, original_umat=source_path,
        transformed_umat=Path(str(summary.get("transformed_source"))), generated_dir=transform_dir,
        ntens=ntens, abaqus_command="abaqus", abaqus_modules="", run_prefix="",
        material_test_mode=batch._material_test_mode(config), run_compile_smoke=False,
        compare_outputs=FORCED_OUTPUTS, material_props=material_props,
        comparison_abs_tolerance=batch._validation_float(config, "absolute_tolerance"),
        comparison_rel_tolerance=batch._validation_float(config, "relative_tolerance"),
        comparison_ddsdde_abs_tolerance=batch._validation_float(config, "ddsdde_absolute_tolerance"),
        comparison_ddsdde_rel_tolerance=batch._validation_float(config, "ddsdde_relative_tolerance"),
    )
    record["validation_dir"] = str(validation_dir)
    record["material_test_mode"] = batch._material_test_mode(config)
    record["status"] = "prepared"
    return record


def _user_file(script: Path) -> str:
    match = re.search(r"user='([^']+)'", script.read_text())
    if not match:
        raise RuntimeError(f"no user= file in {script}")
    return match.group(1)


def run_pair(validation_dir: Path, case: str) -> Dict[str, object]:
    """Both jobs, one at a time, under claudeP names; outputs renamed for the extractor."""
    from umat_oti.validation.abaqus_runner import (
        AbaqusRunResult, completed_despite_teardown_abort,
    )
    from umat_oti.validation.job_builder import update_validation_report

    runs = {}
    for side, expected, script, key in (
            ("orig", "original_umat_validation", "run_original_abaqus.sh", "original_run_status"),
            ("oti", "otis_umat_validation", "run_otis_abaqus.sh", "transformed_run_status")):
        job = f"{JOB_PREFIX}_{_safe(case)}_{side}"
        user = _user_file(validation_dir / script)
        command = ["abaqus", f"job={job}", f"input={expected}.inp", f"user={user}", "double=both",
                   "interactive"]
        started = time.perf_counter()
        process = subprocess.run(command, cwd=validation_dir, capture_output=True, text=True, timeout=3600)
        seconds = time.perf_counter() - started
        (validation_dir / f"{side}_abaqus_stdout.log").write_text(process.stdout)
        (validation_dir / f"{side}_abaqus_stderr.log").write_text(process.stderr)
        for produced in validation_dir.glob(f"{job}.*"):
            produced.rename(validation_dir / (expected + produced.name[len(job):]))
        if process.returncode == 0:
            status, message = "completed", ""
        else:
            teardown = completed_despite_teardown_abort(validation_dir, expected)
            status = "completed_after_teardown_abort" if teardown else "failed"
            message = teardown
        result = AbaqusRunResult(status, command, validation_dir, process.returncode,
                                 validation_dir / f"{side}_abaqus_stdout.log",
                                 validation_dir / f"{side}_abaqus_stderr.log", message,
                                 process.stdout[-2000:], process.stderr[-2000:])
        update_validation_report(validation_dir, {key: result.to_json()})
        sta = validation_dir / f"{expected}.sta"
        runs[side] = {"job": job, "status": status, "returncode": process.returncode,
                      "seconds": round(seconds, 1),
                      "sta_completed": sta.is_file()
                      and "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" in sta.read_text(errors="replace")}
    return runs


def compare(validation_dir: Path) -> Dict[str, object]:
    from umat_oti.validation.abaqus_runner import extract_results
    from umat_oti.validation.compare_results import compare_validation_results

    extraction = extract_results(validation_dir, "abaqus", "", "")
    compare_validation_results(validation_dir)
    report = json.loads((validation_dir / "comparison_report.json").read_text())
    dd = report.get("ddsdde_comparison") or {}
    st = report.get("stress_comparison") or {}
    sv = report.get("state_variable_comparison") or {}
    return {"extraction": extraction.status, "overall_pass": report.get("pass"),
            "overall_status": report.get("status"),
            "ddsdde_status": dd.get("status"), "ddsdde_pass": dd.get("pass"),
            "ddsdde_max_abs": dd.get("max_abs_difference"), "ddsdde_max_rel": dd.get("max_rel_difference"),
            "ddsdde_compared_increments": dd.get("compared_increment_count"),
            "stress_pass": st.get("pass"), "stress_max_abs": st.get("max_abs_difference"),
            "stress_max_rel": st.get("max_rel_difference"),
            "statev_status": sv.get("status"), "statev_pass": sv.get("pass"),
            "errors": report.get("errors", [])}


def classify(row: Dict[str, object]) -> str:
    if not row.get("compared"):
        return "not compared"
    rel = row["compared"].get("ddsdde_max_rel")
    if rel is None:
        return "DDSDDE unavailable"
    return "exact" if float(rel) == 0.0 else "differs"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--case", action="append", dest="cases")
    parser.add_argument("--abaqus", action="store_true", help="run the Abaqus jobs (one at a time)")
    parser.add_argument("--variants", action="store_true",
                        help="also run the documented input corrections (separate rows)")
    parser.add_argument("--out", type=Path, default=default_out())
    parser.add_argument("--work", type=Path, default=default_work() / "claim4")
    args = parser.parse_args(argv)
    origins = require_worktree_imports()
    configs = sorted((umat_repo() / "benchmarks").glob("*.json"))
    if args.cases:
        configs = [c for c in configs if c.stem in args.cases]
    work = args.work.resolve(); work.mkdir(parents=True, exist_ok=True)
    rows: List[Dict[str, object]] = []
    jobs = [(config, None) for config in configs]
    if args.variants:
        jobs += [(config, "variant") for config in configs if config.stem in VARIANTS]
    for config, variant in jobs:
        label = config.stem + ("__variant" if variant else "")
        print(f"[claim4] {label}", flush=True)
        try:
            source_config = variant_config(config, work) if variant else config
            row = prepare(source_config, work,
                          VARIANTS[config.stem].get("props") if variant else None)
        except Exception as exc:  # noqa: BLE001 - recorded, never skipped
            row = {"case": config.stem, "status": "error", "error": f"{type(exc).__name__}: {exc}"[:500],
                   "transform_success": False, "blockers": []}
        row["case"] = config.stem
        row["variant"] = VARIANTS[config.stem]["why"] if variant else None
        if args.abaqus and row["status"] == "prepared":
            if not shutil.which("abaqus"):
                row["status"] = "abaqus_not_on_path"
            else:
                row["runs"] = run_pair(Path(row["validation_dir"]), label)
                row["compared"] = compare(Path(row["validation_dir"]))
                row["status"] = "compared"
        row["classification"] = classify(row)
        row["slide"] = SLIDE.get(config.stem)
        c = row.get("compared") or {}
        print(f"          {row['status']}: DDSDDE {format_e(c.get('ddsdde_max_abs'), 3)} / "
              f"{format_e(c.get('ddsdde_max_rel'), 2)}  pass={c.get('overall_pass')}  slide {row['slide']}",
              flush=True)
        rows.append(row)
    variant_rows = [r for r in rows if r.get("variant")]
    rows_main = [r for r in rows if not r.get("variant")]
    compared = [r for r in rows_main if r.get("compared")]
    on_slide = [r for r in compared if r["case"] in SLIDE]
    summary = {
        "cases": len(rows_main), "compared": len(compared),
        "overall_pass": sum(1 for r in compared if r["compared"]["overall_pass"]),
        "ddsdde_exact": sum(1 for r in compared if r["classification"] == "exact"),
        "ddsdde_differs": sum(1 for r in compared if r["classification"] == "differs"),
        "slide_cases_compared": len(on_slide),
        "slide_cases_pass": sum(1 for r in on_slide if r["compared"]["overall_pass"]),
        "slide_cases_exact": sum(1 for r in on_slide if r["classification"] == "exact"),
        "not_compared": [r["case"] for r in rows_main if not r.get("compared")],
        "variants": {r["case"]: {"status": r["status"], "classification": r["classification"],
                                 "ddsdde_max_abs": (r.get("compared") or {}).get("ddsdde_max_abs"),
                                 "ddsdde_max_rel": (r.get("compared") or {}).get("ddsdde_max_rel"),
                                 "overall_pass": (r.get("compared") or {}).get("overall_pass")}
                     for r in variant_rows},
        "slide": {"verified": "18/18", "exact": 12, "differ": 6},
    }
    payload = {"claim": "slide 8", "summary": summary, "rows": rows,
               "procedure": "run_config_transform + build_validation_workspace (DDSDDE forced) + "
                            "paired Abaqus jobs + extract_results + compare_validation_results",
               "provenance": provenance(), "imports": origins}
    out = args.out.resolve()
    write_json(out / "claim4_benchmark_ddsdde.json", payload)
    print(json.dumps(summary, indent=1))
    print(f"wrote {out / 'claim4_benchmark_ddsdde.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
