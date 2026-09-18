#!/usr/bin/env python3
"""Run every presentation-claim reproduction and write one summary.

    python presentation/run_all.py                 # offline parts of all claims
    python presentation/run_all.py --abaqus        # plus the Abaqus jobs (claims 3, 4)
    python presentation/run_all.py --quick         # reduced cases, a few minutes

Environment (see docs/PRESENTATION_CLAIMS.md): PYTHONPATH must put this
repository, the UMAT repository's ``src`` and the OTILib build first;
``UMAT_OTI_REPO`` names the UMAT repository. Large intermediate files go to
``--work`` (default $PRESENTATION_WORK or /tmp), the small JSON/CSV results to
``--out`` (default presentation/results/, git-ignored).

The summary ``presentation_summary.json`` lists, per claim, the slide's value
and the measured one. ``presentation/expected/presentation_summary.json`` is the
committed copy from the run documented in docs/PRESENTATION_CLAIMS.md.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from presentation import (  # noqa: E402
    claim1_sensitivity_sweep, claim2_flowrule_jacobian, claim3_cp_residual_c3d8,
    claim4_benchmark_ddsdde, claim5_constitutive_jacobians,
)
from presentation.common import default_out, default_work, provenance, require_worktree_imports, write_json  # noqa: E402


def _load(out: Path, name: str):
    path = out / name
    return json.loads(path.read_text()) if path.is_file() else None


def summarise(out: Path) -> dict:
    summary: dict = {}
    c1 = _load(out, "claim1_sensitivity_sweep.json")
    if c1:
        s = c1["summary"]
        summary["claim1_slide13"] = {
            "slide": {"models": 18, "directions": 76, "worst": 1.6e-7, "limit": 1e-5},
            "measured": {"models_attempted": s["models_attempted"],
                         "models_both_below_1e-5": s["models_both_below_1e-5"],
                         "directions": s["parameter_directions_in_passing_models"],
                         "worst": s["worst_among_passing"],
                         "strict_max_over_path_models_below_1e-5": s["models_both_below_1e-5_strict_max_over_path"],
                         "umat_sweep_tool_funnel": (c1.get("umat_sweep_tool", {}).get("_run") or {}).get("funnel")},
            "per_model": {r["model"]: {"p1": r["p1_dsigma_dp"], "p2": r["p2_dsigmavm_dp"],
                                       "slide_p1": r["slide_p1"], "slide_p2": r["slide_p2"]}
                          for r in s["rows"]}}
    c2 = _load(out, "claim2_flowrule_jacobian.json")
    if c2:
        per_compiler = {}
        for compiler, result in c2["results"].items():
            a, t = result["accuracy"]["summary"], result["timing"]["ratios"]
            per_compiler[compiler] = {
                "oti_vs_analytical_max_rel": a["oti_vs_analytical_max_rel"],
                "fd_centred_best_over_oti": a["ratio_fd_centred_best_to_oti_smooth"],
                "fd_centred_plateau_over_oti": a.get("ratio_fd_centred_plateau_to_oti_smooth"),
                "fd_forward_best_over_oti": a["ratio_fd_forward_best_to_oti_all"],
                "time_oti_over_analytical": t["oti_over_analytical"],
                "time_fd_centred_known_step_over_oti": t["fd_centred_known_step_over_oti"],
                "time_fd_forward_known_step_over_oti": t["fd_forward_known_step_over_oti"],
                "time_fd_centred_with_step_search_over_oti": t["fd_centred_step_search_over_oti"]}
        summary["claim2_slides26_27"] = {"slide": c2["slide"], "measured": per_compiler,
                                         "sources": c2["sources"]}
    c3 = _load(out, "claim3_cp_residual_c3d8.json")
    if c3:
        summary["claim3_slides28_32"] = {
            "slide": {"nrmse_limit": 1e-8, "mesh": "same values on a larger mesh"},
            "measured": {
                "sigma_vm_final_MPa": c3["single_element"]["sigma_vm_final_MPa"],
                "nrmse_oti_vs_analytic_worst": c3["analytic_chain_rule"]["nrmse_worst"],
                "fd_best_case_nrmse_worst": c3["central_fd"]["best_case_nrmse_worst"],
                "mesh_vs_single_max_rel": c3["larger_mesh"]["dvm_mesh_vs_single_max_rel"],
                "regime_checks_refined_path": c3["weighted_sensitivities_refined_path"]["checks"],
                "cost_kernel_x_plain": c3["cost"]["compiled_kernel"]["normalised_to_plain"],
                "cost_python_x_plain": c3["cost"]["python_analysis"]["normalised_to_plain"],
                "abaqus": {k: c3["abaqus"].get(k) for k in ("status", "nominal", "larger_mesh",
                                                            "fd_of_abaqus_best_step_vs_oti_worst")}}}
    c4 = _load(out, "claim4_benchmark_ddsdde.json")
    if c4:
        summary["claim4_slide8"] = {
            "slide": c4["summary"]["slide"],
            "measured": {k: v for k, v in c4["summary"].items() if k != "slide"},
            "rows": {r["case"] + ("__variant" if r.get("variant") else ""):
                     {"status": r["status"], "classification": r["classification"],
                      "ddsdde_max_abs": (r.get("compared") or {}).get("ddsdde_max_abs"),
                      "ddsdde_max_rel": (r.get("compared") or {}).get("ddsdde_max_rel"),
                      "pass": (r.get("compared") or {}).get("overall_pass"), "slide": r.get("slide")}
                     for r in c4["rows"]}}
    c5 = _load(out, "claim5_constitutive_jacobians.json")
    if c5:
        summary["claim5_slide25"] = {"slide": c5["summary"]["slide"],
                                     "measured": {k: v for k, v in c5["summary"].items()
                                                  if k not in ("slide", "slide_table")}}
    return summary


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--abaqus", action="store_true")
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--only", action="append", choices=["1", "2", "3", "4", "5"])
    parser.add_argument("--summary-only", action="store_true",
                        help="only rebuild presentation_summary.json from the results already in --out")
    parser.add_argument("--out", type=Path, default=default_out())
    parser.add_argument("--work", type=Path, default=default_work())
    args = parser.parse_args(argv)
    require_worktree_imports()
    out, work = args.out.resolve(), args.work.resolve()
    common = ["--out", str(out)]
    plan = {
        "1": (claim1_sensitivity_sweep.main,
              common + ["--work", str(work / "claim1")]
              + (["--model", "m1_elastic", "--model", "m3_j2", "--model", "m5_cpflow", "--increments", "60"]
                 if args.quick else [])),
        "2": (claim2_flowrule_jacobian.main,
              common + ["--work", str(work / "claim2")] + (["--calls", "5000", "--repetitions", "2"] if args.quick else [])),
        "3": (claim3_cp_residual_c3d8.main,
              common + ["--work", str(work / "claim3")]
              + (["--refine", "2", "--mesh", "2", "--timing-reps", "2", "--kernel-reps", "20"] if args.quick else [])
              + (["--abaqus", "--abaqus-work", str(work / "claim3_abaqus")] if args.abaqus else [])),
        "4": (claim4_benchmark_ddsdde.main,
              common + ["--work", str(work / "claim4"), "--variants"] + (["--abaqus"] if args.abaqus else [])),
        "5": (claim5_constitutive_jacobians.main, common + ["--work", str(work / "claim5")]),
    }
    status = {}
    for key in ([] if args.summary_only else (args.only or ["1", "2", "3", "4", "5"])):
        function, arguments = plan[key]
        started = time.perf_counter()
        print(f"\n===== claim {key}: {' '.join(arguments)}", flush=True)
        try:
            code = function(arguments)
        except SystemExit as exc:
            code = exc.code if isinstance(exc.code, int) else 1
        status[key] = {"exit": code, "seconds": round(time.perf_counter() - started, 1)}
    payload = {"status": status, "abaqus": args.abaqus, "quick": args.quick,
               "claims": summarise(out), "provenance": provenance()}
    write_json(out / "presentation_summary.json", payload)
    print(json.dumps(status, indent=1))
    print(f"wrote {out / 'presentation_summary.json'}")
    return 0 if all(v["exit"] == 0 for v in status.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
