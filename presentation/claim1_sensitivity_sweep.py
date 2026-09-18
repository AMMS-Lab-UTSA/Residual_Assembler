#!/usr/bin/env python3
"""Slide 13: parameter sensitivities of the 20-model collection, both programs.

What the slide says
    "18 material models, 76 parameter directions verified against finite
    differences, worst-case agreement 1.6e-7, both levels stay below 1e-5."
    Columns: Dsigma_Dp (HYPAD vs fd) and DsigmaVm_Dp (HYPAD vs fd).

What produced the slide (recovered from RA origin/cross-platform-hardening,
results/build_sweep_tables.py + results/figures/sweep_error_tables.json)
    20 materials were run; the 18 on the slide are the ones below 1e-5 on both
    programs (sweep_drucker_prager P1 6.7e-3 and sweep_perzyna_linear P2 1.5e-4
    were dropped; their 8 directions make 84 - 8 = 76). Program 1 was the OTI
    DSIGMA_DP of the compiled provider marched over a 150-increment uniaxial
    strain ramp to 3 %, compared with centred FD of the ORIGINAL UMAT, reported
    as RMSE(OTI - FD) / max|FD| per parameter. Program 2 was d(sigma_vM)/dp at a
    single material point (the docstring said "single C3D8" but the code marched
    one integration point: under prescribed homogeneous strain du/dp = 0), at
    the LAST increment of a 120-increment ramp, |OTI - FD| / |FD|.

What this script measures (current code, nothing copied)
    For every model in the UMAT repository's own model list:

    Program 1a  the UMAT sweep tool ``tools/run_parameter_sensitivity_sweep.py``
                (source transform + OTI material-point driver vs centred FD of the
                separately compiled ORIGINAL), all 20 models, its own funnel.
    Program 1b  ``umat-oti-provider build`` -> compiled object -> its DSIGMA_DP
                through UMAT_OTI_EVAL over the complete loading path (incoming
                derivatives carried for path-dependent models; stateless models
                evaluated at the total strain from the virgin state, as the old
                eval_point did), vs centred FD of the ORIGINAL UMAT bundled in the
                same object, re-marched over the whole path for every perturbation.
    Program 2   the Residual Assembler replay: ``residual_core.replay.path_material.
                PathMaterial`` is what links the compiled object and replays it;
                d(sigma_vM)/dp = (d sigma_vM / d sigma) : d sigma/dp at every
                increment, vs centred FD of sigma_vM of the ORIGINAL re-marched.
    Reference   step ladder 1e-3..1e-7 (relative), plateau-selected per parameter
                without looking at OTI; at an increment where a perturbed run
                switches between the elastic and inelastic branch (a kink at an
                increment boundary) the second-order one-sided stencil on the
                nominal side is used, and those increments are listed.
    Cross-checks the Program-1 transform driver's DSIGMA_DP (UMAT sweep tool)
                against the compiled provider's UMAT_OTI_MARCH on the identical
                path (independent builds and drivers), and MARCH vs EVAL.

Both the slide metric (RMSE / final increment) and a stricter max-over-path
metric are written. Output: ``claim1_sensitivity_sweep.json`` and ``.csv``.

    python presentation/claim1_sensitivity_sweep.py            # all 20 models
    python presentation/claim1_sensitivity_sweep.py --model m3_j2 --model m1_elastic
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from presentation.common import (  # noqa: E402
    FD_LADDER, build_provider, default_out, default_work, format_e,
    fresh_dir, plateau, provenance, require_worktree_imports, umat_repo, write_json,
)

TOLERANCE = 1.0e-5          # the slide's "both levels stay below 1e-5"
FLOOR = 1.0e-3              # inactive-parameter floor used by the slide's metric

#: The slide's table (slide 13), kept only to print it beside the measurement.
SLIDE = {
    "m1_elastic": (1.0e-10, 1.6e-10), "m2_cubic": (1.1e-10, 6.3e-11),
    "m3_j2": (7.7e-11, 4.5e-08), "m5_cpflow": (1.4e-10, 2.8e-10),
    "m6_fcc": (1.6e-07, 2.2e-08), "sweep_aniso_ortho": (1.4e-10, 1.4e-10),
    "sweep_damage_elastic": (2.8e-10, 4.9e-10), "sweep_eco": (8.6e-11, 1.6e-10),
    "sweep_j2_bilinear": (1.1e-10, 8.5e-09), "sweep_j2_combined": (1.1e-10, 1.1e-08),
    "sweep_j2_kinematic": (1.1e-10, 6.9e-08), "sweep_lame_elastic": (3.2e-10, 2.8e-10),
    "sweep_maxwell_ve": (3.0e-10, 2.2e-10), "sweep_mooney_small": (3.7e-10, 6.8e-08),
    "sweep_real_ECL_TEMP": (1.1e-10, 5.5e-08), "sweep_real_PCO": (8.0e-11, 4.5e-08),
    "sweep_thermoelastic": (8.6e-11, 1.6e-10), "sweep_transiso": (5.4e-09, 1.5e-08),
}
SLIDE_SUMMARY = {"models": 18, "directions": 76, "worst": 1.6e-7}

PHYSICS = {
    "m1_elastic": "Isotropic elasticity", "m2_cubic": "Cubic anisotropic elasticity",
    "m3_j2": "J2 plasticity, linear hardening", "m5_cpflow": "Viscoplastic flow (thermally activated)",
    "m6_fcc": "FCC crystal, 12 slip systems", "sweep_aniso_ortho": "Orthotropic elasticity",
    "sweep_damage_elastic": "Elasticity + scalar damage", "sweep_drucker_prager": "Drucker-Prager plasticity",
    "sweep_eco": "Cosserat elasticity (ECO)", "sweep_j2_bilinear": "J2 plasticity, bilinear",
    "sweep_j2_combined": "J2 plasticity, combined", "sweep_j2_kinematic": "J2 plasticity, kinematic",
    "sweep_lame_elastic": "Isotropic elasticity (Lame)", "sweep_maxwell_ve": "Maxwell viscoelasticity",
    "sweep_mooney_small": "Neo-Hookean (small strain)", "sweep_perzyna_linear": "Perzyna viscoplasticity",
    "sweep_real_ECL_TEMP": "Thermo-elasticity (ECL_TEMP)", "sweep_real_PCO": "Couple-stress plast. (PCO)",
    "sweep_thermoelastic": "Thermo-elasticity", "sweep_transiso": "Transversely isotropic",
}


def loading_path(n: int = 150, emax: float = 0.03):
    """The slide's Program-1 path: monotonic uniaxial strain to 3 %, dt = 1/n."""
    increments = [[emax / n, 0.0, 0.0, 0.0, 0.0, 0.0] for _ in range(n)]
    return increments, [1.0 / n] * n


def model_list() -> List[str]:
    tool = umat_repo() / "tools" / "run_parameter_sensitivity_sweep.py"
    out = subprocess.run([sys.executable, str(tool), "--list"], capture_output=True,
                         text=True, check=True)
    return [line.strip() for line in out.stdout.splitlines() if line.strip()]


# --------------------------------------------------------------------------- #
# Program 1a: the UMAT repository's own sweep tool
# --------------------------------------------------------------------------- #
def run_umat_sweep(models: List[str], work: Path) -> Dict[str, dict]:
    tool = umat_repo() / "tools" / "run_parameter_sensitivity_sweep.py"
    sweep_work = fresh_dir(work / "umat_sweep_work")
    results = sweep_work / "results"
    command = [sys.executable, str(tool), "--work-dir", str(sweep_work),
               "--results-dir", str(results)]
    for model in models:
        command += ["--model", model]
    started = time.perf_counter()
    completed = subprocess.run(command, capture_output=True, text=True, cwd=umat_repo())
    seconds = time.perf_counter() - started
    if completed.returncode != 0:
        return {"_run": {"returncode": completed.returncode, "seconds": seconds,
                         "error": (completed.stdout + completed.stderr)[-2000:]}}
    round_payload = json.loads((results / "parameter_sensitivity_round.json").read_text())
    rows: Dict[str, List[dict]] = {}
    with (results / "table6_comparison_rows.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            rows.setdefault(row["model"], []).append(row)
    out: Dict[str, dict] = {"_run": {"returncode": 0, "seconds": round(seconds, 1),
                                     "command": "python tools/run_parameter_sensitivity_sweep.py "
                                                + " ".join(f"--model {m}" for m in models)
                                                + " --work-dir <work> --results-dir <work>/results",
                                     "funnel": round_payload["funnel"],
                                     "failure_taxonomy": round_payload["failure_taxonomy"]}}
    for record in round_payload["models"]:
        model = record["model"]
        verified = record["stages"].get("derivatives_verified", {})
        stress_rows = [r for r in rows.get(model, []) if r["array"] == "DSIGMA_DP"]
        by_param: Dict[str, list] = {}
        for r in stress_rows:
            by_param.setdefault(r["parameter"], []).append(
                (int(r["increment"]), int(r["component"]), float(r["oti"]), float(r["reference"])))
        metric = _slide_metric_from_rows(by_param)
        out[model] = {
            "furthest_stage": record["furthest_stage"],
            "verdict": verified.get("status"),
            "parameters": record.get("parameters", []),
            "directions_verified": record.get("verified_parameter_directions", []),
            "rows": verified.get("rows"), "rows_agreeing": verified.get("rows_agreeing"),
            "rows_unresolved": verified.get("rows_reference_unresolved"),
            "tool_worst_relative_error": verified.get("worst_relative_error"),
            "slide_metric_rmse_rel_per_parameter": metric,
            "slide_metric_worst": max(metric.values()) if metric else None,
            "driver_dsigma": {p: [(i, c, o) for i, c, o, _ in v] for p, v in by_param.items()},
            "path": _sweep_path(model),
        }
    return out


def _slide_metric_from_rows(by_param: Dict[str, list]) -> Dict[str, float]:
    if not by_param:
        return {}
    ref_max = {p: max(abs(v[3]) for v in rows) for p, rows in by_param.items()}
    gscale = max(ref_max.values()) or 1.0
    out = {}
    for p, rows in by_param.items():
        err = np.array([v[2] - v[3] for v in rows])
        out[p] = float(np.sqrt(np.mean(err ** 2)) / max(ref_max[p], FLOOR * gscale))
    return out


def _sweep_path(model: str):
    contract = json.loads((umat_repo() / "parameter_sensitivity" / "contracts" / f"{model}.json")
                          .read_text())
    driver = contract["material_point_driver"]
    return {"dstran_per_increment": driver["dstran_per_increment"],
            "n_increments": int(driver["n_increments"]), "dtime": 1.0}


# --------------------------------------------------------------------------- #
# Program 1b + Program 2 through the compiled provider
# --------------------------------------------------------------------------- #
def run_model(model: str, work: Path, n_increments: int, sweep: dict | None) -> dict:
    from residual_core.replay.path_material import PathMaterial, dmises_dstress, mises

    v2_path = umat_repo() / "parameter_sensitivity" / "models" / model / "contract_v2.json"
    v2 = json.loads(v2_path.read_text())
    props = [float(v) for v in v2["validation"]["props_values"]]
    record: dict = {"model": model, "physics": PHYSICS.get(model, model),
                    "props": props, "loading_path": {
                        "kind": "uniaxial strain ramp", "increments": n_increments,
                        "total_strain_11": 0.03, "dt": 1.0 / n_increments}}
    build = build_provider(v2_path, work / model / "provider")
    record["provider_build"] = {k: v for k, v in build.items() if k != "contract"}
    if build["returncode"] != 0:
        record["status"] = "provider_build_failed"
        return record
    contract = build["contract"]
    material = PathMaterial(build["object"], contract, workdir=str(work / model / "link"))
    names = list(material.params)
    pidx = [i - 1 for i in material.props_index]
    path_dep = bool(contract["history"]["path_dependent"])
    record.update(parameters=names, props_index=material.props_index,
                  path_dependent=path_dep, eval_arguments=len(contract["symbols"]["oti_eval_signature"]),
                  has_march=material.has_march)
    increments, dts = loading_path(n_increments)
    nt, ns, npar = material.ntens, material.nstatev, material.nparam

    # ---- the provider's own whole-path export (UMAT_OTI_MARCH), cross-check only
    started = time.perf_counter()
    march = material.march_fast(props, increments, dts)
    record["seconds_march"] = time.perf_counter() - started
    march_dsig = np.array([m["dsigma_dp"] for m in march])          # (N, 6, P)

    # ---- DSIGMA_DP through the replay ABI (UMAT_OTI_EVAL, per increment) ------
    started = time.perf_counter()
    if path_dep:
        replay = material.march_oti(props, increments, dts)
        stress = np.array([r["stress"] for r in replay])
        dsig = np.array([r["dsigma_dp"] for r in replay])
        route = "EVAL per increment, DSIGMA_DP_IN/DSTATEV_DP_IN carried (PathMaterial.march_oti)"
    else:
        total = np.cumsum(np.asarray(increments), axis=0)
        stress, dsig = [], []
        for strain, dt in zip(total, dts):
            s, _st, _dd, ds, _dsv = material._oti_step(
                props, np.zeros(nt), np.zeros(ns), np.zeros((nt, npar)),
                np.zeros((ns, npar)), strain, dt)
            stress.append(s); dsig.append(ds)
        stress, dsig = np.array(stress), np.array(dsig)
        route = ("EVAL at the total strain from the virgin state (stateless 16-argument ABI "
                 "has no incoming-derivative carry)")
    record["seconds_replay"] = time.perf_counter() - started
    record["replay_route"] = route
    g = np.array([dmises_dstress(s) for s in stress])
    dvm = np.einsum("ni,nip->np", g, dsig)                           # (N, P)

    # ---- primal parity: replay primal vs ORIGINAL ---------------------------
    nominal = material.march_regular(props, increments, dts)
    sig_orig = np.array([r["stress"] for r in nominal])
    scale = max(float(np.max(np.abs(sig_orig))), 1e-30)
    record["primal_parity_scaled"] = float(np.max(np.abs(stress - sig_orig)) / scale)
    record["sigma_vm_final"] = mises(sig_orig[-1])
    nominal_branch = _branches(nominal, ns)
    record["inelastic_increments"] = int(sum(nominal_branch))

    # ---- independent reference: FD of the ORIGINAL UMAT, whole path ----------
    def original(prop_values):
        run = material.march_regular(prop_values, increments, dts)
        s = np.array([r["stress"] for r in run])
        vm = np.array([mises(row) for row in s])
        return np.concatenate([s, vm[:, None]], axis=1), _branches(run, ns)

    fd_s, fd_vm, fd_meta = {}, {}, {}
    started = time.perf_counter()
    for k, name in enumerate(names):
        ladder, flags = branch_aware_ladder(original, props, pidx[k], nominal_branch)
        pick_s = plateau({h: a[:, :6] for h, a in ladder.items()})
        pick_vm = plateau({h: a[:, 6] for h, a in ladder.items()})
        fd_s[name], fd_vm[name] = pick_s["value"], pick_vm["value"]
        fd_meta[name] = {"stress_step": pick_s["step"], "stress_uncertainty": pick_s["uncertainty_scaled"],
                         "vm_step": pick_vm["step"], "vm_uncertainty": pick_vm["uncertainty_scaled"],
                         "one_sided_increments": flags}
    record["seconds_fd"] = time.perf_counter() - started
    record["fd_reference"] = {
        "ladder": list(FD_LADDER), "selection": "plateau of consecutive steps",
        "stencil": ("centred (f(p+h)-f(p-h))/2h; at an increment where a perturbed run takes a "
                    "different elastic/inelastic branch than the nominal run (a kink at the "
                    "increment boundary, where the two one-sided derivatives differ), the "
                    "second-order one-sided stencil (-3f(p)+4f(p+-h)-f(p+-2h))/(+-2h) on the side "
                    "that stays on the nominal branch"),
        "per_parameter": fd_meta, "original_runs": 4 * len(FD_LADDER) * len(names) + 1}
    p1 = dsig
    p2_vm = dvm

    # ---- metrics -------------------------------------------------------------
    g1 = max(float(np.max(np.abs(fd_s[n]))) for n in names) or 1.0
    g2p = max(float(np.max(np.abs(fd_vm[n]))) for n in names) or 1.0
    g2f = max(abs(float(fd_vm[n][-1])) for n in names) or 1.0
    p1_rmse, p1_max, p2_final, p2_path, p2_rmse = {}, {}, {}, {}, {}
    for k, name in enumerate(names):
        d1 = max(float(np.max(np.abs(fd_s[name]))), FLOOR * g1)
        err1 = p1[:, :, k] - fd_s[name]
        p1_rmse[name] = float(np.sqrt(np.mean(err1 ** 2)) / d1)
        p1_max[name] = float(np.max(np.abs(err1)) / d1)
        err2 = p2_vm[:, k] - fd_vm[name]
        p2_final[name] = float(abs(err2[-1]) / max(abs(float(fd_vm[name][-1])), FLOOR * g2f))
        d2 = max(float(np.max(np.abs(fd_vm[name]))), FLOOR * g2p)
        p2_path[name] = float(np.max(np.abs(err2)) / d2)
        p2_rmse[name] = float(np.sqrt(np.mean(err2 ** 2)) / d2)
    record["program1_provider"] = {"metric": "RMSE_(path,components)(OTI-FD)/max|FD| (slide metric)",
                                   "source": "DSIGMA_DP of the compiled Program-1 provider through UMAT_OTI_EVAL",
                                   "per_parameter": p1_rmse, "worst": max(p1_rmse.values()),
                                   "max_abs_over_path_rel": p1_max, "worst_max_rel": max(p1_max.values())}
    record["program2_replay"] = {"metric": "|OTI-FD|/|FD| of d sigma_vM/dp at the final increment (slide metric)",
                                 "per_parameter_final": p2_final, "worst": max(p2_final.values()),
                                 "per_parameter_max_over_path": p2_path, "worst_max_over_path": max(p2_path.values()),
                                 "per_parameter_rmse_over_path": p2_rmse}
    # cross-check 1: the provider's MARCH export vs the EVAL replay (same object,
    # two entry points). MARCH runs at TEMP=293.15 like the Program-1 driver; the
    # replay shims pass TEMP=0, so a temperature-dependent model differs here by
    # construction (reported, not counted as a derivative error).
    record["crosscheck_march_vs_eval_replay"] = float(
        np.max(np.abs(march_dsig - p1)) / max(float(np.max(np.abs(p1))), 1e-30))
    record["reads_temperature"] = _reads_temperature(model)
    # cross-check 2: the Program-1 transform driver (UMAT sweep tool) vs the
    # compiled provider's MARCH on the identical path and temperature
    record["crosscheck_transform_driver_vs_replay"] = _driver_vs_replay(
        material, props, names, sweep, path_dep)
    record["status"] = "measured"
    return record


def _branches(run, nstatev: int) -> List[bool]:
    """Per increment: did the ORIGINAL update change its state (inelastic step)?"""
    previous = np.zeros(nstatev)
    out = []
    for row in run:
        state = np.asarray(row["statev"], float)
        out.append(bool(nstatev and np.any(state != previous)))
        previous = state
    return out


def branch_aware_ladder(original, props, index, nominal_branch):
    """Centred differences over the step ladder, one-sided where a stencil straddles a kink.

    Returns ({step: (N, 7) estimate}, {step: [flagged increments (1-based)]}).
    """
    base = float(props[index])
    scale = abs(base) if base != 0.0 else 1.0
    f0 = original(list(props))[0]
    ladder, flags = {}, {}
    nominal = np.asarray(nominal_branch)
    for step in FD_LADDER:
        h = step * scale
        runs = {}
        for mult in (-2, -1, 1, 2):
            prop = list(props); prop[index] = base + mult * h
            runs[mult] = original(prop)
        estimate = (runs[1][0] - runs[-1][0]) / (2.0 * h)
        plus_ok = (np.asarray(runs[1][1]) == nominal) & (np.asarray(runs[2][1]) == nominal)
        minus_ok = (np.asarray(runs[-1][1]) == nominal) & (np.asarray(runs[-2][1]) == nominal)
        straddle = ~((np.asarray(runs[1][1]) == nominal) & (np.asarray(runs[-1][1]) == nominal))
        flagged = []
        for n in np.flatnonzero(straddle):
            if plus_ok[n]:
                estimate[n] = (-3.0 * f0[n] + 4.0 * runs[1][0][n] - runs[2][0][n]) / (2.0 * h)
                flagged.append((int(n) + 1, "one-sided +"))
            elif minus_ok[n]:
                estimate[n] = (3.0 * f0[n] - 4.0 * runs[-1][0][n] + runs[-2][0][n]) / (2.0 * h)
                flagged.append((int(n) + 1, "one-sided -"))
            else:
                flagged.append((int(n) + 1, "unresolved (both sides leave the nominal branch)"))
        ladder[float(step)] = estimate
        if flagged:
            flags[f"{step:g}"] = flagged
    return ladder, flags


def _reads_temperature(model: str) -> bool:
    import re
    source = (umat_repo() / "parameter_sensitivity" / "models" / model / "umat.for").read_text(errors="replace")
    body = [line for line in source.splitlines() if line[:1] not in ("C", "c", "*", "!")]
    return any(re.search(r"=\s*.*\bTEMP\b", line) for line in body)


def _driver_vs_replay(material, props, names, sweep, path_dep):
    if not sweep or not sweep.get("driver_dsigma"):
        return None
    path = sweep["path"]
    increments = [list(path["dstran_per_increment"])] * path["n_increments"]
    dts = [path["dtime"]] * path["n_increments"]
    replay = material.march_fast(props, increments, dts)
    worst, scale = 0.0, 0.0
    for name, values in sweep["driver_dsigma"].items():
        k = names.index(name)
        for increment, component, driver_value in values:
            replay_value = float(replay[increment - 1]["dsigma_dp"][component - 1, k])
            worst = max(worst, abs(replay_value - driver_value))
            scale = max(scale, abs(driver_value))
    return {"max_abs_difference_scaled": worst / (scale or 1.0),
            "path": f"{path['n_increments']} increments of {path['dstran_per_increment']}",
            "compared": "UMAT sweep OTI driver DSIGMA_DP vs provider MARCH DSIGMA_DP"}


# --------------------------------------------------------------------------- #
def summarise(records: List[dict], sweep: dict) -> dict:
    measured = [r for r in records if r.get("status") == "measured"]
    rows = []
    for r in records:
        s = sweep.get(r["model"], {}) if sweep else {}
        p1 = r.get("program1_provider", {}).get("worst")
        p2 = r.get("program2_replay", {}).get("worst")
        both = p1 is not None and p2 is not None and p1 < TOLERANCE and p2 < TOLERANCE
        rows.append({
            "model": r["model"], "physics": r.get("physics"),
            "parameters": len(r.get("parameters", [])),
            "p1_dsigma_dp": p1, "p2_dsigmavm_dp": p2,
            "p1_max_over_path": r.get("program1_provider", {}).get("worst_max_rel"),
            "p2_max_over_path": r.get("program2_replay", {}).get("worst_max_over_path"),
            "p1_umat_sweep_tool": s.get("slide_metric_worst"),
            "umat_sweep_verdict": s.get("verdict"),
            "crosscheck_march_vs_replay": r.get("crosscheck_march_vs_eval_replay"),
            "crosscheck_driver_vs_replay": (r.get("crosscheck_transform_driver_vs_replay") or {}).get(
                "max_abs_difference_scaled"),
            "both_below_1e-5": both,
            "slide_p1": SLIDE.get(r["model"], (None, None))[0],
            "slide_p2": SLIDE.get(r["model"], (None, None))[1],
            "on_slide": r["model"] in SLIDE,
            "status": r.get("status"),
        })
    passing = [row for row in rows if row["both_below_1e-5"]]
    worst_pass = max((max(row["p1_dsigma_dp"], row["p2_dsigmavm_dp"]) for row in passing), default=None)
    worst_all = max((max(row["p1_dsigma_dp"], row["p2_dsigmavm_dp"]) for row in rows
                     if row["p1_dsigma_dp"] is not None), default=None)
    strict = [row for row in rows if row["p1_max_over_path"] is not None
              and row["p1_max_over_path"] < TOLERANCE and row["p2_max_over_path"] < TOLERANCE]
    return {
        "models_attempted": len(records),
        "models_measured": len(measured),
        "models_both_below_1e-5": len(passing),
        "parameter_directions_in_passing_models": sum(row["parameters"] for row in passing),
        "parameter_directions_total": sum(row["parameters"] for row in rows),
        "worst_among_passing": worst_pass,
        "worst_over_all_models": worst_all,
        "models_both_below_1e-5_strict_max_over_path": len(strict),
        "failing_models": [row["model"] for row in rows if not row["both_below_1e-5"]],
        "slide": SLIDE_SUMMARY,
        "rows": rows,
    }


def write_csv(path: Path, rows: List[dict]) -> None:
    columns = ["model", "physics", "parameters", "p1_dsigma_dp", "p2_dsigmavm_dp",
               "p1_max_over_path", "p2_max_over_path", "p1_umat_sweep_tool", "umat_sweep_verdict",
               "crosscheck_march_vs_replay", "crosscheck_driver_vs_replay", "both_below_1e-5",
               "slide_p1", "slide_p2"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--model", action="append", dest="models")
    parser.add_argument("--increments", type=int, default=150)
    parser.add_argument("--out", type=Path, default=default_out())
    parser.add_argument("--work", type=Path, default=default_work() / "claim1")
    parser.add_argument("--skip-umat-sweep", action="store_true",
                        help="do not run the UMAT repository's sweep tool (Program 1a)")
    args = parser.parse_args(argv)
    origins = require_worktree_imports()
    models = args.models or model_list()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)

    sweep = {} if args.skip_umat_sweep else run_umat_sweep(models, work)
    records = []
    for model in models:
        print(f"[claim1] {model}", flush=True)
        try:
            records.append(run_model(model, work, args.increments, sweep.get(model)))
        except Exception as exc:  # noqa: BLE001 - recorded as a failure, never skipped
            records.append({"model": model, "status": "error", "error": f"{type(exc).__name__}: {exc}"})
        r = records[-1]
        if r.get("status") == "measured":
            print(f"          P1 {format_e(r['program1_provider']['worst'])}  "
                  f"P2 {format_e(r['program2_replay']['worst'])}  "
                  f"(path max P1 {format_e(r['program1_provider']['worst_max_rel'])}, "
                  f"P2 {format_e(r['program2_replay']['worst_max_over_path'])})", flush=True)
        else:
            print(f"          {r.get('status')}: {r.get('error', '')[:200]}", flush=True)

    summary = summarise(records, sweep)
    for record in records:
        # the per-row driver values are only needed for the cross-check
        (sweep.get(record["model"]) or {}).pop("driver_dsigma", None)
    payload = {"claim": "slide 13", "summary": summary,
               "umat_sweep_tool": sweep, "models": records,
               "provenance": provenance(), "imports": origins}
    out = args.out.resolve()
    write_json(out / "claim1_sensitivity_sweep.json", payload)
    write_csv(out / "claim1_sensitivity_sweep.csv", summary["rows"])
    print(f"\nmodels attempted {summary['models_attempted']}, both programs < 1e-5: "
          f"{summary['models_both_below_1e-5']} ({summary['parameter_directions_in_passing_models']} "
          f"directions); worst among them {format_e(summary['worst_among_passing'])}; "
          f"failing: {summary['failing_models']}")
    print(f"slide: {SLIDE_SUMMARY}")
    print(f"wrote {out / 'claim1_sensitivity_sweep.json'}")
    return 0 if summary["models_measured"] == summary["models_attempted"] else 1


if __name__ == "__main__":
    sys.exit(main())
