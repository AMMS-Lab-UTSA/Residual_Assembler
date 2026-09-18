#!/usr/bin/env python3
"""Slides 28-32: residual-method sensitivities of the crystal-plasticity flow model.

What the slides say
    Single C3D8 (8 IPs), bottom face fixed, top face sheared to 5 % strain, the
    thermally activated flow model with six parameters (tau0, dG, p, q, gam0, H):
    all six sensitivities of sigma_vM in one OTI run (slide 29); weighted
    sensitivities per regime (elastic ~ 0; tau0 grows at yield; dG decreases and q
    grows in plastic flow; H small); OTI vs an analytical chain rule, NRMSE < 1e-8
    (slide 30); the same values on a larger mesh (slide 31); plain UMAT vs HYPAD vs
    best-case central FD cost (slide 32).

What produced them (recovered from RA origin/cross-platform-hardening)
    results/cp_residual_sensitivities.py (one integration point, E12: 0 -> 0.05 in
    50 increments of dt = 1, OTI vs FD), results/cp_mesh_residual.py (4x4x4, affine
    field prescribed on the whole boundary, mesh vs single 1.0e-14) and
    results/cp_analysis_timing.py (a different model, m6_fcc). No analytic chain
    rule was implemented there; the old reference was FD only.

What this script does (current packages)
    * ``umat-oti-provider build`` of parameter_sensitivity/models/m5_cpflow.
    * One C3D8 (unit cube), bottom face y=0 fixed, top face y=1: u1 = gamma*H,
      u2 = 0, u3 free; gamma 0 -> 0.1 (tensor E12 0 -> 0.05) in 50 increments of
      dt = 1 (the old study's path). Nonlinear solve with the ORIGINAL UMAT, then
      the residual method (presentation/c3d8_residual.py): K du/dp = -dR/dp with
      the provider's carried DSIGMA_DP, all six parameters in one enriched run.
    * References: (a) hand-derived chain rule of the same update
      (presentation/cpflow_analytic.py), NRMSE = RMSE / max|reference| per
      parameter; (b) centred FD of the ORIGINAL UMAT, the whole C3D8 analysis
      re-solved per perturbation, over a step ladder.
    * Weighted sensitivities |p dsigma_vM/dp| and their shares per regime, on the
      old path and on a 10x finer path at the same strain rate (regimes resolved).
    * 4x4x4 mesh under the study's homogeneous shear (affine field on every
      boundary node, interior free): element-average sensitivities vs one element.
    * Cost: (a) the Python analysis driver, plain vs enriched vs 13 re-solves;
      (b) the compiled kernel (presentation/cp_timing_driver.f90 linked with the
      same provider object): plain UMAT vs UMAT_OTI_EVAL vs UMAT_OTI_MARCH.
    * ``--abaqus``: the same single C3D8 in Abaqus with the unchanged umat.for
      (nominal + perturbed jobs, one at a time, double-precision .fil output).

    python presentation/claim3_cp_residual_c3d8.py [--abaqus]
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from presentation import c3d8_residual as C  # noqa: E402
from presentation import cpflow_analytic as A  # noqa: E402
from presentation.common import (  # noqa: E402
    PRESENTATION, build_provider, default_out, default_work, format_e, fresh_dir, provenance,
    require_worktree_imports, umat_repo, write_json,
)

MODEL = "m5_cpflow"
GAMMA = 0.1            # engineering shear, i.e. tensor E12 = 0.05 ("5 % strain")
LADDER = (1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7)
SLIDE = {"nrmse_limit": 1e-8, "mesh": "4x4x4 (old study)", "mesh_vs_single_old": 1.0473501463248156e-14}


def material():
    from residual_core.replay.path_material import PathMaterial
    v2 = umat_repo() / "parameter_sensitivity" / "models" / MODEL / "contract_v2.json"
    props = [float(v) for v in json.loads(v2.read_text())["validation"]["props_values"]]
    return v2, props, PathMaterial


def element_series(problem, solution):
    vm = np.array([C.volume_average(problem, C.mises_field(s)) for s in solution.stress])
    out = {"vm": vm}
    if solution.dstress_dp is not None:
        out["dvm"] = np.array([C.volume_average(problem, C.dmises_field(s, d))
                               for s, d in zip(solution.stress, solution.dstress_dp)])
    return out


def nrmse(values: np.ndarray, reference: np.ndarray) -> np.ndarray:
    return np.sqrt(np.mean((values - reference) ** 2, axis=0)) / np.max(np.abs(reference), axis=0)


def regimes(eqp: np.ndarray, dgamma: float, names, weights: np.ndarray, vm: np.ndarray):
    """Classify increments by the plastic share of the equivalent strain increment."""
    frac = np.diff(np.concatenate([[0.0], eqp])) / (dgamma / np.sqrt(3.0))
    steady = float(np.max(frac))
    labels = np.where(frac < 0.05, "elastic", np.where(frac >= 0.95 * steady, "plastic_flow", "yield_transition"))
    share = 100.0 * weights / np.maximum(weights.sum(axis=1, keepdims=True), 1e-300)
    report = {"plastic_fraction_definition": "d eqp / (d gamma_12 / sqrt 3); elastic < 0.05, "
                                             "plastic flow >= 0.95 x its steady value, transition between",
              "steady_plastic_fraction": steady, "increments": {}}
    for label in ("elastic", "yield_transition", "plastic_flow"):
        idx = np.flatnonzero(labels == label)
        if not idx.size:
            report["increments"][label] = {"count": 0}
            continue
        first, last = idx[0], idx[-1]
        report["increments"][label] = {
            "count": int(idx.size), "first": int(first) + 1, "last": int(last) + 1,
            "share_first_percent": dict(zip(names, share[first].round(3).tolist())),
            "share_last_percent": dict(zip(names, share[last].round(3).tolist())),
            "sum_weighted_over_vm_max": float(np.max(weights[idx].sum(axis=1) / np.maximum(vm[idx], 1e-300))),
        }
    k = {n: i for i, n in enumerate(names)}
    el = report["increments"]["elastic"]; tr = report["increments"]["yield_transition"]
    pl = report["increments"]["plastic_flow"]
    checks = {}
    if el.get("count"):
        checks["elastic: sum|p dvm/dp| / sigma_vM (max over elastic increments)"] = el["sum_weighted_over_vm_max"]
    if tr.get("count"):
        start = max(tr["first"] - 2, 0)
        checks["transition: tau0 share at start -> end (%)"] = [float(share[start, k["tau0"]]),
                                                                  float(share[tr["last"] - 1, k["tau0"]])]
        checks["transition: tau0 share rises"] = bool(share[tr["last"] - 1, k["tau0"]] > share[start, k["tau0"]])
    if pl.get("count"):
        a, b = pl["first"] - 1, pl["last"] - 1
        for name in ("dG", "q"):
            checks[f"plastic flow: {name} share first -> last (%)"] = [float(share[a, k[name]]), float(share[b, k[name]])]
        checks["plastic flow: dG share decreases"] = bool(share[b, k["dG"]] < share[a, k["dG"]] - 1e-9)
        checks["plastic flow: q share grows"] = bool(share[b, k["q"]] > share[a, k["q"]] + 1e-9)
    total = weights.sum(axis=1)
    checks["sum|p dvm/dp| at increment 1 / at the end"] = float(total[0] / total[-1])
    if el.get("count"):
        checks["sum|p dvm/dp| at the end of the elastic regime / at the end"] = float(total[el["last"] - 1] / total[-1])
    mean_plastic = share[labels == "plastic_flow"].mean(axis=0) if pl.get("count") else share.mean(axis=0)
    ranking = sorted(zip(names, mean_plastic.tolist()), key=lambda item: -item[1])
    checks["ranking by mean share in plastic flow (%)"] = [[n, round(v, 2)] for n, v in ranking]
    checks["three dominant parameters"] = [n for n, _ in ranking[:3]]
    checks["slide: tau0, dG, q dominate"] = sorted(n for n, _ in ranking[:3]) == sorted(["tau0", "dG", "q"])
    checks["H share: max over the path (%)"] = float(np.max(share[:, k["H"]]))
    checks["H share: at the end (%)"] = float(share[-1, k["H"]])
    checks["whole path: dG share first -> last (%)"] = [float(share[0, k["dG"]]), float(share[-1, k["dG"]])]
    checks["whole path: q share first -> last (%)"] = [float(share[0, k["q"]]), float(share[-1, k["q"]])]
    report["checks"] = checks
    report["labels"] = labels.tolist()
    report["share_percent"] = share
    return report


def strain_path_of(problem, solution):
    """The strain increments integration point 1 of element 1 actually saw."""
    previous = [np.zeros_like(solution.u[0])] + solution.u[:-1]
    return np.array([C.batch_strain(problem.ops, (u - up)[problem.edofs])[0, 0]
                     for u, up in zip(solution.u, previous)])


def run_single(mat, props, names, n_increments, dt):
    problem = C.simple_shear_problem(1, GAMMA, n_increments, dt)
    started = time.perf_counter()
    nominal = C.solve(mat, props, problem, mode="regular")
    t_plain = time.perf_counter() - started
    started = time.perf_counter()
    enriched = C.solve(mat, props, problem, mode="oti", sensitivities=True)
    t_enriched = time.perf_counter() - started
    series = element_series(problem, enriched)
    nominal_vm = element_series(problem, nominal)["vm"]
    # the analytic chain rule on the strain path the element actually saw
    ref = A.march(props, strain_path_of(problem, enriched), [dt] * n_increments, names)
    return problem, nominal, enriched, series, nominal_vm, ref, t_plain, t_enriched


def fd_reference(mat, props, pidx, problem, names):
    """Centred FD of sigma_vM of the ORIGINAL UMAT, the whole C3D8 analysis re-solved."""
    out = {}
    for k, name in enumerate(names):
        base = props[pidx[k]]
        per_step = {}
        for step in LADDER:
            h = step * abs(base)
            vals = []
            for sign in (+1, -1):
                prop = list(props); prop[pidx[k]] = base + sign * h
                sol = C.solve(mat, prop, problem, mode="regular")
                vals.append(element_series(problem, sol)["vm"])
            per_step[step] = (vals[0] - vals[1]) / (2.0 * h)
        out[name] = per_step
    return out


def kernel_timing(build: dict, props, work: Path, n_increments: int, dt: float, reps: int):
    exe_dir = fresh_dir(work / "kernel_timing")
    exe = exe_dir / "cp_timing"
    subprocess.run(["gfortran", "-O1", "-fPIC", str(PRESENTATION / "cp_timing_driver.f90"),
                    build["object"], "-o", str(exe)], check=True, capture_output=True, text=True)
    stdin = f"{reps} 8 {n_increments} 8 6 {dt}\n{' '.join(repr(p) for p in props)}\n0 0 0 {GAMMA / n_increments!r} 0 0\n"
    samples: Dict[str, List[float]] = {}
    for _ in range(5):
        out = subprocess.run([str(exe)], input=stdin, capture_output=True, text=True, check=True).stdout
        for line in out.splitlines():
            key, value = line.split()
            if key != "check":
                samples.setdefault(key, []).append(float(value))
    med = {k: statistics.median(v) for k, v in samples.items()}
    return {"compile": "gfortran -O1 -fPIC cp_timing_driver.f90 <provider object> (object built -O1 -fPIC "
                       "-fcheck=bounds by umat-oti-provider; ORIGINAL and OTI parts share those flags)",
            "integration_points": 8, "increments": n_increments, "repetitions_per_sample": reps,
            "samples": 5, "seconds_per_analysis_median": med, "samples_seconds": samples,
            "normalised_to_plain": {"plain_umat": 1.0,
                                    "hypad_eval_carry": med["eval"] / med["plain"],
                                    "hypad_march": med["march"] / med["plain"],
                                    "central_fd_13_runs": 13.0}}


def python_timing(mat, props, n_increments, dt, reps: int):
    problem = C.simple_shear_problem(1, GAMMA, n_increments, dt)
    C.solve(mat, props, problem, mode="regular")
    C.solve(mat, props, problem, mode="oti", sensitivities=True)
    plain, enriched = [], []
    for _ in range(reps):
        t = time.perf_counter(); C.solve(mat, props, problem, mode="regular"); plain.append(time.perf_counter() - t)
        t = time.perf_counter(); C.solve(mat, props, problem, mode="oti", sensitivities=True)
        enriched.append(time.perf_counter() - t)
    p, e = statistics.median(plain), statistics.median(enriched)
    return {"repetitions": reps, "warm_up": 1, "plain_seconds_median": p, "hypad_seconds_median": e,
            "central_fd_seconds_estimate": 13 * p,
            "normalised_to_plain": {"plain_umat": 1.0, "hypad": e / p, "central_fd_13_runs": 13.0},
            "note": "Python driver: element assembly and ctypes calls dominate both runs"}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--increments", type=int, default=50)
    parser.add_argument("--dt", type=float, default=1.0)
    parser.add_argument("--refine", type=int, default=10, help="refinement factor for the regime path")
    parser.add_argument("--mesh", type=int, default=4)
    parser.add_argument("--timing-reps", type=int, default=10)
    parser.add_argument("--kernel-reps", type=int, default=2000)
    parser.add_argument("--abaqus", action="store_true", help="also run the Abaqus jobs (one at a time)")
    parser.add_argument("--abaqus-work", type=Path, default=None)
    parser.add_argument("--out", type=Path, default=default_out())
    parser.add_argument("--work", type=Path, default=default_work() / "claim3")
    args = parser.parse_args(argv)
    origins = require_worktree_imports()
    work = args.work.resolve(); work.mkdir(parents=True, exist_ok=True)
    v2, props, PathMaterial = material()
    build = build_provider(v2, work / "provider")
    if build["returncode"]:
        raise SystemExit(build["error"])
    mat = PathMaterial(build["object"], build["contract"], workdir=str(work / "link"))
    names = list(mat.params)
    pidx = [i - 1 for i in mat.props_index]
    pvals = np.array([props[i] for i in pidx])
    record: Dict[str, object] = {"claim": "slides 28-32", "model": MODEL, "props": props,
                                 "parameters": names, "props_index": mat.props_index,
                                 "provider_build": {k: v for k, v in build.items() if k != "contract"}}

    # ---- single C3D8, the study's path ----------------------------------------
    n, dt = args.increments, args.dt
    problem, nominal, enriched, series, nominal_vm, ref, t_plain, t_enriched = run_single(mat, props, names, n, dt)
    dvm = series["dvm"]
    record["single_element"] = {
        "mesh": "1 C3D8, 8 integration points (Abaqus Gauss order, B-bar)",
        "boundary_conditions": "y=0 face: u=0; y=1 face: u1=gamma*H, u2=0; u3 of the top face free",
        "free_dofs": int(problem.free.size),
        "loading": f"gamma_12 0 -> {GAMMA} (E12 0 -> {GAMMA / 2}) in {n} increments, dt={dt}",
        "sigma_vm_final_MPa": float(series["vm"][-1]),
        "enriched_vs_plain_primal_max_rel": float(np.max(np.abs(series["vm"] - nominal_vm)) / np.max(nominal_vm)),
        "newton_iterations_max": int(max(enriched.iterations)),
        "max_du_dp": enriched.max_du_dp,
        "max_strain_sensitivity_increment": enriched.max_strain_sensitivity_increment,
        "sensitivities_in_one_run": len(names),
        "dvm_dp_final": dict(zip(names, dvm[-1].tolist())),
    }
    # ---- analytic chain rule -----------------------------------------------
    primal_gap = float(np.max(np.abs(ref["vm"] - series["vm"])) / np.max(ref["vm"]))
    err = nrmse(dvm, ref["dvm"])
    record["analytic_chain_rule"] = {
        "reference": "presentation/cpflow_analytic.py (hand-derived implicit-function chain rule)",
        "branches": sorted(set(ref["branch"])),
        "primal_sigma_vm_max_rel": primal_gap,
        "nrmse_definition": "sqrt(mean_n (OTI_n - ref_n)^2) / max_n |ref_n| per parameter",
        "nrmse_per_parameter": dict(zip(names, err.tolist())),
        "nrmse_worst": float(err.max()), "slide_limit": SLIDE["nrmse_limit"],
        "passes_slide_limit": bool(err.max() < SLIDE["nrmse_limit"]),
    }
    # the analytic reference checked on its own: centred FD of its own primal
    self_check = {}
    for k, name in enumerate(names):
        errors = {}
        for step in (1e-4, 1e-5, 1e-6):
            h = step * abs(props[pidx[k]])
            plus = list(props); plus[pidx[k]] += h
            minus = list(props); minus[pidx[k]] -= h
            fd_a = (A.march(plus, strain_path_of(problem, enriched), [dt] * n, names)["vm"]
                    - A.march(minus, strain_path_of(problem, enriched), [dt] * n, names)["vm"]) / (2 * h)
            errors[f"{step:g}"] = float(nrmse(fd_a[:, None], ref["dvm"][:, k:k + 1])[0])
        self_check[name] = min(errors.values())
    record["analytic_chain_rule"]["self_check_vs_fd_of_analytic_primal_best_nrmse"] = self_check
    print(f"[claim3] sigma_vM final {series['vm'][-1]:.4f} MPa; OTI vs analytic NRMSE worst "
          f"{format_e(err.max())}; max|du/dp| {enriched.max_du_dp:.1e}", flush=True)

    # ---- FD of the ORIGINAL, whole analysis re-solved ------------------------
    started = time.perf_counter()
    fd = fd_reference(mat, props, pidx, problem, names)
    fd_seconds = time.perf_counter() - started
    fd_report = {}
    for k, name in enumerate(names):
        errors = {f"{s:g}": float(nrmse(fd[name][s][:, None], ref["dvm"][:, k:k + 1])[0]) for s in LADDER}
        best = min(errors, key=errors.get)
        fd_report[name] = {"nrmse_vs_analytic_per_step": errors, "best_step": best,
                           "best_nrmse": errors[best],
                           "oti_nrmse": float(err[k]),
                           "oti_vs_fd_best_nrmse": float(nrmse(dvm[:, k:k + 1], fd[name][float(best)][:, None])[0])}
    record["central_fd"] = {"per_parameter": fd_report, "ladder": list(LADDER),
                            "analyses_per_parameter_per_step": 2, "seconds_whole_ladder": fd_seconds,
                            "best_case_nrmse_worst": max(v["best_nrmse"] for v in fd_report.values()),
                            "note": "best step chosen with the analytic value (the 'h optimised' best case)"}

    # ---- weighted sensitivities per regime ----------------------------------
    dgam = GAMMA / n
    weights = np.abs(dvm) * pvals[None, :]
    record["weighted_sensitivities_old_path"] = regimes(ref["eqp"], dgam, names, weights, series["vm"])
    n_fine, dt_fine = n * args.refine, dt / args.refine
    fine = C.simple_shear_problem(1, GAMMA, n_fine, dt_fine)
    fine_sol = C.solve(mat, props, fine, mode="oti", sensitivities=True)
    fine_series = element_series(fine, fine_sol)
    fine_ref = A.march(props, [[0, 0, 0, GAMMA / n_fine, 0, 0]] * n_fine, [dt_fine] * n_fine, names)
    fine_weights = np.abs(fine_series["dvm"]) * pvals[None, :]
    record["weighted_sensitivities_refined_path"] = regimes(
        fine_ref["eqp"], GAMMA / n_fine, names, fine_weights, fine_series["vm"])
    record["weighted_sensitivities_refined_path"]["path"] = f"{n_fine} increments, dt={dt_fine} (same strain rate)"
    record["weighted_sensitivities_refined_path"]["nrmse_vs_analytic_worst"] = float(
        nrmse(fine_series["dvm"], fine_ref["dvm"]).max())

    # ---- larger mesh ----------------------------------------------------------
    started = time.perf_counter()
    mesh_problem = C.simple_shear_problem(args.mesh, GAMMA, n, dt, all_boundary=True)
    mesh_sol = C.solve(mat, props, mesh_problem, mode="oti", sensitivities=True)
    mesh_series = element_series(mesh_problem, mesh_sol)
    record["larger_mesh"] = {
        "mesh": f"{args.mesh}x{args.mesh}x{args.mesh} C3D8 ({len(mesh_problem.conn)} elements, "
                f"{mesh_problem.ndof} dof, {mesh_problem.free.size} free)",
        "boundary_conditions": "affine simple shear u1 = gamma*y prescribed on every boundary node, interior free "
                               "(the study's controlled-deformation test)",
        "max_du_dp": mesh_sol.max_du_dp,
        "dvm_mesh_vs_single_max_rel": float(np.max(np.abs(mesh_series["dvm"] - dvm)) / np.max(np.abs(dvm))),
        "vm_mesh_vs_single_max_rel": float(np.max(np.abs(mesh_series["vm"] - series["vm"])) / np.max(series["vm"])),
        "seconds": time.perf_counter() - started,
        "slide_old_value": SLIDE["mesh_vs_single_old"],
    }
    print(f"[claim3] {record['larger_mesh']['mesh']}: sensitivities vs single element "
          f"{format_e(record['larger_mesh']['dvm_mesh_vs_single_max_rel'])}", flush=True)

    # ---- cost -----------------------------------------------------------------
    record["cost"] = {
        "python_analysis": python_timing(mat, props, n, dt, args.timing_reps),
        "compiled_kernel": kernel_timing(build, props, work, n, dt, args.kernel_reps),
        "fd_accuracy_best_case_vs_hypad": {
            "hypad_nrmse_worst": float(err.max()),
            "fd_best_step_nrmse_worst": record["central_fd"]["best_case_nrmse_worst"]},
        "load_average": os.getloadavg(),
    }
    ck = record["cost"]["compiled_kernel"]["normalised_to_plain"]
    py = record["cost"]["python_analysis"]["normalised_to_plain"]
    print(f"[claim3] cost (x plain): kernel HYPAD EVAL {ck['hypad_eval_carry']:.2f}, MARCH "
          f"{ck['hypad_march']:.2f}; python analysis HYPAD {py['hypad']:.2f}; central FD 13", flush=True)

    # ---- Abaqus -----------------------------------------------------------------
    if args.abaqus:
        from presentation import claim3_abaqus
        abaqus_work = (args.abaqus_work or Path(os.environ.get(
            "PRESENTATION_ABAQUS_WORK", work / "abaqus"))).resolve()
        record["abaqus"] = claim3_abaqus.run(
            props=props, names=names, pidx=pidx, n_increments=n, dt=dt, gamma=GAMMA,
            work=abaqus_work, python_vm=series["vm"], oti_dvm=dvm, analytic_dvm=ref["dvm"])
    else:
        record["abaqus"] = {"status": "not_run (pass --abaqus)"}

    record["series"] = {"gamma_12": (GAMMA * np.arange(1, n + 1) / n).tolist(),
                        "sigma_vm": series["vm"].tolist(), "dvm_dp_oti": dvm.tolist(),
                        "dvm_dp_analytic": ref["dvm"].tolist(), "eqp": ref["eqp"].tolist()}
    record["provenance"] = provenance(); record["imports"] = origins
    out = args.out.resolve()
    write_json(out / "claim3_cp_residual_c3d8.json", record)
    print(f"wrote {out / 'claim3_cp_residual_c3d8.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
