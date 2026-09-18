#!/usr/bin/env python3
"""Slides 26-27: HYPAD vs the hand-coded analytical Jacobian of a crystal-plasticity routine.

What the slides say
    "The OTI-extracted derivatives match the hand-coded analytical derivatives
    with a max rel error of 10E-15. The FD maximum error is about 10E7 times
    larger than the OTI maximum error. FD takes about 9.00x more time than OTI
    (getting the correct perturbation)."

Which routine and study (recovered, see docs/PRESENTATION_CLAIMS.md)
    The "crystal plasticity UMAT" box is the stand-alone flow-rule subroutine
    ``computeFlowRule`` of the UTSA constitutive driver (12 FCC slip systems,
    thermally activated slip, T = 1123 K), kept in the git repository
    https://github.com/santiagarcia/OTI_computeflowrule (local clone
    ~/Downloads/driver_utsa, HEAD fde918b). ``computeFlowRule.f90`` carries the
    hand-coded analytical derivatives d(dgamma)/d(tau, backstress, thermal SSD,
    cross-slip, athermal SSD); ``computeFlowRule_otis.f90`` is its OTI version
    (module OTIM6N1, five seeded directions). Those sources are NOT part of this
    repository (their licence is not stated here); this script takes them from
    ``--source-dir``, records their SHA-256, and builds them. The benchmark that
    produced the slide numbers compared OTI with analytical only; no finite-
    difference script or timing record for the FD statements exists anywhere on
    this machine, so the FD part is measured here from scratch.

What is measured
    * OTI vs analytical: max over active slip systems and variables of
      |OTI - A| / |A|.
    * Centred, forward and backward FD of the ANALYTICAL routine's primal output
      (one slip system perturbed at a time), over a 7-step ladder
      h = s * max(|z|, |tau_alpha|), s = 1e-2 ... 1e-8; the error of every step
      against the analytical value, the best step (the "h optimised" best case)
      and the step a plateau rule picks without knowing the answer.
    * Timing, same compiler and flags for every method, R repetitions of N calls
      after a warm-up call: analytical (with its Jacobian), OTI, primal only,
      centred FD with a known step (10 primal calls), forward FD with a known
      step (6 calls), centred FD including the search over the 7-step ladder.
      Built with ifort -O2 (the flags of the original benchmark) and gfortran -O2.

    python presentation/claim2_flowrule_jacobian.py --source-dir ~/Downloads/driver_utsa
"""

from __future__ import annotations

import argparse
import csv
import io
import os
import shutil
import statistics
import subprocess
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from presentation.common import (  # noqa: E402
    PRESENTATION, default_out, default_work, format_e, fresh_dir, provenance, sha256, write_json,
)

SOURCES = ("kinds.f90", "master_parameters.f90", "real_utils.f90", "otim6n1.f90",
           "computeFlowRule.f90", "computeFlowRule_otis.f90")
VARIABLES = ("resolved shear stress", "back stress", "thermal SSD resistance",
             "thermal cross-slip resistance", "athermal SSD resistance")
DRIVER = PRESENTATION / "claim2_flowrule_driver.f90"
SLIDE = {"oti_vs_analytical_max_rel": 1e-15, "fd_over_oti_error_ratio": 1e7,
         "fd_over_oti_time": 9.0}
COMPILERS = {"ifort": ["-O2"], "gfortran": ["-O2", "-ffree-line-length-none"]}


def default_source_dir() -> Path:
    return Path(os.environ.get("CLAIM2_FLOWRULE_DIR", Path.home() / "Downloads" / "driver_utsa"))


def source_identity(source_dir: Path) -> Dict[str, object]:
    missing = [name for name in SOURCES if not (source_dir / name).is_file()]
    if missing:
        raise SystemExit(f"flow-rule sources missing from {source_dir}: {missing}")
    commit = ""
    try:
        commit = subprocess.check_output(["git", "-C", str(source_dir), "rev-parse", "HEAD"],
                                         text=True, stderr=subprocess.DEVNULL).strip()
    except (OSError, subprocess.CalledProcessError):
        pass
    return {"directory": str(source_dir), "git_commit": commit,
            "sha256": {name: sha256(source_dir / name) for name in SOURCES}}


def build(compiler: str, source_dir: Path, work: Path) -> Path:
    flags = COMPILERS[compiler]
    directory = fresh_dir(work / compiler)
    objects = []
    for name in SOURCES:
        obj = directory / (Path(name).stem + ".o")
        subprocess.run([compiler, *flags, "-c", str(source_dir / name), "-o", str(obj)],
                       cwd=directory, check=True, capture_output=True, text=True)
        objects.append(str(obj))
    driver_obj = directory / "driver.o"
    subprocess.run([compiler, *flags, "-c", str(DRIVER), "-o", str(driver_obj)],
                   cwd=directory, check=True, capture_output=True, text=True)
    exe = directory / "claim2"
    subprocess.run([compiler, *flags, str(driver_obj), *objects, "-o", str(exe)],
                   cwd=directory, check=True, capture_output=True, text=True)
    return exe


def read_values(exe: Path) -> Dict[str, dict]:
    out = subprocess.run([str(exe), "values"], capture_output=True, text=True, check=True).stdout
    data: Dict[str, dict] = {}
    for row in csv.DictReader(io.StringIO(out)):
        key = row["kind"].strip()
        variable, step = int(row["variable"]), float(row["step"])
        system = int(row["slip_system"]) - 1
        data.setdefault(key, {}).setdefault((variable, step), np.zeros(12))[system] = float(row["value"])
    return data


def accuracy(data: Dict[str, dict]) -> Dict[str, object]:
    primal_a = data["primal_analytical"][(0, 0.0)]
    primal_o = data["primal_oti"][(0, 0.0)]
    active = primal_a != 0.0
    steps = sorted({s for (_v, s) in data["fd_centred"]}, reverse=True)
    per_variable = {}
    oti_worst = 0.0
    for v in range(1, 6):
        A = data["analytical"][(v, 0.0)]
        O = data["oti"][(v, 0.0)]
        mask = active & (A != 0.0)
        rel = np.abs(O[mask] - A[mask]) / np.abs(A[mask])
        oti_rel = float(rel.max()) if rel.size else 0.0
        oti_worst = max(oti_worst, oti_rel)
        entry = {"variable": VARIABLES[v - 1], "active_systems": int(mask.sum()),
                 "oti_vs_analytical_max_rel": oti_rel,
                 "oti_vs_analytical_max_abs": float(np.max(np.abs(O - A)))}
        for kind in ("fd_centred", "fd_forward", "fd_backward"):
            errs = {}
            for s in steps:
                F = data[kind][(v, s)]
                errs[s] = float(np.max(np.abs(F[mask] - A[mask]) / np.abs(A[mask]))) if mask.any() else 0.0
            best = min(errs, key=errs.get)
            # a plateau rule that never looks at the analytical value
            diffs = []
            for coarse, fine in zip(steps, steps[1:]):
                Fc, Ff = data[kind][(v, coarse)][mask], data[kind][(v, fine)][mask]
                size = np.maximum(np.abs(Fc), np.abs(Ff))
                gap = np.where(size > 0.0, np.abs(Fc - Ff) / np.where(size > 0.0, size, 1.0), 0.0)
                diffs.append((float(np.max(gap)) if gap.size else 0.0, coarse))
            plateau_step = min(diffs)[1]
            entry[kind] = {"max_rel_error_per_step": {f"{s:g}": e for s, e in errs.items()},
                           "best_step": best, "best_error": errs[best],
                           "plateau_step": plateau_step, "plateau_error": errs[plateau_step]}
        decoupled = max(float(np.max(np.abs(data["fd_centred_all_systems"][(v, s)][mask]
                                            - data["fd_centred"][(v, s)][mask])
                                     / np.abs(A[mask]))) if mask.any() else 0.0 for s in steps)
        entry["all_systems_vs_one_system_fd_max_rel"] = decoupled
        per_variable[v] = entry
    primal_rel = float(np.max(np.abs(primal_o[active] - primal_a[active]) / np.abs(primal_a[active])))
    smooth = [1, 2, 3, 4]
    summary = {
        "active_slip_systems": int(active.sum()),
        "primal_oti_vs_analytical_max_rel": primal_rel,
        "oti_vs_analytical_max_rel": oti_worst,
        "fd_centred_best_step_max_rel_all_variables": max(per_variable[v]["fd_centred"]["best_error"] for v in range(1, 6)),
        "fd_centred_best_step_max_rel_smooth_variables": max(per_variable[v]["fd_centred"]["best_error"] for v in smooth),
        "fd_forward_best_step_max_rel_all_variables": max(per_variable[v]["fd_forward"]["best_error"] for v in range(1, 6)),
        "fd_centred_plateau_step_max_rel_smooth_variables": max(per_variable[v]["fd_centred"]["plateau_error"] for v in smooth),
    }
    denominator = max(oti_worst, np.finfo(float).eps * 1e-3)
    summary["ratio_fd_centred_best_to_oti_smooth"] = summary["fd_centred_best_step_max_rel_smooth_variables"] / denominator
    summary["ratio_fd_forward_best_to_oti_all"] = summary["fd_forward_best_step_max_rel_all_variables"] / denominator
    summary["ratio_fd_centred_plateau_to_oti_smooth"] = (
        summary["fd_centred_plateau_step_max_rel_smooth_variables"] / denominator)
    summary["kink_note"] = (
        "the athermal SSD resistance is 0 for every slip system and enters through "
        "sqrt(ssd^2 + gnd^2) with gnd = 0, i.e. |ssd|: the hand-coded and the OTI "
        "derivative are the one-sided (ssd > 0) slope, the centred difference averages "
        "the two one-sided slopes; forward FD resolves the one-sided value")
    return {"per_variable": per_variable, "summary": summary}


def timing(exe: Path, calls: int, repetitions: int) -> Dict[str, object]:
    out = subprocess.run([str(exe), "timing", str(calls), str(repetitions)],
                         capture_output=True, text=True, check=True).stdout
    samples: Dict[str, List[float]] = {}
    for line in out.splitlines():
        method, _rep, value = [part.strip() for part in line.split(",")]
        if method != "checksum":
            samples.setdefault(method, []).append(float(value))
    median = {k: statistics.median(v) for k, v in samples.items()}
    ratios = {
        "oti_over_analytical": median["oti"] / median["analytical"],
        "primal_over_analytical": median["primal"] / median["analytical"],
        "fd_centred_known_step_over_oti": median["fd_centred_known_step"] / median["oti"],
        "fd_forward_known_step_over_oti": median["fd_forward_known_step"] / median["oti"],
        "fd_centred_step_search_over_oti": median["fd_centred_step_search"] / median["oti"],
    }
    return {"calls": calls, "repetitions": repetitions,
            "per_call_seconds_median": median, "per_call_seconds_samples": samples,
            "ratios": ratios}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source-dir", type=Path, default=default_source_dir())
    parser.add_argument("--out", type=Path, default=default_out())
    parser.add_argument("--work", type=Path, default=default_work() / "claim2")
    parser.add_argument("--calls", type=int, default=200000)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--compiler", action="append", dest="compilers")
    args = parser.parse_args(argv)
    source_dir = args.source_dir.expanduser().resolve()
    identity = source_identity(source_dir)
    compilers = args.compilers or [c for c in COMPILERS if shutil.which(c)]
    results = {}
    for compiler in compilers:
        print(f"[claim2] building with {compiler} {' '.join(COMPILERS[compiler])}", flush=True)
        exe = build(compiler, source_dir, args.work.resolve())
        acc = accuracy(read_values(exe))
        tim = timing(exe, args.calls, args.repetitions)
        results[compiler] = {"flags": COMPILERS[compiler], "accuracy": acc, "timing": tim,
                             "version": subprocess.run([compiler, "--version"], capture_output=True,
                                                       text=True).stdout.splitlines()[0]}
        s, r = acc["summary"], tim["ratios"]
        print(f"  OTI vs analytical max rel {format_e(s['oti_vs_analytical_max_rel'])}; "
              f"centred FD best step max rel (smooth variables) "
              f"{format_e(s['fd_centred_best_step_max_rel_smooth_variables'])} "
              f"(ratio {format_e(s['ratio_fd_centred_best_to_oti_smooth'])}); forward FD all "
              f"{format_e(s['fd_forward_best_step_max_rel_all_variables'])}", flush=True)
        print(f"  time: OTI/analytical {r['oti_over_analytical']:.2f}, FD centred known step/OTI "
              f"{r['fd_centred_known_step_over_oti']:.2f}, with step search "
              f"{r['fd_centred_step_search_over_oti']:.2f}", flush=True)
    load = os.getloadavg()
    payload = {"claim": "slides 26-27", "slide": SLIDE, "sources": identity,
               "timing_protocol": {
                   "method": "wall clock (system_clock) around N calls, after one untimed call; "
                             "R repetitions; median reported",
                   "same_flags_for_every_method": True, "threads": 1,
                   "load_average_at_end": load},
               "results": results, "provenance": provenance()}
    out = args.out.resolve()
    write_json(out / "claim2_flowrule_jacobian.json", payload)
    print(f"wrote {out / 'claim2_flowrule_jacobian.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
