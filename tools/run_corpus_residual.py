#!/usr/bin/env python3
"""Drive verified corpus UMATs through RA's residual assembly and sensitivities.

Every case runs in its own subprocess: a UMAT that calls XIT/STOP, or prints on
every call, takes down or floods only its own process, whose stdout is kept in
``<out>/logs/<key>.log``. Records (one per UMAT x feature x problem) are
appended to ``<out>/records.jsonl``; evidence JSON under ``<out>/evidence``.

    PYTHONPATH=/home/ammslab3/softwarex_work/final-umat/src \\
    UMAT_OTI_REPO=/home/ammslab3/softwarex_work/final-umat \\
    python tools/run_corpus_residual.py --keys c7bf17b21519e33da0b7bbb1 \\
        --out /home/ammslab3/softwarex_work/corpus_campaign/batches/B1/noether

    --control sweep_j2_bilinear   a labelled NON-corpus plastic control
    --quick                       single element, 3 increments
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def _one(args) -> int:
    from residual_core.corpus.runner import control_case, run_case
    from residual_core.corpus.sources import load_case
    case = control_case(args.one[len("control_"):]) if args.one.startswith("control_") \
        else load_case(args.one)
    features = tuple(args.features.split(","))
    records = run_case(case, Path(args.out), quick=args.quick, features=features,
                       kinematics=None if args.kinematics == "default" or not case.finite
                       else args.kinematics)
    with open(args.records_out, "w", encoding="utf-8") as stream:
        for record in records:
            stream.write(json.dumps(record, default=str) + "\n")
    return 0


def _fmt(value):
    return "-" if value is None else ("%.1e" % value if isinstance(value, float) else str(value))


def summarise_ledger(path: Path) -> str:
    """Markdown: one row per (UMAT, problem), the three features side by side."""
    rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    table = {}
    for r in rows:
        key = (r["key"], r.get("problem", "-"))
        table.setdefault(key, {"source": r.get("source_id", ""), "family": r.get("family", ""),
                               "kin": r.get("kinematics", "")})[r["feature"]] = r
    out = ["| key | source | family | kin | problem | residual_sens (max_rel, n) | "
           "global_sens (max_rel, n) | QoI | assembly (K_exact, Newton order, RA K author / OTI DDSDDE) |",
           "|---|---|---|---|---|---|---|---|---|"]
    for (key, problem), cells in table.items():
        def cell(feature):
            r = cells.get(feature)
            if r is None:
                return "-"
            if r["status"] in ("verified", "failed") and r.get("counts"):
                c = r["counts"]
                extra = "".join(", %d %s" % (c[k], k) for k in ("nonsmooth", "unresolved") if c.get(k))
                return "%s (%s, %d/%d%s)" % (r["status"], _fmt(r.get("max_rel")), c["verified"],
                                              r.get("comparisons", 0), extra)
            return "%s%s" % (r["status"], " [%s]" % r.get("failure_class") if r.get("failure_class") else "")
        g = cells.get("global_sens", {})
        a = cells.get("assembly_consistency", {}).get("checks") or {}
        assembly = "-" if not a else "%s; %s, %s; %s / %s" % (
            cells["assembly_consistency"]["status"], a.get("K_exact_oti"),
            _fmt(a.get("newton_order")), _fmt(a.get("K_ra_ddsdde_max_rel")),
            _fmt(a.get("K_ra_with_oti_ddsdde_max_rel")))
        out.append("| %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
            key[:12], cells["source"][:48], cells["family"], cells["kin"], problem,
            cell("residual_sens"), cell("global_sens"),
            "%s (%s)" % (g.get("qoi_status", "-"), _fmt(g.get("qoi_max_rel"))), assembly))
    return "\n".join(out)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--keys", nargs="*", default=[])
    parser.add_argument("--control", nargs="*", default=[])
    parser.add_argument("--all-verified", action="store_true")
    parser.add_argument("--out", required=True)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--features", default="residual_sens,global_sens,assembly")
    parser.add_argument("--timeout", type=float, default=3600.0)
    parser.add_argument("--kinematics", default="default",
                        help="finite-strain contract: default (= nlgeom: Abaqus NLGEOM, C3D8 "
                             "selectively reduced, Hughes-Winget DROT, pre-rotated STRESS/STRAN), "
                             "nlgeom_mean, nlgeom_full, or b1_legacy (DROT=I, full integration). "
                             "Small-strain cases ignore it.")
    parser.add_argument("--summarise", action="store_true",
                        help="print a markdown table of <out>/records.jsonl and exit")
    parser.add_argument("--one", help=argparse.SUPPRESS)
    parser.add_argument("--records-out", help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if args.one:
        return _one(args)
    if args.summarise:
        print(summarise_ledger(Path(args.out) / "records.jsonl"))
        return 0
    out = Path(args.out)
    (out / "logs").mkdir(parents=True, exist_ok=True)
    keys = list(args.keys) + ["control_" + c for c in args.control]
    if args.all_verified:
        from residual_core.corpus.sources import verified_keys
        keys += [k for k in verified_keys() if k not in keys]
    ledger = out / "records.jsonl"
    for key in keys:
        partial = out / "logs" / ("%s.records.jsonl" % key)
        log = out / "logs" / ("%s.log" % key)
        command = [sys.executable, __file__, "--one", key, "--out", str(out),
                   "--records-out", str(partial), "--features", args.features,
                   "--kinematics", args.kinematics]
        if args.quick:
            command.append("--quick")
        started = time.time()
        with open(log, "w", encoding="utf-8") as stream:
            try:
                code = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT,
                                      timeout=args.timeout, env=os.environ.copy()).returncode
            except subprocess.TimeoutExpired:
                code = "timeout"
        elapsed = time.time() - started
        if code == 0 and partial.is_file():
            rows = [json.loads(line) for line in partial.read_text().splitlines() if line.strip()]
        else:
            rows = [{"schema": "ra-corpus-residual/2", "key": key, "feature": feature,
                     "status": "failed", "failure_class": "case_process_died",
                     "reason": "subprocess exit %s; see %s" % (code, log)}
                    for feature in ("residual_sens", "global_sens")]
        for row in rows:
            row["seconds_case"] = round(elapsed, 1)
        with open(ledger, "a", encoding="utf-8") as stream:
            for row in rows:
                stream.write(json.dumps(row, default=str) + "\n")
        print("%s  %s  (%.0fs)" % (key, ", ".join("%s/%s=%s" % (r.get("problem", "-"),
              r["feature"], r["status"]) for r in rows), elapsed), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
