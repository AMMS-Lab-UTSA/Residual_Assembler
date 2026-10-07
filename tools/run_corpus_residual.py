#!/usr/bin/env python3
"""Drive verified corpus UMATs through RA's residual assembly and sensitivities.

Every case runs in its own subprocess: a UMAT that calls XIT/STOP, or prints on
every call, takes down or floods only its own process, whose stdout is kept in
``<out>/logs/<key>.log``. Records (one per UMAT x feature x problem) are
appended to ``<out>/records.jsonl``; evidence JSON under ``<out>/evidence``.

    export CORPUS_WORKSPACE=/path/to/workspace   # holds final-umat/, corpus_run/,
                                                 # discovery_cache/ (required)
    export UMAT_OTI_REPO=$CORPUS_WORKSPACE/final-umat   # optional; this is the default
    PYTHONPATH=$UMAT_OTI_REPO/src \\
    python tools/run_corpus_residual.py --keys c7bf17b21519e33da0b7bbb1 \\
        --out $CORPUS_WORKSPACE/corpus_campaign/batches/B1/noether

    --routine-verified            every routine-level-verified source (109 of 238)
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


def fold_cell(records):
    """One (source, feature) cell from its per-problem records, by the manifest's
    rule: failed if any problem failed; verified if at least min(2, n) of the n
    EVALUATED problems (not unsupported) verified. States "verified on k of n
    evaluated problems (of r run)", flags a cell that rests on 1 of 3, and carries
    the plant test and the slot coverage:

    * power: a verified problem has power if the same comparisons, judged with the
      analytic derivative scaled by 1 + 1e-4, FAIL. The cell is power-checked if at
      least min(2, n) of its verified problems have power; otherwise its text says
      it is NOT power-checked and on which problems (the verdict on a derivative the
      FD cannot see says nothing).
    * slots: "verified on 7 of 8 slots; PROPS(5) integer, non-differentiable,
      excluded" where the provider was built without a slot."""
    per = {r.get("problem"): r["status"] for r in records}
    verified = [p for p, v in per.items() if v == "verified"]
    evaluated = [p for p, v in per.items() if v != "unsupported"]
    if any(v == "failed" for v in per.values()):
        cell = "failed"
    elif verified and len(verified) >= min(2, len(evaluated)):
        cell = "verified"
    elif set(per.values()) == {"unsupported"}:
        cell = "unsupported"
    else:
        cell = "not_attempted"
    run = len(per)
    thin = cell == "verified" and len(verified) == 1 and run >= 3
    powered = [r.get("problem") for r in records if r["status"] == "verified"
               and (r.get("plant") or {}).get("power")]
    unpowered = [p for p in verified if p not in powered]
    power_checked = cell == "verified" and len(powered) >= min(2, len(evaluated))
    text = "%s: verified on %d of %d evaluated problems (of %d run)" % (
        cell, len(verified), len(evaluated), run)
    if thin:
        text += "; rests on 1 of %d problems" % run
    if cell == "verified":
        text += "; power-checked on %d of %d" % (len(powered), len(verified))
        if not power_checked:
            text += "; NOT POWER-CHECKED (no power at 1e-4 on: %s)" % ", ".join(unpowered)
    unseeded = next((r.get("unseeded_slots") for r in records if r.get("unseeded_slots")), None)
    slots = None
    if unseeded:
        seeded = next(len(r["wrt"].split("[")[1].split("]")[0].split(",")) for r in records
                      if r.get("wrt") and "[" in r["wrt"])
        slots = "verified on %d of %d slots; %s integer, non-differentiable, excluded" % (
            seeded, seeded + len(unseeded),
            ", ".join("PROPS(%d)" % u["props_index"] for u in unseeded))
        text += "; " + slots
    return {"cell": cell, "k": len(verified), "n": len(evaluated), "run": run,
            "rests_on_1_of_3": thin, "power_checked": power_checked,
            "powered_problems": powered, "unpowered_problems": unpowered,
            "slots": slots, "text": text}


def summarise_sources(path: Path) -> str:
    """Per-source cells for residual_sens and global_sens, and the two ways the
    counts must always be quoted: on at least 2 of 3 problems, and on any problem."""
    groups = {}
    for line in path.read_text().splitlines():
        if line.strip():
            r = json.loads(line)
            if r["feature"] in ("residual_sens", "global_sens"):
                groups.setdefault((r["key"], r.get("source_id", "")), {}).setdefault(
                    r["feature"], []).append(r)
    out = ["| key | source | residual_sens | global_sens |", "|---|---|---|---|"]
    tally = {"residual_sens": [0, 0, 0, 0], "global_sens": [0, 0, 0, 0]}
    for (key, source), cells in sorted(groups.items(), key=lambda kv: kv[0][1]):
        folded = {f: fold_cell(cells[f]) for f in tally if f in cells}
        for f, c in folded.items():
            if c["cell"] == "verified":
                tally[f][1] += 1
                tally[f][3] += 1 if c["power_checked"] else 0
                if not c["rests_on_1_of_3"]:
                    tally[f][0] += 1
                    tally[f][2] += 1 if c["power_checked"] else 0
        out.append("| %s | %s | %s | %s |" % (key[:12], source[:60],
                   folded.get("residual_sens", {}).get("text", "-"),
                   folded.get("global_sens", {}).get("text", "-")))
    out.append("")
    for f, (two, anyp, two_pc, any_pc) in tally.items():
        out.append("%s: verified on at least 2 of 3 problems: %d (power-checked: %d); "
                   "on any problem: %d (power-checked: %d) (of %d sources)"
                   % (f, two, two_pc, anyp, any_pc, len(groups)))
    return "\n".join(out)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--keys", nargs="*", default=[])
    parser.add_argument("--control", nargs="*", default=[])
    parser.add_argument("--all-verified", action="store_true",
                        help="every registry record in terminal state fully_verified")
    parser.add_argument("--routine-verified", action="store_true",
                        help="every routine-level-verified source (D-8: eligible, DDSDDE "
                             "verified in the corpus manifest) -- the 109 of 238")
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
    parser.add_argument("--summarise-sources", action="store_true",
                        help="print per-source cells (verified on k of n problems, 1-of-3 flag) "
                             "of <out>/records.jsonl and exit")
    parser.add_argument("--one", help=argparse.SUPPRESS)
    parser.add_argument("--records-out", help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if args.one:
        return _one(args)
    if args.summarise_sources:
        print(summarise_sources(Path(args.out) / "records.jsonl"))
        return 0
    if args.summarise:
        print(summarise_ledger(Path(args.out) / "records.jsonl"))
        return 0
    from residual_core.corpus.sources import CorpusWorkspaceUnset, paths
    try:
        if args.keys or args.all_verified or args.routine_verified:
            paths().registry  # corpus cases need the workspace; say so before any subprocess
        if args.control:
            paths().umat_repo
    except CorpusWorkspaceUnset as exc:
        print("run_corpus_residual: %s" % exc, file=sys.stderr)
        return 2
    out = Path(args.out)
    (out / "logs").mkdir(parents=True, exist_ok=True)
    keys = list(args.keys) + ["control_" + c for c in args.control]
    if args.all_verified:
        from residual_core.corpus.sources import verified_keys
        keys += [k for k in verified_keys() if k not in keys]
    if args.routine_verified:
        from residual_core.corpus.sources import routine_verified_keys
        keys += [k for k in routine_verified_keys() if k not in keys]
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
