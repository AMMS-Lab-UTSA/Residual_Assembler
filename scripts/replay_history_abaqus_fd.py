#!/usr/bin/env python3
"""Abaqus finite-difference reruns of a deck for the history-replay evidence.

For every provider parameter (names and PROPS slots from the completed
contract) and every relative step h, the deck is rerun with that PROPS entry
scaled by (1 +/- h). ``--tight`` inserts ``*Controls, parameters=field`` with
R_n = 1e-9 and C_n = 1e-8 after the ``*Static`` data line, so every rerun is
equilibrated far below the single-precision resolution of the ODB (Abaqus's
defaults 5e-3 / 1e-2 make the difference quotient noisy at ~1e-3). Jobs run
one at a time; success is judged from the .sta file; each ODB is exported
with ``residual_core/replay/odb_export_npz.py`` and then deleted.

    python scripts/replay_history_abaqus_fd.py --deck Analysis.inp --umat umat.for \\
        --contract umat_m3_j2_oti.json --steps 5e-3,2e-3,1e-3 --tight --prefix claudeC_j2cant \\
        --work /scratch/j2_fd_tight

Writes ``manifest.tsv`` (job, parameter, sign, step, value, status, seconds).
"""
import argparse
import json
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXPORTER = ROOT / "residual_core" / "replay" / "odb_export_npz.py"


def perturbed(text, slot, factor, tight):
    lines = text.splitlines()
    out, index = [], 0
    while index < len(lines):
        line = lines[index]
        out.append(line)
        keyword = line.strip().lower()
        if keyword.startswith("*user material"):
            data = []
            index += 1
            while index < len(lines) and not lines[index].lstrip().startswith("*"):
                data += [float(x) for x in lines[index].split(",") if x.strip()]
                index += 1
            data[slot - 1] *= factor
            out += [", ".join(repr(v) for v in data[s:s + 8]) for s in range(0, len(data), 8)]
            continue
        if keyword.startswith("*static"):
            index += 1
            out.append(lines[index])
            if tight:
                out += ["*Controls, parameters=field, field=global", "1e-9, 1e-8"]
        index += 1
    result = "\n".join(out) + "\n"
    if tight and result.lower().count("parameters=field") != 1:
        raise SystemExit("the deck already sets *Controls, parameters=field")
    return result


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--deck", required=True, type=Path)
    ap.add_argument("--umat", required=True, type=Path)
    ap.add_argument("--contract", required=True, type=Path)
    ap.add_argument("--steps", default="5e-3,2e-3,1e-3")
    ap.add_argument("--parameters", help="comma list (default: every contract parameter)")
    ap.add_argument("--tight", action="store_true")
    ap.add_argument("--nominal", action="store_true", help="also rerun the unperturbed deck")
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--work", required=True, type=Path)
    ap.add_argument("--abaqus", default="abaqus")
    ap.add_argument("--cpus", default="4")
    a = ap.parse_args(argv)
    a.work.mkdir(parents=True, exist_ok=True)
    contract = json.loads(a.contract.read_text())
    slots = {p["name"]: p["props_index"] for p in contract["parameters"]}
    names = a.parameters.split(",") if a.parameters else list(slots)
    text = a.deck.read_text()
    manifest = a.work / "manifest.tsv"
    if not manifest.is_file():
        manifest.write_text("job\tparameter\tsign\tstep\tfactor\tstatus\tseconds\n")
    done = {line.split("\t")[0] for line in manifest.read_text().splitlines()[1:]}
    jobs = [("%s_nominal" % a.prefix, "nominal", "0", "0", 1.0, None)] if a.nominal else []
    for h in a.steps.split(","):
        for name in names:
            for sign, factor in (("p", 1 + float(h)), ("m", 1 - float(h))):
                jobs.append(("%s_%s_%s%s" % (a.prefix, name, sign, h.replace(".", "p").replace("-", "m")),
                             name, sign, h, factor, slots[name]))
    for job, name, sign, h, factor, slot in jobs:
        if job in done:
            continue
        deck = perturbed(text, slot or 1, factor if slot else 1.0, a.tight)
        (a.work / (job + ".inp")).write_text(deck)
        started = time.time()
        subprocess.run([a.abaqus, "job=" + job, "input=" + job + ".inp", "user=" + str(a.umat.resolve()),
                        "double=both", "cpus=" + a.cpus, "interactive"], cwd=a.work,
                       capture_output=True, text=True)
        sta = a.work / (job + ".sta")
        status = "analysis_failed"
        if sta.is_file() and "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" in sta.read_text():
            export = subprocess.run([a.abaqus, "python", str(EXPORTER), "--", job + ".odb", job + "_fields.npz"],
                                    cwd=a.work, capture_output=True, text=True)
            status = "completed" if export.returncode == 0 else "export_failed"
        for suffix in (".odb", ".com", ".prt", ".odb_f", ".023", ".mdl", ".stt", ".res", ".sim", ".dat"):
            path = a.work / (job + suffix)
            if path.exists():
                path.unlink()
        with manifest.open("a") as stream:
            stream.write("%s\t%s\t%s\t%s\t%r\t%s\t%d\n" % (job, name, sign, h, factor, status,
                                                           time.time() - started))
        print(job, status, "%.0f s" % (time.time() - started), flush=True)


if __name__ == "__main__":
    main()
