#!/usr/bin/env python3
"""Regenerate the committed history-replay example from real Abaqus runs.

Needs a licensed Abaqus (2021 or newer) with a Fortran compiler for user
subroutines, and the UMAT-OTI repository for the ORIGINAL m3_j2 UMAT source.
Runs one job at a time in ``--work`` (never inside the repository):

1. writes the deck (12 x 4 x 2 C3D8 J2 beam, root clamped, tip pushed
   0.08 mm in 10 fixed increments, tight Abaqus convergence controls);
2. runs the nominal analysis and exports its ODB with
   ``residual_core/replay/odb_export_npz.py`` (Abaqus Python);
3. reruns the analysis with each PROPS entry scaled by (1 +/- h), h = 1e-2,
   5e-3, 2e-3 and 1e-3, exports each, and forms central differences of
   selected outputs (an existing export in --work is reused);
4. writes into ``examples/replay_history/j2_beam/``: ``Analysis.inp``,
   ``fields.npz`` (float32: exactly the single-precision ODB values) and
   ``abaqus_fd.json``.

    python examples/replay_history/regenerate_example.py --work /scratch/dir

A job counts as completed only when its .sta file says THE ANALYSIS HAS
COMPLETED SUCCESSFULLY (Abaqus 2021.HF5 may abort in teardown afterwards).
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE))
from make_beam_deck import deck  # noqa: E402

N, PUSH, STEPS = (12, 4, 2), 0.08, 10
PROPS = [200000.0, 0.3, 250.0, 2000.0]
NAMES = ["E", "nu", "SIGY0", "H"]
STEPS_FD = (0.01, 0.005, 0.002, 0.001)
PREFIX = "beam_j2"


def run_job(work, job, text, umat, abaqus):
    export = work / (job + "_fields.npz")
    if export.is_file() and (work / (job + ".inp")).is_file() and (work / (job + ".inp")).read_text() == text:
        return export, 0.0
    (work / (job + ".inp")).write_text(text)
    started = time.time()
    subprocess.run([abaqus, "job=" + job, "input=" + job + ".inp", "user=" + str(umat), "double=both",
                    "cpus=1", "interactive"], cwd=work, capture_output=True, text=True)
    sta = work / (job + ".sta")
    if not sta.is_file() or "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" not in sta.read_text():
        raise SystemExit("%s did not complete (see %s)" % (job, sta))
    export = work / (job + "_fields.npz")
    completed = subprocess.run([abaqus, "python", str(REPO / "residual_core/replay/odb_export_npz.py"), "--",
                                str(work / (job + ".odb")), str(export)], cwd=work, capture_output=True, text=True)
    if completed.returncode or not export.is_file():
        raise SystemExit("export of %s failed:\n%s%s" % (job, completed.stdout, completed.stderr))
    return export, time.time() - started


def outputs(npz):
    d = np.load(npz)
    coords, conn = d["coords"], d["conn"]
    labels = list(d["node_labels"])
    nx, ny, nz = N
    tip = np.isclose(coords[:, 0], nx)
    mid_top = labels.index(1 + nx // 2 + (nx + 1) * ny)            # x = L/2, y = H, z = 0
    tip_top = labels.index(1 + nx + (nx + 1) * ny)                  # x = L, y = H, z = 0
    return {"tip_RF2": d["RF"][:, tip, 1].sum(axis=1),
            "midtop_U2": d["U"][:, mid_top, 1], "tiptop_U1": d["U"][:, tip_top, 0],
            "e1_ip1_S11": d["S"][:, 0, 0, 0], "e1_ip1_S12": d["S"][:, 0, 0, 3],
            "e1_ip1_SDV1": d["SDV"][:, 0, 0, 0]}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", required=True, type=Path)
    ap.add_argument("--umat", type=Path, help="ORIGINAL m3_j2 umat.for (default: $UMAT_OTI_REPO/...)")
    ap.add_argument("--abaqus", default="abaqus")
    a = ap.parse_args(argv)
    work = a.work.resolve()
    if REPO in work.parents or work == REPO:
        raise SystemExit("--work must be outside the repository")
    work.mkdir(parents=True, exist_ok=True)
    umat = a.umat or Path(os.environ.get("UMAT_OTI_REPO", "")) / "parameter_sensitivity/models/m3_j2/umat.for"
    if not umat.is_file():
        raise SystemExit("ORIGINAL UMAT not found: %s (set UMAT_OTI_REPO or --umat)" % umat)
    shutil.copy2(umat, work / "umat.for")
    target = HERE / "j2_beam"
    target.mkdir(exist_ok=True)
    nominal_text = deck(N, PUSH, STEPS, PROPS, 1, tight=True)
    export, seconds = run_job(work, PREFIX + "_nominal", nominal_text, work / "umat.for", a.abaqus)
    print("nominal: %.1f s" % seconds)
    (target / "Analysis.inp").write_text(nominal_text)
    data = dict(np.load(export))
    compact = {}
    for key, value in data.items():
        if value.dtype == np.float64:
            single = value.astype(np.float32)
            if not np.array_equal(single.astype(np.float64), value):
                raise SystemExit("%s is not exactly single precision; refusing to truncate" % key)
            value = single
        compact[key] = value
    np.savez_compressed(target / "fields.npz", **compact)
    base = outputs(export)
    reference = {"props": PROPS, "parameters": NAMES, "steps": list(STEPS_FD),
                 "note": ("central differences (q(p(1+h)) - q(p(1-h))) / (2 h p) of Abaqus reruns with "
                          "tight controls (R_n 1e-9, C_n 1e-8); ODB field output is single precision"),
                 "nominal": {key: value.tolist() for key, value in base.items()}, "fd": {}}
    for index, name in enumerate(NAMES):
        reference["fd"][name] = {}
        for h in STEPS_FD:
            runs = []
            for sign in (+1, -1):
                props = list(PROPS)
                props[index] = PROPS[index] * (1 + sign * h)
                job = "%s_%s_%s%s" % (PREFIX, name, "p" if sign > 0 else "m", str(h).replace(".", "p"))
                export_fd, seconds = run_job(work, job, deck(N, PUSH, STEPS, props, 1, tight=True),
                                             work / "umat.for", a.abaqus)
                print("%s: %.1f s" % (job, seconds))
                runs.append(outputs(export_fd))
            reference["fd"][name][str(h)] = {key: ((runs[0][key] - runs[1][key]) / (2 * h * PROPS[index])).tolist()
                                             for key in base}
    (target / "abaqus_fd.json").write_text(json.dumps(reference, indent=1) + "\n")
    print("wrote", target)


if __name__ == "__main__":
    main()
