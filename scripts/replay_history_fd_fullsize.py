#!/usr/bin/env python3
"""Whole-model ORIGINAL-UMAT finite differences at full model size (evidence).

For the given parameters, re-equilibrates the whole model in Python with the
ORIGINAL UMAT (the regular path of the provider object) at p(1 +/- h) over the
increments of the ODB export, for a ladder of h, and compares with the OTI
sensitivities of the same discrete problem (``history_verify.whole_model_fd``).
Run one process per parameter to use several cores.

    python scripts/replay_history_fd_fullsize.py --deck Analysis.inp --fields fields.npz \\
        --object umat_oti.obj --parameters E --steps 1e-3,1e-4,1e-5 --out fd_E.json
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from residual_core.replay.history import HistoryEngine  # noqa: E402
from residual_core.replay.history_inputs import load_recorded_fields, read_history_model  # noqa: E402
from residual_core.replay.history_material import HistoryMaterial  # noqa: E402
from residual_core.replay.history_verify import summarize_fd, whole_model_fd  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--deck", required=True, type=Path)
    ap.add_argument("--fields", required=True, type=Path)
    ap.add_argument("--object", required=True, type=Path)
    ap.add_argument("--parameters", required=True)
    ap.add_argument("--steps", default="1e-3,1e-4,1e-5")
    ap.add_argument("--rtol", type=float, default=1e-13)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args(argv)
    started = time.perf_counter()
    model = read_history_model(a.deck)
    times = load_recorded_fields(a.fields, model).time
    contract = json.loads(a.object.with_suffix(".json").read_text())
    material = HistoryMaterial(a.object, contract, str(a.out.with_suffix("")) + "_link")
    engine = HistoryEngine(model, material)
    result = whole_model_fd(engine, times, parameters=a.parameters.split(","),
                            steps=tuple(float(x) for x in a.steps.split(",")), rtol=a.rtol)
    result.pop("reference_result")
    result["summary"] = summarize_fd(result)
    result["wall_s"] = time.perf_counter() - started
    a.out.write_text(json.dumps(result, indent=1) + "\n")
    print(json.dumps(result["summary"], indent=1))


if __name__ == "__main__":
    main()
