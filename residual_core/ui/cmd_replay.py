"""Public entry point for the bounded compiled J2 replay workflow."""

import json
from pathlib import Path
import subprocess
import sys

import numpy as np


def register(subparsers):
    parser = subparsers.add_parser("replay", help="compiled J2 C3D8 replay with total-history sensitivities")
    parser.add_argument("record", type=Path)
    parser.add_argument("--object", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--solve", action="store_true", help="generate a synthetic equilibrated record from the model")
    parser.add_argument("--verify", action="store_true", help="independent ORIGINAL whole-path and whole-model FD sweep")
    parser.set_defaults(func=run)


def run(args):
    from ..replay.connected import solve_history
    from ..replay.path_material import PathMaterial
    from ..replay.record import ReplayRecord, RecordError
    from ..replay.verification import verify_connected

    public = args.out / "public"
    private = args.out / "private"
    public.mkdir(parents=True, exist_ok=True)
    private.mkdir(parents=True, exist_ok=True)
    try:
        contract = json.loads(args.contract.read_text())
        material = PathMaterial(str(args.object.resolve()), contract, str(private / "link"))
        record = ReplayRecord.load(str(args.record))
        result = solve_history(record, material, replay=not args.solve)
        verification = verify_connected(record, material, result) if args.verify else {"passed": None, "status": "not_run"}
        summary = {"status": "completed", "verified": bool(args.verify),
                   "sensitivity_semantics": result["sensitivity_semantics"],
                   "parameters": result["parameters"], "increments": len(result["increments"]),
                   "integration_points": len(result["element_ids"]) * 8,
                   "provenance": result["provenance"], "verification": verification,
                   "max_equilibrium_error": max(row["equilibrium_error"] for row in result["increments"])}
        (private / "record.json").write_text(json.dumps(result.pop("record"), indent=2, allow_nan=False) + "\n")
        (private / "result.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
        (public / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
        print(f"replay: {summary['increments']} increments, {summary['integration_points']} IPs; "
              f"total-history du/dp; verified={summary['verified']}")
        print(f"summary: {public / 'summary.json'}")
        return 0
    except (OSError, ValueError, RuntimeError, RecordError, AssertionError, KeyError,
            np.linalg.LinAlgError, subprocess.CalledProcessError) as error:
        message = str(error)
        if isinstance(error, subprocess.CalledProcessError):
            message += "\n" + (error.stderr or "")
        (public / "summary.json").write_text(json.dumps({"status": "failed", "verified": False, "error": message}, indent=2) + "\n")
        print(f"replay failed: {message}", file=sys.stderr)
        return 2