"""run_abaqus_validation.py — run an Abaqus job for a verification-zoo case.

Ready to run on a machine with Abaqus + Intel Fortran. Skips cleanly (exit 0,
prints "Abaqus not available: validation pending") when Abaqus is absent, so it
never fails the offline suite.

Usage (with Abaqus):
    python scripts/run_abaqus_validation.py --job my_job --inp model.inp --user umat.for
    python scripts/run_abaqus_validation.py --job my_job --inp model.inp --user uel.for --unsymm

Steps performed (Step 1 of the roadmap):
    1. locate the abaqus executable
    2. launch: abaqus interactive double job=<job> input=<inp> user=<sub>
    3. report the job status; the ODB is then consumed by extract_odb_fields.py
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(__file__))
from _abaqus_env import abaqus_available, note_pending  # noqa: E402


def build_command(args) -> list:
    exe = os.environ.get("ABAQUS_CMD") or "abaqus"
    cmd = [exe, "interactive", "double", "job=%s" % args.job, "input=%s" % args.inp]
    if args.user:
        cmd.append("user=%s" % args.user)
    if args.unsymm:
        cmd.append("unsymm=yes")
    if args.cpus and args.cpus > 1:
        cmd.append("cpus=%d" % args.cpus)
    if args.extra:
        cmd.extend(args.extra)
    return cmd


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Run an Abaqus job (skips w/o Abaqus)")
    p.add_argument("--job", required=True)
    p.add_argument("--inp", required=True)
    p.add_argument("--user", help="UMAT/UEL source (.for/.f)")
    p.add_argument("--unsymm", action="store_true", help="unsymmetric solver (F-bar UELs)")
    p.add_argument("--cpus", type=int, default=1)
    p.add_argument("--extra", nargs="*", help="extra abaqus args")
    p.add_argument("--dry-run", action="store_true", help="print the command and exit")
    args = p.parse_args(argv)

    if not os.path.exists(args.inp):
        note_pending("input file not found: %s" % args.inp)

    cmd = build_command(args)
    if args.dry_run:
        print("would run:", " ".join(cmd))
        return 0
    if not abaqus_available():
        note_pending("would run: %s" % " ".join(cmd))

    print("running:", " ".join(cmd))
    try:
        rc = subprocess.call(cmd)
    except FileNotFoundError:
        note_pending("abaqus executable not launchable")
        return 0
    if rc != 0:
        print("Abaqus job returned non-zero exit code %d" % rc, file=sys.stderr)
    else:
        print("Abaqus job '%s' completed; ODB: %s.odb" % (args.job, args.job))
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
