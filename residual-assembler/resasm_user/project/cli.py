"""``resasm project ...`` -- the project-oriented CLI (M5-A).

    resasm project create  DIR --inp model.inp --umat umat.for [--name N]
    resasm project inspect DIR
    resasm project configure DIR --param C11:1 --param C12:2 --param C44:3
    resasm project prepare DIR
    resasm project run     DIR          # run -> export -> assemble -> validate -> report
    resasm project validate DIR
    resasm project report  DIR
    resasm project status  DIR

Every subcommand is a thin wrapper over the reusable action engine; all the work
lives in idempotent, resumable actions (see resasm_user/project/*).
"""

from __future__ import annotations

import argparse
import sys
from typing import List

from . import (Project, ProjectError, State, run_action, run_pipeline,
               set_parameter_selection)


def _print_results(results) -> bool:
    ok = True
    for r in results:
        mark = {"complete": "OK ", "skipped": "-- ", "failed": "!! "}.get(r.status, "?")
        print("  [%s] %-24s %s" % (mark, r.action,
                                   r.message or (r.error or "")))
        ok = ok and r.ok
    return ok


def _status(project: Project) -> None:
    s = project.summary()
    print("project: %s   state: %s" % (s["name"], s["state"]))
    print("  root      : %s" % s["root"])
    print("  inputs    : %s" % ", ".join("%s=%s" % kv for kv in s["inputs"].items()))
    print("  parameters: %s" % (", ".join(s["parameters"]) or "(none configured)"))
    if s["actions"]:
        print("  actions   :")
        for k, v in s["actions"].items():
            print("      %-24s %s" % (k, v))


def main(argv: List[str] = None) -> int:
    p = argparse.ArgumentParser(prog="resasm project",
                                description="Project-oriented sensitivity workflow")
    sub = p.add_subparsers(dest="cmd", required=True)

    c = sub.add_parser("create", help="create a project from an .inp (+ UMAT)")
    c.add_argument("dir")
    c.add_argument("--inp", required=True)
    c.add_argument("--umat")
    c.add_argument("--name")

    for name, help_ in (("inspect", "scan model, UMAT and toolchain"),
                        ("prepare", "generate OTI UMAT, SDV layout, env, patched .inp"),
                        ("run", "run Abaqus -> export -> assemble -> validate -> report"),
                        ("validate", "run validation checks"),
                        ("report", "write the results report"),
                        ("status", "show project state and action status")):
        s = sub.add_parser(name, help=help_)
        s.add_argument("dir")

    cfg = sub.add_parser("configure", help="select parameters (NAME:PROPS_INDEX)")
    cfg.add_argument("dir")
    cfg.add_argument("--param", action="append", default=[], metavar="NAME:INDEX",
                     help="repeatable, e.g. --param C11:1 --param C12:2")

    args = p.parse_args(argv)

    try:
        if args.cmd == "create":
            proj = Project.create(args.dir, name=(args.name or _basename(args.dir)),
                                  inp=args.inp, umat=args.umat)
            print("created project '%s' at %s" % (proj.data["name"], proj.root))
            return 0

        proj = Project.load(args.dir)

        if args.cmd == "status":
            _status(proj)
            return 0

        if args.cmd == "inspect":
            ok = _print_results(run_pipeline(proj, State.INSPECTED))
            _status(proj)
            return 0 if ok else 1

        if args.cmd == "configure":
            sel = []
            for tok in args.param:
                if ":" not in tok:
                    print("ERROR: --param must be NAME:INDEX, got %r" % tok,
                          file=sys.stderr)
                    return 2
                nm, ix = tok.split(":", 1)
                sel.append({"name": nm.strip(), "index": int(ix)})
            set_parameter_selection(proj, sel)
            ok = _print_results([run_action(proj, "configure_parameters")])
            return 0 if ok else 1

        if args.cmd == "prepare":
            ok = _print_results(run_pipeline(proj, State.PREPARED))
            return 0 if ok else 1

        if args.cmd == "run":
            ok = _print_results(run_pipeline(proj, State.REPORTED))
            _status(proj)
            return 0 if ok else 1

        if args.cmd == "validate":
            ok = _print_results([run_action(proj, "run_validation")])
            return 0 if ok else 1

        if args.cmd == "report":
            ok = _print_results([run_action(proj, "generate_report")])
            return 0 if ok else 1

    except ProjectError as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        return 2
    return 0


def _basename(path: str) -> str:
    import os
    return os.path.basename(os.path.normpath(path)) or "project"


if __name__ == "__main__":
    raise SystemExit(main())
