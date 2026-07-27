#!/usr/bin/env python3
"""Render a terminal still of the guided sensitivity-job flow (init -> validate
-> run) with actionable validation errors, for the deck / docs.

    python results/demos/still_guided_job.py    # -> results/figures/guided_job_cli.png
"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from termcast import still, prompt

RA = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
OUT = os.path.join(RA, "results", "figures", "guided_job_cli.png")
P = prompt(path="~/study")

LINES = [
    P + [("resasm-job init job.json", "cmd")],
    [("wrote a starter sensitivity job -> ", "fg"), ("job.json", "blue")],
    [("  edit the material / analysis paths, parameters and outputs, then validate", "gray")],
    [],
    P + [("resasm-job validate job.json", "cmd")],
    [("3 problem(s) in job.json:", "fg")],
    [("  - material.oti_umat not found: material/umat_oti.obj", "red")],
    [("  - parameter ", "red"), ("'YOung'", "yellow"),
     (" is not in the material contract; available: E, nu", "red")],
    [("  - outputs[1].component ", "red"), ("'RF9'", "yellow"),
     (" invalid for reaction_force; use RF1, RF2, RF3", "red")],
    [],
    [("# every problem at once, in plain language — never a raw KeyError or Fortran dump", "comment")],
    P + [("resasm-job validate job.json", "cmd")],
    [("OK: job.json is a valid sensitivity job", "green")],
    [],
    P + [("resasm-job run job.json", "cmd")],
    [("  replay OTI UMAT  ", "fg"), ("ok", "green"),
     ("   assemble + solve  ", "fg"), ("ok", "green"),
     ("   -> ", "fg"), ("sensitivity_out/results.json", "blue")],
]


def main():
    still(OUT, "resasm — guided sensitivity job (validate before you run)", LINES)
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
