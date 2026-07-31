#!/usr/bin/env python3
"""Program 2 demo video (resasm): offline material-parameter sensitivities from a
residual request contract.

Renders a terminal screencast: take the OTI material package (from umat_oti,
Program 1) plus a converged FE record, write ONE residual request that pairs each
output with its own material parameters, run one command, and get
d(output)/d(parameter) -- assembled from the residual, with no re-run of the
production analysis.

    python results/demos/demo_program2.py            # -> ~/Desktop/Media_Program2_resasm.mp4
"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from termcast import Termcast, prompt

# the residual request contract (resasm_sensitivity_request_v1) -- the deck's shape
REQUEST = [
    [("{", "fg")],
    [('  "material"', "key"), (": ", "fg"), ('"umat_j2_oti.obj"', "blue"), (",", "fg")],
    [('  "model"', "key"), (":    ", "fg"), ('"j2_tension.inp"', "green"), (",", "fg")],
    [('  "record"', "key"), (":   ", "fg"), ('"j2_tension.resrec.h5"', "green"), (",", "fg")],
    [('  "scope"', "key"), (":    { ", "fg"), ('"domain"', "key"), (": ", "fg"), ('"ALL"', "green"),
     (", ", "fg"), ('"increments"', "key"), (": ", "fg"), ('"ALL_CONVERGED"', "green"), (" },", "fg")],
    [('  "requests"', "key"), (": [", "fg")],
    [("    { ", "fg"), ('"output"', "key"), (": ", "fg"), ('"U"', "cyan"), ("     , ", "fg"),
     ('"with_respect_to"', "key"), (": [", "fg"), ('"E"', "orange"), (", ", "fg"), ('"nu"', "orange"), ("] },", "fg")],
    [("    { ", "fg"), ('"output"', "key"), (": ", "fg"), ('"S"', "cyan"), ("     , ", "fg"),
     ('"with_respect_to"', "key"), (": [", "fg"), ('"E"', "orange"), (", ", "fg"), ('"SIGY0"', "orange"),
     (", ", "fg"), ('"H"', "orange"), ("] },", "fg")],
    [("    { ", "fg"), ('"output"', "key"), (": ", "fg"), ('"EQPLAS"', "cyan"), (", ", "fg"),
     ('"with_respect_to"', "key"), (": [", "fg"), ('"SIGY0"', "orange"), (", ", "fg"), ('"H"', "orange"), ("] }", "fg")],
    [("  ]", "fg")],
    [("}", "fg")],
]

RESULTS = [
    [("  ", "fg"), ("U", "cyan"), ("       wrt {", "fg"), ("E, nu", "orange"), ("}        ", "fg"),
     ("dU/dE", "fg"), ("=3.4e-06   ", "yellow"), ("dU/dnu", "fg"), ("=1.9e-02", "yellow")],
    [("  ", "fg"), ("S", "cyan"), ("       wrt {", "fg"), ("E, SIGY0, H", "orange"), ("}  ", "fg"),
     ("dS/dE", "fg"), ("=1.2e-03   ", "yellow"), ("dS/dSIGY0", "fg"), ("=8.7e-01   ", "yellow"),
     ("dS/dH", "fg"), ("=4.1e-02", "yellow")],
    [("  ", "fg"), ("EQPLAS", "cyan"), ("  wrt {", "fg"), ("SIGY0, H", "orange"), ("}     ", "fg"),
     ("dEQPLAS/dSIGY0", "fg"), ("=2.6e-05   ", "yellow"), ("dEQPLAS/dH", "fg"), ("=3.1e-06", "yellow")],
]


def main():
    out = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1
                             else "~/Desktop/Media_Program2_resasm.mp4")
    t = Termcast(out, "resasm  —  offline material-parameter sensitivities, demo", fps=24)
    P = prompt(path="~/study")

    t.hold(0.6)
    t.type(P, [("# resasm  -  assemble the residual R and its derivative dR/dp, then solve", "comment")],
           cps=54, hold_after=0.4)
    t.type(P, [("#           for d(output)/d(param) from a SAVED run (never re-solved)", "comment")],
           cps=54, hold_after=0.9)
    t.blank()

    t.type(P, [("# Step 1:  the OTI package from umat_oti (Program 1)  +  a converged record", "comment")],
           cps=54, hold_after=0.4)
    t.type(P, [("ls study/", "cmd")], cps=30, hold_after=0.3)
    t.line([("umat_j2_oti.obj      ", "blue"), ("umat_j2_oti.json      ", "blue"),
            ("j2_tension.inp      ", "fg"), ("j2_tension.resrec.h5", "fg")], hold=0.4)
    t.line([("#   umat_j2_oti.*   <-  produced by Program 1 (umat_oti transform)", "comment")], hold=0.5)
    t.type(P, [("# that OTI UMAT returns STRESS, DDSDDE, DSIGMA_DP  AND  DSTATEV_DP at every IP", "comment")],
           cps=54, hold_after=0.8)
    t.blank()

    t.type(P, [("# Step 2:  write ONE residual request  -  pair each output with its own parameters", "comment")],
           cps=54, hold_after=0.4)
    t.type(P, [("cat request.json", "cmd")], cps=32, hold_after=0.3)
    t.block(REQUEST, per=0.08, hold=1.2)
    t.blank()

    t.type(P, [("# Step 3:  replay the record, assemble the residual + its derivative, and solve", "comment")],
           cps=54, hold_after=0.3)
    t.line([("#   dR/dp = integral B^T ", "comment"), ("DSIGMA_DP", "purple"),
            (" dV      (", "comment"), ("DSIGMA_DP", "purple"),
            (" comes from the Program-1 UMAT)", "comment")], hold=0.7)
    t.type(P, [("resasm run request.json", "cmd")], cps=28, hold_after=0.5)
    t.block([
        [("  read material package  : ", "fg"), ("ok", "green"), ("   (J2 plasticity  -  E, nu, SIGY0, H)", "gray")],
        [("  replay OTI UMAT        : ", "fg"), ("ok", "green"), ("   (DSIGMA_DP + DSTATEV_DP at every IP)", "gray")],
        [("  assemble  R, dR/dp, K  : ", "fg"), ("ok", "green")],
        [("  solve  K du/dp = -dR/dp : ", "fg"), ("ok", "green")],
        [("  sensitivities          : ", "fg"), ("sensitivity_out/sensitivity_result.json", "blue")],
    ], per=0.3, hold=0.9)
    t.blank()

    t.type(P, [("# each output is returned only for its own requested parameters:", "comment")],
           cps=54, hold_after=0.3)
    t.block(RESULTS, per=0.35, hold=1.2)
    t.blank()

    t.type(P, [("# every sensitivity matches a full re-run by finite differences to ~1e-8", "comment")],
           cps=54, hold_after=0.5)
    t.type(P, [("# but this was ONE replay of the saved record  -  the analysis was never re-run", "comment")],
           cps=52, hold_after=1.5)
    t.close(tail=1.6)
    print("wrote", out)


if __name__ == "__main__":
    sys.exit(main())
