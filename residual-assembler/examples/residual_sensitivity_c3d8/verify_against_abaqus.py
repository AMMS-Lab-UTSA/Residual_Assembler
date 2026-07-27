# -*- coding: utf-8 -*-
"""Verify the residual-method sensitivity engine against a live Abaqus solve.

The engine returns the EXACT analytic derivative; Abaqus is finite-differenced
across two solves, so a match to O(dE^2) (~1e-5 with dE/E=1e-3) confirms the
engine is right (the residual is the finite-difference truncation, not engine error).

Two cases, covering both output kinds:
  A  load control       -> nonzero du/dE           (the sensitivity solve matters)
  B  displacement control -> nonzero d sigma11 / dE (the stress-output chain rule)

Usage (needs Abaqus for the runs):
  python verify_against_abaqus.py gen        # write Case A inps + engine prediction
  <run Abaqus at E-,E+ ; extract fields>     # see run_and_compare() / README
  python verify_against_abaqus.py compare
  python verify_against_abaqus.py genB / compareB   # Case B
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

from examples.residual_sensitivity_c3d8 import sensitivity_engine as eng

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.join(HERE, "abaqus_runs")
os.makedirs(RUN, exist_ok=True)
E0, DELTA = 210000.0, 1e-3

INP_A = """*Heading
 elastic C3D8 load control, E={E}
*Node
1, 0.0,0.0,0.0
2, 1.0,0.0,0.0
3, 1.0,1.0,0.0
4, 0.0,1.0,0.0
5, 0.0,0.0,1.0
6, 1.0,0.0,1.0
7, 1.0,1.0,1.0
8, 0.0,1.0,1.0
*Element, type=C3D8, elset=EALL
1, 1,2,3,4,5,6,7,8
*Nset, nset=NX0
1,4,5,8
*Nset, nset=NX1
2,3,6,7
*Nset, nset=NY0
1,2,5,6
*Nset, nset=NZ0
1,2,3,4
*Solid Section, elset=EALL, material=STEEL
*Material, name=STEEL
*Elastic
{E}, 0.3
*Boundary
NX0, XSYMM
NY0, YSYMM
NZ0, ZSYMM
*Step, name=Step-1, nlgeom=NO
*Static
1.0,1.0,1e-05,1.0
{load}
*Output, field
*Node Output
U, RF
*Element Output, directions=YES
S
*End Step
"""
LOAD_A = "*Cload\nNX1, 1, 25.0"
LOAD_B = "*Boundary\nNX1, 1, 1, 0.001"


def _write_inps(prefix, load):
    Es = {"m": E0 * (1 - DELTA), "0": E0, "p": E0 * (1 + DELTA)}
    for tag, E in Es.items():
        with open(os.path.join(RUN, "%s_%s.inp" % (prefix, tag)), "w") as fh:
            fh.write(INP_A.format(E=E, load=load))
    return Es


def gen():
    Es = _write_inps("elasA", LOAD_A)
    fixed = eng.symmetry_bcs()
    fext = np.zeros(24)
    for n in eng.X1:
        fext[eng.dof(n, 0)] = 100.0 / 4
    res = eng.sensitivity_dE(E0, fixed, fext)
    json.dump({"Em": Es["m"], "Ep": Es["p"], "dudE_engine": res["dudE"].tolist()},
              open(os.path.join(RUN, "engineA_pred.json"), "w"), indent=2)
    print("Case A: wrote 3 inps + engineA_pred.json to", RUN)


def genB():
    Es = _write_inps("elasB", LOAD_B)
    fixed = eng.symmetry_bcs()
    for n in eng.X1:
        fixed[eng.dof(n, 0)] = 1.0e-3
    res = eng.sensitivity_dE(E0, fixed, np.zeros(24))
    json.dump({"Em": Es["m"], "Ep": Es["p"], "dS_engine": res["dsig_dE_total"].tolist()},
              open(os.path.join(RUN, "engineB_pred.json"), "w"), indent=2)
    print("Case B: wrote 3 inps + engineB_pred.json to", RUN)


def _disp(fp):
    d = json.load(open(fp))["displacements"]
    U = np.zeros(24)
    for ns, comps in d.items():
        n = int(ns) - 1
        for c, v in enumerate(comps):
            U[eng.dof(n, c)] = v
    return U


def _stress(fp):
    return np.asarray(json.load(open(fp))["stress_ip"]["1"], float)


def compare():
    p = json.load(open(os.path.join(RUN, "engineA_pred.json")))
    dudE_fd = (_disp(os.path.join(RUN, "fieldsA_p.json"))
               - _disp(os.path.join(RUN, "fieldsA_m.json"))) / (p["Ep"] - p["Em"])
    eng_v = np.asarray(p["dudE_engine"])
    print("Case A  du/dE : engine (exact)  vs  Abaqus central-difference")
    worst = 0.0
    for i in range(24):
        if abs(eng_v[i]) < 1e-14 and abs(dudE_fd[i]) < 1e-14:
            continue
        rel = abs(eng_v[i] - dudE_fd[i]) / max(abs(dudE_fd[i]), 1e-30)
        worst = max(worst, rel)
        print("  node %d comp %d : engine %+.6e  Abaqus %+.6e  rel %.2e"
              % (i // 3 + 1, i % 3 + 1, eng_v[i], dudE_fd[i], rel))
    print("  worst rel err (nonzero DOFs): %.2e  => %s"
          % (worst, "MATCH" if worst < 5e-3 else "CHECK"))


def compareB():
    p = json.load(open(os.path.join(RUN, "engineB_pred.json")))
    dS_fd = (_stress(os.path.join(RUN, "fieldsB_p.json"))
             - _stress(os.path.join(RUN, "fieldsB_m.json"))) / (p["Ep"] - p["Em"])
    eng_v = np.asarray(p["dS_engine"])
    print("Case B  dSigma11/dE : engine (exact)  vs  Abaqus central-difference")
    worst = 0.0
    for ip in range(8):
        rel = abs(eng_v[ip, 0] - dS_fd[ip, 0]) / max(abs(dS_fd[ip, 0]), 1e-30)
        worst = max(worst, rel)
        print("  IP %d : engine %+.6e  Abaqus %+.6e  rel %.2e"
              % (ip + 1, eng_v[ip, 0], dS_fd[ip, 0], rel))
    print("  worst rel err: %.2e  => %s" % (worst, "MATCH" if worst < 5e-3 else "CHECK"))


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "compare"
    {"gen": gen, "compare": compare, "genB": genB, "compareB": compareB}[cmd]()
