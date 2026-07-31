#!/usr/bin/env python3
"""Program 1 demo video (umat_oti): take a NORMAL UMAT and make it ALSO return
its parameter derivatives (DSIGMA_DP, DSTATEV_DP), driven by the established
transform contract (resasm_umat_transform_v2).

    python results/demos/demo_program1.py            # -> ~/Desktop/Media_Program1_umat_oti.mp4
"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(__file__))
from termcast import Termcast, prompt

# --- the ORIGINAL UMAT the transformer consumes (3D isotropic elasticity,
#     E = PROPS(1), nu = PROPS(2); ntens=6, nprops=2, nstatev=0) --------------
UMAT = r"""      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,
     1 DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,
     2 CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     3 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C     3D isotropic linear elasticity.   E = PROPS(1),   nu = PROPS(2).
      INCLUDE 'ABA_PARAM.INC'
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS),
     1 DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),
     2 TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3),DROT(3,3),
     3 DFGRD0(3,3),DFGRD1(3,3)
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0)
      DIMENSION DDS(6,6), DSTRESS(6)
C
      E   = PROPS(1)
      ENU = PROPS(2)
      ELAM = E*ENU/((ONE+ENU)*(ONE-TWO*ENU))
      EG   = E/(TWO*(ONE+ENU))
C
C     ELASTICITY MATRIX
      DO K1 = 1,3
        DO K2 = 1,3
          DDS(K2,K1) = ELAM
        END DO
        DDS(K1,K1) = ELAM + TWO*EG
      END DO
      DDS(4,4) = EG
      DDS(5,5) = EG
      DDS(6,6) = EG
C
C     STRESS INCREMENT AND UPDATE
      DO K1 = 1,NTENS
        DSTRESS(K1) = ZERO
        DO K2 = 1,NTENS
          DSTRESS(K1) = DSTRESS(K1) + DDS(K1,K2)*DSTRAN(K2)
        END DO
        STRESS(K1) = STRESS(K1) + DSTRESS(K1)
      END DO
C
C     CONSISTENT TANGENT  DDSDDE
      DO K1 = 1,NTENS
        DO K2 = 1,NTENS
          DDSDDE(K1,K2) = DDS(K1,K2)
        END DO
      END DO
      RETURN
      END
""".split("\n")

# UMAT interface variables get distinct colors
ROLE = {"STRESS": "cyan", "STATEV": "cyan", "DDSDDE": "purple", "DSTRESS": "blue",
        "DSTRAN": "orange", "PROPS": "orange"}
_TOK = re.compile(r"[A-Za-z_][A-Za-z0-9_]*|\d+\.?\d*[Dd]?\d*|\s+|[^\w\s]")


def hl(text):
    if text[:1] in "Cc*!":
        return [(text, "comment")]
    spans = []
    for tok in _TOK.findall(text):
        if tok in ROLE:
            spans.append((tok, ROLE[tok]))
        elif re.match(r"\d", tok):
            spans.append((tok, "yellow"))
        else:
            spans.append((tok, "fg"))
    return spans


# fill the transform contract straight from the UMAT interface
MAP = [
    ("dimensions", "ntens 6  nprops 2  nstatev 0", "yellow", "the sizes the UMAT declares"),
    ("parameters", "E = PROPS(1),  nu = PROPS(2)", "orange", "which PROPS to differentiate  (the seed)"),
    ("derivative", "STRESS  ->  DSIGMA_DP", "cyan", "differentiate the stress update"),
    ("history", "STATEV  ->  DSTATEV_DP", "purple", "differentiate the state  (none here: nstatev 0)"),
]

# the ESTABLISHED transform contract (resasm_umat_transform_v2)
CONTRACT = [
    [("{", "fg")],
    [('  "schema"', "key"), (": ", "fg"), ('"resasm_umat_transform_v2"', "green"), (",", "fg")],
    [('  "source"', "key"), (": { ", "fg"), ('"entry_point"', "key"), (": ", "fg"), ('"UMAT"', "green"),
     (", ", "fg"), ('"main_file"', "key"), (": ", "fg"), ('"umat.for"', "green"), (" },", "fg")],
    [('  "kinematics"', "key"), (": ", "fg"), ('"small_strain"', "green"), (",", "fg")],
    [('  "dimensions"', "key"), (": { ", "fg"), ('"ntens"', "key"), (": ", "fg"), ("6", "yellow"),
     (", ", "fg"), ('"nprops"', "key"), (": ", "fg"), ("2", "yellow"), (", ", "fg"),
     ('"nstatev"', "key"), (": ", "fg"), ("0", "yellow"), (" },", "fg")],
    [('  "parameters"', "key"), (": [ {", "fg"), ('"name"', "key"), (": ", "fg"), ('"E"', "orange"),
     (", ", "fg"), ('"props_index"', "key"), (": ", "fg"), ("1", "yellow"), ("},", "fg"),
     (" {", "fg"), ('"name"', "key"), (": ", "fg"), ('"nu"', "orange"), (", ", "fg"),
     ('"props_index"', "key"), (": ", "fg"), ("2", "yellow"), ("} ],", "fg")],
    [('  "derivative"', "key"), (": { ", "fg"), ('"response"', "key"), (": ", "fg"), ('"STRESS"', "cyan"),
     (", ", "fg"), ('"export"', "key"), (": ", "fg"), ('"DSIGMA_DP"', "purple"), (" },", "fg")],
    [('  "history"', "key"), (":    { ", "fg"), ('"state"', "key"), (": ", "fg"), ('"STATEV"', "cyan"),
     (", ", "fg"), ('"export"', "key"), (": ", "fg"), ('"DSTATEV_DP"', "purple"), (" },", "fg")],
    [('  "output"', "key"), (": { ", "fg"), ('"object"', "key"), (": ", "fg"), ('"umat_elastic_oti.obj"', "blue"),
     (", ", "fg"), ('"contract"', "key"), (": ", "fg"), ('"umat_elastic_oti.json"', "blue"), (" }", "fg")],
    [("}", "fg")],
]


def main():
    out = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1
                             else "~/Desktop/Media_Program1_umat_oti.mp4")
    t = Termcast(out, "umat-oti  —  make a UMAT return its parameter derivatives, demo", fps=24)
    P = prompt()

    t.hold(0.6)
    t.type(P, [("# umat_oti  -  make a normal UMAT ALSO return dSTRESS/dp and dSTATEV/dp  (OTI)", "comment")],
           cps=54, hold_after=0.9)
    t.blank()
    t.type(P, [("# Step 1:  read the original UMAT  (its interface variables are highlighted)", "comment")],
           cps=54, hold_after=0.4)
    t.type(P, [("cat -n umat.for", "cmd")], cps=32, hold_after=0.3)
    rows = [[("%4d  " % (i + 1), "gray")] + hl(line) for i, line in enumerate(UMAT) if line != ""]
    t.block(rows, per=0.05, hold=0.9)
    t.type(P, [("# a normal UMAT returns  STRESS,  STATEV  and the tangent  DDSDDE", "comment")],
           cps=54, hold_after=0.9)
    t.blank()

    t.type(P, [("# Step 2:  fill the transform contract from that interface", "comment")], cps=54, hold_after=0.4)
    for label, var, ck, desc in MAP:
        t.line([("#  ", "comment"), ("%-11s" % label, "comment"), ("%-30s" % var, ck),
                (desc, "comment")], hold=0.16)
    t.hold(0.9)
    t.type(P, [("cat transform.json", "cmd")], cps=34, hold_after=0.3)
    t.block(CONTRACT, per=0.09, hold=1.1)
    t.blank()

    t.type(P, [("# Step 3:  one command transforms it", "comment")], cps=54, hold_after=0.4)
    t.type(P, [("umat-oti transform transform.json", "cmd")], cps=28, hold_after=0.5)
    t.block([
        [("  parse UMAT + trace parameter dependencies : ", "fg"), ("ok", "green")],
        [("  promote to OTI, seed PROPS, extract coeffs : ", "fg"), ("ok", "green")],
        [("  object            : ", "fg"), ("umat_elastic_oti.obj", "blue")],
        [("  completed contract : ", "fg"), ("umat_elastic_oti.json", "blue")],
    ], per=0.28, hold=0.9)
    t.blank()

    t.type(P, [("# the OTI UMAT returns all the original outputs, PLUS the derivatives:", "comment")],
           cps=54, hold_after=0.3)
    t.line([("#   a normal UMAT :        ", "comment"), ("STRESS", "cyan"), (" , ", "comment"),
            ("STATEV", "cyan"), (" , ", "comment"), ("DDSDDE", "purple")])
    t.line([("#   the OTI UMAT ALSO :    ", "comment"), ("DSIGMA_DP", "purple"), (" = dSTRESS/dp", "comment"),
            ("    ", "comment"), ("DSTATEV_DP", "purple"), (" = dSTATEV/dp", "comment")], hold=0.6)
    t.line([("      DSIGMA_DP(", "fg"), ("i", "yellow"), (",", "fg"), ("p", "yellow"), (") = ", "fg"),
            ("GETIM", "purple"), ("( ", "fg"), ("STRESS_OTI", "cyan"), ("(", "fg"), ("i", "yellow"),
            ("), ", "fg"), ("p", "yellow"), (" )", "fg")], hold=0.8)
    t.type(P, [("# same inputs, same original outputs  -  just the exact parameter derivatives", "comment")],
           cps=52, hold_after=1.4)
    t.close(tail=1.6)
    print("wrote", out)


if __name__ == "__main__":
    sys.exit(main())
