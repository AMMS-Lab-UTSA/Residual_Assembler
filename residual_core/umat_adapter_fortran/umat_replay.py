#!/usr/bin/env python3
"""
umat_replay.py -- Python 3 orchestrator for the standalone UMAT material-point
replay adapter (Verification Mode 2).

WHAT IT DOES
------------
Marches the Grilli crystal-plasticity UMAT through a deformation-gradient
history ONE INCREMENT AT A TIME (plasticity is history-dependent; you must
never jump to the final step), by driving the compiled Fortran binary
`umat_driver` (see build.sh / build.bat).  It can:

  * --dry-run
        Fabricate a trivial single-IP history (identity, then a small uniaxial
        ramp), feed it through the driver, and print STRESS/STATEV per
        increment.  Proves the end-to-end plumbing with a known answer
        (identity deformation -> zero stress).  Uses the mock UMAT if the real
        driver is not built.

  * --fields fields.json --element E --ip I
        Reconstruct the deformation-gradient history AT one integration point
        from the ODB fields (nodal U + C3D8 shape-function gradients, per
        residual_core/CONTRACT.md sections 2 & 5), replay it, then compare the
        replayed STRESS against the ODB `S` and STATEV against the ODB `SDV`,
        printing per-increment absolute/relative error tables.

DRIVER BINARY
-------------
Looked up (in order):  --driver PATH  ->  build/umat_driver[.exe] (REAL UMAT)
->  build/umat_driver_mock[.exe] (mock).  Use --mock to force the mock, or
--no-run to only emit the driver input (plumbing/inspection, no binary needed).

IMPORTANT HONESTY NOTES
-----------------------
* The real `umat_driver` requires the real UMAT object, which needs Intel ifort
  + Abaqus (MKL) to compile (gfortran cannot build the unmodified sources -- see
  README.md).  Until then, only the MOCK material is available and the
  STRESS/SDV comparison is against the mock, not crystal plasticity.
* Exact replay needs the deformation gradient at EVERY solver increment.  ODB
  field output written every N increments (the example uses frequency=200) is
  too sparse for an exact plasticity replay; frame-to-frame marching is then a
  sub-stepping approximation.  For an exact Mode-2 check, request field output
  every increment (frequency=1) or use fixed time increments.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys

try:
    import numpy as np
except ImportError:
    sys.stderr.write("ERROR: numpy is required for umat_replay.py\n")
    raise

# --------------------------------------------------------------------------
# C3D8 kinematics.  Prefer the shared residual_core module; fall back to an
# inline implementation consistent with CONTRACT.md sections 2 & 4.
# --------------------------------------------------------------------------
try:
    # repo root on path -> import the shared kernel (post-refactor location)
    sys.path.insert(0, os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))))
    from residual_core.formulations.c3d8_kernel import ABAQUS_C3D8_NODES as _NODES  # type: ignore
    from residual_core.formulations.c3d8_kernel import ABAQUS_C3D8_GAUSS as _GAUSS  # type: ignore
    from residual_core.formulations.c3d8_kernel import shape_grad_natural as _shape_grad_natural  # type: ignore
    _GAUSS_PTS = np.asarray(_GAUSS[0], dtype=float)
    _NODES = np.asarray(_NODES, dtype=float)
    _HAVE_C3D8 = True
except Exception:
    _HAVE_C3D8 = False
    # Node natural coords in Abaqus C3D8 connectivity order (CONTRACT sec 2)
    _NODES = np.array([
        [-1, -1, -1], [+1, -1, -1], [+1, +1, -1], [-1, +1, -1],
        [-1, -1, +1], [+1, -1, +1], [+1, +1, +1], [-1, +1, +1],
    ], dtype=float)
    _g = 1.0 / math.sqrt(3.0)
    # Abaqus C3D8 integration-point order (CONTRACT sec 2)
    _GAUSS_PTS = np.array([
        [-_g, -_g, -_g], [+_g, -_g, -_g], [-_g, +_g, -_g], [+_g, +_g, -_g],
        [-_g, -_g, +_g], [+_g, -_g, +_g], [-_g, +_g, +_g], [+_g, +_g, +_g],
    ], dtype=float)

    def _shape_grad_natural(xi):
        """dN_a/dxi_k, shape (8,3), trilinear C3D8."""
        xi = np.asarray(xi, dtype=float)
        g = np.empty((8, 3), dtype=float)
        for a in range(8):
            xa = _NODES[a]
            for k in range(3):
                other = [j for j in range(3) if j != k]
                g[a, k] = 0.125 * xa[k] * (1.0 + xi[other[0]] * xa[other[0]]) \
                                        * (1.0 + xi[other[1]] * xa[other[1]])
        return g


def deformation_gradient(Xe, Ue, xi):
    """
    F(3,3) at natural coord xi for a C3D8 element.
      Xe : (8,3) reference nodal coords
      Ue : (8,3) nodal displacements
      xi : (3,)  natural coord of the integration point
    F = I + dU/dX,  dU_i/dX_j = sum_a U_a[i] * dN_a/dX_j.
    """
    Xe = np.asarray(Xe, dtype=float)
    Ue = np.asarray(Ue, dtype=float)
    dN_dxi = np.asarray(_shape_grad_natural(xi), dtype=float)   # (8,3)
    J = Xe.T @ dN_dxi                       # J[i,j] = dX_i/dxi_j
    Jinv = np.linalg.inv(J)
    dN_dX = dN_dxi @ Jinv                    # (8,3): dN_a/dX_j
    H = Ue.T @ dN_dX                         # H[i,j] = dU_i/dX_j
    return np.eye(3) + H


# --------------------------------------------------------------------------
# Driver input / output
# --------------------------------------------------------------------------
def _fmt_row(vals):
    return " ".join("%.16e" % v for v in vals)


def write_driver_input(path, nstatv, nprops, ntens, props, statev0, increments):
    """
    Emit the plain-text, list-directed input consumed by umat_driver.
    `increments` is a list of dicts with keys:
      noel, npt, kstep, kinc, dtime, temp, dtemp, time (2,), dfgrd0 (3x3), dfgrd1 (3x3)
    """
    props = list(props)
    statev0 = list(statev0)
    assert len(props) == nprops, "props length %d != nprops %d" % (len(props), nprops)
    assert len(statev0) == nstatv, "statev0 length %d != nstatv %d" % (len(statev0), nstatv)
    with open(path, "w") as f:
        f.write("%d %d %d %d\n" % (nstatv, nprops, ntens, len(increments)))
        f.write(_fmt_row(props) + "\n")
        f.write(_fmt_row(statev0) + "\n")
        for inc in increments:
            f.write("%d %d %d %d\n" % (inc["noel"], inc["npt"], inc["kstep"], inc["kinc"]))
            f.write("%.16e %.16e %.16e %.16e %.16e\n" % (
                inc["dtime"], inc["temp"], inc["dtemp"], inc["time"][0], inc["time"][1]))
            F0 = np.asarray(inc["dfgrd0"], dtype=float).reshape(3, 3)
            F1 = np.asarray(inc["dfgrd1"], dtype=float).reshape(3, 3)
            f.write(_fmt_row(F0.reshape(-1)) + "\n")
            f.write(_fmt_row(F1.reshape(-1)) + "\n")


def parse_driver_output(path):
    """Parse umat_driver output into a list of records."""
    with open(path) as f:
        toks = []
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            toks.append(s)
    # first non-comment line: nstatv ntens ninc
    hdr = toks[0].split()
    nstatv, ntens, ninc = int(hdr[0]), int(hdr[1]), int(hdr[2])
    recs = []
    idx = 1
    # values may be wrapped across lines; flatten remaining tokens into a stream
    stream = " ".join(toks[idx:]).split()
    p = 0

    def take(n):
        nonlocal p
        vals = [float(x) for x in stream[p:p + n]]
        p += n
        return vals

    for _ in range(ninc):
        kinc = int(float(stream[p])); p += 1
        stress = take(ntens)
        statev = take(nstatv)
        ddsdde = take(ntens * ntens)
        pnewdt = take(1)[0]
        recs.append(dict(kinc=kinc, stress=stress, statev=statev,
                         ddsdde=ddsdde, pnewdt=pnewdt))
    return recs


# --------------------------------------------------------------------------
# Driver discovery / invocation
# --------------------------------------------------------------------------
def find_driver(explicit, prefer_mock):
    here = os.path.dirname(os.path.abspath(__file__))
    bd = os.path.join(here, "build")
    ext = ".exe" if os.name == "nt" else ""
    real = os.path.join(bd, "umat_driver" + ext)
    mock = os.path.join(bd, "umat_driver_mock" + ext)
    if explicit:
        return explicit if os.path.exists(explicit) else None
    order = [mock, real] if prefer_mock else [real, mock]
    for c in order:
        if os.path.exists(c):
            return c
    return None


def run_driver(driver, infile, outfile):
    subprocess.run([driver, infile, outfile], check=True)
    return parse_driver_output(outfile)


# --------------------------------------------------------------------------
# fields.json helpers
# --------------------------------------------------------------------------
def load_fields(path):
    with open(path) as f:
        return json.load(f)


# History / indexing sanity thresholds (CONTRACT sec 6).
_NELEMENTS = 18315   # nElements in mycommon.f -> /UMPS/ common leading dimension
_JUMP_TOL = 0.5      # max per-increment ||DFGRD1 - DFGRD0||_F before warning


def element_history(fields, eid, ip):
    """
    Build the DFGRD history for element `eid`, integration point `ip` (1..8),
    plus the reference S / SDV series to compare against.
    Returns (increments, ref_S, ref_SDV, frame_labels).
    """
    nodes = {int(k): np.asarray(v, dtype=float) for k, v in fields["nodes"].items()}
    conn = [int(n) for n in fields["elements"][str(eid)]]
    Xe = np.array([nodes[n] for n in conn], dtype=float)      # (8,3) ref coords
    xi = _GAUSS_PTS[ip - 1]

    frames = fields["frames"]
    # sort frames by time to be safe
    frames = sorted(frames, key=lambda fr: fr.get("time", 0.0))

    # F at each output frame (relative to the undeformed reference)
    F_series, times, S_series, SDV_series, labels = [], [], [], [], []
    for fr in frames:
        U = fr["U"]
        Ue = np.array([U[str(n)] for n in conn], dtype=float)  # (8,3)
        F = deformation_gradient(Xe, Ue, xi)
        F_series.append(F)
        times.append(float(fr.get("time", 0.0)))
        labels.append(fr.get("frame", len(labels)))
        S = fr.get("S", {}).get(str(eid))
        S_series.append(np.asarray(S[ip - 1], dtype=float) if S is not None else None)
        SDV = fr.get("SDV", {}).get(str(eid)) if "SDV" in fr else None
        SDV_series.append(np.asarray(SDV[ip - 1], dtype=float) if SDV is not None else None)

    # prepend an undeformed reference state at t=0 (F=I) if the first frame is
    # not already the reference, so increment 1 marches I -> F(frame0).
    incs = []
    F_prev = np.eye(3)
    t_prev = 0.0
    kstep = frames[0].get("step", "Step-1") if frames else "Step-1"
    for k, F in enumerate(F_series):
        t = times[k]
        dt = t - t_prev if (t - t_prev) > 0 else 1.0
        # Abaqus TIME(1)=step time and TIME(2)=total time are BOTH the values at
        # the START of the increment; dt (DTIME) carries the increment length.
        incs.append(dict(noel=int(eid), npt=int(ip), kstep=1, kinc=k + 1,
                         dtime=dt, temp=293.0, dtemp=0.0, time=(t_prev, t_prev),
                         dfgrd0=F_prev, dfgrd1=F))
        F_prev = F
        t_prev = t

    # --- history-jump guards (CONTRACT sec 6: plasticity needs full history) ---
    if len(labels) < 2:
        sys.stderr.write(
            "WARNING: replay history has %d output frame(s) (< 2); a single "
            "increment cannot capture the deformation history that history-"
            "dependent plasticity requires (CONTRACT sec 6).\n" % len(labels))
    max_jump = 0.0
    for inc in incs:
        dF = np.asarray(inc["dfgrd1"], dtype=float) - np.asarray(inc["dfgrd0"], dtype=float)
        max_jump = max(max_jump, float(np.linalg.norm(dF)))   # Frobenius
    if max_jump > _JUMP_TOL:
        sys.stderr.write(
            "WARNING: max per-increment ||DFGRD1 - DFGRD0|| (Frobenius) = %.3e "
            "exceeds %.2f; increments are too coarse for an exact history-"
            "dependent plasticity replay (CONTRACT sec 6).\n"
            % (max_jump, _JUMP_TOL))
    if int(eid) > _NELEMENTS:
        sys.stderr.write(
            "WARNING: element label %d > nElements (%d); would index the static "
            "/UMPS/ common out of range (the real UMAT indexes kFp(noel,...)).\n"
            % (int(eid), _NELEMENTS))

    return incs, S_series, SDV_series, labels


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------
_STRESS_COMPS = ["S11", "S22", "S33", "S12", "S13", "S23"]


def print_stress_table(recs, ref_S, labels):
    tiny = 1e-12
    print("\n==== replayed STRESS vs ODB S (per output frame) ====")
    for k, rec in enumerate(recs):
        if k >= len(ref_S):
            break
        s_rep = np.asarray(rec["stress"], dtype=float)
        s_odb = ref_S[k]
        tag = "frame %s (kinc=%d)" % (labels[k], rec["kinc"])
        if s_odb is None:
            print("  %s : no ODB S; replayed = [%s]" %
                  (tag, ", ".join("%.4e" % v for v in s_rep)))
            continue
        aerr = np.abs(s_rep - s_odb)
        rerr = aerr / np.maximum(np.abs(s_odb), tiny)
        print("  %s" % tag)
        for i, c in enumerate(_STRESS_COMPS[:len(s_rep)]):
            print("    %-4s replay=% .6e  odb=% .6e  abs=%.3e  rel=%.3e"
                  % (c, s_rep[i], s_odb[i], aerr[i], rerr[i]))
        print("    max |abs|=%.3e  max rel=%.3e" % (aerr.max(), rerr.max()))


def print_sdv_table(recs, ref_SDV, labels, sdv_labels):
    tiny = 1e-12
    if all(s is None for s in ref_SDV):
        print("\n(no SDV in fields.json; skipping STATEV vs SDV comparison)")
        return
    print("\n==== replayed STATEV vs ODB SDV ====")
    if sdv_labels is None:
        print("  (no 'sdv_labels' in fields.json; assuming SDV list = STATEV[1..len])")
    for k, rec in enumerate(recs):
        if k >= len(ref_SDV) or ref_SDV[k] is None:
            continue
        sdv_odb = ref_SDV[k]
        statev = np.asarray(rec["statev"], dtype=float)
        if sdv_labels is not None:
            idx = [int(l) - 1 for l in sdv_labels]
        else:
            idx = list(range(len(sdv_odb)))
        print("  frame %s (kinc=%d)" % (labels[k], rec["kinc"]))
        for j, comp in enumerate(idx):
            if comp < 0 or comp >= len(statev):
                continue
            rep = statev[comp]
            odb = sdv_odb[j]
            aerr = abs(rep - odb)
            rerr = aerr / max(abs(odb), tiny)
            lab = ("SDV%d" % sdv_labels[j]) if sdv_labels is not None else ("SDV%d" % (comp + 1))
            print("    %-8s replay=% .6e  odb=% .6e  abs=%.3e  rel=%.3e"
                  % (lab, rep, odb, aerr, rerr))


# --------------------------------------------------------------------------
# Modes
# --------------------------------------------------------------------------
# Default material constants: HCP alpha-uranium example (Compression111.inp).
DEFAULT_PROPS = [0.0, 0.89931, -0.4373, 0.0, 0.26422, 0.54336,
                 -0.79684, 0.34846, 0.71661, 0.6042, 1.0]


def parse_props(arg):
    if not arg:
        return list(DEFAULT_PROPS)
    return [float(x) for x in arg.replace(",", " ").split()]


def dry_run(args):
    nstatv = args.nstatv
    props = parse_props(args.props)
    nprops = len(props)
    ntens = 6
    statev0 = [0.0] * nstatv

    # synthetic single-IP history: increment 1 is identity (known answer:
    # STRESS == 0), then a small uniaxial stretch ramp along z.
    incs = []
    n = args.increments
    F_prev = np.eye(3)
    for k in range(1, n + 1):
        ezz = 0.0 if k == 1 else 1.0e-3 * (k - 1)
        F1 = np.diag([1.0, 1.0, 1.0 + ezz])
        # TIME(1)/TIME(2) are both the increment-START time; DTIME is the length.
        incs.append(dict(noel=1, npt=1, kstep=1, kinc=k, dtime=0.01,
                         temp=293.0, dtemp=0.0,
                         time=(0.01 * (k - 1), 0.01 * (k - 1)),
                         dfgrd0=F_prev, dfgrd1=F1))
        F_prev = F1

    here = os.path.dirname(os.path.abspath(__file__))
    work = os.path.join(here, "build")
    os.makedirs(work, exist_ok=True)
    infile = os.path.join(work, "replay_dryrun_in.txt")
    outfile = os.path.join(work, "replay_dryrun_out.txt")
    write_driver_input(infile, nstatv, nprops, ntens, props, statev0, incs)
    print("wrote driver input : %s" % infile)

    driver = find_driver(args.driver, prefer_mock=args.mock or not args.driver)
    if args.no_run or driver is None:
        if driver is None and not args.no_run:
            print("\nNO DRIVER BINARY FOUND under build/ (run build.sh / build.bat).")
            print("The driver input above is ready; nothing was executed.")
        else:
            print("--no-run: not executing any binary.")
        return 0

    is_mock = "mock" in os.path.basename(driver)
    print("driver             : %s%s" % (driver, "   [MOCK material]" if is_mock else "   [REAL UMAT]"))
    recs = run_driver(driver, infile, outfile)
    print("wrote driver output: %s\n" % outfile)

    print("==== per-increment replay (dry-run) ====")
    for rec in recs:
        s = rec["stress"]
        st = rec["statev"]
        # STATEV(35)=cumulative slip-like scalar; STATEV(48:53)=Cauchy stress
        print("  kinc=%d  STRESS=[%s]  STATEV35=%.3e  STATEV48-50=[% .3e % .3e % .3e]  pnewdt=%.2f"
              % (rec["kinc"], ", ".join("% .4e" % v for v in s),
                 st[34] if len(st) > 34 else float("nan"),
                 st[47] if len(st) > 47 else float("nan"),
                 st[48] if len(st) > 48 else float("nan"),
                 st[49] if len(st) > 49 else float("nan"),
                 rec["pnewdt"]))

    if is_mock:
        s1 = np.asarray(recs[0]["stress"], dtype=float)
        ok = np.allclose(s1, 0.0, atol=1e-9)
        print("\nsanity (identity increment -> zero stress): %s (|S|max=%.2e)"
              % ("PASS" if ok else "FAIL", np.abs(s1).max()))
        print("NOTE: results are from the MOCK elastic UMAT, not crystal plasticity.")
    return 0


def fields_run(args):
    fields = load_fields(args.fields)
    eid, ip = args.element, args.ip
    props = parse_props(args.props)
    nprops = len(props)
    nstatv = args.nstatv
    ntens = 6
    statev0 = [0.0] * nstatv
    sdv_labels = fields.get("sdv_labels")  # optional extension emitted by extractor

    incs, ref_S, ref_SDV, labels = element_history(fields, eid, ip)
    print("element %d  IP %d : %d output frame(s) -> %d increment(s)"
          % (eid, ip, len(labels), len(incs)))
    if not _HAVE_C3D8:
        print("(using inline C3D8 kinematics; residual_core/c3d8_residual not importable)")

    here = os.path.dirname(os.path.abspath(__file__))
    work = os.path.join(here, "build")
    os.makedirs(work, exist_ok=True)
    infile = os.path.join(work, "replay_e%d_ip%d_in.txt" % (eid, ip))
    outfile = os.path.join(work, "replay_e%d_ip%d_out.txt" % (eid, ip))
    write_driver_input(infile, nstatv, nprops, ntens, props, statev0, incs)
    print("wrote driver input : %s" % infile)

    driver = find_driver(args.driver, prefer_mock=args.mock)
    if args.no_run or driver is None:
        if driver is None and not args.no_run:
            print("\nNO DRIVER BINARY FOUND under build/ (run build.sh / build.bat).")
            print("Driver input is ready for inspection; nothing executed.")
        return 0
    is_mock = "mock" in os.path.basename(driver)
    print("driver             : %s%s" % (driver, "   [MOCK]" if is_mock else "   [REAL UMAT]"))
    recs = run_driver(driver, infile, outfile)
    print("wrote driver output: %s" % outfile)

    print_stress_table(recs, ref_S, labels)
    print_sdv_table(recs, ref_SDV, labels, sdv_labels)
    if is_mock:
        print("\nNOTE: comparison is against the MOCK elastic UMAT; build the real")
        print("      umat_driver (ifort+Abaqus) for a crystal-plasticity comparison.")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description="Replay the Grilli UMAT at one IP, "
                                             "increment by increment.")
    ap.add_argument("--fields", help="fields.json (CONTRACT sec 5); omit for --dry-run")
    ap.add_argument("--dry-run", action="store_true",
                    help="fabricate an identity+ramp single-IP history")
    ap.add_argument("--element", type=int, default=1, help="element label (fields mode)")
    ap.add_argument("--ip", type=int, default=1, help="integration point 1..8 (fields mode)")
    ap.add_argument("--props", default=None,
                    help="comma/space list of the material constants "
                         "(default: HCP Compression111 constants)")
    ap.add_argument("--nstatv", type=int, default=125, help="number of SDVs (default 125)")
    ap.add_argument("--increments", type=int, default=5, help="dry-run increment count")
    ap.add_argument("--driver", default=None, help="path to a umat_driver binary")
    ap.add_argument("--mock", action="store_true", help="force the mock driver binary")
    ap.add_argument("--no-run", action="store_true",
                    help="only write the driver input; do not execute a binary")
    args = ap.parse_args(argv)

    if args.dry_run or not args.fields:
        return dry_run(args)
    return fields_run(args)


if __name__ == "__main__":
    sys.exit(main())
