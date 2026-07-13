"""
run_checks.py -- Step-4 / Step-7 verification driver for c3d8_residual.

Runs the three Abaqus-INDEPENDENT correctness checks from CONTRACT section 4,
prints a results table, writes the numeric report to `last_report.txt`, and
exits nonzero on any failure.

    1. Divergence-theorem patch test  (real Compression111 element + distorted)
    2. Uniform uniaxial sanity         (unit cube)
    3. Tangent finite-difference check:
         3a  linear small-strain analytic tangent K = sum B0^T D B0 detJ0 w
         3b  finite-strain force FD vs exact fixed-sigma Jacobian

No Abaqus is required or invoked.
"""

import os
import sys
import datetime
import numpy as np

# make residual_core importable regardless of CWD / test-dir depth
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = _HERE
while _ROOT != os.path.dirname(_ROOT) and not os.path.isdir(os.path.join(_ROOT, "residual_core")):
    _ROOT = os.path.dirname(_ROOT)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from residual_core.formulations.c3d8_kernel import (   # noqa: E402
    COMPRESSION111_ELEM1_XE,
    element_internal_force_small_strain,
    element_internal_force_finite_strain,
    element_tangent,
    force_tangent_fixed_sigma,
    surface_traction_nodal_forces,
    linear_stress_target,
    shape_functions,
    b_matrix_reference,
    isotropic_D,
    unit_cube_Xe,
    ABAQUS_C3D8_GAUSS,
)

REPORT_PATH = os.path.join(_HERE, "last_report.txt")


def check1_divergence():
    """Divergence-theorem patch test. Returns list of row dicts."""
    rng = np.random.default_rng(12345)
    rows = []
    geoms = {
        "Compression111_elem1": COMPRESSION111_ELEM1_XE.copy(),
        "distorted_hex": COMPRESSION111_ELEM1_XE + rng.uniform(-0.6, 0.6, (8, 3)),
    }
    stresses = {
        "random_symmetric": rng.uniform(-100, 100, size=6),
        "uniaxial_s33_250": np.array([0, 0, 250.0, 0, 0, 0]),
    }
    for gname, Xe in geoms.items():
        for sname, s0 in stresses.items():
            f_vol = element_internal_force_small_strain(Xe, np.tile(s0, (8, 1)))
            f_surf = surface_traction_nodal_forces(Xe, s0)
            denom = max(np.linalg.norm(f_surf), 1e-30)
            rel = float(np.linalg.norm(f_vol - f_surf) / denom)
            equil = float(np.linalg.norm(f_vol.reshape(8, 3).sum(axis=0)))
            ok = (rel < 1e-8) and (equil < 1e-8 * max(denom, 1.0))
            rows.append(dict(name="1 divergence: %s / %s" % (gname, sname),
                             metric="rel=%.3e |sumF|=%.3e" % (rel, equil),
                             tol="rel<1e-8", ok=ok))
    return rows


def check1b_linear_stress():
    """LINEAR (spatially varying) stress patch test: catches per-IP weight and
    sigma<->point pairing errors that the uniform patch test is blind to.
    Affine cube -> exact quadrature -> machine precision. Does NOT probe the
    Abaqus IP *ordering* (needs an ODB)."""
    rng = np.random.default_rng(2024)
    S0 = rng.uniform(-50, 50, (3, 3)); S0 = 0.5 * (S0 + S0.T)
    grads = [0.5 * (G + G.T) for G in rng.uniform(-8, 8, (3, 3, 3))]

    def sigma_fn(x):
        M = S0 + grads[0] * x[0] + grads[1] * x[1] + grads[2] * x[2]
        return np.array([M[0, 0], M[1, 1], M[2, 2], M[0, 1], M[0, 2], M[1, 2]])

    Xe = COMPRESSION111_ELEM1_XE.copy()
    pts = ABAQUS_C3D8_GAUSS.points
    sigma_ip = np.array([sigma_fn(shape_functions(pts[k]) @ Xe) for k in range(8)])
    f_asm = element_internal_force_small_strain(Xe, sigma_ip)
    f_tgt = linear_stress_target(Xe, sigma_fn)
    rel = float(np.linalg.norm(f_asm - f_tgt) / max(np.linalg.norm(f_tgt), 1e-30))
    ok = rel < 1e-9
    return [dict(name="1b linear-stress patch (per-IP weights+pairing)",
                 metric="rel=%.3e" % rel, tol="rel<1e-9", ok=ok)]


def check2_uniaxial():
    Xe = unit_cube_Xe()
    S33 = 300.0
    s0 = np.array([0, 0, S33, 0, 0, 0.0])
    f = element_internal_force_small_strain(Xe, np.tile(s0, (8, 1))).reshape(8, 3)
    top, bot = [4, 5, 6, 7], [0, 1, 2, 3]
    expect = S33 * 1.0 / 4.0
    err = max(
        float(np.max(np.abs(f[top, 2] - expect))),
        float(np.max(np.abs(f[bot, 2] + expect))),
        float(np.max(np.abs(f[:, 0]))),
        float(np.max(np.abs(f[:, 1]))),
    )
    ok = err < 1e-10
    return [dict(name="2 uniaxial unit cube (+z=+S33A/4, sides 0)",
                 metric="max_err=%.3e" % err, tol="<1e-10", ok=ok)]


def check3a_linear_tangent():
    Xe = COMPRESSION111_ELEM1_XE.copy()
    D = isotropic_D(E=200000.0, nu=0.3)
    pts, wts = ABAQUS_C3D8_GAUSS.points, ABAQUS_C3D8_GAUSS.weights
    B0s = [b_matrix_reference(Xe, pts[k])[0] for k in range(8)]

    def residual(U):
        sig = np.array([D @ (B0s[k] @ U) for k in range(8)])
        return element_internal_force_small_strain(Xe, sig)

    K = element_tangent(Xe, np.zeros(24), D, mode='small')
    scale = np.linalg.norm(Xe) / np.sqrt(8)
    h = 1e-6 * max(scale, 1.0)
    Kfd = np.zeros((24, 24))
    for j in range(24):
        Up = np.zeros(24); Up[j] = h
        Um = np.zeros(24); Um[j] = -h
        Kfd[:, j] = (residual(Up) - residual(Um)) / (2 * h)
    abs_err = float(np.max(np.abs(Kfd - K)))
    fro = float(np.linalg.norm(Kfd - K) / max(np.linalg.norm(K), 1e-30))
    max_entry = abs_err
    ok = fro < 1e-6
    row = dict(name="3a linear tangent FD vs analytic K",
               metric="abs=%.3e rel_Fro=%.3e max=%.3e" % (abs_err, fro, max_entry),
               tol="rel<1e-6", ok=ok)
    return [row], dict(abs_err=abs_err, rel_fro=fro, max_entry=max_entry)



def check3b_finite_strain(U_amp):
    rng = np.random.default_rng(7)
    Xe = COMPRESSION111_ELEM1_XE.copy()
    sigma_ip = np.tile(rng.uniform(-50, 50, size=6), (8, 1))
    U0 = rng.uniform(-U_amp, U_amp, size=24) if U_amp > 0 else np.zeros(24)

    def force(U):
        return element_internal_force_finite_strain(Xe, U, sigma_ip)

    scale = np.linalg.norm(Xe) / np.sqrt(8)
    h = 1e-6 * max(scale, 1.0)
    Kfd = np.zeros((24, 24))
    for j in range(24):
        Up = U0.copy(); Up[j] += h
        Um = U0.copy(); Um[j] -= h
        Kfd[:, j] = (force(Up) - force(Um)) / (2 * h)
    Kan = force_tangent_fixed_sigma(Xe, U0, sigma_ip)
    abs_err = float(np.max(np.abs(Kfd - Kan)))
    fro = float(np.linalg.norm(Kfd - Kan) / max(np.linalg.norm(Kan), 1e-30))
    ok = fro < 1e-6
    row = dict(name="3b finite-strain force FD vs fixed-sigma J (U_amp=%.2f)" % U_amp,
               metric="abs=%.3e rel_Fro=%.3e" % (abs_err, fro),
               tol="rel<1e-6", ok=ok)
    return [row], dict(U_amp=U_amp, abs_err=abs_err, rel_fro=fro)


def main():
    rows = []
    rows += check1_divergence()
    rows += check1b_linear_stress()
    rows += check2_uniaxial()
    r3a, d3a = check3a_linear_tangent()
    rows += r3a
    r3b0, d3b0 = check3b_finite_strain(0.0)
    r3b1, d3b1 = check3b_finite_strain(0.5)
    rows += r3b0 + r3b1

    # ---- print table ----
    namew = max(len(r["name"]) for r in rows)
    metw = max(len(r["metric"]) for r in rows)
    line = "-" * (namew + metw + 28)
    print(line)
    print("%-*s  %-*s  %-10s  %s" % (namew, "CHECK", metw, "METRIC", "TOL", "RESULT"))
    print(line)
    for r in rows:
        print("%-*s  %-*s  %-10s  %s"
              % (namew, r["name"], metw, r["metric"], r["tol"],
                 "PASS" if r["ok"] else "FAIL"))
    print(line)
    all_ok = all(r["ok"] for r in rows)
    print("OVERALL: %s" % ("ALL PASS" if all_ok else "FAILURE"))

    # ---- write numeric report ----
    with open(REPORT_PATH, "w") as fh:
        fh.write("c3d8_residual verification report\n")
        fh.write("generated: %s\n" % datetime.datetime.now().isoformat(timespec="seconds"))
        fh.write("numpy: %s\n\n" % np.__version__)
        fh.write("%-58s %-32s %-10s %s\n" % ("CHECK", "METRIC", "TOL", "RESULT"))
        fh.write("%s\n" % ("-" * 110))
        for r in rows:
            fh.write("%-58s %-32s %-10s %s\n"
                     % (r["name"], r["metric"], r["tol"], "PASS" if r["ok"] else "FAIL"))
        fh.write("\n")
        fh.write("TANGENT FINITE-DIFFERENCE NUMERIC DETAIL\n")
        fh.write("  3a linear small-strain tangent:\n")
        fh.write("      absolute error (max entry) : %.6e\n" % d3a["abs_err"])
        fh.write("      relative error (Frobenius) : %.6e   (assert < 1e-6)\n" % d3a["rel_fro"])
        fh.write("      max single-entry error     : %.6e\n" % d3a["max_entry"])
        fh.write("  3b finite-strain force FD vs exact fixed-sigma Jacobian:\n")
        for d in (d3b0, d3b1):
            fh.write("      U_amp=%.2f  abs=%.6e  rel_Fro=%.6e   (assert < 1e-6)\n"
                     % (d["U_amp"], d["abs_err"], d["rel_fro"]))
        fh.write("\nOVERALL: %s\n" % ("ALL PASS" if all_ok else "FAILURE"))
    print("wrote %s" % REPORT_PATH)

    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
