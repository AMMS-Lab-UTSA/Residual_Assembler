#!/usr/bin/env python3
"""PROGRAM 1 (JHU-side) -- build the elastic OTI provider vertical slice.

source/config -> OTI/provider build -> loadable binary + manifest -> primal parity

Given the provider source that implements the agreed C ABI, this:
  1. compiles it to a loadable shared object (.so);
  2. records source/model hash, build id, ABI version, parameter/state layout;
  3. verifies the OTI material's REAL part matches the regular material response
     (the matched-twin parity check JHU owns);
  4. emits material_manifest.json (resasm_material_package_v1) and a build report;
  5. packages ONLY the distributable artifacts into ./dist/.

This is the ONLY place that knows the constitutive law. Nothing here does FE
mesh handling, residual assembly or sensitivity solves -- those belong to the
collaborator tool (Program 2, the Residual_Assembler repo).

    python build_provider.py            # build + validate + package into ./dist
"""

import ctypes
import hashlib
import json
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CONTRACT = os.path.join(HERE, os.pardir, "contract")
_SRC_ROOT = os.path.join(HERE, os.pardir, os.pardir, "src")
if _SRC_ROOT not in sys.path:
    sys.path.insert(0, _SRC_ROOT)
from umat_oti.runtime import load_shared_library, shared_library_suffix   # noqa: E402

SRC = os.path.join(HERE, "elastic_reference.f90")
DIST = os.path.join(HERE, "dist")
MODEL_ID = "reference_elastic_isotropic"
SO_NAME = "libmat_elastic_oti" + shared_library_suffix()


class MatDesc(ctypes.Structure):
    _fields_ = [("abi_version", ctypes.c_int), ("ntens", ctypes.c_int),
                ("nprops", ctypes.c_int), ("nstatev", ctypes.c_int),
                ("kinematics", ctypes.c_int), ("order", ctypes.c_int)]


def _sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()[:16]


def _regular_D(E, nu):
    """The REGULAR elastic material (analytic), Abaqus Voigt, engineering shear."""
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    D = np.zeros((6, 6))
    D[:3, :3] = lam
    for i in range(3):
        D[i, i] += 2 * mu
    for i in range(3, 6):
        D[i, i] = mu
    return D


def compile_so():
    os.makedirs(DIST, exist_ok=True)
    so = os.path.join(DIST, SO_NAME)
    subprocess.check_call([os.environ.get("FC") or "gfortran",
                           "-shared", "-fPIC", "-O2", SRC, "-o", so])
    return so


def parity_check(so):
    """Verify the OTI binary's REAL stress + tangent match the regular material.
    (For a real model the regular and OTI objects are separate builds from the
    same source; here they share the analytic law so parity is exact.)"""
    lib = load_shared_library(so).lib
    lib.mat_eval_v1.restype = ctypes.c_int
    lib.mat_eval_v1.argtypes = [ctypes.POINTER(MatDesc), ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int), ctypes.c_int, ctypes.POINTER(ctypes.c_double),
        ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
        ctypes.c_void_p, ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int)]
    desc = MatDesc(1, 6, 2, 0, 0, 1)
    rng = np.random.default_rng(0)
    max_sig = max_tan = 0.0
    for _ in range(50):
        E, nu = rng.uniform(1e4, 3e5), rng.uniform(0.05, 0.45)
        eps = rng.uniform(-1e-2, 1e-2, 6)
        props = (ctypes.c_double * 2)(E, nu)
        seed = (ctypes.c_int * 2)(1, 2)
        kin = (ctypes.c_double * 6)(*eps)
        tdata = (ctypes.c_double * 4)(0, 1, 0, 0)
        sig = (ctypes.c_double * 6)(); dsig = (ctypes.c_double * 12)(); dds = (ctypes.c_double * 36)()
        st = ctypes.c_int(0); nul = ctypes.c_void_p(0)
        rc = lib.mat_eval_v1(ctypes.byref(desc), props, seed, 2, kin, nul, nul, nul, tdata,
                             sig, dsig, nul, nul, dds, ctypes.byref(st))
        assert rc == 0, "provider returned rc=%d" % rc
        D = _regular_D(E, nu)
        max_sig = max(max_sig, float(np.max(np.abs(np.array(sig) - D @ eps))))
        max_tan = max(max_tan, float(np.max(np.abs(np.array(dds).reshape(6, 6) - D))))
    return {"stress_real_vs_regular_max_abs": max_sig,
            "tangent_vs_regular_max_abs": max_tan, "samples": 50}


def main():
    t0 = time.time()
    contract_ver = json.load(open(os.path.join(CONTRACT, "CONTRACT_VERSION.json")))
    src_hash = _sha(SRC)
    so = compile_so()
    so_hash = _sha(so)
    parity = parity_check(so)
    build_id = "elastic-ref-%s" % src_hash

    manifest = {
        "schema": "resasm_material_package_v1",
        "model_id": MODEL_ID,
        "kinematics": "small_strain",
        "ntens": 6, "nprops": 2, "nstatev": 0,
        "parameters": [
            {"name": "E", "index": 1, "oti_direction": 1, "units": "MPa"},
            {"name": "nu", "index": 2, "oti_direction": 2, "units": "-"},
        ],
        "state_layout": [],
        "outputs": ["stress", "consistent_tangent"],
        "sensitivity_capabilities": {"runtime_parameter_seeding": True,
                                     "maximum_directions": 2, "order": 1},
        "abi": {"symbol": "mat_eval_v1", "version": 1},
        "contract_version": contract_ver["combined_hash"],
        "binaries": {
            # regular twin: same source, its own (analytic) build identity
            "regular": {"hash": src_hash, "build_id": build_id, "note": "same source as OTI"},
            "oti": {"path": SO_NAME, "hash": so_hash, "build_id": build_id},
        },
        "provenance": {
            "generated_by": "oti_provider/elastic/build_provider.py",
            "source_hash": src_hash,
            "note": "REFERENCE elastic provider -- analytic derivatives, NOT a real OTI toolchain build",
        },
    }
    json.dump(manifest, open(os.path.join(DIST, "material_manifest.json"), "w"), indent=2)

    report = {
        "schema": "resasm_provider_build_report_v1",
        "model_id": MODEL_ID, "build_id": build_id,
        "source_hash": src_hash, "oti_binary_hash": so_hash,
        "contract_version": contract_ver["combined_hash"],
        "abi_version": 1,
        "parity": parity,
        "parity_passed": bool(parity["stress_real_vs_regular_max_abs"] < 1e-9
                              and parity["tangent_vs_regular_max_abs"] < 1e-9),
        "build_seconds": time.time() - t0,
    }
    json.dump(report, open(os.path.join(DIST, "build_report.json"), "w"), indent=2)

    print("=" * 68)
    print(" Program 1 provider build  |  model=%s" % MODEL_ID)
    print("=" * 68)
    print(" distributable artifacts -> %s" % DIST)
    print("   %s        (opaque OTI binary)" % SO_NAME)
    print("   material_manifest.json          (contract %s)" % contract_ver["combined_hash"])
    print("   build_report.json")
    print(" source_hash=%s  oti_hash=%s  build_id=%s" % (src_hash, so_hash, build_id))
    print(" real-part parity vs regular: stress %.2e  tangent %.2e  -> %s"
          % (parity["stress_real_vs_regular_max_abs"], parity["tangent_vs_regular_max_abs"],
             "PASS" if report["parity_passed"] else "FAIL"))
    print(" NOTE: JHU distributes ONLY ./dist to collaborators. Source never leaves.")
    return 0 if report["parity_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
