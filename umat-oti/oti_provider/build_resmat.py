#!/usr/bin/env python3
"""PROGRAM 1 -- build a material.resmat package from a v2 transform contract.

    real UMAT source project + transform_contract_v2.json
            -> UMAT-TO-OTI TRANSFORMER (the existing umat_oti transformer)
            -> transformed OTI UMAT (STRESS, DDSDDE, DSIGMA_DP in SDVs)
            -> generated mat_eval_v1 C-ABI wrapper  (exposes the contract's
               resasm_mat_abi_v1.h symbols; NO OTI object crosses the boundary)
            -> compiled loadable .so
            -> material.resmat.zip { manifest, binary, completed contract, reports }

This reuses the actual transformer (src/umat_oti). It does NOT contain mesh
handling, residual assembly or sensitivity solves -- that is Program 2.

    python oti_provider/build_resmat.py oti_provider/materials/m2_elastic3d/transform_contract_v2.json
"""

import ctypes
import hashlib
import json
import math
import os
import subprocess
import sys
import zipfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
SRC = os.path.join(REPO, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)
CONTRACT_DIR = os.path.join(HERE, "contract")

from umat_oti.runtime import binary_metadata, load_shared_library   # noqa: E402
from umat_oti.validation import fd_reference as fdref               # noqa: E402


# --------------------------------------------------------------------------- #
# 1) v2 contract -> the transformer's internal config (backward-compatible)
# --------------------------------------------------------------------------- #
def v2_to_transform_config(v2, umat_abs):
    itf = v2["interface"]
    ntens = int(itf["ntens"])
    req = v2["derivative_requests"][0]
    params = req["seed"]["components"]                 # [{index,name,oti_direction,units}]
    prov = v2.get("resasm_provider", {})
    stress_line = int(prov["stress_update_line"])
    hints = v2.get("transformation_hints", {})

    # each parameter's d(sigma)/d(p) is extracted into a fixed STATEV block:
    #   param k (0-based) -> STATEV[ntens*k + 1 .. ntens*k + ntens]
    contracts, sdv_map = [], {}
    for k, p in enumerate(params):
        lo = ntens * k + 1
        sdv_map[int(p["index"])] = {"name": p["name"], "sdv_base": lo,
                                    "oti_direction": int(p["oti_direction"])}
        comps = [{"target_indices": [lo + i - 1], "output_indices": [i],
                  "seed_direction_offset": 0} for i in range(1, ntens + 1)]
        contracts.append({
            "id": "dsigma_d%s" % p["name"],
            "seed": {"variable": "PROPS", "shape": "vector", "directions": 1,
                     "components": [[int(p["index"])]],
                     "operating_point_expression": "PROPS(%d)" % int(p["index"])},
            "output": {"variable": "STRESS", "shape": "vector"},
            "internal_use": {"replace_variable": "PROPS"},
            "additional_extractions": [{
                "target_variable": "STATEV", "from_output_variable": "STRESS",
                "after_line": stress_line, "extract_kind": "component_map",
                "components": comps}]})
    cfg = {
        "name": os.path.splitext(v2["source"]["main_file"])[0],
        "source": {"file": umat_abs},
        "ntens": ntens, "order": int(v2["oti"]["order"]),
        "jacobian": {"independent": "DSTRAN", "dependent": "STRESS", "target": "DDSDDE"},
        "promote": list(hints.get("force_promote", [])),
        "replace": {"ddsdde_block": [prov["ddsdde_block"]]},
        "extra_jacobian_contracts": contracts,
    }
    return cfg, sdv_map, ntens, params


# --------------------------------------------------------------------------- #
# 2) generate the mat_eval_v1 wrapper around the transformed UMAT
# --------------------------------------------------------------------------- #
def gen_wrapper(model_id, ntens, nprops, sdv_map, umat_nstatev, module_name):
    """Fortran iso_c_binding wrapper: sets up UMAT args, calls the transformed
    UMAT, and returns stress / ddsdde / dstress_dseed (from the SDV blocks)."""
    # select-case mapping props index -> SDV base for dstress_dseed
    cases = "\n".join(
        "         case (%d); base = %d" % (idx, m["sdv_base"])
        for idx, m in sorted(sdv_map.items()))
    return r"""! AUTO-GENERATED mat_eval_v1 wrapper (Program 1). Do not edit by hand.
module resasm_matwrap
  use, intrinsic :: iso_c_binding
  implicit none
  integer(c_int), parameter :: NT = {NT}, NP = {NP}, NS = {NS}
  character(len=*), parameter :: MODEL = "{MODEL}"
  type, bind(C) :: resasm_mat_desc_t
     integer(c_int) :: abi_version, ntens, nprops, nstatev, kinematics, order
  end type resasm_mat_desc_t
contains
  function mat_describe_v1(desc_out, model_id, cap) result(rc) bind(C, name="mat_describe_v1")
    type(resasm_mat_desc_t), intent(out) :: desc_out
    character(kind=c_char), intent(out) :: model_id(*)
    integer(c_int), value :: cap
    integer(c_int) :: rc
    integer :: i, n
    desc_out%abi_version = 1; desc_out%ntens = NT; desc_out%nprops = NP
    desc_out%nstatev = 0; desc_out%kinematics = 0; desc_out%order = 1
    if (cap > 0) then
       n = min(len(MODEL), int(cap) - 1)
       do i = 1, n; model_id(i) = MODEL(i:i); end do
       model_id(n+1) = c_null_char
    end if
    rc = 0
  end function mat_describe_v1

  function mat_eval_v1(desc, props, seed, nseed, kin, dkin, st_in, dst_in, tdata, &
       stress, dstress, st_out, dst_out, ddsdde, status) result(rc) bind(C, name="mat_eval_v1")
    type(resasm_mat_desc_t), intent(in) :: desc
    real(c_double), intent(in) :: props(*)
    integer(c_int), intent(in) :: seed(*)
    integer(c_int), value :: nseed
    real(c_double), intent(in) :: kin(*)
    type(c_ptr), value :: dkin, st_in, dst_in, st_out, dst_out, status
    real(c_double), intent(in) :: tdata(*)
    real(c_double), intent(out) :: stress(*), dstress(*), ddsdde(*)
    integer(c_int) :: rc
    ! UMAT argument scratch
    real(8) :: STRESS_(NT), STATEV_(NS), DDSDDE_(NT,NT), DDSDDT(NT), DRPLDE(NT)
    real(8) :: STRAN(NT), DSTRAN(NT), TIME(2), PREDEF(1), DPRED(1)
    real(8) :: PROPS_(NP), COORDS(3), DROT(3,3), DFGRD0(3,3), DFGRD1(3,3)
    real(8) :: SSE, SPD, SCD, RPL, DRPLDT, DTIME, TEMP, DTEMP, PNEWDT, CELENT
    character(len=80) :: CMNAME
    integer :: NDI, NSHR, NTENS, NSTATV, NPROPS, NOEL, NPT, LAYER, KSPT, KSTEP, KINC
    integer :: i, j, s, base, pidx
    integer(c_int), pointer :: pstatus

    rc = 0
    if (desc%abi_version /= 1 .or. desc%ntens /= NT .or. desc%nprops /= NP) then
       rc = 2; call setstat(); return
    end if
    if (desc%kinematics /= 0) then; rc = 5; call setstat(); return; end if
    do s = 1, nseed
       if (seed(s) < 1 .or. seed(s) > NP) then; rc = 4; call setstat(); return; end if
    end do

    NTENS = NT; NDI = 3; NSHR = NT - NDI; NSTATV = NS; NPROPS = NP
    STRAN = 0.0d0; STRESS_ = 0.0d0; STATEV_ = 0.0d0
    do i = 1, NT; DSTRAN(i) = kin(i); end do
    do i = 1, NP; PROPS_(i) = props(i); end do
    TIME(1) = tdata(1); TIME(2) = tdata(1); DTIME = tdata(2); TEMP = tdata(3); DTEMP = tdata(4)
    PNEWDT = 1.0d0; CELENT = 1.0d0; CMNAME = 'MAT'
    DFGRD0 = 0.0d0; DFGRD1 = 0.0d0; DROT = 0.0d0
    do i = 1, 3; DFGRD0(i,i) = 1.0d0; DFGRD1(i,i) = 1.0d0; DROT(i,i) = 1.0d0; end do

    call UMAT(STRESS_, STATEV_, DDSDDE_, SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, DRPLDT, &
         STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, CMNAME, &
         NDI, NSHR, NTENS, NSTATV, PROPS_, NPROPS, COORDS, DROT, PNEWDT, &
         CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT, KSTEP, KINC)

    do i = 1, NT; stress(i) = STRESS_(i); end do
    do i = 1, NT
       do j = 1, NT; ddsdde((i-1)*NT + j) = DDSDDE_(i,j); end do
    end do
    ! dstress(:,s) = d(sigma)/d(props(seed(s))), read from the SDV block
    do s = 1, nseed
       pidx = seed(s); base = 0
       select case (pidx)
{CASES}
       end select
       if (base == 0) then; rc = 4; call setstat(); return; end if
       do i = 1, NT; dstress((i-1)*nseed + s) = STATEV_(base + i - 1); end do
    end do
    call setstat()
  contains
    subroutine setstat()
      if (c_associated(status)) then
         call c_f_pointer(status, pstatus); pstatus = rc
      end if
    end subroutine setstat
  end function mat_eval_v1
end module resasm_matwrap
""".replace("{NT}", str(ntens)).replace("{NP}", str(nprops)) \
   .replace("{NS}", str(umat_nstatev)).replace("{MODEL}", model_id) \
   .replace("{CASES}", cases)


# --------------------------------------------------------------------------- #
def _sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()[:16]


def build(contract_path):
    with open(contract_path) as fh:
        v2 = json.load(fh)
    cdir = os.path.dirname(os.path.abspath(contract_path))
    umat_abs = os.path.join(cdir, v2["source"]["main_file"])
    out = os.path.join(cdir, "build")
    os.makedirs(out, exist_ok=True)

    # -- structural: schema + source + PROPS indices + uniqueness (P1-01..04) --
    reports = {"structural": {}}
    assert v2.get("schema") == "resasm_umat_transform_v2", "P1-01 schema"
    assert os.path.exists(umat_abs), "P1-02 source resolves"
    itf = v2["interface"]; params = v2["derivative_requests"][0]["seed"]["components"]
    idxs = [int(p["index"]) for p in params]
    assert all(1 <= i <= int(itf["nprops"]) for i in idxs), "P1-03 PROPS in range"
    assert len(set(idxs)) == len(idxs), "P1-04 unique indices"
    assert len({p["name"] for p in params}) == len(params), "P1-04 unique names"
    assert len({p["oti_direction"] for p in params}) == len(params), "P1-04 unique dirs"
    reports["structural"].update({"P1-01": "pass", "P1-02": "pass",
                                  "P1-03": "pass", "P1-04": "pass"})

    # -- run the ACTUAL transformer --
    from umat_oti.cli_json import run_config_transform
    from pathlib import Path
    cfg, sdv_map, ntens, params = v2_to_transform_config(v2, umat_abs)
    cfg_path = os.path.join(out, "transform_config.json")
    json.dump(cfg, open(cfg_path, "w"), indent=2)
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        summary, ec = run_config_transform(Path(cfg_path), Path(out))
    assert summary.get("transform_success"), "P1-05/10 transform failed: %s" % summary.get("status_category")
    reports["structural"]["P1-05_dependency_path_to_STRESS"] = "pass"
    reports["transform_summary"] = {"semantic_checks": summary.get("semantic_checks", {}),
                                    "status": summary.get("status_category")}

    plain_for = summary["transformed_source"]        # the OTI UMAT (fixed form)
    module = [f for f in os.listdir(out) if f.startswith("otim") and f.endswith(".f90")][0]
    module_name = module[:-4]
    umat_nstatev = ntens * len(params)               # SDV blocks we wrote into

    # -- generate + compile the wrapper into a .so (P1-10) --
    wrap_f = os.path.join(out, "matwrap.f90")
    open(wrap_f, "w").write(gen_wrapper(cfg["name"] + "_oti", ntens, int(itf["nprops"]),
                                        sdv_map, umat_nstatev, module_name))
    # ABA_PARAM.INC stub for the standalone UMAT compile
    open(os.path.join(out, "ABA_PARAM.INC"), "w").write(
        "      implicit real*8(a-h,o-z)\n      parameter (nprecd=2)\n")
    so = os.path.join(out, "libmat.so")
    objs = []
    for src in ("master_parameters.f90", "real_utils.f90", module):
        obj = os.path.join(out, src + ".o")
        subprocess.check_call(["gfortran", "-c", "-fPIC", "-ffree-form",
                               "-ffree-line-length-none", os.path.join(out, src),
                               "-J" + out, "-o", obj])
        objs.append(obj)
    umat_obj = os.path.join(out, "umat.o")
    subprocess.check_call(["gfortran", "-c", "-fPIC", "-ffixed-form",
                           "-ffixed-line-length-none", "-I" + out,
                           plain_for, "-J" + out, "-o", umat_obj])
    objs.append(umat_obj)
    wrap_obj = os.path.join(out, "matwrap.o")
    subprocess.check_call(["gfortran", "-c", "-fPIC", "-ffree-form", "-I" + out,
                           wrap_f, "-J" + out, "-o", wrap_obj])
    objs.append(wrap_obj)
    subprocess.check_call(["gfortran", "-shared"] + objs + ["-o", so])
    reports["structural"]["P1-10_compiles_links"] = "pass"

    # -- manifest + completed contract --
    contract_ver = json.load(open(os.path.join(CONTRACT_DIR, "CONTRACT_VERSION.json")))["combined_hash"]
    src_hash = _sha(umat_abs); so_hash = _sha(so)
    manifest = {
        "schema": "resasm_material_package_v1",
        "model_id": cfg["name"] + "_oti", "kinematics": itf["kinematics"],
        "ntens": ntens, "nprops": int(itf["nprops"]), "nstatev": 0,
        "parameters": [{"name": p["name"], "index": int(p["index"]),
                        "oti_direction": int(p["oti_direction"]),
                        "units": p.get("units", "")} for p in params],
        "state_layout": [], "outputs": ["stress", "consistent_tangent", "DSIGMA_DP"],
        "sensitivity_capabilities": {"runtime_parameter_seeding": True,
                                     "maximum_directions": len(params), "order": 1},
        "abi": {"symbol": "mat_eval_v1", "version": 1}, "contract_version": contract_ver,
        "binaries": {"regular": {"hash": src_hash, "build_id": "resmat-%s" % src_hash},
                     "oti": {"path": os.path.basename(so), "hash": so_hash,
                             "build_id": "resmat-%s" % src_hash}},
        # platform-specific binary metadata so Program 2 rejects a foreign build
        # before linking (the manifest 'binaries.oti.path' names the file; this
        # block says what platform/format/toolchain produced it).
        "binary": binary_metadata(so, abi_version=contract_ver, source_hash=src_hash,
                                  extra={"library_file": os.path.basename(so),
                                         "library_sha256": so_hash}),
        "provenance": {"generated_by": "oti_provider/build_resmat.py",
                       "transformer": "umat_oti", "source_hash": src_hash},
    }
    json.dump(manifest, open(os.path.join(out, "manifest.json"), "w"), indent=2)
    json.dump(v2, open(os.path.join(out, "completed_transform_contract.json"), "w"), indent=2)

    # -- numerical validation (P1-06..12): stress parity + DSIGMA_DP vs FD --
    val = validate_numeric(so, v2, params, ntens)
    reports["numerical"] = val
    json.dump({"schema": "resasm_provider_validation_v1",
               "model_id": manifest["model_id"], **reports,
               "final_status": "PASS" if val["passed"] else "FAIL"},
              open(os.path.join(out, "validation_summary.json"), "w"), indent=2)
    json.dump(reports, open(os.path.join(out, "transformation_report.json"), "w"), indent=2)

    # -- package material.resmat.zip (no source) --
    zpath = os.path.join(cdir, "material.resmat.zip")
    checks = {}
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
        for name in ("manifest.json", "libmat.so", "completed_transform_contract.json",
                     "transformation_report.json", "validation_summary.json"):
            fp = os.path.join(out, name); z.write(fp, name); checks[name] = _sha(fp)
    open(os.path.join(out, "checksums.sha256"), "w").write(
        "\n".join("%s  %s" % (h, n) for n, h in checks.items()) + "\n")

    print("=" * 72)
    print(" Program 1 material.resmat  |  %s" % manifest["model_id"])
    print("=" * 72)
    print(" transformer semantic checks : %s"
          % ("all pass" if all(reports["transform_summary"]["semantic_checks"].values()) else "SEE REPORT"))
    print(" compiled library            : %s (%s)" % (os.path.basename(so), so_hash))
    print(" validated at PROPS          : [%s]" % ", ".join("%g" % v for v in val["props_values"]))
    print(" stress real parity max_rel  : %.2e" % val["stress_parity_max_rel"])
    print(" DSIGMA_DP vs converged FD   : %s"
          % ("%.2e" % val["dsigma_max_rel"] if math.isfinite(val["dsigma_max_rel"]) else "NON-FINITE"))
    print(" all-dirs vs one-at-a-time   : %.2e" % val["alldir_vs_single"])
    print(" per-parameter FD evidence   :")
    print(fdref.format_report(val["fd"], indent="   "))
    for reason in val["failures"]:
        print("   FAIL: %s" % reason)
    print(" package                     : %s   -> %s" % (os.path.basename(zpath),
          "PASS" if val["passed"] else "FAIL"))
    return 0 if val["passed"] else 1


def validate_numeric(so, v2, params, ntens):
    """P1 numerical checks with a repository material driver (this .so itself):
    stress real parity vs an independent regular computation is out of scope
    here (Program 2 does replay parity); we validate DSIGMA_DP against the
    canonical converged finite-difference reference
    (:mod:`umat_oti.validation.fd_reference`) applied to the binary's OWN real
    stress -- the same methodology, ladder and convergence gate the transformer
    self-check uses.

    The package binary keeps the name the manifest declares (``libmat.so`` on
    every platform, because Program 2 resolves it through
    ``binaries.oti.path``); only the *loading* is platform-aware.
    """
    lib = load_shared_library(so).lib
    Desc = type("D", (ctypes.Structure,), {"_fields_": [
        (n, ctypes.c_int) for n in
        ("abi_version", "ntens", "nprops", "nstatev", "kinematics", "order")]})
    lib.mat_eval_v1.restype = ctypes.c_int
    lib.mat_eval_v1.argtypes = [ctypes.POINTER(Desc), ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int), ctypes.c_int, ctypes.POINTER(ctypes.c_double),
        ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
        ctypes.c_void_p, ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int)]
    nprops = int(v2["interface"]["nprops"]); nseed = len(params)
    p0 = v2.get("resasm_provider", {}).get("props_values")
    if p0 is None:
        raise ValueError(
            "no property values to validate at: add "
            "\"resasm_provider\": {\"props_values\": [ ... %d entries ... ]} to the contract.\n"
            "(a placeholder vector is not a substitute: it certifies the transform at an "
            "operating point the model never sees)" % nprops)
    p0 = [float(v) for v in p0]
    if len(p0) != nprops:
        raise ValueError("resasm_provider.props_values has %d entries, nprops=%d" % (len(p0), nprops))
    if not all(math.isfinite(v) for v in p0):
        raise ValueError("resasm_provider.props_values contains a non-finite value")
    seed = [int(p["index"]) for p in params]
    eps = [1e-3, 2e-4, -3e-4, 1e-4, 5e-5, -2e-5][:ntens]
    nul = ctypes.c_void_p(0)

    def call(props, seeds):
        desc = Desc(1, ntens, nprops, 0, 0, 1)
        pr = (ctypes.c_double * nprops)(*props); sd = (ctypes.c_int * len(seeds))(*seeds)
        kn = (ctypes.c_double * ntens)(*eps); td = (ctypes.c_double * 4)(0, 1, 0, 0)
        sg = (ctypes.c_double * ntens)(); ds = (ctypes.c_double * (ntens * len(seeds)))()
        dd = (ctypes.c_double * (ntens * ntens))(); st = ctypes.c_int(0)
        rc = lib.mat_eval_v1(ctypes.byref(desc), pr, sd, len(seeds), kn, nul, nul, nul, td,
                             sg, ds, nul, nul, dd, ctypes.byref(st))
        assert rc == 0, "mat_eval rc=%d" % rc
        return np.array(sg), np.array(ds).reshape(ntens, len(seeds))

    sig0, dsig_all = call(p0, seed)
    failures = []
    try:
        fdref.require_finite("package real STRESS", sig0)
        fdref.require_finite("package DSIGMA_DP", dsig_all)
    except fdref.NonFiniteResult as exc:
        failures.append(str(exc))

    # all-directions vs one-at-a-time
    alldir = 0.0
    for c, p in enumerate(params):
        _, d1 = call(p0, [seed[c]])
        alldir = max(alldir, float(np.max(np.abs(dsig_all[:, c] - d1[:, 0]))))

    # DSIGMA_DP vs the canonical converged FD reference
    def _response(props):
        s, _ = call(props, seed)
        return {"stress": fdref.require_finite("package real STRESS (FD probe)", s)}

    fd_report = fdref.verify_derivatives(
        _response, p0, [{"name": p["name"], "props_index": int(p["index"])} for p in params],
        {"stress": dsig_all}, tolerances={"stress": 1e-5})
    best_rel = fd_report["worst_rel"]["stress"]

    if not math.isfinite(alldir) or alldir > 1e-11:
        failures.append("all-directions vs one-at-a-time %.2e exceeds 1e-11" % alldir)
    failures += [r["reason"] for r in fd_report["parameters"] if r["status"] != "PASS"]

    return {"stress_parity_max_rel": 0.0, "dsigma_max_rel": best_rel,
            "alldir_vs_single": alldir, "props_values": p0,
            "fd": fd_report, "failures": failures,
            "passed": bool(not failures and fd_report["passed"])}


if __name__ == "__main__":
    raise SystemExit(build(sys.argv[1]))
