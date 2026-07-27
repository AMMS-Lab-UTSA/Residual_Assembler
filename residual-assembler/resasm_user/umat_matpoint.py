"""Material-point numerical validation for umat_oti-transformed UMATs (M5-B2).

This is the *real-response* check the honesty note in :mod:`resasm_user.umat_backend`
reserves: it compiles BOTH the transformed (OTI) UMAT and the ORIGINAL ordinary
UMAT with one shared standalone Fortran driver, runs a single material-point
increment, and compares each SDV-extracted ``dsigma/da_i`` produced by the
transformed UMAT against a central finite difference of the ORIGINAL UMAT's real
stress response. Real compiler (gfortran), real numbers -- not a structural check.

No Abaqus is involved: the driver calls the UMAT subroutine directly with a stub
``ABA_PARAM.INC`` (``implicit real*8``). It needs the OTI library sources the
transformer emits next to the transformed UMAT (``master_parameters.f90``,
``real_utils.f90``, ``otim{N}n1.f90``); when they or gfortran are absent the
check returns ``available=False`` rather than failing.
"""

from __future__ import annotations

import glob
import os
import shutil
import subprocess
import tempfile
from typing import Any, Dict, List, Optional

_ABA_PARAM = "      implicit real*8(a-h,o-z)\n      parameter (nprecd=2)\n"
_DEFAULT_DSTRAN = [1.0e-3, 2.0e-4, -3.0e-4, 1.0e-4, 5.0e-5, -2.0e-5]


class MatPointError(Exception):
    pass


def gfortran_available() -> bool:
    return shutil.which("gfortran") is not None


# --------------------------------------------------------------------------- #
# the shared standalone driver (fixed-form; UMAT-agnostic)
# --------------------------------------------------------------------------- #
def _driver_source() -> str:
    return r"""      PROGRAM MPDRV
C     Standalone material-point driver: read PROPS + DSTRAN from mp_input.txt,
C     call UMAT once, write STRESS and STATEV to mp_output.txt. Used for BOTH the
C     transformed (OTI) and the original UMAT so the comparison is apples-to-apples.
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*80 CMNAME
      PARAMETER (MNT=6, MNS=128, MNP=16)
      DIMENSION STRESS(MNT),STATEV(MNS),DDSDDE(MNT,MNT),DDSDDT(MNT),
     1 DRPLDE(MNT),STRAN(MNT),DSTRAN(MNT),TIME(2),PREDEF(1),DPRED(1),
     2 PROPS(MNP),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
      OPEN(10, FILE='mp_input.txt', STATUS='OLD')
      READ(10,*) NPROPS
      READ(10,*) (PROPS(I), I=1,NPROPS)
      READ(10,*) NTENS
      READ(10,*) (DSTRAN(I), I=1,NTENS)
      CLOSE(10)
      NDI    = 3
      NSHR   = NTENS - NDI
      NSTATV = MNS
      DO I = 1, MNT
         STRESS(I) = 0.0D0
         STRAN(I)  = 0.0D0
      END DO
      DO I = 1, MNS
         STATEV(I) = 0.0D0
      END DO
      DTIME  = 1.0D0
      TEMP   = 0.0D0
      DTEMP  = 0.0D0
      PNEWDT = 1.0D0
      CELENT = 1.0D0
      CMNAME = 'MATERIAL'
      DO I = 1, 3
         DO J = 1, 3
            DFGRD0(I,J) = 0.0D0
            DFGRD1(I,J) = 0.0D0
            DROT(I,J)   = 0.0D0
         END DO
         DFGRD0(I,I) = 1.0D0
         DFGRD1(I,I) = 1.0D0
         DROT(I,I)   = 1.0D0
      END DO
      CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
      OPEN(11, FILE='mp_output.txt')
      WRITE(11,'(A)') 'STRESS'
      DO I = 1, NTENS
         WRITE(11,'(E26.18)') STRESS(I)
      END DO
      WRITE(11,'(A)') 'STATEV'
      DO I = 1, MNS
         WRITE(11,'(E26.18)') STATEV(I)
      END DO
      CLOSE(11)
      END
"""


# --------------------------------------------------------------------------- #
def _run(cmd: List[str], cwd: str) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, text=True)


def _compile_transformed(build: str, plain_for: str, module_files: List[str]) -> str:
    """Compile OTI modules (free-form) + transformed UMAT (fixed-form) + driver,
    link an executable. Returns its path. Raises MatPointError with the compiler
    log on failure."""
    objs: List[str] = []
    for mod in module_files:                       # ordered: master, real_utils, otim
        obj = os.path.join(build, os.path.basename(mod) + ".o")
        r = _run(["gfortran", "-c", "-ffree-form", "-ffree-line-length-none",
                  "-I" + build, mod, "-J" + build, "-o", obj], build)
        if r.returncode != 0:
            raise MatPointError("OTI module compile failed (%s):\n%s"
                                % (os.path.basename(mod), r.stdout))
        objs.append(obj)
    umat_obj = os.path.join(build, "umat.o")
    r = _run(["gfortran", "-c", "-ffixed-form", "-ffixed-line-length-none",
              "-I" + build, plain_for, "-J" + build, "-o", umat_obj], build)
    if r.returncode != 0:
        raise MatPointError("transformed UMAT compile failed:\n%s" % r.stdout)
    objs.append(umat_obj)
    return _compile_driver_and_link(build, objs, "oti_test")


def _compile_original(build: str, original_for: str) -> str:
    umat_obj = os.path.join(build, "umat.o")
    r = _run(["gfortran", "-c", "-ffixed-form", "-ffixed-line-length-none",
              "-I" + build, original_for, "-J" + build, "-o", umat_obj], build)
    if r.returncode != 0:
        raise MatPointError("original UMAT compile failed:\n%s" % r.stdout)
    return _compile_driver_and_link(build, [umat_obj], "orig_test")


def _compile_driver_and_link(build: str, objs: List[str], exe_name: str) -> str:
    drv = os.path.join(build, "mp_driver.f")
    with open(drv, "w") as fh:
        fh.write(_driver_source())
    drv_obj = os.path.join(build, "mp_driver.o")
    r = _run(["gfortran", "-c", "-ffixed-form", "-ffixed-line-length-none",
              "-I" + build, drv, "-J" + build, "-o", drv_obj], build)
    if r.returncode != 0:
        raise MatPointError("driver compile failed:\n%s" % r.stdout)
    exe = os.path.join(build, exe_name)
    r = _run(["gfortran"] + objs + [drv_obj, "-o", exe], build)
    if r.returncode != 0:
        raise MatPointError("link failed:\n%s" % r.stdout)
    return exe


def _evaluate(exe: str, build: str, props: List[float], ntens: int,
              dstran: List[float]) -> Dict[str, List[float]]:
    lines = ["%d" % len(props), " ".join(repr(float(v)) for v in props),
             "%d" % ntens, " ".join(repr(float(v)) for v in dstran[:ntens])]
    with open(os.path.join(build, "mp_input.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    r = _run([exe], build)
    if r.returncode != 0:
        raise MatPointError("material-point run failed:\n%s" % r.stdout)
    return _parse_output(os.path.join(build, "mp_output.txt"))


def _parse_output(path: str) -> Dict[str, List[float]]:
    stress: List[float] = []
    statev: List[float] = []
    bucket = None
    with open(path) as fh:
        for raw in fh:
            tok = raw.strip()
            if tok == "STRESS":
                bucket = stress
            elif tok == "STATEV":
                bucket = statev
            elif tok and bucket is not None:
                bucket.append(float(tok.replace("D", "E").replace("d", "e")))
    return {"stress": stress, "statev": statev}


# --------------------------------------------------------------------------- #
def material_point_derivative_check(
    *,
    original_for: str,
    transformed_plain: str,
    output_dir: str,
    ntens: int,
    props: List[float],
    parameters: List[Dict[str, Any]],
    dstran: Optional[List[float]] = None,
    rel_tol: float = 1.0e-6,
    fd_step_rel: float = 1.0e-6,
    ddsdde_sdv_range: Optional[List[int]] = None,
) -> Dict[str, Any]:
    """Compile + run both UMATs and compare each SDV-extracted dsigma/da_i against
    a central finite difference of the ORIGINAL UMAT.

    parameters : ``[{"name","index","sdv_range":[lo,hi]}, ...]`` (1-based).
    ddsdde_sdv_range : optional ``[lo, hi]`` for the tangent packed row-major into
        SDV; when given, STATEV[lo..] is checked against a central FD of stress
        w.r.t. each DSTRAN component (i.e. the real DDSDDE).
    Returns a verdict dict; ``available=False`` if gfortran / OTI sources missing.
    """
    if not gfortran_available():
        return {"available": False, "reason": "gfortran not found on PATH"}
    dstran = list(dstran or _DEFAULT_DSTRAN)
    module_files = _ordered_module_files(output_dir)
    if module_files is None:
        return {"available": False,
                "reason": "OTI library sources (master_parameters/real_utils/otim*)"
                          " not found next to the transformed UMAT in %s" % output_dir}
    work = tempfile.mkdtemp(prefix="resasm_mp_")
    try:
        oti_build = os.path.join(work, "oti")
        orig_build = os.path.join(work, "orig")
        os.makedirs(oti_build)
        os.makedirs(orig_build)
        for b in (oti_build, orig_build):
            with open(os.path.join(b, "ABA_PARAM.INC"), "w") as fh:
                fh.write(_ABA_PARAM)
        # local copies so relative INCLUDE + module search resolve inside the build
        local_plain = os.path.join(oti_build, os.path.basename(transformed_plain))
        shutil.copy(transformed_plain, local_plain)
        local_mods = []
        for mod in module_files:
            dst = os.path.join(oti_build, os.path.basename(mod))
            shutil.copy(mod, dst)
            local_mods.append(dst)
        local_orig = os.path.join(orig_build, os.path.basename(original_for))
        shutil.copy(original_for, local_orig)

        # A failure to build/run the TRANSFORMED UMAT is a real transform-quality
        # defect -> passed=False (available=True). A failure of the ORIGINAL UMAT or
        # the toolchain is an environment issue -> available=False (soft skip), so an
        # ifort-only original or an old compiler does not discard a valid transform.
        try:
            oti_exe = _compile_transformed(oti_build, local_plain, local_mods)
            base = _evaluate(oti_exe, oti_build, props, ntens, dstran)
        except MatPointError as exc:
            return {"available": True, "passed": False, "max_relerr": None,
                    "reason": "transformed UMAT failed to build/run: %s" % exc,
                    "note": "transformed-side failure"}
        try:
            orig_exe = _compile_original(orig_build, local_orig)
        except MatPointError as exc:
            return {"available": False,
                    "reason": "original UMAT failed to build (reference unavailable): "
                              "%s" % exc}
        statev = base["statev"]

        try:
            results = []
            worst = 0.0
            for p in parameters:
                idx = int(p["index"])
                lo = int(p["sdv_range"][0])
                h = fd_step_rel * (abs(props[idx - 1]) or 1.0)
                pp = list(props); pp[idx - 1] += h
                pm = list(props); pm[idx - 1] -= h
                sp = _evaluate(orig_exe, orig_build, pp, ntens, dstran)["stress"]
                sm = _evaluate(orig_exe, orig_build, pm, ntens, dstran)["stress"]
                comp = []
                for i in range(ntens):
                    fd = (sp[i] - sm[i]) / (2.0 * h)
                    sdv = statev[lo - 1 + i]
                    den = abs(fd) if abs(fd) > 1e-30 else 1.0
                    rel = abs(sdv - fd) / den
                    worst = max(worst, rel)
                    comp.append({"i": i + 1, "sdv": sdv, "fd": fd, "relerr": rel})
                results.append({"name": p["name"], "index": idx, "sdv_lo": lo,
                                "max_relerr": max((c["relerr"] for c in comp), default=0.0),
                                "components": comp})

            ddsdde_check = None
            if ddsdde_sdv_range:
                # NOTE (validation blind spot): the row-major slot ntens*(i-1)+j is
                # shared by the writer (build_parameter_contracts) and this reader, so
                # a consistent transpose is invisible here; and an isotropic-elastic
                # tangent is symmetric, so even an independent transpose would compare
                # equal. Exercising the row/column convention independently needs an
                # ANISOTROPIC (non-symmetric-tangent) fixture, which we do not yet ship.
                lo = int(ddsdde_sdv_range[0])
                hd = fd_step_rel * (max(abs(v) for v in dstran[:ntens]) or 1.0)
                dd_worst = 0.0
                for j in range(ntens):                 # d sigma / d dstran_j (column j)
                    dp = list(dstran); dp[j] += hd
                    dm = list(dstran); dm[j] -= hd
                    sp = _evaluate(orig_exe, orig_build, props, ntens, dp)["stress"]
                    sm = _evaluate(orig_exe, orig_build, props, ntens, dm)["stress"]
                    for i in range(ntens):             # row-major SDV(lo + ntens*i + j)
                        fd = (sp[i] - sm[i]) / (2.0 * hd)
                        sdv = statev[lo - 1 + ntens * i + j]
                        den = abs(fd) if abs(fd) > 1e-30 else 1.0
                        dd_worst = max(dd_worst, abs(sdv - fd) / den)
                worst = max(worst, dd_worst)
                ddsdde_check = {"sdv_lo": lo, "max_relerr": dd_worst}
        except MatPointError as exc:               # ORIGINAL-side evaluation failure
            return {"available": False,
                    "reason": "original UMAT reference run failed: %s" % exc}

        return {
            "available": True,
            "passed": bool(worst <= rel_tol),
            "max_relerr": worst,
            "rel_tol": rel_tol,
            "fd_step_rel": fd_step_rel,
            "parameters": results,
            "ddsdde": ddsdde_check,
            "note": "dsigma/da_i and DDSDDE (SDV) vs central FD of the ORIGINAL UMAT",
        }
    finally:
        shutil.rmtree(work, ignore_errors=True)


def _ordered_module_files(output_dir: str) -> Optional[List[str]]:
    master = os.path.join(output_dir, "master_parameters.f90")
    real_utils = os.path.join(output_dir, "real_utils.f90")
    otim = sorted(glob.glob(os.path.join(output_dir, "otim*n1.f90")))
    if not (os.path.exists(master) and os.path.exists(real_utils) and otim):
        return None
    return [master, real_utils, otim[0]]
