"""Program 2 path-dependent material replay: march a distributed OTI .obj over a
recorded strain path, carrying the physical state AND its parameter sensitivity.

For a path-dependent material the residual method must replay increment by
increment because the stress update depends on the previous stress/state. This
module links the opaque ``umat_<name>_oti.obj`` into a small ``.so`` exposing two
C entry points and marches them in Python:

  oti_path : calls UMAT_OTI_EVAL with the incoming DSIGMA_DP_IN / DSTATEV_DP_IN
             and returns the TOTAL DSIGMA_DP / DSTATEV_DP for the increment
             (the eval chains the state-transition Jacobians internally).
  reg_path : calls the REGULAR UMAT so an independent finite-difference reference
             can be marched over the same path (no hard-coded constitutive law).

No source transformation and no OTI promotion happen here -- this only links and
replays a binary the collaborator received, exactly like objlink.
"""
from __future__ import annotations

import ctypes
import hashlib
import os
import subprocess
import tempfile
from typing import Any, List, Mapping, Sequence, Tuple

import numpy as np

from residual_core.runtime import (
    LibraryLoadError,
    check_binary_compatibility,
    load_shared_library,
    prefer_static_runtime,
    shared_library_suffix,
    static_fortran_link_flags,
)

_D = ctypes.c_double
_IP = ctypes.POINTER(_D)


def _sha(path: str) -> str:
    return hashlib.sha256(open(path, "rb").read()).hexdigest()[:16]


_PATH_SHIM = r"""
      subroutine oti_path(props,np,strin,statin,ns,dsigin,dstvin,
     1 dstran,dt,nt,nprm,strout,ddsdde,statout,dsigout,dstvout)
     2 bind(C,name="oti_path")
      use iso_c_binding
      real(c_double),intent(in)::props(*),strin(*),statin(*),dstran(*)
      real(c_double),intent(in)::dsigin(*),dstvin(*)
      real(c_double),value::dt
      integer(c_int),value::np,ns,nt,nprm
      real(c_double),intent(out)::strout(*),ddsdde(*),statout(*),dsigout(*),dstvout(*)
      real*8 STRESS(nt),STATEV(max(ns,1)),DDSDDE_(nt,nt),STRAN(nt),DSTRN(nt)
      real*8 TIME(2),PROPS_(np),DSIGMA_DP(nt,nprm),DSTATEV_DP(max(ns,1),nprm)
      real*8 DSIGMA_DP_IN(nt,nprm),DSTATEV_DP_IN(max(ns,1),nprm)
      integer i,j
      do i=1,nt; STRESS(i)=strin(i); STRAN(i)=0.d0; DSTRN(i)=dstran(i); end do
      do i=1,max(ns,1); STATEV(i)=0.d0; end do
      do i=1,ns; STATEV(i)=statin(i); end do
      do i=1,np; PROPS_(i)=props(i); end do
      do j=1,nprm; do i=1,nt; DSIGMA_DP_IN(i,j)=dsigin((i-1)*nprm+j); end do; end do
      do j=1,nprm; do i=1,ns; DSTATEV_DP_IN(i,j)=dstvin((i-1)*nprm+j); end do; end do
      TIME(1)=0.d0; TIME(2)=0.d0
      call UMAT_OTI_EVAL(STRESS,STATEV,DDSDDE_,STRAN,DSTRN,TIME,dt,0.d0,0.d0,
     1 PROPS_,np,nt,ns,nprm,DSIGMA_DP,DSTATEV_DP,DSIGMA_DP_IN,DSTATEV_DP_IN)
      do i=1,nt; strout(i)=STRESS(i); end do
      do i=1,nt; do j=1,nt; ddsdde((i-1)*nt+j)=DDSDDE_(i,j); end do; end do
      do i=1,ns; statout(i)=STATEV(i); end do
      do j=1,nprm; do i=1,nt; dsigout((i-1)*nprm+j)=DSIGMA_DP(i,j); end do; end do
      do j=1,nprm; do i=1,ns; dstvout((i-1)*nprm+j)=DSTATEV_DP(i,j); end do; end do
      end subroutine

      subroutine reg_path(props,np,strin,statin,ns,dstran,dt,nt,
     1 strout,ddsdde,statout) bind(C,name="reg_path")
      use iso_c_binding
      real(c_double),intent(in)::props(*),strin(*),statin(*),dstran(*)
      real(c_double),value::dt
      integer(c_int),value::np,ns,nt
      real(c_double),intent(out)::strout(*),ddsdde(*),statout(*)
      real*8 STRESS(nt),STATEV(max(ns,1)),DDSDDE_(nt,nt),STRAN(nt),DSTRN(nt)
      real*8 TIME(2),PREDEF(1),DPRED(1),PROPS_(np),COORDS(3),DROT(3,3)
      real*8 DFGRD0(3,3),DFGRD1(3,3),DDSDDT(nt),DRPLDE(nt)
      real*8 SSE,SPD,SCD,RPL,DRPLDT
      character*80 CMNAME
      integer i,j
      SSE=0.d0; SPD=0.d0; SCD=0.d0; RPL=0.d0; DRPLDT=0.d0
      do i=1,nt; STRESS(i)=strin(i); STRAN(i)=0.d0; DSTRN(i)=dstran(i)
        DDSDDT(i)=0.d0; DRPLDE(i)=0.d0; end do
      do i=1,max(ns,1); STATEV(i)=0.d0; end do
      do i=1,ns; STATEV(i)=statin(i); end do
      do i=1,np; PROPS_(i)=props(i); end do
      do i=1,3; COORDS(i)=0.d0; do j=1,3; DROT(i,j)=0.d0
        DFGRD0(i,j)=0.d0; DFGRD1(i,j)=0.d0; end do
        DROT(i,i)=1.d0; DFGRD0(i,i)=1.d0; DFGRD1(i,i)=1.d0; end do
      TIME(1)=0.d0; TIME(2)=0.d0
      call UMAT(STRESS,STATEV,DDSDDE_,SSE,SPD,SCD,RPL,DDSDDT,DRPLDE,DRPLDT,
     1 STRAN,DSTRN,TIME,dt,0.d0,0.d0,PREDEF,DPRED,CMNAME,
     2 3,nt-3,nt,max(ns,1),PROPS_,np,COORDS,DROT,1.d0,1.d0,DFGRD0,DFGRD1,
     3 1,1,1,1,1,1)
      do i=1,nt; strout(i)=STRESS(i); end do
      do i=1,nt; do j=1,nt; ddsdde((i-1)*nt+j)=DDSDDE_(i,j); end do; end do
      do i=1,ns; statout(i)=STATEV(i); end do
      end subroutine
"""

# march_c is only linkable against objects that expose UMAT_OTI_MARCH (built with
# the state-alive whole-path march). Appended to the shim only when the completed
# contract declares a "march" symbol, so elastic/no-march objects still link.
_MARCH_SHIM = r"""
      subroutine march_c(props,np,path,npath,dt,nt,ns,nprm,dsig,strout,
     1 ddout) bind(C,name="march_c")
      use iso_c_binding
      real(c_double),intent(in)::props(*),path(*),dt(*)
      integer(c_int),value::np,npath,nt,ns,nprm
      real(c_double),intent(out)::dsig(*),strout(*),ddout(*)
      call UMAT_OTI_MARCH(props,np,path,npath,dt,nt,ns,nprm,dsig,strout,
     1 ddout)
      end subroutine
"""


class PathMaterial:
    """Linked path-dependent OTI material (single integration point marcher)."""

    def __init__(self, obj_path: str, contract: Mapping[str, Any], workdir: str = None):
        if contract.get("schema") != "resasm_umat_oti_contract_v1":
            raise ValueError("expected completed contract schema resasm_umat_oti_contract_v1")
        declared = contract.get("object", {}).get("sha256")
        if declared and _sha(obj_path) != declared:
            raise ValueError("object hash != contract object.sha256")
        # reject a platform-incompatible object before linking (clear diagnostic
        # instead of an opaque ld error on the wrong-format object).
        check_binary_compatibility(obj_path, contract.get("binary"))
        d = contract["dimensions"]
        self.ntens = int(d["ntens"]); self.nprops = int(d["nprops"])
        self.nstatev = int(d["nstatev"]); self.nparam = int(d["nparam"])
        self.params = [p["name"] for p in contract["parameters"]]
        self.props_index = [int(p["props_index"]) for p in contract["parameters"]]
        self.work = workdir or tempfile.mkdtemp(prefix="resasm_path_")
        # the completed contract stores "march": null for stateless materials, so
        # a truthy test is required -- "march" in contract is True even when the
        # object exposes no march symbol.
        self._want_march = bool(contract.get("march"))
        try:
            so = self._link(obj_path, self._want_march)
            self._handle = load_shared_library(so)
        except (LibraryLoadError, OSError, subprocess.CalledProcessError):
            # the object doesn't expose umat_oti_march_ (an older/partial build):
            # linking or loading the march shim fails. Fall back to the no-march
            # shim -- the per-increment oti_path/reg_path marches still work, only
            # the fast whole-path march is lost. If the no-march shim also fails,
            # that error propagates (it is a genuine build/runtime problem).
            if not self._want_march:
                raise
            so = self._link(obj_path, False)
            self._handle = load_shared_library(so)
            self._want_march = False
        self.lib = self._handle.lib
        self.lib.oti_path.argtypes = [_IP, ctypes.c_int, _IP, _IP, ctypes.c_int, _IP, _IP,
                                      _IP, _D, ctypes.c_int, ctypes.c_int, _IP, _IP, _IP, _IP, _IP]
        self.lib.reg_path.argtypes = [_IP, ctypes.c_int, _IP, _IP, ctypes.c_int, _IP,
                                      _D, ctypes.c_int, _IP, _IP, _IP]
        self.has_march = hasattr(self.lib, "march_c")
        if self.has_march:
            self.lib.march_c.argtypes = [_IP, ctypes.c_int, _IP, ctypes.c_int, _IP,
                                         ctypes.c_int, ctypes.c_int, ctypes.c_int, _IP, _IP, _IP]

    def _link(self, obj_path: str, with_march: bool = True) -> str:
        os.makedirs(self.work, exist_ok=True)
        src = os.path.join(self.work, "path_shim.for")
        shim = _PATH_SHIM + (_MARCH_SHIM if with_march else "")
        open(src, "w").write(shim)
        o = os.path.join(self.work, "path_shim.o")
        fc = os.environ.get("FC") or "gfortran"
        # capture output: a link failure here is caught by __init__ and retried
        # without the march shim, so the linker's noise is not surfaced unless the
        # retry also fails.
        subprocess.run([fc, "-c", "-fPIC", "-ffixed-form", "-ffixed-line-length-none",
                        src, "-o", o], check=True, capture_output=True, text=True)
        so = os.path.join(self.work, "path_material" + shared_library_suffix())
        link = [fc, "-shared", o, obj_path, "-o", so]
        if prefer_static_runtime():
            link = [fc, "-shared"] + static_fortran_link_flags() + [o, obj_path, "-o", so]
        subprocess.run(link, check=True, capture_output=True, text=True)
        return so

    # --- single increment -------------------------------------------------
    def _oti_step(self, props, strin, statin, dsigin, dstvin, dstran, dt):
        nt, ns, nprm = self.ntens, self.nstatev, self.nparam
        pr = (_D * self.nprops)(*props); si = (_D * nt)(*strin)
        st = (_D * max(ns, 1))(*(list(statin) + [0.0] * (max(ns, 1) - ns)))
        de = (_D * nt)(*dstran)
        dsi = (_D * (nt * nprm))(*np.asarray(dsigin).flatten())
        dsv = (_D * (ns * nprm))(*np.asarray(dstvin).flatten())
        so = (_D * nt)(); dd = (_D * (nt * nt))(); sto = (_D * max(ns, 1))()
        dso = (_D * (nt * nprm))(); dsvo = (_D * (ns * nprm))()
        self.lib.oti_path(pr, self.nprops, si, st, ns, dsi, dsv, de, dt, nt, nprm,
                          so, dd, sto, dso, dsvo)
        return (np.array(so), np.array(sto)[:ns], np.array(dd).reshape(nt, nt),
                np.array(dso).reshape(nt, nprm), np.array(dsvo).reshape(ns, nprm))

    def _reg_step(self, props, strin, statin, dstran, dt):
        nt, ns = self.ntens, self.nstatev
        pr = (_D * self.nprops)(*props); si = (_D * nt)(*strin)
        st = (_D * max(ns, 1))(*(list(statin) + [0.0] * (max(ns, 1) - ns)))
        de = (_D * nt)(*dstran)
        so = (_D * nt)(); dd = (_D * (nt * nt))(); sto = (_D * max(ns, 1))()
        self.lib.reg_path(pr, self.nprops, si, st, ns, de, dt, nt, so, dd, sto)
        return np.array(so), np.array(sto)[:ns], np.array(dd).reshape(nt, nt)

    # --- full-path marches ------------------------------------------------
    def march_oti(self, props, dstran_path, dt_path):
        """March the OTI eval; return per-increment dict with stress, statev,
        dsigma_dp (nt,nparam), dstatev_dp (nstatev,nparam) -- TOTAL derivatives."""
        nt, ns, nprm = self.ntens, self.nstatev, self.nparam
        s = np.zeros(nt); st = np.zeros(ns)
        dsig = np.zeros((nt, nprm)); dstv = np.zeros((ns, nprm))
        out = []
        for de, dt in zip(dstran_path, dt_path):
            s, st, dd, dsig, dstv = self._oti_step(props, s, st, dsig, dstv, de, dt)
            out.append({"stress": s.copy(), "statev": st.copy(), "ddsdde": dd.copy(),
                        "dsigma_dp": dsig.copy(), "dstatev_dp": dstv.copy()})
        return out

    def march_fast(self, props, dstran_path, dt_path):
        """EFFICIENT state-alive whole-path march (UMAT_OTI_MARCH, param-only OTI):
        the whole path runs in ONE Fortran call with the OTI state carried across
        increments. Returns per-increment {dsigma_dp (nt,nparam), ddsdde (nt,nt)}.
        ~2x faster than march_oti and gives the same DSIGMA_DP to machine precision."""
        if not getattr(self, "has_march", False):
            raise RuntimeError("this .obj has no UMAT_OTI_MARCH symbol (rebuild with the march)")
        nt, ns, nprm = self.ntens, self.nstatev, self.nparam
        npath = len(dstran_path)
        path = np.asarray(dstran_path, float).T                 # (nt, npath)
        pr = (_D * self.nprops)(*props)
        pa = (_D * (nt * npath))(*path.flatten("F"))
        dtc = (_D * npath)(*dt_path)
        ds = (_D * (nt * nprm * npath))(); so = (_D * nt)(); dd = (_D * (nt * nt * npath))()
        self.lib.march_c(pr, self.nprops, pa, npath, dtc, nt, ns, nprm, ds, so, dd)
        dsig = np.array(ds).reshape(nt, nprm, npath, order="F")
        ddm = np.array(dd).reshape(nt, nt, npath, order="F")
        return [{"dsigma_dp": dsig[:, :, k], "ddsdde": ddm[:, :, k]} for k in range(npath)]

    def march_regular(self, props, dstran_path, dt_path):
        """March the REGULAR UMAT; return per-increment stress + statev arrays."""
        nt, ns = self.ntens, self.nstatev
        s = np.zeros(nt); st = np.zeros(ns); out = []
        for de, dt in zip(dstran_path, dt_path):
            s, st, _ = self._reg_step(props, s, st, de, dt)
            out.append({"stress": s.copy(), "statev": st.copy()})
        return out


def mises(sig: Sequence[float]) -> float:
    s = np.asarray(sig)
    return float(np.sqrt(0.5 * ((s[0] - s[1]) ** 2 + (s[1] - s[2]) ** 2 + (s[2] - s[0]) ** 2)
                         + 3.0 * (s[3] ** 2 + s[4] ** 2 + s[5] ** 2)))


def dmises_dstress(sig: Sequence[float]) -> np.ndarray:
    """d(sigma_vM)/d(sigma) in Voigt (11,22,33,12,13,23), engineering-shear stress."""
    s = np.asarray(sig, float); sm = mises(s)
    if sm < 1e-30:
        return np.zeros(6)
    hydro = (s[0] + s[1] + s[2]) / 3.0
    dev = np.array([s[0] - hydro, s[1] - hydro, s[2] - hydro, s[3], s[4], s[5]])
    g = np.empty(6)
    g[:3] = 1.5 * dev[:3] / sm
    g[3:] = 3.0 * dev[3:] / sm
    return g
