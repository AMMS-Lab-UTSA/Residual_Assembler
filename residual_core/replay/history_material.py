"""Batched access to a compiled provider object for whole-model history replay.

The collaborator receives ``OTI_UMAT.obj`` plus its completed contract
(``resasm_umat_oti_contract_v1``). This module links that object with a small
Fortran shim into a private shared library and exposes two batched calls, one
ctypes call per increment for every integration point of the mesh:

``total``
    ``UMAT_OTI_EVAL_TOTAL`` at every point. Parameter direction ``j`` carries
    the derivative of every history input (incoming STRESS, STATEV, the strain
    at the start of the increment STRAN and the strain increment DSTRAN) and
    the unit PROPS seed, so the returned DSIGMA_DP/DSTATEV_DP are TOTAL
    first-order derivatives through the update. It also returns DDSDDE and
    dSTATEV/dDSTRAN from the unit DSTRAN directions.

``regular``
    the ORIGINAL (real-valued) UMAT compiled into the same object. It never
    touches OTI arithmetic, so a solution equilibrated with it is an
    independent reference for finite differences.

Nothing here knows which constitutive model is inside: dimensions, the
parameter order and the PROPS slots come from the contract. The object's hash
is checked against the contract before linking. Array layout on the Python
side is C order with the point axis first: stress (n, NTENS), state
(n, NSTATV), DSIGMA_DP (n, NTENS, NPARAM), DSTATEV_DP (n, NSTATV, NPARAM),
DDSDDE (n, NTENS, NTENS) with ``[i, a, b] = dSTRESS_a/dDSTRAN_b``,
dSTATEV/dDSTRAN (n, NSTATV, NTENS).
"""
from __future__ import annotations

import ctypes
import hashlib
import os
import subprocess
import tempfile
from typing import Any, Mapping

import numpy as np

from residual_core.runtime import (
    check_binary_compatibility,
    load_shared_library,
    prefer_static_runtime,
    shared_library_suffix,
    static_fortran_link_flags,
)

TOTAL_SYMBOL = "umat_oti_eval_total_"
TOTAL_SIGNATURE_LENGTH = 28

_SHIM = r"""
subroutine hist_total(n, nprops, ntens, nstatv, nparam, props, stress, statev, &
    stran, dstran, time2, dtime, temp, dtemp, dsig_in, dstv_in, stran_dp, &
    dstran_dp, coords, celent, noel, npt, kstep, kinc, ddsdde, dsig, dstv, &
    dstv_de, pnewdt) bind(C, name="hist_total")
  use iso_c_binding
  implicit none
  integer(c_int), value :: n, nprops, ntens, nstatv, nparam, kstep, kinc
  real(c_double), intent(in) :: props(nprops), time2(2)
  real(c_double), value :: dtime, temp, dtemp
  real(c_double), intent(inout) :: stress(ntens, n), statev(max(nstatv, 1), n)
  real(c_double), intent(in) :: stran(ntens, n), dstran(ntens, n), coords(3, n), celent(n)
  real(c_double), intent(in) :: dsig_in(nparam, ntens, n), dstv_in(nparam, max(nstatv, 1), n)
  real(c_double), intent(in) :: stran_dp(nparam, ntens, n), dstran_dp(nparam, ntens, n)
  integer(c_int), intent(in) :: noel(n), npt(n)
  real(c_double), intent(out) :: ddsdde(ntens, ntens, n), dsig(nparam, ntens, n)
  real(c_double), intent(out) :: dstv(nparam, max(nstatv, 1), n)
  real(c_double), intent(out) :: dstv_de(ntens, max(nstatv, 1), n), pnewdt(n)
  real(8) :: s(ntens), sv(max(nstatv, 1)), dd(ntens, ntens), ds(ntens, nparam)
  real(8) :: dv(max(nstatv, 1), nparam), dsi(ntens, nparam), dvi(max(nstatv, 1), nparam)
  real(8) :: e0(ntens, nparam), de(ntens, nparam), dve(max(nstatv, 1), ntens), pn
  integer :: i
  do i = 1, n
    s = stress(:, i)
    sv = 0.0d0
    if (nstatv > 0) sv(1:nstatv) = statev(1:nstatv, i)
    dsi = transpose(dsig_in(:, :, i))
    dvi = 0.0d0
    if (nstatv > 0) dvi(1:nstatv, :) = transpose(dstv_in(:, 1:nstatv, i))
    e0 = transpose(stran_dp(:, :, i))
    de = transpose(dstran_dp(:, :, i))
    call umat_oti_eval_total(s, sv, dd, stran(:, i), dstran(:, i), time2, dtime, &
        temp, dtemp, props, nprops, ntens, nstatv, nparam, ds, dv, dsi, dvi, &
        e0, de, dve, coords(:, i), celent(i), noel(i), npt(i), kstep, kinc, pn)
    stress(:, i) = s
    if (nstatv > 0) statev(1:nstatv, i) = sv(1:nstatv)
    ddsdde(:, :, i) = transpose(dd)
    dsig(:, :, i) = transpose(ds)
    dstv(:, :, i) = 0.0d0
    dstv_de(:, :, i) = 0.0d0
    if (nstatv > 0) then
      dstv(:, 1:nstatv, i) = transpose(dv(1:nstatv, :))
      dstv_de(:, 1:nstatv, i) = transpose(dve(1:nstatv, :))
    end if
    pnewdt(i) = pn
  end do
end subroutine hist_total

subroutine hist_regular(n, nprops, ntens, nstatv, props, stress, statev, stran, &
    dstran, time2, dtime, temp, dtemp, coords, celent, noel, npt, kstep, kinc, &
    ddsdde, pnewdt) bind(C, name="hist_regular")
  use iso_c_binding
  implicit none
  integer(c_int), value :: n, nprops, ntens, nstatv, kstep, kinc
  real(c_double), intent(in) :: props(nprops), time2(2)
  real(c_double), value :: dtime, temp, dtemp
  real(c_double), intent(inout) :: stress(ntens, n), statev(max(nstatv, 1), n)
  real(c_double), intent(in) :: stran(ntens, n), dstran(ntens, n), coords(3, n), celent(n)
  integer(c_int), intent(in) :: noel(n), npt(n)
  real(c_double), intent(out) :: ddsdde(ntens, ntens, n), pnewdt(n)
  real(8) :: s(ntens), sv(max(nstatv, 1)), dd(ntens, ntens), pr(nprops)
  real(8) :: sse, spd, scd, rpl, ddsddt(ntens), drplde(ntens), drpldt
  real(8) :: predef(1), dpred(1), drot(3, 3), dfgrd0(3, 3), dfgrd1(3, 3), pn, cel
  real(8) :: st(ntens), dst(ntens), tm(2), co(3), dtm, tp, dtp
  character(len=80) :: cmname
  integer :: i, k
  cmname = 'MATERIAL'
  do i = 1, n
    s = stress(:, i)
    sv = 0.0d0
    if (nstatv > 0) sv(1:nstatv) = statev(1:nstatv, i)
    pr = props
    st = stran(:, i); dst = dstran(:, i); tm = time2; co = coords(:, i)
    dtm = dtime; tp = temp; dtp = dtemp; cel = celent(i)
    sse = 0.0d0; spd = 0.0d0; scd = 0.0d0; rpl = 0.0d0; drpldt = 0.0d0
    ddsddt = 0.0d0; drplde = 0.0d0; predef = 0.0d0; dpred = 0.0d0; dd = 0.0d0
    drot = 0.0d0; dfgrd0 = 0.0d0; dfgrd1 = 0.0d0
    do k = 1, 3
      drot(k, k) = 1.0d0; dfgrd0(k, k) = 1.0d0; dfgrd1(k, k) = 1.0d0
    end do
    pn = 1.0d0
    call umat(s, sv, dd, sse, spd, scd, rpl, ddsddt, drplde, drpldt, st, dst, tm, &
        dtm, tp, dtp, predef, dpred, cmname, 3, ntens - 3, ntens, max(nstatv, 1), &
        pr, nprops, co, drot, pn, cel, dfgrd0, dfgrd1, noel(i), npt(i), 1, 1, kstep, kinc)
    stress(:, i) = s
    if (nstatv > 0) statev(1:nstatv, i) = sv(1:nstatv)
    ddsdde(:, :, i) = transpose(dd)
    pnewdt(i) = pn
  end do
end subroutine hist_regular
"""


class ProviderError(ValueError):
    """The provider object or its contract cannot serve a history replay."""


def _digest(path: str) -> str:
    with open(path, "rb") as stream:
        return hashlib.sha256(stream.read()).hexdigest()


def _pointer(array: np.ndarray):
    return array.ctypes.data_as(ctypes.c_void_p)


class HistoryMaterial:
    """A provider object linked for batched total-derivative evaluation."""

    def __init__(self, obj_path: str, contract: Mapping[str, Any], workdir: str | None = None,
                 *, compiler: str | None = None):
        if contract.get("schema") != "resasm_umat_oti_contract_v1":
            raise ProviderError("expected the completed provider contract resasm_umat_oti_contract_v1")
        if contract.get("kinematics") != "small_strain":
            raise ProviderError("history replay supports small_strain providers only")
        obj_path = os.path.abspath(str(obj_path))
        full = contract.get("object", {}).get("sha256_full")
        if not full:
            raise ProviderError("contract object.sha256_full is required to bind the object")
        if _digest(obj_path) != full:
            raise ProviderError("OTI object sha256 does not match contract object.sha256_full")
        symbols = contract.get("symbols", {})
        if symbols.get("oti_eval_total") != TOTAL_SYMBOL or len(
                symbols.get("oti_eval_total_signature", [])) != TOTAL_SIGNATURE_LENGTH:
            raise ProviderError(
                "provider lacks UMAT_OTI_EVAL_TOTAL (the total-derivative entry point); rebuild it "
                "with `umat-oti-provider build <contract_v2.json> --out DIR` from a UMAT-OTI version "
                "that emits oti_eval_total")
        check_binary_compatibility(obj_path, contract.get("binary"))
        dimensions = contract["dimensions"]
        self.ntens = int(dimensions["ntens"])
        self.nprops = int(dimensions["nprops"])
        self.nstatev = int(dimensions["nstatev"])
        self.nparam = int(dimensions["nparam"])
        if self.ntens != 6:
            raise ProviderError("history replay requires 3D NTENS=6")
        parameters = contract.get("parameters", [])
        if len(parameters) != self.nparam or [p.get("oti_direction") for p in parameters] != list(
                range(1, self.nparam + 1)):
            raise ProviderError("contract parameters must list OTI directions 1..NPARAM in order")
        self.params = [str(p["name"]) for p in parameters]
        self.props_index = [int(p["props_index"]) for p in parameters]
        if len(set(self.props_index)) != self.nparam or not all(
                1 <= slot <= self.nprops for slot in self.props_index):
            raise ProviderError("parameter PROPS indices must be unique and within NPROPS")
        self.contract = dict(contract)
        self.object_path = obj_path
        self.object_sha256 = full
        self.work = workdir or tempfile.mkdtemp(prefix="resasm_history_")
        os.makedirs(self.work, exist_ok=True)
        self._compiler = compiler or os.environ.get("FC") or "gfortran"
        self._handle = load_shared_library(self._link())
        self.lib = self._handle.lib
        self.lib.hist_total.restype = None
        self.lib.hist_regular.restype = None

    def _link(self) -> str:
        source = os.path.join(self.work, "history_shim.f90")
        with open(source, "w") as stream:
            stream.write(_SHIM)
        shim = os.path.join(self.work, "history_shim.o")
        library = os.path.join(self.work, "history_material" + shared_library_suffix())
        run = dict(check=True, capture_output=True, text=True)
        try:
            subprocess.run([self._compiler, "-O2", "-fPIC", "-ffree-form", "-c", source, "-o", shim], **run)
            link = [self._compiler, "-shared"]
            if prefer_static_runtime():
                link += static_fortran_link_flags()
            subprocess.run(link + [shim, self.object_path, "-o", library], **run)
        except subprocess.CalledProcessError as error:
            raise ProviderError("linking the provider object failed:\n%s\n%s" % (
                " ".join(error.cmd), (error.stdout or "") + (error.stderr or ""))) from error
        return library

    # ------------------------------------------------------------------ utils
    def props_with(self, props, values=None):
        """PROPS with the named parameters replaced (for perturbed references)."""
        props = np.array(props, dtype=float)
        if props.shape != (self.nprops,):
            raise ProviderError("PROPS length %d does not match contract NPROPS %d" % (len(props), self.nprops))
        for name, value in (values or {}).items():
            props[self.props_index[self.params.index(name)] - 1] = value
        return props

    def parameter_values(self, props):
        props = np.asarray(props, dtype=float)
        return np.array([props[slot - 1] for slot in self.props_index])

    def _common(self, props, stress, state, stran, dstran, coords, celent, noel, npt):
        n = stress.shape[0]
        props = np.ascontiguousarray(props, dtype=float)
        if props.shape != (self.nprops,) or not np.all(np.isfinite(props)):
            raise ProviderError("PROPS must be %d finite reals" % self.nprops)
        arrays = {
            "stress": (stress, (n, self.ntens)), "state": (state, (n, self.nstatev)),
            "stran": (stran, (n, self.ntens)), "dstran": (dstran, (n, self.ntens)),
            "coords": (coords, (n, 3)), "celent": (celent, (n,)),
        }
        out = {}
        for name, (value, shape) in arrays.items():
            value = np.ascontiguousarray(value, dtype=float)
            if value.shape != shape or not np.all(np.isfinite(value)):
                raise ProviderError("%s must be finite with shape %s, got %s" % (name, shape, value.shape))
            out[name] = value.copy() if name in ("stress", "state") else value
        if self.nstatev == 0:
            out["state_buffer"] = np.zeros((n, 1))
        else:
            out["state_buffer"] = out["state"]
        out["noel"] = np.ascontiguousarray(noel, dtype=np.int32)
        out["npt"] = np.ascontiguousarray(npt, dtype=np.int32)
        if out["noel"].shape != (n,) or out["npt"].shape != (n,):
            raise ProviderError("NOEL/NPT must have one entry per point")
        return n, props, out

    # ------------------------------------------------------------- the calls
    def total(self, props, stress, state, stran, dstran, *, dstress_in, dstate_in,
              stran_dp, dstran_dp, time, dtime, coords, celent, noel, npt,
              kstep=1, kinc=1, temp=0.0, dtemp=0.0):
        """Batched UMAT_OTI_EVAL_TOTAL; returns a dict of updated arrays."""
        n, props, a = self._common(props, stress, state, stran, dstran, coords, celent, noel, npt)
        nt, ns, npar = self.ntens, self.nstatev, self.nparam
        seeds = {}
        for name, value, shape in (("dstress_in", dstress_in, (n, nt, npar)),
                                   ("dstate_in", dstate_in, (n, ns, npar)),
                                   ("stran_dp", stran_dp, (n, nt, npar)),
                                   ("dstran_dp", dstran_dp, (n, nt, npar))):
            value = np.asarray(value, dtype=float)
            if value.shape != shape or not np.all(np.isfinite(value)):
                raise ProviderError("%s must be finite with shape %s" % (name, shape))
            seeds[name] = np.ascontiguousarray(value)
        dstate_in = seeds["dstate_in"] if ns else np.zeros((n, 1, npar))
        time2 = np.ascontiguousarray(time, dtype=float)
        if time2.shape != (2,) or not np.isfinite(dtime) or dtime <= 0:
            raise ProviderError("TIME must be two reals and DTIME finite and positive")
        ddsdde = np.empty((n, nt, nt))
        dsig = np.empty((n, nt, npar))
        dstv = np.empty((n, max(ns, 1), npar))
        dstv_de = np.empty((n, max(ns, 1), nt))
        pnewdt = np.empty(n)
        ci, cd = ctypes.c_int, ctypes.c_double
        self.lib.hist_total(
            ci(n), ci(self.nprops), ci(nt), ci(ns), ci(npar), _pointer(props),
            _pointer(a["stress"]), _pointer(a["state_buffer"]), _pointer(a["stran"]),
            _pointer(a["dstran"]), _pointer(time2), cd(dtime), cd(temp), cd(dtemp),
            _pointer(seeds["dstress_in"]), _pointer(np.ascontiguousarray(dstate_in)),
            _pointer(seeds["stran_dp"]), _pointer(seeds["dstran_dp"]), _pointer(a["coords"]),
            _pointer(a["celent"]), _pointer(a["noel"]), _pointer(a["npt"]), ci(kstep), ci(kinc),
            _pointer(ddsdde), _pointer(dsig), _pointer(dstv), _pointer(dstv_de), _pointer(pnewdt))
        self._check(pnewdt, "OTI")
        return {"stress": a["stress"], "state": a["state_buffer"][:, :ns] if ns else np.zeros((n, 0)),
                "ddsdde": ddsdde, "dstress_dp": dsig, "dstate_dp": dstv[:, :ns, :],
                "dstate_ddstran": dstv_de[:, :ns, :]}

    def regular(self, props, stress, state, stran, dstran, *, time, dtime, coords, celent,
                noel, npt, kstep=1, kinc=1, temp=0.0, dtemp=0.0):
        """Batched ORIGINAL UMAT; returns stress, state, DDSDDE (as the UMAT returns it)."""
        n, props, a = self._common(props, stress, state, stran, dstran, coords, celent, noel, npt)
        nt, ns = self.ntens, self.nstatev
        time2 = np.ascontiguousarray(time, dtype=float)
        if time2.shape != (2,) or not np.isfinite(dtime) or dtime <= 0:
            raise ProviderError("TIME must be two reals and DTIME finite and positive")
        ddsdde = np.empty((n, nt, nt))
        pnewdt = np.empty(n)
        ci, cd = ctypes.c_int, ctypes.c_double
        self.lib.hist_regular(
            ci(n), ci(self.nprops), ci(nt), ci(ns), _pointer(props), _pointer(a["stress"]),
            _pointer(a["state_buffer"]), _pointer(a["stran"]), _pointer(a["dstran"]),
            _pointer(time2), cd(dtime), cd(temp), cd(dtemp), _pointer(a["coords"]),
            _pointer(a["celent"]), _pointer(a["noel"]), _pointer(a["npt"]), ci(kstep), ci(kinc),
            _pointer(ddsdde), _pointer(pnewdt))
        self._check(pnewdt, "ORIGINAL")
        return {"stress": a["stress"], "state": a["state_buffer"][:, :ns] if ns else np.zeros((n, 0)),
                "ddsdde": ddsdde}

    @staticmethod
    def _check(pnewdt, which):
        if not np.all(np.isfinite(pnewdt)) or np.any(pnewdt < 1.0):
            bad = int(np.sum(~np.isfinite(pnewdt) | (pnewdt < 1.0)))
            raise ProviderError("%s UMAT requested a time cut-back (PNEWDT<1) at %d points; the "
                                "recorded increment cannot be replayed as saved" % (which, bad))
