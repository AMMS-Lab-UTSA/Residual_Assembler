"""A verified corpus UMAT, built live as a parameter-seeded OTI provider.

The Residual Assembler is asked here to assemble residuals from UMATs it did
not choose: the corpus cases the UMAT-OTI pipeline verified in Abaqus. This
module gets one of those routines into the process, twice over, from ONE
compile of ONE source:

``regular``
    the ORIGINAL routine, compiled unchanged (bar the declaration-only
    adaptation below, recorded when made). It never sees OTI arithmetic, so
    finite differences of it are an independent reference.

``total``
    the PROPS-seeded OTI lift of the same routine, through the provider's
    ``UMAT_OTI_EVAL_TOTAL`` (small strain, F = I) or ``UMAT_OTI_EVAL_TOTAL_F``
    (deformation-gradient driven) entry point. Parameter direction ``j``
    carries the unit PROPS seed plus whatever incoming seeds the caller gives
    (stress, state, STRAN, DSTRAN, DFGRD0, DFGRD1), so first-order outputs are
    linear in those seeds -- which is what the assembler's local, total and
    tangent computations all use.

The build is ``umat_oti.provider.build.build_provider`` from the UMAT
repository, imported from wherever ``PYTHONPATH`` points; the module refuses
to run against an installed ``umat_oti`` that is not the checkout the caller
named (the workspace venv carries a stale one).

Source adaptation (stated, never silent)
----------------------------------------
The provider's wrappers pass a scalar KSTEP in the 4th-from-last UMAT slot.
A source that declares that dummy as ``JSTEP(4)`` (the Abaqus 2021+ spelling)
then fails to compile against the lifted module's explicit interface with
"Rank mismatch in argument 'jstep'". When -- and only when -- ``JSTEP`` occurs
nowhere but in the argument list and its one declaration, the copy built here
renames the dummy to ``KSTEP`` and moves the ``(4)`` array onto an unused
local. Since the routine never reads JSTEP the arithmetic is unchanged; the
before/after digests and the edit are recorded in the build record. Any other
use of JSTEP is refused as ``unsupported`` with the reason.
"""
from __future__ import annotations

import ctypes
import hashlib
import json
import os
import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

from residual_core.runtime import load_shared_library

from .sources import CorpusCase

__all__ = ["ProviderBuildFailed", "CorpusProvider", "adapt_source", "build_provider_for"]


class ProviderBuildFailed(RuntimeError):
    """The live OTI provider could not be produced for this source."""

    def __init__(self, failure_class: str, detail: str):
        self.failure_class = failure_class
        self.detail = detail
        super().__init__("%s: %s" % (failure_class, detail))


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


# ---------------------------------------------------------------------------
# declaration-only source adaptation
# ---------------------------------------------------------------------------
_JSTEP = re.compile(r"\bjstep\b", re.IGNORECASE)


def _is_comment(line: str, fixed: bool) -> bool:
    stripped = line.lstrip()
    if stripped.startswith("!"):
        return True
    return fixed and bool(line) and line[0] in "cC*!"


def _code_part(line: str, fixed: bool) -> str:
    if _is_comment(line, fixed):
        return ""
    return line.split("!", 1)[0]


def adapt_source(text: str, fixed: bool) -> tuple[str, List[str]]:
    """Return (source, notes). See the module docstring for the one edit."""
    lines = text.splitlines()
    hits = [i for i, line in enumerate(lines) if _JSTEP.search(_code_part(line, fixed))]
    if not hits:
        return text, []
    occurrences = sum(len(_JSTEP.findall(_code_part(lines[i], fixed))) for i in hits)
    declared = [i for i in hits if re.search(r"\bjstep\s*\(\s*4\s*\)", _code_part(lines[i], fixed),
                                             re.IGNORECASE)]
    if occurrences != 2 or len(declared) != 1:
        raise ProviderBuildFailed(
            "jstep_array_in_use",
            "the source names JSTEP %d times (expected 2: the argument list and one "
            "JSTEP(4) declaration); the provider wrappers pass a scalar KSTEP there and "
            "a routine that reads JSTEP cannot be adapted without changing it" % occurrences)
    decl = declared[0]
    header = [i for i in hits if i != decl][0]
    # The statement holding the declaration: walk back over continuation lines.
    start = decl
    while start > 0 and _continued(lines[start], fixed):
        start -= 1
    keyword = _code_part(lines[start], fixed).strip().lower()
    out = list(lines)
    out[header] = _JSTEP.sub("KSTEP", lines[header], count=1)
    if keyword.startswith("dimension") or re.match(r"^\d*\s*dimension", keyword):
        # a DIMENSION list: keep the (4) on a local nobody reads; KSTEP stays
        # implicitly INTEGER as the Abaqus include makes it
        out[decl] = re.sub(r"\bjstep\s*\(\s*4\s*\)", "JSTEPUNUSED(4)", lines[decl],
                           flags=re.IGNORECASE)
        how = "JSTEP(4) in a DIMENSION list renamed to the unused local JSTEPUNUSED(4)"
    else:
        out[decl] = re.sub(r"\bjstep\s*\(\s*4\s*\)", "kstep", lines[decl], flags=re.IGNORECASE)
        how = "JSTEP(4) in a typed declaration replaced by the scalar kstep"
    notes = [
        "declaration-only adaptation: dummy argument JSTEP renamed KSTEP at line %d; %s "
        "(line %d). JSTEP occurs nowhere else, so no executable statement changed. Reason: "
        "umat_oti.provider wrappers pass a scalar KSTEP in that slot and gfortran 9 refuses "
        "the rank mismatch against the lifted module's explicit interface."
        % (header + 1, how, decl + 1)]
    return "\n".join(out) + ("\n" if text.endswith("\n") else ""), notes


def _continued(line: str, fixed: bool) -> bool:
    """Is ``line`` a continuation of the previous one?"""
    if fixed:
        return len(line) > 5 and line[5] not in (" ", "0", "\t") and not _is_comment(line, True)
    return False


def _previous_continues(lines, index):  # pragma: no cover - free form helper
    return lines[index - 1].rstrip().endswith("&")


def source_reads_argument(path: Path, name: str, fixed: bool) -> bool:
    """Does any executable statement name the dummy ``name`` beyond the argument
    list and one declaration? (Conservative: any third occurrence counts.)

    Used for DROT: UMAT_OTI_EVAL_TOTAL_F passes DROT with NO derivative seed,
    so a routine that reads DROT loses d DROT/d u and d DROT/d p from every
    derivative -- its derivative claims are refused, not reported.
    """
    pattern = re.compile(r"\b%s\b" % re.escape(name), re.IGNORECASE)
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    count = sum(len(pattern.findall(_code_part(line, fixed))) for line in text.splitlines())
    return count > 2


# ---------------------------------------------------------------------------
# build
# ---------------------------------------------------------------------------
def _umat_oti_build_module(umat_repo: Optional[Path]):
    try:
        from umat_oti.provider import build as provider_build
    except ImportError as error:  # pragma: no cover - environment
        raise ProviderBuildFailed("umat_oti_unavailable", str(error))
    origin = Path(provider_build.__file__).resolve()
    if umat_repo is not None and Path(umat_repo).resolve() not in origin.parents:
        raise ProviderBuildFailed(
            "stale_umat_oti",
            "umat_oti resolves to %s, not the checkout %s; run with "
            "PYTHONPATH=%s/src" % (origin, umat_repo, umat_repo))
    return provider_build


def _umat_tree_state(module_file: str) -> Dict[str, Any]:
    """HEAD and a digest of the uncommitted diff of the umat_oti that built this.

    Other agents edit that tree while this runs; an object built from a dirty
    tree is reproducible only with the diff, so its digest is recorded.
    """
    repo = Path(module_file).resolve().parents[3]
    try:
        head = subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"], capture_output=True,
                              text=True, check=True).stdout.strip()
        diff = subprocess.run(["git", "-C", str(repo), "diff", "HEAD", "--", "src"],
                              capture_output=True, text=True, check=True).stdout
        files = subprocess.run(["git", "-C", str(repo), "diff", "HEAD", "--name-only", "--", "src"],
                               capture_output=True, text=True, check=True).stdout.split()
    except (OSError, subprocess.CalledProcessError) as error:
        return {"repo": str(repo), "error": str(error)}
    return {"repo": str(repo), "head": head, "dirty_files": files,
            "diff_sha256": _sha(diff.encode("utf-8")) if diff else None}


def build_provider_for(case: CorpusCase, out_dir: Path, *, umat_repo: Optional[Path] = None,
                       parameters: Optional[List[int]] = None) -> Dict[str, Any]:
    """Build (or reuse) the provider object for ``case`` under ``out_dir``.

    Returns the build record: paths, digests, adaptation notes. ``parameters``
    are 1-based PROPS slots to seed; default every slot.
    """
    out_dir = Path(out_dir)
    src_dir = out_dir / "src"
    build_dir = out_dir / "build"
    src_dir.mkdir(parents=True, exist_ok=True)
    original = case.source_path.read_bytes()
    text = original.decode("utf-8", errors="replace")
    fixed = case.source_form == "fixed"
    adapted, notes = adapt_source(text, fixed)
    name = case.source_path.name.replace(" ", "_")
    (src_dir / name).write_text(adapted, encoding="utf-8")
    slots = list(parameters or range(1, len(case.props) + 1))
    kinematics = "finite_strain" if case.finite else "small_strain"
    stem = "corpus_%s" % case.key[:12]
    contract = {
        "schema": "resasm_umat_transform_v2",
        "source": {"entry_point": "UMAT", "main_file": name},
        "kinematics": kinematics,
        "dimensions": {"ntens": case.ntens, "nprops": len(case.props), "nstatev": case.nstatv},
        "parameters": [{"name": "P%d" % slot, "props_index": slot} for slot in slots],
        "derivative": {"response": "STRESS", "export": "DSIGMA_DP"},
        "history": {"state": "STATEV", "export": "DSTATEV_DP", "propagate": True,
                    "path_dependent": True},
        "output": {"object": stem + "_oti.obj", "contract": stem + "_oti.json"},
    }
    contract_path = src_dir / "contract_v2.json"
    contract_text = json.dumps(contract, indent=2) + "\n"
    record_path = out_dir / "build_record.json"
    fingerprint = _sha((adapted + contract_text).encode("utf-8"))
    if record_path.is_file():
        try:
            previous = json.loads(record_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            previous = {}
        if (previous.get("fingerprint") == fingerprint
                and Path(previous.get("object", "")).is_file()
                and _sha(Path(previous["object"]).read_bytes()) == previous.get("object_sha256")):
            previous["reused"] = True
            return previous
    contract_path.write_text(contract_text, encoding="utf-8")
    module = _umat_oti_build_module(umat_repo)
    try:
        result = module.build_provider(contract_path, build_dir)
    except Exception as error:  # ProviderBuildError, CalledProcessError, ...
        message = str(error)
        failure = "provider_build_failed"
        undefined = re.findall(r"undefined reference to `([^']+)'", message)
        if "Rank mismatch" in message:
            failure = "provider_rank_mismatch"
        elif "NTENS=6" in message:
            failure = "provider_ntens_unsupported"
        elif "Cannot open module file" in message:
            failure = "provider_companion_module_missing"
        elif undefined:
            failure = "provider_unresolved_symbol:" + ",".join(sorted(set(undefined)))
        elif "includes/extra sources" in message:
            failure = "provider_include_unsupported"
        raise ProviderBuildFailed(failure, message[-4000:])
    obj = Path(result["object"])
    completed = json.loads(Path(result["contract"]).read_text(encoding="utf-8"))
    record = {
        "fingerprint": fingerprint,
        "object": str(obj),
        "object_sha256": _sha(obj.read_bytes()),
        "contract": result["contract"],
        "completed_contract": completed,
        "source_sha256_original": _sha(original),
        "source_sha256_built": _sha(adapted.encode("utf-8")),
        "source_adapted": bool(notes),
        "adaptation_notes": notes,
        "parameters": contract["parameters"],
        "kinematics": kinematics,
        "umat_oti_module": module.__file__,
        "umat_tree": _umat_tree_state(module.__file__),
        "reused": False,
    }
    record_path.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    return record


# ---------------------------------------------------------------------------
# the linked library
# ---------------------------------------------------------------------------
_SDVINI_SHIM = r"""
subroutine corpus_sdvini(n, nstatv, statev, coords, noel, npt) bind(C, name="corpus_sdvini")
  use iso_c_binding
  implicit none
  integer(c_int), value :: n, nstatv
  real(c_double), intent(inout) :: statev(nstatv, n)
  real(c_double), intent(in) :: coords(3, n)
  integer(c_int), intent(in) :: noel(n), npt(n)
  real(8) :: sv(nstatv), co(3)
  integer :: i, ns, ncrds, el, ip, layer, kspt
  do i = 1, n
    sv = statev(:, i); co = coords(:, i); ns = nstatv; ncrds = 3
    el = noel(i); ip = npt(i); layer = 1; kspt = 1
    call sdvini(sv, co, ns, ncrds, el, ip, layer, kspt)
    statev(:, i) = sv
  end do
end subroutine corpus_sdvini
"""

_SHIM = r"""
subroutine corpus_regular(n, nprops, ntens, nstatv, props, stress, statev, stran, dstran, &
    time2, dtime, temp, dtemp, coords, celent, noel, npt, kstep, kinc, dfgrd0, dfgrd1, drot, &
    ddsdde, pnewdt) bind(C, name="corpus_regular")
  use iso_c_binding
  implicit none
  integer(c_int), value :: n, nprops, ntens, nstatv, kstep, kinc
  real(c_double), intent(in) :: props(nprops), time2(2)
  real(c_double), value :: dtime, temp, dtemp
  real(c_double), intent(inout) :: stress(ntens, n), statev(nstatv, n)
  real(c_double), intent(in) :: stran(ntens, n), dstran(ntens, n), coords(3, n), celent(n)
  real(c_double), intent(in) :: dfgrd0(3, 3, n), dfgrd1(3, 3, n), drot(3, 3, n)
  integer(c_int), intent(in) :: noel(n), npt(n)
  real(c_double), intent(out) :: ddsdde(ntens, ntens, n), pnewdt(n)
  real(8) :: s(ntens), sv(nstatv), dd(ntens, ntens), pr(nprops)
  real(8) :: sse, spd, scd, rpl, ddsddt(ntens), drplde(ntens), drpldt
  real(8) :: predef(1), dpred(1), dr(3, 3), f0(3, 3), f1(3, 3), pn, cel
  real(8) :: st(ntens), dst(ntens), tm(2), co(3), dtm, tp, dtp
  character(len=80) :: cmname
  integer :: i, ndi, nshr, ns, np, nt, el, ip, layer, kspt, ks, ki
  cmname = 'MATERIAL'
  do i = 1, n
    s = stress(:, i); sv = statev(:, i); pr = props
    st = stran(:, i); dst = dstran(:, i); tm = time2; co = coords(:, i)
    dtm = dtime; tp = temp; dtp = dtemp; cel = celent(i)
    sse = 0.0d0; spd = 0.0d0; scd = 0.0d0; rpl = 0.0d0; drpldt = 0.0d0
    ddsddt = 0.0d0; drplde = 0.0d0; predef = 0.0d0; dpred = 0.0d0; dd = 0.0d0
    dr = drot(:, :, i); f0 = dfgrd0(:, :, i); f1 = dfgrd1(:, :, i)
    pn = 1.0d0; ndi = 3; nshr = ntens - 3; nt = ntens; ns = nstatv; np = nprops
    el = noel(i); ip = npt(i); layer = 1; kspt = 1; ks = kstep; ki = kinc
    call umat(s, sv, dd, sse, spd, scd, rpl, ddsddt, drplde, drpldt, st, dst, tm, &
        dtm, tp, dtp, predef, dpred, cmname, ndi, nshr, nt, ns, &
        pr, np, co, dr, pn, cel, f0, f1, el, ip, layer, kspt, ks, ki)
    stress(:, i) = s; statev(:, i) = sv; ddsdde(:, :, i) = dd; pnewdt(i) = pn
  end do
end subroutine corpus_regular

subroutine corpus_total(n, nprops, ntens, nstatv, nparam, finite, props, stress, statev, &
    stran, dstran, time2, dtime, temp, dtemp, dsig_in, dstv_in, stran_dp, dstran_dp, &
    coords, celent, noel, npt, kstep, kinc, dfgrd0, dfgrd1, drot, dfgrd0_dp, dfgrd1_dp, &
    ddsdde, dsig, dstv, dstv_de, pnewdt) bind(C, name="corpus_total")
  use iso_c_binding
  implicit none
  integer(c_int), value :: n, nprops, ntens, nstatv, nparam, finite, kstep, kinc
  real(c_double), intent(in) :: props(nprops), time2(2)
  real(c_double), value :: dtime, temp, dtemp
  real(c_double), intent(inout) :: stress(ntens, n), statev(nstatv, n)
  real(c_double), intent(in) :: stran(ntens, n), dstran(ntens, n), coords(3, n), celent(n)
  real(c_double), intent(in) :: dsig_in(ntens, nparam, n), dstv_in(nstatv, nparam, n)
  real(c_double), intent(in) :: stran_dp(ntens, nparam, n), dstran_dp(ntens, nparam, n)
  real(c_double), intent(in) :: dfgrd0(3, 3, n), dfgrd1(3, 3, n), drot(3, 3, n)
  real(c_double), intent(in) :: dfgrd0_dp(3, 3, nparam, n), dfgrd1_dp(3, 3, nparam, n)
  integer(c_int), intent(in) :: noel(n), npt(n)
  real(c_double), intent(out) :: ddsdde(ntens, ntens, n), dsig(ntens, nparam, n)
  real(c_double), intent(out) :: dstv(nstatv, nparam, n), dstv_de(nstatv, ntens, n), pnewdt(n)
  real(8) :: s(ntens), sv(nstatv), dd(ntens, ntens), ds(ntens, nparam), dv(nstatv, nparam)
  real(8) :: dsi(ntens, nparam), dvi(nstatv, nparam), e0(ntens, nparam), de(ntens, nparam)
  real(8) :: dve(nstatv, ntens), pn, pr(nprops), st(ntens), dst(ntens), tm(2), co(3), cel
  real(8) :: f0(3, 3), f1(3, 3), dr(3, 3), f0p(3, 3, nparam), f1p(3, 3, nparam)
  integer :: i, el, ip, ks, ki, np, nt, ns, npar
  do i = 1, n
    s = stress(:, i); sv = statev(:, i); pr = props
    st = stran(:, i); dst = dstran(:, i); tm = time2; co = coords(:, i); cel = celent(i)
    dsi = dsig_in(:, :, i); dvi = dstv_in(:, :, i)
    e0 = stran_dp(:, :, i); de = dstran_dp(:, :, i)
    f0 = dfgrd0(:, :, i); f1 = dfgrd1(:, :, i); dr = drot(:, :, i)
    f0p = dfgrd0_dp(:, :, :, i); f1p = dfgrd1_dp(:, :, :, i)
    el = noel(i); ip = npt(i); ks = kstep; ki = kinc; np = nprops; nt = ntens; ns = nstatv
    npar = nparam
    if (finite /= 0) then
      call umat_oti_eval_total_f(s, sv, dd, st, dst, tm, dtime, temp, dtemp, pr, np, nt, ns, &
          npar, ds, dv, dsi, dvi, e0, de, dve, co, cel, el, ip, ks, ki, pn, f0, f1, dr, f0p, f1p)
    else
      call umat_oti_eval_total(s, sv, dd, st, dst, tm, dtime, temp, dtemp, pr, np, nt, ns, &
          npar, ds, dv, dsi, dvi, e0, de, dve, co, cel, el, ip, ks, ki, pn)
    end if
    stress(:, i) = s; statev(:, i) = sv; ddsdde(:, :, i) = dd
    dsig(:, :, i) = ds; dstv(:, :, i) = dv; dstv_de(:, :, i) = dve; pnewdt(i) = pn
  end do
end subroutine corpus_total

subroutine corpus_flush() bind(C, name="corpus_flush")
  flush(6)
end subroutine corpus_flush
"""


def _ptr(array: np.ndarray):
    return array.ctypes.data_as(ctypes.c_void_p)


class MaterialCallError(RuntimeError):
    """A UMAT call returned something an increment cannot be built on."""


@dataclass
class CorpusProvider:
    """The linked provider: ``regular`` and ``total`` over a batch of points.

    Python-side layouts (C order, point axis first):
      stress (n, NTENS); state (n, NSTATV); F (n, 3, 3) with F[q, i, J];
      DDSDDE (n, NTENS, NTENS) with [q, a, b] = d sigma_a / d eps_b;
      dsig (n, NTENS, NPARAM); dstate (n, NSTATV, NPARAM);
      dF (n, NPARAM, 3, 3) with dF[q, j, i, J] = d F_iJ / d p_j.
    """

    record: Dict[str, Any]
    case: CorpusCase
    workdir: Path
    lib: Any = field(init=False, repr=False)
    nparam: int = field(init=False)
    slots: List[int] = field(init=False)

    def __post_init__(self):
        self.workdir = Path(self.workdir)
        self.workdir.mkdir(parents=True, exist_ok=True)
        self.slots = [int(p["props_index"]) for p in self.record["parameters"]]
        self.nparam = len(self.slots)
        self.ntens = int(self.case.ntens)
        self.nstatv = int(self.case.nstatv)
        self.nprops = len(self.case.props)
        self.finite = bool(self.case.finite)
        symbols = subprocess.run(["nm", "-g", self.record["object"]], capture_output=True,
                                 text=True).stdout
        self.has_sdvini = bool(re.search(r"\sT\s+sdvini_\s*$", symbols, re.MULTILINE))
        self.reads_drot = source_reads_argument(self.case.source_path, "drot",
                                                self.case.source_form == "fixed")
        library = self.workdir / ("corpus_provider_%s%s.so" % (self.record["object_sha256"][:16],
                                                               "_sdvini" if self.has_sdvini else ""))
        if not library.is_file():
            shim = self.workdir / "corpus_shim.f90"
            shim.write_text(_SHIM + (_SDVINI_SHIM if self.has_sdvini else ""), encoding="utf-8")
            run = dict(check=True, capture_output=True, text=True)
            try:
                subprocess.run(["gfortran", "-O2", "-fPIC", "-ffree-form", "-c", str(shim),
                                "-o", str(self.workdir / "corpus_shim.o")], **run)
                subprocess.run(["gfortran", "-shared", str(self.workdir / "corpus_shim.o"),
                                self.record["object"], "-o", str(library)], **run)
            except subprocess.CalledProcessError as error:
                raise ProviderBuildFailed("shim_link_failed", (error.stdout or "") + (error.stderr or ""))
        self._handle = load_shared_library(str(library))
        self.lib = self._handle.lib
        for symbol in ("corpus_regular", "corpus_total", "corpus_flush"):
            getattr(self.lib, symbol).restype = None
        self.library_path = str(library)

    def flush(self):
        self.lib.corpus_flush()

    def sdvini(self, state, coords, noel, npt):
        """The routine's own SDVINI (when the source defines one), as Abaqus
        calls it when the deck asks for user-subroutine initial state."""
        state = np.ascontiguousarray(state, dtype=float).copy()
        n = state.shape[0]
        self.lib.corpus_sdvini.restype = None
        self.lib.corpus_sdvini(ctypes.c_int(n), ctypes.c_int(self.nstatv), _ptr(state),
                               _ptr(np.ascontiguousarray(coords, dtype=float)),
                               _ptr(np.ascontiguousarray(noel, dtype=np.int32)),
                               _ptr(np.ascontiguousarray(npt, dtype=np.int32)))
        return state

    # ------------------------------------------------------------------
    def _points(self, stress, state, stran, dstran, F0, F1, coords, celent, noel, npt,
                drot=None):
        n = stress.shape[0]
        def c(a, shape):
            a = np.ascontiguousarray(a, dtype=float)
            if a.shape != shape:
                raise ValueError("expected %s, got %s" % (shape, a.shape))
            return a
        out = dict(
            stress=c(stress, (n, self.ntens)).copy(), state=c(state, (n, self.nstatv)).copy(),
            stran=c(stran, (n, self.ntens)), dstran=c(dstran, (n, self.ntens)),
            F0=np.ascontiguousarray(np.transpose(c(F0, (n, 3, 3)), (0, 2, 1))),
            F1=np.ascontiguousarray(np.transpose(c(F1, (n, 3, 3)), (0, 2, 1))),
            drot=np.ascontiguousarray(np.tile(np.eye(3), (n, 1, 1)) if drot is None else
                                      np.transpose(c(drot, (n, 3, 3)), (0, 2, 1))),
            coords=c(coords, (n, 3)), celent=c(celent, (n,)),
            noel=np.ascontiguousarray(noel, dtype=np.int32),
            npt=np.ascontiguousarray(npt, dtype=np.int32))
        return n, out

    def props_array(self, props):
        props = np.ascontiguousarray(props, dtype=float)
        if props.shape != (self.nprops,) or not np.all(np.isfinite(props)):
            raise ValueError("PROPS must be %d finite reals" % self.nprops)
        return props

    def regular(self, props, stress, state, stran, dstran, F0, F1, *, time, dtime, coords,
                celent, noel, npt, kstep=1, kinc=1, temp=0.0, dtemp=0.0, drot=None):
        props = self.props_array(props)
        n, a = self._points(stress, state, stran, dstran, F0, F1, coords, celent, noel, npt, drot)
        ddsdde = np.empty((n, self.ntens, self.ntens))
        pnewdt = np.empty(n)
        ci, cd = ctypes.c_int, ctypes.c_double
        self.lib.corpus_regular(
            ci(n), ci(self.nprops), ci(self.ntens), ci(self.nstatv), _ptr(props),
            _ptr(a["stress"]), _ptr(a["state"]), _ptr(a["stran"]), _ptr(a["dstran"]),
            _ptr(np.ascontiguousarray(time, dtype=float)), cd(dtime), cd(temp), cd(dtemp),
            _ptr(a["coords"]), _ptr(a["celent"]), _ptr(a["noel"]), _ptr(a["npt"]),
            ci(kstep), ci(kinc), _ptr(a["F0"]), _ptr(a["F1"]), _ptr(a["drot"]),
            _ptr(ddsdde), _ptr(pnewdt))
        self._check(pnewdt, a, "ORIGINAL")
        return {"stress": a["stress"], "state": a["state"],
                "ddsdde": np.ascontiguousarray(np.transpose(ddsdde, (0, 2, 1))),
                "pnewdt": pnewdt}

    def total(self, props, stress, state, stran, dstran, F0, F1, *, time, dtime, coords,
              celent, noel, npt, dstress_in=None, dstate_in=None, stran_dp=None,
              dstran_dp=None, dF0_dp=None, dF1_dp=None, kstep=1, kinc=1, temp=0.0,
              dtemp=0.0, drot=None):
        props = self.props_array(props)
        n, a = self._points(stress, state, stran, dstran, F0, F1, coords, celent, noel, npt, drot)
        nt, ns, npar = self.ntens, self.nstatv, self.nparam
        def seed(value, shape):
            if value is None:
                return np.zeros(shape)
            value = np.asarray(value, dtype=float)
            if value.shape != shape or not np.all(np.isfinite(value)):
                raise ValueError("seed must be finite with shape %s, got %s" % (shape, value.shape))
            return value
        # Fortran (ntens, nparam, n) == C (n, nparam, ntens)
        dsi = np.ascontiguousarray(np.transpose(seed(dstress_in, (n, nt, npar)), (0, 2, 1)))
        dvi = np.ascontiguousarray(np.transpose(seed(dstate_in, (n, ns, npar)), (0, 2, 1)))
        e0 = np.ascontiguousarray(np.transpose(seed(stran_dp, (n, nt, npar)), (0, 2, 1)))
        de = np.ascontiguousarray(np.transpose(seed(dstran_dp, (n, nt, npar)), (0, 2, 1)))
        # Fortran (3,3,nparam,n) == C (n, nparam, 3[col J], 3[row i])
        f0p = np.ascontiguousarray(np.transpose(seed(dF0_dp, (n, npar, 3, 3)), (0, 1, 3, 2)))
        f1p = np.ascontiguousarray(np.transpose(seed(dF1_dp, (n, npar, 3, 3)), (0, 1, 3, 2)))
        ddsdde = np.empty((n, nt, nt))
        dsig = np.empty((n, npar, nt))
        dstv = np.empty((n, npar, ns))
        dstv_de = np.empty((n, nt, ns))
        pnewdt = np.empty(n)
        ci, cd = ctypes.c_int, ctypes.c_double
        self.lib.corpus_total(
            ci(n), ci(self.nprops), ci(nt), ci(ns), ci(npar), ci(1 if self.finite else 0),
            _ptr(props), _ptr(a["stress"]), _ptr(a["state"]), _ptr(a["stran"]),
            _ptr(a["dstran"]), _ptr(np.ascontiguousarray(time, dtype=float)), cd(dtime),
            cd(temp), cd(dtemp), _ptr(dsi), _ptr(dvi), _ptr(e0), _ptr(de), _ptr(a["coords"]),
            _ptr(a["celent"]), _ptr(a["noel"]), _ptr(a["npt"]), ci(kstep), ci(kinc),
            _ptr(a["F0"]), _ptr(a["F1"]), _ptr(a["drot"]), _ptr(f0p), _ptr(f1p),
            _ptr(ddsdde), _ptr(dsig), _ptr(dstv), _ptr(dstv_de), _ptr(pnewdt))
        self._check(pnewdt, a, "OTI")
        return {"stress": a["stress"], "state": a["state"],
                "ddsdde": np.ascontiguousarray(np.transpose(ddsdde, (0, 2, 1))),
                "dstress_dp": np.ascontiguousarray(np.transpose(dsig, (0, 2, 1))),
                "dstate_dp": np.ascontiguousarray(np.transpose(dstv, (0, 2, 1))),
                "dstate_ddstran": np.ascontiguousarray(np.transpose(dstv_de, (0, 2, 1))),
                "pnewdt": pnewdt}

    def _check(self, pnewdt, a, which):
        bad = (~np.isfinite(pnewdt)) | (pnewdt < 1.0)
        if np.any(bad):
            raise MaterialCallError("%s UMAT requested a time cut-back (PNEWDT<1) at %d points"
                                    % (which, int(bad.sum())))
        # STATEV entries the verification run found undefined in the original
        # (sources.load_case) may hold anything, NaN included; whether the
        # routine reads them is the hidden-state probe's question, not this one's
        undefined = [i - 1 for i in (getattr(self.case, "extra", None) or {}).get(
            "undefined_statev", ()) if 1 <= i <= a["state"].shape[-1]]
        for name in ("stress", "state"):
            value = a[name]
            if name == "state" and undefined:
                value = np.delete(value, undefined, axis=-1)
            if not np.all(np.isfinite(value)):
                raise MaterialCallError("%s UMAT returned a non-finite %s" % (which, name))
