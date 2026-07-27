"""Generic Abaqus-UMAT material backend.

Presents an arbitrary Abaqus UMAT (``STRESS, DDSDDE, STATEV, DFGRD0/1, PROPS,
...``) to the formulation layer as a :class:`~residual_core.materials.base.Material`.
The UMAT contract *is* the Material contract (see materials/base.py), so this
class is a thin bridge that:

    * receives ``kinematics = {'F0':F0, 'F1':F1, 'element':eid, 'ip':ip}`` from a
      finite-strain formulation (kinematic_input='deformation_gradient'),
    * carries the previous per-IP STATEV slab in ``state_prev``,
    * carries PROPS in ``binding.constants``,
    * calls a UMAT through one of two configurable backends, and
    * returns ``(cauchy_stress(6,), ddsdde(6,6), statev_new(nstatv,), diagnostics)``.

Two backends
------------
``backend='python'`` (offline-testable)
    Call an *injected* Python callable that has the UMAT signature. NO subprocess
    is spawned. Used for mock/testing UMATs so the whole material-update ->
    residual -> tangent pipeline is exercisable without ifort/Abaqus. The
    callable is invoked as ``umat_fn(**umat_kwargs)`` with the canonical Abaqus
    argument names (lowercase) and must return ``(stress, ddsdde, statev,
    pnewdt)`` (a shorter tuple, or ``None`` for in-place mutation of the passed
    arrays, is also accepted). See ``mock_isotropic_elastic_umat`` below for a
    reference implementation.

``backend='fortran'`` (real UMAT)
    Shell out, ONE increment per call, to the compiled ``umat_driver`` binary
    from ``residual_core/umat_adapter_fortran/`` -- reusing that package's
    verified input/output plumbing (``umat_replay.py``: ``write_driver_input``,
    ``find_driver``, ``run_driver``). STATEV is carried forward across calls via
    ``state_prev`` -> the driver's ``STATEV0``. If no suitable binary is built,
    a clear, actionable error naming the ifort+Abaqus requirement is raised --
    the adapter never fabricates a stress.

    HONEST CAVEAT: the real Grilli UMAT keeps *additional* history in a static
    ``/UMPS/`` COMMON block (kFp, kcurlFp, kSigma0, ...) that persists between
    increments only within ONE process. A per-call subprocess re-zeros that
    COMMON block every call, so this backend faithfully carries STATEV but NOT
    the COMMON-block history. For an exact real-UMAT single-IP history replay,
    drive the whole loading path through the one-process
    ``umat_adapter_fortran/umat_replay.py`` instead.
"""

from __future__ import annotations

import importlib.util
import os
import tempfile
from typing import Any, Callable, Dict, Optional

import numpy as np

from .base import Material, MaterialBinding


# ---------------------------------------------------------------------------
# Reuse the verified Fortran-driver plumbing WITHOUT importing it as a package
# (umat_adapter_fortran has no __init__.py) and WITHOUT modifying it.
# ---------------------------------------------------------------------------
_UMAT_REPLAY = None


def _load_umat_replay():
    """Import ``umat_adapter_fortran/umat_replay.py`` by file path (cached).

    Reuses its driver helpers (find_driver / write_driver_input / run_driver)
    so the fortran backend shares the exact byte-for-byte driver I/O contract.
    Loading it runs only module-level code (path setup + a c3d8_kernel import);
    its argparse ``main()`` is guarded by ``__name__ == '__main__'``.
    """
    global _UMAT_REPLAY
    if _UMAT_REPLAY is not None:
        return _UMAT_REPLAY
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.normpath(
        os.path.join(here, "..", "umat_adapter_fortran", "umat_replay.py"))
    if not os.path.exists(path):
        raise RuntimeError(
            "cannot locate the Fortran driver bridge umat_replay.py "
            "(expected at %s)" % path)
    spec = importlib.util.spec_from_file_location("umat_replay_bridge", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    _UMAT_REPLAY = mod
    return mod


def _time2(t) -> np.ndarray:
    """Coerce whatever the formulation passes as ``time`` into TIME(2)
    = [step time, total time]."""
    if t is None:
        return np.array([0.0, 0.0])
    a = np.atleast_1d(np.asarray(t, dtype=float)).ravel()
    if a.size == 0:
        return np.array([0.0, 0.0])
    if a.size == 1:
        return np.array([a[0], a[0]])
    return np.array([a[0], a[1]])


class UmatAdapter(Material):
    """Generic Abaqus UMAT presented as a Material.

    Parameters
    ----------
    n_state_vars : int
        NSTATV. May be overridden per element by ``binding.n_state_vars``.
    backend : {'fortran', 'python'}
        'python' calls ``umat_fn`` in-process (offline-testable); 'fortran'
        shells out to the compiled ``umat_driver``.
    umat_fn : callable, optional
        Required for backend='python'. Called ``umat_fn(**umat_kwargs)`` with
        the canonical Abaqus UMAT argument names; returns ``(stress, ddsdde,
        statev, pnewdt)`` (shorter tuple or None -> in-place also accepted).
    driver : str, optional
        Explicit path to a ``umat_driver`` binary (backend='fortran').
    cmname : str
        CMNAME passed to the UMAT.
    default_temp : float
        Temperature used when none is supplied via fields/options.
    rotation_prop_offset, rotation_statev_offset, fp_statev_offset : int or None
        Drive :meth:`init_state`'s generic first-increment seeding (rotation
        matrix from PROPS into STATEV, plastic Fp=I into STATEV). Set any to
        ``None`` to disable that piece. Defaults follow the common Abaqus
        crystal-plasticity convention (PROPS[1:10] rotation, STATEV[0:9]
        rotation, STATEV[80:89] Fp). Override for a different UMAT.
    allow_mock_driver : bool
        If False, backend='fortran' refuses a ``*_mock`` driver binary (it is
        not the intended real UMAT) and raises the ifort+Abaqus error instead.
    prefer_mock : bool
        Passed to the driver-discovery order.
    keep_files : bool
        Keep the per-call driver scratch files (backend='fortran') for
        inspection instead of deleting them.
    """

    name = "umat"
    stress_measure = "cauchy"
    tangent_measure = "ddsdde"
    kinematic_input = "deformation_gradient"

    def __init__(self,
                 n_state_vars: int = 0,
                 backend: str = "fortran",
                 umat_fn: Optional[Callable[..., Any]] = None,
                 driver: Optional[str] = None,
                 cmname: str = "UMATMATERIAL",
                 default_temp: float = 293.0,
                 rotation_prop_offset: Optional[int] = 1,
                 rotation_statev_offset: Optional[int] = 0,
                 fp_statev_offset: Optional[int] = 80,
                 allow_mock_driver: bool = True,
                 prefer_mock: bool = False,
                 keep_files: bool = False,
                 name: Optional[str] = None):
        if backend not in ("fortran", "python"):
            raise ValueError("backend must be 'fortran' or 'python', got %r"
                             % (backend,))
        self.n_state_vars = int(n_state_vars)
        self.backend = backend
        self.umat_fn = umat_fn
        self.driver = driver
        self.cmname = str(cmname)
        self.default_temp = float(default_temp)
        self.rotation_prop_offset = rotation_prop_offset
        self.rotation_statev_offset = rotation_statev_offset
        self.fp_statev_offset = fp_statev_offset
        self.allow_mock_driver = bool(allow_mock_driver)
        self.prefer_mock = bool(prefer_mock)
        self.keep_files = bool(keep_files)
        if name is not None:
            self.name = str(name)

    # ------------------------------------------------------------------
    # Material contract
    # ------------------------------------------------------------------
    def evaluate(self, kinematics: Dict[str, Any], state_prev: np.ndarray,
                 binding: MaterialBinding, time, dtime: float,
                 fields: Optional[dict], options: Optional[dict]):
        options = options or {}
        fields = fields or {}

        F0 = np.asarray(kinematics["F0"], dtype=float).reshape(3, 3)
        F1 = np.asarray(kinematics["F1"], dtype=float).reshape(3, 3)
        noel = int(kinematics.get("element", 1) or 1)
        npt = int(kinematics.get("ip", 1) or 1)

        props = np.asarray(list(binding.constants), dtype=float) \
            if (binding is not None and binding.constants) else np.zeros(0)
        nprops = int(props.size)

        nstatv = int((binding.n_state_vars if binding is not None else 0)
                     or self.n_state_vars)
        statev = self._normalize_statev(state_prev, nstatv)
        nstatv = int(statev.size)

        ntens = 6
        time_arr = _time2(time)
        dtime = float(dtime) if dtime else 0.0
        temp = float(fields.get("temp", options.get("temp", self.default_temp)))
        dtemp = float(fields.get("dtemp", options.get("dtemp", 0.0)))
        kstep = int(options.get("kstep", 1))
        kinc = options.get("kinc")
        if kinc is None:
            # First increment of the first step fires the UMAT init block; a
            # zero step-time at increment start is the offline heuristic for it.
            kinc = 1 if time_arr[0] <= 0.0 else 2
        kinc = int(kinc)
        coords = np.asarray(kinematics.get("coords", np.zeros(3)),
                            dtype=float).ravel()[:3]
        if coords.size < 3:
            coords = np.zeros(3)
        celent = float(options.get("celent", 1.0))

        umat_kwargs = dict(
            stress=np.zeros(ntens),
            statev=statev.copy(),
            ddsdde=np.zeros((ntens, ntens)),
            sse=0.0, spd=0.0, scd=0.0,
            rpl=0.0, ddsddt=np.zeros(ntens), drplde=np.zeros(ntens), drpldt=0.0,
            stran=np.zeros(ntens), dstran=np.zeros(ntens),
            time=time_arr.copy(), dtime=dtime, temp=temp, dtemp=dtemp,
            predef=np.zeros(1), dpred=np.zeros(1), cmname=self.cmname,
            ndi=3, nshr=3, ntens=ntens, nstatv=nstatv,
            props=props.copy(), nprops=nprops,
            coords=coords.copy(), drot=np.eye(3), pnewdt=1.0, celent=celent,
            dfgrd0=F0.copy(), dfgrd1=F1.copy(),
            noel=noel, npt=npt, layer=1, kspt=1, kstep=kstep, kinc=kinc,
        )

        if self.backend == "python":
            stress, ddsdde, statev_new, pnewdt, extra = self._eval_python(umat_kwargs)
        else:
            stress, ddsdde, statev_new, pnewdt, extra = self._eval_fortran(umat_kwargs)

        stress = np.asarray(stress, dtype=float).ravel()[:ntens]
        ddsdde = np.asarray(ddsdde, dtype=float).reshape(ntens, ntens)
        statev_new = np.asarray(statev_new, dtype=float).ravel()

        diag = {"material": self.name, "backend": self.backend,
                "pnewdt": float(pnewdt), "kinc": kinc, "kstep": kstep,
                "element": noel, "ip": npt}
        diag.update(extra)
        return stress, ddsdde, statev_new, diag

    # ------------------------------------------------------------------
    # backend: python (in-process, no subprocess)
    # ------------------------------------------------------------------
    def _eval_python(self, kw: Dict[str, Any]):
        if self.umat_fn is None:
            raise ValueError(
                "UmatAdapter(backend='python') requires an injected `umat_fn` "
                "(a callable with the Abaqus UMAT signature). None was given.")
        # references we can read back if the callable mutates in place
        stress_ref = kw["stress"]
        ddsdde_ref = kw["ddsdde"]
        statev_ref = kw["statev"]
        ret = self.umat_fn(**kw)
        pnewdt = float(kw.get("pnewdt", 1.0))
        if ret is None:
            return stress_ref, ddsdde_ref, statev_ref, pnewdt, {"binary": None}
        if not isinstance(ret, (tuple, list)):
            raise TypeError("python umat_fn must return a tuple "
                            "(stress, ddsdde, statev, pnewdt) or None; got %r"
                            % type(ret))
        stress = ret[0] if len(ret) >= 1 else stress_ref
        ddsdde = ret[1] if len(ret) >= 2 else ddsdde_ref
        statev = ret[2] if len(ret) >= 3 else statev_ref
        if len(ret) >= 4 and ret[3] is not None:
            pnewdt = float(ret[3])
        return stress, ddsdde, statev, pnewdt, {"binary": None}

    # ------------------------------------------------------------------
    # backend: fortran (one increment per subprocess call)
    # ------------------------------------------------------------------
    def _eval_fortran(self, kw: Dict[str, Any]):
        R = _load_umat_replay()
        driver = R.find_driver(self.driver, prefer_mock=self.prefer_mock)
        is_mock = bool(driver) and ("mock" in os.path.basename(driver))
        if driver is None or (is_mock and not self.allow_mock_driver):
            raise RuntimeError(self._no_driver_message(found_mock=is_mock))

        ntens = int(kw["ntens"])
        nstatv = int(kw["nstatv"])
        props = list(np.asarray(kw["props"], dtype=float).ravel())
        statev0 = list(np.asarray(kw["statev"], dtype=float).ravel())
        time_arr = np.asarray(kw["time"], dtype=float).ravel()
        inc = dict(noel=int(kw["noel"]), npt=int(kw["npt"]),
                   kstep=int(kw["kstep"]), kinc=int(kw["kinc"]),
                   dtime=float(kw["dtime"]), temp=float(kw["temp"]),
                   dtemp=float(kw["dtemp"]),
                   time=(float(time_arr[0]), float(time_arr[1])),
                   dfgrd0=np.asarray(kw["dfgrd0"], dtype=float),
                   dfgrd1=np.asarray(kw["dfgrd1"], dtype=float))

        workdir = tempfile.mkdtemp(prefix="umat_adapter_")
        infile = os.path.join(workdir, "inc_in.txt")
        outfile = os.path.join(workdir, "inc_out.txt")
        try:
            R.write_driver_input(infile, nstatv, len(props), ntens,
                                 props, statev0, [inc])
            recs = R.run_driver(driver, infile, outfile)
            rec = recs[0]
            stress = np.asarray(rec["stress"], dtype=float)
            statev = np.asarray(rec["statev"], dtype=float)
            ddsdde = np.asarray(rec["ddsdde"], dtype=float).reshape(ntens, ntens)
            pnewdt = float(rec["pnewdt"])
        finally:
            if not self.keep_files:
                for p in (infile, outfile):
                    try:
                        os.remove(p)
                    except OSError:
                        pass
                try:
                    os.rmdir(workdir)
                except OSError:
                    pass
        return stress, ddsdde, statev, pnewdt, {
            "binary": driver, "is_mock_binary": is_mock}

    def _no_driver_message(self, found_mock: bool) -> str:
        here = os.path.dirname(os.path.abspath(__file__))
        build = os.path.normpath(
            os.path.join(here, "..", "umat_adapter_fortran", "build"))
        lead = ("Only the isotropic-elastic MOCK driver binary was found, which "
                "is NOT the intended real UMAT."
                if found_mock else
                "No compiled 'umat_driver' binary was found.")
        return (
            "UmatAdapter(backend='fortran') could not run: %s\n"
            "Expected a real driver at: %s/umat_driver[.exe].\n"
            "The real crystal-plasticity UMAT object only compiles with Intel "
            "ifort + Abaqus (MKL); gfortran cannot build the unmodified sources "
            "(two ifort-specific constructs -- see "
            "residual_core/umat_adapter_fortran/README.md 'Build status'). "
            "Build it inside the Abaqus environment "
            "(e.g. `abaqus make library=umat.for`, then link umat_driver.o), "
            "pass driver=<path> to a real binary, or use backend='python' with "
            "an injected umat_fn for offline testing. "
            "The adapter does NOT fabricate a stress." % (lead, build))

    # ------------------------------------------------------------------
    # helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _normalize_statev(state_prev, nstatv: int) -> np.ndarray:
        if state_prev is None:
            return np.zeros(int(nstatv), dtype=float)
        s = np.asarray(state_prev, dtype=float).ravel()
        if not nstatv:
            return s
        if s.size == nstatv:
            return s.copy()
        out = np.zeros(int(nstatv), dtype=float)
        m = min(int(nstatv), s.size)
        out[:m] = s[:m]
        return out

    def init_state(self, binding: MaterialBinding, coords=None,
                   n_ip: int = 8) -> np.ndarray:
        """Generic first-increment STATEV seeding, mirroring a UMAT's
        ``kinc<=1`` init block where knowable from PROPS:

          * rotation matrix PROPS[rotation_prop_offset : +9] -> STATEV[
            rotation_statev_offset : +9],
          * plastic deformation gradient Fp = I -> STATEV diagonal at
            fp_statev_offset (+0, +4, +8).

        Any piece whose offset is ``None`` (or that does not fit) is skipped.
        Override for UMAT-specific extras (see CrystalPlasticityAdapter).
        """
        n = int((binding.n_state_vars if binding is not None else 0)
                or self.n_state_vars)
        state = np.zeros((int(n_ip), n), dtype=float)
        if n == 0:
            return state
        c = np.asarray(list(binding.constants), dtype=float) \
            if (binding is not None and binding.constants) else np.zeros(0)

        ro, so = self.rotation_prop_offset, self.rotation_statev_offset
        if (ro is not None and so is not None
                and c.size >= ro + 9 and n >= so + 9):
            state[:, so:so + 9] = c[ro:ro + 9]

        fo = self.fp_statev_offset
        if fo is not None and n >= fo + 9:
            for d in (0, 4, 8):
                state[:, fo + d] = 1.0
        return state


# ---------------------------------------------------------------------------
# Reference / testing mock UMAT (python backend)
# ---------------------------------------------------------------------------
def mock_isotropic_elastic_umat(*, dfgrd1, ddsdde, statev, ntens, props,
                                stress, **kwargs):
    """A trivial, offline UMAT with the Abaqus signature, for testing wiring.

    Small-strain isotropic linear elasticity read from the total deformation
    gradient (NOT crystal plasticity):

        eps = sym(F1) - I          (engineering shear in Voigt)
        sigma = D(E, nu) : eps,    DDSDDE = D

    PROPS convention here: ``[E, nu]`` (defaults 100000, 0.3). Identity
    deformation -> zero stress. STATEV is passed through unchanged except that
    the Cauchy stress is mirrored into slots 48:53 (0-based 47:53) when present,
    matching the real-UMAT convention. Returns ``(stress, ddsdde, statev,
    pnewdt)``.
    """
    from ..core.voigt import isotropic_D   # neutral helper, not a formulation import

    p = np.asarray(props, dtype=float).ravel()
    E = float(p[0]) if p.size >= 1 else 100000.0
    nu = float(p[1]) if p.size >= 2 else 0.3
    D = isotropic_D(E, nu)

    F1 = np.asarray(dfgrd1, dtype=float).reshape(3, 3)
    eps = np.array([
        F1[0, 0] - 1.0,
        F1[1, 1] - 1.0,
        F1[2, 2] - 1.0,
        F1[0, 1] + F1[1, 0],          # engineering shear gamma12 = 2*e12
        F1[0, 2] + F1[2, 0],          # gamma13
        F1[1, 2] + F1[2, 1],          # gamma23
    ], dtype=float)
    sigma = D @ eps

    sv = np.asarray(statev, dtype=float).ravel().copy()
    if sv.size >= 53:
        sv[47:53] = sigma
    out_stress = np.asarray(stress, dtype=float).ravel().copy()
    out_stress[:ntens] = sigma[:ntens]
    out_ddsdde = np.asarray(ddsdde, dtype=float).reshape(ntens, ntens).copy()
    out_ddsdde[:, :] = D[:ntens, :ntens]
    return out_stress, out_ddsdde, sv, 1.0


# ---------------------------------------------------------------------------
# Offline self-test: prove the Material<->Formulation contract composes for the
# finite-strain / UMAT path using the python backend (no Abaqus needed).
#   run:  python -m residual_core.materials.umat_adapter
# ---------------------------------------------------------------------------
def _selftest() -> bool:
    from ..formulations.solid_c3d8_finite_strain import SolidC3D8FiniteStrain
    from ..formulations.c3d8_kernel import COMPRESSION111_ELEM1_XE

    ok = True
    E, nu = 100000.0, 0.3
    nstatv = 60  # exercise state carrying (>53 so slots 48:53 are written)
    mat = UmatAdapter(n_state_vars=nstatv, backend="python",
                      umat_fn=mock_isotropic_elastic_umat)
    binding = MaterialBinding(material=mat, constants=[E, nu],
                              n_state_vars=nstatv, name="mock_umat")
    form = SolidC3D8FiniteStrain()
    Xe = COMPRESSION111_ELEM1_XE.copy()
    mstate = mat.init_state(binding, n_ip=8)

    # --- (a) F = I (zero displacement) -> stress/residual ~ 0 --------------
    u0 = np.zeros(24)
    r0, K0, s0, d0 = form.eval_element(
        1, "C3D8", Xe, u0, {"dofs_prev": np.zeros(24)}, mstate, binding,
        time=(0.0, 0.0), dtime=1.0, fields=None,
        options={"compute_tangent": True})
    r0max = float(np.abs(r0).max())
    zero_ok = (K0 is not None and K0.shape == (24, 24)
               and np.all(np.isfinite(r0)) and np.all(np.isfinite(K0))
               and r0max < 1e-8)
    ok = ok and zero_ok
    print("  [%s] F=I: |r|max=%.3e (expect ~0), K shape=%s"
          % ("PASS" if zero_ok else "FAIL", r0max,
             None if K0 is None else K0.shape))

    # --- (b) small stretch -> finite residual + (24,24) tangent -----------
    u = np.zeros((8, 3))
    u[:, 2] = 1e-3 * Xe[:, 2]          # ~0.1% uniaxial stretch along z
    u = u.reshape(-1)
    r, K, s_new, diag = form.eval_element(
        2, "C3D8", Xe, u, {"dofs_prev": np.zeros(24)}, mstate, binding,
        time=(0.0, 0.0), dtime=1.0, fields=None,
        options={"compute_tangent": True})
    finite_ok = (r.shape == (24,) and np.all(np.isfinite(r))
                 and K is not None and K.shape == (24, 24)
                 and np.all(np.isfinite(K)) and float(np.abs(r).max()) > 0.0
                 and s_new is not None and s_new.shape == (8, nstatv))
    ok = ok and finite_ok
    # spot-check the mirrored Cauchy stress in STATEV[47:53]
    sigma_ip0 = s_new[0, 47:53] if s_new is not None else np.zeros(6)
    print("  [%s] stretch: |r|max=%.3e (finite), K shape=%s, state shape=%s"
          % ("PASS" if finite_ok else "FAIL", float(np.abs(r).max()),
             None if K is None else K.shape,
             None if s_new is None else s_new.shape))
    print("        STATEV[48:53] (Cauchy) IP1 = %s"
          % np.array2string(sigma_ip0, precision=3))

    print("  OVERALL: %s" % ("PASS" if ok else "FAIL"))
    return ok


if __name__ == "__main__":
    import sys
    print("UmatAdapter offline self-test (python backend + finite-strain C3D8)")
    sys.exit(0 if _selftest() else 1)
