"""OTILib adapter — the production hypercomplex backend seam.

OTILib (Order Truncated Imaginary numbers; Aristizabal et al.) is the real HYPAD
algebra for **arbitrary-order, multi-parameter** residual sensitivity. Its Python
bindings are the package ``pyoti`` built from source at
https://github.com/mauriaristi/otilib.git (branch ``master``, **GPLv3**). This
module hides all OTILib-specific naming behind a clean internal API
(``OtiContext``) so the rest of the framework never imports ``pyoti`` directly and
no GPLv3 code is vendored into the core.

Dual1 (``algebra/dual1.py``) was only a first-order smoke test — OTILib is the
production backend. See ``docs/otilib_integration.md`` and
``docs/otilib_api_inventory.md``.

Availability / detection
------------------------
The genuine OTI library must be installed separately (conda + CMake source build;
Windows requires WSL). It is **NOT** the unrelated PyPI package also named
``pyoti`` — do not ``pip install pyoti``. Detection tries, in order:

  A. explicit environment variables ``OTILIB_ROOT`` / ``PYOTI_PATH`` (python
     package roots) and records ``OTILIB_INCLUDE_DIR`` / ``OTILIB_LIBRARY_DIR``;
  B. importing the installed source-build package ``pyoti.sparse``;
  C. known local external checkouts ``external/gpl/otilib`` or
     ``third_party/gpl/otilib`` (clearly marked GPLv3).

A module is accepted only if it exposes the genuine OTI surface (``e`` +
``number`` + ``get_im``), which rejects the PyPI squat. When absent,
``otilib_available()`` returns False and callers skip cleanly — there is **no**
silent fall back to Dual1.

Canonical direction representation
----------------------------------
Directions are exponent multi-indices ``kappa`` (length = number of bases). For
``m=3``: order 1 -> [1,0,0],[0,1,0],[0,0,1]; order 2 -> [2,0,0],[1,1,0],[1,0,1],
[0,2,0],[0,1,1],[0,0,2]. Each direction has a derivative *recovery factor*
``prod_i factorial(kappa_i)``. OTILib addresses the same direction by the flat,
1-based basis-index list with multiplicity (``kappa=(2,0,1)`` -> ``[1,1,3]``),
which is what ``oti.e(...)`` / ``get_im(...)`` accept. This matches
``core/sensitivity_package.py``.
"""

from __future__ import annotations

import math
import os
import sys
from itertools import combinations_with_replacement
from typing import Any, Dict, List, Optional, Sequence, Tuple

# --------------------------------------------------------------------------- #
# backend detection: import the GENUINE OTI library, verify its API surface
# --------------------------------------------------------------------------- #
_OTI = None
_OTI_ERROR = ""
_OTI_API = ""

_ENV_ROOT_VARS = ("OTILIB_ROOT", "PYOTI_PATH")
_ENV_HINT_VARS = ("OTILIB_INCLUDE_DIR", "OTILIB_LIBRARY_DIR")
_LOCAL_ROOTS = ("external/gpl/otilib", "third_party/gpl/otilib")

_MISSING_MSG = (
    "OTILib backend requested but genuine OTILib was not found.\n"
    "Do not install PyPI pyoti (an UNRELATED library).\n"
    "Install from https://github.com/mauriaristi/otilib.git (branch master, "
    "GPLv3): build with conda + CMake (Windows requires WSL) — see "
    "docs/otilib_integration.md, or run scripts/setup_otilib.sh. You may also set "
    "OTILIB_ROOT or PYOTI_PATH to an existing build.")


def _repo_root() -> str:
    # residual_core/algebra/otilib_adapter.py -> repo root is three levels up.
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.dirname(os.path.dirname(here))


def _candidate_python_dirs() -> List[str]:
    """Directories that may contain an importable ``pyoti`` package: env-var
    roots (A) and local GPLv3 checkouts (C). Each root is probed at its top level
    and at the ``src/python`` / ``build`` subdirectories used by the OTILib build."""
    roots: List[str] = []
    for var in _ENV_ROOT_VARS:
        v = os.environ.get(var)
        if v:
            roots.append(v)
    base = _repo_root()
    for rel in _LOCAL_ROOTS:
        roots.append(os.path.join(base, *rel.split("/")))
    dirs: List[str] = []
    for r in roots:
        for sub in ("", os.path.join("src", "python"), "build"):
            d = os.path.join(r, sub) if sub else r
            if d and os.path.isdir(d):
                dirs.append(d)
    return dirs


def _accept(mod) -> bool:
    """Genuine OTI signature (rejects the unrelated PyPI squat)."""
    return (hasattr(mod, "e") and hasattr(mod, "number")
            and hasattr(mod, "get_im"))


def _probe_backend():
    """Import the real OTI library and confirm it exposes the OTI API (not the
    unrelated PyPI package of the same name). Sets module globals."""
    global _OTI, _OTI_ERROR, _OTI_API
    if _OTI is not None or _OTI_ERROR:
        return
    # (A + C) make env/local source-build python dirs importable.
    for d in _candidate_python_dirs():
        if d not in sys.path:
            sys.path.insert(0, d)
    # (B) import the genuine module; prefer the flexible dynamic-sparse impl.
    for modname in ("pyoti.sparse", "pyoti.static", "pyoti"):
        try:
            mod = __import__(modname, fromlist=["*"])
        except Exception:
            continue
        if modname == "pyoti" and not _accept(mod):
            mod = getattr(mod, "sparse", mod)
        if _accept(mod):
            _OTI = mod
            _OTI_API = getattr(mod, "__name__", modname)
            return
    _OTI_ERROR = _MISSING_MSG


def otilib_available() -> bool:
    _probe_backend()
    return _OTI is not None


def otilib_status() -> Dict[str, Any]:
    _probe_backend()
    hints = {v: os.environ.get(v) for v in (_ENV_ROOT_VARS + _ENV_HINT_VARS)
             if os.environ.get(v)}
    return {"available": _OTI is not None, "api_module": _OTI_API,
            "repo": "https://github.com/mauriaristi/otilib.git",
            "license": "GPLv3", "env": hints, "error": _OTI_ERROR}


# --------------------------------------------------------------------------- #
# direction enumeration (backend-independent; canonical exponent multi-indices)
# --------------------------------------------------------------------------- #
def num_coefficients_total(m: int, nt: int) -> int:
    """N = C(m + nt, m)."""
    return 1 if m == 0 else math.comb(m + nt, m)


def num_coefficients_order(m: int, p: int) -> int:
    """N^(p) = C(p + m - 1, p)."""
    if m == 0:
        return 1 if p == 0 else 0
    return math.comb(p + m - 1, p)


def order_directions(m: int, p: int) -> List[Dict[str, Any]]:
    """Ordered order-p directions: each {exponents, label, recovery_factor}."""
    out: List[Dict[str, Any]] = []
    if m <= 0 or p <= 0:
        return out
    for combo in combinations_with_replacement(range(1, m + 1), p):
        exps = [0] * m
        for b in combo:
            exps[b - 1] += 1
        parts = []
        for i, k in enumerate(exps, start=1):
            if k == 1:
                parts.append("e%d" % i)
            elif k > 1:
                parts.append("e%d^%d" % (i, k))
        factor = 1
        for k in exps:
            factor *= math.factorial(k)
        out.append({"exponents": tuple(exps), "label": "*".join(parts),
                    "recovery_factor": factor})
    return out


# --------------------------------------------------------------------------- #
# OtiContext — clean internal API over OTILib
# --------------------------------------------------------------------------- #
class OtiContext:
    """A hypercomplex working context with ``m`` bases and truncation order
    ``nt``. Wraps OTILib and exposes a framework-neutral API.

    All methods raise ``OtiUnavailableError`` if OTILib is not installed; callers
    should gate on ``otilib_available()`` and skip cleanly.
    """

    def __init__(self, num_bases: int, order: int):
        if not otilib_available():
            raise OtiUnavailableError(_OTI_ERROR)
        self.m = int(num_bases)
        self.nt = int(order)
        self._oti = _OTI

    # ---- scalar construction / seeding ----------------------------------
    def scalar(self, real_value: float):
        """A pure-real OTI scalar (all imaginary coefficients zero), carrying the
        context truncation order ``nt``."""
        return self._oti.number(float(real_value), order=self.nt)

    def seed(self, real_value: float, basis_index: int):
        """``real_value + e_{basis_index}`` — perturb one parameter along its own
        imaginary basis direction (1-based index), truncated at order ``nt``."""
        return float(real_value) + self._basis(basis_index)

    def _basis(self, i: int):
        return self._oti.e(int(i), order=self.nt)

    # ---- extraction ------------------------------------------------------
    def real(self, x) -> float:
        r = getattr(x, "real", None)
        if r is not None:
            return float(r)
        return float(x)

    def coeff(self, x, exponents: Sequence[int]) -> float:
        """The *raw* OTI coefficient of ``x`` along the direction with the given
        exponent multi-index (e.g. (2,0,0) -> e1^2). Multiply by the direction's
        recovery factor to obtain the true partial derivative, or use ``deriv``."""
        exps = tuple(int(e) for e in exponents)
        if all(e == 0 for e in exps):
            return self.real(x)
        humdir = self._exponents_to_index_list(exps)
        get_im = getattr(x, "get_im", None)
        if callable(get_im):
            return float(get_im(humdir))
        fn = getattr(self._oti, "get_im", None)
        if callable(fn):
            return float(fn(humdir, x))
        raise OtiApiError("OTILib build exposes no get_im; see "
                          "docs/otilib_api_inventory.md")

    def deriv(self, x, exponents: Sequence[int]) -> float:
        """The true partial derivative along ``exponents`` (OTILib ``get_deriv``
        already applies the recovery factor prod_i kappa_i!)."""
        exps = tuple(int(e) for e in exponents)
        if all(e == 0 for e in exps):
            return self.real(x)
        humdir = self._exponents_to_index_list(exps)
        get_deriv = getattr(x, "get_deriv", None)
        if callable(get_deriv):
            return float(get_deriv(humdir))
        fn = getattr(self._oti, "get_deriv", None)
        if callable(fn):
            return float(fn(humdir, x))
        # last resort: raw coefficient times our own recovery factor
        factor = 1
        for k in exps:
            factor *= math.factorial(k)
        return factor * self.coeff(x, exps)

    def set_coeff(self, x, exponents: Sequence[int], value: float):
        """Set the OTI coefficient of ``x`` along ``exponents`` to ``value``.
        Returns the (possibly new) scalar. Implemented additively via
        ``oti.e(...)`` so it works on any OTILib build; the delta ``value -
        current`` makes it a true set, not an accumulate."""
        exps = tuple(int(e) for e in exponents)
        if all(e == 0 for e in exps):
            try:
                x.real = float(value)
                return x
            except Exception:
                return x + (float(value) - self.real(x))
        current = self.coeff(x, exps)
        return x + (float(value) - current) * self._direction_number(exps)

    def zero_like(self, x):
        return self.scalar(0.0)

    def from_coefficients(self, real: float, coeff_map: Dict[Tuple[int, ...], float]):
        """Build an OTI scalar from a real part and a {exponents: value} map."""
        x = self.scalar(real)
        for exps, val in coeff_map.items():
            x = x + float(val) * self._direction_number(tuple(exps))
        return x

    # ---- direction bookkeeping ------------------------------------------
    def order_directions(self, p: int) -> List[Dict[str, Any]]:
        return order_directions(self.m, p)

    def num_coefficients_total(self) -> int:
        return num_coefficients_total(self.m, self.nt)

    def num_coefficients_order(self, p: int) -> int:
        return num_coefficients_order(self.m, p)

    # ---- internal helpers -----------------------------------------------
    def _exponents_to_index_list(self, exps: Tuple[int, ...]) -> List[int]:
        """OTILib addresses a direction by the flat, 1-based list of basis indices
        with multiplicity, e.g. (2,0,1) -> [1,1,3]. This is what ``oti.e`` and
        ``get_im`` accept."""
        idx: List[int] = []
        for i, k in enumerate(exps, start=1):
            idx.extend([i] * k)
        return idx

    def _direction_number(self, exps: Tuple[int, ...]):
        """The pure OTI number e_1^{k1} e_2^{k2} ... via a single ``oti.e`` call on
        the flat basis-index list."""
        humdir = self._exponents_to_index_list(exps)
        return self._oti.e(humdir, order=self.nt)


# --------------------------------------------------------------------------- #
class OtiUnavailableError(RuntimeError):
    pass


class OtiApiError(RuntimeError):
    pass
