"""ctypes binding to the versioned material C ABI (``contract/resasm_mat_abi_v1.h``).

This is THE adapter between the Python residual tool and a compiled material
binary. It is deliberately identical for the reference provider and the real
JHU OTI binary -- only the shared-library path differs. No OTI/dual object ever
crosses this boundary: the binary returns the real response and, in separate
arrays, the first-order derivative coefficients for each seeded parameter.

Everything here is defensive about the ABI contract (dimensions, layouts,
NULL-able optional pointers, explicit failure) so a mis-shaped call fails loudly
instead of reading garbage.
"""

from __future__ import annotations

import ctypes
from typing import Dict, List, Optional, Sequence

import numpy as np

from residual_core.runtime import LibraryLoadError, load_shared_library

ABI_VERSION = 1
KIN_SMALL_STRAIN = 0
KIN_FINITE_STRAIN = 1

# status/return codes mirrored from resasm_mat_abi.h
_ERR_NAMES = {
    0: "OK", 1: "ERR_ABI_VERSION", 2: "ERR_DIMS", 3: "ERR_NULL",
    4: "ERR_SEED", 5: "ERR_KINEMATICS", 6: "ERR_STATE",
    7: "ERR_CONVERGENCE", 8: "ERR_INTERNAL",
}


class MatEvalError(RuntimeError):
    """A material-binary call returned a non-zero status, or the request was
    malformed before the call. Carries the numeric code and its name."""

    def __init__(self, code: int, context: str = ""):
        self.code = int(code)
        self.name = _ERR_NAMES.get(self.code, "ERR_%d" % self.code)
        msg = "material ABI call failed: %s (code %d)" % (self.name, self.code)
        if context:
            msg += " -- " + context
        super().__init__(msg)


class MatDesc(ctypes.Structure):
    """Mirror of ``resasm_mat_desc_t``."""
    _fields_ = [
        ("abi_version", ctypes.c_int),
        ("ntens", ctypes.c_int),
        ("nprops", ctypes.c_int),
        ("nstatev", ctypes.c_int),
        ("kinematics", ctypes.c_int),
        ("order", ctypes.c_int),
    ]


_c_double_p = ctypes.POINTER(ctypes.c_double)
_c_int_p = ctypes.POINTER(ctypes.c_int)


def _arr_d(seq) -> ctypes.Array:
    a = np.ascontiguousarray(seq, dtype=np.float64).ravel()
    return (ctypes.c_double * a.size)(*a.tolist())


def _arr_i(seq) -> ctypes.Array:
    a = np.ascontiguousarray(seq, dtype=np.int32).ravel()
    return (ctypes.c_int * a.size)(*a.tolist())


class MaterialABI:
    """A loaded material binary, callable through the v1 C ABI.

    Parameters
    ----------
    so_path : path to the compiled shared object exposing ``mat_eval_v1`` and
        ``mat_describe_v1``.
    """

    def __init__(self, so_path: str):
        self.so_path = str(so_path)
        try:
            self._handle = load_shared_library(self.so_path)
            self._lib = self._handle.lib
        except (LibraryLoadError, OSError) as exc:   # pragma: no cover - env dependent
            raise MatEvalError(3, "cannot load material binary %r: %s"
                               % (self.so_path, exc))
        self._bind()
        self.desc, self.model_id = self._describe()

    # -- symbol binding --------------------------------------------------- #
    def _bind(self) -> None:
        for sym in ("mat_eval_v1", "mat_describe_v1"):
            if not hasattr(self._lib, sym):
                raise MatEvalError(
                    8, "binary %r does not export %r (not a v%d ABI provider?)"
                       % (self.so_path, sym, ABI_VERSION))
        self._lib.mat_describe_v1.restype = ctypes.c_int
        self._lib.mat_describe_v1.argtypes = [
            ctypes.POINTER(MatDesc), ctypes.c_char_p, ctypes.c_int]
        self._lib.mat_eval_v1.restype = ctypes.c_int
        self._lib.mat_eval_v1.argtypes = [
            ctypes.POINTER(MatDesc),      # desc
            _c_double_p, _c_int_p, ctypes.c_int,   # props, seed, nseed
            _c_double_p,                  # kin
            ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,  # dkin, state_in, dstate_in
            _c_double_p,                  # time_data
            _c_double_p, _c_double_p,     # stress, dstress_dseed
            ctypes.c_void_p, ctypes.c_void_p,  # state_out, dstate_out
            _c_double_p,                  # ddsdde
            _c_int_p,                     # status
        ]

    def _describe(self):
        desc = MatDesc()
        buf = ctypes.create_string_buffer(128)
        rc = self._lib.mat_describe_v1(ctypes.byref(desc), buf, 128)
        if rc != 0:
            raise MatEvalError(rc, "mat_describe_v1")
        if desc.abi_version != ABI_VERSION:
            raise MatEvalError(1, "binary reports ABI v%d, tool speaks v%d"
                               % (desc.abi_version, ABI_VERSION))
        return desc, buf.value.decode("ascii", "replace")

    # -- metadata --------------------------------------------------------- #
    def metadata(self) -> Dict[str, object]:
        d = self.desc
        return {
            "model_id": self.model_id, "abi_version": d.abi_version,
            "ntens": d.ntens, "nprops": d.nprops, "nstatev": d.nstatev,
            "kinematics": d.kinematics, "order": d.order,
        }

    # -- evaluation ------------------------------------------------------- #
    def eval_point(
        self,
        props: Sequence[float],
        seed_indices: Sequence[int],
        kin: Sequence[float],
        *,
        state_in: Optional[Sequence[float]] = None,
        dstate_dseed_in: Optional[Sequence[float]] = None,
        dkin_dseed: Optional[Sequence[float]] = None,
        time_data: Sequence[float] = (0.0, 1.0, 0.0, 0.0),
    ) -> Dict[str, np.ndarray]:
        """Evaluate one material point for one increment.

        Returns ``{stress (ntens,), dstress_dseed (ntens, nseed), ddsdde
        (ntens, ntens), state (nstatev,), dstate_dseed (nstatev, nseed)}``.
        Raises :class:`MatEvalError` on a malformed request or a non-zero
        status from the binary.
        """
        d = self.desc
        ntens, nprops, nstatev = d.ntens, d.nprops, d.nstatev
        seed = list(int(s) for s in seed_indices)
        nseed = len(seed)
        if nseed < 1:
            raise MatEvalError(4, "at least one seed parameter is required")
        for s in seed:
            if s < 1 or s > nprops:
                raise MatEvalError(4, "seed index %d out of [1, %d]" % (s, nprops))
        if len(props) != nprops:
            raise MatEvalError(2, "props length %d != nprops %d" % (len(props), nprops))
        kin_len = ntens if d.kinematics == KIN_SMALL_STRAIN else 18
        if len(kin) != kin_len:
            raise MatEvalError(
                2, "kinematics length %d != expected %d for kinematics=%d"
                   % (len(kin), kin_len, d.kinematics))

        props_c = _arr_d(props)
        seed_c = _arr_i(seed)
        kin_c = _arr_d(kin)
        # pad/truncate to exactly 4 as a Python list BEFORE building the ctypes
        # array (slicing a ctypes array returns a list, which cannot be cast)
        _td = list(time_data)[:4] + [0.0] * max(0, 4 - len(time_data))
        time_c = _arr_d(_td)
        stress_c = (ctypes.c_double * ntens)()
        dstress_c = (ctypes.c_double * (ntens * nseed))()
        ddsdde_c = (ctypes.c_double * (ntens * ntens))()
        status_c = ctypes.c_int(0)

        def _optin(vals, n):
            if vals is None:
                return ctypes.c_void_p(0)
            if len(vals) != n:
                raise MatEvalError(2, "optional input length %d != %d" % (len(vals), n))
            return ctypes.cast(_arr_d(vals), ctypes.c_void_p)

        dkin_c = _optin(dkin_dseed, kin_len * nseed)
        state_in_c = _optin(state_in, nstatev)
        dstate_in_c = _optin(dstate_dseed_in, nstatev * nseed)
        # state outputs only if the model carries state
        if nstatev > 0:
            state_out_arr = (ctypes.c_double * nstatev)()
            dstate_out_arr = (ctypes.c_double * (nstatev * nseed))()
            state_out_c = ctypes.cast(state_out_arr, ctypes.c_void_p)
            dstate_out_c = ctypes.cast(dstate_out_arr, ctypes.c_void_p)
        else:
            state_out_arr = dstate_out_arr = None
            state_out_c = dstate_out_c = ctypes.c_void_p(0)

        rc = self._lib.mat_eval_v1(
            ctypes.byref(self.desc),
            ctypes.cast(props_c, _c_double_p),
            ctypes.cast(seed_c, _c_int_p), ctypes.c_int(nseed),
            ctypes.cast(kin_c, _c_double_p),
            dkin_c, state_in_c, dstate_in_c,
            ctypes.cast(time_c, _c_double_p),
            ctypes.cast(stress_c, _c_double_p),
            ctypes.cast(dstress_c, _c_double_p),
            state_out_c, dstate_out_c,
            ctypes.cast(ddsdde_c, _c_double_p),
            ctypes.byref(status_c),
        )
        # The RETURN CODE is the sole authority on validity (ABI rule 6); status
        # is informational detail and does not by itself invalidate the outputs.
        if rc != 0:
            raise MatEvalError(rc, "mat_eval_v1 (status detail %d)" % status_c.value)

        out = {
            "stress": np.array(stress_c, dtype=float),
            "dstress_dseed": np.array(dstress_c, dtype=float).reshape(ntens, nseed),
            "ddsdde": np.array(ddsdde_c, dtype=float).reshape(ntens, ntens),
        }
        if nstatev > 0:
            out["state"] = np.array(state_out_arr, dtype=float)
            out["dstate_dseed"] = np.array(dstate_out_arr, dtype=float).reshape(nstatev, nseed)
        else:
            out["state"] = np.zeros(0)
            out["dstate_dseed"] = np.zeros((0, nseed))
        return out
