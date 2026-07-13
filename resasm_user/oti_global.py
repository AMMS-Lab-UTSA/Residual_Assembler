"""Order-by-order sensitivity generation + solve for a *global* residual.

This is the HYPAD residual method (Aristizabal et al.) applied to a user-supplied
global residual, rather than to element assembly:

    seed all parameters   a_i* = a_i + e_i    (one OTI number carries every one)
    u* = converged real u
    for p = 1..q:
        R* = residual(u*, a*)               # evaluated with OTI scalars
        R^(p) = order-p coefficients of R*
        solve  T U^(p) = -R^(p)             # on the free DOFs
        inject U^(p) into u*                # BEFORE order p+1

Python path uses OTILib (`OtiContext`) directly. Black-box path delegates R^(p)
generation to the user's executable and only does the linear solve here.
"""

from __future__ import annotations

from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np


def build_direction_map(m: int, order: int) -> Dict[int, List[List[int]]]:
    from residual_core.algebra.otilib_adapter import order_directions
    return {p: [list(d["exponents"]) for d in order_directions(m, p)]
            for p in range(1, order + 1)}


def direction_labels(m: int, order: int) -> Dict[int, List[Dict[str, Any]]]:
    from residual_core.algebra.otilib_adapter import order_directions
    return {p: order_directions(m, p) for p in range(1, order + 1)}


def _free_mask(free_mask, ndof) -> np.ndarray:
    if free_mask is None:
        return np.ones(ndof, dtype=bool)
    return np.asarray(free_mask, dtype=bool)


# --------------------------------------------------------------------------- #
def solve_python(residual_fn: Callable, u, parameters: Dict[str, float],
                 order: int, T: np.ndarray, free_mask=None,
                 state=None, time=(0.0, 0.0), dtime=0.0
                 ) -> Tuple[Dict[int, np.ndarray], Dict[int, np.ndarray], Dict[str, Any]]:
    """OTILib order-by-order for a global residual callable. Requires OTILib."""
    from residual_core.algebra.otilib_adapter import (OtiContext, otilib_available,
                                                       otilib_status)
    if not otilib_available():
        raise OtiUnavailable(otilib_status()["error"])

    names = list(parameters.keys())
    m = len(names)
    ndof = int(np.asarray(u, float).size)
    ctx = OtiContext(num_bases=m, order=order)
    seeded = {n: ctx.seed(float(parameters[n]), i + 1) for i, n in enumerate(names)}
    u_star = [ctx.scalar(float(np.asarray(u, float).ravel()[j])) for j in range(ndof)]

    free = _free_mask(free_mask, ndof)
    T = np.asarray(T, float)
    Tff = T[np.ix_(free, free)]

    R_orders: Dict[int, np.ndarray] = {}
    U_orders: Dict[int, np.ndarray] = {}
    for p in range(1, order + 1):
        R_oti = list(residual_fn(u_star, seeded, state, time))
        if len(R_oti) != ndof:
            raise ValueError(
                "residual returned length %d but problem.unknowns is %d"
                % (len(R_oti), ndof))
        dirs = ctx.order_directions(p)
        R_p = np.zeros((ndof, len(dirs)))
        for col, d in enumerate(dirs):
            for j in range(ndof):
                R_p[j, col] = ctx.coeff(R_oti[j], d["exponents"])
        U_p = np.zeros((ndof, len(dirs)))
        U_p[free, :] = np.linalg.solve(Tff, -R_p[free, :])
        for col, d in enumerate(dirs):
            for j in range(ndof):
                u_star[j] = ctx.set_coeff(u_star[j], d["exponents"], U_p[j, col])
        R_orders[p] = R_p
        U_orders[p] = U_p

    diag = {"engine": "otilib", "basis_count": m, "truncation_order": order,
            "backend": otilib_status()}
    return R_orders, U_orders, diag


# --------------------------------------------------------------------------- #
def solve_executable(exe, u, parameters: Dict[str, float], order: int,
                     T: Optional[np.ndarray], free_mask=None,
                     state=None, time=(0.0, 0.0), dtime=0.0
                     ) -> Tuple[Dict[int, np.ndarray], Dict[int, np.ndarray],
                                np.ndarray, Dict[str, Any]]:
    """Black-box order-by-order. The executable returns R^(p); the linear solve
    is done here. The tangent may come pre-loaded (``T``) or from the response."""
    names = list(parameters.keys())
    m = len(names)
    ndof = int(np.asarray(u, float).size)
    full_dm = build_direction_map(m, order)

    R_orders: Dict[int, np.ndarray] = {}
    U_orders: Dict[int, np.ndarray] = {}
    u_star_coeffs: Dict[int, Any] = {}
    Tused = None if T is None else np.asarray(T, float)
    diag: Dict[str, Any] = {"engine": "executable"}

    for p in range(1, order + 1):
        resp = exe.eval_rhs(u, parameters, order=p, direction_map={p: full_dm[p]},
                            u_star_coefficients=u_star_coeffs, state=state,
                            time=time, dtime=dtime)
        if Tused is None and resp.get("tangent") is not None:
            Tused = np.asarray(resp["tangent"], float)
        R_p = resp["orders"].get(p)
        if R_p is None:
            raise ValueError(
                "executable did not return R^(%d) (response had orders %s)"
                % (p, sorted(resp["orders"].keys())))
        R_p = np.asarray(R_p, float)
        if Tused is None:
            raise ValueError(
                "no tangent available: provide tangent.file or return `tangent` "
                "in the black-box response")
        free = _free_mask(free_mask, ndof)
        Tff = Tused[np.ix_(free, free)]
        U_p = np.zeros_like(R_p)
        U_p[free, :] = np.linalg.solve(Tff, -R_p[free, :])
        u_star_coeffs[p] = U_p.tolist()
        R_orders[p] = R_p
        U_orders[p] = U_p
        if resp.get("diagnostics"):
            diag.setdefault("executable_diagnostics", {})[p] = resp["diagnostics"]

    return R_orders, U_orders, Tused, diag


class OtiUnavailable(RuntimeError):
    pass
