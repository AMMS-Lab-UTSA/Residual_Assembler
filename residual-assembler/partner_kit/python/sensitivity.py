"""Core kit logic: real assembly, Newton solve, tangent acquisition, RHS
generation (order 1) and the sensitivity solve. All runs locally on the partner
machine against their provider; nothing here needs the partner's source or mesh.
"""

from __future__ import annotations

from typing import Any, Dict, Optional, Tuple

import numpy as np

from .hypercomplex import imag_part, make_seed
from .residual_provider import (ResidualProvider, GlobalResidualProvider,
                                ElementResidualProvider, BlackBoxResidualProvider)


# --------------------------------------------------------------------------- #
# real residual + tangent assembly (works for element and global providers)
# --------------------------------------------------------------------------- #
def assemble_real(provider: ResidualProvider, U, state, values, time, dtime,
                  want_tangent: bool = False) -> Tuple[np.ndarray, Optional[np.ndarray]]:
    U = np.asarray(U, float)
    ndof = provider.ndof
    R = np.zeros(ndof)
    T = np.zeros((ndof, ndof)) if want_tangent else None
    contributed_T = False

    if isinstance(provider, GlobalResidualProvider):
        r = provider.eval_global_residual(U, state, values, time, dtime)
        R[:] = [float(x) for x in r]
    elif isinstance(provider, ElementResidualProvider):
        for eid in provider.elements():
            dofs = np.asarray(provider.element_dof_map(eid), int)
            u_e = U[dofs]
            st = provider.element_state(eid)
            r_e, k_e, _s = provider.eval_element_residual(
                eid, u_e, st, values, time, dtime)
            for j, val in enumerate(r_e):
                R[dofs[j]] += float(val)
            if want_tangent and k_e is not None:
                T[np.ix_(dofs, dofs)] += np.asarray(k_e, float)
                contributed_T = True
    else:
        raise TypeError("unknown provider level: %s" % provider.kind())

    if want_tangent and not contributed_T and isinstance(provider, GlobalResidualProvider):
        T = None  # global providers supply T via get_tangent/apply_tangent
    if want_tangent and T is not None and not contributed_T and \
            isinstance(provider, ElementResidualProvider):
        T = None
    return R, T


def acquire_tangent(provider: ResidualProvider, U, state, values, time, dtime,
                    external: Optional[np.ndarray] = None,
                    reconstruct_limit: int = 4000):
    """Return (T_dense_or_None, source_str). Tries, in order: external file,
    element-assembled T, provider.get_tangent, matrix-free reconstruction."""
    if external is not None:
        return np.asarray(external, float), "external-file"
    # element-assembled
    _R, T = assemble_real(provider, U, state, values, time, dtime, want_tangent=True)
    if T is not None:
        return T, "backend-assembled"
    # dense provider tangent
    Tg = provider.get_tangent(U, state, values, time, dtime)
    if Tg is not None:
        return np.asarray(Tg, float), "provider-get_tangent"
    # matrix-free reconstruction (small problems only)
    ndof = provider.ndof
    if ndof <= reconstruct_limit:
        cols = []
        ok = True
        for i in range(ndof):
            e = np.zeros(ndof); e[i] = 1.0
            tx = provider.apply_tangent(U, state, values, e, time, dtime)
            if tx is None:
                ok = False
                break
            cols.append(np.asarray(tx, float))
        if ok:
            return np.column_stack(cols), "matrix-free-reconstructed"
    return None, "unavailable"


# --------------------------------------------------------------------------- #
# real Newton solve (step 1) — optional convenience
# --------------------------------------------------------------------------- #
def newton_solve(provider: ResidualProvider, u0=None, state=None, time=(0.0, 0.0),
                 dtime=0.0, max_iter: int = 100, tol: float = 1e-11):
    values = provider.parameter_values()
    ndof = provider.ndof
    U = np.ones(ndof) if u0 is None else np.asarray(u0, float).copy()
    free = provider.free_mask()
    if free is None:
        free = np.ones(ndof, dtype=bool)
    for _ in range(max_iter):
        R, _T = assemble_real(provider, U, state, values, time, dtime, want_tangent=False)
        rn = float(np.linalg.norm(R[free])) if free.any() else float(np.linalg.norm(R))
        if rn < tol:
            break
        T, _src = acquire_tangent(provider, U, state, values, time, dtime)
        if T is None:
            raise RuntimeError("Newton needs a tangent; provider supplied none")
        Tff = T[np.ix_(free, free)]
        dU = np.linalg.solve(Tff, -R[free])
        U = U.copy()
        U[free] += dU
    return U


# --------------------------------------------------------------------------- #
# RHS generation: R^(1) = dR/dp
# --------------------------------------------------------------------------- #
def generate_rhs(provider: ResidualProvider, U, order: int = 1, state=None,
                 time=(0.0, 0.0), dtime=0.0, algebra: str = "dual1"):
    """Return (R_p [ndof x m], diagnostics, raw_response_or_None)."""
    U = np.asarray(U, float)
    params = list(provider.parameters)
    m = len(params)
    ndof = provider.ndof
    R = np.zeros((ndof, m))
    diag: Dict[str, Any] = {"algebra": algebra, "provider": provider.name,
                            "order": order, "hypercomplex_ready": True, "notes": []}
    values = provider.parameter_values()

    # black-box: the partner's executable does the overloading internally
    if isinstance(provider, BlackBoxResidualProvider):
        seeds = {p: i + 1 for i, p in enumerate(params)}
        resp = provider.eval_rhs(order, U, values, seeds, time, dtime)
        R = provider._runner.residual_coefficients(resp)
        diag["source"] = "blackbox-executable"
        diag["response_diagnostics"] = resp.get("diagnostics", {})
        return R, diag, resp

    if order != 1:
        raise NotImplementedError("the kit generates order 1 only (dual1); higher "
                                  "orders need an OTI/HYPAD scalar or a black-box "
                                  "provider that returns R^(p).")

    for col, p in enumerate(params):
        seeded = dict(values)
        seeded[p] = make_seed(values[p], algebra, order)
        try:
            if isinstance(provider, GlobalResidualProvider):
                r = provider.eval_global_residual(U, state, seeded, time, dtime)
                for i, val in enumerate(r):
                    R[i, col] += imag_part(val)
            elif isinstance(provider, ElementResidualProvider):
                for eid in provider.elements():
                    dofs = np.asarray(provider.element_dof_map(eid), int)
                    u_e = U[dofs]
                    st = provider.element_state(eid)
                    r_e, _k, _s = provider.eval_element_residual(
                        eid, u_e, st, seeded, time, dtime)
                    for j, val in enumerate(r_e):
                        R[dofs[j], col] += imag_part(val)
            else:
                raise TypeError("unknown provider level")
        except (TypeError, ValueError) as exc:
            diag["hypercomplex_ready"] = False
            diag["notes"].append("parameter %r not hypercomplex-ready: %s" % (p, exc))
    diag["source"] = "local-%s" % algebra
    return R, diag, None


def solve_sensitivity(T, R_p, free_mask=None):
    """Solve T U^(p) = -R^(p) on the free partition (zeros at prescribed DOFs)."""
    T = np.asarray(T, float)
    R_p = np.asarray(R_p, float)
    rhs = -R_p
    U = np.zeros_like(rhs)
    if free_mask is None:
        return np.linalg.solve(T, rhs)
    free = np.asarray(free_mask, bool)
    Tff = T[np.ix_(free, free)]
    U[free, :] = np.linalg.solve(Tff, rhs[free, :])
    return U
