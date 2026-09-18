"""Hand-derived chain-rule sensitivities of the m5_cpflow update (slide 30).

An independent re-implementation, in NumPy, of the thermally activated flow
model in ``parameter_sensitivity/models/m5_cpflow/umat.for`` together with its
exact derivative by the implicit-function theorem. It shares no code with the
OTI transformation or with the compiled provider: the primal is re-coded from
the equations, and every derivative below was derived by hand.

Update (Voigt 11,22,33,12,13,23, engineering shear):

    sigma_tr = sigma_n + C : d_eps                  (C from E, nu only)
    s_tr     = sigma_vM(sigma_tr),  n = dev(sigma_tr) / s_tr
    F(dg)    = dg - dt * gdot(x) = 0,
    gdot     = gam0 * exp(-dG * (1 - x^p)^q),
    x        = (s_tr - 3 G dg) / (tau0 + H (eqp_n + dg))
    sigma    = sigma_tr - 3 G dg n,       eqp = eqp_n + dg

(the "previous stress minus the elastic correction produced by plastic slip").
With theta any of (tau0, dG, p, q, gam0, H):

    d dg / d theta = -(dF/dtheta + dF/ds_tr ds_tr/dtheta + dF/deqp_n deqp_n/dtheta) / (dF/d dg)
    d sigma / d theta = d sigma_n/d theta - 3G (d dg/d theta n + dg d n/d theta)
    d sigma_vM / d theta = (d sigma_vM / d sigma) : d sigma / d theta     (projection)

The Fortran source runs a fixed 60-step Newton with clamps; this module runs the
same iteration for the primal and refuses to differentiate an increment on
which a clamp is active (the derivative there is the clamp's, not the flow
rule's), reporting it instead.
"""

from __future__ import annotations

from typing import Dict, List, Sequence

import numpy as np

#: PROPS layout of m5_cpflow (1-based): E, nu, tau0, dG, p, q, gam0, H.
PARAMETER_PROPS = {"tau0": 3, "dG": 4, "p": 5, "q": 6, "gam0": 7, "H": 8}


def elastic_matrix(E: float, nu: float) -> np.ndarray:
    ebulk3 = E / (1.0 - 2.0 * nu)
    eg2 = E / (1.0 + nu)
    eg = eg2 / 2.0
    elam = (ebulk3 - eg2) / 3.0
    C = np.zeros((6, 6))
    C[:3, :3] = elam
    for i in range(3):
        C[i, i] = eg2 + elam
    for i in range(3, 6):
        C[i, i] = eg
    return C


def mises(s: np.ndarray) -> float:
    return float(np.sqrt(0.5 * ((s[0] - s[1]) ** 2 + (s[1] - s[2]) ** 2 + (s[2] - s[0]) ** 2)
                         + 3.0 * (s[3] ** 2 + s[4] ** 2 + s[5] ** 2)))


def dmises(s: np.ndarray) -> np.ndarray:
    vm = mises(s)
    if vm < 1e-30:
        return np.zeros(6)
    h = (s[0] + s[1] + s[2]) / 3.0
    g = np.empty(6)
    g[:3] = 1.5 * (s[:3] - h) / vm
    g[3:] = 3.0 * s[3:] / vm
    return g


def update(props: Sequence[float], stress_n: np.ndarray, eqp_n: float, deps: np.ndarray,
           dt: float, dstress_n: np.ndarray, deqp_n: np.ndarray, names: Sequence[str]):
    """One increment. ``dstress_n`` (6, P) and ``deqp_n`` (P,) are w.r.t. ``names``."""
    E, nu, tau0, dG, pexp, qexp, gam0, hard = [float(v) for v in props]
    eg3 = 3.0 * (E / (1.0 + nu) / 2.0)
    C = elastic_matrix(E, nu)
    s_tr = stress_n + C @ deps
    ds_tr = dstress_n.copy()                              # C does not depend on theta
    smises = mises(s_tr)
    info = {"branch": "elastic_zero_stress", "clamp": None}
    if smises <= 1e-12:
        return s_tr, eqp_n, ds_tr, deqp_n.copy(), 0.0, info
    shydro = s_tr[:3].mean()
    flow = np.concatenate([(s_tr[:3] - shydro) / smises, s_tr[3:] / smises])

    # --- primal: the same fixed 60-step Newton as the Fortran source --------
    dg = 0.0
    clamp = None
    for _ in range(60):
        sbar = smises - eg3 * dg
        sres = tau0 + hard * (eqp_n + dg)
        x = sbar / sres
        if x <= 0.0:
            x = 1e-8
        if x >= 1.0:
            gdot, dgdot = gam0, 0.0
        else:
            br = 1.0 - x ** pexp
            r = br ** qexp
            gdot = gam0 * np.exp(-dG * r)
            dxdd = (-eg3 * sres - sbar * hard) / (sres * sres)
            dgdot = gdot * (-dG) * (qexp * br ** (qexp - 1.0)) * (-pexp * x ** (pexp - 1.0) * dxdd)
        f = dg - dt * gdot
        df = 1.0 - dt * dgdot
        dg = dg - f / df
        clamp = None
        if dg < 0.0:
            dg, clamp = 0.0, "lower"
        if dg > 0.999 * smises / eg3:
            dg, clamp = 0.999 * smises / eg3, "upper"
    syield = smises - eg3 * dg
    stress = np.concatenate([flow[:3] * syield + shydro, flow[3:] * syield])
    eqp = eqp_n + dg

    # --- derivatives -------------------------------------------------------
    names = list(names)
    npar = len(names)
    gm = dmises(s_tr)
    dsm = gm @ ds_tr                                          # (P,)
    dshydro = ds_tr[:3].mean(axis=0)                          # (P,)
    dev_ds = np.vstack([ds_tr[:3] - dshydro[None, :], ds_tr[3:]])
    dflow = dev_ds / smises - np.outer(flow, dsm) / smises    # (6, P)

    sbar = smises - eg3 * dg
    sres = tau0 + hard * (eqp_n + dg)
    x_raw = sbar / sres
    unit = {name: np.eye(npar)[k] for k, name in enumerate(names)}
    zero = np.zeros(npar)
    if clamp == "upper":
        ddg = 0.999 * dsm / eg3
        info = {"branch": "clamp_upper", "clamp": clamp}
    elif clamp == "lower":
        ddg = zero.copy()
        info = {"branch": "clamp_lower", "clamp": clamp}
    elif x_raw >= 1.0:
        ddg = dt * unit.get("gam0", zero)                     # dg = dt*gam0
        info = {"branch": "saturated_x_ge_1", "clamp": None}
    else:
        x_clamped = x_raw <= 0.0
        x = 1e-8 if x_clamped else x_raw
        br = 1.0 - x ** pexp
        r = br ** qexp
        gdot = gam0 * np.exp(-dG * r)
        dgdot_dr = -dG * gdot
        dr_dx = qexp * br ** (qexp - 1.0) * (-pexp * x ** (pexp - 1.0))
        dr_dq = r * np.log(br)
        dr_dp = qexp * br ** (qexp - 1.0) * (-(x ** pexp) * np.log(x))
        if x_clamped:
            dx_dd = dx_ds = dx_de = dx_tau = dx_h = 0.0
        else:
            dx_dd = (-eg3 * sres - sbar * hard) / sres ** 2
            dx_ds = 1.0 / sres
            dx_de = -sbar * hard / sres ** 2
            dx_tau = -sbar / sres ** 2
            dx_h = -sbar * (eqp_n + dg) / sres ** 2
        f_dd = 1.0 - dt * dgdot_dr * dr_dx * dx_dd
        dx = dx_ds * dsm + dx_de * deqp_n + dx_tau * unit.get("tau0", zero) + dx_h * unit.get("H", zero)
        dr = dr_dx * dx + dr_dq * unit.get("q", zero) + dr_dp * unit.get("p", zero)
        dgd = dgdot_dr * dr - r * gdot * unit.get("dG", zero) + (gdot / gam0) * unit.get("gam0", zero)
        f_theta = -dt * dgd
        ddg = -f_theta / f_dd
        info = {"branch": "x_clamped_1e-8" if x_clamped else "flow", "clamp": None,
                "x": x_raw, "gdot": gdot}
    dsyield = dsm - eg3 * ddg
    dstress = np.vstack([dflow[:3] * syield + np.outer(flow[:3], dsyield) + dshydro[None, :],
                         dflow[3:] * syield + np.outer(flow[3:], dsyield)])
    deqp = deqp_n + ddg
    return stress, eqp, dstress, deqp, dg, info


def march(props: Sequence[float], increments: Sequence[Sequence[float]], dts: Sequence[float],
          names: Sequence[str]) -> Dict[str, object]:
    """Whole path from the virgin state. Returns stress, eqp, derivatives, sigma_vM and its derivative."""
    npar = len(names)
    stress = np.zeros(6); eqp = 0.0
    dstress = np.zeros((6, npar)); deqp = np.zeros(npar)
    out: Dict[str, List] = {"stress": [], "eqp": [], "dstress": [], "deqp": [], "dg": [],
                            "vm": [], "dvm": [], "branch": []}
    for deps, dt in zip(increments, dts):
        stress, eqp, dstress, deqp, dg, info = update(
            props, stress, eqp, np.asarray(deps, float), dt, dstress, deqp, names)
        out["stress"].append(stress.copy()); out["eqp"].append(eqp)
        out["dstress"].append(dstress.copy()); out["deqp"].append(deqp.copy())
        out["dg"].append(dg); out["vm"].append(mises(stress))
        out["dvm"].append(dmises(stress) @ dstress); out["branch"].append(info["branch"])
    return {k: (np.array(v) if k != "branch" else v) for k, v in out.items()}
