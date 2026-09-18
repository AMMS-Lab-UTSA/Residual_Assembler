"""Independent references for a history replay.

Nothing here uses the OTI derivatives it checks:

``tangent_check``
    At sampled integration points of chosen increments, central differences
    of the ORIGINAL UMAT (the regular path compiled into the provider object)
    with respect to DSTRAN, at the replayed start-of-increment state, compared
    with the provider's DDSDDE. Two step sizes give the plateau.

``whole_model_fd``
    For every parameter, the whole model is re-equilibrated in Python with the
    ORIGINAL UMAT at p(1 +/- h) over the same increments (Newton to a tight
    residual; the Newton matrix may be the provider tangent because the
    converged state depends only on the residual), and central differences of
    displacements, reactions, stresses and state variables are formed for a
    ladder of relative steps h. For every parameter, field and increment the
    adjacent pair of steps with the smallest spread is the plateau (larger h:
    truncation error, smaller h: Newton/roundoff error); the OTI derivative is
    compared with the finer estimate of that pair. The OTI sensitivities of
    the same discrete problem come from an OTI solve at the nominal parameters.
"""
from __future__ import annotations

import time as _time
from typing import Dict, List, Sequence

import numpy as np

from .history import HistoryEngine, HistoryResult, run_history, von_mises, von_mises_gradient


def _increment_inputs(engine: HistoryEngine, result: HistoryResult, number: int):
    previous = result.increments[number - 2] if number > 1 else None
    npts = engine.ne * engine.nq
    stress = previous.stress.reshape(npts, 6) if previous else np.zeros((npts, 6))
    state = previous.state.reshape(npts, -1) if previous else np.zeros((npts, engine.material.nstatev))
    u_prev = previous.u if previous else np.zeros(engine.ndof)
    u_now = result.increments[number - 1].u
    strain_prev = engine.flat(engine.strain(u_prev))
    dstran = engine.flat(engine.strain(u_now)) - strain_prev
    t_prev = result.times[number - 1]
    dt = result.times[number] - t_prev
    return stress, state, strain_prev, dstran, np.array([t_prev, t_prev]), dt


def tangent_check(engine: HistoryEngine, result: HistoryResult, increments: Sequence[int],
                  sample: int = 12, steps=(2e-9, 1e-9)) -> Dict:
    """DDSDDE of the provider vs central FD of the ORIGINAL UMAT at sampled points."""
    material = engine.material
    worst_error, worst_spread, checked = 0.0, 0.0, 0
    for number in increments:
        stress, state, strain_prev, dstran, time2, dt = _increment_inputs(engine, result, number)
        vm = von_mises(result.increments[number - 1].stress).reshape(-1)
        order = np.argsort(vm)
        picks = np.unique(np.concatenate([order[-sample // 2:], order[:: max(1, len(order) // (sample // 2))][: sample // 2]]))
        sub = lambda a: a[picks]
        common = dict(time=time2, dtime=dt, coords=sub(engine.ip_coords.reshape(-1, 3)),
                      celent=sub(engine.celent), noel=sub(engine.noel), npt=sub(engine.npt), kinc=number)
        npar, ns = material.nparam, material.nstatev
        zero = np.zeros((len(picks), 6, npar))
        tangent = material.total(result.props, sub(stress), sub(state), sub(strain_prev), sub(dstran),
                                 dstress_in=zero, dstate_in=np.zeros((len(picks), ns, npar)),
                                 stran_dp=zero, dstran_dp=zero, **common)["ddsdde"]
        estimates = []
        for h in steps:
            column = np.empty_like(tangent)
            for axis in range(6):
                shift = np.zeros_like(sub(dstran))
                shift[:, axis] = h
                plus = material.regular(result.props, sub(stress), sub(state), sub(strain_prev),
                                        sub(dstran) + shift, **common)["stress"]
                minus = material.regular(result.props, sub(stress), sub(state), sub(strain_prev),
                                         sub(dstran) - shift, **common)["stress"]
                column[:, :, axis] = (plus - minus) / (2 * h)
            estimates.append(column)
        scale = np.abs(estimates[1]).max(axis=(1, 2), keepdims=True)
        worst_error = max(worst_error, float((np.abs(tangent - estimates[1]) / scale).max()))
        worst_spread = max(worst_spread, float((np.abs(estimates[0] - estimates[1]) / scale).max()))
        checked += len(picks)
    return {"points_checked": checked, "increments": list(increments), "steps": list(steps),
            "max_relative_error": worst_error, "fd_plateau_spread": worst_spread,
            "reference": "central FD of the ORIGINAL UMAT w.r.t. DSTRAN at the replayed state"}


def _snapshot(engine: HistoryEngine, result: HistoryResult):
    u = np.stack([r.u for r in result.increments])
    rf = np.stack([r.reaction for r in result.increments])
    stress = np.stack([r.stress for r in result.increments])
    state = np.stack([r.state for r in result.increments])
    return {"U": u, "RF": rf, "S": stress, "SDV": state, "MISES": von_mises(stress)}


def whole_model_fd(engine: HistoryEngine, times, *, parameters: Sequence[str] = None,
                   steps=(1e-3, 3e-4, 1e-4, 3e-5, 1e-5), rtol: float = 1e-13, reference: HistoryResult = None) -> Dict:
    """Central FD of the whole equilibrated model with the ORIGINAL UMAT (see module doc)."""
    material, model = engine.material, engine.model
    parameters = list(parameters or material.params)
    started = _time.perf_counter()
    if reference is None:
        reference = run_history(engine, times=times, sensitivities=True, rtol=rtol)
    nominal = material.parameter_values(model.props)
    base_props = np.array(model.props, dtype=float)
    fd: Dict[str, Dict[float, Dict[str, np.ndarray]]] = {}
    for name in parameters:
        j = material.params.index(name)
        fd[name] = {}
        for h in steps:
            runs = []
            for sign in (+1, -1):
                props = material.props_with(base_props, {name: nominal[j] * (1 + sign * h)})
                solved = run_history(engine, times=times, props=props, sensitivities=False,
                                     material_path="regular", newton_tangent="oti", rtol=rtol)
                runs.append(_snapshot(engine, solved))
            delta = 2 * h * nominal[j]
            fd[name][h] = {key: (runs[0][key] - runs[1][key]) / delta for key in runs[0]}
    oti = {
        "U": np.stack([r.du for r in reference.increments]),
        "RF": np.stack([r.dreaction for r in reference.increments]),
        "S": np.stack([r.dstress for r in reference.increments]),
        "SDV": np.stack([r.dstate for r in reference.increments]),
    }
    stress = np.stack([r.stress for r in reference.increments])
    oti["MISES"] = np.einsum("neqa,neqam->neqm", von_mises_gradient(stress), oti["S"])
    constrained = engine.constrained
    values = _snapshot(engine, reference)
    comparison = {}
    for name in parameters:
        j = material.params.index(name)
        rows = {}
        for key in ("U", "RF", "S", "SDV", "MISES"):
            ours, q = oti[key][..., j], values[key]
            ladder = [fd[name][h][key] for h in steps]
            if key == "RF":
                ours, q = ours[:, constrained], q[:, constrained]
                ladder = [estimate[:, constrained] for estimate in ladder]
            per_increment = []
            for n in range(ours.shape[0]):
                natural = max(np.abs(q[n]).max(), 1e-300) / abs(nominal[j])     # |q|/|p|
                pairs = []
                for k in range(len(steps) - 1):
                    coarse, fine = ladder[k][n], ladder[k + 1][n]
                    scale = max(np.abs(fine).max(), 1e-300)
                    pairs.append((float(np.abs(coarse - fine).max() / scale), k, fine, scale))
                spread, k, fine, scale = min(pairs, key=lambda item: item[0])
                per_increment.append({
                    "error": float(np.abs(ours[n] - fine).max() / scale), "spread": spread,
                    "steps": [steps[k], steps[k + 1]],
                    "weighted_fd": float(np.abs(fine).max() / natural),
                    "weighted_oti": float(np.abs(ours[n]).max() / natural),
                    "weighted_error": float(np.abs(ours[n] - fine).max() / natural)})
            rows[key] = per_increment
        comparison[name] = rows
    return {"parameters": parameters, "steps": list(steps), "rtol": rtol, "comparison": comparison,
            "seconds": _time.perf_counter() - started, "reference_result": reference,
            "definition": ("per increment and field, over the adjacent step pair with the smallest "
                           "spread: error = max|OTI - FD(h_fine)| / max|FD(h_fine)|, "
                           "spread = max|FD(h_coarse) - FD(h_fine)| / max|FD(h_fine)|, weighted_* = "
                           "max|.| |p| / max|q| (the derivative on the scale of the field itself); "
                           "FD from ORIGINAL-UMAT equilibria at p(1 +/- h)")}


def summarize_fd(result: Dict, *, zero_level: float = 1e-6) -> Dict[str, Dict[str, Dict[str, float]]]:
    """Worst error and spread per parameter and field.

    An increment whose FD derivative is below ``zero_level`` on the field's
    own scale (|p| max|dq/dp| / max|q|) is a zero reference (e.g. du/dE under
    pure displacement control while elastic, or a yield-stress derivative
    before yield): there OTI and FD are compared on the field scale
    (|OTI - FD| |p| / max|q|) instead of relative to FD noise.
    """
    summary = {}
    for name, rows in result["comparison"].items():
        summary[name] = {}
        for key, values in rows.items():
            live = [v for v in values if v["weighted_fd"] > zero_level]
            zeros = [v for v in values if v["weighted_fd"] <= zero_level]
            summary[name][key] = {
                "nonzero_increments": len(live), "zero_increments": len(zeros),
                "max_error": max((v["error"] for v in live), default=0.0),
                "max_spread": max((v["spread"] for v in live), default=0.0),
                "zero_reference_max_weighted_oti": max((v["weighted_oti"] for v in zeros), default=0.0),
                "zero_reference_max_weighted_error": max((v["weighted_error"] for v in zeros), default=0.0)}
    return summary
