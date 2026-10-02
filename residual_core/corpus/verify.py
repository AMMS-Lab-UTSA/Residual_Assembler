"""Finite-difference adjudication, entry by entry (B2 rule; Vera B1 review B/C/D).

A comparison is one analytic array ``A`` (the OTI-derived derivative) against
centred differences ``D(h) = (q(p+h) - q(p-h)) / 2h`` of an INDEPENDENT
evaluation (the original routine) over a ladder of steps. Every entry ``e`` is
judged on its own; no column norm, array maximum or path maximum enters the
tolerance of an entry.

1. FD-only plateau. The reference is judged from the FD sequence alone, never
   from its agreement with A. For entry e, a plateau is a run of >= ``PLATEAU``
   consecutive steps whose successive estimates differ by at most
   ``max(rtol |median|, noise_e(h))``, with the round-off model
   ``noise_e(h) = 8 eps_eval (smax_e + |D_e(h)| h) / h`` (``smax_e`` the
   magnitude that cancels in the evaluated entry, ``eps_eval`` the precision the
   ORIGINAL computes in). The longest plateau is taken (ties: larger steps).
   ``D_e`` = its median, ``u_e`` = its spread ``max |D_e(h) - D_e|``.
2. Structural zero: ``|D_e| <= 1e-12 max|D|``; judged by
   ``|A_e| <= atol + 1e-12 max|D| + 2 u_e``. Zero within resolution: no
   plateau, but ``|D_e(h)| <= noise_e(h)`` at EVERY step (the FD sees nothing);
   judged by ``|A_e| <= atol + min_h noise_e(h)`` (any larger value would have
   shown in the difference at the largest step).
3. Resolution: ``u_e > 1e-3 |D_e|`` (the FD cannot resolve the entry) ->
   UNRESOLVED, never verified. No plateau -> UNRESOLVED.
4. Entry passes iff ``|A_e - D_e| <= atol + rtol |D_e| + 2 u_e``.
5. A failing or unresolved entry is examined for a kink: one-sided differences
   D+ and D- of a smooth q differ by O(h); if ``|D+ - D-|`` exceeds
   ``10 rtol |D|`` at h=1e-3, has not fallen 10x by h=1e-5, and stands 10x
   clear of the round-off bound at both steps, the entry is NONSMOOTH (kept
   out of verified and failed).
6. Comparison verdict: ``failed`` if any entry failed; else ``nonsmooth`` if
   any significant entry is nonsmooth; else ``unresolved`` if any SIGNIFICANT
   entry (|D_e| > 1e-6 max|D|) is unresolved, or no entry was resolved; else
   ``verified``. A comparison is never verified on partial coverage of its
   significant entries; the counts are kept in the result.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence

import numpy as np

__all__ = ["STEPS", "PLATEAU", "EPS64", "EPS32", "adjudicate", "summarise"]

STEPS: Sequence[float] = (1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 1e-7)
PLATEAU = 3
EPS64 = float(np.finfo(np.float64).eps)
EPS32 = float(np.finfo(np.float32).eps)
SIGNIFICANT = 1e-6
UNRESOLVED_SPREAD = 1e-3
ZERO = 1e-12


def _plateaus(est: np.ndarray, noise: np.ndarray, rtol: float):
    """est, noise (k, E) -> (start, end) index per entry of the chosen plateau
    (-1 where none), plateau median D and spread u.

    A step can belong to a plateau only if its round-off bound is below
    ``UNRESOLVED_SPREAD |D(h)|`` (a step that cannot resolve the entry to 1e-3
    cannot certify it). Among the admissible windows of >= PLATEAU steps the
    longest is taken, ties broken by the smaller spread. An entry whose
    differences are exactly 0 at every step is a structural zero (D = u = 0).
    """
    k, E = est.shape
    finite = np.isfinite(est)
    admissible = finite & (noise <= UNRESOLVED_SPREAD * np.abs(np.nan_to_num(est)))
    start = np.full(E, -1)
    end = np.full(E, -1)
    D = np.full(E, np.nan)
    u = np.full(E, np.nan)
    exact_zero = np.all(finite, axis=0) & np.all(est == 0.0, axis=0)
    D[exact_zero], u[exact_zero] = 0.0, 0.0
    start[exact_zero], end[exact_zero] = 0, k - 1
    found = exact_zero.copy()
    for length in range(k, PLATEAU - 1, -1):
        best_u = np.full(E, np.inf)
        best = np.full(E, -1)
        for a in range(0, k - length + 1):
            b = a + length - 1
            seg = np.nan_to_num(est[a:b + 1])
            med = np.median(seg, axis=0)
            lim = np.maximum(rtol * np.abs(med), noise[a + 1:b + 1])
            ok = (np.all(admissible[a:b + 1], axis=0)
                  & np.all(np.abs(np.diff(seg, axis=0)) <= lim, axis=0) & ~found)
            spread = np.max(np.abs(seg - med), axis=0)
            better = ok & (spread < best_u)
            best_u[better], best[better] = spread[better], a
        hit = best >= 0
        for e in np.flatnonzero(hit):
            a = best[e]
            seg = est[a:a + length, e]
            start[e], end[e] = a, a + length - 1
            D[e] = np.median(seg)
            u[e] = np.max(np.abs(seg - D[e]))
        found |= hit
        if found.all():
            break
    return start, end, D, u


def adjudicate(analytic, plus: Dict[float, np.ndarray], minus: Dict[float, np.ndarray],
               base, steps: Dict[float, float], *, rtol: float, atol: float,
               value_scale=None, eps_eval: float = EPS64,
               step_level: Optional[Dict[float, float]] = None) -> dict:
    """``plus/minus[s]`` evaluated at ``p +/- steps[s]`` (s the relative step;
    a missing or None entry -- e.g. a re-solve that failed -- breaks the ladder
    there). ``value_scale``: per-entry magnitude cancelling in the evaluated
    quantity (default |base|). ``step_level[s]``: relative error level of the
    evaluations at step s when it exceeds eps_eval (e.g. the free residual a
    nonlinear re-solve actually achieved); the round-off bound at that step
    becomes 8 max(eps_eval, level) (...)/h."""
    A = np.asarray(analytic, dtype=float).ravel()
    order = [s for s in sorted(steps, reverse=True)]
    hs = np.array([steps[s] for s in order])
    base = np.asarray(base, dtype=float).ravel()
    smax = np.abs(base) if value_scale is None else np.broadcast_to(
        np.asarray(value_scale, dtype=float).ravel() if np.ndim(value_scale) else value_scale,
        A.shape).astype(float)
    k, E = len(order), A.size
    est = np.full((k, E), np.nan)
    dplus = np.full((k, E), np.nan)
    dminus = np.full((k, E), np.nan)
    for i, s in enumerate(order):
        if plus.get(s) is None or minus.get(s) is None:
            continue
        qp = np.asarray(plus[s], dtype=float).ravel()
        qm = np.asarray(minus[s], dtype=float).ravel()
        est[i] = (qp - qm) / (2.0 * hs[i])
        dplus[i] = (qp - base) / hs[i]
        dminus[i] = (base - qm) / hs[i]
    levels = np.array([max(eps_eval, (step_level or {}).get(s, 0.0)) for s in order])
    noise = 8.0 * levels[:, None] * (smax[None, :] + np.abs(np.nan_to_num(est)) * hs[:, None]) / hs[:, None]
    start, end, D, u = _plateaus(est, noise, rtol)
    # zero within resolution: no plateau, but every difference lies inside its
    # round-off bound -- the FD sees nothing; judged against the smallest bound
    finite_all = np.all(np.isfinite(est), axis=0)
    inside = finite_all & np.all(np.abs(np.nan_to_num(est)) <= noise, axis=0)
    fd_zero = np.isnan(D) & inside
    zero_bound = np.min(noise, axis=0)
    D[fd_zero], u[fd_zero] = 0.0, 0.0
    start[fd_zero], end[fd_zero] = 0, len(order) - 1
    has = ~np.isnan(D)
    Dmax = float(np.nanmax(np.abs(D))) if has.any() else 0.0
    zero = has & (np.abs(D) <= ZERO * Dmax)
    resolved = has & ~zero & (u <= UNRESOLVED_SPREAD * np.abs(D))
    verdict = np.full(E, "unresolved", dtype=object)
    tol_zero = atol + ZERO * Dmax + 2.0 * np.nan_to_num(u) + np.where(fd_zero, zero_bound, 0.0)
    verdict[zero] = np.where(np.abs(A[zero]) <= tol_zero[zero], "pass", "fail")
    tol = atol + rtol * np.abs(np.nan_to_num(D)) + 2.0 * np.nan_to_num(u)
    err = np.abs(A - np.nan_to_num(D))
    verdict[resolved] = np.where(err[resolved] <= tol[resolved], "pass", "fail")
    # kinks: failing / unresolved entries whose one-sided asymmetry does not shrink
    by = {s: i for i, s in enumerate(order)}
    if 1e-3 in by and 1e-5 in by:
        i3, i5 = by[1e-3], by[1e-5]
        ref = np.where(has, np.abs(np.nan_to_num(D)), np.abs(np.nan_to_num(est[i3])))
        big = np.abs(dplus[i3] - dminus[i3])
        small = np.abs(dplus[i5] - dminus[i5])
        # a kink's one-sided jump does not shrink with h; round-off grows as 1/h,
        # so both asymmetries must also stand clear of the round-off bound
        kink = ((big > 10.0 * rtol * np.maximum(ref, atol / max(rtol, 1e-300)))
                & (small > 0.1 * big)
                & (big > 10.0 * noise[i3]) & (small > 10.0 * noise[i5]))
        kink &= np.isin(verdict, ["fail", "unresolved"])
        verdict[kink] = "nonsmooth"
    absD = np.abs(np.nan_to_num(D, nan=0.0))
    estmax = np.nanmax(np.abs(est), axis=0) if np.isfinite(est).any() else np.zeros(E)
    size = np.maximum(np.where(has, absD, np.nan_to_num(estmax)), np.abs(A))
    scale_all = float(np.max(size)) if E else 0.0
    significant = size > SIGNIFICANT * scale_all
    counts = {v: int(np.sum(verdict == v)) for v in ("pass", "fail", "unresolved", "nonsmooth")}
    counts["structural_zero"] = int(np.sum(zero & ~fd_zero))
    counts["zero_within_resolution"] = int(np.sum(fd_zero))
    counts["entries"] = int(E)
    counts["significant"] = int(np.sum(significant))
    counts["significant_unresolved"] = int(np.sum(significant & (verdict == "unresolved")))
    if counts["fail"]:
        status = "failed"
    elif np.any(significant & (verdict == "nonsmooth")):
        status = "nonsmooth"
    elif counts["significant_unresolved"] or not np.any(resolved | zero):
        status = "unresolved"
    else:
        status = "verified"
    judged = resolved | zero
    rel = np.where(judged, err / np.maximum(np.abs(np.nan_to_num(D)), atol / max(rtol, 1e-300)), 0.0)
    # error / tolerance per judged entry (the manifest's max_error, tolerance 1)
    ratio = np.where(zero, np.abs(A) / tol_zero, np.where(resolved, err / tol, 0.0))
    worst = np.argsort(-np.where(judged, rel, -1.0))[:3]
    lengths = np.where(has, end - start + 1, 0)
    return {
        "verdict": status, "counts": counts,
        "max_error_ratio": float(np.max(ratio)) if E and judged.any() else None,
        # plateau statistics over NON-ZERO resolved entries only; a zero entry's
        # "plateau" (all steps inside the round-off bound, or a flat zero run)
        # says nothing about FD convergence and is reported on its own
        "min_plateau_observed": int(np.min(lengths[resolved])) if resolved.any() else None,
        "min_plateau_zero_entries": int(np.min(lengths[zero])) if zero.any() else None,
        "max_abs": float(np.max(np.where(judged, err, 0.0))) if E else 0.0,
        "max_rel": float(np.max(rel)) if E else 0.0,
        "resolution_max": float(np.max(np.where(resolved, np.nan_to_num(u) / np.maximum(absD, 1e-300), 0.0)))
        if E else 0.0,
        "plateau_min_steps": int(np.min(lengths[resolved])) if resolved.any() else 0,
        "plateau_steps_seen": sorted({(order[start[e]], order[end[e]]) for e in np.flatnonzero(resolved)})[:6],
        "worst": [{"entry": int(e), "analytic": float(A[e]), "fd": float(D[e]) if has[e] else None,
                   "spread": float(u[e]) if has[e] else None, "verdict": str(verdict[e])}
                  for e in worst if E],
        "analytic_scale": float(np.max(np.abs(A))) if E else 0.0,
        "fd_scale": Dmax, "eps_eval": eps_eval, "rtol": rtol, "atol": atol,
        "steps": order,
    }


def summarise(comparisons: List[dict]) -> dict:
    """Roll comparisons up into one record status (the campaign's status set).

    ``failed``: any comparison failed. ``verified``: none failed, none
    unresolved, and the verified comparisons outnumber the nonsmooth ones (a
    nonsmooth comparison sits on a kink of the material response, where the
    derivative does not exist; they are listed, never counted verified).
    Otherwise ``unsupported`` with ``reason_class`` ``fd_reference_unresolved``
    (the FD reference could not resolve a significant entry at this precision)
    or ``nonsmooth_dominant`` -- never verified on partial coverage.
    """
    counts = {"verified": 0, "failed": 0, "nonsmooth": 0, "unresolved": 0}
    for c in comparisons:
        counts[c["verdict"]] += 1
    judged = [c for c in comparisons if c["verdict"] in ("verified", "failed")]
    reason = None
    if not comparisons:
        status = "not_attempted"
    elif counts["failed"]:
        status = "failed"
    elif counts["unresolved"]:
        status, reason = "unsupported", "fd_reference_unresolved"
    elif counts["verified"] > counts["nonsmooth"]:
        status = "verified"
    else:
        status, reason = "unsupported", "nonsmooth_dominant"
    ratios = [c.get("max_error_ratio") for c in judged if c.get("max_error_ratio") is not None]
    plateaus = [c.get("min_plateau_observed") for c in judged
                if c.get("min_plateau_observed") is not None]
    return {"status": status, "reason_class": reason, "counts": counts,
            "comparisons": len(comparisons),
            "max_error_ratio": max(ratios) if ratios else None,
            "min_plateau_observed": min(plateaus) if plateaus else None,
            "max_abs": max((c["max_abs"] for c in judged if c.get("max_abs") is not None), default=None),
            "max_rel": max((c["max_rel"] for c in judged if c.get("max_rel") is not None), default=None),
            "resolution_max": max((c.get("resolution_max") or 0.0 for c in judged), default=None)}
