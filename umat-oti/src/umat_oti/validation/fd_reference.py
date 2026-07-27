"""The canonical finite-difference reference for material-parameter derivatives.

This is the single implementation behind *every* OTI-vs-FD claim Program 1
makes: the transformer's own build-time self-check
(``oti_provider/umat_transform.py``) and the figure/table scripts that feed the
deck (``residual-assembler/results/program1_material_validation.py``).  Both
import from here so a number in a slide and a number in a build log mean the
same thing.

Method
------
For parameter :math:`p` with base value :math:`p_0`, the reference is a
**parameter-scaled centered difference**

.. math::  F(h) = \\frac{f(p_0 + h|p_0|) - f(p_0 - h|p_0|)}{2h|p_0|}

evaluated over a ladder of *relative* steps from ``1e-3`` down to ``1e-8``.
Scaling by :math:`|p_0|` is what makes one ladder usable across parameters that
differ by ten orders of magnitude (``C11 = 168000`` next to ``gd0 = 0.001``).

Centered differences carry truncation error :math:`O(h^2)` and roundoff error
:math:`O(\\varepsilon/h)`, so accuracy improves as :math:`h` falls, bottoms out,
and then degrades.  The step is chosen by **self-consistency of the FD sequence
alone** -- the relative change between successive estimates,

.. math::  \\delta_k = \\frac{\\|F(h_k) - F(h_{k-1})\\|_\\infty}
                             {\\max(\\|F(h_k)\\|_\\infty, \\|F(h_{k-1})\\|_\\infty)}

is minimised over the ladder.  :func:`richardson` then extrapolates the chosen
pair as an independent check on that plateau.

**The selection never looks at the OTI value.**  That matters: picking the step
that happens to agree best with the quantity under test would make the check
circular, and would have hidden the real story on the FCC crystal model, whose
FD reference simply has not converged at ``1e-4``.  Only after a converged
reference is established is OTI compared against it.

Any non-finite number -- in the response, in an FD estimate, in the OTI value --
is an unconditional failure.  It is never smoothed into a ``0.0``.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

#: Relative perturbations tried for every parameter, coarse -> fine.
DEFAULT_LADDER: tuple[float, ...] = (
    1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 3e-7, 1e-7, 3e-8, 1e-8,
)

#: Denominator floor for relative measures.
TINY = 1e-30

#: A relative successive change below this between two *different* step sizes is
#: cancellation of a near-zero derivative, not convergence.  Central-difference
#: roundoff is ~eps/h, so genuine agreement between adjacent rungs bottoms out
#: far above this; a value this small means the estimates rounded identically.
_ROUNDOFF_FLOOR = 1e-12


class NonFiniteResult(ValueError):
    """A response, derivative or FD estimate contained a NaN or an infinity."""


def require_finite(label: str, values) -> np.ndarray:
    """Return ``values`` as an array, or raise :class:`NonFiniteResult`."""
    arr = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(arr)):
        bad = int(np.count_nonzero(~np.isfinite(arr)))
        raise NonFiniteResult("%s contains %d non-finite value(s) of %d" % (label, bad, arr.size))
    return arr


def relative_discrepancy(candidate, reference) -> float:
    """Max-norm *discrepancy* of *candidate* against the scale of *reference*.

    Deliberately not called "relative error": the reference here is a
    finite-difference estimate, not an analytical truth, so this is the relative
    disagreement between two approximations of the same derivative -- an
    OTI-vs-FD discrepancy -- bounded below by the FD reference's own truncation
    uncertainty (see :attr:`FDStudy.convergence_metric`). "error" is reserved for
    comparison against an analytical derivative, which this framework does not
    have.
    """
    c = np.asarray(candidate, dtype=float)
    r = np.asarray(reference, dtype=float)
    if not (np.all(np.isfinite(c)) and np.all(np.isfinite(r))):
        return float("nan")
    return float(np.max(np.abs(c - r)) / max(float(np.max(np.abs(r))), TINY))


#: Back-compatible alias. Prefer :func:`relative_discrepancy`.
relative_error = relative_discrepancy


def richardson(coarse, fine, ratio: float):
    """Richardson extrapolation of two centered differences.

    With error :math:`O(h^2)` and ``ratio = h_coarse / h_fine``:
    :math:`F = (r^2 F_\\text{fine} - F_\\text{coarse}) / (r^2 - 1)`.
    """
    r2 = float(ratio) ** 2
    return (r2 * np.asarray(fine, dtype=float) - np.asarray(coarse, dtype=float)) / (r2 - 1.0)


@dataclass
class FDStep:
    """One rung of the ladder."""

    h_rel: float
    h_abs: float
    estimate: np.ndarray | None
    finite: bool
    delta_prev: float | None = None       # self-consistency vs the coarser rung

    def to_json(self) -> dict:
        return {
            "h_rel": self.h_rel,
            "h_abs": self.h_abs,
            "finite": self.finite,
            "delta_vs_previous": self.delta_prev,
        }


@dataclass
class FDStudy:
    """The FD convergence history for one (parameter, response-block) pair."""

    block: str
    parameter: str
    props_index: int
    base_value: float
    steps: list[FDStep] = field(default_factory=list)
    selected: int | None = None
    converged: bool = False
    convergence_delta: float = float("nan")
    richardson_gap: float | None = None
    note: str = ""

    @property
    def convergence_metric(self) -> tuple[float, str]:
        """How converged the *selected* estimate is, and how that was judged.

        The successive change :math:`\\delta` locates the plateau but overstates
        the residual error: it is the gap between two rungs, not the gap between
        the chosen rung and the limit.  When Richardson extrapolation of the
        chosen pair is available its gap is the sharper -- and standard --
        estimate of that residual, so it is what certifies convergence.
        """
        if self.richardson_gap is not None and math.isfinite(self.richardson_gap):
            return self.richardson_gap, "richardson"
        return self.convergence_delta, "successive"

    @property
    def selected_step(self) -> FDStep | None:
        return None if self.selected is None else self.steps[self.selected]

    @property
    def reference(self) -> np.ndarray | None:
        step = self.selected_step
        return None if step is None else step.estimate

    @property
    def selected_h_rel(self) -> float | None:
        step = self.selected_step
        return None if step is None else step.h_rel

    def to_json(self) -> dict:
        metric, basis = self.convergence_metric
        return {
            "block": self.block,
            "parameter": self.parameter,
            "props_index": self.props_index,
            "base_value": self.base_value,
            "h_tested": [s.h_rel for s in self.steps],
            "fd_self_consistency": [s.delta_prev for s in self.steps],
            "selected_h_rel": self.selected_h_rel,
            "converged": self.converged,
            "convergence_delta": self.convergence_delta,
            "richardson_gap": self.richardson_gap,
            "convergence_metric": metric,
            "convergence_basis": basis,
            "note": self.note,
        }


def _norm(a) -> float:
    return float(np.max(np.abs(a))) if a is not None and a.size else 0.0


def parameter_studies(
    evaluate,
    props,
    props_index: int,
    *,
    parameter: str = "",
    ladder: tuple[float, ...] = DEFAULT_LADDER,
    scale: float | None = None,
) -> dict[str, FDStudy]:
    """Run the ladder for one parameter and build a study per response block.

    ``evaluate(props) -> dict[str, ndarray]`` is called twice per rung; every
    block it returns is differenced from the *same* pair of evaluations, so
    adding blocks (stress, state, ...) costs nothing extra.

    Selection happens here; the *convergence verdict* is applied separately by
    :func:`apply_convergence_gate`, because the threshold belongs to the
    tolerance the block is being certified against.

    ``scale`` overrides the perturbation scale.  It defaults to :math:`|p_0|`,
    which is what makes one relative ladder work across parameters spanning
    orders of magnitude.  A caller differencing a *vector* whose components may
    legitimately be zero (a strain increment, say) passes the vector's own scale
    instead, so a zero component is not perturbed by a full unit.
    """
    props = list(props)
    j = int(props_index) - 1
    p0 = float(props[j])
    scale = float(scale) if scale else (abs(p0) or 1.0)

    per_block: dict[str, list[FDStep]] = {}
    for h_rel in ladder:
        h_abs = h_rel * scale
        pp = list(props); pp[j] = p0 + h_abs
        pm = list(props); pm[j] = p0 - h_abs
        try:
            up = evaluate(pp)
            um = evaluate(pm)
        except NonFiniteResult:
            up = um = None

        keys = sorted(up) if up is not None else sorted(per_block) or []
        for key in keys:
            block = per_block.setdefault(key, [])
            if up is None or um is None:
                block.append(FDStep(h_rel, h_abs, None, False))
                continue
            a = np.asarray(up[key], dtype=float)
            b = np.asarray(um[key], dtype=float)
            est = (a - b) / (2.0 * h_abs)
            ok = bool(np.all(np.isfinite(est)))
            block.append(FDStep(h_rel, h_abs, est if ok else None, ok))

    studies: dict[str, FDStudy] = {}
    for key, steps in per_block.items():
        study = FDStudy(block=key, parameter=parameter, props_index=int(props_index), base_value=p0,
                        steps=steps)
        _select(study)
        studies[key] = study
    return studies


def apply_convergence_gate(study: FDStudy, conv_tol: float) -> None:
    """Decide whether *study*'s reference is converged enough to certify at ``conv_tol``."""
    if study.selected is None:
        study.converged = False
        return
    metric, basis = study.convergence_metric
    study.converged = bool(math.isfinite(metric) and metric <= conv_tol)
    if not study.converged:
        study.note = ("FD reference not converged: %s estimate %.2e > %.2e"
                      % (basis, metric, conv_tol))
    elif study.note.startswith("FD reference not converged"):
        study.note = ""


def _select(study: FDStudy) -> None:
    """Pick the plateau of the FD sequence. Never consults the OTI value.

    Selection is by **two-sided self-consistency**: a rung is a plateau only if
    it agrees with its coarser neighbour *and* the next finer rung agrees with
    it.  Scoring each rung by ``max(delta_here, delta_next)`` and minimising that
    is what makes the choice robust in the roundoff regime.  A single-sided
    "minimum successive change" rule can lock onto an isolated cancellation dip
    -- two fine-step estimates of a near-zero derivative that happen to differ by
    ~1e-16 -- and certify a step where roundoff, not truncation, dominates.
    Requiring the agreement to persist across two refinements rejects that.
    """
    steps = study.steps
    finite = [i for i, s in enumerate(steps) if s.finite]
    if not finite:
        study.note = "every finite-difference evaluation was non-finite"
        return

    # successive change between adjacent finite rungs (coarse -> fine)
    deltas: dict[int, float] = {}
    prev = None
    for i in finite:
        if prev is not None:
            a, b = steps[prev].estimate, steps[i].estimate
            denom = max(_norm(a), _norm(b), TINY)
            delta = float(np.max(np.abs(b - a))) / denom
            steps[i].delta_prev = delta
            deltas[i] = delta
        prev = i

    if not deltas:                           # only one usable rung
        study.selected = finite[0]
        study.convergence_delta = float("nan")
        study.note = "only one finite rung; convergence could not be established"
        return

    # two-sided score: agreement with the coarser neighbour AND the finer one.
    # deltas[i] pairs rung i with its coarser neighbour; deltas[next_i] pairs the
    # next finer rung with i.  A genuine plateau makes both small, so an isolated
    # cancellation dip (delta tiny at one rung, large at its neighbour) is
    # rejected by the max.
    #
    # The finest rung has no finer neighbour to corroborate it.  It may still be
    # the right choice -- for a stiff elastic constant the central difference is
    # genuinely still converging at 1e-8 -- so it self-corroborates, UNLESS its
    # successive change is below the roundoff floor.  A relative change at or
    # near zero between two *different* step sizes cannot come from a real
    # derivative (roundoff alone is ~eps/h); it is total cancellation of a
    # near-zero derivative, and must not be read as perfect convergence (the
    # J2-vs-nu failure mode).
    order = sorted(deltas)                    # finite rungs that have a delta
    finest = order[-1]
    best_i, best_score = None, float("inf")
    for pos, i in enumerate(order):
        if pos + 1 < len(order):
            look_ahead = deltas[order[pos + 1]]
        else:                                 # the finest rung: self-corroborate
            look_ahead = deltas[i] if deltas[i] >= _ROUNDOFF_FLOOR else float("inf")
        score = max(deltas[i], look_ahead)
        if score < best_score:
            best_i, best_score = i, score

    if best_i is None:                        # only rung was a cancelled finest
        best_i = finest
        best_score = deltas[finest]

    study.selected = best_i
    study.convergence_delta = best_score

    # independent check on the plateau: Richardson-extrapolate the chosen pair
    coarser = [i for i in finite if i < best_i]
    if coarser:
        prev_i = coarser[-1]
        ratio = steps[prev_i].h_rel / steps[best_i].h_rel
        if ratio > 1.0:
            extrap = richardson(steps[prev_i].estimate, steps[best_i].estimate, ratio)
            study.richardson_gap = relative_error(steps[best_i].estimate, extrap)
    return


@dataclass
class ParameterVerdict:
    """OTI compared against a converged FD reference, for one parameter/block.

    ``oti_vs_fd`` is an **OTI-FD relative discrepancy**, not an error against an
    analytical truth. ``within_fd_uncertainty`` is True when that discrepancy is
    no larger than the FD reference's own convergence uncertainty (its Richardson
    gap), i.e. OTI and FD agree to the limit of what the FD reference can
    resolve; in that regime one cannot say which is closer to the true
    derivative.
    """

    block: str
    parameter: str
    study: FDStudy
    oti_vs_fd: float = float("nan")
    fd_uncertainty: float = float("nan")
    within_fd_uncertainty: bool = False
    passed: bool = False
    reason: str = ""

    def to_json(self) -> dict:
        d = self.study.to_json()
        d.update({"oti_vs_fd_rel_discrepancy": self.oti_vs_fd,
                  "fd_reference_uncertainty": self.fd_uncertainty,
                  "within_fd_uncertainty": self.within_fd_uncertainty,
                  "status": "PASS" if self.passed else "FAIL",
                  "reason": self.reason})
        return d


def compare(study: FDStudy, oti_value, tol: float) -> ParameterVerdict:
    """Compare an OTI derivative against an already-selected FD reference.

    Reports the OTI-FD relative discrepancy and whether it lies within the FD
    reference's own uncertainty. It never asserts OTI is *more accurate* than FD
    from a small discrepancy alone -- only that they agree to within the FD
    reference's resolution.
    """
    verdict = ParameterVerdict(block=study.block, parameter=study.parameter, study=study)
    metric, _ = study.convergence_metric
    verdict.fd_uncertainty = metric
    ref = study.reference
    if ref is None:
        verdict.reason = study.note or "no usable finite-difference reference"
        return verdict

    oti = np.asarray(oti_value, dtype=float)
    if not np.all(np.isfinite(oti)):
        verdict.oti_vs_fd = float("nan")
        verdict.reason = "OTI derivative is non-finite"
        return verdict
    if oti.shape != ref.shape:
        verdict.reason = "shape mismatch: OTI %s vs FD %s" % (oti.shape, ref.shape)
        return verdict

    verdict.oti_vs_fd = relative_discrepancy(oti, ref)
    if not math.isfinite(verdict.oti_vs_fd):
        verdict.reason = "OTI-FD comparison produced a non-finite value"
        return verdict
    # "comparable to the Richardson estimate" -> within the FD reference's own
    # uncertainty (allow a small factor, since the gap is itself an estimate).
    verdict.within_fd_uncertainty = bool(math.isfinite(metric)
                                         and verdict.oti_vs_fd <= 3.0 * metric)
    if not study.converged:
        verdict.reason = ("FD reference not converged (%s); OTI-FD discrepancy %.2e is not certifiable"
                          % (study.note, verdict.oti_vs_fd))
        return verdict
    if verdict.oti_vs_fd > tol:
        verdict.reason = ("OTI-FD discrepancy %.2e exceeds tolerance %.2e at converged h=%.0e"
                          % (verdict.oti_vs_fd, tol, study.selected_h_rel))
        return verdict

    verdict.passed = True
    if verdict.within_fd_uncertainty:
        verdict.reason = ("OTI agrees with FD within the FD reference's uncertainty "
                          "(discrepancy %.2e <= ~3x Richardson gap %.2e at h=%.0e)"
                          % (verdict.oti_vs_fd, metric, study.selected_h_rel))
    else:
        verdict.reason = ("converged FD at h=%.0e (Richardson gap %.2e); OTI-FD discrepancy %.2e "
                          "<= tol %.0e" % (study.selected_h_rel, metric, verdict.oti_vs_fd, tol))
    return verdict


def verify_derivatives(
    evaluate,
    props,
    parameters,
    oti_values,
    *,
    tolerances,
    ladder: tuple[float, ...] = DEFAULT_LADDER,
    conv_tol=None,
) -> dict:
    """Full OTI-vs-FD verification over every parameter and response block.

    Parameters
    ----------
    evaluate:
        ``evaluate(props) -> {block: ndarray}``; the *response* whose parameter
        derivative is under test (stress, state, ...).
    props:
        the base property vector -- the physical operating point.
    parameters:
        ``[{"name": ..., "props_index": ...}, ...]``, in OTI direction order.
    oti_values:
        ``{block: ndarray}`` with the parameter as the **last** axis, so
        ``oti_values[block][..., c]`` matches ``evaluate``'s block shape.
    tolerances:
        ``{block: tol}`` -- the OTI-vs-FD tolerance per block. Unchanged from
        the values Program 1 has always used.
    conv_tol:
        ``{block: tol}`` (or a scalar) for FD self-consistency.  Defaults to a
        tenth of the comparison tolerance: the reference has to be an order of
        magnitude better than what it certifies.
    """
    if conv_tol is None:
        conv_tol = {k: 0.1 * v for k, v in tolerances.items()}
    elif not isinstance(conv_tol, dict):
        conv_tol = {k: float(conv_tol) for k in tolerances}

    verdicts: list[ParameterVerdict] = []
    for c, spec in enumerate(parameters):
        name = spec["name"]
        studies = parameter_studies(evaluate, props, int(spec["props_index"]),
                                    parameter=name, ladder=ladder)
        for block, tol in tolerances.items():
            study = studies.get(block)
            if study is None:
                verdicts.append(ParameterVerdict(
                    block, name, FDStudy(block, name, int(spec["props_index"]), 0.0),
                    reason="response block %r was never produced" % block))
                continue
            apply_convergence_gate(study, conv_tol[block])
            verdicts.append(compare(study, np.asarray(oti_values[block])[..., c], tol))

    worst = {b: 0.0 for b in tolerances}
    for v in verdicts:
        e = v.oti_vs_fd
        if not math.isfinite(e):
            worst[v.block] = float("nan")
        elif math.isfinite(worst[v.block]):
            worst[v.block] = max(worst[v.block], e)

    return {
        "method": "parameter-scaled centered differences, plateau-selected, Richardson-checked",
        "ladder": list(ladder),
        "tolerances": dict(tolerances),
        "convergence_tolerances": dict(conv_tol),
        "worst_rel": worst,
        "passed": all(v.passed for v in verdicts) and bool(verdicts),
        "parameters": [v.to_json() for v in verdicts],
    }


def format_report(report: dict, *, indent: str = " ") -> str:
    """The per-parameter evidence table."""
    rows = report["parameters"]
    if not rows:
        return indent + "(no parameters verified)"
    head = ("%-8s %-7s %13s %8s %10s %13s  %s"
            % ("param", "block", "h tested", "sel. h", "FD unc.", "OTI-FD discr.", "status"))
    lines = [indent + head, indent + "-" * len(head)]
    for r in rows:
        hs = r.get("h_tested") or []
        span = "%.0e..%.0e" % (hs[0], hs[-1]) if hs else "-"
        sel = "%.0e" % r["selected_h_rel"] if r.get("selected_h_rel") else "-"
        metric = r.get("convergence_metric")
        conv_s = ("%.1e%s" % (metric, "R" if r.get("convergence_basis") == "richardson" else "d")
                  if metric is not None and math.isfinite(metric) else "-")
        err = r.get("oti_vs_fd_rel_discrepancy", r.get("oti_vs_fd_rel"))
        within = r.get("within_fd_uncertainty")
        err_s = ("%.2e%s" % (err, "*" if within else " ")
                 if err is not None and math.isfinite(err) else "NON-FINITE")
        lines.append(indent + "%-8s %-7s %13s %8s %10s %13s  %s"
                     % (r["parameter"][:8], r["block"][:7], span, sel, conv_s, err_s, r["status"]))
        if r["status"] != "PASS":
            lines.append(indent + "    reason: %s" % r["reason"])
    lines.append(indent + "(FD unc. = FD reference uncertainty: Richardson gap [R] or successive "
                          "change [d] at the selected step;")
    lines.append(indent + " OTI-FD discr. = OTI-vs-FD relative discrepancy; * = within the FD "
                          "reference's own uncertainty)")
    return "\n".join(lines)
