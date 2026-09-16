"""Residual derivatives of a verified fixture, against finite differences.

``python -m residual_core.core.fixture_residual_check --fixture F --out R``

The check the interface's "Check residual derivatives" button runs, and the
report its screen draws. For one frozen fixture it answers, with the numbers:

**May this fixture be used at all?** Only when every one of the six
acceptance gates it carries is explicitly ``True`` and it is not a named
attribution subject (:data:`ATTRIBUTION_SUBJECTS`). Anything else is refused
before a single residual is assembled, and the refusal says which gate.

**dR/du** -- the element residual of a C3D8 unit cube (the element the
verification ran on), ``R(u) = sum_k w_k detJ_k B_k^T sigma_k(u)``, and its
Jacobian ``K = sum_k w_k detJ_k B_k^T D_k B_k`` built from the DDSDDE the
converted build reported. Differenced column by column over a sweep of step
sizes. Where this repository has its own model of the material
(:mod:`residual_core.materials.j2_radial_return` for the bundled J2) the
stress inside ``R`` is RECOMPUTED from the fixture's previous state and the
displaced strain by an independent implementation, so the difference is a
true derivative of the residual -- including at the increment the material
yields in, where a secant of the recorded history is not a derivative of
anything. Where it does not, the recorded linearization ``sigma_n + D_n B du``
is differenced instead, and the report names that as an assembly-consistency
check rather than a constitutive one.

**dR/dq** and **dR/dp** -- where the model is available: the analytic side is
``d sigma / d q`` (closed form) and ``d sigma / d p`` (forward-mode dual
numbers) through this repository's own assembly (``element_dR_dq``,
``integrate_dRe_dp_single``); the difference re-evaluates the residual at a
perturbed state or parameter with the independent model. Where no model is
available these are reported NOT ESTABLISHED, with what would establish them.

**Which layer?** :func:`~residual_core.core.where_it_went_wrong.diagnose` is
run on the same fixture and its findings are carried in the report, so a
failure above arrives with the first layer that did not hold.

Unknown is never a pass: a check that could not run has status
``not_established`` and the summary counts it apart from the ones that held.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Callable, Optional, Sequence

import numpy as np

from residual_core.core.finite_difference import sweep_steps
from residual_core.formulations.c3d8_kernel import (
    ABAQUS_C3D8_GAUSS, b_matrix_reference, element_internal_force_small_strain,
    element_tangent)

__all__ = ["ATTRIBUTION_SUBJECTS", "BUNDLED_J2_SHA256", "STEPS", "eligibility",
           "check_fixture", "main"]

#: Committed fixtures kept as the SUBJECT of an attribution test and never as
#: something to assemble from. The viscoelastic case is a withdrawn
#: primal_mismatch_explained entry: its two builds disagreed.
ATTRIBUTION_SUBJECTS = frozenset({"umat_viscoelastic--c6ae96a734.json"})

#: The bytes of ``UMATs/UMATs/generic_ps/j2_props.f``. Identity by content,
#: not by name: a file called j2_props.f with other bytes is another material.
BUNDLED_J2_SHA256 = ("4362e759c4642cb2f74e20c50b86a9472b1e80d6a7fef663bf8faa3d"
                     "fa8de1c8")

STEPS = (1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8)
TOLERANCE = 1e-6

UNIT_CUBE = np.array([[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.],
                      [0., 0., 1.], [1., 0., 1.], [1., 1., 1.], [0., 1., 1.]])

HOLDS, FAILS, NOT_ESTABLISHED = "holds", "fails", "not_established"


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def eligibility(fixture, path: Optional[Path] = None) -> dict:
    """Whether this fixture may be assembled from, and exactly why not."""
    name = Path(path or getattr(fixture, "path", "") or "").name
    not_true = list(getattr(fixture, "gates_not_true", ()))
    reasons = []
    if name in ATTRIBUTION_SUBJECTS:
        reasons.append(f"{name} is kept only as an attribution subject: it "
                       f"comes from a withdrawn primal_mismatch_explained "
                       f"entry whose two builds disagreed")
    for gate in not_true:
        value = getattr(fixture, "evidence", {}).get(gate, "absent")
        reasons.append(f"acceptance gate {gate} is "
                       f"{'not established' if value in (None, 'absent') else json.dumps(value)}"
                       f", not true")
    return {"eligible": not reasons, "reasons": reasons,
            "gates": {g: getattr(fixture, "evidence", {}).get(g, "absent")
                      for g in _gates()}}


def _gates() -> tuple:
    from residual_core.materials.verified_fixture import EVIDENCE_GATES
    return EVIDENCE_GATES


def _homogeneous(strain: Sequence[float]) -> np.ndarray:
    """Nodal displacements of the unit cube producing ``strain`` everywhere."""
    e = np.asarray(strain, dtype=float)
    grad = np.array([[e[0], e[3] / 2.0, e[4] / 2.0],
                     [e[3] / 2.0, e[1], e[5] / 2.0],
                     [e[4] / 2.0, e[5] / 2.0, e[2]]])
    return (UNIT_CUBE @ grad.T).reshape(-1)


def _points():
    out = []
    for point, weight in zip(ABAQUS_C3D8_GAUSS.points,
                             ABAQUS_C3D8_GAUSS.weights):
        B, detJ = b_matrix_reference(UNIT_CUBE, point)
        out.append((B, detJ * weight))
    return out


def _is_bundled_j2(fixture) -> bool:
    return (getattr(fixture, "source_sha256", "") == BUNDLED_J2_SHA256
            and len(fixture.props) == 4 and fixture.ntens == 6
            and not fixture.finite_strain)


def _reference_update():
    """The OTHER implementation, from the UMAT repository, or None."""
    try:
        from umat_oti.validation.j2_reference import (J2Parameters, J2State,
                                                      integrate_increment)
    except Exception:                               # noqa: BLE001
        return None

    def evaluate(props, stress_n, eqplas_n, dstrain):
        params = J2Parameters(E=float(props[0]), nu=float(props[1]),
                              SIGY0=float(props[2]), H=float(props[3]))
        state = J2State(stress=tuple(float(s) for s in stress_n),
                        statev=(float(eqplas_n),))
        result = integrate_increment(params, state,
                                     tuple(float(e) for e in dstrain))
        return np.asarray(result.stress, dtype=float)

    return evaluate


def _record(name: str, increment, status: str, detail: str, sweep=None,
            **measured) -> dict:
    record = {"name": name, "increment": increment, "status": status,
              "detail": detail, **measured}
    record.setdefault("step", measured.get("step"))
    if sweep is not None:
        record["sweep"] = [{"step": s, "relative": e}
                           for s, e in zip(sweep.steps, sweep.errors)]
        record["best"] = sweep.best
        record["plateau"] = sweep.plateau
        record["plateau_span"] = list(sweep.plateau_span)
    return record


def _verdict(sweep) -> str:
    return HOLDS if sweep.converged else FAILS


def _check_increment(fixture, index: int, model: Optional[Callable]) -> list:
    """dR/du, dR/dq and dR/dp at converted record ``index`` (>= 1)."""
    previous, current = fixture.converted[index - 1], fixture.converted[index]
    increment = (f"step {current.step} increment {current.increment}"
                 if current.step is not None else current.increment)
    checks = []
    points = _points()
    base_u = _homogeneous(current.dstrain)
    D = np.asarray(current.tangent, dtype=float)
    K = element_tangent(UNIT_CUBE, None, np.broadcast_to(D, (8, 6, 6)),
                        mode="small")
    q_n = float(previous.state[0]) if previous.state.size else 0.0
    scale_u = max(float(np.max(np.abs(base_u))), 1e-12)

    if model is not None:
        props = [float(p) for p in fixture.props]
        stress_prev = np.asarray(previous.stress, dtype=float)

        def residual_u(u):
            # Each integration point is updated from ITS OWN strain, B_k u. A
            # perturbation of one nodal degree of freedom strains the points
            # differently; taking one point's strain for all eight was a flat
            # 0.82 disagreement in the first version of this check.
            stress = np.array([model(props, stress_prev, q_n, B @ u)
                               for B, _w in points])
            return element_internal_force_small_strain(UNIT_CUBE, stress)

        recomputed = model(props, stress_prev, q_n, current.dstrain)
        primal = float(np.max(np.abs(recomputed - current.stress))
                       / max(float(np.max(np.abs(current.stress))), 1e-30))
        checks.append(_record(
            "stress recomputed by the independent model", increment,
            HOLDS if primal <= 1e-8 else FAILS,
            f"umat_oti.validation.j2_reference from the fixture's previous "
            f"stress and EQPLAS and this increment's strain reproduces the "
            f"converted build's stress to {primal:.3e}", relative=primal))
        what = ("stress recomputed at every perturbed displacement by the "
                "independent J2 model, so this is a derivative of the "
                "residual and not a secant of the recorded history")
    else:
        stress_n = np.asarray(current.stress, dtype=float)

        def residual_u(u):
            stress = np.array([stress_n + D @ (B @ (u - base_u))
                               for B, _w in points])
            return element_internal_force_small_strain(UNIT_CUBE, stress)
        what = ("no model of this material is carried here, so the recorded "
                "linearization sigma_n + D_n B du is differenced: this checks "
                "shape functions, quadrature, B, Voigt mapping and the element "
                "Jacobian against the element residual, not the constitutive "
                "tangent")

    def difference_u(h):
        columns = []
        for dof in range(24):
            plus, minus = base_u.copy(), base_u.copy()
            plus[dof] += h
            minus[dof] -= h
            columns.append((residual_u(plus) - residual_u(minus)) / (2.0 * h))
        return np.array(columns).T

    sweep = sweep_steps(K, difference_u, steps=STEPS, tolerance=TOLERANCE,
                        scale=scale_u)
    checks.append(_record(
        "dR/du", increment, _verdict(sweep),
        f"K from the converted build's DDSDDE against a centred difference "
        f"of the element residual: {sweep.verdict()}; {what}", sweep,
        analytic_source="converted build's DDSDDE (Abaqus, OTI)",
        worst_components=sweep.worst_components()))

    if model is None:
        for name in ("dR/dq", "dR/dp"):
            checks.append(_record(
                name, increment, NOT_ESTABLISHED,
                "this fixture carries no d sigma/d" + name[-1] + " and this "
                "repository has no model of the material to take one from",
                would_establish="a fixture carrying the transformed UMAT's "
                                "own DSIGMA_DQ / DSIGMA_DP, or a model of "
                                "this material"))
        return checks

    from residual_core.core.state_sensitivity import element_dR_dq
    from residual_core.core.umat_oti_sensitivity import integrate_dRe_dp_single
    from residual_core.materials import j2_radial_return as j2

    strain = current.dstrain
    dq = j2.dsigma_dq(props, stress_prev, q_n, strain)
    per_point = np.tile(dq.reshape(1, 6, 1), (8, 1, 1))
    dRdq_columns = element_dR_dq(UNIT_CUBE, per_point)
    analytic_q = dRdq_columns.sum(axis=1)       # same q moved at every point

    def residual_q(q):
        stress = model(props, stress_prev, q, strain)
        return element_internal_force_small_strain(UNIT_CUBE,
                                                   np.tile(stress, (8, 1)))

    if not np.any(dq):
        checks.append(_record(
            "dR/dq", increment, NOT_ESTABLISHED,
            "d sigma/d EQPLAS is exactly zero here (the trial state is "
            "elastic), so a comparison could not fail and is not reported as "
            "a pass", would_establish="an increment whose trial state yields"))
    else:
        sweep = sweep_steps(
            analytic_q, lambda h: (residual_q(q_n + h) - residual_q(q_n - h))
            / (2.0 * h), steps=STEPS, tolerance=TOLERANCE,
            scale=max(abs(q_n), 1e-4))
        checks.append(_record(
            "dR/dq", increment, _verdict(sweep),
            f"sum over points of element_dR_dq with closed-form d sigma/d "
            f"EQPLAS against a centred difference of the residual in EQPLAS "
            f"at every point: {sweep.verdict()}", sweep,
            analytic_source="closed-form radial return in this repository",
            eqplas_n=q_n))

    dp = j2.dsigma_dp(props, stress_prev, q_n, strain)
    for k, name in enumerate(j2.PARAMETERS):
        analytic_p = integrate_dRe_dp_single(
            [B for B, _w in points], [w for _B, w in points],
            [dp[:, k] for _ in points])

        def residual_p(value, k=k):
            perturbed = list(props)
            perturbed[k] = value
            stress = model(perturbed, stress_prev, q_n, strain)
            return element_internal_force_small_strain(
                UNIT_CUBE, np.tile(stress, (8, 1)))

        if not np.any(analytic_p):
            checks.append(_record(
                f"dR/d{name}", increment, NOT_ESTABLISHED,
                f"d sigma/d {name} is exactly zero at this increment, so a "
                f"comparison could not fail",
                would_establish=f"an increment where the stress depends on "
                                f"{name}"))
            continue
        sweep = sweep_steps(
            analytic_p, lambda h: (residual_p(props[k] + h)
                                   - residual_p(props[k] - h)) / (2.0 * h),
            steps=STEPS, tolerance=TOLERANCE, scale=abs(props[k]))
        checks.append(_record(
            f"dR/d{name}", increment, _verdict(sweep),
            f"integrate_dRe_dp_single with d sigma/d {name} by dual numbers "
            f"against a centred difference of the residual with the "
            f"independent model: {sweep.verdict()}", sweep,
            analytic_source="forward-mode dual numbers through this "
                            "repository's radial return"))
    return checks


def _increments_to_check(fixture) -> list:
    """The first increment the state leaves zero in, and the last one.

    Chosen by the state, not by the size of a tangent change: the lead's
    record of a failed filter shows a smoothly nonlinear material changes its
    tangent every increment. The yield increment is where a secant check
    fails and a derivative of the residual must still hold.
    """
    records = fixture.converted
    chosen = []
    for index in range(1, len(records)):
        before, after = records[index - 1].state, records[index].state
        if before.size and np.all(before == 0.0) and np.any(after != 0.0):
            chosen.append(index)
            break
    # The first increment of every later step: a reversal or a hold starts a
    # new regime, and a derivative that holds only on the first ramp has not
    # been shown to hold on the path.
    for index in range(1, len(records)):
        if (records[index].step is not None
                and records[index].step != records[index - 1].step
                and index not in chosen):
            chosen.append(index)
    last = len(records) - 1
    if last >= 1 and last not in chosen:
        chosen.append(last)
    return sorted(chosen)


def _independent_force(stress_at_ip) -> np.ndarray:
    """B^T sigma summed point by point, written out rather than called."""
    out = np.zeros(24)
    for (B, weight), stress in zip(_points(), stress_at_ip):
        out = out + B.T @ np.asarray(stress, dtype=float) * weight
    return out


def _diagnosis(fixture):
    """The layer diagnosis, with every probe this fixture makes possible."""
    from residual_core.core.where_it_went_wrong import (
        ParameterSensitivityProbe, StateSensitivityProbe, diagnose)

    kwargs = {"assemble": lambda s: element_internal_force_small_strain(
                  UNIT_CUBE, s),
              "reference": _independent_force}
    if _is_bundled_j2(fixture) and len(fixture.converted) >= 2:
        from residual_core.core.state_sensitivity import element_dR_dq
        from residual_core.core.umat_oti_sensitivity import (
            integrate_dRe_dp_single)
        from residual_core.materials import j2_radial_return as j2

        props = [float(p) for p in fixture.props]
        kwargs["material_update"] = (
            lambda stress, state, dstrain: j2.update(
                props, [float(v) for v in stress], float(state[0]),
                [float(e) for e in dstrain])[0])
        previous, current = fixture.converted[-2], fixture.converted[-1]
        q_n, strain = float(previous.state[0]), current.dstrain
        stress_prev = [float(v) for v in previous.stress]

        def residual(values, q):
            stress = np.array(j2.update(values, stress_prev, q, strain)[0],
                              dtype=float)
            return element_internal_force_small_strain(
                UNIT_CUBE, np.tile(stress, (8, 1)))

        dq = j2.dsigma_dq(props, stress_prev, q_n, strain)
        if np.any(dq):
            kwargs["state_probe"] = StateSensitivityProbe(
                analytic=element_dR_dq(
                    UNIT_CUBE, np.tile(dq.reshape(1, 6, 1), (8, 1, 1))
                ).sum(axis=1),
                difference_at=lambda h: (residual(props, q_n + h)
                                         - residual(props, q_n - h)) / (2 * h),
                scale=max(abs(q_n), 1e-4), label="dR/dEQPLAS",
                tolerance=TOLERANCE)
        dp = j2.dsigma_dp(props, stress_prev, q_n, strain)
        k = j2.PARAMETERS.index("H")
        if np.any(dp[:, k]):
            def moved(h):
                plus, minus = list(props), list(props)
                plus[k] += h
                minus[k] -= h
                return (residual(plus, q_n) - residual(minus, q_n)) / (2 * h)
            kwargs["parameter_probe"] = ParameterSensitivityProbe(
                analytic=integrate_dRe_dp_single(
                    [B for B, _w in _points()], [w for _B, w in _points()],
                    [dp[:, k] for _ in range(8)]),
                difference_at=moved, scale=abs(props[k]), label="dR/dH",
                tolerance=TOLERANCE)
    return diagnose(fixture, **kwargs)


def check_fixture(path: Path, *, fingerprint: Optional[str] = None,
                  accept_fingerprint_of_producer: bool = False) -> dict:
    """The whole report for one fixture file. Never raises for a bad fixture."""
    from residual_core.materials.verified_fixture import (
        CURRENT_TRANSFORM_FINGERPRINT, FixtureError, load)

    path = Path(path)
    report = {"schema": "residual-assembler/fixture-residual-check/1",
              "fixture": str(path), "fixture_sha256": _sha256(path),
              "recorded_generation": CURRENT_TRANSFORM_FINGERPRINT,
              "checks": [], "refused": ""}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        report["refused"] = f"the fixture could not be read: {error}"
        return report
    produced_at = str(payload.get("transform_fingerprint") or "")
    report["fixture_fingerprint"] = produced_at
    accept = fingerprint or CURRENT_TRANSFORM_FINGERPRINT
    if accept_fingerprint_of_producer and produced_at:
        accept = produced_at
        report["fingerprint_note"] = (
            f"read at the fingerprint that produced it ({produced_at}); this "
            f"repository's recorded generation is "
            f"{CURRENT_TRANSFORM_FINGERPRINT}. This is a fresh fixture checked "
            f"as evidence about the build that made it, NOT a committed "
            f"regression baseline")
    try:
        fixture = load(path, fingerprint=accept)
    except FixtureError as error:
        report["refused"] = str(error)
        return report
    report["source_id"] = fixture.source_id
    report["element_type"] = fixture.element_type
    decision = eligibility(fixture, path)
    report["eligibility"] = decision
    if not decision["eligible"]:
        report["refused"] = "; ".join(decision["reasons"])
        return report
    if not (fixture.ntens == 6 and fixture.element_type.upper() == "C3D8"
            and not fixture.finite_strain):
        report["checks"].append(_record(
            "dR/du", None, NOT_ESTABLISHED,
            f"{fixture.element_type} with NTENS={fixture.ntens}, "
            f"{fixture.kinematics}: this check assembles a small-strain C3D8 "
            f"unit cube only, so it has nothing mechanically valid to say "
            f"about this fixture",
            would_establish="a small-strain C3D8 fixture, or the kernel for "
                            "this element and kinematics"))
    else:
        model = None
        if _is_bundled_j2(fixture):
            model = _reference_update()
        report["material_model"] = ("umat_oti.validation.j2_reference "
                                    "(independent) and residual_core."
                                    "materials.j2_radial_return (analytic)"
                                    if model else "none")
        for index in _increments_to_check(fixture):
            report["checks"].extend(_check_increment(fixture, index, model))

    try:
        found = _diagnosis(fixture)
        report["diagnosis"] = {
            "ok": found.ok, "complete": found.complete, "blame": found.blame,
            "not_established": found.not_established,
            "findings": [dict(f.as_dict()) for f in found.findings]}
    except Exception as error:                      # noqa: BLE001
        report["diagnosis"] = {"error": f"{type(error).__name__}: {error}"}

    held = [c for c in report["checks"] if c["status"] == HOLDS]
    failed = [c for c in report["checks"] if c["status"] == FAILS]
    unknown = [c for c in report["checks"] if c["status"] == NOT_ESTABLISHED]
    headline = (f"{len(held)} held, {len(failed)} failed, {len(unknown)} not "
                f"established")
    if failed:
        headline = ("FAILED: " + ", ".join(f"{c['name']} at increment "
                                           f"{c['increment']}" for c in failed)
                    + f" ({headline})")
    report["summary"] = {"held": len(held), "failed": len(failed),
                         "not_established": len(unknown),
                         "headline": headline,
                         "all_applicable_held": bool(held) and not failed}
    return report


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--fixture", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--accept-fingerprint-of-producer",
                        action="store_true",
                        help="read a fresh fixture at the fingerprint that "
                             "produced it, and say so in the report")
    args = parser.parse_args(argv)
    report = check_fixture(
        args.fixture,
        accept_fingerprint_of_producer=args.accept_fingerprint_of_producer)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=1, default=float),
                        encoding="utf-8")
    print(report.get("refused") or report["summary"]["headline"])
    if report.get("refused"):
        return 3
    return 0 if report["summary"]["all_applicable_held"] else 2


if __name__ == "__main__":
    sys.exit(main())
