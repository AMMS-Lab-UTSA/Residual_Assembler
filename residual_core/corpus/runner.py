"""Run one corpus UMAT through the residual / sensitivity / consistency checks.

Produces JSONL records, one per (UMAT, feature, problem), whose fields match
the manifest feature cells ``residual_sens`` and ``global_sens``::

    status, reference, max_abs, max_rel, tolerance, what, wrt, held_fixed, evidence

plus an ``assembly_consistency`` record per problem (patch test, reaction
equilibrium, K against a difference of R, Newton rate, OTI/original primal
parity). Evidence -- every step of every ladder -- goes to a JSON file whose
path the record names.
"""
from __future__ import annotations

import json
import os
import time
from dataclasses import asdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

from . import mesh as M
from .engine import KINEMATICS, Assembly, Driver, NewtonFailed, default_kinematics, run_history
from .provider import (CorpusProvider, MaterialCallError, ProviderBuildFailed,
                       build_provider_for)
from .sources import CorpusCase, load_case, paths
from .verify import EPS32, EPS64, PLATEAU, STEPS, adjudicate, summarise
from . import hidden_state
from . import parameters as P

__all__ = ["SCHEMA", "control_case", "default_problems", "run_case", "unsupported_records"]

SCHEMA = "ra-corpus-residual/2"
RTOL_LOCAL = 1e-6
RTOL_GLOBAL = 1e-6
ATOL_FACTOR = 1e-9
RTOL_TANGENT = 1e-6
PATCH_TOL = 1e-9
NEWTON_TOL = 1e-12
OBJECTIVITY_TOL = 1e-9
TANGENT_STEPS = (1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 1e-7, 1e-8)

WHAT = {
    "residual_sens": "assembled nodal residual vector R(u_n) of the C3D8 mesh, every DOF "
                     "(free rows and the reaction rows), at every increment n of the path, "
                     "under the stated kinematics contract",
    "global_sens": "converged nodal displacements u_n at every increment, the reaction QoI "
                   "(sum of R over the loaded DOFs) and a free-node displacement QoI",
}
HELD = {
    "residual_sens": "u_n and u_{n-1} (hence DSTRAN, DFGRD0, DFGRD1, DROT), incoming "
                     "STRESS/STATEV/STRAN at every point (before the host's pre-rotation), "
                     "TIME/DTIME/KINC, geometry, COORDS (the positions handed to the routine; "
                     "see coords_contract) -- local, one increment",
    "global_sens": "prescribed displacements as functions of the load factor, load schedule "
                   "and time, geometry, COORDS (see coords_contract), initial STATEV -- total: history propagated through "
                   "STRESS, STATEV, STRAN, F0 increment to increment; equilibrium re-converged",
}
REFERENCE = {
    "residual_sens": "centred FD of R assembled from the ORIGINAL UMAT (separately compiled, no "
                     "OTI) at p+/-h, incoming state restored identically (h=0 replays "
                     "bit-exact); relative steps %s; entrywise FD-only plateau >=%d steps, "
                     "tolerance atol + rtol|D_e| + 2 u_e" % (list(STEPS), PLATEAU),
    "global_sens": "centred FD of complete nonlinear re-solves of the whole load history with "
                   "the ORIGINAL UMAT at p+/-h (Newton to 1e-12 relative free residual, or to the "
                   "residual's round-off floor confirmed by one further Newton step; "
                   "nominal re-run bit-exact); relative steps %s; entrywise FD-only plateau "
                   ">=%d steps, tolerance atol + rtol|D_e| + 2 u_e" % (list(STEPS), PLATEAU),
}


# ---------------------------------------------------------------------------
def control_case(model: str) -> CorpusCase:
    """A NON-corpus control from final-umat/parameter_sensitivity/models.

    The 44 verified corpus cases contain no rate-independent plasticity; the
    plastic regression is therefore a labelled control, never counted in the
    corpus denominators.
    """
    root = paths().umat_repo / "parameter_sensitivity" / "models" / model
    contract = json.loads((root / "contract_v2.json").read_text(encoding="utf-8"))
    dims = contract["dimensions"]
    return CorpusCase(
        key="control_" + model, source_id="final-umat/parameter_sensitivity/models/%s/%s"
        % (model, contract["source"]["main_file"]), source_path=root / contract["source"]["main_file"],
        source_form="fixed", family="control: " + model, terminal_state="control_not_corpus",
        verification_fingerprint="", kinematics="small strain", ntens=int(dims["ntens"]),
        nstatv=int(dims["nstatev"]), props=[float(v) for v in contract["validation"]["props_values"]],
        material_provenance="contract_v2.json validation.props_values", element_type="C3D8")


#: halvings of one planned increment the analytic solve may take (1/64 of it)
CUTBACKS = 6


def default_problems(case: CorpusCase, *, quick: bool = False) -> List[M.Problem]:
    amplitude = 0.1 if case.finite else 0.01
    if case.key.startswith("control_"):
        amplitude = 0.01
    single = M.uniaxial(M.brick((1, 1, 1)), amplitude)
    if quick:
        single.path = ((1.0, 2), (0.5, 1))
        problems = [single]
    else:
        problems = [single, M.clamped_shear(M.brick((2, 2, 2)), amplitude),
                    M.clamped_tension(M.brick((2, 2, 2)), amplitude)]
    total = case.extra.get("loading_period_total")
    if total:
        # the probe spans the analysis time the verification run's loading spans
        # (a routine that normalises its clock -- e.g. growth by TIME/TOTALT --
        # is driven over the history its author designed, not 3x past it)
        for problem in problems:
            problem.segment_period = float(total) / len(problem.path)
            problem.description += "; analysis time %g (the verification run's loading)" % total
    element = M.author_element(case)
    if element is not None:
        # the verification run recorded one element of the author's mesh because
        # the routine caches its point's position: COORDS come from there
        for problem in problems:
            problem.material_element = element
            problem.description += ("; COORDS: probe points mapped (trilinear) into the "
                                    "author's element (%s)" % case.extra.get("node_provenance", "")[:60])
    return problems


def _scale(p):
    return abs(p) if p != 0.0 else 1.0


def case_props_map(case: CorpusCase) -> List[dict]:
    """The source's PROPS assignments and their declared domains (cached)."""
    if "props_map" not in case.extra:
        try:
            text = Path(case.source_path).read_text(encoding="utf-8", errors="replace")
        except OSError:
            text = ""
        case.extra["props_map"] = P.props_map(text, case.source_form)
    return case.extra["props_map"]


def _ladder(case: CorpusCase, slot: int, value: float):
    """The FD steps of PROPS(slot) whose p +/- h stays inside the slot's declared
    domain; the dropped ones are tagged ``outside_parameter_domain``. A nominal
    value already outside its declared domain says the declaration does not fit
    this source: the whole ladder is kept and that is stated."""
    hs = _steps(value)
    domain = P.domain_of(case_props_map(case), slot)
    note = {"parameter_domain": None if domain is None else
            {"open_interval": [domain[0], domain[1]], "basis": domain[2]}}
    if domain is None:
        return hs, dict(note, dropped_steps=[])
    if not P.inside(value, domain):
        note["parameter_domain"]["nominal_outside"] = True
        return hs, dict(note, dropped_steps=[])
    kept = {s: h for s, h in hs.items() if P.inside(value + h, domain) and P.inside(value - h, domain)}
    dropped = [{"step": s, "reason": "outside_parameter_domain"} for s in hs if s not in kept]
    return kept, dict(note, dropped_steps=dropped)


def _steps(value):
    return {s: s * _scale(value) for s in STEPS}


def parity_tolerance(case: CorpusCase) -> float:
    """OTI-vs-original primal parity bound: 1e-12, or 64 eps32 for a source that
    computes in single precision (the lift carries those locals in double)."""
    return 64.0 * EPS32 if single_precision_declarations(case) else 1e-12


def eps_eval_for(case: CorpusCase) -> float:
    """The precision the ORIGINAL computes in, for the FD round-off model."""
    return EPS32 if single_precision_declarations(case) else EPS64


def _same_bits(a: dict, b: dict, undefined_statev=()) -> List[str]:
    """Names of the outputs that are not bit-identical; STATEV entries in
    ``undefined_statev`` (declared undefined in the original AND shown unread by
    the hidden-state probe) are left out."""
    if undefined_statev:
        a, b = hidden_state._masked(a, undefined_statev), hidden_state._masked(b, undefined_statev)
    bad = []
    for k in sorted(set(a) | set(b)):
        x, y = a.get(k), b.get(k)
        if x is None or y is None:
            continue
        x, y = np.ascontiguousarray(x), np.ascontiguousarray(y)
        if x.shape != y.shape or x.tobytes() != y.tobytes():
            bad.append(k)
    return bad


class HiddenStateTrip(RuntimeError):
    """The ORIGINAL routine's outputs depend on more than its arguments."""


# ---------------------------------------------------------------------------
def _regular_at(drv: Driver, props, run, n, schedule):
    inc = schedule[n]
    U_prev = run.U[n - 1] if n else np.zeros(drv.asm.ndof)
    return drv.regular(props, U_prev, run.U[n], run.incoming[n], inc, n + 1)


def residual_sensitivity(provider, case, problem, analytic, kinematics):
    asm = Assembly(problem, provider.finite, kinematics)
    drv = Driver(provider, asm)
    schedule = analytic.schedule or problem.schedule()
    props = np.asarray(case.props, dtype=float)
    eps = eps_eval_for(case)
    first = [_regular_at(drv, props, analytic, n, schedule) for n in range(len(analytic.U))]
    base = [asm.residual(analytic.U[n], o["stress"]) for n, o in enumerate(first)]
    vscale = [asm.residual_scale_entrywise(analytic.U[n], o["stress"]) for n, o in enumerate(first)]
    force = max(float(np.max(np.abs(v))) for v in vscale) or 1.0
    comparisons, replays = [], 0
    for j, slot in enumerate(provider.slots):
        value = props[slot - 1]
        hs, ladder = _ladder(case, slot, value)
        atol = ATOL_FACTOR * force / _scale(value)
        for n in range(len(analytic.U)):
            plus, minus = {}, {}
            for s, h in hs.items():
                pp, pm = props.copy(), props.copy()
                pp[slot - 1] += h
                pm[slot - 1] -= h
                plus[s] = asm.residual(analytic.U[n], _regular_at(drv, pp, analytic, n, schedule)["stress"])
                minus[s] = asm.residual(analytic.U[n], _regular_at(drv, pm, analytic, n, schedule)["stress"])
            # h = 0 replay after the perturbed calls: every output bit-identical
            again = _regular_at(drv, props, analytic, n, schedule)
            replays += 1
            bad = _same_bits(first[n], again, case.extra.get("undefined_statev") or ())
            if bad:
                raise HiddenStateTrip("h=0 replay of the ORIGINAL at increment %d after the P%d "
                                      "ladder differs in %s" % (n + 1, slot, bad))
            c = adjudicate(analytic.local_dR_dp[n][:, j], plus, minus, base[n], hs,
                           rtol=RTOL_LOCAL, atol=atol, value_scale=vscale[n], eps_eval=eps)
            c.update(parameter="P%d" % slot, increment=n + 1, **ladder)
            comparisons.append(c)
    return comparisons, {"force_scale": force, "h0_replays_bit_exact": replays, "eps_eval": eps}


#: largest max|V_ref - V_analytic| / max|V| at which the nominal reference is
#: taken to sit on the analytic solve's equilibrium
BRANCH_TOL = 1e-6
#: iteration matrices tried, in order, for a reference solve of the ORIGINAL
REFERENCE_MATRICES = ("auto", "exact", "ra_then_exact")
#: (iteration matrix, cut-backs) tried, in order, for the analytic (OTI) solve;
#: the first is the B2 behaviour
ANALYTIC_STRATEGIES = (("exact", 0), ("ra_then_exact", 0), ("ra_stalled_exact", None),
                       ("exact", None))
#: the reference matrix that steers the original as the analytic matrix steered OTI
#: (an analytic solve with the exact matrix keeps the B2 reference order)
_SAME_STEERING = {"ra_then_exact": "ra_then_exact", "ra_stalled_exact": "auto"}


def analytic_solve(provider, problem, props, kinematics, detail=None):
    """The OTI history with sensitivities, under the first strategy of
    ``ANALYTIC_STRATEGIES`` that converges (``None`` cut-backs = ``CUTBACKS``).
    Failed attempts are listed in ``detail['solve_attempts']``; the last error
    is raised when none converges."""
    attempts, error = [], None
    for matrix, cut in ANALYTIC_STRATEGIES:
        cut = CUTBACKS if cut is None else cut
        try:
            run = run_history(provider, problem, props, material="oti", sensitivities=True,
                              kinematics=kinematics, newton_matrix=matrix, cutbacks=cut)
        except (NewtonFailed, MaterialCallError, np.linalg.LinAlgError) as failed:
            attempts.append({"newton_matrix": matrix, "cutbacks": cut,
                             "error": "%s: %s" % (type(failed).__name__, failed)})
            error = failed
            continue
        if detail is not None:
            detail["solve_attempts"] = attempts
        return run
    if detail is not None:
        detail["solve_attempts"] = attempts
    raise error


def common_schedule(provider, problem, props, kinematics, analytic, detail=None):
    """One set of increments for the analytic solve and every reference solve.

    The analytic solve cut back; the ORIGINAL at the nominal PROPS is solved on
    those increments, itself allowed to cut any of them back (the two builds
    agree to round-off, but at a hard increment Newton may converge for one and
    not the other). If it had to, the analytic solve is repeated on the finer
    increments without cut-backs of its own, and that run is returned; when it
    does not converge there, the first analytic run is kept (its reference
    solves then fail, an unsupported global cell)."""
    first = _SAME_STEERING.get(analytic.newton_matrix_arg)
    try:
        ref = reference_solve(provider, problem, props, kinematics, analytic.schedule,
                              first=first, cutbacks=CUTBACKS)
    except (NewtonFailed, MaterialCallError, np.linalg.LinAlgError) as error:
        if detail is not None:
            detail["common_schedule"] = "reference not solved: %s: %s" % (
                type(error).__name__, error)
        return analytic
    if len(ref.schedule) == len(analytic.schedule):
        return analytic
    try:
        again = run_history(provider, problem, props, material="oti", sensitivities=True,
                            kinematics=kinematics, newton_matrix=analytic.newton_matrix_arg,
                            schedule=ref.schedule)
    except (NewtonFailed, MaterialCallError, np.linalg.LinAlgError) as error:
        if detail is not None:
            detail["common_schedule"] = "analytic not solved on the reference's %d increments: " \
                "%s: %s" % (len(ref.schedule), type(error).__name__, error)
        return analytic
    again.cutbacks = analytic.cutbacks + [dict(c, by="reference") for c in ref.cutbacks]
    if detail is not None:
        detail["common_schedule"] = "analytic re-solved on the reference's %d increments" \
            % len(ref.schedule)
    return again


def steered_analytic(provider, problem, props, kinematics, reference, detail=None):
    """The analytic solve started, increment by increment, from the ORIGINAL's
    converged solution on the ORIGINAL's increments; None if it does not
    converge there.

    Where the equilibrium is not unique, the two builds -- equal to round-off --
    may each converge only into a basin of their own (0e56ba3c: K_ff indefinite
    at increment 9). Started on the original's solution, the analytic solve
    still solves R = 0 with its OWN residual; accepted, it sits on the
    reference's equilibrium, and the comparisons there are meaningful (the
    branch gate in global_sensitivity checks that it stayed there)."""
    for matrix in ("exact", "ra_then_exact"):
        try:
            run = run_history(provider, problem, props, material="oti", sensitivities=True,
                              kinematics=kinematics, newton_matrix=matrix,
                              schedule=reference.schedule, start_from=reference.V)
        except (NewtonFailed, MaterialCallError, np.linalg.LinAlgError) as error:
            if detail is not None:
                detail.setdefault("steering_attempts", []).append(
                    {"newton_matrix": matrix, "error": "%s: %s" % (type(error).__name__, error)})
            continue
        if detail is not None:
            detail["steered_onto_reference"] = {
                "newton_matrix": matrix, "increments": len(reference.schedule),
                "reference_newton_matrix": reference.newton_matrix_arg,
                "max_V_difference": max(float(np.max(np.abs(a - b)))
                                        for a, b in zip(run.V, reference.V))}
        run.reference_steering = reference.newton_matrix_arg
        return run
    return None


def reference_solve(provider, problem, props, kinematics, schedule, first=None, **kw):
    """Solve the ORIGINAL on exactly ``schedule``. The iteration matrix only
    steers Newton (the converged state solves R_original = 0 whichever matrix
    led there), so when the original's own DDSDDE matrix does not converge the
    exact and the mixed matrices are tried; ``first`` is tried before them (the
    steering of the analytic solve). The matrix used is recorded."""
    error = None
    order = ([first] if first else []) + [m for m in REFERENCE_MATRICES if m != first]
    for matrix in order:
        try:
            return run_history(provider, problem, props, material="regular",
                               kinematics=kinematics, schedule=schedule, newton_matrix=matrix,
                               **kw)
        except (NewtonFailed, MaterialCallError, np.linalg.LinAlgError) as failed:
            error = error or failed
            if isinstance(failed, MaterialCallError) and "PNEWDT" not in str(failed) \
                    and "inverted" not in str(failed):
                raise                       # the routine itself cannot be evaluated
    raise error


def _state_rel(a, b) -> float:
    """max relative difference of two IncrementStates (stress, state, strain)."""
    worst = 0.0
    for x, y in ((a.stress, b.stress), (a.state, b.state), (a.stran, b.stran)):
        x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
        if x.size:
            worst = max(worst, float(np.max(np.abs(x - y))) / max(float(np.max(np.abs(x))), 1e-300))
    return worst


def branch_departures(case, analytic, base, branch) -> set:
    """0-based increments whose comparisons are not judged because the
    reference left the analytic equilibrium (Vera R-1).

    Applies only when primal parity holds on this problem: if the OTI build's
    stress differs from the original's beyond the parity tolerance, a different
    solution may be a primal defect, and the comparisons stay judged (and fail).
    From the first departure on, every increment is excused -- the history
    differs -- unless both its solution and its incoming state (stress, STATEV,
    strain) match again to BRANCH_TOL."""
    parity = max(analytic.primal_parity) if getattr(analytic, "primal_parity", None) else 0.0
    if parity > parity_tolerance(case):
        return set()
    departed = [n for n, d in enumerate(branch) if d > BRANCH_TOL]
    if not departed:
        return set()
    out = set()
    for n in range(departed[0], len(branch)):
        rejoined = (branch[n] <= BRANCH_TOL
                    and _state_rel(analytic.incoming[n], base.incoming[n]) <= BRANCH_TOL)
        if not rejoined:
            out.add(n)
    return out


def global_sensitivity(provider, case, problem, analytic, kinematics):
    props = np.asarray(case.props, dtype=float)
    eps = eps_eval_for(case)
    # every reference solve uses the increments the analytic solve used (cut-backs
    # included), so the histories compare increment by increment
    schedule = analytic.schedule or None
    first = getattr(analytic, "reference_steering", None) or \
        _SAME_STEERING.get(analytic.newton_matrix_arg)
    base = reference_solve(provider, problem, props, kinematics, schedule, first=first,
                           keep_outputs=True)
    V0 = np.array(base.V)
    Q0 = base.qoi(problem)
    force = max(base.force_scale) or 1.0
    uscale = float(np.max(np.abs(V0))) or 1.0
    # the reference must sit on the equilibrium the analytic solve found: where
    # the equilibrium is not unique (K_ff indefinite past a bifurcation), the
    # original's Newton can land on another branch, and a difference around it
    # says nothing about the analytic derivative (B9: SeaShell, |dV| = 0.66)
    branch = [float(np.max(np.abs(V0[n] - analytic.V[n]))) / uscale for n in range(len(V0))]
    elsewhere = branch_departures(case, analytic, base, branch)
    dQ = analytic.qoi_dp(problem)                      # (ninc, nq, npar)
    comparisons, matrices, failures = [], set(), []
    achieved = [max(h[-1] for h in base.newton)]       # final relative free residuals
    for j, slot in enumerate(provider.slots):
        value = props[slot - 1]
        hs, ladder = _ladder(case, slot, value)
        Vp, Vm, Qp, Qm, levels = {}, {}, {}, {}, {}
        for s, h in hs.items():
            pp, pm = props.copy(), props.copy()
            pp[slot - 1] += h
            pm[slot - 1] -= h
            try:
                rp = reference_solve(provider, problem, pp, kinematics, schedule, first=first)
                rm = reference_solve(provider, problem, pm, kinematics, schedule, first=first)
            except (NewtonFailed, MaterialCallError, np.linalg.LinAlgError) as error:
                # the ORIGINAL cannot be re-solved at p +/- h (e.g. a viscosity of 0
                # perturbed negative): that step is missing from the ladder
                failures.append({"parameter": "P%d" % slot, "step": s,
                                 "reason": "%s: %s" % (type(error).__name__, error)})
                continue
            matrices.update(rp.newton_matrix + rm.newton_matrix)
            achieved.extend(max(h[-1] for h in r.newton) for r in (rp, rm))
            # round-off model of a converged re-solve: the free residual the two
            # solves ACHIEVED at this step and increment (relative to the force
            # scale), not the 1e-12 tolerance they sit below. A model that is
            # too small only rejects steps -- the tolerance uses the MEASURED
            # plateau spread.
            levels[s] = [max(a[-1], b[-1]) for a, b in zip(rp.newton, rm.newton)]
            Vp[s], Vm[s] = np.array(rp.V), np.array(rm.V)
            Qp[s], Qm[s] = rp.qoi(problem), rm.qoi(problem)
        if not Vp:
            comparisons.append({"verdict": "unresolved", "parameter": "P%d" % slot,
                                "value": float(value), "max_abs": None, "max_rel": None,
                                "reason": "no step of the ladder could be re-solved"})
            continue
        for n in range(len(analytic.V)):
            if n in elsewhere:
                for quantity in ["u"] + [q["name"] for q in problem.qois]:
                    comparisons.append({
                        "verdict": "unresolved", "parameter": "P%d" % slot, "increment": n + 1,
                        "quantity": quantity, "max_abs": None, "max_rel": None,
                        "reason": "reference_on_another_equilibrium: the ORIGINAL's nominal "
                                  "solve differs from the analytic solution by %.3g of |V| "
                                  "at this increment" % branch[n]})
                continue
            c = adjudicate(analytic.du_dp[n][:, j], {s: Vp[s][n] for s in Vp},
                           {s: Vm[s][n] for s in Vm}, V0[n], hs, rtol=RTOL_GLOBAL,
                           atol=ATOL_FACTOR * uscale / _scale(value),
                           value_scale=uscale, eps_eval=eps,
                           step_level={s: levels[s][n] for s in levels})
            c.update(parameter="P%d" % slot, increment=n + 1, quantity="u", **ladder)
            comparisons.append(c)
            for k, q in enumerate(problem.qois):
                qscale = force if q["kind"] == "reaction" else uscale
                c = adjudicate(np.array([dQ[n, k, j]]), {s: Qp[s][n, k:k + 1] for s in Qp},
                               {s: Qm[s][n, k:k + 1] for s in Qm}, Q0[n, k:k + 1], hs,
                               rtol=RTOL_GLOBAL, atol=ATOL_FACTOR * qscale / _scale(value),
                               value_scale=force if q["kind"] == "reaction" else uscale,
                               eps_eval=eps, step_level={s: levels[s][n] for s in levels})
                c.update(parameter="P%d" % slot, increment=n + 1, quantity=q["name"], **ladder)
                comparisons.append(c)
    # nominal re-run after every perturbed history: all outputs bit-identical
    again = run_history(provider, problem, props, material="regular", kinematics=kinematics,
                        keep_outputs=True, schedule=schedule,
                        newton_matrix=base.newton_matrix_arg)
    for n, (o1, o2) in enumerate(zip(base.outputs, again.outputs)):
        bad = _same_bits(o1, o2, case.extra.get("undefined_statev") or ()) + (["U"] if base.U[n].tobytes() != again.U[n].tobytes() else [])
        if bad:
            raise HiddenStateTrip("nominal ORIGINAL history re-run after the perturbed re-solves "
                                  "differs at increment %d in %s" % (n + 1, bad))
    return comparisons, {"force_scale": force, "u_scale": uscale, "eps_eval": eps,
                         "reference_vs_analytic_V_rel": branch,
                         "reference_on_another_equilibrium": sorted(n + 1 for n in elsewhere),
                         "resolve_residual_achieved_max": max(achieved),
                         "reference_newton_matrices": sorted(matrices),
                         "reference_resolve_failures": failures,
                         "nominal_rerun_bit_exact": True,
                         "reference_newton_iterations_max": max(len(h) for h in base.newton)}


def check_increments(problem, schedule=None) -> List[int]:
    """0-based increments where K is checked: the end of every path segment."""
    if schedule:
        return [n for n, inc in enumerate(schedule) if inc.get("segment_end")]
    out, count = [], 0
    for _, n in problem.path:
        count += n
        out.append(count - 1)
    return sorted(set(out))


def tangent_check(provider, case, problem, analytic, kinematics, increments=None):
    """K (exact, OTI) and the DDSDDE tangents against a centred FD of R(u)
    (ORIGINAL routine) at converged states, incoming state of that increment held."""
    asm = Assembly(problem, provider.finite, kinematics)
    drv = Driver(provider, asm)
    schedule = analytic.schedule or problem.schedule()
    props = np.asarray(case.props, dtype=float)
    # FD in u: no parameter-derived single-precision constant moves with u; the
    # tolerance uses the measured plateau spread, so the model only decides
    # which steps can resolve an entry, never whether it passes
    eps = EPS64
    results = []
    for n in (increments if increments is not None else check_increments(problem, schedule)):
        inc = schedule[n]
        U_prev = analytic.U[n - 1] if n else np.zeros(asm.ndof)
        U = analytic.U[n]
        incoming = analytic.incoming[n]
        local = drv.oti(props, U_prev, U, incoming, inc, n + 1)
        K = drv.exact_tangent(props, U_prev, U, incoming, inc, n + 1, local)
        reg = drv.regular(props, U_prev, U, incoming, inc, n + 1)
        K_ra = asm.tangent_ra(U, reg["stress"], reg["ddsdde"])
        K_ra_oti = asm.tangent_ra(U, local["stress"], local["ddsdde"])
        uscale = float(np.max(np.abs(U))) or 1.0

        def R(Ut):
            return asm.residual(Ut, drv.regular(props, U_prev, Ut, incoming, inc, n + 1)["stress"])

        base = R(U)
        vscale = np.repeat(asm.residual_scale_entrywise(U, reg["stress"])[:, None], asm.ndof, axis=1)
        hs = {s: s * uscale for s in TANGENT_STEPS}
        plus = {s: np.zeros((asm.ndof, asm.ndof)) for s in hs}
        minus = {s: np.zeros((asm.ndof, asm.ndof)) for s in hs}
        for s, h in hs.items():
            for k in range(asm.ndof):
                e = np.zeros(asm.ndof)
                e[k] = h
                try:
                    plus[s][:, k] = R(U + e)
                    minus[s][:, k] = R(U - e)
                except (MaterialCallError, ValueError):
                    plus[s] = minus[s] = None
                    break
        Kscale = float(np.max(np.abs(K))) or 1.0
        basem = np.repeat(base[:, None], asm.ndof, axis=1)
        out = {"increment": n + 1}
        for name, matrix in (("K_exact_oti", K), ("K_ddsdde_author", K_ra),
                             ("K_ddsdde_oti", K_ra_oti)):
            c = adjudicate(matrix, plus, minus, basem, hs, rtol=RTOL_TANGENT,
                           atol=1e-12 * Kscale, value_scale=vscale, eps_eval=eps)
            c.pop("worst", None)
            out[name] = c
        results.append(out)
    return results


def newton_rates(run) -> dict:
    """Observed convergence order of the analytic (OTI, exact-K) Newton solves."""
    orders, worst_iterations = [], 0
    starts = getattr(run, "exact_from", None) or [0] * len(run.newton)
    for history, start in zip(run.newton, starts):
        # only the iterations taken with the exact matrix (ra_then_exact starts
        # with the DDSDDE one)
        history = history[start:]
        worst_iterations = max(worst_iterations, len(history))
        for a, b in zip(history[1:-1], history[2:]):
            if 1e-14 < b and a < 1e-2 and a > 1e-13:
                orders.append(np.log(b) / np.log(a))
    median = float(np.median(orders)) if orders else None
    ok = (median is None and worst_iterations <= 3) or (median is not None and median >= 1.8)
    return {"pass": bool(ok), "median_order": median, "orders_measured": len(orders),
            "max_iterations": worst_iterations,
            "note": "orders from successive scaled free-residual pairs r_k<1e-2 and r_k+1>1e-14; "
                    "a path solved in <=3 iterations per increment is linear-exact"}


def objectivity_check(provider, case, problem, kinematics) -> dict:
    """Superposed rigid rotation (finite strain only): the same path with a hold
    segment, rotated by Q(t) during the hold and after it, against the same path
    unrotated. Frame displacements, frame reactions, sigma' = Q sigma Q^T, and
    the total derivatives dV/dp and dQoI/dp must agree to round-off."""
    from residual_core.formulations import c3d8_nlgeom as nl
    plain, rotated = M.with_superposed_rotation(problem)
    a = run_history(provider, plain, case.props, material="oti", sensitivities=True,
                    kinematics=kinematics)
    b = run_history(provider, rotated, case.props, material="oti", sensitivities=True,
                    kinematics=kinematics)
    fs = max(a.force_scale) or 1.0
    us = float(np.max(np.abs(np.array(a.V)))) or 1.0
    sig = 0.0
    for n, (oa, ob) in enumerate(zip(a.outgoing, b.outgoing)):
        Q = rotated.rotations[n]
        Sa = nl.stress_voigt_to_tensor(oa.stress)
        Sb = nl.stress_voigt_to_tensor(ob.stress)
        sig = max(sig, float(np.max(np.abs(Q @ Sa @ Q.T - Sb))) / max(float(np.max(np.abs(Sa))), 1e-300))
    mags = max(float(np.max(np.abs(np.linalg.norm(x.reshape(-1, 3), axis=1)
                                   - np.linalg.norm(y.reshape(-1, 3), axis=1))))
               for x, y in zip(a.reactions, b.reactions)) / fs
    dq_a, dq_b = a.qoi_dp(plain), b.qoi_dp(rotated)
    # yardstick per parameter: max(|dV/dp_j|, |V|/|p_j|) -- a derivative that is
    # zero in exact arithmetic (displacement control, a stiffness that scales the
    # whole response) is round-off on both sides, and round-off over round-off
    # is no measure of objectivity
    pscale = np.array([_scale(case.props[slot - 1]) for slot in provider.slots])
    dv_scale = np.maximum(np.max(np.abs(np.array(a.du_dp)), axis=(0, 1)), us / pscale)
    dv = max(float(np.max(np.abs(x - y) / dv_scale)) for x, y in zip(a.du_dp, b.du_dp))
    metrics = {
        "reaction_frame_max_rel": max(float(np.max(np.abs(x - y))) for x, y in zip(a.reactions, b.reactions)) / fs,
        "reaction_magnitude_max_rel": mags,
        "V_max_rel": max(float(np.max(np.abs(x - y))) for x, y in zip(a.V, b.V)) / us,
        "sigma_QsigmaQT_max_rel": sig,
        "dQoI_dp_max_rel": float(np.max(np.abs(dq_a - dq_b))) / max(float(np.max(np.abs(dq_a))), 1e-300),
        "dV_dp_max_rel": dv,
    }
    return {"pass": bool(all(v <= OBJECTIVITY_TOL for v in metrics.values())), **metrics,
            "tolerance": OBJECTIVITY_TOL, "rotation": "axis (1,2,3), 60 deg over a 3-increment hold "
            "after the first segment, kept to the end", "increments": len(a.V)}


def kinematics_delta(provider, case, problem) -> dict:
    """Primal differences between the B1 convention and the NLGEOM contract on
    the same path (original routine): what the change of kinematics alone does."""
    out = {}
    runs = {}
    for name in ("b1_legacy", "nlgeom_full", "nlgeom_centroid", "nlgeom"):
        try:
            runs[name] = run_history(provider, problem, case.props, material="regular",
                                     kinematics=name)
        except (NewtonFailed, MaterialCallError, np.linalg.LinAlgError) as error:
            out[name] = "%s: %s" % (type(error).__name__, error)
    ref = runs.get("nlgeom")
    if ref is None:
        return out
    q0 = ref.qoi(problem)
    s0 = max(float(np.max(np.abs(o.stress))) for o in ref.outgoing) or 1.0
    for name, r in runs.items():
        if name == "nlgeom":
            continue
        q = r.qoi(problem)
        out[name + "_vs_nlgeom"] = {
            "qoi_max_rel": [float(np.max(np.abs(q[:, k] - q0[:, k])) / (np.max(np.abs(q0[:, k])) or 1.0))
                            for k in range(q.shape[1])],
            "stress_max_rel": max(float(np.max(np.abs(a.stress - b.stress)))
                                  for a, b in zip(r.outgoing, ref.outgoing)) / s0}
    return out


def patch_check(provider, case, kinematics) -> dict:
    gradient = np.array([[1.0, 0.3, -0.2], [0.1, -0.4, 0.25], [-0.15, 0.2, 0.5]])
    gradient *= 0.05 if provider.finite else 1e-3
    mesh = M.brick((3, 3, 3), distort=0.15)
    problem = M.patch_test(mesh, gradient)
    problem.material_element = M.author_element(case)
    run = run_history(provider, problem, case.props, material="oti", kinematics=kinematics)
    U = run.U[-1].reshape(-1, 3)
    exact = mesh.coords @ gradient.T
    interior = [a for a, x in enumerate(mesh.coords)
                if all(1e-9 < x[d] < mesh.size[d] - 1e-9 for d in range(3))]
    err = float(np.max(np.abs(U[interior] - exact[interior]))) / float(np.max(np.abs(exact)))
    stress = run.outgoing[-1].stress
    spread = float(np.max(np.abs(stress - stress.mean(axis=0))))
    sscale = float(np.max(np.abs(stress))) or 1.0
    passed = bool(err <= PATCH_TOL and spread / sscale <= PATCH_TOL)
    homogeneous = material_is_homogeneous(provider, case, gradient)
    return {"pass": passed if homogeneous else None,
            "applicable": homogeneous,
            "why_not_applicable": None if homogeneous else
            "the ORIGINAL returns different stress for the same deformation at different "
            "COORDS/NOEL/NPT (position-dependent material), so the exact solution is not "
            "the linear field and the patch test does not apply",
            "interior_nodes": len(interior), "interior_displacement_rel_error": err,
            "stress_spread_rel": spread / sscale, "tolerance": PATCH_TOL,
            "mesh": "3x3x3 C3D8, interior nodes moved by 0.15 h", "gradient": gradient.tolist()}


def material_is_homogeneous(provider, case, gradient) -> bool:
    """Same deformation, same incoming state, two positions/labels: same stress?"""
    F = np.tile(np.eye(3) + gradient, (2, 1, 1)) if provider.finite else np.tile(np.eye(3), (2, 1, 1))
    eps = np.array([gradient[0, 0], gradient[1, 1], gradient[2, 2], gradient[0, 1] + gradient[1, 0],
                    gradient[0, 2] + gradient[2, 0], gradient[1, 2] + gradient[2, 1]])
    state = np.zeros((2, provider.nstatv))
    coords = np.array([[0.11, 0.23, 0.37], [0.83, 0.61, 0.71]])
    element = M.author_element(case)
    if element is not None:       # the positions the probes hand this routine
        coords = M.material_coordinates(coords, M.brick((1, 1, 1)), element)
    noel, npt = np.array([1, 7]), np.array([1, 5])
    initial = case.extra.get("initial_statev") or []
    if initial:
        state[:, :len(initial)] = np.asarray(initial, dtype=float)[:provider.nstatv]
    if case.extra.get("initial_state_from_user_subroutine") and provider.has_sdvini:
        state = provider.sdvini(state, coords, noel, npt)
    try:
        out = provider.regular(case.props, np.zeros((2, provider.ntens)), state,
                               np.zeros((2, provider.ntens)), np.tile(eps, (2, 1)),
                               np.tile(np.eye(3), (2, 1, 1)), F, time=np.zeros(2), dtime=0.5,
                               coords=coords, celent=np.ones(2), noel=noel, npt=npt, kinc=1)
    except MaterialCallError:
        return True
    s = out["stress"]
    return bool(np.max(np.abs(s[0] - s[1])) <= 1e-12 * max(float(np.max(np.abs(s))), 1e-300))


def equilibrium_check(run) -> dict:
    worst = 0.0
    for R, scale in zip(run.reactions, run.force_scale):
        total = R.reshape(-1, 3).sum(axis=0)
        scale = scale or 1.0
        worst = max(worst, float(np.max(np.abs(total))) / scale)
    return {"pass": worst <= 1e-12, "max_net_force_rel": worst,
            "statement": "sum over all nodes of the assembled R (= sum of reactions, no "
                         "external load) vanishes per direction at every converged increment, "
                         "relative to max_dof sum_e |f_e|"}


# ---------------------------------------------------------------------------


def single_precision_declarations(case: CorpusCase) -> List[str]:
    """Lines declaring default-kind REAL variables (4-byte) in the source.

    A routine computing in single precision is piecewise constant in its
    parameters at ~1e-7 relative, so a finite difference of the ORIGINAL cannot
    resolve a derivative to 1e-6 there, and its primal differs from a lift
    carried in double. Reported as context, never used to change a verdict.
    """
    import re
    hits = []
    pattern = re.compile(r"^\s*real\s*(::|\s)\s*[a-z]", re.IGNORECASE)
    for number, line in enumerate(case.source_path.read_text(errors="replace").splitlines(), 1):
        if line[:1] in "cC*!" or line.lstrip().startswith("!"):
            continue
        if pattern.match(line) and not re.match(r"^\s*real\s*\*", line, re.IGNORECASE):
            hits.append("%d: %s" % (number, line.strip()))
    return hits


def build_fields(provider_record: dict) -> dict:
    """The BUILD every derivative cell of this record verified (Vera B1 E)."""
    tree = provider_record.get("umat_tree") or {}
    return {"build": {"kind": "provider", "sha256": provider_record.get("object_sha256") or "",
                      "fingerprint": provider_record.get("fingerprint") or ""},
            "build_detail": "umat_oti.provider.build_provider: ONE compile of the source "
                            "holding the PROPS-seeded OTI lift (analytic side) and the original "
                            "routine (reference side, symbol umat)",
            "object_sha256": provider_record.get("object_sha256"),
            "build_fingerprint": provider_record.get("fingerprint"),
            "source_sha256_built": provider_record.get("source_sha256_built"),
            "umat_head": tree.get("head"), "umat_src_diff_sha256": tree.get("diff_sha256")}


SCOPE = {"residual_sens": "local", "global_sens": "total"}
QUANTITY = {"residual_sens": "dR/dp: assembled nodal residual R(u_n) of the C3D8 mesh, every DOF, "
                             "every increment, differentiated w.r.t. each PROPS slot",
            "global_sens": "du_n/dp (every DOF, every increment) and dQoI/dp for the reaction and "
                           "displacement QoIs, differentiated w.r.t. each PROPS slot"}
_RA_ROOT = Path(__file__).resolve().parents[2]


def _roots():
    """Named roots for manifest locators: this repository always; the workspace
    folders only when ``CORPUS_WORKSPACE`` (and so the corpus) is configured."""
    roots = [("ra", _RA_ROOT)]
    where = paths()
    if where.workspace is not None:
        roots += [("campaign", where.workspace / "corpus_campaign"),
                  ("corpus_run", where.workspace / "corpus_run"),
                  ("umat", where.umat_repo)]
    return [(name, str(base.resolve())) for name, base in roots]


def locator(path) -> str:
    """``<root>:<relative>`` for the manifest contract (absolute if under no root)."""
    text = str(Path(path).resolve()) if path else ""
    for root, base in _roots():
        if text.startswith(base + "/"):
            return "%s:%s" % (root, text[len(base) + 1:])
    return text


def contract_fields(feature: str, summary: dict, provider_record: dict) -> dict:
    """The manifest merge contract (scout B2) for a derivative record."""
    if feature not in SCOPE:
        return {}
    status = summary.get("status")
    decided = status in ("verified", "failed")
    out = {"reference": "fd",
           "reference_kind": "fd_resolve" if feature == "global_sens" else "fd",
           "fd_steps": list(STEPS), "tolerance": 1.0,
           "max_error": summary.get("max_error_ratio") if decided else None,
           "tolerance_rule_id": "entrywise/1",
           "rtol": RTOL_LOCAL if feature == "residual_sens" else RTOL_GLOBAL,
           "min_plateau_observed": summary.get("min_plateau_observed"),
           "plateau_basis": "fd_only", "scope": SCOPE[feature],
           "quantity": QUANTITY[feature],
           "build": {"kind": "provider", "sha256": provider_record.get("object_sha256") or "",
                     "fingerprint": provider_record.get("fingerprint") or ""}}
    return out


def _record(case, problem, feature, provider_record, summary, evidence, kin, extra=None):
    rec = {"schema": SCHEMA, "key": case.key, "source_id": case.source_id,
           "family": case.family, "kinematics": case.kinematics,
           "kinematics_contract": kin.name, "kinematics_detail": kin.describe(),
           "problem": problem.name,
           "problem_description": problem.description, "feature": feature,
           "status": summary["status"], "reason_class": summary.get("reason_class"),
           "reference": REFERENCE.get(feature, ""),
           "max_abs": summary.get("max_abs"), "max_rel": summary.get("max_rel"),
           "resolution_max": summary.get("resolution_max"),
           "tolerance": summary.get("tolerance"), "what": WHAT.get(feature, ""),
           "wrt": "PROPS slots %s (all seeded)" % [p["props_index"] for p in
                                                   provider_record["parameters"]],
           "held_fixed": HELD.get(feature, ""), "evidence": locator(evidence),
           "material_path": "live OTI provider (umat_oti.provider.build_provider), "
                            "original routine from the same compile for references",
           "source_adapted": provider_record.get("source_adapted", False),
           "adaptation_notes": provider_record.get("adaptation_notes", []),
           "counts": summary.get("counts"), "comparisons": summary.get("comparisons")}
    rec.update(build_fields(provider_record))
    rec["reference_text"] = rec["reference"]
    rec.update(contract_fields(feature, summary, provider_record))
    if rec["status"] != "verified" and not rec.get("reason"):
        rec["reason"] = summary.get("reason") or (
            "%s: %s" % (summary.get("reason_class"), summary.get("counts"))
            if summary.get("reason_class") else "comparison counts %s" % (summary.get("counts"),))
    single = single_precision_declarations(case)
    if single:
        rec["single_precision_declarations"] = single[:5]
        rec["fd_roundoff_model"] = "float32 (the original computes in single precision)"
    if feature in ("residual_sens", "global_sens"):
        rec["props_map"] = case_props_map(case)
        rec["coords_contract"] = coords_contract(case, problem)
    rec.update(extra or {})
    return rec


def coords_contract(case: CorpusCase, problem) -> dict:
    """Where COORDS come from and where the mechanics happen (disclosure)."""
    element = getattr(problem, "material_element", None)
    if element is None:
        return {"coords": "probe_reference_positions", "mechanics": "probe brick",
                "mixed_geometry": False,
                "statement": "COORDS are the probe integration points' reference positions"}
    return {"coords": "author_element_mapped", "mechanics": "probe brick",
            "mixed_geometry": True,
            "author_element_corners": np.asarray(element).tolist(),
            "provenance": case.extra.get("node_provenance", ""),
            "statement": "COORDS are the probe integration points mapped (trilinear) into one "
                         "element of the author's mesh; the mechanics -- geometry, strains, "
                         "residual -- stay on the probe brick. COORDS do not move with u."}


def unsupported_records(case: CorpusCase, failure: ProviderBuildFailed, problems, kin,
                        status="unsupported") -> List[dict]:
    out = []
    for problem in problems:
        for feature in ("residual_sens", "global_sens"):
            out.append({"schema": SCHEMA, "key": case.key, "source_id": case.source_id,
                        "family": case.family, "kinematics": case.kinematics,
                        "kinematics_contract": kin.name,
                        "problem": problem.name, "feature": feature,
                        "status": status, "failure_class": failure.failure_class,
                        "reason": failure.detail[-800:],
                        "material_path": "none: the live OTI provider could not be built; "
                                         "fixture playback (residual_core.materials."
                                         "fixture_playback) would supply PRIMAL stress/tangent "
                                         "at one frozen increment only and carries no "
                                         "parameter derivative, so neither feature can be "
                                         "assessed through it",
                        "reference": REFERENCE[feature], "what": WHAT[feature],
                        "held_fixed": HELD[feature], "max_abs": None, "max_rel": None,
                        "tolerance": None, "evidence": None})
    return out


def _history_calls(provider, problem, run, kinematics) -> List[dict]:
    """Keyword arguments of every ORIGINAL call of a converged history (probe input)."""
    asm = Assembly(problem, provider.finite, kinematics)
    drv = Driver(provider, asm)
    calls = []
    for n, inc in enumerate(run.schedule or problem.schedule()):
        U_prev = run.U[n - 1] if n else np.zeros(asm.ndof)
        inp = asm.inputs(U_prev, run.U[n], run.incoming[n])
        kw = drv._kw(inc, n + 1)
        calls.append(dict(props=np.asarray(run.props), stress=inp["stress_in"],
                          state=run.incoming[n].state, stran=inp["stran_in"],
                          dstran=inp["dstran"], F0=inp["F0"], F1=inp["F1"], drot=inp["drot"], **kw))
    return calls


def run_case(case: CorpusCase, out_root: Path, *, problems: Optional[Sequence[M.Problem]] = None,
             quick: bool = False, features=("residual_sens", "global_sens", "assembly"),
             kinematics=None) -> List[dict]:
    out_root = Path(out_root)
    kin = KINEMATICS[kinematics] if isinstance(kinematics, str) else \
        (kinematics or default_kinematics(case.finite))
    problems = list(problems or default_problems(case, quick=quick))
    provider_dir = out_root / "providers" / case.key
    try:
        record = build_provider_for(case, provider_dir, umat_repo=paths().umat_repo)
        provider = CorpusProvider(record, case, provider_dir / "lib")
    except ProviderBuildFailed as failure:
        return unsupported_records(case, failure, problems, kin)
    if provider.finite and kin.rotation == "abaqus" and provider.reads_drot:
        return unsupported_records(case, ProviderBuildFailed(
            "drot_read_without_derivative_seed",
            "the source reads DROT; UMAT_OTI_EVAL_TOTAL_F passes DROT without a derivative "
            "seed, so d DROT/d u and d DROT/d p would be missing from K and from every "
            "derivative"), problems, kin)
    evidence_dir = out_root / "evidence" / case.key / kin.name
    evidence_dir.mkdir(parents=True, exist_ok=True)
    records = []
    patch = None
    if "assembly" in features:
        try:
            patch = patch_check(provider, case, kin)
        except (NewtonFailed, MaterialCallError, np.linalg.LinAlgError) as error:
            patch = {"pass": False, "error": "%s: %s" % (type(error).__name__, error)}
    probe_result = None
    for problem in problems:
        started = time.time()
        evidence = evidence_dir / ("%s.json" % problem.name)
        detail = {"case": {k: (str(v) if isinstance(v, Path) else v)
                           for k, v in asdict(case).items()},
                  "kinematics": {"name": kin.name, "detail": kin.describe()},
                  "problem": {"name": problem.name, "description": problem.description,
                              "path": problem.path, "qois": problem.qois,
                              "ndof": problem.mesh.ndof, "free": int(problem.free.size),
                              "coords_from_author_element": problem.material_element},
                  "provider": {k: v for k, v in record.items() if k != "completed_contract"}}
        try:
            # the exact matrix on the planned increments first (B2); if that does
            # not converge, other steering (the converged state is the same
            # R = 0) and increment cut-backs as Abaqus would
            analytic = analytic_solve(provider, problem, case.props, kin, detail)
        except (NewtonFailed, MaterialCallError, np.linalg.LinAlgError) as error:
            reason = "%s: %s" % (type(error).__name__, error)
            analytic = None
            try:
                ref = reference_solve(provider, problem, case.props, kin, None, cutbacks=CUTBACKS)
                analytic = steered_analytic(provider, problem, case.props, kin, ref, detail)
                status, failure = "failed", "oti_primal_solve_failed_original_succeeds"
            except (NewtonFailed, MaterialCallError, np.linalg.LinAlgError) as again:
                status, failure = "unsupported", "generic_problem_inadmissible"
                reason += " | ORIGINAL on the same path: %s: %s" % (type(again).__name__, again)
        if analytic is None:
            for feature in ("residual_sens", "global_sens"):
                records.append(_record(case, problem, feature, record,
                                       {"status": status}, evidence, kin,
                                       {"failure_class": failure, "reason": reason}))
            detail["error"] = reason
            evidence.write_text(json.dumps(detail, indent=1, default=_json_default), encoding="utf-8")
            continue
        if analytic.cutbacks:
            analytic = common_schedule(provider, problem, case.props, kin, analytic, detail)
        # uninitialised-memory probe, once per source, on the first solved history
        # every problem is probed, on its own solved history (Vera, B8 review)
        probe_result = hidden_state.probe(record, case, provider_dir,
                                          _history_calls(provider, problem, analytic, kin),
                                          undefined_statev=case.extra.get("undefined_statev") or ())
        detail["hidden_state_probe"] = probe_result
        if probe_result["status"] == "trip":
            evidence.write_text(json.dumps(detail, indent=1, default=_json_default), encoding="utf-8")
            return _not_attempted(case, problems, record, kin, evidence,
                                  "uninitialised-memory probe: %s" % probe_result["statement"],
                                  probe_result)
        defined = [k for k in range(case.nstatv)
                   if k + 1 not in (case.extra.get("undefined_statev") or ())]
        moved = [float(np.max(np.abs(o.state[:, defined] - i.state[:, defined])))
                 if o.state.size and defined else 0.0
                 for i, o in zip(analytic.incoming, analytic.outgoing)]
        activity = {"increments_with_state_change": int(sum(m > 0.0 for m in moved)),
                    "max_state_change": max(moved) if moved else 0.0,
                    "max_abs_dstate_dp": max(float(np.max(np.abs(d))) if d.size else 0.0
                                             for d in analytic.dstate_dp)}
        parity = max(analytic.primal_parity)
        parity_note = {} if parity <= 1e-12 else {
            "failure_class_if_failed": "oti_primal_differs_from_original",
            "primal_parity_max": parity}
        solve = {"cutbacks": len(analytic.cutbacks),
                 "increments_planned": len(problem.schedule()),
                 "newton_matrix": sorted(set(analytic.newton_matrix)),
                 "converged_by_roundoff_floor": analytic.converged_by.count("roundoff_floor")}
        detail["analytic"] = {"solve": solve, "cutback_log": analytic.cutbacks,
                              "converged_by": analytic.converged_by,
                              "roundoff_floor": analytic.roundoff_floor,
                              "newton": analytic.newton, "newton_matrix": analytic.newton_matrix,
                              "backtracks": analytic.backtracks, "activity": activity,
                              "primal_parity_max": parity,
                              "sensitivity_equilibrium_max": max(analytic.sensitivity_equilibrium),
                              "seconds": analytic.seconds}
        common = {"increments": len(analytic.U), "activity": activity, "solve": solve,
                  "hidden_state_probe": (probe_result or {}).get("status"),
                  "primal_parity_max": parity}
        try:
            if "residual_sens" in features:
                comps, info = residual_sensitivity(provider, case, problem, analytic, kin)
                summary = summarise(comps)
                summary["tolerance"] = {"rtol": RTOL_LOCAL, "atol": "%g * max(sum_e|f_e|) / |p|"
                                        % ATOL_FACTOR, "entry": "atol + rtol|D_e| + 2 u_e",
                                        "plateau_steps": PLATEAU, "eps_eval": info["eps_eval"]}
                detail["residual_sens"] = {"summary": summary, "info": info, "comparisons": comps}
                extra = dict(common, h0_replays_bit_exact=info["h0_replays_bit_exact"])
                if summary["status"] == "failed" and parity_note:
                    extra.update(failure_class=parity_note["failure_class_if_failed"])
                records.append(_record(case, problem, "residual_sens", record, summary, evidence,
                                       kin, extra))
            if "global_sens" in features:
                try:
                    comps, info = global_sensitivity(provider, case, problem, analytic, kin)
                    summary = summarise(comps)
                except (NewtonFailed, MaterialCallError, np.linalg.LinAlgError) as error:
                    # the ORIGINAL could not be solved at the nominal PROPS on the
                    # analytic solve's increments: there is no reference, which
                    # says nothing about the derivative (not a failed comparison)
                    comps, info = [], {"error": "%s: %s" % (type(error).__name__, error)}
                    summary = {"status": "unsupported", "reason": info["error"],
                               "reason_class": "reference_solve_failed"}
                summary["tolerance"] = {"rtol": RTOL_GLOBAL, "atol": "%g * scale / |p|" % ATOL_FACTOR,
                                        "entry": "atol + rtol|D_e| + 2 u_e",
                                        "plateau_steps": PLATEAU, "eps_eval": info.get("eps_eval")}
                qoi = summarise([c for c in comps if c.get("quantity") != "u"])
                detail["global_sens"] = {"summary": summary, "qoi_summary": qoi, "info": info,
                                         "comparisons": comps}
                records.append(_record(case, problem, "global_sens", record, summary, evidence, kin, dict(
                    common, qoi=[q["name"] for q in problem.qois],
                    qoi_status=qoi["status"], qoi_max_rel=qoi.get("max_rel"),
                    reference_newton=info.get("reference_newton_matrices"),
                    reference_resolve_failures=len(info.get("reference_resolve_failures", [])),
                    reference_on_another_equilibrium=info.get("reference_on_another_equilibrium"),
                    nominal_rerun_bit_exact=info.get("nominal_rerun_bit_exact"),
                    sensitivity_equilibrium_max=max(analytic.sensitivity_equilibrium),
                    **({"failure_class": parity_note["failure_class_if_failed"]}
                       if summary["status"] == "failed" and parity_note else {}))))
        except HiddenStateTrip as trip:
            evidence.write_text(json.dumps(detail, indent=1, default=_json_default), encoding="utf-8")
            return _not_attempted(case, problems, record, kin, evidence, str(trip), probe_result)
        if "assembly" in features:
            records.append(_assembly_record(provider, case, problem, analytic, kin, record,
                                            evidence, detail, patch))
        steered = detail.get("steered_onto_reference")
        if steered:
            # the analytic solve only converged started on the original's solution:
            # say so on every record of this problem, and what it came to
            for rec in records:
                if rec.get("problem") == problem.name and rec["feature"] in (
                        "residual_sens", "global_sens"):
                    rec["steered_onto_reference"] = steered
                    if rec["status"] != "verified":
                        verified_n = (rec.get("counts") or {}).get("verified", 0)
                        rec["reason"] = ("%s; steered onto the original's equilibrium "
                                         "(steered_onto_reference); %s, %d verified comparisons"
                                         % (rec.get("reason") or "", rec["status"], verified_n))
        detail["seconds"] = time.time() - started
        evidence.write_text(json.dumps(detail, indent=1, default=_json_default), encoding="utf-8")
    provider.flush()
    return records


def _not_attempted(case, problems, record, kin, evidence, reason, probe_result):
    """A hidden-state trip: the whole source is not attempted (Vera B1 A)."""
    out = []
    for problem in problems:
        for feature in ("residual_sens", "global_sens", "assembly_consistency"):
            out.append(_record(case, problem, feature, record, {"status": "not_attempted"},
                               evidence, kin, {"failure_class": "hidden_state_trip",
                                               "reason": reason,
                                               "hidden_state_probe": (probe_result or {}).get("status")}))
    return out


def _assembly_record(provider, case, problem, analytic, kin, record, evidence, detail, patch):
    checks = {"equilibrium": equilibrium_check(analytic),
              "newton": newton_rates(analytic),
              "primal_parity": {"pass": max(analytic.primal_parity) <= parity_tolerance(case),
                                "max_rel": max(analytic.primal_parity),
                                "tolerance": parity_tolerance(case)},
              "patch_test": patch}
    try:
        checks["tangent"] = tangent_check(provider, case, problem, analytic, kin)
    except (NewtonFailed, MaterialCallError, np.linalg.LinAlgError) as error:
        checks["tangent"] = [{"error": str(error)}]
    if provider.finite:
        try:
            checks["objectivity"] = objectivity_check(provider, case, problem, kin)
        except (NewtonFailed, MaterialCallError, np.linalg.LinAlgError) as error:
            checks["objectivity"] = {"pass": False, "error": "%s: %s" % (type(error).__name__, error)}
        checks["kinematics_delta"] = kinematics_delta(provider, case, problem)
    tangents = checks["tangent"]
    exact = [t.get("K_exact_oti", {}).get("verdict") for t in tangents]
    exact_ok = bool(exact) and all(v == "verified" for v in exact)
    patch_ok = bool(patch and (patch.get("pass") or patch.get("applicable") is False))
    parts = {"equilibrium": checks["equilibrium"]["pass"], "newton": checks["newton"]["pass"],
             "primal_parity": checks["primal_parity"]["pass"], "patch": patch_ok,
             "K_exact": exact_ok}
    if provider.finite:
        parts["objectivity"] = bool(checks["objectivity"].get("pass"))
    status = "verified" if all(parts.values()) else "failed"
    detail["assembly"] = checks
    worst = lambda name, key: max((t.get(name, {}).get(key) or 0.0 for t in tangents), default=None)
    verdicts = lambda name: [t.get(name, {}).get("verdict") for t in tangents]
    return _record(case, problem, "assembly_consistency", record,
                   {"status": status, "max_abs": worst("K_exact_oti", "max_abs"),
                    "max_rel": worst("K_exact_oti", "max_rel"),
                    "tolerance": {"tangent_rtol": RTOL_TANGENT, "patch": PATCH_TOL,
                                  "equilibrium": 1e-12, "objectivity": OBJECTIVITY_TOL}},
                   evidence, kin, {
        "what": "assembled R: patch test, net force, K vs FD of R(u) at the end of every path "
                "segment, Newton order, OTI-vs-original primal parity, objectivity under a "
                "superposed rigid rotation (finite strain)",
        "wrt": "nodal displacements u (tangent check)",
        "held_fixed": "incoming state of each checked increment",
        "reference": "patch: exact linear field; K: centred FD of R assembled from the "
                     "ORIGINAL UMAT, relative steps %s of max|u|, same entrywise rule; "
                     "objectivity: the unrotated run" % (list(TANGENT_STEPS),),
        "parts": parts,
        "checks": {"equilibrium": checks["equilibrium"]["pass"],
                   "newton_order": checks["newton"]["median_order"],
                   "newton_max_iterations": checks["newton"]["max_iterations"],
                   "primal_parity": checks["primal_parity"]["max_rel"],
                   "patch": (patch and patch.get("pass")) if (patch or {}).get(
                       "applicable", True) else "not_applicable",
                   "K_checked_increments": [t.get("increment") for t in tangents],
                   "K_exact_oti": exact,
                   "K_ddsdde_author": verdicts("K_ddsdde_author"),
                   "K_ddsdde_author_max_rel": worst("K_ddsdde_author", "max_rel"),
                   "K_ddsdde_oti": verdicts("K_ddsdde_oti"),
                   "K_ddsdde_oti_max_rel": worst("K_ddsdde_oti", "max_rel"),
                   "objectivity": {k: v for k, v in (checks.get("objectivity") or {}).items()
                                   if k not in ("rotation",)},
                   "kinematics_delta": checks.get("kinematics_delta")}})


def _json_default(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, set):
        return sorted(value)
    return str(value)
