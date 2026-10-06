"""Why the corpus residual run could not solve most routine-verified sources,
and the fixes (B8). Each test fails without its fix.

1. Newton stalled at the residual's round-off floor (nearly incompressible
   neo-Hookean, bulk/shear ~ 5e3: relative free residual ~4e-12 against a
   1e-12 tolerance) -- ``engine.ROUNDOFF_FLOOR_FACTOR``.
2. The probe brick sat at [0,1]^3 while the routine caches its point's
   position (COORDS -> STATEV) and grows by it; the verification run had
   recorded one element of the author's mesh -- ``mesh.material_coordinates``.
3. The probe ran 3x past the analysis time the verification run's loading
   spans (a clock normalised by TOTALT) -- ``loading_period_total``.
4. The hidden-state probe tripped on a STATEV entry the verification run had
   already found undefined in the original; it is now left out only if the
   routine provably never reads it -- ``hidden_state.probe(undefined_statev=)``.
5. The exact iteration matrix alone diverged between clamps where the DDSDDE
   matrix converges -- ``newton_matrix='ra_then_exact'``.
6. The objectivity check divided round-off by round-off when dV/dp is zero in
   exact arithmetic -- yardstick max(|dV/dp|, |V|/|p|).
7. No increment control: an increment Abaqus would cut back failed the whole
   path -- ``run_history(cutbacks=)``; the re-solves at p +/- h replay the
   increments the analytic solve took (``schedule=``).
8. A reference solve of the ORIGINAL failed only because its own DDSDDE
   iteration matrix stalled -- ``runner.reference_solve`` tries the exact and
   the mixed matrices on the same increments.
9. The FD ladder ended at h/|p| = 1e-7, before a parameter near a singular
   limit (nu = 0.4995) shows a plateau -- three smaller steps.
10. A reference that cannot be solved was booked as a FAILED derivative; it
   is ``unsupported`` (reason_class reference_solve_failed).
11. The analytic solve had one strategy; ``runner.analytic_solve`` tries the
   exact matrix first (B2), then other steering and cut-backs.
12. Analytic and original builds agree to round-off, yet at a hard increment
   Newton may converge for one and not the other: ``runner.common_schedule``
   gives both one set of increments.

B9 (Vera's B8 review):
13. A ladder step whose p +/- h leaves the slot's declared domain (Poisson's
   ratio in (-1, 1/2), moduli > 0, read from the source's NAME = PROPS(i)) is
   dropped and tagged ``outside_parameter_domain``; the zero-within-resolution
   bound is taken over the kept steps only (f7be16bc, nu = 0.499).
14. A floor acceptance needs one confirming Newton step; a step that halves
   the residual and stays under the floor is accepted, as is one that does
   not halve it; the better state is kept.
15. The unread check writes 0, -1.2345e30, +1.2345e30 and NaN; every problem
   is probed; records name COORDS as held fixed and disclose the mixed geometry.
16. Past a bifurcation (K_ff indefinite) the ORIGINAL's nominal solve can land
   on another equilibrium than the analytic solve; a difference around it is
   no reference for the analytic derivative. Such increments are unresolved
   (reference_on_another_equilibrium), never judged (B9: SeaShell, |dV| 0.66),
   but only while primal parity holds, and from the first departure on unless
   solution and incoming state both match again (Vera R-1).
17. Where only the ORIGINAL converges (0e56ba3c: K_ff indefinite at increment
   9), the analytic solve is started on the original's solution and, if it
   converges there with its own residual, judged on that path.
18. Labels only (Vera B9): a record whose analytic solve was steered says so and
   what it came to; problems excused by a branch departure have their own
   reason class, distinct from fd_reference_unresolved; the per-source tally
   quotes "verified on k of n problems" and flags a cell that rests on 1 of 3.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from residual_core.corpus import hidden_state  # noqa: E402
from residual_core.corpus import mesh as M  # noqa: E402

#: identity growth (G = I): nearly incompressible neo-Hookean, nu = 0.4999,
#: STATEV(7) undefined in the original (an unassigned local is stored there)
FLAT_SOURCE = ("Jeff97__Programming-Plane-Strain-Plates-through-Growth-Under-Body-Forces/"
               "Examples-In-Section-3/Flat/Th01/PureGrowth.for",
               "4cea0702a6151771f58c45cee174c7ac4ae3192884cc505da240f323a77c982f")
#: growth G11 = 1 + (pi/2 - pi Y - 1) t / TOTALT with Y = COORDS(2) cached in STATEV(2)
ARC_DOWN_SOURCE = ("Jeff97__Programming-Plane-Strain-Plates-through-Growth-Under-Body-Forces/"
                   "Examples-In-Section-3/ArcDown/Th01/PureGrowth.for",
                   "c5b21e39d4b10ebac1d3ea8ae592ff3df0dd4fc556f61b349ffcf43a97890cb8")
Z2_SOURCE = ("Jeff97__Realization-of-planar-and-surface-conformal-mappings/"
             "Mesh_Convergence_test/2D/Z2/10/Growth-Z2.for",
             "b41f018c4396eab01dd9306aff57bd3978f79df9a3929aa5028319688af26cb1")
ARC_SOURCE = ("Jeff97__Programming-Plane-Strain-Plates-through-Growth-Under-Body-Forces/"
              "Examples-In-Section-3/ArcUp/Th002/PureGrowth.for",
              "f7a30a6be5719b7ca68c4e6aa1197e6a7536c6fa49b64b125c3622c03a88b37a")


# --------------------------------------------------------------------------- #
# pure python
# --------------------------------------------------------------------------- #
def test_probe_points_map_into_the_authors_element_corner_to_corner():
    element = np.array([[-1.03527617, -3.86370325, 1.439394], [-1.06115806, -3.96029592, 1.439394],
                        [-0.91233581, -3.99720454, 1.439394], [-0.89008373, -3.89971161, 1.439394],
                        [-1.03527617, -3.86370325, 1.5151515], [-1.06115806, -3.96029592, 1.5151515],
                        [-0.91233581, -3.99720454, 1.5151515], [-0.89008373, -3.89971161, 1.5151515]])
    for divisions in ((1, 1, 1), (2, 2, 2)):
        mesh = M.brick(divisions, size=(2.0, 1.0, 0.5))
        corners = [mesh.node_at(p) for p in [(0, 0, 0), (2, 0, 0), (2, 1, 0), (0, 1, 0),
                                              (0, 0, .5), (2, 0, .5), (2, 1, .5), (0, 1, .5)]]
        mapped = M.material_coordinates(mesh.coords, mesh, element)
        assert np.allclose(mapped[corners], element, atol=1e-14)
        centre = M.material_coordinates(np.array([[1.0, 0.5, 0.25]]), mesh, element)
        assert np.allclose(centre[0], element.mean(axis=0), atol=1e-14)


def _case(**extra):
    from residual_core.corpus.sources import CorpusCase
    return CorpusCase(key="f" * 24, source_id="x", source_path=Path("x.f"), source_form="fixed",
                      family="growth", terminal_state="fully_verified",
                      verification_fingerprint="", kinematics="finite", ntens=6, nstatv=9,
                      props=[1.0e6], material_provenance="", element_type="C3D8H",
                      extra=dict({"initial_statev": [], "undefined_statev": [],
                                  "node_coordinates": []}, **extra))


def test_the_probe_spans_the_analysis_time_of_the_verification_runs_loading():
    from residual_core.corpus.runner import default_problems
    for total in (1.0, 10.0, 2.2):
        for problem in default_problems(_case(loading_period_total=total)):
            last = problem.schedule()[-1]
            assert last["time_start"] + last["dtime"] == pytest.approx(total, rel=1e-12)
            assert "analysis time" in problem.description
    # no stated loading: one unit of time per path segment, as before
    problem = default_problems(_case(loading_period_total=None))[0]
    last = problem.schedule()[-1]
    assert last["time_start"] + last["dtime"] == pytest.approx(len(problem.path))


def test_an_authors_element_sets_coords_on_every_problem():
    from residual_core.corpus.runner import default_problems
    rows = [[k + 1] + list(map(float, c)) for k, c in enumerate(
        [(0, .01, 0), (.0026, .01, 0), (.0026, .0125, 0), (0, .0125, 0),
         (0, .01, .001), (.0026, .01, .001), (.0026, .0125, .001), (0, .0125, .001)])]
    for problem in default_problems(_case(node_coordinates=rows)):
        assert problem.material_element is not None
        from residual_core.corpus.engine import Assembly
        coords = Assembly(problem, True).coords
        assert coords[:, 1].min() > 0.01 and coords[:, 1].max() < 0.0125
    assert all(p.material_element is None for p in default_problems(_case()))


def test_a_cut_back_halves_load_factor_and_time_and_keeps_the_segment_end_last():
    from residual_core.corpus.engine import _halves
    inc = {"lambda_start": 0.4, "lambda": 1.0, "dtime": 0.2, "time_start": 1.0,
           "segment": 1, "segment_end": True}
    a, b = _halves(inc)
    assert (a["lambda_start"], a["lambda"], b["lambda_start"], b["lambda"]) == (0.4, 0.7, 0.7, 1.0)
    assert (a["time_start"], a["dtime"], b["time_start"], b["dtime"]) == (1.0, 0.1, 1.1, 0.1)
    assert (a["segment_end"], b["segment_end"]) == (False, True)
    assert a["cutback_depth"] == b["cutback_depth"] == 1
    assert [x["cutback_depth"] for x in _halves(a)] == [2, 2]


def test_a_poisson_ratio_near_one_half_is_resolved_by_the_ladder():
    """R ~ K(nu) with K = 1/(1 - 2 nu): the derivative a central difference must
    resolve at nu = 0.4995. On the old ladder (down to 1e-7) no 3-step plateau."""
    from residual_core.corpus.verify import STEPS, adjudicate
    nu = 0.4995

    def R(v):
        return np.array([1.0 / (1.0 - 2.0 * v), 0.3 / (1.0 - 2.0 * v) + v])

    exact = np.array([2.0 / (1.0 - 2.0 * nu) ** 2, 0.6 / (1.0 - 2.0 * nu) ** 2 + 1.0])

    def judge(steps):
        hs = {s: s * nu for s in steps}
        return adjudicate(exact, {s: R(nu + h) for s, h in hs.items()},
                          {s: R(nu - h) for s, h in hs.items()}, R(nu), hs,
                          rtol=1e-6, atol=0.0)["verdict"]
    old = [s for s in STEPS if s >= 1e-7]
    assert judge(old) == "unresolved"
    assert judge(STEPS) == "verified"
    # a wrong derivative is not rescued by the extra steps
    hs = {s: s * nu for s in STEPS}
    assert adjudicate(exact * (1 + 1e-4), {s: R(nu + h) for s, h in hs.items()},
                      {s: R(nu - h) for s, h in hs.items()}, R(nu), hs,
                      rtol=1e-6, atol=0.0)["verdict"] == "failed"


class _Run:
    def __init__(self, schedule, matrix="exact", cutbacks=()):
        self.schedule, self.newton_matrix_arg, self.cutbacks = schedule, matrix, list(cutbacks)


def test_the_analytic_solve_falls_back_in_order_and_says_so(monkeypatch):
    from residual_core.corpus import runner
    from residual_core.corpus.engine import NewtonFailed
    tried = []

    def fake(provider, problem, props, **kw):
        tried.append((kw["newton_matrix"], kw["cutbacks"]))
        if len(tried) < 3:
            raise NewtonFailed("increment 1: stand-in")
        return _Run([1, 2], kw["newton_matrix"])
    monkeypatch.setattr(runner, "run_history", fake)
    detail = {}
    run = runner.analytic_solve(None, None, [1.0], None, detail)
    assert tried == [("exact", 0), ("ra_then_exact", 0), ("ra_stalled_exact", runner.CUTBACKS)]
    assert run.newton_matrix_arg == "ra_stalled_exact" and len(detail["solve_attempts"]) == 2


def test_analytic_and_reference_solves_share_the_finer_increments(monkeypatch):
    from residual_core.corpus import runner
    coarse, fine = [{"i": 1}, {"i": 2}], [{"i": 1}, {"i": 1.5}, {"i": 2}]
    analytic = _Run(coarse, "ra_stalled_exact", [{"increment": 1}])
    calls = []

    def fake_reference(provider, problem, props, kin, schedule, first=None, **kw):
        calls.append(("reference", schedule, first, kw.get("cutbacks")))
        return _Run(fine, cutbacks=[{"increment": 2}])

    def fake_history(provider, problem, props, **kw):
        calls.append(("analytic", kw["schedule"], kw["newton_matrix"], kw.get("cutbacks", 0)))
        return _Run(kw["schedule"], kw["newton_matrix"])
    monkeypatch.setattr(runner, "reference_solve", fake_reference)
    monkeypatch.setattr(runner, "run_history", fake_history)
    detail = {}
    run = runner.common_schedule(None, None, [1.0], None, analytic, detail)
    assert calls[0] == ("reference", coarse, "auto", runner.CUTBACKS)
    assert calls[1] == ("analytic", fine, "ra_stalled_exact", 0)
    assert run.schedule == fine and len(run.cutbacks) == 2
    assert "re-solved" in detail["common_schedule"]
    # the reference needed no cut-back of its own: the analytic run is kept
    calls.clear()
    monkeypatch.setattr(runner, "reference_solve", lambda *a, **k: _Run(coarse))
    assert runner.common_schedule(None, None, [1.0], None, analytic) is analytic


def test_the_routine_verified_sources_are_read_from_the_manifest(tmp_path, monkeypatch):
    import json
    from residual_core.corpus.sources import paths, routine_verified_keys
    rows = []
    for key, eligible, ddsdde in (("a" * 24, "verified", "verified"), ("b" * 24, "verified", "not_attempted"),
                                  ("c" * 24, "blocked", "blocked"), ("d" * 24, "verified", "verified")):
        rows.append({"registry": {"key": key, "terminal_state": "tangent_not_verified"},
                     "pipeline": {"eligible": {"status": eligible}},
                     "features": {"ddsdde": {"status": ddsdde}}})
    target = tmp_path / "final-umat" / "paper_results" / "corpus" / "manifest"
    target.mkdir(parents=True)
    (target / "corpus_manifest.json").write_text(json.dumps({"rows": rows}))
    monkeypatch.setenv("CORPUS_WORKSPACE", str(tmp_path))
    monkeypatch.delenv("UMAT_OTI_REPO", raising=False)
    assert routine_verified_keys(paths()) == ["a" * 24, "d" * 24]


POISSON_SOURCE = """      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
C     PROPS(2) - NU
        EMOD=PROPS(1)
        ENU=PROPS(2)
        ETA = PROPS(3) ! a viscosity, no declared domain
      END
"""


HELPER_FIRST = """      SUBROUTINE HELPER(PROPS, N)
      DIMENSION PROPS(N)
      XNU = PROPS(1)
      END SUBROUTINE HELPER
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
        EMOD=PROPS(1)
      ENDIF
      END
      REAL FUNCTION G(PROPS)
      ANU = PROPS(2)
      END
"""


def test_the_props_map_reads_only_the_entry_routine():
    """Vera hardening: a helper's own PROPS dummy is another array."""
    from residual_core.corpus import parameters as P
    mapping = P.props_map(HELPER_FIRST, "fixed")
    assert [(m["index"], m["name"]) for m in mapping] == [(1, "EMOD")]
    assert P.domain_of(mapping, 1)[:2] == (0.0, float("inf"))


def test_the_props_map_reads_declared_domains_from_the_source():
    from residual_core.corpus import parameters as P
    mapping = P.props_map(POISSON_SOURCE, "fixed")
    assert [(m["index"], m["name"]) for m in mapping] == [(1, "EMOD"), (2, "ENU"), (3, "ETA")]
    assert P.domain_of(mapping, 2)[:2] == (-1.0, 0.5)
    assert P.domain_of(mapping, 1)[:2] == (0.0, float("inf"))
    assert P.domain_of(mapping, 3) is None


def test_a_step_across_nu_one_half_is_dropped_and_the_zero_bound_uses_the_kept_steps(tmp_path):
    """f7be16bc: nu = 0.499, the steps 0.1 ... 3e-3 cross nu = 1/2 where the
    routine evaluates another material. Here q = c nu below 1/2 and its mirror
    above: every difference sits inside its round-off bound, and with the full
    ladder the zero bound comes from a crossing step and FAILS a right A = c."""
    from residual_core.corpus import runner
    from residual_core.corpus.verify import adjudicate
    source = tmp_path / "u.for"
    source.write_text(POISSON_SOURCE)
    case = _case()
    case.source_path, case.props = source, [906512.0, 0.499, 1.0]
    nu, c, big = 0.499, 1000.0, 1e16

    def q(v):
        return np.array([c * v if v < 0.5 else c * (1.0 - v)])
    full = runner._steps(nu)
    kept, note = runner._ladder(case, 2, nu)
    assert set(full) - set(kept) == {0.1, 0.03, 0.01, 0.003}
    assert {d["reason"] for d in note["dropped_steps"]} == {"outside_parameter_domain"}
    assert note["parameter_domain"]["open_interval"] == [-1.0, 0.5]

    def judge(hs):
        return adjudicate(np.array([c]), {s: q(nu + h) for s, h in hs.items()},
                          {s: q(nu - h) for s, h in hs.items()}, q(nu), hs, rtol=1e-6,
                          atol=0.0, value_scale=big)
    assert judge(full)["verdict"] == "failed"
    kept_verdict = judge(kept)
    assert kept_verdict["verdict"] != "failed"
    assert kept_verdict["counts"]["zero_within_resolution"] == 1
    # E has a domain (> 0) that no step of its ladder leaves: nothing dropped
    assert runner._ladder(case, 1, 906512.0)[1]["dropped_steps"] == []
    # no declared domain: the whole ladder
    assert runner._ladder(case, 3, 1.0)[0] == runner._steps(1.0)


def test_the_held_fixed_texts_name_coords():
    from residual_core.corpus import runner
    for feature in ("residual_sens", "global_sens"):
        assert "COORDS" in runner.HELD[feature] and "coords_contract" in runner.HELD[feature]


def test_branch_excused_problems_have_their_own_reason_class():
    from residual_core.corpus.verify import summarise
    excused = {"verdict": "unresolved", "reason": "reference_on_another_equilibrium: ..."}
    fd = {"verdict": "unresolved"}
    ok = {"verdict": "verified"}
    assert summarise([ok, excused, excused])["reason_class"] == "reference_on_another_equilibrium"
    assert summarise([ok, excused, fd])["reason_class"] == "fd_reference_unresolved"
    assert summarise([ok, fd])["reason_class"] == "fd_reference_unresolved"
    assert summarise([ok, ok])["status"] == "verified"


def test_the_per_source_tally_says_k_of_n_and_flags_one_of_three():
    import importlib.util
    spec = importlib.util.spec_from_file_location("run_corpus_residual",
                                                  ROOT / "tools" / "run_corpus_residual.py")
    tool = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(tool)
    def rec(problem, status, feature="residual_sens"):
        return {"key": "k" * 24, "source_id": "src", "feature": feature, "problem": problem,
                "status": status}
    three = [rec(p, s) for p, s in (("a", "verified"), ("b", "unsupported"), ("c", "unsupported"))]
    cell = tool.fold_cell(three)
    assert (cell["cell"], cell["k"], cell["n"], cell["run"]) == ("verified", 1, 1, 3)
    assert cell["rests_on_1_of_3"] is True
    assert "verified on 1 of 1 evaluated problems (of 3 run)" in cell["text"]
    assert "rests on 1 of 3" in cell["text"]
    two = tool.fold_cell([rec("a", "verified"), rec("b", "verified"), rec("c", "unsupported")])
    assert (two["k"], two["n"], two["rests_on_1_of_3"]) == (2, 2, False)
    assert tool.fold_cell([rec("a", "verified"), rec("b", "failed"), rec("c", "verified")])["cell"] == "failed"


class _Branch:
    """A reference run of a linear model V_n = p a_n (+ offset_n); its incoming
    state at increment n is V_{n-1} (+ state_offsets)."""

    def __init__(self, p, offsets, state_offsets=(0.0, 0.0)):
        from types import SimpleNamespace
        a = np.array([[1.0, 2.0], [3.0, 4.0]])
        self.V = [p * a[n] + offsets[n] for n in range(2)]
        self.U = [v.copy() for v in self.V]
        self.incoming = [SimpleNamespace(stress=np.array([1.0 + state_offsets[n]]),
                                         state=np.zeros(0), stran=np.zeros(0)) for n in range(2)]
        self.force_scale, self.newton = [1.0, 1.0], [[1e-14], [1e-14]]
        self.newton_matrix, self.newton_matrix_arg = ["ra", "ra"], "auto"
        self.outputs = [{"stress": v.copy()} for v in self.V]

    def qoi(self, problem):
        return np.zeros((2, 0))


def _global(monkeypatch, tmp_path, analytic_offsets, state_offsets=(0.0, 0.0), parity=0.0,
            du_scale=1.0):
    from types import SimpleNamespace
    from residual_core.corpus import runner
    source = tmp_path / "u.for"
    source.write_text("      EMOD=PROPS(1)\n")
    case = _case()
    case.source_path, case.props = source, [2.0]
    monkeypatch.setattr(runner, "reference_solve",
                        lambda prov, prob, props, kin, sched, first=None, **kw:
                        _Branch(props[0], [0.0, 0.0]))
    monkeypatch.setattr(runner, "run_history", lambda *a, **k: _Branch(2.0, [0.0, 0.0]))
    mine = _Branch(2.0, analytic_offsets, state_offsets)
    analytic = SimpleNamespace(schedule=[{}, {}], newton_matrix_arg="exact", V=mine.V,
                               incoming=mine.incoming, primal_parity=[parity, parity],
                               du_dp=[du_scale * np.array([[1.0], [2.0]]),
                                      du_scale * np.array([[3.0], [4.0]])],
                               qoi_dp=lambda problem: np.zeros((2, 0, 1)))
    comps, info = runner.global_sensitivity(SimpleNamespace(slots=[1]), case,
                                            SimpleNamespace(qois=[]), analytic, None)
    return {c["increment"]: c for c in comps}, info


def test_a_reference_on_another_equilibrium_is_not_judged(monkeypatch, tmp_path):
    by_inc, info = _global(monkeypatch, tmp_path, [0.0, 0.5])      # inc 2: another branch
    assert by_inc[1]["verdict"] == "verified"
    assert by_inc[2]["verdict"] == "unresolved"
    assert by_inc[2]["reason"].startswith("reference_on_another_equilibrium")
    assert info["reference_on_another_equilibrium"] == [2]
    # without the offset both increments are judged (and verified)
    by_inc, info = _global(monkeypatch, tmp_path, [0.0, 0.0])
    assert {c["verdict"] for c in by_inc.values()} == {"verified"}
    assert not info["reference_on_another_equilibrium"]


def test_the_branch_gate_never_excuses_a_primal_defect(monkeypatch, tmp_path):
    """R-1: with primal parity broken, a different solution may be the defect
    itself; the comparisons stay judged, and a wrong derivative stays FAILED."""
    by_inc, info = _global(monkeypatch, tmp_path, [0.0, 0.5], parity=1e-6, du_scale=1.1)
    assert not info["reference_on_another_equilibrium"]
    assert by_inc[2]["verdict"] == "failed" and by_inc[1]["verdict"] == "failed"


def test_after_a_departure_only_a_full_rejoin_is_judged(monkeypatch, tmp_path):
    """R-1: from the first departure on, increments are excused unless both the
    solution and the incoming state match again."""
    by_inc, info = _global(monkeypatch, tmp_path, [0.5, 0.0], state_offsets=(0.0, 0.3))
    assert info["reference_on_another_equilibrium"] == [1, 2]     # inc 2: history differs
    by_inc, info = _global(monkeypatch, tmp_path, [0.5, 0.0])
    assert info["reference_on_another_equilibrium"] == [1]        # inc 2 fully rejoined
    assert by_inc[2]["verdict"] == "verified"


class _FakeLib:
    """``regular`` of a one-point routine: STRESS from F, STATEV(7) garbage
    (differs per build); ``reads`` makes STRESS depend on incoming STATEV(7)."""

    def __init__(self, garbage, reads):
        self.garbage, self.reads = garbage, reads

    def regular(self, **kw):
        state = np.array(kw["state"], dtype=float, copy=True)
        stress = np.array([[1.0, 2.0, 3.0, 0.0, 0.0, 0.0]])
        if self.reads:
            stress = stress + 1e-3 * np.nan_to_num(state[:, 6:7], nan=7.0)
        state[:, 6] = self.garbage
        return {"stress": stress, "state": state, "ddsdde": np.eye(6)[None], "pnewdt": np.ones(1)}


@pytest.mark.parametrize("reads", [False, True])
def test_a_declared_undefined_statev_is_excused_only_if_never_read(monkeypatch, reads):
    libs = {"snan": _FakeLib(np.nan, reads), "zero": _FakeLib(0.0, reads)}
    monkeypatch.setattr(hidden_state, "build_probe",
                        lambda record, case, out, which: (libs[which], []))
    calls = [{"state": np.zeros((1, 9))}, {"state": np.full((1, 9), 0.5)}]
    case = _case()
    # without the declaration: the garbage entry trips the probe (as in B2)
    assert hidden_state.probe({}, case, Path("."), calls)["status"] == "trip"
    result = hidden_state.probe({}, case, Path("."), calls, undefined_statev=[7])
    if reads:
        assert result["status"] == "trip" and result["declared_undefined_read"]
        assert "READ" in result["statement"]
    else:
        assert result["status"] == "clean", result
        assert result["calls_differing_only_in_declared_undefined"] == len(calls)
        assert "never read" in result["statement"]
        assert "nan" in result["statement"] and "1.2345e+30" in result["statement"]


# --------------------------------------------------------------------------- #
# live corpus UMATs
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="module")
def built(tmp_path_factory):
    from test_corpus_residual import _live
    _live()
    from residual_core.corpus.provider import CorpusProvider, build_provider_for
    from residual_core.corpus.sources import key_for_source, load_case, paths
    root = tmp_path_factory.mktemp("solver_robustness")
    out = {}
    for name, source in (("flat", FLAT_SOURCE), ("arc", ARC_SOURCE),
                         ("arc_down", ARC_DOWN_SOURCE), ("z2", Z2_SOURCE)):
        case = load_case(key_for_source(*source))
        record = build_provider_for(case, root / case.key, umat_repo=paths().umat_repo)
        out[name] = (case, CorpusProvider(record, case, root / case.key / "lib"), record, root)
    return out


@pytest.mark.integration
@pytest.mark.fortran
@pytest.mark.slow
def test_newton_converges_at_the_roundoff_floor_and_stalls_without_it(built, monkeypatch):
    from residual_core.corpus import engine
    from residual_core.corpus.runner import default_problems
    case, provider, _, _ = built["flat"]
    problem = default_problems(case, quick=True)[0]
    run = engine.run_history(provider, problem, case.props, material="regular")
    assert "roundoff_floor" in run.converged_by
    for how, h, floor, scale in zip(run.converged_by, run.newton, run.roundoff_floor,
                                    run.force_scale):
        if how == "roundoff_floor":       # accepted after one confirming Newton step:
            # under the floor again (halved or not) or not halved
            first = next(i for i, r in enumerate(h) if r * scale <= floor * 1.001)
            assert len(h) >= first + 2, "the confirming step is in the history"
            assert h[first + 1] >= 0.5 * h[first] or h[first + 1] * scale <= floor * 1.001
            assert h[-1] == min(h[first], h[first + 1])
    worst = max(h[-1] for h in run.newton)
    assert 1e-12 < worst < 1e-9          # above the tolerance, at the floor
    monkeypatch.setattr(engine, "ROUNDOFF_FLOOR_FACTOR", 0.0)
    with pytest.raises(engine.NewtonFailed):
        engine.run_history(provider, problem, case.props, material="regular")


@pytest.mark.integration
@pytest.mark.fortran
@pytest.mark.slow
def test_a_position_caching_growth_law_solves_only_inside_the_authors_element(built):
    from residual_core.corpus.engine import NewtonFailed, run_history
    from residual_core.corpus.provider import MaterialCallError
    from residual_core.corpus.runner import default_problems
    case, provider, _, _ = built["arc"]
    assert M.author_element(case) is not None
    problem = default_problems(case, quick=True)[0]
    run = run_history(provider, problem, case.props, material="regular")
    grown = np.array([o.state[:, 7] for o in run.outgoing])        # STATEV(8) = G11
    assert grown.max() > 1.0 + 1e-3, "the growth law must act"
    problem.material_element = None      # COORDS on the probe brick, Y in [0, 1]
    with pytest.raises((NewtonFailed, MaterialCallError)):
        run_history(provider, problem, case.props, material="regular")


@pytest.mark.integration
@pytest.mark.fortran
@pytest.mark.slow
def test_the_declared_undefined_statev_of_a_corpus_routine_is_unread(built):
    from residual_core.corpus.engine import run_history
    from residual_core.corpus.runner import _history_calls, default_problems
    case, provider, record, root = built["flat"]
    assert case.extra["undefined_statev"] == [7]
    problem = default_problems(case, quick=True)[0]
    run = run_history(provider, problem, case.props, material="regular")
    calls = _history_calls(provider, problem, run, None)
    assert hidden_state.probe(record, case, root / "p0", calls)["status"] == "trip"
    result = hidden_state.probe(record, case, root / "p1", calls, undefined_statev=[7])
    assert result["status"] == "clean", result


@pytest.mark.integration
@pytest.mark.fortran
@pytest.mark.slow
def test_the_ddsdde_matrix_steers_where_the_exact_one_alone_diverges(built):
    from residual_core.corpus.engine import NewtonFailed, run_history
    from residual_core.corpus.runner import default_problems
    case, provider, _, _ = built["flat"]
    problem = [p for p in default_problems(case) if p.name.startswith("clamped_tension")][0]
    with pytest.raises(NewtonFailed):
        run_history(provider, problem, case.props, material="oti", newton_matrix="exact")
    run = run_history(provider, problem, case.props, material="oti", newton_matrix="ra_then_exact",
                      sensitivities=True)
    assert len(run.U) == len(problem.schedule())
    assert max(run.sensitivity_equilibrium) < 1e-8


@pytest.mark.integration
@pytest.mark.fortran
@pytest.mark.slow
def test_objectivity_of_a_derivative_that_is_zero_in_exact_arithmetic(built):
    """dV/dC0 = 0 under displacement control when C0 scales the whole response:
    the check must not divide round-off by round-off (B8 measured 5.0)."""
    from residual_core.corpus.runner import default_problems, objectivity_check
    case, provider, _, _ = built["flat"]
    problem = default_problems(case, quick=True)[0]
    result = objectivity_check(provider, case, problem, "nlgeom")
    assert result["pass"], result


@pytest.mark.integration
@pytest.mark.fortran
@pytest.mark.slow
def test_an_increment_cut_back_solves_the_path_and_the_re_solves_replay_it(built):
    from residual_core.corpus.engine import NewtonFailed, run_history
    from residual_core.corpus.runner import check_increments, default_problems
    case, provider, _, _ = built["arc"]
    problem = [p for p in default_problems(case) if p.name.startswith("clamped_tension")][0]
    planned = len(problem.schedule())
    with pytest.raises(NewtonFailed):
        run_history(provider, problem, case.props, material="regular")
    run = run_history(provider, problem, case.props, material="regular", cutbacks=6)
    assert run.cutbacks and len(run.schedule) > planned
    assert len(check_increments(problem, run.schedule)) == len(problem.path)
    assert run.schedule[-1]["lambda"] == problem.path[-1][0]
    # a perturbed re-solve on exactly those increments, no cut-back of its own
    props = np.array(case.props) * (1.0 + 1e-4)
    again = run_history(provider, problem, props, material="regular", schedule=run.schedule)
    assert [i["lambda"] for i in again.schedule] == [i["lambda"] for i in run.schedule]


@pytest.mark.integration
@pytest.mark.fortran
@pytest.mark.slow
def test_a_reference_solve_is_not_lost_to_the_originals_own_iteration_matrix(built):
    from residual_core.corpus.engine import NewtonFailed, run_history
    from residual_core.corpus.runner import default_problems, reference_solve
    case, provider, _, _ = built["arc_down"]
    problem = [p for p in default_problems(case) if p.name.startswith("clamped_shear")][0]
    with pytest.raises(NewtonFailed):
        run_history(provider, problem, case.props, material="regular")
    run = reference_solve(provider, problem, case.props, None, None)
    assert run.newton_matrix_arg != "auto" and len(run.U) == len(problem.schedule())
    again = run_history(provider, problem, case.props, material="regular",
                        newton_matrix=run.newton_matrix_arg, schedule=run.schedule)
    assert all(a.tobytes() == b.tobytes() for a, b in zip(run.U, again.U))


@pytest.mark.integration
@pytest.mark.fortran
@pytest.mark.slow
def test_a_reference_that_cannot_be_solved_is_not_a_failed_derivative(built, monkeypatch):
    from residual_core.corpus import runner
    from residual_core.corpus.engine import NewtonFailed
    case, _, _, root = built["flat"]

    def no_reference(*args, **kw):
        raise NewtonFailed("increment 1: stand-in for a reference that does not converge")
    monkeypatch.setattr(runner, "global_sensitivity", no_reference)
    probed = []
    real_probe = runner.hidden_state.probe

    def spy(record, case, out, calls, **kw):
        probed.append(len(calls))
        return real_probe(record, case, out, calls, **kw)
    monkeypatch.setattr(runner.hidden_state, "probe", spy)
    problems = runner.default_problems(case, quick=True)
    second = M.clamped_shear(M.brick((1, 1, 1)), 0.1, path=((1.0, 2),))
    second.material_element = problems[0].material_element
    records = runner.run_case(case, root / "run", problems=problems + [second],
                              features=("global_sens",))
    assert len(probed) == 2, "every problem is probed"
    records = [r for r in records if r["feature"] == "global_sens"]
    assert len(records) == 2
    for record in records:
        assert record["status"] == "unsupported"
        assert record["reason_class"] == "reference_solve_failed"
        assert "COORDS" in record["held_fixed"]
        assert record["coords_contract"]["mixed_geometry"] is True
        assert record["coords_contract"]["mechanics"] == "probe brick"
        assert {"index": 1, "name": "C0"}.items() <= record["props_map"][0].items()


@pytest.mark.integration
@pytest.mark.fortran
@pytest.mark.slow
def test_the_analytic_solve_is_steered_onto_the_originals_equilibrium(built):
    from residual_core.corpus.engine import NewtonFailed, run_history
    from residual_core.corpus.runner import default_problems, reference_solve, steered_analytic
    case, provider, _, _ = built["z2"]
    problem = [p for p in default_problems(case) if p.name.startswith("clamped_shear")][0]
    ref = reference_solve(provider, problem, case.props, None, None, cutbacks=6)
    with pytest.raises(NewtonFailed):        # unsteered, on the very same increments
        run_history(provider, problem, case.props, material="oti", newton_matrix="exact",
                    schedule=ref.schedule)
    detail = {}
    run = steered_analytic(provider, problem, case.props, None, ref, detail)
    assert run is not None and len(run.U) == len(ref.U)
    assert detail["steered_onto_reference"]["max_V_difference"] <= 1e-6 * float(
        np.max(np.abs(np.array(ref.V))))
    assert run.reference_steering == ref.newton_matrix_arg
    assert len(run.du_dp) == len(ref.U)


@pytest.mark.integration
@pytest.mark.fortran
@pytest.mark.slow
def test_a_steered_record_says_so_and_what_it_came_to(built):
    from residual_core.corpus import runner
    case, _, _, root = built["z2"]
    problem = [p for p in runner.default_problems(case) if p.name.startswith("clamped_shear")][0]
    records = runner.run_case(case, root / "run_steered", problems=[problem],
                              features=("residual_sens", "global_sens"))
    for record in [r for r in records if r["feature"] in ("residual_sens", "global_sens")]:
        assert record["steered_onto_reference"]["max_V_difference"] < 1e-6
        if record["status"] != "verified":
            assert ("steered onto the original's equilibrium (steered_onto_reference); "
                    "%s, %d verified comparisons" % (record["status"],
                                                      record["counts"]["verified"])
                    ) in record["reason"]
