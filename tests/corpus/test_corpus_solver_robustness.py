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
                         ("arc_down", ARC_DOWN_SOURCE)):
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
    problem = runner.default_problems(case, quick=True)[0]
    records = runner.run_case(case, root / "run", problems=[problem], features=("global_sens",))
    (record,) = [r for r in records if r["feature"] == "global_sens"]
    assert record["status"] == "unsupported"
    assert record["reason_class"] == "reference_solve_failed"
