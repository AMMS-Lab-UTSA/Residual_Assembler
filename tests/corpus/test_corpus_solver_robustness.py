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
19. The FD re-solves at p +/- h are started from the NOMINAL solution of each
   increment (B12). Started cold, Newton at p +/- h did not converge, even at
   h = 1e-8, on paths whose nominal solve needed hundreds of backtracks: 113 of
   157 unresolved global records had missing ladder steps. Same equilibrium
   followed from the nominal one, same tolerance, same residual.
20. A slot whose value flows through an INTEGER variable has no derivative; the
   provider is built without it, and every record names the slot left out
   (B12: Worlthen enhanced/simplified curing, PROPS(5)).
21. The host speaks six components; a routine of NTENS = 4 (plane strain) reads
   and writes the first four, the others are zero (B12). The provider maps in
   both directions, for the original and for the OTI build.
22. The plant test (Vera, B12): every comparison is judged again with the analytic
   derivative scaled by 1 + 1e-4; a verified cell has power only if that is caught.
   The cell text carries "power-checked p of k", names the problems without power,
   and says "verified on 7 of 8 slots; PROPS(5) integer, non-differentiable,
   excluded" where a slot was left out.
18. Labels only (Vera B9): a record whose analytic solve was steered says so and
   what it came to; problems excused by a branch departure have their own
   reason class, distinct from fd_reference_unresolved; the per-source tally
   quotes "verified on k of n problems" and flags a cell that rests on 1 of 3.
"""
from __future__ import annotations

import json
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
UP_BODY_SOURCE = ("Jeff97__Programming-Plane-Strain-Plates-through-Growth-Under-Body-Forces/"
                  "Examples-In-Section-3/ArcUp/Th005/BodyForce-Growth-2Stages.for",
                  "1bf3ad8383f4b28f021b00667334cadd851f1f8b2b4b219c081a76e74f0719c1")
SCALLOP_SOURCE = ("Jeff97__General-shape-control-of-shell/Abaqus_Files/2Dto2D/From-2D-to-2D-Scallop.for",
                  "52ae0b5f7a057fed224b64ef0eb4739147b606332ea49432d3ed4577e7d1e7ab")
CURING_SOURCE = ("Worlthen__20220314-abqus-simulation/abaqus/enhanced/enhanced_curing.for",
                 "6f0ff5e8b246f8c508cb2924b69e99c7fa290b0c745a9679be4aa145811b7be9")
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


def test_the_cell_text_carries_power_and_the_slot_coverage():
    import importlib.util
    spec = importlib.util.spec_from_file_location("run_corpus_residual",
                                                  ROOT / "tools" / "run_corpus_residual.py")
    tool = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(tool)
    wrt = "PROPS slots [1, 2, 3, 4, 6, 7, 8] (seeded; NOT seeded, non-differentiable: [5])"
    unseeded = [{"props_index": 5, "reason": "non_differentiable_integer_parameter_path"}]

    def rec(problem, status, power, **extra):
        return dict({"key": "k" * 24, "source_id": "src", "feature": "global_sens",
                     "problem": problem, "status": status, "plant": {"power": power}}, **extra)
    full = tool.fold_cell([rec("a", "verified", True), rec("b", "verified", True),
                           rec("c", "verified", True)])
    assert full["power_checked"] and "power-checked on 3 of 3" in full["text"]
    thin = tool.fold_cell([rec("a", "verified", True), rec("b", "verified", False),
                           rec("c", "unsupported", False)])
    assert thin["cell"] == "verified" and not thin["power_checked"]
    assert thin["unpowered_problems"] == ["b"] and "NOT POWER-CHECKED" in thin["text"]
    # a cell whose only verification has no power says so in the text itself
    blind = tool.fold_cell([rec("a", "verified", False), rec("b", "unsupported", False),
                            rec("c", "unsupported", False)])
    assert "NOT POWER-CHECKED (no power at 1e-4 on: a)" in blind["text"]
    slots = tool.fold_cell([rec(p, "verified", True, wrt=wrt, unseeded_slots=unseeded)
                            for p in ("a", "b", "c")])
    assert slots["slots"] == ("verified on 7 of 8 slots; PROPS(5) integer, "
                              "non-differentiable, excluded")
    assert slots["slots"] in slots["text"]


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

    def __init__(self, p, offsets, state_offsets=(0.0, 0.0), gain=1.0):
        from types import SimpleNamespace
        a = gain * np.array([[1.0, 2.0], [3.0, 4.0]])
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
            du_scale=1.0, gain=1.0):
    from types import SimpleNamespace
    from residual_core.corpus import runner
    source = tmp_path / "u.for"
    source.write_text("      EMOD=PROPS(1)\n")
    case = _case()
    case.source_path, case.props = source, [2.0]
    monkeypatch.setattr(runner, "reference_solve",
                        lambda prov, prob, props, kin, sched, first=None, **kw:
                        _Branch(props[0], [0.0, 0.0], gain=gain))
    monkeypatch.setattr(runner, "run_history", lambda *a, **k: _Branch(2.0, [0.0, 0.0], gain=gain))
    mine = _Branch(2.0, analytic_offsets, state_offsets, gain=gain)
    analytic = SimpleNamespace(schedule=[{}, {}], newton_matrix_arg="exact", V=mine.V,
                               incoming=mine.incoming, primal_parity=[parity, parity],
                               du_dp=[du_scale * gain * np.array([[1.0], [2.0]]),
                                      du_scale * gain * np.array([[3.0], [4.0]])],
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


def test_every_perturbed_resolve_is_started_from_the_nominal_solution(monkeypatch, tmp_path):
    from types import SimpleNamespace
    from residual_core.corpus import runner
    calls = []
    source = tmp_path / "u.for"
    source.write_text("      EMOD=PROPS(1)\n")
    case = _case()
    case.source_path, case.props = source, [2.0]
    nominal = _Branch(2.0, [0.0, 0.0])

    def fake(prov, prob, props, kin, sched, first=None, **kw):
        calls.append((props[0], kw.get("start_from")))
        return _Branch(props[0], [0.0, 0.0])
    monkeypatch.setattr(runner, "reference_solve", fake)
    monkeypatch.setattr(runner, "run_history", lambda *a, **k: _Branch(2.0, [0.0, 0.0]))
    analytic = SimpleNamespace(schedule=[{}, {}], newton_matrix_arg="exact", V=nominal.V,
                               incoming=nominal.incoming, primal_parity=[0.0, 0.0],
                               du_dp=[np.array([[1.0], [2.0]]), np.array([[3.0], [4.0]])],
                               qoi_dp=lambda problem: np.zeros((2, 0, 1)))
    runner.global_sensitivity(SimpleNamespace(slots=[1]), case, SimpleNamespace(qois=[]),
                              analytic, None)
    nominal_call, perturbed = calls[0], calls[1:]
    assert nominal_call[1] is None and len(perturbed) >= 2
    assert all(start is not None and all(np.array_equal(a, b) for a, b in zip(start, nominal.V))
               for value, start in perturbed)


def test_a_planted_error_is_caught_where_the_derivative_is_resolved(monkeypatch, tmp_path):
    """The plant test (Vera, B12): the analytic scaled by 1 + 1e-4 must FAIL."""
    from residual_core.corpus import runner
    by_inc, info = _global(monkeypatch, tmp_path, [0.0, 0.0])
    assert {c["verdict"] for c in by_inc.values()} == {"verified"}
    plant = runner.plant_summary(info["planted"])
    assert plant["power"] and plant["status"] == "failed"
    assert plant["factor"] == pytest.approx(1.0001) and plant["u"]["caught"] == plant["u"]["comparisons"] > 0


def test_a_verdict_on_a_derivative_the_fd_cannot_see_has_no_power(monkeypatch, tmp_path):
    """A response that does not depend on the parameter: the analytic derivative is
    zero, verified trivially, and scaling zero plants nothing. No power."""
    from residual_core.corpus import runner
    by_inc, info = _global(monkeypatch, tmp_path, [0.0, 0.0], gain=0.0)
    assert {c["verdict"] for c in by_inc.values()} == {"verified"}
    plant = runner.plant_summary(info["planted"])
    assert not plant["power"] and plant["caught"] == 0 and plant["status"] != "failed"


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
                         ("arc_down", ARC_DOWN_SOURCE), ("z2", Z2_SOURCE),
                         ("up_body", UP_BODY_SOURCE), ("scallop", SCALLOP_SOURCE),
                         ("curing", CURING_SOURCE)):
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


def _resolve_pair(built, name, h):
    """Nominal solve, and the solve at p(1+h) of slot 1 cold and from the nominal solution."""
    from residual_core.corpus import runner
    from residual_core.corpus.engine import NewtonFailed
    case, provider, _, _ = built[name]
    problem = [p for p in runner.default_problems(case) if p.name.startswith("clamped_tension")][0]
    analytic = runner.analytic_solve(provider, problem, case.props, None, {})
    if analytic.cutbacks:
        analytic = runner.common_schedule(provider, problem, case.props, None, analytic, {})
    first = getattr(analytic, "reference_steering", None) or \
        runner._SAME_STEERING.get(analytic.newton_matrix_arg)
    props = np.asarray(case.props, dtype=float)
    nominal = runner.reference_solve(provider, problem, props, None, analytic.schedule, first=first)
    perturbed = props.copy()
    perturbed[provider.slots[0] - 1] *= 1.0 + h
    cold = warm = None
    try:
        cold = runner.reference_solve(provider, problem, perturbed, None, analytic.schedule,
                                      first=first)
    except NewtonFailed:
        pass
    warm = runner.reference_solve(provider, problem, perturbed, None, analytic.schedule,
                                  first=first, start_from=nominal.V)
    return nominal, cold, warm


@pytest.mark.integration
@pytest.mark.fortran
@pytest.mark.slow
def test_a_resolve_that_does_not_converge_cold_converges_from_the_nominal_solution(built):
    """Scallop, clamped tension: the nominal solve needs 500+ backtracks; at
    p(1 + 1e-6) a cold Newton does not converge. From the nominal solution it does."""
    nominal, cold, warm = _resolve_pair(built, "scallop", 1e-6)
    assert sum(nominal.backtracks) > 100
    assert cold is None, "this problem is the one a cold re-solve cannot do"
    assert len(warm.V) == len(nominal.V)


@pytest.mark.integration
@pytest.mark.fortran
@pytest.mark.slow
def test_where_both_converge_the_warm_resolve_is_the_cold_one(built):
    """Starting from the nominal solution changes the basin Newton starts in, not
    the equilibrium it finds: where the cold re-solve converges too, the two
    agree to round-off (ArcUp Th005 BodyForce, clamped tension, h = 1e-6)."""
    nominal, cold, warm = _resolve_pair(built, "up_body", 1e-6)
    assert cold is not None
    scale = float(np.max(np.abs(np.array(nominal.V))))
    worst = max(float(np.max(np.abs(a - b))) for a, b in zip(cold.V, warm.V)) / scale
    assert worst < 1e-9
    assert max(float(np.max(np.abs(a - b))) for a, b in zip(warm.V, nominal.V)) / scale > 1e-12, \
        "the perturbation must move the solution; otherwise nothing was compared"


def test_a_provider_is_built_without_the_slot_that_flows_through_an_integer(monkeypatch, tmp_path):
    from residual_core.corpus import provider as P
    source = tmp_path / "u.for"
    source.write_text("      SUBROUTINE UMAT\n      END\n")
    case = _case()
    case.source_path, case.props = source, [1.0, 2.0, 3.0, 4.0, 5.0]
    seen = []

    class Module:
        __file__ = str(tmp_path / "build.py")

        @staticmethod
        def build_provider(contract_path, build_dir):
            contract = json.loads(Path(contract_path).read_text())
            slots = [p["props_index"] for p in contract["parameters"]]
            seen.append(slots)
            if 5 in slots:
                raise RuntimeError("non_differentiable_integer_parameter_path: PROPS(5) flows "
                                   "through INTEGER variable K; integer conversion is "
                                   "non-differentiable.")
            build_dir = Path(build_dir)
            build_dir.mkdir(parents=True, exist_ok=True)
            obj = build_dir / "x.obj"
            obj.write_bytes(b"obj")
            out = build_dir / "x.json"
            out.write_text("{}")
            return {"object": str(obj), "contract": str(out)}
    monkeypatch.setattr(P, "_umat_oti_build_module", lambda repo: Module)
    monkeypatch.setattr(P, "_umat_tree_state", lambda f: {})
    record = P.build_provider_for(case, tmp_path / "out", umat_repo=tmp_path)
    assert seen == [[1, 2, 3, 4, 5], [1, 2, 3, 4]]
    assert [p["props_index"] for p in record["parameters"]] == [1, 2, 3, 4]
    assert [u["props_index"] for u in record["unseeded_slots"]] == [5]
    assert record["unseeded_slots"][0]["reason"] == "non_differentiable_integer_parameter_path"
    # any other build failure is still raised, and the last slot is never dropped
    class Other(Module):
        @staticmethod
        def build_provider(contract_path, build_dir):
            raise RuntimeError("non_differentiable_integer_parameter_path: PROPS(1) flows ...")
    case.props = [1.0]
    monkeypatch.setattr(P, "_umat_oti_build_module", lambda repo: Other)
    with pytest.raises(P.ProviderBuildFailed):
        P.build_provider_for(case, tmp_path / "out2", umat_repo=tmp_path)


@pytest.mark.integration
@pytest.mark.fortran
@pytest.mark.slow
def test_a_source_with_an_integer_slot_verifies_on_the_differentiable_ones(built, tmp_path):
    from residual_core.corpus import runner
    case, provider, record, _ = built["curing"]
    assert [u["props_index"] for u in record["unseeded_slots"]] == [5]
    assert 5 not in provider.slots and len(provider.slots) == len(case.props) - 1
    problem = runner.default_problems(case, quick=True)[0]
    records = runner.run_case(case, tmp_path, problems=[problem],
                              features=("residual_sens", "global_sens"))
    for r in [r for r in records if r["feature"] in ("residual_sens", "global_sens")]:
        assert r["status"] == "verified", r
        assert "NOT seeded, non-differentiable: [5]" in r["wrt"]
        assert r["unseeded_slots"][0]["props_index"] == 5


# --------------------------------------------------------------------------- #
# NTENS-generic provider: the host's six components, the routine's NTENS
# --------------------------------------------------------------------------- #
ELASTIC_UMAT = """      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS),
     1 DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),TIME(2),
     2 PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3),DROT(3,3),
     3 DFGRD0(3,3),DFGRD1(3,3)
      E=PROPS(1)
      XNU=PROPS(2)
      XLAM=E*XNU/((1.D0+XNU)*(1.D0-2.D0*XNU))
      G=E/(2.D0*(1.D0+XNU))
      DO I=1,NTENS
        DO J=1,NTENS
          DDSDDE(I,J)=0.D0
        ENDDO
      ENDDO
      DO I=1,NDI
        DO J=1,NDI
          DDSDDE(I,J)=XLAM
        ENDDO
        DDSDDE(I,I)=XLAM+2.D0*G
      ENDDO
      DO I=NDI+1,NTENS
        DDSDDE(I,I)=G
      ENDDO
      DO I=1,NTENS
        DO J=1,NTENS
          STRESS(I)=STRESS(I)+DDSDDE(I,J)*DSTRAN(J)
        ENDDO
      ENDDO
      RETURN
      END
"""


def _regular_only_provider(tmp_path, ntens):
    import subprocess
    from residual_core.corpus.hidden_state import _REGULAR_SHIM, _ProbeProvider
    (tmp_path / "umat.f").write_text(ELASTIC_UMAT)
    (tmp_path / "shim.f90").write_text(_REGULAR_SHIM)
    run = dict(cwd=tmp_path, check=True, capture_output=True, text=True)
    subprocess.run(["gfortran", "-O0", "-fPIC", "-std=legacy", "-ffixed-form",
                    "-ffixed-line-length-none", "-c", "umat.f", "-o", "umat.o"], **run)
    subprocess.run(["gfortran", "-O0", "-fPIC", "-ffree-form", "-c", "shim.f90", "-o", "shim.o"], **run)
    subprocess.run(["gfortran", "-shared", "shim.o", "umat.o", "-o", "lib.so"], **run)
    case = _case()
    case.ntens, case.nstatv, case.props = ntens, 1, [210.0e3, 0.3]
    record = {"parameters": [{"props_index": 1}, {"props_index": 2}]}
    return _ProbeProvider(record, case, tmp_path, tmp_path / "lib.so"), case


def test_components_map_between_the_host_and_a_routine_of_fewer():
    from residual_core.corpus.provider import ProviderBuildFailed, to_host, to_material
    host = np.arange(1.0, 7.0)
    assert list(to_material(host, 4, 0)) == [1, 2, 3, 4]
    assert list(to_material(host, 6, 0)) == [1, 2, 3, 4, 5, 6]
    assert list(to_host(np.array([1.0, 2.0, 3.0, 4.0]), 4, 0)) == [1, 2, 3, 4, 0, 0]
    matrix = np.arange(16.0).reshape(1, 4, 4)
    padded = to_host(to_host(matrix, 4, 1), 4, 2)
    assert padded.shape == (1, 6, 6) and np.array_equal(padded[0, :4, :4], matrix[0])
    assert not padded[0, 4:, :].any() and not padded[0, :, 4:].any()
    for unsupported in (3, 2, 5):
        with pytest.raises(ProviderBuildFailed):
            to_material(host, unsupported, 0)


@pytest.mark.integration
@pytest.mark.fortran
def test_a_plane_strain_routine_sees_four_components_and_returns_the_host_six(tmp_path):
    """The same elastic law as NTENS = 6 and NTENS = 4: the leading components and
    the leading block of the tangent agree; the shear components a plane-strain
    routine does not have are zero for the host."""
    for name in ("d3", "pe"):
        (tmp_path / name).mkdir()
    three_d, _ = _regular_only_provider(tmp_path / "d3", 6)
    plane, _ = _regular_only_provider(tmp_path / "pe", 4)
    rng = np.random.default_rng(0)
    n = 3
    dstran = rng.normal(size=(n, 6)) * 1e-3
    flat = dstran.copy()
    flat[:, 4:] = 0.0                      # the 3D routine must see a plane-strain increment
    common = dict(time=np.zeros(2), dtime=1.0, coords=np.zeros((n, 3)), celent=np.ones(n),
                  noel=np.arange(1, n + 1), npt=np.ones(n, dtype=int))
    args = (np.array([210.0e3, 0.3]), np.zeros((n, 6)), np.zeros((n, 1)), np.zeros((n, 6)))
    eye = np.tile(np.eye(3), (n, 1, 1))
    a = three_d.regular(*args, flat, eye, eye, **common)
    b = plane.regular(*args, dstran, eye, eye, **common)    # e13, e23 present: not seen
    assert np.allclose(b["stress"][:, :4], a["stress"][:, :4], rtol=1e-14, atol=0.0)
    assert not b["stress"][:, 4:].any()
    assert np.allclose(b["ddsdde"][:, :4, :4], a["ddsdde"][:, :4, :4], rtol=1e-14, atol=0.0)
    assert not b["ddsdde"][:, 4:, :].any() and not b["ddsdde"][:, :, 4:].any()
    assert plane.ntens == three_d.ntens == 6 and plane.material_ntens == 4


def test_the_oti_call_maps_seeds_in_and_results_out(monkeypatch):
    """The OTI entry with NTENS = 4: seeds reach the library with four components and
    the outputs come back padded to six. A stand-in library records and fills them."""
    from residual_core.corpus import provider as P
    seen = {}

    class Lib:
        def corpus_total(self, *args):
            n, nprops, nt, ns, npar = (a.value for a in args[:5])
            arrays = [a for a in args if isinstance(a, np.ndarray)]
            dsi = arrays[6]
            seen.update(nt=nt, dsi_shape=dsi.shape, stress_shape=arrays[1].shape)
            ddsdde, dsig, dstv, dstv_de, pnewdt = arrays[-5:]
            ddsdde[:] = 1.0
            dsig[:] = 2.0
            dstv[:] = 3.0
            dstv_de[:] = 4.0
            pnewdt[:] = 1.0
    monkeypatch.setattr(P, "_ptr", lambda a: a)
    provider = P.CorpusProvider.__new__(P.CorpusProvider)
    provider.material_ntens, provider.ntens, provider.nstatv, provider.nparam = 4, 6, 2, 3
    provider.nprops, provider.finite, provider.lib = 2, True, Lib()
    provider.case = _case()
    n = 2
    out = provider.total(np.array([1.0, 0.3]), np.zeros((n, 6)), np.zeros((n, 2)), np.zeros((n, 6)),
                         np.zeros((n, 6)), np.tile(np.eye(3), (n, 1, 1)), np.tile(np.eye(3), (n, 1, 1)),
                         time=np.zeros(2), dtime=1.0, coords=np.zeros((n, 3)), celent=np.ones(n),
                         noel=np.arange(1, n + 1), npt=np.ones(n, dtype=int),
                         dstress_in=np.ones((n, 6, 3)))
    assert seen["nt"] == 4 and seen["dsi_shape"] == (n, 3, 4) and seen["stress_shape"] == (n, 4)
    assert out["stress"].shape == (n, 6) and out["ddsdde"].shape == (n, 6, 6)
    assert out["dstress_dp"].shape == (n, 6, 3) and out["dstate_ddstran"].shape == (n, 2, 6)
    assert (out["dstress_dp"][:, :4] == 2.0).all() and not out["dstress_dp"][:, 4:].any()
    assert (out["ddsdde"][:, :4, :4] == 1.0).all() and not out["ddsdde"][:, 4:].any()
    assert (out["dstate_ddstran"][:, :, :4] == 4.0).all() and not out["dstate_ddstran"][:, :, 4:].any()


@pytest.mark.integration
@pytest.mark.fortran
@pytest.mark.slow
def test_a_live_verified_record_carries_a_plant_that_is_caught(built):
    from residual_core.corpus import runner
    case, _, _, root = built["flat"]
    problem = runner.default_problems(case, quick=True)[0]
    records = runner.run_case(case, root / "run_plant", problems=[problem],
                              features=("residual_sens", "global_sens"))
    for record in [r for r in records if r["feature"] in ("residual_sens", "global_sens")]:
        assert record["status"] == "verified"
        plant = record["plant"]
        assert plant["power"] and plant["status"] == "failed" and plant["caught"] > 0
        assert plant["factor"] == pytest.approx(1.0 + 1e-4)
