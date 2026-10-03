"""A verified corpus UMAT drives RA's assembled residual and its derivatives.

The materials are not written here. One is a corpus case the UMAT pipeline
verified in Abaqus (CAEAssistant's isotropic elasticity); the plastic one is a
labelled NON-corpus control (the 44 verified corpus cases contain no
rate-independent plasticity). Each is compiled live, twice from one source:
the original routine, and its PROPS-seeded OTI lift.

What is checked, and against what:

* dR/dp at fixed u_n and fixed incoming state, assembled from the OTI stress
  sensitivity, against a centred difference of R assembled from the ORIGINAL
  routine at p +/- h with the incoming state restored -- over a ladder of
  steps, verified only on a >=3-step plateau spanning a decade;
* du/dp through the whole load history (K du/dp = -dR/dp, history carried),
  against complete nonlinear re-solves of the path at p +/- h;
* the assembly itself: net force, primal parity of the two builds, K against a
  difference of R(u), patch test.

The pure-Python pieces (the source adaptation, the ladder adjudication and its
nonsmooth flag) are tested without any corpus data.
"""
from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from residual_core.corpus.provider import ProviderBuildFailed, adapt_source  # noqa: E402
from residual_core.corpus.verify import adjudicate, summarise  # noqa: E402

#: Materials are named by source id and the sha256 of the acquired file. Their
#: registry keys derive from the transform fingerprint and change at every freeze.
ELASTIC_SOURCE = ("CAEAssistant-Group__UMAT-Abaqus-Isotropic-Elasticity-Isothermal-Suboutine/ISOTROPIC-ELASTICITY.for",
                  "27cda337a45be53431149436fbad42bd19da0bb5ec2a68f6370c17b608326a86")


# --------------------------------------------------------------------------- #
# pure python
# --------------------------------------------------------------------------- #
FIXED_HEADER = """      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION STRESS(NTENS),
     4 JSTEP(4)
      STRESS(1)=PROPS(1)
      RETURN
      END
"""


@pytest.mark.unit
def test_an_unused_jstep_array_is_adapted_and_said_so():
    text, notes = adapt_source(FIXED_HEADER, fixed=True)
    assert "KSPT,KSTEP,KINC)" in text
    assert "JSTEPUNUSED(4)" in text
    assert "STRESS(1)=PROPS(1)" in text, "no executable statement may change"
    assert len(notes) == 1 and "declaration-only" in notes[0]


@pytest.mark.unit
def test_a_jstep_the_routine_reads_is_refused():
    used = FIXED_HEADER.replace("STRESS(1)=PROPS(1)", "STRESS(1)=PROPS(1)*JSTEP(2)")
    with pytest.raises(ProviderBuildFailed) as raised:
        adapt_source(used, fixed=True)
    assert raised.value.failure_class == "jstep_array_in_use"


@pytest.mark.unit
def test_a_source_without_jstep_is_untouched():
    plain = FIXED_HEADER.replace("JSTEP", "KSTEP").replace("KSTEP(4)", "KTEMP(4)")
    text, notes = adapt_source(plain, fixed=True)
    assert text == plain and notes == []


def _ladder(f, p, analytic):
    steps = {s: s * abs(p) for s in (1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 1e-7)}
    plus = {s: np.atleast_1d(f(p + h)) for s, h in steps.items()}
    minus = {s: np.atleast_1d(f(p - h)) for s, h in steps.items()}
    return adjudicate(np.atleast_1d(analytic), plus, minus, np.atleast_1d(f(p)), steps,
                      rtol=1e-6, atol=1e-12)


@pytest.mark.unit
def test_a_right_derivative_of_a_smooth_function_verifies_on_a_plateau():
    c = _ladder(np.sin, 0.7, np.cos(0.7))
    assert c["verdict"] == "verified"
    assert c["min_plateau_observed"] >= 3
    assert c["max_error_ratio"] <= 1.0


@pytest.mark.unit
def test_a_wrong_derivative_of_a_smooth_function_fails():
    c = _ladder(np.sin, 0.7, 1.001 * np.cos(0.7))
    assert c["verdict"] == "failed"
    assert summarise([c])["status"] == "failed"


@pytest.mark.unit
def test_a_kink_is_flagged_nonsmooth_rather_than_failed_or_verified():
    def kinked(x):        # slope 1 below 1, slope 3 above: a yield point at p
        return np.where(x < 1.0, x, 1.0 + 3.0 * (x - 1.0))
    c = _ladder(kinked, 1.0, 1.0)
    assert c["verdict"] == "nonsmooth"
    rolled = summarise([c])
    assert rolled["status"] == "unsupported" and rolled["counts"]["verified"] == 0


def _vector(analytic, f, p):
    from residual_core.corpus.verify import STEPS
    steps = {s: s * abs(p) for s in STEPS}
    plus = {s: f(p + h) for s, h in steps.items()}
    minus = {s: f(p - h) for s, h in steps.items()}
    return adjudicate(analytic, plus, minus, f(p), steps, rtol=1e-6, atol=1e-30)


def _column(x):
    # one large entry and two tiny ones (1e-7 and 1e-9 of the column)
    return np.array([np.sin(x), 1e-7 * np.cos(x), 1e-9 * np.exp(x)])


_TRUE = np.array([np.cos(0.7), -1e-7 * np.sin(0.7), 1e-9 * np.exp(0.7)])


@pytest.mark.unit
@pytest.mark.parametrize("entry,factor", [(1, 1.005), (1, 0.0), (2, 10.0), (0, 1.0 + 3e-6)])
def test_an_error_on_any_entry_fails_however_small_the_entry(entry, factor):
    """Vera B1 B/C: a column-norm tolerance passes +0.5% at 1e-7 of the column,
    x10 at 1e-9, -100% on a small entry; the entrywise rule fails all of them."""
    wrong = _TRUE.copy()
    wrong[entry] *= factor
    c = _vector(wrong, _column, 0.7)
    assert c["verdict"] == "failed", c["counts"]
    assert c["max_error_ratio"] > 1.0


@pytest.mark.unit
def test_the_right_column_verifies_entry_by_entry():
    c = _vector(_TRUE, _column, 0.7)
    assert c["verdict"] == "verified" and c["counts"]["fail"] == 0
    assert c["max_error_ratio"] <= 1.0 and c["min_plateau_observed"] >= 3


@pytest.mark.unit
def test_a_single_precision_routine_is_resolved_only_with_its_own_roundoff_model():
    """FD of a float32 computation cannot reach 1e-6 in double's round-off
    model (unresolved, never verified); with the float32 model the spread u_e
    enters the tolerance and a 1e-4 error is still caught."""
    from residual_core.corpus.verify import EPS32, STEPS
    f = lambda x: np.float64(np.float32(3.0 * x)) * np.array([1.0, 2.0])
    steps = {s: s * 0.7 for s in STEPS}
    plus = {s: f(0.7 + h) for s, h in steps.items()}
    minus = {s: f(0.7 - h) for s, h in steps.items()}
    right, wrong = np.array([3.0, 6.0]), np.array([3.0, 6.0 * (1 + 1e-4)])
    assert adjudicate(right, plus, minus, f(0.7), steps, rtol=1e-6, atol=1e-30)["verdict"] == "unresolved"
    ok = adjudicate(right, plus, minus, f(0.7), steps, rtol=1e-6, atol=1e-30, eps_eval=EPS32)
    assert ok["verdict"] == "verified" and ok["resolution_max"] < 1e-3
    bad = adjudicate(wrong, plus, minus, f(0.7), steps, rtol=1e-6, atol=1e-30, eps_eval=EPS32)
    assert bad["verdict"] == "failed"


@pytest.mark.unit
def test_a_derivative_the_difference_cannot_see_is_judged_against_the_roundoff_bound():
    f = lambda x: np.array([1.0 + 1e-17 * np.sin(x), np.sin(x)])
    assert _vector(np.array([0.0, np.cos(0.7)]), f, 0.7)["verdict"] == "verified"
    assert _vector(np.array([1e-6, np.cos(0.7)]), f, 0.7)["verdict"] == "failed"


@pytest.mark.unit
def test_partial_coverage_is_never_verified():
    from residual_core.corpus.verify import summarise as roll
    ok = {"verdict": "verified", "max_abs": 0.0, "max_rel": 0.0, "max_error_ratio": 0.1,
          "min_plateau_observed": 3}
    assert roll([ok, {"verdict": "unresolved", "max_abs": None, "max_rel": None}])["status"] \
        == "unsupported"
    assert roll([ok, {"verdict": "nonsmooth", "max_abs": None, "max_rel": None}, ok])["status"] \
        == "verified"
    assert roll([ok, {"verdict": "nonsmooth", "max_abs": None, "max_rel": None}])["status"] \
        == "unsupported"


# --------------------------------------------------------------------------- #
# live corpus UMAT
# --------------------------------------------------------------------------- #
def _live():
    if shutil.which("gfortran") is None:
        pytest.skip("gfortran not on PATH (environmental blocker)")
    from residual_core.corpus.sources import paths
    where = paths()
    if where.workspace is None:
        pytest.skip("CORPUS_WORKSPACE is unset; set it to the workspace folder that holds "
                    "final-umat/, corpus_run/ and discovery_cache/ (environmental blocker)")
    if not where.available():
        pytest.skip("corpus registry / verification records / acquisition cache not present "
                    "under CORPUS_WORKSPACE=%s" % where.workspace)
    try:
        import umat_oti.provider.build as build
    except ImportError:
        pytest.skip("umat_oti not importable; set PYTHONPATH=<final-umat>/src")
    if where.umat_repo.resolve() not in Path(build.__file__).resolve().parents:
        pytest.skip("umat_oti resolves to %s, not UMAT_OTI_REPO=%s (stale checkout)"
                    % (build.__file__, where.umat_repo))
    return where


@pytest.fixture(scope="module")
def out_root(tmp_path_factory):
    return tmp_path_factory.mktemp("corpus_residual")


@pytest.mark.integration
@pytest.mark.fortran
def test_a_verified_elastic_umat_gives_a_residual_derivative_a_difference_agrees_with(out_root):
    _live()
    from residual_core.corpus.runner import run_case
    from residual_core.corpus.sources import key_for_source, load_case
    case = load_case(key_for_source(*ELASTIC_SOURCE))
    assert case.terminal_state == "fully_verified"
    records = run_case(case, out_root, quick=True)
    by_feature = {r["feature"]: r for r in records}
    residual = by_feature["residual_sens"]
    assert residual["status"] == "verified", residual
    assert residual["counts"]["failed"] == 0 and residual["counts"]["nonsmooth"] == 0
    assert residual["counts"]["verified"] == residual["comparisons"] > 0
    assert residual["max_rel"] < 1e-6
    assert residual["reference"] == "fd" and "ORIGINAL" in residual["reference_text"]
    assert "u_n" in residual["held_fixed"] and residual["scope"] == "local"
    assert residual["build"]["kind"] == "provider" and residual["build"]["sha256"]
    assert residual["tolerance"] == 1.0 and residual["max_error"] <= 1.0
    assert residual["plateau_basis"] == "fd_only" and residual["min_plateau_observed"] >= 3
    assert by_feature["global_sens"]["status"] == "verified", by_feature["global_sens"]
    assembly = by_feature["assembly_consistency"]
    assert assembly["status"] == "verified", assembly["checks"]
    assert assembly["checks"]["primal_parity"] < 1e-12
    assert residual["h0_replays_bit_exact"] > 0
    assert Path(residual["evidence"]).is_file()
    # the derivative is not trivially zero: dR/dE carries the stress / E
    import json
    detail = json.loads(Path(residual["evidence"]).read_text())
    scales = [c["analytic_scale"] for c in detail["residual_sens"]["comparisons"]]
    assert max(scales) > 0.0


@pytest.mark.integration
@pytest.mark.fortran
@pytest.mark.slow
def test_a_plastic_control_through_load_unload_reload_holds_on_the_whole_path(out_root):
    _live()
    from residual_core.corpus import mesh as M
    from residual_core.corpus.runner import control_case, run_case
    case = control_case("m3_j2")
    problem = M.uniaxial(M.brick((1, 1, 1)), 0.01)          # load, unload, reload
    records = run_case(case, out_root, problems=[problem])
    by_feature = {r["feature"]: r for r in records}
    for feature in ("residual_sens", "global_sens", "assembly_consistency"):
        assert by_feature[feature]["status"] == "verified", by_feature[feature]
    activity = by_feature["global_sens"]["activity"]
    assert activity["increments_with_state_change"] >= 3, "plasticity must actually happen"
    assert activity["max_abs_dstate_dp"] > 0.0, "the history derivative must be carried"
    assert by_feature["global_sens"]["qoi_status"] == "verified"


# mholla growth, finite strain, 6 state variables
AREA_STRETCH_SOURCE = ("mholla__growth/umats/umat_area_stretch.f",
                       "869907d36a7efc5b7594c2fd35eb9c8e802ef839ac8749a6665e772fd23240e8")
# irfancn umat_elastic: integrates the incoming stress
HYPOELASTIC_SOURCE = ("irfancn__Abaqus-UMAT-elastic/umat_elastic.for",
                      "5e500764da04aaabba279113de51697b41641f62f5912ba5e44144e4fff31572")


@pytest.mark.integration
@pytest.mark.fortran
@pytest.mark.slow
def test_a_finite_strain_corpus_umat_with_state_verifies_under_the_abaqus_contract(out_root):
    """Finite strain WITH history (B2): Abaqus NLGEOM kinematics, state carried,
    local and total derivatives against the original, K, objectivity."""
    _live()
    from residual_core.corpus.runner import run_case
    from residual_core.corpus.sources import key_for_source, load_case
    case = load_case(key_for_source(*AREA_STRETCH_SOURCE))
    assert case.finite and case.nstatv > 0
    records = run_case(case, out_root, quick=True)
    by_feature = {r["feature"]: r for r in records}
    for feature in ("residual_sens", "global_sens", "assembly_consistency"):
        assert by_feature[feature]["status"] == "verified", by_feature[feature]
        assert by_feature[feature]["kinematics_contract"] == "nlgeom"
    assert by_feature["global_sens"]["activity"]["increments_with_state_change"] > 0
    assert by_feature["assembly_consistency"]["parts"]["objectivity"] is True
    assert by_feature["global_sens"]["scope"] == "total"


@pytest.mark.integration
@pytest.mark.fortran
@pytest.mark.slow
def test_the_hypoelastic_corpus_umat_is_objective_only_with_drot(out_root):
    _live()
    from residual_core.corpus import mesh as M
    from residual_core.corpus.provider import CorpusProvider, build_provider_for
    from residual_core.corpus.runner import objectivity_check
    from residual_core.corpus.sources import key_for_source, load_case, paths
    case = load_case(key_for_source(*HYPOELASTIC_SOURCE))
    record = build_provider_for(case, out_root / "providers" / case.key,
                                umat_repo=paths().umat_repo)
    provider = CorpusProvider(record, case, out_root / "providers" / case.key / "lib")
    problem = M.clamped_shear(M.brick((1, 1, 2), size=(1.0, 1.0, 2.0)), 0.1)
    problem.path = ((1.0, 2), (0.5, 1))
    assert objectivity_check(provider, case, problem, "nlgeom")["pass"]
    legacy = objectivity_check(provider, case, problem, "b1_legacy")
    assert not legacy["pass"] and legacy["sigma_QsigmaQT_max_rel"] > 1e-2


def test_zero_entries_are_reported_apart_from_the_plateau_statistics():
    """B3 review (Vera, LOW): a zero entry's all-steps 'plateau' is not evidence
    of FD convergence; the plateau minimum is taken over non-zero resolved
    entries and the zero entries are reported on their own."""
    from residual_core.corpus.verify import STEPS
    p = np.array([0.7, -1.3])
    f = lambda q: np.array([np.sin(q[0]) * 2.0, 0.0])        # entry 1 identically zero
    base = f(p)
    steps = {s: s for s in STEPS}
    plus = {s: f(p + np.array([s, 0.0])) for s in STEPS}
    minus = {s: f(p - np.array([s, 0.0])) for s in STEPS}
    c = adjudicate(np.array([2.0 * np.cos(0.7), 0.0]), plus, minus, base, steps,
                   rtol=1e-6, atol=1e-12)
    assert c["verdict"] == "verified"
    assert c["min_plateau_zero_entries"] == len(STEPS)
    assert 3 <= c["min_plateau_observed"] < len(STEPS)
