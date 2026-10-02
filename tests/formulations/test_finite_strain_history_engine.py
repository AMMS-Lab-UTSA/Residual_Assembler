"""Finite strain WITH history through RA's residual assembly, end to end.

A non-corpus UMAT written here: isotropic HYPOELASTIC (it integrates the
pre-rotated stress, so it is objective only if the host applies DROT exactly
as Abaqus does) with a softening state variable

    STATEV(1) <- STATEV(1) + DSTRAN : DSTRAN        (rotation invariant)
    STRESS    <- STRESS + C(E, nu) : DSTRAN / (1 + ALPHA * STATEV(1))

so the response depends on its whole history and d STATEV/d p is non-zero.
Compiled live twice from one source (original + PROPS-seeded OTI lift) by
umat_oti.provider; every reference below is the ORIGINAL routine.

Checked, under the Abaqus NLGEOM contract (``nlgeom``: mean-dilatation Fbar,
Hughes-Winget DROT, pre-rotated STRESS/STRAN, B-bar residual):

(a) K vs a centred difference of R(u) at the end of every path segment;
(b) Newton converges quadratically with that K;
(c) dR/dp at fixed u_n and incoming state vs FD of R from the original;
(d) total du/dp and dQoI/dp through the history vs complete re-solves at p+/-h;
(e) a superposed rigid rotation leaves frame displacements, reactions (and
    their magnitudes) and du/dp unchanged and gives sigma' = Q sigma Q^T --
    while the B1 convention (DROT = I, no pre-rotation) visibly does not.
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

UMAT = """      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
      INCLUDE 'ABA_PARAM.INC'
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
      EMOD = PROPS(1)
      ENU = PROPS(2)
      ALPHA = PROPS(3)
      ALAM = EMOD*ENU/((1.D0+ENU)*(1.D0-2.D0*ENU))
      AMU = EMOD/(2.D0*(1.D0+ENU))
      DSQ = 0.D0
      DO I = 1, NDI
        DSQ = DSQ + DSTRAN(I)**2
      END DO
      DO I = NDI+1, NTENS
        DSQ = DSQ + 0.5D0*DSTRAN(I)**2
      END DO
      STATEV(1) = STATEV(1) + DSQ
      FAC = 1.D0/(1.D0 + ALPHA*STATEV(1))
      DO I = 1, NTENS
        DO J = 1, NTENS
          DDSDDE(I,J) = 0.D0
        END DO
      END DO
      DO I = 1, NDI
        DO J = 1, NDI
          DDSDDE(I,J) = FAC*ALAM
        END DO
        DDSDDE(I,I) = FAC*(ALAM + 2.D0*AMU)
      END DO
      DO I = NDI+1, NTENS
        DDSDDE(I,I) = FAC*AMU
      END DO
      DO I = 1, NTENS
        DO J = 1, NTENS
          STRESS(I) = STRESS(I) + DDSDDE(I,J)*DSTRAN(J)
        END DO
      END DO
      RETURN
      END
"""

PROPS = [1000.0, 0.3, 40.0]


def _toolchain():
    if shutil.which("gfortran") is None:
        pytest.skip("gfortran not on PATH (environmental blocker)")
    try:
        import umat_oti.provider.build  # noqa: F401
    except ImportError:
        pytest.skip("umat_oti not importable; set PYTHONPATH=<final-umat>/src")


@pytest.fixture(scope="module")
def provider(tmp_path_factory):
    _toolchain()
    from residual_core.corpus.provider import CorpusProvider, build_provider_for
    from residual_core.corpus.sources import CorpusCase
    root = tmp_path_factory.mktemp("finite_history")
    source = root / "hypo_softening.f"
    source.write_text(UMAT, encoding="utf-8")
    case = CorpusCase(key="test_hypo_softening", source_id="tests/formulations (not corpus)",
                      source_path=source, source_form="fixed", family="test",
                      terminal_state="test", verification_fingerprint="",
                      kinematics="finite", ntens=6, nstatv=1, props=list(PROPS),
                      material_provenance="test", element_type="C3D8")
    record = build_provider_for(case, root / "provider")
    return CorpusProvider(record, case, root / "provider" / "lib"), case


def _problem():
    from residual_core.corpus import mesh as M
    problem = M.clamped_shear(M.brick((1, 1, 2), size=(1.0, 1.0, 2.0)), 0.15)
    problem.path = ((1.0, 3), (0.4, 2), (1.2, 2))
    return problem


@pytest.fixture(scope="module")
def analytic(provider):
    from residual_core.corpus.engine import run_history
    prov, case = provider
    return run_history(prov, _problem(), case.props, material="oti", sensitivities=True,
                       kinematics="nlgeom")


@pytest.mark.integration
@pytest.mark.fortran
def test_history_is_carried_and_differentiated(provider, analytic):
    moved = [float(np.max(np.abs(o.state - i.state)))
             for i, o in zip(analytic.incoming, analytic.outgoing)]
    assert min(moved) > 0.0, "the state must evolve in every increment"
    assert max(float(np.max(np.abs(d))) for d in analytic.dstate_dp) > 0.0
    assert max(analytic.primal_parity) <= 1e-12
    assert max(analytic.sensitivity_equilibrium) <= 1e-10


@pytest.mark.integration
@pytest.mark.fortran
def test_K_is_the_derivative_of_R_at_converged_states(provider, analytic):
    from residual_core.corpus.runner import tangent_check
    prov, case = provider
    for t in tangent_check(prov, case, _problem(), analytic, "nlgeom"):
        assert t["K_exact_oti"]["verdict"] == "verified", t["K_exact_oti"]["counts"]
        assert t["K_exact_oti"]["max_error_ratio"] <= 1.0


@pytest.mark.integration
@pytest.mark.fortran
def test_newton_converges_quadratically(analytic):
    from residual_core.corpus.runner import newton_rates
    rates = newton_rates(analytic)
    assert rates["pass"] and rates["median_order"] >= 1.8, rates


@pytest.mark.integration
@pytest.mark.fortran
def test_the_local_residual_sensitivity_matches_a_difference_of_the_original(provider, analytic):
    from residual_core.corpus.runner import residual_sensitivity
    from residual_core.corpus.verify import summarise
    prov, case = provider
    comps, info = residual_sensitivity(prov, case, _problem(), analytic, "nlgeom")
    rolled = summarise(comps)
    assert rolled["status"] == "verified", rolled
    assert rolled["max_error_ratio"] <= 1.0 and rolled["min_plateau_observed"] >= 3
    assert info["h0_replays_bit_exact"] == len(comps)


@pytest.mark.integration
@pytest.mark.fortran
def test_the_total_sensitivity_matches_resolves_of_the_whole_history(provider, analytic):
    from residual_core.corpus.runner import global_sensitivity
    from residual_core.corpus.verify import summarise
    prov, case = provider
    comps, info = global_sensitivity(prov, case, _problem(), analytic, "nlgeom")
    rolled = summarise(comps)
    assert rolled["status"] == "verified", rolled
    assert info["nominal_rerun_bit_exact"] and not info["reference_resolve_failures"]
    qoi = summarise([c for c in comps if c.get("quantity") != "u"])
    assert qoi["status"] == "verified"


@pytest.mark.integration
@pytest.mark.fortran
def test_a_superposed_rotation_changes_nothing_under_the_abaqus_contract(provider):
    from residual_core.corpus.runner import objectivity_check
    prov, case = provider
    result = objectivity_check(prov, case, _problem(), "nlgeom")
    assert result["pass"], result
    assert result["reaction_magnitude_max_rel"] <= 1e-9
    assert result["dV_dp_max_rel"] <= 1e-9


@pytest.mark.integration
@pytest.mark.fortran
def test_without_drot_the_hypoelastic_law_is_not_objective(provider):
    """The B1 convention (DROT = I, stress not pre-rotated) is what this batch
    replaces: for a law that integrates the incoming stress it is wrong at O(1)."""
    from residual_core.corpus.runner import objectivity_check
    prov, case = provider
    result = objectivity_check(prov, case, _problem(), "b1_legacy")
    assert not result["pass"]
    assert result["sigma_QsigmaQT_max_rel"] > 1e-2
