"""dR/du, assembled from a verified UMAT's tangent, against a difference of R.

The material tangent in these fixtures was verified in the UMAT repository:
the original routine ran in Abaqus, the OTI-converted build ran the same deck,
their stress and state histories agreed over the whole path, and the converted
build's DDSDDE matched a centred difference of the ORIGINAL's stress over a
plateau of step sizes. None of that says anything about what THIS repository
does with it. Integrating a correct tangent with the wrong quadrature weights,
the wrong Voigt order or the wrong DOF ordering produces a stiffness that is
not the derivative of the force the same code assembles, and every check that
compares the assembly against a hand-written form of the same assembly shares
whatever error is in both.

So this differences the thing itself. The residual is assembled at ``u``, at
``u + h e_j`` and at ``u - h e_j`` for every one of the 24 degrees of freedom,
and the centred difference is compared column by column with the assembled
stiffness.

Over a PLATEAU of step sizes, not at one. A centred difference has two error
terms that move in opposite directions -- truncation falling as ``h^2`` and
cancellation rising as ``eps/h`` -- so a single number taken anywhere on that
curve is indistinguishable from a number taken at the one step where two
different matrices happen to cross. What settles it is a run of step sizes
over which the disagreement stops moving. This is the same rule the UMAT
pipeline holds its own differences to, and it is held here for the same
reason.

Two derivatives are differenced, because the repository computes two and they
are not the same object:

* the small-strain material stiffness ``sum B0^T D B0 detJ0 w``, which is the
  derivative of ``sum B0^T sigma(u) detJ0 w`` when the stress follows the
  tangent; and
* the exact Jacobian of the finite-strain force at FIXED Cauchy stress, which
  is what ``force_tangent_fixed_sigma`` claims to be and is deliberately NOT
  the conventional geometric stiffness inside ``element_tangent``. That
  distinction is checked too, because the cheapest way to make a failing
  finite-difference pass is to quietly make one of them into the other.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from residual_core.formulations.c3d8_kernel import (  # noqa: E402
    ABAQUS_C3D8_GAUSS, b_matrix_reference, b_matrix_spatial,
    element_internal_force_finite_strain, element_internal_force_small_strain,
    element_tangent, force_tangent_fixed_sigma)
from residual_core.materials.verified_fixture import (  # noqa: E402
    tangent_convention)

sys.path.insert(0, str(Path(__file__).resolve().parent))
from verified_fixtures import (FIXTURES, all_fixtures,  # noqa: E402
                               hex_fixtures)

#: A unit cube, which is the geometry the corpus verification decks are
#: generated on, so the element being differenced here is the element the
#: numbers were produced in.
UNIT_CUBE = np.array([[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.],
                      [0., 0., 1.], [1., 0., 1.], [1., 1., 1.], [0., 1., 1.]])

#: Relative step sizes, spanning five decades. Enough that truncation
#: dominates one end and cancellation the other, so a plateau between them is
#: a plateau and not the whole range.
STEPS = (1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8)

#: How close a differenced column has to sit to the assembled one, inside the
#: plateau. Measured: the four frozen fixtures reach 1.5e-15 to 2.8e-15 on the
#: small-strain assembly and 3.0e-13 to 4.0e-13 on the finite-strain one, so
#: this leaves five decades of headroom and the test fails on a wrong assembly
#: rather than on a machine's rounding.
AGREES = 1e-7

#: How many step sizes have to agree, and over what span. One step cannot
#: separate a truncation error from a cancellation one; two adjacent ones
#: barely can. Three spanning a decade is the same bar the UMAT pipeline's
#: tangent verdict uses.
PLATEAU_STEPS = 3


def fixtures():
    """The hexahedra whose run passed every gate and that carry a tangent."""
    usable = [f for f in hex_fixtures()
              if f.converted and f.converted[-1].tangent is not None]
    if not usable:
        pytest.skip("no fixture carries a six-component tangent")
    return usable


def strain_at_each_point(displacement, reference=True, Xe=UNIT_CUBE):
    """B u at every integration point, in the fixture's own Voigt order."""
    out = []
    for point in ABAQUS_C3D8_GAUSS.points:
        if reference:
            B, _detJ = b_matrix_reference(Xe, point)
        else:                                      # pragma: no cover
            B, _detJ = b_matrix_spatial(Xe, point)
        out.append(B @ np.asarray(displacement, dtype=float))
    return out


def centred_jacobian(force, step: float, size: int = 24) -> np.ndarray:
    """d(force)/d(displacement) by a centred difference, column by column."""
    jacobian = np.zeros((size, size))
    for dof in range(size):
        plus, minus = np.zeros(size), np.zeros(size)
        plus[dof], minus[dof] = step, -step
        jacobian[:, dof] = (np.asarray(force(plus), dtype=float)
                            - np.asarray(force(minus), dtype=float)) / (2.0 * step)
    return jacobian


def sweep(analytic: np.ndarray, force, scale: float) -> list:
    """The disagreement at every step size, as the pipeline records one."""
    size = float(np.max(np.abs(analytic))) or 1.0
    rows = []
    for relative in STEPS:
        step = relative * scale
        numerical = centred_jacobian(force, step)
        rows.append({"step": relative,
                     "relative": float(np.max(np.abs(analytic - numerical))) / size})
    return rows


def plateau(rows: list, tolerance: float = AGREES) -> list:
    """The longest run of consecutive steps that all agree within tolerance."""
    best: list = []
    run: list = []
    for row in rows:
        if row["relative"] <= tolerance:
            run.append(row)
            best = max(best, run, key=len)
        else:
            run = []
    return best


def describe(rows: list) -> str:
    return ", ".join(f"h={r['step']:g}: {r['relative']:.3e}" for r in rows)


# ---------------------------------------------------------------------------
# the small-strain residual
# ---------------------------------------------------------------------------
def test_the_small_strain_stiffness_is_the_derivative_of_the_small_strain_force():
    """Differenced, not re-derived. ``sum B0^T D B0 detJ0 w`` against a centred
    difference of ``sum B0^T sigma(u) detJ0 w``, over every fixture that
    carries a six-component verified tangent.

    Measured on the four frozen fixtures: 2.06e-15, 2.83e-15, 2.03e-15 and
    1.52e-15 at their best step, and inside 1e-7 at all six step sizes from
    1e-3 to 1e-8.

    What that does and does not establish is worth being exact about. With the
    stress following the tangent, ``sum B0^T sigma(u)`` is LINEAR in ``u``, so
    a centred difference of it has no truncation error at all and the flat
    sweep is arithmetic rather than convergence. What it checks is therefore
    the thing that can actually be wrong here -- quadrature weights, Jacobian
    determinants, Voigt order, DOF ordering, the transpose -- since any of
    those makes the two matrices differ at EVERY step. The convergence
    question is asked where there is something to converge:
    ``test_the_fixed_stress_jacobian_is_the_derivative_of_the_finite_strain_force``
    and its displaced sibling, whose force is genuinely nonlinear in ``u``.
    """
    for fixture in fixtures():
        D = fixture.converted[-1].tangent
        base = np.asarray(fixture.original[-1].stress, dtype=float)
        K = element_tangent(UNIT_CUBE, np.zeros(24), np.tile(D, (8, 1, 1)),
                            mode="small")

        def force(displacement, D=D, base=base):
            stress = np.array([base + D @ strain for strain
                               in strain_at_each_point(displacement)])
            return element_internal_force_small_strain(UNIT_CUBE, stress)

        rows = sweep(K, force, scale=1.0)
        held = plateau(rows)
        assert len(held) >= PLATEAU_STEPS, (
            f"{fixture.path.name}: the assembled stiffness is not the "
            f"derivative of the assembled force over a plateau -- "
            f"{describe(rows)}")
        assert held[-1]["step"] <= held[0]["step"] / 10.0, (
            f"{fixture.path.name}: the agreement spans less than a decade of "
            f"step size, so it is a crossing rather than a plateau -- "
            f"{describe(rows)}")


def test_the_agreement_is_a_plateau_and_not_one_lucky_step():
    """A step far off the plateau disagrees, which is what makes the plateau
    evidence rather than a coincidence.

    Measured on the frozen fixtures: at a step of 1e-14 the centred difference
    of a stiffness of order 3.5e+04 loses its significant digits to
    cancellation and misses by 2.1e-04 to 2.8e-04, while the same comparison
    inside the plateau sits at 1.5e-15 to 2.8e-15 -- eleven orders apart. A
    difference that gave the same answer at every step would not be resolving
    anything, and its agreement would mean nothing.
    """
    fixture = fixtures()[0]
    D = fixture.converted[-1].tangent
    base = np.asarray(fixture.original[-1].stress, dtype=float)
    K = element_tangent(UNIT_CUBE, np.zeros(24), np.tile(D, (8, 1, 1)),
                        mode="small")

    def force(displacement):
        stress = np.array([base + D @ strain for strain
                           in strain_at_each_point(displacement)])
        return element_internal_force_small_strain(UNIT_CUBE, stress)

    size = float(np.max(np.abs(K)))
    cancelled = float(np.max(np.abs(K - centred_jacobian(force, 1e-14)))) / size
    on_plateau = min(row["relative"] for row in sweep(K, force, scale=1.0))
    assert cancelled > on_plateau * 1e3, (
        f"a step of 1e-14 gave {cancelled:.3e} and the plateau gives "
        f"{on_plateau:.3e}; if one step size is as good as any other, the "
        f"difference is not resolving anything and the agreement means "
        f"nothing")


# ---------------------------------------------------------------------------
# the finite-strain residual
# ---------------------------------------------------------------------------
def test_the_fixed_stress_jacobian_is_the_derivative_of_the_finite_strain_force():
    """The exact Jacobian of the code path as written, differenced.

    ``element_internal_force_finite_strain`` integrates over the CURRENT
    configuration, so its derivative at fixed Cauchy stress carries both the
    change of the spatial gradient and the change of the current volume.
    ``force_tangent_fixed_sigma`` claims to be that derivative; this
    differences the force to find out, using the Cauchy stress a verified UMAT
    actually returned rather than a stress invented to make it work.

    This force IS nonlinear in ``u`` -- the spatial gradient and the current
    volume both move with it -- so the sweep has a real error curve. Measured
    on the four fixtures: 3.0e-13 to 4.0e-13 at h=1e-3, rising linearly with
    1/h to 3.3e-08 to 4.6e-08 at h=1e-8, which is the cancellation term
    eps/h and is what a centred difference of a smooth function looks like.
    """
    for fixture in fixtures():
        stress = np.tile(np.asarray(fixture.original[-1].stress, dtype=float),
                         (8, 1))
        K = force_tangent_fixed_sigma(UNIT_CUBE, np.zeros(24), stress)

        def force(displacement, stress=stress):
            return element_internal_force_finite_strain(
                UNIT_CUBE, displacement, stress)

        rows = sweep(K, force, scale=1.0)
        held = plateau(rows)
        assert len(held) >= PLATEAU_STEPS, (
            f"{fixture.path.name}: the fixed-stress Jacobian is not the "
            f"derivative of the finite-strain force -- {describe(rows)}")


def test_the_finite_strain_derivative_holds_where_the_element_has_actually_moved():
    """Differenced about a displaced, sheared configuration rather than about
    ``u = 0``.

    At ``u = 0`` the current configuration is the reference one, so every term
    that depends on the geometry having moved is being evaluated at the one
    point where it cannot be wrong. This puts a 2% stretch and a 1% shear on
    the element first and differences there, which is where the spatial
    gradient, the current volume and the non-symmetry of the fixed-stress
    Jacobian all actually do something.

    Measured on the four frozen fixtures: the assembled Jacobian agrees with
    the centred difference to 2.7e-13 to 4.2e-13 relative at its best step and
    holds inside 1e-7 at all six, and the Jacobian at the displaced
    configuration differs from the one at ``u = 0`` by 1.1e-02 to 2.6e-02 --
    so the displacement is not a no-op and this is not the previous check
    under another name.
    """
    stretch = np.array([0.02, -0.006, -0.006])
    displaced = np.concatenate(
        [node * stretch + np.array([0.01 * node[1], 0.0, 0.0])
         for node in UNIT_CUBE])
    for fixture in fixtures():
        stress = np.tile(np.asarray(fixture.original[-1].stress, dtype=float),
                         (8, 1))
        K = force_tangent_fixed_sigma(UNIT_CUBE, displaced, stress)
        rest = force_tangent_fixed_sigma(UNIT_CUBE, np.zeros(24), stress)
        size = float(np.max(np.abs(K))) or 1.0
        assert float(np.max(np.abs(K - rest))) / size > 1e-3, (
            f"{fixture.path.name}: moving the element changed nothing, so "
            f"this is the undisplaced check under another name")

        def force(step, stress=stress):
            return element_internal_force_finite_strain(
                UNIT_CUBE, displaced + step, stress)

        rows = sweep(K, force, scale=1.0)
        held = plateau(rows)
        assert len(held) >= PLATEAU_STEPS, (
            f"{fixture.path.name}: the fixed-stress Jacobian is not the "
            f"derivative of the finite-strain force at a displaced "
            f"configuration -- {describe(rows)}")


def test_the_geometric_term_is_not_the_fixed_stress_jacobian():
    """They are different objects and the repository says so; this is the
    number behind the sentence.

    ``element_tangent``'s geometric term is the conventional symmetric
    initial-stress stiffness used with an objective stress rate;
    ``force_tangent_fixed_sigma``'s is the exact, generally non-symmetric
    derivative of the force routine as written. The cheapest way to make a
    failing finite-difference pass is to make one of them into the other, so
    the difference between them is asserted rather than left implicit.
    """
    for fixture in fixtures():
        stress = np.tile(np.asarray(fixture.original[-1].stress, dtype=float),
                         (8, 1))
        if float(np.max(np.abs(stress))) == 0.0:   # pragma: no cover
            continue
        exact = force_tangent_fixed_sigma(UNIT_CUBE, np.zeros(24), stress)
        D = fixture.converted[-1].tangent
        with_geometry = element_tangent(UNIT_CUBE, np.zeros(24),
                                        np.tile(D, (8, 1, 1)), sigma_ip=stress,
                                        mode="finite")
        material = element_tangent(UNIT_CUBE, np.zeros(24),
                                   np.tile(D, (8, 1, 1)), mode="small")
        geometric = with_geometry - material
        assert np.allclose(geometric, geometric.T, atol=1e-8 * max(
            float(np.max(np.abs(geometric))), 1.0)), (
            f"{fixture.path.name}: the geometric term is documented as "
            f"symmetric and is not")
        difference = float(np.max(np.abs(geometric - exact)))
        assert difference > 1e-12 * float(np.max(np.abs(exact))), (
            f"{fixture.path.name}: the geometric term and the fixed-stress "
            f"Jacobian have become the same matrix; they differ by exactly "
            f"the objective-rate terms and a finite difference of the force "
            f"routine as written targets the second")


# ---------------------------------------------------------------------------
# and what the fixtures are allowed to contain
# ---------------------------------------------------------------------------
def test_every_frozen_fixture_carries_a_completely_finite_history():
    """Asserted here as well as at freeze time, because a claim checked only
    by the tool that makes it is a claim nobody checked.

    Measured over the four frozen fixtures: 4 fixtures, 6 increments each on
    both sides, 48 histories of stress, state, strain, strain increment and
    tangent, every value finite.
    """
    for fixture in fixtures():
        for side in (fixture.original, fixture.converted):
            for record in side:
                for name, values in (("stress", record.stress),
                                     ("state", record.state),
                                     ("strain", record.strain),
                                     ("dstrain", record.dstrain)):
                    assert np.all(np.isfinite(values)), (
                        f"{fixture.path.name} increment {record.increment} "
                        f"{name}")
                if record.tangent is not None:
                    assert np.all(np.isfinite(record.tangent)), (
                        f"{fixture.path.name} increment {record.increment} "
                        f"tangent")


def test_a_fixture_with_a_nan_in_it_is_refused_rather_than_loaded(tmp_path: Path):
    """The one check nobody may skip. A NaN frozen into a baseline does not
    fail: it propagates, and every later comparison is made against it."""
    import json

    from residual_core.materials.verified_fixture import FixtureError, load

    payload = json.loads(fixtures()[0].path.read_text(encoding="utf-8"))
    payload["original"][2]["stress"][3] = float("nan")
    bad = tmp_path / "with_a_nan.json"
    bad.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(FixtureError) as raised:
        load(bad)
    assert "stress[3]" in str(raised.value)
    assert "finite" in str(raised.value)


def test_a_fixture_frozen_from_an_incomplete_run_is_refused(tmp_path: Path):
    """An increment short of a material point did not produce the state a
    later comparison would compare, so it may not become the thing later runs
    are compared against."""
    import json

    from residual_core.materials.verified_fixture import FixtureError, load

    payload = json.loads(fixtures()[0].path.read_text(encoding="utf-8"))
    payload["finite_history"]["history_grouping"]["original"][
        "first_incomplete_increment"] = {"increment": 6, "points": 7,
                                         "expected_points": 8, "time": 0.5}
    bad = tmp_path / "incomplete.json"
    bad.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(FixtureError) as raised:
        load(bad)
    assert "first incomplete increment is 6" in str(raised.value)


def test_a_fixture_with_no_finiteness_evidence_at_all_is_refused(tmp_path: Path):
    """A fixture that says nothing about the run it was frozen from is not a
    fixture whose run was fine; it is a fixture nobody can check."""
    import json

    from residual_core.materials.verified_fixture import FixtureError, load

    payload = json.loads(fixtures()[0].path.read_text(encoding="utf-8"))
    payload.pop("finite_history")
    bad = tmp_path / "no_evidence.json"
    bad.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(FixtureError) as raised:
        load(bad)
    assert "finite_history" in str(raised.value)


def test_the_fixtures_come_from_more_than_one_author():
    """One material passing says the bridge works for one material.

    The frozen set spans nine independently written sources from nine
    repositories: CAEAssistant-Group's isotropic elasticity,
    BristolCompositesInstitute's abaci UMAT, irfancn's finite-strain elastic
    and viscoelastic UMATs, AlexanderJFDR's and Sina-Taghizadeh's neo-Hookean
    models, keisuke58's viscous biofilm, mholla's isotropic-growth UMAT, and
    awhelanUCD's plane-stress Lemaitre damage model.
    """
    sources = {fixture.source_id for fixture in all_fixtures()}
    assert len(sources) >= 8, sorted(sources)
    repositories = {fixture.repository for fixture in all_fixtures()}
    assert len(repositories) >= 8, sorted(repositories)


def test_the_frozen_set_is_named_as_a_selection_rather_than_as_the_corpus():
    """Nine fixtures out of fifty-five verified cases, chosen for coverage.

    Measured on 2026-09-14 by exporting every one of the 55 cases that reached
    ``verified`` in the pass11 store (fingerprint b0d27ee53c630500; the
    exporter refused none of them) and putting each through ``diagnose`` with
    this repository's C3D8 assembly:

    * 1 carries NTENS=3 (CPS4 plane stress), which the C3D8 kernel cannot take;
    * of the remaining 54, **44 have no failing layer**, 6 fail at
      ``constitutive_derivative`` and 4 at ``transformation``;
    * 25 leave ``constitutive_derivative`` NOT ESTABLISHED rather than passing
      or failing it -- almost all of them growth models whose tangent's own
      prediction is thousands of times the stress change it should predict,
      because the stress answers an elastic strain the fixture does not carry;
    * ``state_sensitivity`` is not established for any of them, because no
      fixture carries a dR/dq to difference.

    So the committed set is a SELECTION for coverage and repository size, not
    what this repository can consume: it spans both kinematics, both tangent
    readings, a moving state, a rate-dependent material, a plane-stress case
    the kernel refuses, and one case whose two builds disagreed. The 35 Jeff97
    growth fixtures left out would add authors and not questions.
    """
    frozen = all_fixtures()
    assert len(frozen) == 9
    kinematics = {fixture.kinematics for fixture in frozen}
    assert kinematics >= {"small strain", "finite"}, (
        f"the set spans both kinematics the pipeline drives, not one: "
        f"{sorted(kinematics)}")
    readings = {tangent_convention(f).reading for f in frozen}
    assert {"material", "jaumann"} <= readings, (
        f"the set spans both readings of DDSDDE, so an assembler that picked "
        f"one would fail here: {sorted(str(r) for r in readings)}")
    assert any(f.ntens != 6 for f in frozen), (
        "a case the C3D8 kernel cannot take, so the refusal is exercised")
    assert any(not f.all_six_gates for f in frozen), (
        "a case the pipeline did not pass on all six gates, so the "
        "attribution is exercised against a known disagreement")
    for fixture in hex_fixtures():
        assert fixture.converted[-1].tangent is not None, fixture.path.name
        assert fixture.verification.get("states_checked"), fixture.path.name
