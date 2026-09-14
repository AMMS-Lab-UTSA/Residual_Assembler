"""Which step a boundary condition belongs to, and why it is not cosmetic.

Abaqus does not apply a deck's boundary blocks all at once. Conditions before
the first ``*STEP`` are initial; each step then carries the previous set
forward and modifies it (``OP=MOD``, the default) or replaces it outright
(``OP=NEW``). A parser that flattens the steps into one list keeps every
condition and loses the only thing that says when each applies -- and because
the result is a dictionary of prescribed values, the last block silently wins.

That is not hypothetical here. Every verification deck this project generates
carries four steps -- uniaxial, simple shear, uniaxial reversed, uniaxial held
-- each opening with ``*BOUNDARY, OP=NEW``. Read flat, a deck's FIRST step is
driven by its FOURTH step's displacements, and the strain that comes out has
nothing to do with the increment anybody asked about.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from residual_core.core import constraints  # noqa: E402
from residual_core.core.dof_manager import DofManager  # noqa: E402
from residual_core.core.model import from_abaqus  # noqa: E402
from residual_core.io.abaqus_inp_parser import parse_inp  # noqa: E402
from verified_fixtures import hex_fixtures  # noqa: E402


@pytest.fixture()
def parsed(tmp_path):
    fixture = hex_fixtures()[0]
    path = tmp_path / "deck.inp"
    path.write_text(fixture.deck, encoding="utf-8")
    return parse_inp(str(path))


def test_the_deck_carries_four_steps_and_each_says_how_long_it_takes(parsed):
    assert len(parsed.steps) == 4
    for step in parsed.steps:
        assert step.procedure == "static"
        assert step.time_period and step.initial_increment
        assert step.fixed_increments, (
            f"step {step.index} does not fix its own increment count, so no "
            f"recorded increment can be placed on its loading path")
    assert [s.fixed_increments for s in parsed.steps] == [10, 10, 10, 5]


def test_every_boundary_knows_which_step_wrote_it(parsed):
    assert len(parsed.boundaries) == 96
    assert {b.step for b in parsed.boundaries} == {1, 2, 3, 4}
    assert {b.op for b in parsed.boundaries} == {"NEW"}


def test_a_step_with_op_new_replaces_what_came_before_it(parsed):
    """Abaqus's own rule, applied rather than approximated."""
    for step in (1, 2, 3, 4):
        active = constraints.active_boundaries(parsed, step)
        assert len(active) == 24, (
            f"step {step} should hold the 24 conditions of its own OP=NEW "
            f"block and no others; it holds {len(active)}")
        assert {b.step for b in active} == {step}


def test_reading_the_steps_flat_gives_the_first_step_the_last_steps_motion(parsed):
    """The failure this exists to prevent, measured rather than described."""
    model = from_abaqus(parsed, formulation_policy=lambda e: None)
    dm = DofManager(model.nodes.keys())
    flat = constraints.dirichlet_dofs(model, dm, step=None)
    first = constraints.dirichlet_dofs(model, dm, step=1)
    second = constraints.dirichlet_dofs(model, dm, step=2)
    last = constraints.dirichlet_dofs(model, dm, step=4)

    assert flat == last, (
        "a flat read is the LAST step's conditions wearing the whole deck's "
        "name")
    wrong = sum(1 for dof, value in second.items()
                if abs(flat.get(dof, 0.0) - value) > 1e-15)
    assert wrong >= 4, (
        f"the shear step and the flat read should disagree about several "
        f"degrees of freedom; they disagree about {wrong}")
    assert first != second, "the uniaxial and shear steps prescribe the same "\
                            "motion, so this deck cannot test anything"


def test_conditions_written_before_any_step_are_initial_and_survive_mod(tmp_path):
    """``OP=MOD`` merges onto what is carried forward; ``OP=NEW`` does not."""
    deck = "\n".join([
        "*NODE", "1, 0.0, 0.0, 0.0", "2, 1.0, 0.0, 0.0",
        "*ELEMENT, TYPE=C3D8, ELSET=E", "1, 1, 2, 2, 1, 1, 2, 2, 1",
        "*BOUNDARY", "1, 1, 1, 0.0",
        "*STEP", "*STATIC", "0.1, 1.0",
        "*BOUNDARY", "2, 2, 2, 0.5",
        "*END STEP",
        "*STEP", "*STATIC", "0.1, 1.0",
        "*BOUNDARY, OP=NEW", "2, 3, 3, 0.25",
        "*END STEP", ""])
    path = tmp_path / "ops.inp"
    path.write_text(deck, encoding="utf-8")
    parsed = parse_inp(str(path))

    initial = [b for b in parsed.boundaries if b.step == 0]
    assert len(initial) == 1 and initial[0].op == "MOD"

    first = constraints.active_boundaries(parsed, 1)
    assert len(first) == 2, ("a MOD step keeps the initial condition and adds "
                             "its own")
    second = constraints.active_boundaries(parsed, 2)
    assert len(second) == 1 and second[0].step == 2, (
        "a NEW step replaces everything carried forward, including the "
        "initial conditions")


def test_a_step_whose_increment_count_is_not_fixed_says_so(tmp_path):
    """A step that can cut back does not record how many increments it took,
    and guessing would put a recorded increment at the wrong place on the
    path."""
    deck = "\n".join([
        "*NODE", "1, 0.0, 0.0, 0.0",
        "*STEP", "*STATIC", "0.1, 1.0, 1e-9, 0.5",
        "*BOUNDARY", "1, 1, 1, 0.0", "*END STEP", ""])
    path = tmp_path / "adaptive.inp"
    path.write_text(deck, encoding="utf-8")
    step = parse_inp(str(path)).steps[0]
    assert step.initial_increment == 0.1 and step.max_increment == 0.5
    assert step.fixed_increments is None


def test_nlgeom_is_read_the_way_abaqus_reads_it(tmp_path):
    deck = "\n".join([
        "*NODE", "1, 0.0, 0.0, 0.0",
        "*STEP, NLGEOM=YES", "*STATIC", "0.1, 1.0", "*END STEP",
        "*STEP, NLGEOM=NO", "*STATIC", "0.1, 1.0", "*END STEP",
        "*STEP", "*STATIC", "0.1, 1.0", "*END STEP", ""])
    path = tmp_path / "nlgeom.inp"
    path.write_text(deck, encoding="utf-8")
    steps = parse_inp(str(path)).steps
    assert [s.nlgeom for s in steps] == [True, False, False], (
        "absent means NO, which is Abaqus's default, and reading an absent "
        "flag as YES would put every small-strain deck on the finite-strain "
        "path")
