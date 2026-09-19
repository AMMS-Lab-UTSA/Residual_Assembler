"""The general deck path refuses a deck state it cannot reproduce.

``resasm requirements | doctor | assemble | verify MODEL.inp`` assemble one
state from one set of loads and supports: every ``*Cload`` of the deck at once
(``core.loads.external_force``) and every ``*Boundary`` as one list
(``core.constraints.dirichlet_dofs(step=None)``), with stress-driven and
small-strain elements integrated over the reference configuration. A deck with
several steps, or a geometrically nonlinear step assembled that way, would
give a residual of a load case and a geometry the analysis never had. The
readiness says so ("one-step deck in scope: no", with the reason), and
assembly and verification refuse with exit 2. ``resasm history`` replays a
deck step by step.

The decks are the shipped one-element example, with a second step or with
NLGEOM=YES; its stress field (uniaxial 300, the exact equilibrium of the
75-per-node load) makes the one-step deck verify.
"""
import json
from pathlib import Path

import pytest

from residual_core.ui import cli

ROOT = Path(__file__).resolve().parents[2]
CLEAN = ROOT / "examples" / "presentation_request" / "Analysis.inp"

SECOND_STEP = ("*Step, name=Unloading, nlgeom=NO\n*Static\n0.5, 1.\n*Cload, op=NEW\n"
               "LOADED, 1, 30.\n*End Step\n")


def _deck(tmp_path, change):
    deck = tmp_path / "Analysis.inp"
    deck.write_text(change(CLEAN.read_text()))
    fields = tmp_path / "fields.json"
    fields.write_text(json.dumps({"stress_ip": {"1": [[300.0, 0, 0, 0, 0, 0]] * 8}}))
    return deck, fields


def run(capsys, *argv):
    code = cli.main([str(a) for a in argv])
    captured = capsys.readouterr()
    return code, captured.out, captured.err


def test_the_one_step_deck_verifies(tmp_path, capsys):
    deck, fields = _deck(tmp_path, lambda text: text)
    code, out, err = run(capsys, "verify", deck, "--fields", fields)
    assert code == 0, out + err
    assert "equilibrium          = PASS" in out
    code, out, _ = run(capsys, "requirements", deck, "--mode", "stress-driven", "--fields", fields)
    assert "  one-step deck in scope: yes" in out


@pytest.mark.parametrize("change, reason", [
    (lambda text: text.rstrip("\n") + "\n" + SECOND_STEP,
     "the deck has 2 steps, and the general assembly would apply the *Cload lines of all "
     "of them at once"),
    (lambda text: text.replace("nlgeom=NO", "nlgeom=YES"),
     "the deck's step is geometrically nonlinear (NLGEOM=YES), but stress-driven assembly "
     "integrates over the reference configuration"),
], ids=["two-steps", "NLGEOM"])
def test_a_deck_state_the_general_path_cannot_place_is_refused(tmp_path, capsys, change, reason):
    """Regression: a two-step deck was assembled with both steps' loads summed
    (75 + 30 per node, although the second step's OP=NEW leaves 30), and an
    NLGEOM=YES deck in the reference configuration; verify printed an
    equilibrium verdict for both."""
    deck, fields = _deck(tmp_path, change)
    code, out, _ = run(capsys, "requirements", deck, "--mode", "stress-driven", "--fields", fields)
    assert code == 0
    assert "  one-step deck in scope: no" in out
    assert "  Why: %s" % reason in out
    code, out, err = run(capsys, "assemble", deck, "--mode", "stress-driven", "--fields", fields)
    assert code == 2 and out == ""
    assert "  Why: %s" % reason in err
    code, out, err = run(capsys, "verify", deck, "--fields", fields)
    assert code == 2
    assert "equilibrium          =" not in out
    assert "  Why: %s" % reason in err
