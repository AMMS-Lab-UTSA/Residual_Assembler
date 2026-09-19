"""The general deck path names every equilibrium keyword it does not apply.

``resasm inspect | assemble | verify MODEL.inp`` read the deck through
``parse_inp`` and ``core.model.from_abaqus``. That path assembles concentrated
loads and Dirichlet conditions only: ``*Dsload`` and ``*Equation`` are parsed
but not applied, ``*Dload`` and ``*Amplitude`` are not read, and an
``AMPLITUDE=`` on a load or boundary line is ignored. Each of those changes
equilibrium, so none may be dropped silently: ``from_abaqus`` records one note
per keyword on ``Model.unapplied_keywords`` and issues a
``DeckKeywordNotApplied`` warning, which reaches the terminal of every command
on this path, and ``resasm inspect`` lists the notes in its summary.

The decks are the shipped one-element example with those keywords added.
"""
import json
import subprocess
import sys
import warnings
from pathlib import Path

import pytest

from residual_core.core.model import (DeckKeywordNotApplied, from_abaqus,
                                      unapplied_deck_keywords)
from residual_core.io.abaqus_inp_parser import parse_inp
from residual_core.ui import cli

ROOT = Path(__file__).resolve().parents[2]
CLEAN = ROOT / "examples" / "presentation_request" / "Analysis.inp"

#: the start of the note for each keyword the deck below adds
EXPECTED = (
    "*Dsload (1 line(s)): parsed but not applied",
    "*Equation (1 constraint(s)): parsed but not applied",
    "AMPLITUDE=RAMP (1 *Cload line(s)): not applied",
    "*Amplitude: present but not applied",
    "*Dload: present but not applied",
)


def _deck_with_unapplied_keywords(directory, equation=True):
    deck = CLEAN.read_text()
    model_data = ("*Amplitude, name=RAMP\n0., 0., 1., 1.\n"
                  "*Surface, type=ELEMENT, name=TOP\nALL, S2\n")
    if equation:
        model_data += "*Equation\n2\n2, 2, 1., 3, 2, -1.\n"
    deck = deck.replace("*Step, name=Loading", model_data + "*Step, name=Loading", 1)
    deck = deck.replace("*Cload\n", "*Cload, amplitude=RAMP\n", 1)
    deck = deck.replace("*Output, field", "*Dsload\nTOP, P, 10.\n"
                        "*Dload\nALL, GRAV, 9.81, 0., 0., -1.\n*Output, field", 1)
    path = Path(directory) / "Analysis.inp"
    path.write_text(deck)
    fields = Path(directory) / "fields.json"
    fields.write_text(json.dumps({"stress_ip": {"1": [[0.0] * 6] * 8}}))
    return path, fields


def _cli(*argv):
    return subprocess.run([sys.executable, "-m", "residual_core.ui.cli", *map(str, argv)],
                          cwd=str(ROOT), capture_output=True, text=True, check=False)


def test_each_unapplied_keyword_is_named_where_the_deck_is_read(tmp_path):
    deck, _ = _deck_with_unapplied_keywords(tmp_path)
    parsed = parse_inp(str(deck))
    with pytest.warns(DeckKeywordNotApplied) as caught:
        model = from_abaqus(parsed)
    assert len(model.unapplied_keywords) == len(EXPECTED)
    for note, expected in zip(model.unapplied_keywords, EXPECTED):
        assert note.startswith(expected), (note, expected)
    message = str(caught[0].message)
    for expected in EXPECTED:
        assert expected in message
    # what the notes say is what the parse holds: nothing is invented
    assert len(parsed.dsloads) == 1 and len(parsed.equations) == 1
    assert [c.amplitude for c in parsed.cloads] == ["RAMP"]


def test_a_deck_whose_keywords_are_all_applied_gets_no_note():
    """Output requests, the heading and *End Step change nothing and are not
    reported; the shipped deck has them and nothing else unapplied."""
    parsed = parse_inp(str(CLEAN))
    assert {k.lower() for k in parsed.unsupported_keywords} >= {"heading", "output"}
    assert unapplied_deck_keywords(parsed) == []
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeckKeywordNotApplied)
        assert from_abaqus(parsed).unapplied_keywords == []


def test_any_other_unread_keyword_is_named_too(tmp_path):
    deck = CLEAN.read_text().replace(
        "*Solid Section", "*Initial Conditions, type=STRESS\nALL, 1., 0., 0., 0., 0., 0.\n"
        "*Solid Section", 1)
    path = tmp_path / "initial.inp"
    path.write_text(deck)
    assert unapplied_deck_keywords(parse_inp(str(path))) == [
        "*Initial Conditions: present but not read, so not applied"]


def test_resasm_inspect_lists_them_in_its_summary(tmp_path, capsys):
    deck, _ = _deck_with_unapplied_keywords(tmp_path)
    with pytest.warns(DeckKeywordNotApplied):
        assert cli.main(["inspect", str(deck)]) == 0
    out = capsys.readouterr().out
    section = out.split("Deck keywords present but NOT applied (the residual omits them):\n",
                        1)[1].split("\n\n", 1)[0]
    assert len(section.splitlines()) == len(EXPECTED)
    for expected in EXPECTED:
        assert "  - " + expected in section


def test_resasm_inspect_of_the_clean_deck_has_no_such_section(capsys):
    assert cli.main(["inspect", str(CLEAN)]) == 0
    assert "NOT applied" not in capsys.readouterr().out


@pytest.mark.parametrize("command", ["inspect", "assemble", "verify"])
def test_every_general_deck_command_prints_the_diagnostic(tmp_path, command):
    deck, fields = _deck_with_unapplied_keywords(tmp_path, equation=command != "verify")
    extra = {"inspect": (), "assemble": ("--mode", "stress-driven", "--fields", fields),
             "verify": ("--fields", fields)}[command]
    result = _cli(command, deck, *extra)
    assert "DeckKeywordNotApplied" in result.stderr, result.stdout + result.stderr
    for expected in EXPECTED:
        if expected.startswith("*Equation") and command == "verify":
            continue
        assert expected in result.stderr
    if command == "assemble":
        # the residual is still assembled, from the loads that ARE applied
        assert result.returncode == 0
        assert "assembled mode 'stress-driven': ndof=24" in result.stdout


def test_verify_refuses_a_deck_whose_equilibrium_keywords_are_not_applied(tmp_path):
    """Regression: verify printed the diagnostic and then an equilibrium verdict
    computed without the dropped loads. It now refuses with exit 2, as it
    already did for *Equation, naming every keyword it would have dropped."""
    deck, fields = _deck_with_unapplied_keywords(tmp_path, equation=False)
    result = _cli("verify", deck, "--fields", fields)
    assert result.returncode == 2, result.stdout + result.stderr
    assert "equilibrium          =" not in result.stdout
    refusal = [line for line in result.stderr.splitlines()
               if line.startswith("Cannot verify equilibrium: ")]
    assert len(refusal) == 1
    for expected in EXPECTED:
        if not expected.startswith("*Equation"):
            assert expected in refusal[0]


def test_verify_still_runs_when_only_a_material_keyword_is_unread(tmp_path):
    """The stress-driven residual takes the stress from the field, so an unread
    material keyword (*Elastic here) cannot change it: verify gives its verdict."""
    deck = tmp_path / "Analysis.inp"
    deck.write_text(CLEAN.read_text().replace("*Solid Section", "*Elastic\n210000., 0.3\n"
                                              "*Solid Section", 1))
    fields = tmp_path / "fields.json"
    fields.write_text(json.dumps({"stress_ip": {"1": [[300.0, 0.0, 0.0, 0.0, 0.0, 0.0]] * 8}}))
    result = _cli("verify", deck, "--fields", fields)
    assert "*Elastic: present but not read, so not applied" in result.stderr
    assert result.returncode == 0, result.stdout + result.stderr
    assert "equilibrium          = PASS" in result.stdout
