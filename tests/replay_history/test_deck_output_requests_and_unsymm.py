"""A deck as Abaqus ran it: UNSYMM and results-file output requests are accepted, not applied.

``*User Material, UNSYMM`` only tells Abaqus to store its own Newton matrix
unsymmetric, and ``*File Format`` / ``*Node File`` / ``*El File`` (and the
printed-output requests) only choose what Abaqus writes. Neither enters
R = F_int - F_ext, so the replay reads the deck the analysis ran -- with them --
into exactly the model it reads without them, and names what it accepted.
"""
import numpy as np
import pytest

from residual_core.replay.history_inputs import UnsupportedFeature

from history_support import beam_model

OUTPUT_LINES = ("*File Format, ascii\n*Node File, frequency=1\nU, RF\n*El File, frequency=1, position=centroidal\n"
                "S, SDV\n*Node Print, nset=TIP, frequency=5\nU\n*El Print, frequency=0\nS\n*End Step")


def _as_abaqus_ran_it(text):
    return text.replace("*User Material, constants=4", "*User Material, constants=4, unsymm").replace(
        "*End Step", OUTPUT_LINES)


def test_the_deck_abaqus_ran_reads_into_the_same_model(tmp_path):
    plain_path, plain = beam_model(tmp_path, name="plain.inp")
    ran_path, ran = beam_model(tmp_path, name="ran.inp", extra=_as_abaqus_ran_it)
    assert "unsymm" in ran_path.read_text() and "*Node File" in ran_path.read_text()
    for name in ("node_ids", "coords", "element_ids", "connectivity", "props", "prescribed_dofs",
                 "prescribed_start", "prescribed_end", "load_end"):
        np.testing.assert_array_equal(getattr(plain, name), getattr(ran, name), err_msg=name)
    assert (plain.node_sets, plain.element_sets, plain.depvar, plain.time_period) == (
        ran.node_sets, ran.element_sets, ran.depvar, ran.time_period)
    assert plain.notes == []
    assert any(note.startswith("*User Material, UNSYMM") for note in ran.notes)
    for keyword in ("*File Format", "*Node File", "*El File", "*Node Print", "*El Print"):
        assert any(note.startswith(keyword + " ") for note in ran.notes), keyword


def test_an_option_that_changes_the_material_is_still_refused(tmp_path):
    with pytest.raises(UnsupportedFeature, match="mechanical"):
        beam_model(tmp_path, extra=lambda t: t.replace("*User Material, constants=4",
                                                       "*User Material, constants=4, type=THERMAL"))
    with pytest.raises(UnsupportedFeature, match="hybrid"):
        beam_model(tmp_path, name="h.inp", extra=lambda t: t.replace(
            "*User Material, constants=4", "*User Material, constants=4, hybrid formulation=total"))
