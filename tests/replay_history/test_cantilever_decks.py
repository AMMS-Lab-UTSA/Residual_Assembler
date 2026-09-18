"""The presentation cantilever decks are the models the slides describe.

``examples/cantilevers/gen_cantilever.py`` rebuilds the J2
(slide 39) and FCC (slide 15) cantilevers from the numbers the slides print.
This pins those numbers -- mesh, DOF, integration points, steps, PROPS, the
clamped root and the prescribed tip -- as the history engine reads them, so a
change to the generator cannot silently turn them into different models.
"""
import sys

import numpy as np
import pytest

from residual_core.replay.history_inputs import read_history_model

from conftest import RA_ROOT

sys.path.insert(0, str(RA_ROOT / "examples" / "cantilevers"))
from gen_cantilever import MODELS, deck  # noqa: E402

SLIDES = {
    # model: elements, nodes, DOF, integration points, push (mm), PROPS
    "j2": (1536, 2499, 7497, 12288, 0.7, [200000.0, 0.3, 250.0, 2000.0]),
    "fcc": (384, 675, 2025, 3072, 0.25,
            [168000.0, 121000.0, 75000.0, 13.0, 55.0, 800.0, 2.0, 1.4, 0.001, 0.05]),
}


@pytest.mark.parametrize("name", sorted(SLIDES))
def test_the_deck_is_the_slide_model(name, tmp_path):
    elements, nodes, dof, points, push, props = SLIDES[name]
    path = tmp_path / f"{name}.inp"
    path.write_text(deck(name, MODELS[name]["props"]))
    model = read_history_model(path)
    assert len(model.element_ids) == elements
    assert len(model.node_ids) == nodes
    assert model.ndof == dof
    assert 8 * len(model.element_ids) == points
    assert list(model.props) == props
    index = {int(node): i for i, node in enumerate(model.node_ids)}
    root = [index[n] for n in model.node_set("ROOT")]
    tip = [index[n] for n in model.node_set("TIP")]
    x = model.coords[:, 0]
    assert np.allclose(x[root], 0.0) and np.allclose(x[tip], x.max())
    assert len(root) == len(tip) == (MODELS[name]["n"][1] + 1) * (MODELS[name]["n"][2] + 1)
    # exactly: every root DOF clamped at zero, and the y DOF of every tip node
    # (global DOF 3*i + 1) taken from 0 to -push over the step
    expected = {3 * i + c: 0.0 for i in root for c in range(3)}
    expected.update({3 * i + 1: -push for i in tip})
    assert sorted(model.prescribed_dofs.tolist()) == sorted(expected)
    assert np.allclose(model.prescribed_start, 0.0)
    for dof_index, value in zip(model.prescribed_dofs.tolist(), model.prescribed_end.tolist()):
        assert value == pytest.approx(expected[dof_index], abs=1e-15)


def test_steps_and_state_variables():
    assert MODELS["j2"]["steps"] == 40 and MODELS["fcc"]["steps"] == 25
    assert MODELS["j2"]["depvar"] == 1 and MODELS["fcc"]["depvar"] == 12
