"""The cantilever central-difference reference records its method as its note.

``examples/cantilevers/fd_reference.py`` writes ``fd_reference.json`` with a
``note`` that tells a reader how the numbers were made. It recorded the
script's usage line there instead of its method paragraph (718ac7c). The
script is run here, as the example's README runs it, on reruns whose exported
fields are exactly proportional to the perturbed parameter, so the central
differences it records are also known exactly.
"""
from __future__ import annotations

import json
import subprocess
import sys

import numpy as np
import pytest

from history_support import RA_ROOT

pytestmark = pytest.mark.regression

SCRIPT = RA_ROOT / "examples" / "cantilevers" / "fd_reference.py"
sys.path.insert(0, str(SCRIPT.parent))
from gen_cantilever import MODELS  # noqa: E402

MODEL = "j2"
STEPS = ("0.02", "0.01", "0.005")


def _fields(scale: float) -> dict:
    """A 2 x 1 x 1 C3D8 bar exported as export_odb.py writes it, times ``scale``."""
    grid = np.array([(x, y, z) for z in (0, 1) for y in (0, 1) for x in (0, 1, 2)], dtype=float)
    label = {tuple(point): index + 1 for index, point in enumerate(grid)}

    def element(x0):
        corners = [(x0, 0, 0), (x0 + 1, 0, 0), (x0 + 1, 1, 0), (x0, 1, 0),
                   (x0, 0, 1), (x0 + 1, 0, 1), (x0 + 1, 1, 1), (x0, 1, 1)]
        return [label[corner] for corner in corners]

    rng = np.random.default_rng(718)
    frames, nodes, elements = 3, len(grid), 2
    return {"coords": grid, "conn": np.array([element(0), element(1)]),
            "U": scale * rng.normal(size=(frames, nodes, 3)),
            "RF": scale * rng.normal(size=(frames, nodes, 3)),
            "S": scale * rng.normal(size=(frames, elements, 8, 6))}


def test_the_note_is_the_method_and_the_differences_are_exact(tmp_path):
    folder = tmp_path / MODEL
    folder.mkdir()
    spec = MODELS[MODEL]
    parameter, p0 = spec["prop_names"][0], spec["props"][0]
    np.savez(folder / f"cantilever_{MODEL}_nominal_fields.npz", **_fields(1.0))
    manifest = ["job\tparameter\tsign\trel_step\tvalue\tstatus\tseconds"]
    for step in STEPS:
        for sign, factor in (("p", 1 + float(step)), ("m", 1 - float(step))):
            job = f"cantilever_{MODEL}_{parameter}_{sign}{step.replace('.', 'p')}"
            # q(p) = q0 p / p0: every output of these reruns is linear in the parameter
            np.savez(folder / f"{job}_fields.npz", **_fields(factor))
            manifest.append(f"{job}\t{parameter}\t{sign}\t{step}\t{p0 * factor!r}\tcompleted\t1")
    (folder / "fd_manifest.tsv").write_text("\n".join(manifest) + "\n")

    completed = subprocess.run([sys.executable, str(SCRIPT), MODEL, str(tmp_path)],
                               capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    result = json.loads((folder / "fd_reference.json").read_text())

    note = result["note"]
    assert note.startswith("For every parameter p and relative step h"), note
    assert "central difference" in note and "(2 h p)" in note
    assert "python fd_reference.py" not in note and "->" not in note

    outputs = result["outputs"]
    assert set(outputs) == {"tip_RF2", "U2_mid_top", "S11_root_ip", "vm_root_max"}
    nominal = _fields(1.0)
    tip = nominal["coords"][:, 0] == 2
    mid_top = np.flatnonzero((nominal["coords"] == (1, 1, 0)).all(axis=1))[0]
    expected = {"tip_RF2": nominal["RF"][:, tip, 1].sum(axis=1) / p0,
                "U2_mid_top": nominal["U"][:, mid_top, 1] / p0,
                "S11_root_ip": nominal["S"][:, 0, 0, 0] / p0}
    for key, value in expected.items():
        for step in STEPS:
            np.testing.assert_allclose(outputs[key][parameter]["by_step_size"][step], value,
                                       rtol=1e-9, atol=1e-12 * np.abs(value).max(), err_msg=f"{key} h={step}")
