"""Generate the minimal end-to-end example assets (audit item 1).

Writes a solver-neutral JSON model (and, where needed, a field export) into each
``examples/minimal_*`` folder so the CLI walkthroughs in their README files are
reproducible with no Abaqus and no low-level Python:

    residual_core/examples/minimal_truss/model.json
    residual_core/examples/minimal_beam/model.json
    residual_core/examples/minimal_mixed/model.json
    residual_core/examples/minimal_c3d8_stress_driven/model.json + fields.json

Run:  python -m residual_core.examples.generate_minimal
"""

from __future__ import annotations

import json
import os

import numpy as np

from ..core.model import Model, Element
from ..io.neutral_model_io import save_model
from ..materials.base import MaterialBinding

HERE = os.path.dirname(__file__)


def _binding(name, section):
    b = MaterialBinding(material=None, name=name)
    b.section = section
    return b


def build_truss() -> Model:
    m = Model(nodes={1: (0.0, 0.0, 0.0), 2: (2.0, 0.0, 0.0), 3: (4.0, 0.0, 0.0)})
    m.elements[1] = Element(1, "T3D2", [1, 2])
    m.elements[2] = Element(2, "T3D2", [2, 3])
    m.element_material = {1: "steel_bar", 2: "steel_bar"}
    m.materials = {"steel_bar": _binding("steel_bar", {"E": 210000.0, "A": 10.0})}
    m.meta = {"title": "two-bar truss (minimal example)"}
    return m


def build_beam() -> Model:
    m = Model(nodes={1: (0.0, 0.0, 0.0), 2: (5.0, 0.0, 0.0)})
    m.elements[1] = Element(1, "B31", [1, 2])
    m.element_material = {1: "beam_section"}
    m.materials = {"beam_section": _binding(
        "beam_section", {"E": 210000.0, "A": 100.0, "Iz": 833.0, "Iy": 833.0,
                         "G": 80000.0, "J": 1666.0})}
    m.meta = {"title": "single cantilever beam (minimal example)"}
    return m


def build_mixed() -> Model:
    # disjoint node sets so each node's DOF set is unambiguous
    m = Model(nodes={1: (0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0),
                     3: (0.0, 1.0, 0.0), 4: (1.0, 1.0, 0.0)})
    m.elements[1] = Element(1, "T3D2", [1, 2])     # bar -> 3 DOF/node
    m.elements[2] = Element(2, "B31", [3, 4])      # beam -> 6 DOF/node
    m.element_material = {1: "bar", 2: "frame"}
    m.materials = {
        "bar": _binding("bar", {"E": 210000.0, "A": 10.0}),
        "frame": _binding("frame", {"E": 210000.0, "A": 100.0, "Iz": 833.0,
                                    "Iy": 833.0, "G": 80000.0, "J": 1666.0}),
    }
    m.meta = {"title": "mixed truss + beam (heterogeneous DOFs)"}
    return m


def build_c3d8() -> Model:
    # single unit-cube C3D8, connectivity in Abaqus node order
    nodes = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
             5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}
    m = Model(nodes={k: tuple(map(float, v)) for k, v in nodes.items()})
    m.elements[1] = Element(1, "C3D8", [1, 2, 3, 4, 5, 6, 7, 8])
    m.element_material = {1: "solid"}
    m.materials = {"solid": _binding("solid", {"E": 210000.0, "nu": 0.3})}
    # reference an external field export (data lives in fields.json, not inlined)
    m.meta = {"title": "single C3D8 stress-driven (minimal example)",
              "field_refs": {"stress_ip": "fields.json"}}
    return m


def _write_fields(path: str) -> None:
    # uniform uniaxial stress S11 = 100 at all 8 integration points
    sigma = np.zeros((8, 6)); sigma[:, 0] = 100.0
    with open(path, "w", encoding="utf-8") as fh:
        json.dump({"stress_ip": {"1": sigma.tolist()}}, fh, indent=2)


def main() -> None:
    specs = [("minimal_truss", build_truss, False),
             ("minimal_beam", build_beam, False),
             ("minimal_mixed", build_mixed, False),
             ("minimal_c3d8_stress_driven", build_c3d8, True)]
    for folder, builder, with_fields in specs:
        d = os.path.join(HERE, folder)
        os.makedirs(d, exist_ok=True)
        save_model(builder(), os.path.join(d, "model.json"))
        if with_fields:
            _write_fields(os.path.join(d, "fields.json"))
        print("wrote %s" % os.path.join(folder, "model.json"))


if __name__ == "__main__":
    main()
