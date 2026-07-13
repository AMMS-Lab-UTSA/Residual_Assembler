"""Runnable examples that demonstrate the framework is model-agnostic.

Each example builds a *small* model directly (no external solver needed),
assembles a residual through the agnostic core, and checks a known result. They
double as living documentation for the high-level API and for adding backends.

Run all:  ``python -m residual_core.ui.examples``
"""

from __future__ import annotations

import numpy as np

from ..core.model import Model, Element
from ..core.dof_manager import DofManager
from ..core.assembler import Assembler
from ..core.state_manager import StateManager
from ..formulations.registry import build_formulation_registry
from ..formulations.truss2 import Truss2
from ..formulations.beam2 import Beam2
from ..materials.base import MaterialBinding


def _assemble(model, forms):
    dm = DofManager.for_model(model, forms)
    asm = Assembler(model, dm, forms, StateManager(2, 0))
    R, _, diag = asm.assemble(np.zeros(dm.ndof))
    return dm, asm, R, diag


def example_truss_axial():
    """A single bar stretched by a unit axial displacement -> F = EA/L."""
    E, A, L = 210e3, 10.0, 2.0
    model = Model(nodes={1: (0, 0, 0), 2: (L, 0, 0)})
    model.elements[1] = Element(1, "T3D2", [1, 2])
    model.element_formulation[1] = "truss2"
    model.element_material[1] = "bar"
    b = MaterialBinding(material=None, name="bar")
    b.section = {"E": E, "A": A}
    model.materials["bar"] = b

    forms = build_formulation_registry()
    dm = DofManager.for_model(model, forms)
    asm = Assembler(model, dm, forms, StateManager(1, 0))
    U = np.zeros(dm.ndof)
    U[dm.element_dofs([2], Truss2.dof_types)[0]] = 1.0   # UX of node 2 = 1
    R, _, _ = asm.assemble(U)
    fx2 = R[dm.element_dofs([2], Truss2.dof_types)[0]]
    print("[truss] internal axial force = %.4f  (expected EA/L = %.4f)"
          % (fx2, E * A / L))
    assert abs(fx2 - E * A / L) < 1e-6
    return R


def example_beam_cantilever():
    """Cantilever with a tip transverse load -> tip deflection PL^3/(3EI)."""
    E, A, L = 210e3, 100.0, 10.0
    I = 1.0 / 12 * 10 * 10 ** 3   # rectangular-ish
    model = Model(nodes={1: (0, 0, 0), 2: (L, 0, 0)})
    model.elements[1] = Element(1, "B31", [1, 2])
    model.element_formulation[1] = "beam2"
    model.element_material[1] = "sec"
    b = MaterialBinding(material=None, name="sec")
    b.section = {"E": E, "A": A, "Iz": I, "Iy": I, "G": E / 2.6, "J": 2 * I}
    model.materials["sec"] = b

    forms = build_formulation_registry()
    dm = DofManager.for_model(model, forms)
    asm = Assembler(model, dm, forms, StateManager(1, 0))
    beam = forms.get("beam2")
    Kg, _L = beam.global_stiffness(model.coords_of([1, 2]), b.section)
    # node1 fully fixed: solve for node2 dofs under transverse tip load P in UY
    all_local = beam.dof_types
    n2 = dm.element_dofs([2], all_local)          # 6 dofs of node2
    P = 1000.0
    Kff = Kg[6:, 6:]
    f = np.zeros(6); f[1] = P                       # UY at node 2
    d = np.linalg.solve(Kff, f)
    tip = d[1]
    print("[beam]  tip deflection = %.6e  (expected PL^3/3EI = %.6e)"
          % (tip, P * L ** 3 / (3 * E * I)))
    assert abs(tip - P * L ** 3 / (3 * E * I)) / (P * L ** 3 / (3 * E * I)) < 1e-9
    return d


def example_mixed_dispatch():
    """Truss + beam elements in one model on disjoint nodes: the physics-blind
    core dispatches each to its own backend and builds heterogeneous DOFs."""
    model = Model(nodes={1: (0, 0, 0), 2: (1, 0, 0),
                         3: (0, 1, 0), 4: (1, 1, 0)})
    model.elements[1] = Element(1, "T3D2", [1, 2])
    model.elements[2] = Element(2, "B31", [3, 4])
    model.element_formulation = {1: "truss2", 2: "beam2"}
    model.element_material = {1: "bar", 2: "sec"}
    bar = MaterialBinding(material=None, name="bar")
    bar.section = {"E": 210e3, "A": 10.0}
    sec = MaterialBinding(material=None, name="sec")
    sec.section = {"E": 210e3, "A": 100.0, "Iz": 833.0, "Iy": 833.0,
                   "G": 80e3, "J": 1666.0}
    model.materials = {"bar": bar, "sec": sec}
    forms = build_formulation_registry()
    dm = DofManager.for_model(model, forms)
    print("[mixed] node1 dofs=%s  node3 dofs=%s  ndof=%d"
          % (dm.dof_types_of(1), dm.dof_types_of(3), dm.ndof))
    assert dm.n_node_dofs(1) == 3 and dm.n_node_dofs(3) == 6
    asm = Assembler(model, dm, forms, StateManager(1, 0))
    R, _, diag = asm.assemble(np.zeros(dm.ndof))
    assert diag["elements"] == 2
    return R


EXAMPLES = {
    "truss_axial": example_truss_axial,
    "beam_cantilever": example_beam_cantilever,
    "mixed_dispatch": example_mixed_dispatch,
}


def run_all():
    for name, fn in EXAMPLES.items():
        print("== %s ==" % name)
        fn()
    print("all examples passed.")


if __name__ == "__main__":
    run_all()
