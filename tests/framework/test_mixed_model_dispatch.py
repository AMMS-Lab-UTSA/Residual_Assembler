"""Mixed-backend dispatch — the core assembler is physics-blind.

Builds a toy model with THREE different formulation families sharing one global
system:

    element 1: truss  (T3D2) -> truss2   backend, 3 DOF/node
    element 2: beam   (B31)  -> beam2    backend, 6 DOF/node (rotations!)
    element 3: solid  (C3D8) -> solid_c3d8_small_strain + elastic, 3 DOF/node

The nodes are disjoint per element, so each node's DOF set is unambiguous and the
heterogeneous DofManager assigns 3, 6, or 3 DOFs per node accordingly. The single
generic Assembler dispatches every element to its bound backend and scatters the
result WITHOUT knowing what any of them are.

Proof of correct dispatch + scatter: the global residual restricted to each
element's own DOFs equals that element evaluated on its own. If the assembler
mixed up DOF maps or physics, this would fail.

Run: python tests/framework/test_mixed_model_dispatch.py
"""

import os
import sys

import numpy as np
import pytest

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from residual_core.core.model import Model, Element
from residual_core.core.dof_manager import DofManager
from residual_core.core.assembler import Assembler
from residual_core.core.state_manager import StateManager
from residual_core.formulations.truss2 import Truss2
from residual_core.formulations.beam2 import Beam2
from residual_core.formulations.solid_c3d8_small_strain import SolidC3D8SmallStrain
from residual_core.materials.elastic_adapter import IsotropicElastic
from residual_core.materials.base import MaterialBinding


def _script_run(fn, *args):
    """Script-mode runner: a test passes unless it raises AssertionError."""
    try:
        fn(*args)
    except AssertionError as exc:
        print("  FAIL: %s" % exc)
        return False
    except pytest.skip.Exception as exc:
        print("  SKIP: %s" % exc)
    return True


def _build_model():
    m = Model()
    m.nodes = {
        1: (0, 0, 0), 2: (1, 0, 0),                  # truss
        3: (0, 2, 0), 4: (1, 2, 0),                  # beam
        5: (0, 0, 5), 6: (1, 0, 5), 7: (1, 1, 5), 8: (0, 1, 5),   # solid bottom
        9: (0, 0, 6), 10: (1, 0, 6), 11: (1, 1, 6), 12: (0, 1, 6),  # solid top
    }
    m.elements = {
        1: Element(1, "T3D2", [1, 2]),
        2: Element(2, "B31", [3, 4]),
        3: Element(3, "C3D8", [5, 6, 7, 8, 9, 10, 11, 12]),
    }
    m.element_formulation = {1: "truss2", 2: "beam2", 3: "solid_c3d8_small_strain"}
    m.element_material = {1: "bar", 2: "beamsec", 3: "steel"}

    bar = MaterialBinding(material=None, name="bar"); bar.section = {"E": 200e3, "A": 8.0}
    beamsec = MaterialBinding(material=None, name="beamsec")
    beamsec.section = {"E": 210e3, "A": 50.0, "Iz": 400.0, "Iy": 400.0, "J": 700.0}
    steel = MaterialBinding(IsotropicElastic(), [200e3, 0.3], 0, "steel")
    m.materials = {"bar": bar, "beamsec": beamsec, "steel": steel}
    return m


def test_dispatch_and_scatter():
    m = _build_model()
    forms = {"truss2": Truss2(), "beam2": Beam2(),
             "solid_c3d8_small_strain": SolidC3D8SmallStrain()}
    dm = DofManager.for_model(m, forms)
    asm = Assembler(m, dm, forms, StateManager(8, 0))

    expected_ndof = 2 * 3 + 2 * 6 + 8 * 3      # truss + beam + solid
    rng = np.random.default_rng(7)
    U = rng.uniform(-1e-3, 1e-3, dm.ndof)

    R, K, diag = asm.assemble(U, compute_tangent=True)

    checks = []
    checks.append(("ndof", dm.ndof == expected_ndof))
    checks.append(("elements", diag["elements"] == 3))
    checks.append(("all three backends fired",
                   diag["formulations"] == {"truss2": 1, "beam2": 1,
                                            "solid_c3d8_small_strain": 1}))

    # node DOF-set correctness (heterogeneous)
    checks.append(("truss node 3 DOFs", dm.n_node_dofs(1) == 3))
    checks.append(("beam node has rotations", dm.n_node_dofs(3) == 6))
    checks.append(("solid node 3 DOFs", dm.n_node_dofs(5) == 3))

    # scatter correctness: R restricted to each element's DOFs == standalone eval
    for eid, el in m.elements.items():
        form = forms[m.element_formulation[eid]]
        edofs = np.array(dm.element_dofs(el.connectivity, form.dof_types), int)
        coords = m.coords_of(el.connectivity)
        props = m.materials[m.element_material[eid]]
        state = np.zeros((8, 0)) if el.etype == "C3D8" else None
        r_e, _, _, _ = form.eval_element(eid, el.etype, coords, U[edofs], {},
                                         state, props, (0, 0), 0.0, None,
                                         {"compute_tangent": False})
        # only this element touches these DOFs (disjoint nodes) -> exact match
        rel = np.linalg.norm(R[edofs] - r_e) / max(np.linalg.norm(r_e), 1e-30)
        checks.append(("scatter %s (elem %d)" % (form.name, eid), rel < 1e-10))

    for label, v in checks:
        print("  %-32s -> %s" % (label, "PASS" if v else "FAIL"))
    failed = [label for label, v in checks if not v]
    assert not failed, "failed checks: " + "; ".join(failed)


def main():
    print("Mixed-backend dispatch verification (physics-blind core)")
    ok = _script_run(test_dispatch_and_scatter)
    print("OVERALL: %s" % ("ALL PASS" if ok else "FAILURE"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
