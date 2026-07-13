"""Mixed / heterogeneous DOF edge cases for the DofManager (audit item 4).

The DOF manager must assign each node the UNION of the DOFs required by the
formulations touching it — never assume all nodes share the same DOFs. Cases:

    (a) a node belonging only to a truss element   -> (UX,UY,UZ)          3 DOF
    (b) a node belonging only to a beam element    -> (UX..RZ)            6 DOF
    (c) a node belonging to BOTH beam and truss    -> union = 6 DOF
    (d) a node belonging to beam AND solid         -> union = 6 DOF
        (solid-only nodes stay at 3 DOF)

Run: python tests/framework/test_dof_manager_mixed.py
"""

import os
import sys

import pytest

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from residual_core.core.model import Model, Element
from residual_core.core.dof_manager import DofManager
from residual_core.formulations.registry import build_formulation_registry

_TRANS = {"UX", "UY", "UZ"}
_ROT = {"RX", "RY", "RZ"}


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


def _build():
    m = Model()
    m.nodes = {i: (float(i), 0.0, 0.0) for i in range(1, 13)}
    #  T3D2  A: [1,3]         node1 truss-only (a); node3 shared with beam (c)
    #  B31   B: [3,2]         node2 beam-only (b)
    #  B31   C: [4,5]         node4 shared with solid (d); node5 beam-only (b)
    #  C3D8  D: [4,6..12]     node4 beam+solid (d); nodes 6..12 solid-only
    m.elements[1] = Element(1, "T3D2", [1, 3])
    m.elements[2] = Element(2, "B31", [3, 2])
    m.elements[3] = Element(3, "B31", [4, 5])
    m.elements[4] = Element(4, "C3D8", [4, 6, 7, 8, 9, 10, 11, 12])
    m.element_formulation = {1: "truss2", 2: "beam2", 3: "beam2",
                             4: "solid_c3d8_small_strain"}
    return m


def test_union_dofs():
    m = _build()
    forms = build_formulation_registry()
    dm = DofManager.for_model(m, forms)
    checks = []

    def dof_set(nid):
        return set(dm.dof_types_of(nid))

    # (a) truss-only
    checks.append(("(a) truss-only node1 = 3 trans",
                   dof_set(1) == _TRANS and dm.n_node_dofs(1) == 3))
    # (b) beam-only
    checks.append(("(b) beam-only node2 = 6 (trans+rot)",
                   dof_set(2) == _TRANS | _ROT and dm.n_node_dofs(2) == 6))
    checks.append(("(b) beam-only node5 = 6",
                   dof_set(5) == _TRANS | _ROT and dm.n_node_dofs(5) == 6))
    # (c) truss + beam -> union 6
    checks.append(("(c) truss+beam node3 = 6 union",
                   dof_set(3) == _TRANS | _ROT and dm.n_node_dofs(3) == 6))
    # (d) beam + solid -> union 6
    checks.append(("(d) beam+solid node4 = 6 union",
                   dof_set(4) == _TRANS | _ROT and dm.n_node_dofs(4) == 6))
    # solid-only stays 3
    checks.append(("solid-only node7 = 3 trans",
                   dof_set(7) == _TRANS and dm.n_node_dofs(7) == 3))

    # global numbering: dense + contiguous, no overlap
    expected_ndof = (1 * 3    # node1 truss
                     + 6      # node2 beam
                     + 6      # node3 union
                     + 6      # node4 union
                     + 6      # node5 beam
                     + 7 * 3)  # nodes 6..12 solid (7 nodes)
    checks.append(("total ndof = %d" % expected_ndof, dm.ndof == expected_ndof))

    # element_dofs gathers exactly each formulation's own DOFs in order
    e_truss = dm.element_dofs(m.elements[1].connectivity, ("UX", "UY", "UZ"))
    checks.append(("truss element_dofs len = 6", len(e_truss) == 6))
    e_beam = dm.element_dofs(m.elements[2].connectivity,
                             ("UX", "UY", "UZ", "RX", "RY", "RZ"))
    checks.append(("beam element_dofs len = 12", len(e_beam) == 12))
    # no duplicate global indices anywhere (disjoint numbering)
    checks.append(("beam element_dofs unique", len(set(e_beam)) == len(e_beam)))

    for label, v in checks:
        print("  %-38s -> %s" % (label, "PASS" if v else "FAIL"))
    failed = [label for label, v in checks if not v]
    assert not failed, "failed checks: " + "; ".join(failed)


def main():
    print("Mixed / heterogeneous DOF edge cases (DofManager)")
    ok = _script_run(test_union_dofs)
    print("OVERALL: %s" % ("ALL PASS" if ok else "FAILURE"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
