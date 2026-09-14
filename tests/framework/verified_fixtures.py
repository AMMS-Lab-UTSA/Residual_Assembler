"""The frozen fixtures, and which of them a given check may use.

Two selections, because they answer different questions and conflating them
silently skips work or silently breaks it:

``all_fixtures()``  everything committed, including the plane-stress CPS4 case
                    whose tensor has three components. Loader rules,
                    provenance and convention checks apply to all of them.
``hex_fixtures()``  the eight-node hexahedra with a six-component tensor whose
                    run passed all six of the pipeline's evidence gates -- the
                    ones this repository's C3D8 kernel can assemble a residual
                    from AND whose numbers the pipeline stands behind. A check
                    that integrates ``B^T sigma`` on a unit cube may only use
                    these; handing it a three-component stress is not a smaller
                    test, it is a shape error dressed as one, and handing it a
                    case whose two builds disagreed is checking assembly
                    against a disagreement.

Every fixture is loaded through the fingerprint gate, so a stale one fails
here rather than quietly becoming the baseline.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from residual_core.materials.verified_fixture import (
    CURRENT_TRANSFORM_FINGERPRINT, load_all)

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "verified"

#: A unit cube, which is the geometry the corpus verification decks are
#: generated on, so the element being checked is the element the numbers were
#: produced in.
UNIT_CUBE = np.array([[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.],
                      [0., 0., 1.], [1., 0., 1.], [1., 1., 1.], [0., 1., 1.]])


def all_fixtures(require_gates: bool = False) -> list:
    if not FIXTURES.is_dir() or not any(FIXTURES.glob("*.json")):
        pytest.skip(
            f"no verified fixtures in {FIXTURES}; export them from the "
            f"current store with "
            f"UMAT_source_transformation/tools/export_residual_fixture.py")
    loaded = load_all(FIXTURES, fingerprint=CURRENT_TRANSFORM_FINGERPRINT)
    if require_gates:
        loaded = [f for f in loaded if f.all_six_gates]
    return loaded


def fixture_named(name: str):
    """One committed fixture by file name, for a check about that case.

    Used where a fixture is the SUBJECT rather than one of a set: the
    viscoelastic case whose two builds disagreed is evidence about the
    attribution machinery, not a case to assemble residuals from.
    """
    for fixture in all_fixtures():
        if fixture.path.name == name:
            return fixture
    pytest.skip(f"{name} is not among the committed fixtures")


def hex_fixtures(require_gates: bool = True) -> list:
    usable = [f for f in all_fixtures(require_gates)
              if f.ntens == 6 and f.element_type.upper().startswith("C3D8")
              and not f.element_type.upper().endswith("H")]
    if not usable:
        pytest.skip("no committed fixture is an eight-node hexahedron with a "
                    "six-component tensor")
    return usable


def with_tangent(fixtures=None) -> list:
    source = hex_fixtures() if fixtures is None else fixtures
    usable = [f for f in source
              if f.converted and f.converted[-1].tangent is not None]
    if not usable:
        pytest.skip("no fixture carries a six-component tangent")
    return usable
