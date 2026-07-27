"""pytest configuration for the residual-assembler test suite.

Puts the repo root on sys.path (the standalone tests already do this themselves,
but this lets pytest import them cleanly too) and honours the opt-in markers:
tests marked `otilib` / `abaqus` are skipped unless RESASM_RUN_OTILIB=1 /
RESASM_RUN_ABAQUS=1 are set, so a plain `pytest` run stays fully offline.
"""
import os
import sys

import pytest

RA = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if RA not in sys.path:
    sys.path.insert(0, RA)

_GATES = {
    "otilib": ("RESASM_RUN_OTILIB", "needs OTILib (set RESASM_RUN_OTILIB=1)"),
    "abaqus": ("RESASM_RUN_ABAQUS", "needs a real Abaqus run (set RESASM_RUN_ABAQUS=1)"),
}


def pytest_collection_modifyitems(config, items):
    for item in items:
        for marker, (env, reason) in _GATES.items():
            if item.get_closest_marker(marker) and os.environ.get(env) != "1":
                item.add_marker(pytest.mark.skip(reason=reason))
