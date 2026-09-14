"""Fixtures for the contract tests.

These tests are about the AGREEMENT between this repository and
UMAT_source_transformation, not about the assembler. They import nothing from
``residual_core`` on purpose: the contract has to be checkable in a checkout
where the assembler does not import, because a contract that can only be
verified when everything else already works is not a contract, it is a
symptom.

The paths themselves live in :mod:`contract_paths`; see there for why they are
not in this file.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from contract_paths import FIXTURES, SCHEMAS, load_schema, umat_repo  # noqa: F401


@pytest.fixture(scope="session")
def schemas_dir() -> Path:
    return SCHEMAS


@pytest.fixture(scope="session")
def fixtures_dir() -> Path:
    return FIXTURES


@pytest.fixture(scope="session")
def lock() -> dict:
    return json.loads((SCHEMAS / "contract_lock.json").read_text())
