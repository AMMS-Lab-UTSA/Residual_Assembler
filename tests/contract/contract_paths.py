"""Where the contract's files are, for the contract tests.

Separate from ``conftest.py`` and named distinctly on purpose. ``conftest`` is
a name pytest gives to a file in every test directory, so importing helpers
"from conftest" resolves to whichever one happens to be first on sys.path --
and this repository is under active development in its other test
directories. A module nothing else will be called cannot be shadowed.

These tests are about the AGREEMENT between this repository and
UMAT_source_transformation, not about the assembler. They import nothing from
``residual_core`` on purpose: the contract has to be checkable in a checkout
where the assembler does not import, because a contract that can only be
verified when everything else already works is not a contract, it is a
symptom.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SCHEMAS = REPO / "schemas"
FIXTURES = REPO / "tests" / "fixtures" / "verified"
SHARED = Path(__file__).resolve().parent / "_shared"

if str(SHARED) not in sys.path:
    sys.path.insert(0, str(SHARED))

#: Directory names the producing repository is checked out under, searched
#: beside this one. Relative, never absolute: an absolute path is true on
#: exactly one computer, and a cross-repository check that silently skips
#: everywhere else is a check that holds nowhere else.
#:
#: The comparison it enables is opportunistic. ``schemas/contract_lock.json``
#: is what holds when the other checkout is absent, and the skip says so
#: rather than passing quietly.
UMAT_REPO_NAMES = ("UMAT_source_transformation", "wt-D-primal")
UMAT_CONTRACT_RELATIVE = Path("src") / "umat_oti" / "contract"


def umat_repo():
    """The producing checkout beside this one, or ``None``."""
    override = os.environ.get("UMAT_OTI_REPO")
    if override:
        candidate = Path(override)
        return candidate if (candidate / UMAT_CONTRACT_RELATIVE).is_dir() else None
    for base in (REPO.parent, REPO.parent.parent):
        for name in UMAT_REPO_NAMES:
            candidate = base / name
            if (candidate / UMAT_CONTRACT_RELATIVE / "schemas").is_dir():
                return candidate
    return None


def load_schema(name: str) -> dict:
    """One shared schema document by stem, e.g. ``umat_contract_v1``."""
    return json.loads((SCHEMAS / f"{name}.schema.json").read_text())
