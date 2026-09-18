"""Shared set-up for the presentation-claim tests (tests/presentation).

Every test here runs the real claim script end to end on a reduced case: real
builds (gfortran), real compiled providers, real finite-difference references.
Nothing is mocked. Outputs go to pytest's tmp directories, never the repository.
"""

import shutil
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


@pytest.fixture(autouse=True)
def _needs_gfortran():
    if not shutil.which("gfortran"):
        pytest.fail("the presentation claims compile Fortran: gfortran must be on PATH")


@pytest.fixture
def out_dirs(tmp_path, monkeypatch):
    out, work = tmp_path / "out", tmp_path / "work"
    monkeypatch.setenv("PRESENTATION_OUT", str(out))
    monkeypatch.setenv("PRESENTATION_WORK", str(work))
    return out, work
