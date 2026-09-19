"""The ODB exporter tells its user to run a file that exists.

The exporter moved to ``residual_core/io/abaqus_odb_export.py``, but its usage
text and the error it prints under an ordinary Python still said ``abaqus
python extract_abaqus_fields.py``, a file the repository no longer has
(4d34c74). Every command they give is checked here against the repository:
the path it names, taken from the repository root as the usage text is
written, must be the exporter itself.
"""
from __future__ import annotations

import ast
import re
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.regression

ROOT = Path(__file__).resolve().parents[2]
EXPORTER = ROOT / "residual_core" / "io" / "abaqus_odb_export.py"
#: ``abaqus python <script>`` and ``abaqus cae noGUI=<script>``
COMMAND = re.compile(r"abaqus\s+(?:python\s+|cae\s+noGUI=)(\S+\.py)\b")


def _named(text: str) -> list:
    named = COMMAND.findall(text)
    assert named, f"no abaqus command found in:\n{text}"
    return named


def test_the_error_under_an_ordinary_python_names_the_exporter():
    # odbAccess exists only in Abaqus' own Python, so this interpreter takes the error path
    completed = subprocess.run([sys.executable, str(EXPORTER), "--", "--odb", "job.odb"],
                               cwd=ROOT, capture_output=True, text=True)
    assert completed.returncode == 2, completed.stdout + completed.stderr
    assert "must be run with Abaqus' Python" in completed.stderr
    for path in _named(completed.stderr):
        assert (ROOT / path).resolve() == EXPORTER, path


def test_every_command_in_the_usage_text_runs_the_exporter():
    usage = ast.get_docstring(ast.parse(EXPORTER.read_text(encoding="utf-8")))
    named = _named(usage)
    assert len(named) >= 5, named          # the usage text gives five commands
    for path in named:
        assert (ROOT / path).resolve() == EXPORTER, path
