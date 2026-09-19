"""The scripts documented as ``abaqus python ...`` must parse under Python 2.7.

Abaqus 2021 ships Python 2.7.15. A script with Python 3-only syntax, or with
non-ASCII text and no encoding line, ends in a SyntaxError there before it
reads a single field (``scripts/extract_odb_fields.py`` did, with an em dash in
its docstring, ``from __future__ import annotations`` and ``-> int``).

The offline check walks each script's syntax tree for what Python 2.7 cannot
parse; it needs no Abaqus. The ``abaqus`` test runs the README's command under
the real interpreter and is deselected by the offline suite.
"""

import ast
import re
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]

#: Every script a user is told to run with ``abaqus python``, and what they import.
ABAQUS_PYTHON_SCRIPTS = (
    "scripts/extract_odb_fields.py",
    "scripts/_abaqus_env.py",
    "residual_core/io/abaqus_odb_export.py",
    "residual_core/replay/odb_export_npz.py",
    "examples/cantilevers/export_odb.py",
    "examples/finite_strain_c3d8/extract_odb.py",
    "tests/cp_c3d8_umat/stress_driven_residual/export_fields.py",
)

_ENCODING_LINE = re.compile(rb"^[ \t\f]*#.*?coding[:=][ \t]*([-\w.]+)")


def _python3_only(tree):
    """Constructs of the syntax tree that Python 2.7 cannot parse."""
    found = []
    for node in ast.walk(tree):
        line = getattr(node, "lineno", "?")
        if isinstance(node, ast.ImportFrom) and node.module == "__future__":
            found += ["line %s: from __future__ import %s" % (line, alias.name)
                      for alias in node.names if alias.name == "annotations"]
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            arguments = node.args
            if node.returns is not None or any(
                    argument.annotation is not None for argument in
                    arguments.posonlyargs + arguments.args + arguments.kwonlyargs):
                found.append("line %s: annotation on %s()" % (line, node.name))
            if arguments.kwonlyargs or arguments.posonlyargs:
                found.append("line %s: keyword- or positional-only arguments" % line)
            if isinstance(node, ast.AsyncFunctionDef):
                found.append("line %s: async def" % line)
        elif isinstance(node, (ast.AnnAssign, ast.JoinedStr, ast.NamedExpr, ast.Nonlocal,
                               ast.YieldFrom, ast.Await)):
            found.append("line %s: %s" % (line, type(node).__name__))
        elif isinstance(node, ast.BinOp) and isinstance(node.op, ast.MatMult):
            found.append("line %s: the @ operator" % line)
    return found


@pytest.mark.unit
@pytest.mark.parametrize("relative", ABAQUS_PYTHON_SCRIPTS)
def test_the_script_parses_under_python27(relative):
    source = (ROOT / relative).read_bytes()
    first_two = source.splitlines()[:2]
    declared = any(_ENCODING_LINE.match(line) for line in first_two)
    non_ascii = [number for number, line in enumerate(source.splitlines(), 1)
                 if any(byte > 127 for byte in line)]
    assert declared or not non_ascii, (
        "%s: non-ASCII text on lines %s and no encoding line" % (relative, non_ascii))
    problems = _python3_only(ast.parse(source))
    assert not problems, "%s: %s" % (relative, "; ".join(problems))


@pytest.mark.abaqus
def test_the_readme_export_command_runs_under_abaqus_python():
    if shutil.which("abaqus") is None:
        pytest.fail("the abaqus marker was selected but no abaqus launcher is on PATH")
    result = subprocess.run(["abaqus", "python", str(ROOT / "scripts/extract_odb_fields.py"), "--help"],
                            capture_output=True, text=True, cwd=str(ROOT))
    assert result.returncode == 0, result.stdout + result.stderr
    assert "--odb" in result.stdout
