"""CI accepts the skip the live corpus tests give when CORPUS_WORKSPACE is unset.

The CI step "Every skip must name a declared external dependency" fails the run
when a test skips for a reason that is not in its ``declared`` tuple. The live
corpus tests need the private corpus workspace (never redistributed), which a
public runner cannot have, so the reason they give when ``CORPUS_WORKSPACE`` is
unset must be declared there. It was not, from 6066af8 until f85ad0a, and the
step failed on 15 such skips.

The reason is not copied into this file. It is taken from the live tests
themselves: ``_live()`` of test_corpus_residual.py (which the solver-robustness
tests call too) is run with ``CORPUS_WORKSPACE`` removed from the environment,
and the message of the skip it raises is the reason. The step's own Python
block is then read out of .github/workflows/ci.yml and run, unchanged except
for the path of the suite output it reads, on a suite output that holds that
skip. A control output with an undeclared reason must still fail the step.
"""
from __future__ import annotations

import ast
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW = REPO_ROOT / ".github" / "workflows" / "ci.yml"
STEP_NAME = "Every skip must name a declared external dependency"


def _step_python_block() -> str:
    """The Python heredoc of the named CI step, dedented, exactly as CI runs it."""
    lines = WORKFLOW.read_text(encoding="utf-8").splitlines()
    starts = [i for i, line in enumerate(lines) if line.strip() == f"- name: {STEP_NAME}"]
    assert len(starts) == 1, f"ci.yml must have exactly one step named {STEP_NAME!r}"
    start = next(i for i in range(starts[0], len(lines)) if lines[i].strip() == "python - <<'PY'")
    end = next(i for i in range(start + 1, len(lines)) if lines[i].strip() == "PY")
    body = lines[start + 1:end]
    indent = min(len(line) - len(line.lstrip()) for line in body if line.strip())
    return "\n".join(line[indent:] for line in body) + "\n"


def _declared(block: str) -> tuple:
    for node in ast.walk(ast.parse(block)):
        if isinstance(node, ast.Assign) and any(getattr(t, "id", None) == "declared" for t in node.targets):
            return ast.literal_eval(node.value)
    raise AssertionError("the CI step no longer assigns a `declared` tuple")


def _reason_when_the_workspace_is_unset(monkeypatch) -> str:
    """The skip message the live corpus tests give with CORPUS_WORKSPACE unset."""
    monkeypatch.delenv("CORPUS_WORKSPACE", raising=False)
    real_which = shutil.which
    # get past the compiler check, which comes first, on a machine without gfortran
    monkeypatch.setattr(shutil, "which",
                        lambda name, *a, **k: "/usr/bin/gfortran" if name == "gfortran" else real_which(name, *a, **k))
    from test_corpus_residual import _live
    with pytest.raises(pytest.skip.Exception) as skipped:
        _live()
    reason = str(skipped.value.msg)
    assert "CORPUS_WORKSPACE" in reason, f"_live() skipped for another reason: {reason!r}"
    return reason


def _run_step(block: str, reasons, tmp_path: Path) -> subprocess.CompletedProcess:
    paths = re.findall(r'open\("([^"]+)"\)', block)
    assert len(paths) == 1, f"the CI step should read one suite output file, found {paths!r}"
    output = tmp_path / "suite.txt"
    output.write_text("".join(f"SKIPPED [1] tests/corpus/x.py:1: {r}\n" for r in reasons)
                      + "1 passed, %d skipped in 0.01s\n" % len(reasons), encoding="utf-8")
    script = tmp_path / "skip_step.py"
    script.write_text(block.replace(paths[0], str(output)), encoding="utf-8")
    return subprocess.run([sys.executable, str(script)], capture_output=True, text=True)


def test_the_declared_skip_list_names_the_corpus_workspace_reason(monkeypatch):
    reason = _reason_when_the_workspace_is_unset(monkeypatch)
    declared = _declared(_step_python_block())
    assert any(entry in reason for entry in declared), (
        f"the CI skip step would reject the live corpus tests' skip {reason!r}; declared = {declared!r}")


def test_the_ci_skip_step_passes_that_skip_and_still_fails_an_undeclared_one(monkeypatch, tmp_path):
    reason = _reason_when_the_workspace_is_unset(monkeypatch)
    block = _step_python_block()
    accepted = _run_step(block, [reason], tmp_path)
    assert accepted.returncode == 0, accepted.stdout + accepted.stderr
    rejected = _run_step(block, [reason, "a reason no runner should accept"], tmp_path)
    assert rejected.returncode == 1, rejected.stdout + rejected.stderr
    assert "a reason no runner should accept" in rejected.stdout
