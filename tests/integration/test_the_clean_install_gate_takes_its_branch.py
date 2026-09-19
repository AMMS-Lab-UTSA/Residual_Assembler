"""The clean-install gate runs on the branch it is told, and refuses any other.

The gate compared each checkout's branch with one hard-coded branch name and
refused everything else, so it could not run on main, on a fresh clone of
another branch, or on a new branch (5e49cc9). The branch is now ``--branch``;
without it any branch is accepted, and with it a checkout on another branch is
still refused.

Each test builds two real git repositories and runs the gate's ``main`` on
them. The gate is given an Abaqus launcher that does not exist, so a gate that
accepts both checkouts stops at its first external requirement, after the
branch check, and never builds a wheel.
"""
from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.regression

SPEC = importlib.util.spec_from_file_location(
    "clean_install_gate", Path(__file__).resolve().parents[2] / "scripts" / "clean_install_gate.py")
gate = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(gate)

MISSING_ABAQUS = "abaqus-that-is-not-installed"
PAST_THE_BRANCH_CHECK = "required external executable not found"


def _git(repository: Path, *arguments: str) -> str:
    return subprocess.run(["git", "-C", str(repository), *arguments], check=True,
                          capture_output=True, text=True).stdout


def _repository(path: Path, branch: str) -> Path:
    """A repository with one commit, a clean tree and ``branch`` checked out."""
    path.mkdir()
    _git(path, "init", "-q")
    _git(path, "symbolic-ref", "HEAD", f"refs/heads/{branch}")
    (path / "README.md").write_text(f"{path.name}\n")
    _git(path, "add", "README.md")
    _git(path, "-c", "user.name=Gate test", "-c", "user.email=gate-test@example.invalid",
         "commit", "-q", "-m", "one commit")
    assert _git(path, "branch", "--show-current").strip() == branch
    return path


def _run(tmp_path: Path, ra_branch: str, umat_branch: str, *options: str) -> tuple:
    ra = _repository(tmp_path / "ra", ra_branch)
    umat = _repository(tmp_path / "umat", umat_branch)
    work = tmp_path / "work"
    status = gate.main(["--ra-repo", str(ra), "--umat-repo", str(umat), "--python", sys.executable,
                        "--odb", str(tmp_path / "unused.odb"), "--abaqus", MISSING_ABAQUS,
                        "--work", str(work), *options])
    return status, json.loads((work / "report.json").read_text()), ra.resolve(), umat.resolve()


@pytest.mark.parametrize("branch", ["main", "release/2027-01"])
def test_the_named_branch_is_accepted(tmp_path, branch):
    status, report, ra, umat = _run(tmp_path, branch, branch, "--branch", branch)
    assert status == 1 and report["passed"] is False
    assert report["error"].startswith(PAST_THE_BRANCH_CHECK), report["error"]
    assert {path: entry["branch"] for path, entry in report["repositories"].items()} == \
        {str(ra): branch, str(umat): branch}


def test_without_a_branch_a_checkout_on_main_is_accepted(tmp_path):
    status, report, ra, umat = _run(tmp_path, "main", "main")
    assert status == 1
    assert report["error"].startswith(PAST_THE_BRANCH_CHECK), report["error"]
    assert set(report["repositories"]) == {str(ra), str(umat)}
    # no branch was named, so nothing is claimed about a published branch
    assert report["final_branch_clean_clone"] is False


def test_a_checkout_on_another_branch_is_refused(tmp_path):
    status, report, ra, umat = _run(tmp_path, "main", "topic", "--branch", "main")
    assert status == 1 and report["passed"] is False
    assert report["error"] == f"expected branch main, found topic at {umat}"
    # refused before anything else ran: no interpreter was started
    assert all(entry["argv"][0] == "git" for entry in report["commands"])
