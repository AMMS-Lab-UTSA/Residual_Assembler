"""What ``run_report.txt`` says about verification and unsupported features.

Three lines of the history report are verdicts, not echoes: ``Tangent
verified``, ``Reference resolved`` and ``Unsupported feature detected``. Each
is produced here by a real ``resasm history`` run on the committed Abaqus
example (``examples/replay_history/j2_beam``) with the m3_j2 provider built
from UMAT-OTI, and each is checked against the numbers the same run writes to
``sensitivity_results.json``:

* without ``--verify`` the tangent is not verified and no reference is run;
* ``--verify tangent`` spot-checks DDSDDE against central differences of the
  ORIGINAL routine and says ``yes`` with the measured error;
* ``--verify fd`` adds whole-model central differences; with the default step
  ladder adjacent steps agree (a plateau: ``yes``), with a ladder of coarse
  steps they do not (``partially``);
* a deck with a feature the replay does not reproduce stops the run, and the
  report written for the failed run names the feature.
"""
import json
import re

import pytest

from residual_core.ui.cli import main as resasm

from history_support import EXAMPLE

BEAM = EXAMPLE / "j2_beam"


def _history(provider, out, *extra, model=BEAM / "Analysis.inp"):
    obj, _, _ = provider
    return resasm(["history", "--model", str(model), "--fields", str(BEAM / "fields.npz"),
                   "--material", str(obj), "--request", str(BEAM / "sensitivity_request.json"),
                   "--out", str(out), *map(str, extra)])


def _report_line(out, field):
    report = (out / "run_report.txt").read_text()
    prefix = field + ": "
    lines = [line[len(prefix):] for line in report.splitlines() if line.startswith(prefix)]
    assert len(lines) == 1, (field, report)
    return lines[0]


def _verification(out):
    return json.loads((out / "sensitivity_results.json").read_text())["metadata"]["verification"]


@pytest.fixture(scope="module")
def provider(provider_factory):
    return provider_factory("m3_j2")


@pytest.fixture(scope="module")
def fd_run(provider, tmp_path_factory):
    """``--verify fd`` with the default step ladder (about 10 s)."""
    out = tmp_path_factory.mktemp("verify-fd") / "results"
    assert _history(provider, out, "--verify", "fd") == 0
    return out


def test_without_verification_nothing_is_claimed_verified(provider, tmp_path):
    out = tmp_path / "results"
    assert _history(provider, out) == 0
    assert _report_line(out, "Tangent available").startswith("yes: DDSDDE")
    assert _report_line(out, "Tangent verified") == "not run"
    assert _report_line(out, "Derivative verified") == "not run"
    assert _report_line(out, "Reference resolved") == "not applicable (no reference was run)"
    assert _verification(out) == {}


def test_tangent_verified_under_verify_tangent(provider, tmp_path):
    out = tmp_path / "results"
    assert _history(provider, out, "--verify", "tangent") == 0
    line = _report_line(out, "Tangent verified")
    check = _verification(out)["tangent"]
    match = re.fullmatch(r"yes: max relative error (\S+) vs central FD of the ORIGINAL UMAT at "
                         r"(\d+) points \(FD plateau spread (\S+)\)", line)
    assert match, line
    assert match.group(1) == "%.2e" % check["max_relative_error"]
    assert int(match.group(2)) == check["points_checked"] > 0
    assert match.group(3) == "%.2e" % check["fd_plateau_spread"]
    assert check["max_relative_error"] < 1e-5          # the threshold for "yes"
    # a tangent check alone runs no whole-model reference
    assert _report_line(out, "Reference resolved") == "not applicable (no reference was run)"
    assert _report_line(out, "Derivative verified") == "not run"


def test_tangent_verified_under_verify_fd(fd_run):
    line = _report_line(fd_run, "Tangent verified")
    check = _verification(fd_run)["tangent"]
    assert line.startswith("yes: max relative error %.2e " % check["max_relative_error"])


def _worst_spread(verification):
    summary = verification["whole_model_fd"]["summary"]
    return max(row["max_spread"] for rows in summary.values() for row in rows.values()
               if row["nonzero_increments"])


def test_reference_resolved_is_a_plateau_with_the_default_ladder(fd_run):
    verification = _verification(fd_run)
    steps = verification["whole_model_fd"]["steps"]
    spread = _worst_spread(verification)
    assert steps == [1e-3, 3e-4, 1e-4, 3e-5, 1e-5]
    assert spread < 1e-4                                # the threshold for "yes"
    assert _report_line(fd_run, "Reference resolved") == (
        "yes: a plateau (adjacent steps of %s agreeing to %.1e) for every nonzero derivative"
        % (steps, spread))
    assert _report_line(fd_run, "Derivative verified").startswith("yes: whole-model central FD")


def test_reference_resolved_is_partial_when_the_steps_do_not_agree(provider, tmp_path):
    """Coarse steps (30 % and 10 % of each parameter) leave no plateau."""
    out = tmp_path / "results"
    assert _history(provider, out, "--verify", "fd", "--fd-steps", "0.3,0.1") == 0
    spread = _worst_spread(_verification(out))
    assert spread >= 1e-4
    assert _report_line(out, "Reference resolved") == \
        "partially: largest plateau spread %.2e" % spread


@pytest.mark.parametrize("change, feature", [
    (lambda text: text.replace("nlgeom=NO", "nlgeom=YES"), "NLGEOM"),
    (lambda text: text.replace("*End Step", "*Dsload\nALL, P, 1.0\n*End Step"), "Dsload"),
], ids=["NLGEOM", "Dsload"])
def test_an_unsupported_feature_is_reported_as_detected(provider, tmp_path, capsys, change,
                                                         feature):
    deck = tmp_path / "Analysis.inp"
    original = (BEAM / "Analysis.inp").read_text()
    changed = change(original)
    assert changed != original
    deck.write_text(changed)
    out = tmp_path / "results"
    assert _history(provider, out, model=deck) == 2
    assert "history replay failed: " in capsys.readouterr().err
    report = (out / "run_report.txt").read_text()
    assert report.startswith("Status: failed\n")
    line = _report_line(out, "Unsupported feature detected")
    assert line.startswith("yes: ") and feature in line, line
    assert "Error: %s" % line[len("yes: "):] in report
    # nothing was computed, and the report does not say otherwise
    assert _report_line(out, "Residual assembled") == "no"
    assert _report_line(out, "Derivative calculated") == "no"
    assert not (out / "sensitivity_results.json").exists()
