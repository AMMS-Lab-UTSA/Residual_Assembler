"""Every worked example documents what a reader needs to reproduce and judge it.

``examples/README.md`` lists the worked examples, each with a WALKTHROUGH.md.
Each walkthrough must carry: its inputs, the command to run, the expected
output, how the result is checked independently, a tolerance or a measured
agreement, the files it writes (at least one of them machine-readable), the
failure diagnostics and the GUI equivalent. These are sections of every
walkthrough, under the headings below; the test fails naming the example and
the missing item.
"""
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
INDEX = ROOT / "examples" / "README.md"

#: required item -> the walkthrough heading that carries it
SECTIONS = {
    "inputs": "Inputs",
    "command": "Run it from the command line",
    "expected output": "Expected output",
    "numerical verification": "How the result is checked independently",
    "files and report": "What it writes",
    "failure diagnostics": "Common problems",
    "GUI equivalent": "Run it in the GUI",
}
MACHINE_READABLE = re.compile(r"[\w./$-]+\.(?:json|csv|npz|npy)\b")
NUMBER = re.compile(r"\b\d(?:\.\d+)?e[-+]?\d+\b")
AGREEMENT = re.compile(r"toleran|agree|within|error|difference", re.IGNORECASE)
COMMAND = re.compile(r"^\s*(?:resasm|python|umat-oti-provider|abaqus|streamlit)\s", re.MULTILINE)


def listed_walkthroughs():
    """The WALKTHROUGH.md files linked from the example table of the index."""
    table = INDEX.read_text().split("## Before you run an example", 1)[0]
    links = re.findall(r"^\| \d+ \| \[[^\]]+\]\(([^)]+/WALKTHROUGH\.md)\)", table, re.MULTILINE)
    return [(INDEX.parent / link).resolve() for link in links]


def sections(path):
    """{heading: body} of the level-2 sections of a walkthrough (fences kept)."""
    out, heading, body, fenced = {}, None, [], False
    for line in path.read_text().splitlines():
        if line.startswith("```"):
            fenced = not fenced
        if not fenced and line.startswith("## "):
            if heading is not None:
                out[heading] = "\n".join(body)
            heading, body = line[3:].strip(), []
        elif heading is not None:
            body.append(line)
    if heading is not None:
        out[heading] = "\n".join(body)
    return out


def _name(path):
    return str(path.relative_to(ROOT))


def test_the_index_lists_seven_walkthroughs_that_exist():
    walkthroughs = listed_walkthroughs()
    assert len(walkthroughs) == 7, walkthroughs
    missing = [_name(path) for path in walkthroughs if not path.is_file()]
    assert not missing, missing


@pytest.mark.parametrize("walkthrough", listed_walkthroughs(), ids=_name)
def test_each_walkthrough_carries_every_required_section(walkthrough):
    found = sections(walkthrough)
    missing = [item for item, heading in SECTIONS.items()
               if not any(title.startswith(heading) for title in found)]
    assert not missing, "%s lacks: %s" % (_name(walkthrough), ", ".join(missing))


def _section(found, heading):
    return next(body for title, body in found.items() if title.startswith(heading))


@pytest.mark.parametrize("walkthrough", listed_walkthroughs(), ids=_name)
def test_each_walkthrough_gives_a_command_files_and_a_measured_agreement(walkthrough):
    found = sections(walkthrough)
    name = _name(walkthrough)
    command = _section(found, SECTIONS["command"])
    assert "```" in command and COMMAND.search(command), \
        "%s: no runnable command in its command-line section" % name
    written = _section(found, SECTIONS["files and report"])
    assert MACHINE_READABLE.search(written), \
        "%s: its files section names no machine-readable output" % name
    judged = (_section(found, SECTIONS["expected output"]) + "\n"
              + _section(found, SECTIONS["numerical verification"]))
    assert NUMBER.search(judged) and AGREEMENT.search(judged), \
        "%s: no tolerance or measured agreement in its output and checks" % name
    gui = _section(found, SECTIONS["GUI equivalent"])
    assert "**" in gui or "no GUI" in gui, "%s: its GUI section names no control" % name
    problems = _section(found, SECTIONS["failure diagnostics"])
    assert problems.count("\n|") >= 3, "%s: its problems section lists no symptom" % name
