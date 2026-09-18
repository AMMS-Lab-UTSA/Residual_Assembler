"""The GUI must not be able to disagree with the CLI.

The Streamlit front end never reimplements a command: every button builds an
argv list and hands it to ``residual_core.ui.cli.main``. These tests pin that
property, because the moment the GUI grows its own copy of an assembly call it
can start reporting a number the CLI would refuse to produce.

Streamlit itself is an optional extra (``pip install -e ".[gui]"``), so the
import-dependent tests SKIP cleanly when it is absent; the source-level checks
run always.

Run: python -m pytest tests/framework/test_gui_is_a_thin_cli_front_end.py
"""

import os
import sys

import pytest

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

_APP_SOURCE = os.path.join(_ROOT, "residual_core", "app", "streamlit_app.py")
_EXAMPLE = os.path.join(_ROOT, "residual_core", "examples",
                        "minimal_c3d8_stress_driven", "model.json")


def _app():
    pytest.importorskip("streamlit", reason="GUI extra not installed: pip install -e '.[gui]'")
    from residual_core.app import streamlit_app

    return streamlit_app


@pytest.mark.unit
def test_the_app_ships_with_the_package():
    """A GUI that is not packaged is a GUI nobody can start."""
    assert os.path.isfile(_APP_SOURCE)
    assert os.path.isfile(os.path.join(_ROOT, "scripts", "app.py"))


@pytest.mark.unit
def test_the_app_does_not_import_the_assembler_directly():
    """No assembly, no solve, no residual construction inside the GUI.

    The front end may import the *registries* (to fill a dropdown with the
    names that exist) and ``otilib_status`` (to tell the truth in the sidebar).
    It must not import the machinery that produces numbers -- that has to be
    reached through ``cli.main``, so the GUI cannot answer a question the CLI
    would refuse. Function-local imports count; this walks the whole AST.
    """
    import ast

    allowed = {
        "residual_core.ui.cli",              # the bridge itself
        "residual_core.core.requirements",   # mode names for a dropdown
        "residual_core.formulations.registry",
        "residual_core.materials.registry",
        "residual_core.algebra.otilib_adapter",  # honest sidebar status
        # The Solve screen's parameter/output/region widgets (slide 18). The
        # module is part of the GUI and is walked below by the same rule; it
        # may read the .inp only to offer the .inp's own sets in a dropdown.
        "residual_core.app.request_screen",
        "residual_core.io.abaqus_inp_parser",
    }
    sources = [_APP_SOURCE, os.path.join(_ROOT, "residual_core", "app", "request_screen.py")]
    imported = set()
    for source in sources:
        with open(source, "r", encoding="utf-8") as fh:
            tree = ast.parse(fh.read())
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.update(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module and not node.level:
                imported.add(node.module)
                # `from residual_core.app import request_screen` names the
                # module in the alias, not in node.module.
                imported.update(f"{node.module}.{alias.name}" for alias in node.names
                                if node.module == "residual_core.app")

    # `from residual_core.ui import cli` records the package, not the module,
    # so every parent of an allowed module is allowed too.
    permitted = set(allowed)
    for name in allowed:
        parts = name.split(".")
        permitted.update(".".join(parts[:i]) for i in range(1, len(parts)))

    offenders = sorted(
        name for name in imported
        if name.split(".")[0] in {"residual_core", "resasm_user", "partner_kit"}
        and name not in permitted)
    assert not offenders, (
        "the GUI imports framework machinery directly: %s -- reach it through "
        "the CLI instead" % ", ".join(offenders))


@pytest.mark.unit
def test_every_button_goes_through_cli_main():
    app = _app()
    result = app._run(["modes"])
    assert result.code == 0
    assert "stress-driven" in result.stdout
    assert result.command.startswith("resasm modes")


@pytest.mark.unit
def test_a_failing_command_reports_its_real_exit_code():
    """Exit codes are data, not decoration: 2 means 'could not run'."""
    app = _app()
    result = app._run(["inspect", os.path.join(_ROOT, "no_such_model.inp")])
    assert result.code != 0
    assert result.stdout.strip() or result.stderr.strip()


@pytest.mark.unit
def test_usage_errors_do_not_kill_the_app():
    """argparse calls sys.exit(); the GUI has to survive that and say so."""
    app = _app()
    result = app._run(["definitely-not-a-command"])
    assert result.code == 2


@pytest.mark.unit
def test_the_working_directory_is_always_restored():
    """A command that fails mid-flight must not leave the process elsewhere."""
    app = _app()
    before = os.getcwd()
    app._run(["inspect", "missing.inp"], cwd=os.path.dirname(_EXAMPLE))
    assert os.getcwd() == before


@pytest.mark.integration
def test_the_gui_assembles_the_same_residual_as_the_cli():
    app = _app()
    fields = os.path.join(os.path.dirname(_EXAMPLE), "fields.json")
    result = app._run(["assemble", _EXAMPLE, "--mode", "stress-driven",
                       "--fields", fields])
    assert result.code == 0
    assert "assembled mode 'stress-driven'" in result.stdout
    assert "ndof=24" in result.stdout


@pytest.mark.unit
def test_the_shipped_examples_are_discoverable():
    """The example dropdown must not be empty in a fresh checkout."""
    app = _app()
    examples = app._shipped_examples()
    assert "examples/minimal_c3d8_stress_driven" in examples
    assert all(path.is_file() for path in examples.values())


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
