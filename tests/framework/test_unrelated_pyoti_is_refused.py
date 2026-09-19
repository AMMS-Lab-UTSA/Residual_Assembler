"""An unrelated package named ``pyoti`` is detected and refused, by name.

OTILib's Python bindings are a package called ``pyoti`` built from source; the
Python Package Index hosts an unrelated package under the same name. A user
who runs ``pip install pyoti`` gets that package, and because it sits in
site-packages it would shadow a genuine build found later on the path.

``residual_core.algebra.otilib_adapter`` accepts a module only if it exposes
the OTI number API (``e``, ``number``, ``get_im``). These tests put a stand-in
for the unrelated package (a ``pyoti`` package with none of that API and no
``pyoti.sparse``) AHEAD of everything else on the path, in a fresh
interpreter, and check that it is refused with a message naming its file,
both by the adapter and by ``resasm check`` of a job that needs OTILib. The
stand-in is the only thing faked; the probe and the command are the real ones.
"""
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from residual_core.algebra import otilib_adapter

ROOT = Path(__file__).resolve().parents[2]

PROBE = ("import json; from residual_core.algebra import otilib_adapter as A; "
         "status = A.otilib_status(); "
         "print(json.dumps({'available': A.otilib_available(), "
         "'error': status['error'], 'api': status['api_module']}))")


def _unrelated_pyoti(tmp_path):
    """A directory holding a package named pyoti that is not OTILib."""
    package = tmp_path / "site" / "pyoti"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text(
        '"""Stand-in for the unrelated PyPI project that is also called pyoti."""\n'
        '__version__ = "0.0.1"\n\n\n'
        'def describe():\n'
        '    return "an unrelated package"\n')
    return package.parent, package / "__init__.py"


def _environment(path_entries):
    """This interpreter's environment with PYTHONPATH replaced and no
    OTILIB_ROOT/PYOTI_PATH, so only the path decides what ``pyoti`` is."""
    env = {key: value for key, value in os.environ.items()
           if key not in ("OTILIB_ROOT", "PYOTI_PATH", "PYTHONPATH")}
    env["PYTHONPATH"] = os.pathsep.join(str(entry) for entry in path_entries)
    return env


def _genuine_build_dirs():
    """Where this environment's genuine OTILib build is, if it names one."""
    return [value for value in (os.environ.get("PYOTI_PATH"), os.environ.get("OTILIB_ROOT"))
            if value and os.path.isdir(value)]


def _probe(path_entries):
    result = subprocess.run([sys.executable, "-c", PROBE], cwd=str(ROOT),
                            env=_environment(path_entries), capture_output=True,
                            text=True, check=False)
    assert result.returncode == 0, result.stdout + result.stderr
    return json.loads(result.stdout)


@pytest.fixture
def no_local_checkout():
    """The adapter also looks in these checkouts; the tests need them absent."""
    local = [ROOT / relative for relative in otilib_adapter._LOCAL_ROOTS]
    present = [str(path) for path in local if path.exists()]
    if present:
        pytest.skip("a local OTILib checkout would be probed first: %s" % present)


def test_an_unrelated_pyoti_ahead_of_the_genuine_build_is_refused(tmp_path, no_local_checkout):
    site, init = _unrelated_pyoti(tmp_path)
    status = _probe([site, ROOT] + _genuine_build_dirs())
    assert status["available"] is False
    assert status["api"] == ""
    error = status["error"]
    assert error.startswith("Refused a module named pyoti that is not OTILib: ")
    assert str(init) in error
    assert "pip uninstall pyoti" in error
    assert "genuine OTILib was not found" in error


def test_without_the_unrelated_package_the_error_names_no_impostor(tmp_path, no_local_checkout):
    """The refusal message is specific: an empty path gives the plain one."""
    status = _probe([ROOT])
    if status["available"]:
        pytest.skip("a genuine OTILib is installed in site-packages")
    assert "Refused a module named pyoti" not in status["error"]
    assert "genuine OTILib was not found" in status["error"]


@pytest.mark.skipif(not otilib_adapter.otilib_available(),
                    reason="needs the genuine OTILib build to show it is still accepted")
def test_the_genuine_build_is_accepted_when_it_comes_first(tmp_path, no_local_checkout):
    site, _ = _unrelated_pyoti(tmp_path)
    genuine = _genuine_build_dirs() or [
        str(Path(otilib_adapter._OTI.__file__).resolve().parents[1])]
    status = _probe(genuine + [site, ROOT])
    assert status["available"] is True, status["error"]
    assert status["api"].startswith("pyoti")


def test_resasm_check_refuses_a_job_when_only_the_unrelated_package_is_found(tmp_path,
                                                                           no_local_checkout):
    site, init = _unrelated_pyoti(tmp_path)
    env = _environment([site, ROOT] + _genuine_build_dirs())
    job = tmp_path / "job"
    created = subprocess.run(
        [sys.executable, "-m", "residual_core.ui.cli", "init", "--template", "python",
         "--out", str(job)], cwd=str(ROOT), env=env, capture_output=True, text=True,
        check=False)
    assert created.returncode == 0, created.stdout + created.stderr
    checked = subprocess.run(
        [sys.executable, "-m", "residual_core.ui.cli", "check", str(job / "resasm.yml")],
        cwd=str(ROOT), env=env, capture_output=True, text=True, check=False)
    assert checked.returncode == 1, checked.stdout + checked.stderr
    assert "Refused a module named pyoti that is not OTILib" in checked.stdout
    assert str(init) in checked.stdout
