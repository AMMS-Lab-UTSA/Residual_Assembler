"""An installed copy of the producing checkout counts as that checkout; another version does not.

The real-store contract tests compare this repository's reader with the
producer's, and must know that the imported ``umat_oti`` is the companion
checkout's code. They refused a ``umat_oti`` installed (not editable) from
that checkout, as CI and users install it, because its files were not under
the checkout, and skipped (924ee8c). The rule that replaced the location test:
an imported package is the checkout's code when every Python file on either
side exists on the other with the same bytes. Those tests run only where the
private corpus store is present; the rule is checked here on a package laid
out as pip lays it out.
"""
from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from contract_paths import runs_the_checkouts_code

pytestmark = pytest.mark.regression

CONTRACT_CODE = "GATES = ('primal', 'tangent')\n"


def _checkout(root: Path) -> Path:
    package = root / "checkout" / "src" / "umat_oti"
    (package / "contract").mkdir(parents=True)
    (package / "__init__.py").write_text('__version__ = "1.1.0"\n')
    (package / "contract" / "__init__.py").write_text(CONTRACT_CODE)
    (package / "contract" / "schema.json").write_text("{}\n")
    return root / "checkout"


def _copy(checkout: Path, destination: Path) -> Path:
    """``__init__.py`` of a copy of the checkout's package at ``destination``."""
    shutil.copytree(checkout / "src" / "umat_oti", destination)
    return destination / "__init__.py"


def test_a_byte_identical_installed_copy_runs_the_checkouts_code(tmp_path):
    checkout = _checkout(tmp_path)
    assert runs_the_checkouts_code(str(checkout / "src" / "umat_oti" / "__init__.py"), checkout)
    installed = _copy(checkout, tmp_path / "env" / "lib" / "python3" / "site-packages" / "umat_oti")
    assert runs_the_checkouts_code(str(installed), checkout)
    # package data is not code: an install may carry it where the source tree does not
    (installed.parent / "contract" / "schema.json").unlink()
    (installed.parent / "README.txt").write_text("installed\n")
    # nor is a byte-compiled cache
    (installed.parent / "__pycache__").mkdir()
    (installed.parent / "__pycache__" / "stale.py").write_text("x = 1\n")
    assert runs_the_checkouts_code(str(installed), checkout)


def test_an_installed_copy_of_another_version_is_refused(tmp_path):
    checkout = _checkout(tmp_path)
    installed = _copy(checkout, tmp_path / "env" / "lib" / "python3" / "site-packages" / "umat_oti")
    contract = installed.parent / "contract" / "__init__.py"
    contract.write_text("GATES = ('primal',)\n")
    assert not runs_the_checkouts_code(str(installed), checkout)
    contract.write_text(CONTRACT_CODE)
    assert runs_the_checkouts_code(str(installed), checkout)
    added = installed.parent / "added.py"
    added.write_text("")
    assert not runs_the_checkouts_code(str(installed), checkout)
    added.unlink()
    contract.unlink()
    assert not runs_the_checkouts_code(str(installed), checkout)
    # a stale copy under the checkout is under it, and still another version
    stale = _copy(checkout, checkout / "build" / "lib" / "umat_oti")
    (stale.parent / "contract" / "__init__.py").write_text("GATES = ()\n")
    assert not runs_the_checkouts_code(str(stale), checkout)
    assert not runs_the_checkouts_code(None, checkout)
