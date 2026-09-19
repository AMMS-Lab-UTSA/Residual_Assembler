"""Shared helpers for the Abaqus validation scripts.

The scripts in this folder are **ready to run on a machine with Abaqus**, but must
**never fail the offline test suite** when Abaqus is missing. This module provides
the graceful-skip contract used by all of them. It is imported by
``extract_odb_fields.py`` under Abaqus Python 2.7 as well as Python 3, so it
stays free of annotations and f-strings.

Skip contract: when Abaqus (or the in-Abaqus ``odbAccess`` module) is not
available, print ``Abaqus not available: validation pending`` and exit 0.
"""

from __future__ import print_function

import os
import sys

try:
    from shutil import which as _which
except ImportError:  # Python 2.7 (Abaqus Python) has no shutil.which
    from distutils.spawn import find_executable as _which

SKIP_MESSAGE = "Abaqus not available: validation pending"


def abaqus_available():
    """True if an ``abaqus`` executable is on PATH or ABAQUS/ABQ env vars point to one."""
    if _which("abaqus"):
        return True
    for var in ("ABAQUS_CMD", "ABAQUS", "ABQ_CMD"):
        p = os.environ.get(var)
        if p and (_which(p) or os.path.exists(p)):
            return True
    return False


def odb_access_available():
    """True only inside the Abaqus Python interpreter (``odbAccess`` importable)."""
    try:
        import odbAccess  # noqa: F401  (Abaqus-only module)
        return True
    except Exception:
        return False


def skip_if_no_abaqus(need_odb=False):
    """Print the pending message and exit 0 if Abaqus / odbAccess is unavailable."""
    ok = odb_access_available() if need_odb else abaqus_available()
    if not ok:
        print(SKIP_MESSAGE)
        sys.exit(0)


def note_pending(reason):
    """Print a specific 'pending' reason without failing (exit 0)."""
    print("%s (%s)" % (SKIP_MESSAGE, reason))
    sys.exit(0)
