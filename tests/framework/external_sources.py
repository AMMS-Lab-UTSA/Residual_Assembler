"""Clear diagnostics for optional external source submodules.

External Abaqus sources are pinned as git submodules under `sources/`. A fresh
clone does not contain their files, so tests that read one must say *exactly*
that -- and say how to fix it -- instead of surfacing a bare `FileNotFoundError`
or an assertion about some unrelated downstream value.

Use `require_external_file()` at the top of any test that reads a submodule file.
"""

import os

import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))

# Which submodule owns which path prefix, and how a user is allowed to get it.
_PERMISSIVE_HINT = (
    "It is a permissive (MIT/BSD-3) submodule pinned in .gitmodules. Fetch it with:\n"
    "    ./scripts/init_permissive_sources.sh\n"
    "(or: git submodule update --init -- {path})"
)
_RESTRICTED_HINT = (
    "It is a reference-only submodule (copyleft or no-license-granted) and is\n"
    "deliberately marked `update = none` in .gitmodules, so setup never fetches\n"
    "it. No offline test may depend on it."
)


def submodule_for(path):
    """Return (submodule_path, tier) owning `path`, or (None, None)."""
    rel = os.path.relpath(os.path.abspath(path), REPO_ROOT).replace(os.sep, "/")
    if not rel.startswith("sources/"):
        return None, None
    parts = rel.split("/")
    if len(parts) < 3:
        return None, None
    tier = parts[1]
    return "/".join(parts[:3]), tier


def require_external_file(path, why):
    """Skip with an actionable message if an external-source file is absent.

    `why` states what the test needs the file for, so the skip line explains the
    lost coverage rather than just naming a path.
    """
    if os.path.exists(path):
        return path

    submodule, tier = submodule_for(path)
    rel = os.path.relpath(path, REPO_ROOT).replace(os.sep, "/")
    if submodule is None:
        pytest.skip("missing test input %s (needed for: %s)" % (rel, why))

    hint = _RESTRICTED_HINT if tier in ("copyleft", "license-unknown") else \
        _PERMISSIVE_HINT.format(path=submodule)
    pytest.skip(
        "external source not initialized: %s\n"
        "Needed for: %s\n"
        "Owning submodule: %s (tier: %s)\n"
        "%s" % (rel, why, submodule, tier, hint)
    )
