"""Versioned interface surface for the residual-assembler.

`versions` is the single source of truth for every JSON-contract schema tag and
the shared C-ABI / material-package contract version, so consumers can pin them
and drift is caught in one place instead of across a dozen string literals.
"""
from . import versions  # noqa: F401

__all__ = ["versions"]
