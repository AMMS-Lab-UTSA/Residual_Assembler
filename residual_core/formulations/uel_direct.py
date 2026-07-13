"""Requested module name ``uel_direct`` — re-exports the UEL/direct adapter.

The implementation lives in ``uel_adapter`` (kept as-is). Mode 3: wrap an
Abaqus-UEL-like callable that returns the element RHS/AMATRX directly, behind the
standard residual Formulation contract (sign convention handled + verified).
"""

from .uel_adapter import UelAdapter

# Canonical short alias
UelDirect = UelAdapter

__all__ = ["UelDirect", "UelAdapter"]
