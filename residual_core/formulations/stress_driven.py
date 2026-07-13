"""Requested module name ``stress_driven`` — re-exports the stress-driven backend.

The implementation lives in ``stress_driven_adapter`` (kept as-is). Mode 1:
assemble the residual from an externally exported integration-point field, with
no material update.
"""

from .stress_driven_adapter import StressDrivenC3D8

# Canonical short alias
StressDriven = StressDrivenC3D8

__all__ = ["StressDriven", "StressDrivenC3D8"]
