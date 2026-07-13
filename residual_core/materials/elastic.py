"""Requested module name ``elastic`` — re-exports the elastic material backend.

Implementation lives in ``elastic_adapter`` (kept as-is). A real, runnable
small-strain isotropic linear-elastic law: the minimal reference material.
"""

from .elastic_adapter import IsotropicElastic

Elastic = IsotropicElastic

__all__ = ["Elastic", "IsotropicElastic"]
