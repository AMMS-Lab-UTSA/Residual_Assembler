"""Requested module name ``crystal_plasticity`` — re-exports the CP material.

Implementation lives in ``crystal_plasticity_adapter`` (kept as-is). Crystal
plasticity is ONE example material backend here, not the organizing principle of
the framework. It is a deformation-gradient-driven UMAT (Cauchy stress + DDSDDE +
STATEV) and needs Intel ifort + Abaqus to compile the Fortran kernel.
"""

from .crystal_plasticity_adapter import CrystalPlasticityAdapter

CrystalPlasticity = CrystalPlasticityAdapter

__all__ = ["CrystalPlasticity", "CrystalPlasticityAdapter"]
