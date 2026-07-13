"""Material registry — name -> Material instance, with capability specs.

Builds a ``core.registry.Registry`` of the constitutive backends the framework
ships. The material system is NOT UMAT-centered: an entry may be a continuum
stress-strain law, a UMAT bridge, a crystal-plasticity return map, or (in future)
a beam section law, cohesive traction-separation law, or thermal/diffusion flux
law. Each declares its capabilities via ``.spec``.
"""

from __future__ import annotations

from ..core.registry import Registry
from .elastic_adapter import IsotropicElastic


def build_material_registry() -> Registry:
    reg = Registry("material")
    reg.register(IsotropicElastic())
    # optional adapters that may need compiled Fortran / Abaqus
    for key, module, cls in (
            ("umat", "residual_core.materials.umat_adapter", "UmatAdapter"),
            ("crystal_plasticity",
             "residual_core.materials.crystal_plasticity_adapter",
             "CrystalPlasticityAdapter")):
        try:
            mod = __import__(module, fromlist=[cls])
            reg.register(getattr(mod, cls)(), name=key)
        except Exception:
            pass
    return reg


def default_material_registry() -> Registry:
    return build_material_registry()
