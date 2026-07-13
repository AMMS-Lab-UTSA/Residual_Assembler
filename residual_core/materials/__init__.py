"""Material backends + a default registry.

Registry maps a material name -> a Material instance. Bindings (constants +
state size) are attached per model material via materials.base.MaterialBinding.
"""

from .base import Material, MaterialBinding
from .elastic_adapter import IsotropicElastic
from .registry import build_material_registry, default_material_registry


def default_materials():
    reg = {"isotropic_elastic": IsotropicElastic()}
    for name, cls in (("umat", "UmatAdapter"),
                      ("crystal_plasticity", "CrystalPlasticityAdapter")):
        try:
            mod = __import__(
                "residual_core.materials.%s" % (
                    "umat_adapter" if name == "umat" else "crystal_plasticity_adapter"),
                fromlist=[cls])
            reg[name] = getattr(mod, cls)()
        except Exception:
            pass
    return reg


__all__ = ["Material", "MaterialBinding", "IsotropicElastic", "default_materials",
           "build_material_registry", "default_material_registry"]
