"""The small-strain C3D8 never drives a finite-strain material (and says why).

``solid_c3d8_small_strain`` hands its material ``{"strain", "dstrain"}`` and
integrates the returned stress with B0 over the reference volume. A material
that expects a deformation gradient (the compressible neo-Hookean solid, the
UMAT adapter) or returns another stress or tangent measure would be fed the
wrong kinematics or read in the wrong measure. The element refuses it with a
ValueError naming the element, the material and the declared measure; the
finite-strain side of the same rule is pinned in
test_finite_strain_hyperelastic.py (test_reject_incompatible_material).
"""
import numpy as np
import pytest

from residual_core.formulations import c3d8_kernel as kernel
from residual_core.formulations.solid_c3d8_small_strain import SolidC3D8SmallStrain
from residual_core.materials.base import MaterialBinding
from residual_core.materials.elastic_adapter import IsotropicElastic
from residual_core.materials.neo_hookean import CompressibleNeoHookean
from residual_core.materials.umat_adapter import UmatAdapter, mock_isotropic_elastic_umat


def small_strain(material, constants, displacement=None):
    coords = kernel.unit_cube_Xe()
    u = np.zeros(24) if displacement is None else displacement
    return SolidC3D8SmallStrain().eval_element(
        1, "C3D8", coords, u, {}, None, MaterialBinding(material, list(constants)),
        (0.0, 1.0), 1.0, None, {"compute_tangent": True})


def test_a_deformation_gradient_material_is_refused_by_name():
    with pytest.raises(ValueError) as refused:
        small_strain(CompressibleNeoHookean(), [2.3, 4.1])
    message = str(refused.value)
    assert message.startswith("solid_c3d8_small_strain: material compressible_neo_hookean "
                              "has kinematic_input='deformation_gradient'; "
                              "requires 'small_strain'")
    assert "solid_c3d8_finite_strain" in message


def test_the_umat_adapter_is_refused_before_it_is_called():
    calls = []

    def umat(**arguments):
        calls.append(arguments)
        return mock_isotropic_elastic_umat(**arguments)

    with pytest.raises(ValueError, match="kinematic_input='deformation_gradient'"):
        small_strain(UmatAdapter(backend="python", umat_fn=umat), [210000.0, 0.3])
    assert calls == []


@pytest.mark.parametrize("attribute,value", [
    ("kinematic_input", "deformation_gradient"),
    ("stress_measure", "pk2"),
    ("stress_measure", "pk1"),
    ("tangent_measure", "kirchhoff_jaumann_over_j"),
    ("tangent_measure", "material"),
])
def test_every_incompatible_declared_measure_is_refused(attribute, value):
    material = IsotropicElastic()
    setattr(material, attribute, value)
    with pytest.raises(ValueError, match="has %s=%r; requires" % (attribute, value)):
        small_strain(material, [210000.0, 0.3])


def test_the_small_strain_material_still_assembles():
    """The guard refuses nothing it should accept: linear elasticity gives
    R = K u exactly, with K the independent kernel stiffness."""
    E, nu = 210000.0, 0.3
    coords = kernel.unit_cube_Xe()
    displacement = 1e-3 * np.sin(np.arange(24.0))
    force, tangent, _, _ = small_strain(IsotropicElastic(), [E, nu], displacement)
    from residual_core.core.voigt import isotropic_D
    reference = kernel.element_tangent(coords, None, isotropic_D(E, nu), mode="small")
    np.testing.assert_allclose(tangent, reference, rtol=1e-12, atol=1e-6)
    np.testing.assert_allclose(force, reference @ displacement, rtol=1e-12, atol=1e-9)


def test_resasm_assemble_reports_the_refusal(tmp_path, capsys):
    from residual_core import ResidualProblem
    from residual_core.core.model import Element, Model
    from residual_core.ui import cli

    coords = kernel.unit_cube_Xe()
    model = Model(nodes={index + 1: tuple(point) for index, point in enumerate(coords)})
    model.elements = {1: Element(1, "C3D8", list(range(1, 9)))}
    model.element_formulation = {1: "solid_c3d8_small_strain"}
    model.element_material = {1: "solid"}
    model.materials = {"solid": MaterialBinding(CompressibleNeoHookean(), [2.3, 4.1],
                                                name="solid")}
    path = tmp_path / "model.json"
    ResidualProblem(model).save_neutral(str(path))
    assert cli.main(["assemble", str(path), "--mode", "material-replay"]) == 2
    err = capsys.readouterr().err
    assert "Cannot assemble in material-replay mode." in err
    assert "solid_c3d8_small_strain: material compressible_neo_hookean has " \
           "kinematic_input='deformation_gradient'" in err
