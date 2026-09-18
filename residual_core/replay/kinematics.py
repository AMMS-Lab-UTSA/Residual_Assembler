"""C3D8 replay operators recovered from the inspected W2 implementation.

The B-bar formulation lives here so recovery does not modify the W1 kernels.
Material arrays must already be extracted real ABI coefficients.
"""

import numpy as np

from ..formulations import c3d8_kernel as _k


def batch_operators(coordinates, integration="selective_reduced"):
    if integration not in ("full", "selective_reduced"):
        raise ValueError(f"unsupported C3D8 integration: {integration}")
    matrices, weights = [], []
    for element in coordinates:
        operators, determinants = zip(*[
            _k.b_matrix_reference(element, point)
            for point in _k.ABAQUS_C3D8_GAUSS.points])
        operators = np.array(operators)
        weight = np.asarray(determinants) * _k.ABAQUS_C3D8_GAUSS.weights
        if not np.all(np.isfinite(weight)) or np.any(weight <= 0):
            raise ValueError("C3D8 requires finite positive Jacobians at every IP")
        if integration == "selective_reduced":
            volume_row = operators[:, :3, :].sum(axis=1)
            average = np.einsum("qi,q->i", volume_row, weight) / weight.sum()
            operators[:, :3, :] += ((average - volume_row) / 3)[:, None, :]
        matrices.append(operators)
        weights.append(weight)
    return np.array(matrices), np.array(weights)


def batch_strain(operators, displacement):
    return np.einsum("eqai,ei->eqa", operators, displacement)


def batch_element_tangent(operators, tangent, weights):
    return np.einsum("eqai,eqab,eqbj,eq->eij", operators, tangent, operators, weights)


def batch_element_force(operators, field, weights):
    field = np.asarray(field)
    if field.dtype.kind not in "fiu":
        raise ValueError("replay assembly requires extracted real coefficients, not live OTI")
    if field.ndim == 3:
        return np.einsum("eqai,eqa,eq->ei", operators, field, weights)
    if field.ndim == 4:
        return np.einsum("eqai,eqam,eq->eim", operators, field, weights)
    raise ValueError("field must have shape (ne,8,6) or (ne,8,6,nparam)")


def element_tangent_bbar(coordinates, tangent):
    operators, weights = batch_operators([coordinates])
    tangent = np.broadcast_to(tangent, (8, 6, 6))
    return batch_element_tangent(operators, tangent[None], weights)[0]


def element_internal_force_bbar(coordinates, stress):
    operators, weights = batch_operators([coordinates])
    return batch_element_force(operators, np.asarray(stress)[None], weights)[0]


def bbar_b_matrices(coordinates):
    operators, weights = batch_operators([coordinates])
    return list(operators[0]), list(weights[0] / _k.ABAQUS_C3D8_GAUSS.weights)