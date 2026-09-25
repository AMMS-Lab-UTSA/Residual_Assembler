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

# ---------------------------------------------------------------------------
# Finite strain.
#
# Everything above is the REFERENCE configuration: batch_operators builds B
# once from b_matrix_reference and it never moves. That is correct for a
# small-strain replay and silently wrong for an Abaqus NLGEOM=YES analysis,
# where the operators belong in the configuration the increment ends in.
#
# Two separate corrections are needed and BOTH have to be applied. Correcting
# only the kinematics is the easy half, and the state it leaves is the one
# hardest to notice, because the stress field then looks right.
# ---------------------------------------------------------------------------

def deformation_gradients(coordinates, displacement, fbar=True):
    """F at every integration point, (ne, 8, 3, 3).

    With ``fbar`` the volumetric part is taken from the element centroid,

        F_bar = (J_centroid / J_ip)^(1/3)  F_ip

    which is what Abaqus's fully integrated first-order bricks hand the UMAT:
    C3D8 uses selectively reduced integration, evaluating the volumetric
    strain at the centroid so the element does not lock in near-incompressible
    plasticity. Passing F_ip instead leaves the deviatoric response correct and
    the volumetric response wrong, which is close to invisible -- measured on a
    160-element cantilever, the shear components agreed with the recorded
    stress to 0.4% while the stress as a whole was 26% out, by a near-uniform
    hydrostatic offset. With F_bar the same comparison is 7.3e-05.
    """
    coordinates = np.asarray(coordinates, dtype=float)
    displacement = np.asarray(displacement, dtype=float)
    points = _k.ABAQUS_C3D8_GAUSS.points
    centre = np.zeros(3)
    out = np.empty(coordinates.shape[:1] + (len(points), 3, 3))
    for e, (element, move) in enumerate(zip(coordinates, displacement)):
        if fbar:
            _, _, grad_c = _k._b_from_coords(element, centre)
            centroid = np.linalg.det(np.eye(3) + move.T @ grad_c)
        for q, point in enumerate(points):
            _, _, grad = _k._b_from_coords(element, point)
            F = np.eye(3) + move.T @ grad
            if fbar:
                F = F * (centroid / np.linalg.det(F)) ** (1.0 / 3.0)
            out[e, q] = F
    return out


def spatial_operators(coordinates, displacement, integration="selective_reduced"):
    """B and w*detJ in the CURRENT configuration, (ne,8,6,24) and (ne,8).

    The companion of :func:`deformation_gradients`. If the material's strain
    measure is the modified one, so is the virtual work:

        internal force = integral( B_bar^T sigma )
        dR/dp          = integral( B_bar^T dsigma/dp )
        K              = integral( B_bar^T c B_bar ) + initial stress

    Keeping plain B in the assembly while the material responds to F_bar drops
    the centroid coupling out of dsigma/du and makes the tangent non-symmetric
    besides. Measured against a central difference of R(u) on the same
    cantilever: 2.78e-01 with plain B, 3.74e-04 with B_bar -- and at the point
    where K was still 28% wrong the stress already matched the recording to
    7e-05 and the equilibrium residual was already inside Abaqus's own
    convergence tolerance. Neither of the checks one would naturally trust
    caught it.
    """
    if integration not in ("full", "selective_reduced"):
        raise ValueError(f"unsupported C3D8 integration: {integration}")
    coordinates = np.asarray(coordinates, dtype=float)
    displacement = np.asarray(displacement, dtype=float)
    current = coordinates + displacement
    matrices, weights = [], []
    for element in current:
        operators, determinants = zip(*[
            _k.b_matrix_spatial(element, point)
            for point in _k.ABAQUS_C3D8_GAUSS.points])
        operators = np.array(operators)
        weight = np.asarray(determinants) * _k.ABAQUS_C3D8_GAUSS.weights
        if not np.all(np.isfinite(weight)) or np.any(weight <= 0):
            raise ValueError("C3D8 requires finite positive Jacobians at every IP; "
                             "a non-positive one means the element has inverted")
        if integration == "selective_reduced":
            # The volumetric row from the element centroid, matching the F_bar
            # above. The mean-dilatation and centroid forms agree exactly for a
            # parallelepiped and differ only for a distorted element.
            centroid, _ = _k.b_matrix_spatial(element, np.zeros(3))
            trace_c = centroid[:3].sum(axis=0)
            volume_row = operators[:, :3, :].sum(axis=1)
            operators[:, :3, :] += ((trace_c - volume_row) / 3)[:, None, :]
        matrices.append(operators)
        weights.append(weight)
    return np.array(matrices), np.array(weights)


def geometric_stiffness(coordinates, displacement, stress, weights=None):
    """The initial-stress term, (ne, 24, 24).

    It comes from differentiating the current-configuration integral, not from
    the strain measure, so B_bar does not enter it: the unmodified spatial
    gradients are the right ones here.
    """
    coordinates = np.asarray(coordinates, dtype=float)
    displacement = np.asarray(displacement, dtype=float)
    stress = np.asarray(stress, dtype=float)
    current = coordinates + displacement
    points, base = _k.ABAQUS_C3D8_GAUSS.points, _k.ABAQUS_C3D8_GAUSS.weights
    identity = np.eye(3)
    out = np.zeros(coordinates.shape[:1] + (24, 24))
    for e, element in enumerate(current):
        for q, point in enumerate(points):
            _, detJ, dNdx = _k._b_from_coords(element, point)
            w = detJ * base[q] if weights is None else weights[e, q]
            G = dNdx @ _k.voigt_to_tensor(stress[e, q]) @ dNdx.T
            out[e] += np.kron(G, identity) * w
    return out
