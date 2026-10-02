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
# The operators below delegate to residual_core.formulations.c3d8_nlgeom,
# which holds the one implementation. What Abaqus 2021 does for a fully
# integrated C3D8 ("selectively reduced") was established against an Abaqus
# run on a DISTORTED element with an inhomogeneous rotating field (B2 probe,
# corpus_campaign/batches/B2/noether/abaqus_probe): MEAN DILATATION --
#
#     Jbar  = sum_q w0_q J_q / sum_q w0_q          (= v / V)
#     F_bar = (Jbar / J_q)^(1/3) F_q               DFGRD0/1 to 4e-16
#     f_a   = sum_q w0_q Jbar [dev sigma_q g_qa + p_q gbar_a],
#     gbar_a = sum_q w0_q J_q g_qa / sum_q w0_q J_q    reactions to 4.7e-8
#
# Commit 30a5ab4 took the volume change from the element CENTROID instead.
# The two coincide where J is affine in the natural coordinates; on the
# probe's element the centroid form is 2.4e-4 off in DFGRD and 2.1e-2 in the
# reactions, and it fails the patch test on a distorted mesh (2e-5) which
# mean dilatation passes to round-off. It stays available as
# ``volume="centroid"`` / ``integration="centroid"``, labelled as the
# comparison, never as the default.
#
# Two corrections are needed and BOTH have to be applied: the kinematics
# (F_bar to the material) and the assembly (B_bar in the virtual work).
# Keeping plain B while the material responds to F_bar drops the volumetric
# coupling out of dsigma/du and makes the tangent non-symmetric besides.
# ---------------------------------------------------------------------------

_VOLUMES = {"mean": "mean_dilatation", "centroid": "centroid"}


def _element(coordinates, integration):
    from ..formulations.c3d8_nlgeom import C3D8Nlgeom
    return C3D8Nlgeom(np.asarray(coordinates, dtype=float), integration)


def deformation_gradients(coordinates, displacement, fbar=True, volume="mean"):
    """F at every integration point, (ne, 8, 3, 3): what the UMAT is handed
    as DFGRD1. With ``fbar`` the volume change is the element's
    (``volume="mean"``: mean dilatation, Abaqus C3D8; ``"centroid"``: the
    30a5ab4 comparison form), ``F_bar = (Jvol / J_q)^(1/3) F_q``."""
    if volume not in _VOLUMES:
        raise ValueError("volume must be one of %s" % sorted(_VOLUMES))
    integration = _VOLUMES[volume] if fbar else "full"
    return _element(coordinates, integration).Fbar(np.asarray(displacement, dtype=float))


def spatial_operators(coordinates, displacement, integration="selective_reduced"):
    """B (or B_bar) in the CURRENT configuration and the integration weights,
    (ne, 8, 6, 24) and (ne, 8), such that the internal force is
    ``sum_q w_q B_q^T sigma_q``.

    ``integration``: ``"selective_reduced"`` / ``"mean_dilatation"`` (Abaqus
    C3D8: volumetric row the current-volume average, weight w0 Jbar),
    ``"centroid"`` (the 30a5ab4 comparison form) or ``"full"``.
    """
    from ..formulations.c3d8_nlgeom import ALIASES, INTEGRATIONS
    if ALIASES.get(integration, integration) not in INTEGRATIONS:
        raise ValueError(f"unsupported C3D8 integration: {integration}")
    element = _element(coordinates, integration)
    operators, weights, _ = element.bbar(np.asarray(displacement, dtype=float))
    return operators, weights


def internal_force(coordinates, displacement, stress, integration="selective_reduced"):
    """(ne, 24) element internal force in the current configuration."""
    return _element(coordinates, integration).force(np.asarray(displacement, dtype=float),
                                                    np.asarray(stress, dtype=float))


def force_tangent_fixed_stress(coordinates, displacement, stress,
                               integration="selective_reduced"):
    """(ne, 24, 24): the EXACT derivative of ``internal_force`` with the Cauchy
    stress held -- spatial gradients, current volume and the volumetric
    average all move with u. The exact K is this plus
    sum_q w Bbar^T (dsigma/du)."""
    element = _element(coordinates, integration)
    ue = np.asarray(displacement, dtype=float)
    return element.dforce_fixed_stress(ue, np.asarray(stress, dtype=float),
                                       element.unit_directions(element.ne))


def geometric_stiffness(coordinates, displacement, stress, weights=None):
    """The conventional initial-stress term, (ne, 24, 24),
    ``sum_q w (g_a . sigma . g_b) delta_ik`` with the unmodified spatial
    gradients -- what an Abaqus-style tangent adds to Bbar^T c Bbar. It is NOT
    the exact fixed-stress derivative of the B-bar force (see
    ``force_tangent_fixed_stress``); the two differ by the objective-rate terms
    and the variation of the volumetric average.
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
