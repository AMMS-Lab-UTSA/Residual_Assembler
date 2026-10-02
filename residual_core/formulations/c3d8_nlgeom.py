"""C3D8 under the Abaqus/Standard NLGEOM contract for a UMAT with state.

What a UMAT is handed in an Abaqus/Standard geometrically nonlinear step, and
what it must return, is fixed by the UMAT interface (Abaqus User Subroutines
Guide, "UMAT", argument list) and by the strain/rotation increment the solver
forms (Abaqus Theory Guide, "Strain measures" / "Rate of deformation and strain
increment" and "Stress rates"). This module reproduces that contract for the
first-order brick and gives the EXACT linearisation of the resulting discrete
residual. The derivation, assumption ledger and sources are in
``corpus_campaign/batches/B2/noether/DESIGN.md``; the short form:

Per integration point q, increment n (configuration n-1 converged, n trial)::

    F_q      = I + grad_X u                       (unmodified, reference grads)
    Fbar_q   = (J_c / J_q)^(1/3) F_q              selectively reduced C3D8:
                                                  volume change from the centroid
    M        = (Fbar0 + Fbar1) / 2                mid-increment configuration
    L        = (Fbar1 - Fbar0) M^-1               = d(Delta u) / d x_{n-1/2}
    DSTRAN   = sym L        (Voigt, engineering shear)
    DROT     = (I - W/2)^-1 (I + W/2),  W = skew L          (Hughes-Winget)
    STRESS_in = DROT sigma_{n-1} DROT^T           Cauchy, pre-rotated by the host
    STRAN_in  = DROT eps_{n-1} DROT^T             (tensor shear inside the rotation)
    DFGRD0 = Fbar0, DFGRD1 = Fbar1

The UMAT returns the Cauchy stress sigma_n. ``integration='full'`` drops the
centroid modification; ``rotation='none'`` passes DROT = I and does not
pre-rotate (the B1 convention, kept to measure what it changes).

Element residual (current configuration, B-bar consistent with Fbar)::

    f_(a,k) = sum_q w_q detJ0_q J_q [ (dev sigma_q . g_qa)_k + p_q g_ca,k ]

with g_qa = dN_a/dx at q, g_ca the same at the centroid, p = tr(sigma)/3. This
is ``sum_q w detJ Bbar_q^T sigma_q`` for the B-bar of
``residual_core.replay.kinematics.spatial_operators``; with ``full`` it is
``sum_q w detJ B_q^T sigma_q`` (``c3d8_kernel.element_internal_force_finite_strain``).

Exact Jacobian::

    dR/du = [d f / d u]_(sigma fixed)  +  sum_q w detJ Bbar_q^T (dsigma_q/dFbar1_q : dFbar1_q/du)

The first term is analytic here (``dforce_fixed_stress``: the variation of the
spatial gradients, the current volume and the centroid gradients). The second
needs the material's TOTAL derivative with respect to Fbar1, through every
input that depends on it (DFGRD1, DSTRAN, DROT, the pre-rotated STRESS and
STRAN): the caller obtains it by directional seeding (``increment_direction``
and ``rotate_direction`` give the seeds for a direction dFbar1).
"""
from __future__ import annotations

from typing import Dict

import numpy as np

from . import c3d8_kernel as _k

__all__ = [
    "GAUSS", "voigt_eng", "stress_voigt_to_tensor", "stress_tensor_to_voigt",
    "strain_voigt_to_tensor", "strain_tensor_to_voigt", "increment", "increment_direction",
    "rotate", "rotate_direction", "C3D8Nlgeom", "InvertedElement", "INTEGRATIONS", "ROTATIONS",
]

GAUSS = _k.ABAQUS_C3D8_GAUSS
INTEGRATIONS = ("full", "mean_dilatation", "centroid")
#: Abaqus calls the fully integrated first-order brick "selectively reduced";
#: what it does is mean dilatation (B2 probe), so the name maps there.
ALIASES = {"selective_reduced": "mean_dilatation"}
ROTATIONS = ("abaqus", "none")
_I = np.eye(3)
_I6 = np.eye(6)


class InvertedElement(ValueError):
    """det F <= 0 at a point of the element: the trial configuration is inadmissible."""


# ---------------------------------------------------------------------------
# Voigt conventions: (11, 22, 33, 12, 13, 23); stress tensor shear, strain
# engineering shear (Abaqus UMAT convention for both).
# ---------------------------------------------------------------------------
def voigt_eng(t: np.ndarray) -> np.ndarray:
    """(..., 3, 3) -> (..., 6): diagonal, then t_ij + t_ji (engineering shear)."""
    t = np.asarray(t)
    return np.stack([t[..., 0, 0], t[..., 1, 1], t[..., 2, 2],
                     t[..., 0, 1] + t[..., 1, 0], t[..., 0, 2] + t[..., 2, 0],
                     t[..., 1, 2] + t[..., 2, 1]], axis=-1)


def _sym_from_voigt(v: np.ndarray, shear_factor: float) -> np.ndarray:
    v = np.asarray(v)
    out = np.empty(v.shape[:-1] + (3, 3), dtype=v.dtype)
    out[..., 0, 0], out[..., 1, 1], out[..., 2, 2] = v[..., 0], v[..., 1], v[..., 2]
    out[..., 0, 1] = out[..., 1, 0] = shear_factor * v[..., 3]
    out[..., 0, 2] = out[..., 2, 0] = shear_factor * v[..., 4]
    out[..., 1, 2] = out[..., 2, 1] = shear_factor * v[..., 5]
    return out


def stress_voigt_to_tensor(s):
    return _sym_from_voigt(s, 1.0)


def strain_voigt_to_tensor(e):
    return _sym_from_voigt(e, 0.5)


def stress_tensor_to_voigt(T):
    T = np.asarray(T)
    return np.stack([T[..., 0, 0], T[..., 1, 1], T[..., 2, 2],
                     0.5 * (T[..., 0, 1] + T[..., 1, 0]), 0.5 * (T[..., 0, 2] + T[..., 2, 0]),
                     0.5 * (T[..., 1, 2] + T[..., 2, 1])], axis=-1)


def strain_tensor_to_voigt(T):
    return voigt_eng(T)


def _T(a):
    return np.swapaxes(a, -1, -2)


# ---------------------------------------------------------------------------
# Increment kinematics (point-wise; any leading shape)
# ---------------------------------------------------------------------------
def increment(F0: np.ndarray, F1: np.ndarray) -> Dict[str, np.ndarray]:
    """DSTRAN and the Hughes-Winget rotation increment from F0 -> F1.

    Returns ``L`` (incremental displacement gradient w.r.t. the mid-increment
    configuration), ``Minv``, ``dstran`` (Voigt, engineering shear), ``W``,
    ``Ainv`` = (I - W/2)^-1 and ``drot``.
    """
    F0 = np.asarray(F0, dtype=float)
    F1 = np.asarray(F1, dtype=float)
    Minv = np.linalg.inv(0.5 * (F0 + F1))
    L = (F1 - F0) @ Minv
    W = 0.5 * (L - _T(L))
    Ainv = np.linalg.inv(_I - 0.5 * W)
    drot = Ainv @ (_I + 0.5 * W)
    return {"L": L, "Minv": Minv, "dstran": voigt_eng(L), "W": W, "Ainv": Ainv, "drot": drot}


def increment_direction(kin: Dict[str, np.ndarray], dF0, dF1) -> Dict[str, np.ndarray]:
    """Directional derivatives of ``increment`` for dF0, dF1 of shape (..., m, 3, 3).

    ``kin`` arrays have shape (..., 3, 3); returned ``dstran`` (..., m, 6) and
    ``ddrot`` (..., m, 3, 3). Exact (no truncation):

        dL    = (dF1 - dF0) M^-1 - L dM M^-1,   dM = (dF0 + dF1)/2
        ddrot = 1/2 A^-1 dW (drot + I),         A = I - W/2
    """
    Minv = kin["Minv"][..., None, :, :]
    L = kin["L"][..., None, :, :]
    dF0 = np.asarray(dF0, dtype=float)
    dF1 = np.asarray(dF1, dtype=float)
    dL = (dF1 - dF0) @ Minv - L @ (0.5 * (dF0 + dF1)) @ Minv
    dW = 0.5 * (dL - _T(dL))
    ddrot = 0.5 * kin["Ainv"][..., None, :, :] @ dW @ (kin["drot"][..., None, :, :] + _I)
    return {"dL": dL, "dstran": voigt_eng(dL), "dW": dW, "ddrot": ddrot}


def rotate(drot, T):
    """drot T drot^T."""
    return drot @ T @ _T(drot)


def rotate_direction(drot, ddrot, T, dT):
    """Derivative of drot T drot^T: ddrot T drot^T + drot dT drot^T + drot T ddrot^T.

    ``drot``/``T`` (..., 3, 3); ``ddrot``/``dT`` (..., m, 3, 3) (either may be None = 0).
    """
    R = drot[..., None, :, :]
    S = T[..., None, :, :]
    out = 0.0
    if ddrot is not None:
        out = out + ddrot @ S @ _T(R) + R @ S @ _T(ddrot)
    if dT is not None:
        out = out + R @ dT @ _T(R)
    return out


# ---------------------------------------------------------------------------
# The element
# ---------------------------------------------------------------------------
class C3D8Nlgeom:
    """A batch of C3D8 elements, current-configuration residual and its exact
    linearisation.

    ``integration``:
      * ``'full'`` -- 2x2x2 Gauss, unmodified F and B;
      * ``'centroid'`` -- the volume change from the element CENTROID,
        Fbar = (J_c/J_q)^(1/3) F_q, volumetric B row from the centroid gradient:
        the operators commit 30a5ab4 put in ``residual_core.replay.kinematics``;
        NOT what Abaqus does on a distorted element -- a labelled comparison;
      * ``'mean_dilatation'`` (alias ``'selective_reduced'``; the default and
        what Abaqus C3D8 does) -- the volume change averaged over the element,
        Jbar = sum_q w0_q J_q / sum_q w0_q (= v / V), volumetric B row
        bbar_a = (1/v) sum_q w J_q g_qa (Nagtegaal-Parks-Rice / Hughes B-bar).
        Equal to the centroid form only where J is affine in the natural
        coordinates (small strain on a parallelepiped); a trilinear motion
        of even a cubic element separates them at second order.

    Arrays: ``Ue`` (ne, 8, 3) element nodal displacements; stresses (ne, 8, 6)
    Cauchy Voigt; parameter/direction axis ``m`` last for stresses
    (ne, 8, m, 6) and second for displacement directions (ne, m, 8, 3).
    """

    def __init__(self, Xe: np.ndarray, integration: str = "mean_dilatation"):
        integration = ALIASES.get(integration, integration)
        if integration not in INTEGRATIONS:
            raise ValueError("integration must be one of %s (or %s)"
                             % (INTEGRATIONS, sorted(ALIASES)))
        self.integration = integration
        self.sri = integration != "full"
        self.mean = integration == "mean_dilatation"
        self.Xe = np.asarray(Xe, dtype=float)
        self.ne = self.Xe.shape[0]
        points = np.vstack([GAUSS.points, np.zeros((1, 3))])       # 8 IPs + centroid
        dNdxi = np.array([_k.shape_grad_natural(p) for p in points])  # (9, 8, 3)
        J0 = np.einsum("eai,qaj->eqij", self.Xe, dNdxi)              # dX/dxi
        self.detJ0 = np.linalg.det(J0)                                # (ne, 9)
        if np.any(self.detJ0 <= 0) or not np.all(np.isfinite(self.detJ0)):
            raise ValueError("C3D8 requires positive reference Jacobians")
        self.G = np.einsum("qak,eqkj->eqaj", dNdxi, np.linalg.inv(J0))  # dN/dX (ne, 9, 8, 3)
        self.w0 = self.detJ0[:, :8] * GAUSS.weights                   # (ne, 8)
        self.N = np.array([_k.shape_functions(p) for p in GAUSS.points])  # (8, 8)

    # ------------------------------------------------------------ kinematics
    def F(self, Ue) -> np.ndarray:
        """(ne, 9, 3, 3) unmodified F at the 8 IPs and the centroid."""
        return _I + np.einsum("eai,eqaj->eqij", np.asarray(Ue, dtype=float), self.G)

    def _Jvol(self, J):
        """(ne, 1) the element's volume-change measure."""
        if self.mean:
            return (np.sum(self.w0 * J[:, :8], axis=1) / np.sum(self.w0, axis=1))[:, None]
        return J[:, 8:9]

    def Fbar(self, Ue) -> np.ndarray:
        """(ne, 8, 3, 3) the deformation gradient handed to the UMAT."""
        F = self.F(Ue)
        if not self.sri:
            return F[:, :8].copy()
        J = np.linalg.det(F)
        if np.any(J <= 0):
            raise InvertedElement("det F <= 0: the element has inverted")
        alpha = (self._Jvol(J) / J[:, :8]) ** (1.0 / 3.0)
        return alpha[..., None, None] * F[:, :8]

    def dFbar(self, Ue, dUe) -> np.ndarray:
        """(ne, 8, m, 3, 3) directional derivative of Fbar for dUe (ne, m, 8, 3)."""
        F = self.F(Ue)
        dF = np.einsum("emai,eqaj->eqmij", np.asarray(dUe, dtype=float), self.G)
        if not self.sri:
            return dF[:, :8]
        J = np.linalg.det(F)
        tr = np.einsum("eqij,eqmji->eqm", np.linalg.inv(F), dF)         # tr(F^-1 dF) = dJ/J
        Jv = self._Jvol(J)
        if self.mean:
            dlogJv = (np.einsum("eq,eqm->em", self.w0 * J[:, :8], tr[:, :8])
                      / np.sum(self.w0 * J[:, :8], axis=1)[:, None])[:, None, :]
        else:
            dlogJv = tr[:, 8:9, :]
        alpha = (Jv / J[:, :8]) ** (1.0 / 3.0)                           # (ne, 8)
        corr = (dlogJv - tr[:, :8, :]) / 3.0                             # (ne, 8, m)
        return alpha[..., None, None, None] * (dF[:, :8] + corr[..., None, None]
                                               * F[:, :8, None])

    def increment(self, Ue0, Ue1) -> Dict[str, np.ndarray]:
        """What Abaqus/Standard hands the UMAT for the increment u0 -> u1, per IP
        (ne, 8, ...). Established against Abaqus 2021 on a distorted C3D8 with an
        inhomogeneous rotating field (B2 probe, all to <= 7e-15):

        * DFGRD0/DFGRD1 = Fbar (``Fbar``; mean dilatation for C3D8);
        * L_q = (F1 - F0) M_q^-1 with the UNMODIFIED F, M = (F0 + F1)/2;
        * DSTRAN = sym L_q with its trace replaced by the element average of
          tr L weighted by the mid-increment volume w0 det M (centroid: tr L_c;
          full integration: unchanged);
        * DROT = (I - W/2)^-1 (I + W/2), W = skew L_q (Hughes-Winget).
        """
        F0, F1 = self.F(Ue0), self.F(Ue1)
        M = 0.5 * (F0 + F1)
        Minv = np.linalg.inv(M)
        L = (F1 - F0) @ Minv                                             # (ne, 9, 3, 3)
        trL = np.einsum("eqii->eq", L)
        if not self.sri:
            trv = trL[:, :8]
            wm = None
        elif self.mean:
            wm = self.w0 * np.linalg.det(M[:, :8])
            trv = np.broadcast_to((np.sum(wm * trL[:, :8], axis=1) / np.sum(wm, axis=1))[:, None],
                                  (self.ne, 8))
        else:
            wm = None
            trv = np.broadcast_to(trL[:, 8:9], (self.ne, 8))
        Lq = L[:, :8]
        Lbar = Lq + ((trv - trL[:, :8]) / 3.0)[..., None, None] * _I
        W = 0.5 * (Lq - _T(Lq))
        Ainv = np.linalg.inv(_I - 0.5 * W)
        return {"Fbar0": self.Fbar(Ue0), "Fbar1": self.Fbar(Ue1), "dstran": voigt_eng(Lbar),
                "drot": Ainv @ (_I + 0.5 * W), "Ainv": Ainv, "L": L, "Minv": Minv, "trL": trL,
                "trv": np.array(trv), "wm": wm, "M": M, "Ue0": np.asarray(Ue0, dtype=float),
                "Ue1": np.asarray(Ue1, dtype=float)}

    def increment_direction(self, inc: Dict[str, np.ndarray], dUe0=None, dUe1=None):
        """Exact directional derivatives of ``increment`` for dUe0/dUe1
        (ne, m, 8, 3) (None = 0): dFbar0, dFbar1 (ne, 8, m, 3, 3), dstran
        (ne, 8, m, 6), ddrot (ne, 8, m, 3, 3)."""
        m = (dUe0 if dUe0 is not None else dUe1).shape[1]
        zero = np.zeros((self.ne, m, 8, 3))
        dUe0 = zero if dUe0 is None else np.asarray(dUe0, dtype=float)
        dUe1 = zero if dUe1 is None else np.asarray(dUe1, dtype=float)
        dF0 = np.einsum("emai,eqaj->eqmij", dUe0, self.G)
        dF1 = np.einsum("emai,eqaj->eqmij", dUe1, self.G)
        Minv = inc["Minv"][:, :, None]
        L = inc["L"][:, :, None]
        dM = 0.5 * (dF0 + dF1)
        dL = (dF1 - dF0) @ Minv - L @ dM @ Minv                          # (ne, 9, m, 3, 3)
        dtrL = np.einsum("eqmii->eqm", dL)
        if not self.sri:
            dtrv = dtrL[:, :8]
        elif self.mean:
            wm = inc["wm"]
            dwm = wm[:, :, None] * np.einsum("eqij,eqmji->eqm", Minv[:, :8, 0], dM[:, :8])
            sw = np.sum(wm, axis=1)[:, None]
            trv = inc["trv"][:, :1]
            dtrv = ((np.einsum("eqm,eq->em", dwm, inc["trL"][:, :8])
                     + np.einsum("eq,eqm->em", wm, dtrL[:, :8])
                     - trv * np.sum(dwm, axis=1)) / sw)[:, None, :]
            dtrv = np.broadcast_to(dtrv, (self.ne, 8, m))
        else:
            dtrv = np.broadcast_to(dtrL[:, 8:9], (self.ne, 8, m))
        dLq = dL[:, :8]
        dLbar = dLq + ((dtrv - dtrL[:, :8]) / 3.0)[..., None, None] * _I
        dW = 0.5 * (dLq - _T(dLq))
        ddrot = 0.5 * inc["Ainv"][:, :, None] @ dW @ (inc["drot"][:, :, None] + _I)
        return {"dFbar0": self.dFbar(inc["Ue0"], dUe0), "dFbar1": self.dFbar(inc["Ue1"], dUe1),
                "dstran": voigt_eng(dLbar), "ddrot": ddrot}

    def geometry(self, Ue):
        """Spatial gradients g (ne, 9, 8, 3) from the unmodified F, integration
        weights of the residual wF (ne, 8), and the volumetric gradient gv
        (ne, 8, 3).

        full: wF = w0 J_q. centroid (the 30a5ab4 comparison form):
        wF = w0 J_q, gv = g at the centroid. mean dilatation (Abaqus C3D8):
        wF = w0 Jbar -- the volume of the Fbar material, whose Kirchhoff stress
        is Jbar sigma -- and gv = sum_q w0 J_q g_q / sum_q w0 J_q (B2 probe:
        Abaqus reactions reproduced to 4.7e-8, the ODB's float32 resolution).
        """
        F = self.F(Ue)
        detF = np.linalg.det(F)
        if np.any(detF <= 0):
            raise InvertedElement("det F <= 0: the element has inverted")
        g = np.einsum("eqak,eqkj->eqaj", self.G, np.linalg.inv(F))
        wJ = self.w0 * detF[:, :8]
        if self.mean:
            gv = np.einsum("eq,eqak->eak", wJ, g[:, :8]) / np.sum(wJ, axis=1)[:, None, None]
            wF = self.w0 * (np.sum(wJ, axis=1) / np.sum(self.w0, axis=1))[:, None]
        else:
            gv = g[:, 8]
            wF = wJ
        return g, wF, gv, wJ

    def coordinates(self, Ue=None) -> np.ndarray:
        """(ne, 8, 3) integration-point positions (reference, or current with Ue)."""
        X = self.Xe if Ue is None else self.Xe + np.asarray(Ue)
        return np.einsum("qa,eai->eqi", self.N, X)

    # ------------------------------------------------------------- residual
    def _force(self, geo, S):
        """S (ne, 8, m, 3, 3) -> (ne, 24, m)."""
        g, wJ, gv, _ = geo
        if self.sri:
            p = np.einsum("eqmii->eqm", S) / 3.0
            dev = S - p[..., None, None] * _I
            f = (np.einsum("eqmkj,eqaj,eq->eakm", dev, g[:, :8], wJ)
                 + np.einsum("eqm,eak,eq->eakm", p, gv, wJ))
        else:
            f = np.einsum("eqmkj,eqaj,eq->eakm", S, g[:, :8], wJ)
        return f.reshape(self.ne, 24, -1)

    def _stress_basis(self):
        """(ne, 8, 48, 6) the 48 real unit stresses (point q, component s) as a
        stack, column q*6 + s."""
        basis = np.zeros((self.ne, 8, 48, 6))
        for q in range(8):
            basis[:, q, 6 * q:6 * q + 6, :] = _I6
        return basis

    @staticmethod
    def _combine(T, sigma):
        """T (ne, 24, [m,] 48) real, sigma (ne, 8, [m,] 6) of any number type:
        sum_b T[..., b] sigma_b with real coefficients, so a Dual1/OTI stress
        keeps its derivative (no cast of the stress to float anywhere)."""
        T = np.asarray(T).astype(object)
        sigma = np.asarray(sigma, dtype=object)
        if sigma.ndim == 3:                                   # (ne, 8, 6)
            return np.array([np.dot(T[e], sigma[e].reshape(48)) for e in range(len(T))],
                            dtype=object)
        if sigma.ndim == 4:                                   # stacked (ne, 8, m, 6)
            s = np.transpose(sigma, (0, 2, 1, 3)).reshape(len(T), sigma.shape[2], 48)
            return np.array([np.dot(T[e], s[e].T) for e in range(len(T))], dtype=object)
        raise ValueError("unsupported stress shape %s" % (sigma.shape,))

    def force(self, Ue, sigma) -> np.ndarray:
        """Element internal force (ne, 24) for Cauchy stress (ne, 8, 6), or
        (ne, 24, m) for a stack (ne, 8, m, 6) -- linear in sigma at fixed u, so
        the same call assembles dR/dp from dsigma/dp.

        Geometry is real. A non-real stress (object array of Dual1 / OTILib
        numbers, from seeded material constants) is assembled by linearity --
        the real force of each unit stress, combined with the stress values --
        so its derivatives are kept rather than cast away."""
        sigma = np.asarray(sigma)
        if sigma.dtype == object:
            return self._combine(self.force(Ue, self._stress_basis()), sigma)
        sigma = np.asarray(sigma, dtype=float)
        stacked = sigma.ndim == 4
        S = stress_voigt_to_tensor(sigma if stacked else sigma[:, :, None, :])
        f = self._force(self.geometry(Ue), S)
        return f if stacked else f[:, :, 0]

    def dforce_fixed_stress(self, Ue, sigma, dUe) -> np.ndarray:
        """(ne, 24, m): derivative of ``force`` along dUe (ne, m, 8, 3), sigma held.
        Linear in sigma: a non-real (object) stress is combined from the real
        unit-stress results, as in ``force``."""
        sigma = np.asarray(sigma)
        if sigma.dtype == object:
            basis = self._stress_basis()
            T = np.stack([self.dforce_fixed_stress(Ue, basis[:, :, b], dUe) for b in range(48)],
                         axis=-1)                                        # (ne, 24, m, 48)
            out = np.empty(T.shape[:3], dtype=object)
            for e in range(self.ne):
                out[e] = np.dot(T[e].astype(object), sigma[e].reshape(48))
            return out
        S = stress_voigt_to_tensor(np.asarray(sigma, dtype=float))     # (ne, 8, 3, 3)
        g, wF, gv, wJ = self.geometry(Ue)
        L = np.einsum("emai,eqaj->eqmij", np.asarray(dUe, dtype=float), g)   # grad_x du
        trL = np.einsum("eqmii->eqm", L[:, :8])
        dwJ = wJ[:, :, None] * trL
        dg = -np.einsum("eqak,eqmkj->eqmaj", g, L)                      # (ne, 9, m, 8, 3)
        if not self.sri:
            f = (np.einsum("eqkj,eqaj,eqm->eakm", S, g[:, :8], dwJ)
                 + np.einsum("eqkj,eqmaj,eq->eakm", S, dg[:, :8], wJ))
            return f.reshape(self.ne, 24, -1)
        if self.mean:
            v = np.sum(wJ, axis=1)
            dlogJbar = np.sum(dwJ, axis=1) / v[:, None]                  # (ne, m)
            dwF = wF[:, :, None] * dlogJbar[:, None, :]
            dgv = ((np.einsum("eqm,eqak->emak", dwJ, g[:, :8])
                    + np.einsum("eq,eqmak->emak", wJ, dg[:, :8])
                    - gv[:, None] * np.sum(dwJ, axis=1)[:, :, None, None]) / v[:, None, None, None])
        else:
            dwF = dwJ
            dgv = dg[:, 8]
        p = np.einsum("eqii->eq", S) / 3.0
        dev = S - p[..., None, None] * _I
        f = (np.einsum("eqkj,eqaj,eqm->eakm", dev, g[:, :8], dwF)
             + np.einsum("eq,eak,eqm->eakm", p, gv, dwF)
             + np.einsum("eqkj,eqmaj,eq->eakm", dev, dg[:, :8], wF)
             + np.einsum("eq,emak,eq->eakm", p, dgv, wF))
        return f.reshape(self.ne, 24, -1)

    @staticmethod
    def unit_directions(ne: int) -> np.ndarray:
        """(ne, 24, 8, 3): the 24 element DOF unit directions."""
        return np.broadcast_to(np.eye(24).reshape(24, 8, 3), (ne, 24, 8, 3))

    def tangent_exact(self, Ue, sigma, dsig_dFbar1) -> np.ndarray:
        """(ne, 24, 24) exact element Jacobian; dsig_dFbar1 (ne, 8, 6, 3, 3) is
        the material's total d sigma / d Fbar1 (all inputs that move with it)."""
        E = self.unit_directions(self.ne)
        dFb = self.dFbar(Ue, E)                                          # (ne, 8, 24, 3, 3)
        dsig = np.einsum("eqsij,eqmij->eqms", np.asarray(dsig_dFbar1), dFb)
        return self.dforce_fixed_stress(Ue, sigma, E) + self.force(Ue, dsig)

    # ---------------------------------------------- the Abaqus-style tangent
    def bbar(self, Ue):
        """(ne, 8, 6, 24) B (full) or B-bar (SRI) in the current configuration, and wJ."""
        g, wJ, gv, _ = self.geometry(Ue)
        B = np.zeros((self.ne, 8, 6, 8, 3))
        gq = g[:, :8]
        B[:, :, 0, :, 0] = gq[..., 0]
        B[:, :, 1, :, 1] = gq[..., 1]
        B[:, :, 2, :, 2] = gq[..., 2]
        B[:, :, 3, :, 0] = gq[..., 1]
        B[:, :, 3, :, 1] = gq[..., 0]
        B[:, :, 4, :, 0] = gq[..., 2]
        B[:, :, 4, :, 2] = gq[..., 0]
        B[:, :, 5, :, 1] = gq[..., 2]
        B[:, :, 5, :, 2] = gq[..., 1]
        if self.sri:
            corr = (gv[:, None] - gq) / 3.0                         # (ne, 8, 8, 3)
            B[:, :, :3] += corr[:, :, None]
        return B.reshape(self.ne, 8, 6, 24), wJ, g

    def tangent_abaqus(self, Ue, sigma, ddsdde) -> np.ndarray:
        """sum w detJ Bbar^T c Bbar + K_geo with c from DDSDDE read as the Jaumann
        rate of Kirchhoff stress over J (UMAT finite-strain convention) and
        K_geo the initial-stress term (delta_ik g_a.sigma.g_b) -- the matrix an
        Abaqus-style element builds from DDSDDE. NOT the exact Jacobian of
        ``force`` in general; scored against a difference of R."""
        B, wJ, g = self.bbar(Ue)
        sigma = np.asarray(sigma, dtype=float)
        D = np.asarray(ddsdde, dtype=float)
        c = np.empty_like(D)
        for e in range(self.ne):
            for q in range(8):
                c[e, q] = _k.kirchhoff_jaumann_to_spatial(D[e, q], sigma[e, q])
        K = np.einsum("eqsi,eqst,eqtj,eq->eij", B, c, B, wJ)
        S = stress_voigt_to_tensor(sigma)
        Gs = np.einsum("eqak,eqkl,eqbl,eq->eab", g[:, :8], S, g[:, :8], wJ)   # (ne, 8, 8)
        K += np.einsum("eab,ij->eaibj", Gs, _I).reshape(self.ne, 24, 24)
        return K
