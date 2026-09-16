"""Small-strain J2 radial return, written once over any number type.

The bundled J2 UMAT (``UMATs/UMATs/generic_ps/j2_props.f`` in the UMAT
repository, PROPS = E, NU, SIGY0, H and STATEV(1) = EQPLAS) is the one
material in the frozen collection for which this repository carries its own
model. It exists here for one purpose: to give a residual-derivative check an
ANALYTIC side for dR/dp and dR/dq that is not the same code as the finite
difference it is compared with.

* :func:`update` is the radial return, written so that every operation works
  on a ``float`` or on a :class:`~residual_core.algebra.dual1.Dual1`. Seeding
  one parameter with a unit dual part gives ``d sigma / d p`` by forward-mode
  differentiation -- exact to rounding, no step size.
* :func:`dsigma_dq` is ``d sigma_(n+1) / d EQPLAS_n`` in closed form.

The finite-difference side of every check that uses this module re-evaluates
the stress with ``umat_oti.validation.j2_reference.integrate_increment``, a
separate implementation in the other repository. Agreement between the two is
therefore evidence about the assembly and about both models, not one model
checked against itself.

What this does NOT establish: that the TRANSFORMED UMAT's own parameter or
state sensitivities are right. A case whose verification did not claim them
has none to check, and a report built on this module says which side each
number came from.
"""
from __future__ import annotations

import math
from typing import Sequence

import numpy as np

from residual_core.algebra.dual1 import Dual1

__all__ = ["PARAMETERS", "update", "dsigma_dq", "dsigma_dp", "sqrt"]

#: PROPS order in the bundled source.
PARAMETERS = ("E", "NU", "SIGY0", "H")


def sqrt(x):
    if isinstance(x, Dual1):
        root = math.sqrt(x.real)
        return Dual1(root, x.imag / (2.0 * root) if root else 0.0)
    return math.sqrt(x)


def _real(x) -> float:
    return x.real if isinstance(x, Dual1) else float(x)


def _stiffness(E, nu):
    lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    mu = E / (2.0 * (1.0 + nu))
    C = [[0.0] * 6 for _ in range(6)]
    for i in range(3):
        for j in range(3):
            C[i][j] = lam
        C[i][i] = lam + 2.0 * mu
    for i in range(3, 6):
        C[i][i] = mu
    return C, mu


def update(props: Sequence, stress_n: Sequence[float], eqplas_n,
           dstrain: Sequence[float]) -> tuple:
    """(stress_(n+1), eqplas_(n+1), yielded) for one increment.

    Engineering shear in ``dstrain``; Voigt order 11, 22, 33, 12, 13, 23.
    """
    E, nu, sigy0, H = props
    C, mu = _stiffness(E, nu)
    trial = [stress_n[i] + sum(C[i][j] * dstrain[j] for j in range(6))
             for i in range(6)]
    pressure = (trial[0] + trial[1] + trial[2]) / 3.0
    dev = [trial[0] - pressure, trial[1] - pressure, trial[2] - pressure,
           trial[3], trial[4], trial[5]]
    q_trial = sqrt(1.5 * (dev[0] * dev[0] + dev[1] * dev[1] + dev[2] * dev[2])
                   + 3.0 * (dev[3] * dev[3] + dev[4] * dev[4]
                            + dev[5] * dev[5]))
    phi = q_trial - (sigy0 + H * eqplas_n)
    if _real(phi) <= 0.0 or _real(q_trial) <= 0.0:
        return trial, eqplas_n, False
    dgamma = phi / (3.0 * mu + H)
    factor = 1.0 - 3.0 * mu * dgamma / q_trial
    stress = [dev[0] * factor + pressure, dev[1] * factor + pressure,
              dev[2] * factor + pressure, dev[3] * factor, dev[4] * factor,
              dev[5] * factor]
    return stress, eqplas_n + dgamma, True


def dsigma_dq(props: Sequence[float], stress_n: Sequence[float],
              eqplas_n: float, dstrain: Sequence[float]) -> np.ndarray:
    """``d sigma_(n+1) / d EQPLAS_n``, closed form; zero in the elastic branch.

    With ``dgamma = (q_trial - SIGY0 - H q_n) / (3 mu + H)`` only ``dgamma``
    depends on ``q_n``: ``d dgamma/d q_n = -H/(3 mu + H)`` and so
    ``d sigma/d q_n = dev_trial * 3 mu H / (q_trial (3 mu + H))``.
    """
    E, nu, sigy0, H = (float(v) for v in props)
    C, mu = _stiffness(E, nu)
    trial = np.array([stress_n[i] + sum(C[i][j] * dstrain[j] for j in range(6))
                      for i in range(6)])
    pressure = trial[:3].mean()
    dev = trial.copy()
    dev[:3] -= pressure
    q_trial = math.sqrt(1.5 * float(np.sum(dev[:3] ** 2))
                        + 3.0 * float(np.sum(dev[3:] ** 2)))
    if q_trial - (sigy0 + H * eqplas_n) <= 0.0 or q_trial <= 0.0:
        return np.zeros(6)
    return dev * (3.0 * mu * H / (q_trial * (3.0 * mu + H)))


def dsigma_dp(props: Sequence[float], stress_n: Sequence[float],
              eqplas_n: float, dstrain: Sequence[float]) -> np.ndarray:
    """``d sigma_(n+1) / d PROPS``, shape (6, 4), by forward-mode dual numbers.

    The incoming state is held fixed: this is the derivative of ONE increment's
    update with respect to the material constants, which is what the residual
    at that increment depends on when the start-of-increment state is given.
    """
    out = np.zeros((6, len(PARAMETERS)))
    for k in range(len(PARAMETERS)):
        seeded = [Dual1(float(v), 1.0 if i == k else 0.0)
                  for i, v in enumerate(props)]
        stress, _q, _y = update(seeded, [float(s) for s in stress_n],
                                float(eqplas_n), [float(e) for e in dstrain])
        out[:, k] = [s.imag if isinstance(s, Dual1) else 0.0 for s in stress]
    return out
