"""Neutral Voigt / elastic-moduli helpers shared across layers.

Lives in ``core`` so both ``formulations`` and ``materials`` can depend on it
DOWNWARD without a materials->formulations layering inversion. Pure numpy, no
element or material knowledge.

Voigt order is Abaqus order (11,22,33,12,13,23) with ENGINEERING shear, matching
CONTRACT.md section 1.
"""

from __future__ import annotations

import numpy as np


def isotropic_D(E, nu):
    """Isotropic linear-elastic 6x6 tangent in Abaqus Voigt order with ENGINEERING
    shear (tau = mu*gamma), so it pairs with the engineering-shear B matrix.

    This is the canonical copy. ``formulations/c3d8_kernel.py`` keeps a private
    duplicate purely so it remains runnable as a bare script; keep the two in sync.
    """
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    D = np.zeros((6, 6), dtype=float)
    D[:3, :3] = lam
    for i in range(3):
        D[i, i] += 2 * mu
    for i in range(3, 6):
        D[i, i] = mu
    return D
