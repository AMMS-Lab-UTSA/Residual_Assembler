"""2-node beam backend — heterogeneous-DOF proof (translations + rotations).

A classic 3D small-strain Euler-Bernoulli frame element. Unlike the solid and
truss backends it carries **six DOFs per node** — three translations and three
rotations — which is exactly what proves the framework's DOF machinery is not
displacement-only:

    dof_types = (UX, UY, UZ, RX, RY, RZ)

Local 12x12 stiffness = axial + St-Venant torsion + Euler-Bernoulli bending about
both local principal axes, rotated to global with the element triad. The element
is linear, so ``F_internal,e = K_e u_e`` and the tangent is constant (Level-4 FD
check is exact).

Section/material data are read from the element ``properties`` (a
MaterialBinding.section dict or a plain dict):
    E   Young's modulus            (required)
    A   cross-section area         (required)
    Iz  2nd moment about local z   (required; ``I`` used for both if only I given)
    Iy  2nd moment about local y   (default = Iz)
    G   shear modulus              (default from nu=0.3)
    J   torsion constant           (default = Iy + Iz)
    x2  orientation reference vec  (optional; default auto)

This is a proof backend, not a production beam/shell element.
"""

from __future__ import annotations

import numpy as np

from .base import Formulation
from .truss2 import _prop


def _local_stiffness(E, G, A, Iy, Iz, J, L):
    """Standard 3D Euler-Bernoulli frame local stiffness (12x12).

    Local DOF order per node: [ux, uy, uz, rx, ry, rz]."""
    K = np.zeros((12, 12))
    EA_L = E * A / L
    GJ_L = G * J / L
    az = 12.0 * E * Iz / L**3
    bz = 6.0 * E * Iz / L**2
    cz = 4.0 * E * Iz / L
    dz = 2.0 * E * Iz / L
    ay = 12.0 * E * Iy / L**3
    by = 6.0 * E * Iy / L**2
    cy = 4.0 * E * Iy / L
    dy = 2.0 * E * Iy / L

    # axial (ux): dofs 0, 6
    K[0, 0] = EA_L;  K[0, 6] = -EA_L
    K[6, 6] = EA_L
    # torsion (rx): dofs 3, 9
    K[3, 3] = GJ_L;  K[3, 9] = -GJ_L
    K[9, 9] = GJ_L
    # bending in x-y plane -> uy (1,7) & rz (5,11), uses Iz
    K[1, 1] = az;  K[1, 5] = bz;   K[1, 7] = -az;  K[1, 11] = bz
    K[5, 5] = cz;  K[5, 7] = -bz;  K[5, 11] = dz
    K[7, 7] = az;  K[7, 11] = -bz
    K[11, 11] = cz
    # bending in x-z plane -> uz (2,8) & ry (4,10), uses Iy (note sign flip)
    K[2, 2] = ay;  K[2, 4] = -by;  K[2, 8] = -ay;  K[2, 10] = -by
    K[4, 4] = cy;  K[4, 8] = by;   K[4, 10] = dy
    K[8, 8] = ay;  K[8, 10] = by
    K[10, 10] = cy
    # symmetrize (only upper triangle filled above)
    K = K + K.T - np.diag(np.diag(K))
    return K


def _triad(coords, ref=None):
    """Local axes (rows) expressed in global coords for the element."""
    coords = np.asarray(coords, dtype=float)
    d = coords[1] - coords[0]
    L = float(np.linalg.norm(d))
    if L <= 0.0:
        raise ValueError("beam2 element has zero length")
    x = d / L
    if ref is None:
        ref = np.array([0.0, 0.0, 1.0])
        if abs(float(x @ ref)) > 0.99:            # nearly parallel -> switch
            ref = np.array([0.0, 1.0, 0.0])
    else:
        ref = np.asarray(ref, dtype=float)
    y = np.cross(ref, x)
    y /= np.linalg.norm(y)
    z = np.cross(x, y)
    R = np.vstack([x, y, z])                        # (3,3) rows = local axes
    return R, L


class Beam2(Formulation):
    name = "beam2"
    element_types = ("B31", "B33", "B31_LIKE", "BEAM2", "FRAME3D")
    dof_types = ("UX", "UY", "UZ", "RX", "RY", "RZ")
    sign_convention = "residual"
    verification_levels = (0, 1, 4)
    supported_modes = ("formulation",)
    required_inputs = ("coords", "connectivity", "dofs", "section: E, A, I")
    optional_inputs = ("G", "J", "Iy", "orientation")
    limitations = ("small-strain Euler-Bernoulli initially",
                   "2 nodes, 6 DOF/node", "no shear deformation (no Timoshenko)")
    verification_tests = ("zero-field", "cantilever-tip", "FD-tangent")
    notes = "Proof backend: exercises rotational DOFs (heterogeneous DOF sets)."
    required_inputs_by_mode = {
        "formulation": ("coords", "connectivity", "dofs", "section: E, A, I")}
    optional_inputs_by_mode = {"formulation": ("G", "J", "Iy", "orientation")}
    material_interface_needed = False
    state_requirements = "none"
    tangent_support = "analytic"
    verification_status = "implemented/simple"

    def n_nodes(self, element_type):
        return 2

    def global_stiffness(self, coords, props):
        E = _prop(props, "E", ("youngs", "young"))
        A = _prop(props, "A", ("area",))
        Iz = _prop(props, "Iz", ("I", "Izz"))
        if E is None or A is None or Iz is None:
            raise ValueError("beam2 needs section properties E, A, Iz (or I).")
        Iy = _prop(props, "Iy", ("Iyy",), default=Iz)
        nu = _prop(props, "nu", (), default=0.3)
        G = _prop(props, "G", ("shear_modulus",), default=E / (2.0 * (1.0 + nu)))
        J = _prop(props, "J", ("torsion",), default=Iy + Iz)
        ref = getattr(props, "section", {}).get("orientation") if hasattr(props, "section") \
            else (props.get("orientation") if isinstance(props, dict) else None)

        R, L = _triad(coords, ref)
        K_local = _local_stiffness(E, G, A, Iy, Iz, J, L)
        # 12x12 transform: block-diag of R over the 4 nodal triads
        T = np.zeros((12, 12))
        for b in range(4):
            T[3*b:3*b+3, 3*b:3*b+3] = R
        K_global = T.T @ K_local @ T
        return K_global, L

    def eval_element(self, element_id, element_type, coords, dofs,
                     solution_state, material_state, properties,
                     time, dtime, fields, options):
        options = options or {}
        compute_tangent = options.get("compute_tangent", True)
        u_e = np.asarray(dofs, dtype=float).reshape(12)
        K, L = self.global_stiffness(coords, properties)
        r = K @ u_e
        tangent = K if compute_tangent else None
        return r, tangent, material_state, {"formulation": self.name, "length": L}
