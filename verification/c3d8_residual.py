"""Small generic C3D8 residual-method driver for the presentation claims.

This is a bounded, dense, small-strain driver used by the crystal-plasticity
claims (slides 28-32). It reuses the replay kinematics of the Residual Assembler
(``residual_core.replay.kinematics``: Abaqus C3D8 2x2x2 Gauss points with the
mean-dilatation B-bar that the fully integrated Abaqus C3D8 uses) and the
compiled-provider ABI (``PathMaterial``). It adds only what those modules do not
have: nonzero prescribed displacements and any provider, not just m3_j2.

Primal: incremental Newton on the free degrees of freedom, every material call
starting from the committed state of the previous increment.

Sensitivity (the residual method): at the converged state of increment n,

    R_n(u_n; p) = sum_e sum_q B^T sigma_n dV - F_ext = 0,
    sigma_n = sigma(sigma_{n-1}, state_{n-1}, B (u_n - u_{n-1}); p),

so, with parameter-independent loads and prescribed displacements,

    K_ff du_n/dp = -[ sum B^T ( S_p - D B du_{n-1}/dp ) dV ]_f ,
    dsigma_n/dp  =  S_p + D B (du_n/dp - du_{n-1}/dp),

where D = DDSDDE (dsigma/dDSTRAN at fixed incoming history) and S_p is the
provider's DSIGMA_DP with the incoming total derivatives carried
(``DSIGMA_DP_IN``/``DSTATEV_DP_IN``). The state derivative is carried the same
way, except that the provider ABI returns no d(state)/d(DSTRAN); its product with
B d(u_n - u_{n-1})/dp is therefore not added. That term is exactly zero whenever
the displacement sensitivity does not change between increments (it is measured
and reported as ``max_strain_sensitivity_increment``); the claims that use this
driver report that bound instead of assuming it.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence

import numpy as np

from residual_core.replay.kinematics import (
    batch_element_force, batch_element_tangent, batch_operators, batch_strain,
)


# --------------------------------------------------------------------------- #
# mesh
# --------------------------------------------------------------------------- #
def cube_mesh(n: int, length: float = 1.0):
    """n x n x n C3D8 mesh of [0, L]^3 in Abaqus node order. Returns (coords, conn)."""
    npc = n + 1
    ident = lambda i, j, k: i + npc * (j + npc * k)  # noqa: E731
    coords = np.array([[length * i / n, length * j / n, length * k / n]
                       for k in range(npc) for j in range(npc) for i in range(npc)], float)
    conn = []
    for k in range(n):
        for j in range(n):
            for i in range(n):
                conn.append([ident(i, j, k), ident(i + 1, j, k), ident(i + 1, j + 1, k),
                             ident(i, j + 1, k), ident(i, j, k + 1), ident(i + 1, j, k + 1),
                             ident(i + 1, j + 1, k + 1), ident(i, j + 1, k + 1)])
    return coords, np.array(conn, int)


@dataclass
class Problem:
    """A C3D8 mesh, prescribed displacements per increment, and time steps."""

    coords: np.ndarray
    conn: np.ndarray
    prescribed: np.ndarray                 # global dof indices (node*3 + comp)
    values: np.ndarray                     # (n_increments, len(prescribed)) total values
    dts: Sequence[float]
    integration: str = "selective_reduced"
    ops: np.ndarray = field(init=False, repr=False)
    weights: np.ndarray = field(init=False, repr=False)
    edofs: np.ndarray = field(init=False, repr=False)

    def __post_init__(self):
        self.ops, self.weights = batch_operators([self.coords[c] for c in self.conn], self.integration)
        self.edofs = np.array([[3 * node + comp for node in element for comp in range(3)]
                               for element in self.conn], int)
        self.ndof = 3 * len(self.coords)
        mask = np.ones(self.ndof, bool)
        mask[self.prescribed] = False
        self.free = np.flatnonzero(mask)
        self.volume = float(self.weights.sum())


def simple_shear_problem(n_elements: int, gamma_max: float, n_increments: int, dt: float,
                         all_boundary: bool = False) -> Problem:
    """Simple shear u1 = gamma * y on the unit cube (engineering gamma_12).

    ``all_boundary=False``: bottom face y=0 fixed, top face y=L given u1 = gamma*L
    and u2 = 0; u3 of the top face and every dof of any other node is free.
    ``all_boundary=True``: the affine field is prescribed on every boundary node
    and only interior nodes are free (the controlled-deformation RVE test).
    """
    coords, conn = cube_mesh(n_elements)
    presc, target = [], []
    tol = 1e-12
    for node, (x, y, z) in enumerate(coords):
        on_boundary = any(abs(c) < tol or abs(c - 1.0) < tol for c in (x, y, z))
        if all_boundary and on_boundary:
            presc += [3 * node, 3 * node + 1, 3 * node + 2]
            target += [y, 0.0, 0.0]
        elif not all_boundary and abs(y) < tol:
            presc += [3 * node, 3 * node + 1, 3 * node + 2]
            target += [0.0, 0.0, 0.0]
        elif not all_boundary and abs(y - 1.0) < tol:
            presc += [3 * node, 3 * node + 1]
            target += [1.0, 0.0]
    target = np.array(target)
    factors = gamma_max * np.arange(1, n_increments + 1) / n_increments
    values = factors[:, None] * target[None, :]
    return Problem(coords, conn, np.array(presc, int), values, [dt] * n_increments)


# --------------------------------------------------------------------------- #
# solver
# --------------------------------------------------------------------------- #
@dataclass
class Solution:
    u: List[np.ndarray]
    stress: List[np.ndarray]               # (ne, 8, 6) per increment
    state: List[np.ndarray]
    iterations: List[int]
    residual: List[float]
    du_dp: Optional[List[np.ndarray]] = None
    dstress_dp: Optional[List[np.ndarray]] = None      # (ne, 8, 6, P)
    max_strain_sensitivity_increment: float = 0.0
    max_du_dp: float = 0.0


def solve(material, props: Sequence[float], problem: Problem, *, mode: str = "regular",
          sensitivities: bool = False, tol: float = 1e-12, max_iter: int = 50) -> Solution:
    """Incremental Newton solve; ``mode`` = 'regular' (ORIGINAL UMAT) or 'oti'."""
    if sensitivities and mode != "oti":
        raise ValueError("sensitivities need the OTI provider (mode='oti')")
    ne = len(problem.conn)
    nt, ns, npar = material.ntens, material.nstatev, material.nparam
    ops, weights, edofs = problem.ops, problem.weights, problem.edofs
    stress = np.zeros((ne, 8, nt)); state = np.zeros((ne, 8, ns))
    dstress = np.zeros((ne, 8, nt, npar)); dstate = np.zeros((ne, 8, ns, npar))
    u_prev = np.zeros(problem.ndof); du_prev = np.zeros((problem.ndof, npar))
    free = problem.free
    out = Solution([], [], [], [], [], [] if sensitivities else None, [] if sensitivities else None)
    zero_ds, zero_dst = np.zeros((nt, npar)), np.zeros((ns, npar))

    def evaluate(u):
        strain = batch_strain(ops, (u - u_prev)[edofs])
        new_s = np.empty_like(stress); new_st = np.empty_like(state)
        tangent = np.empty((ne, 8, nt, nt))
        part_s = np.empty_like(dstress) if mode == "oti" else None
        part_st = np.empty_like(dstate) if mode == "oti" else None
        for e in range(ne):
            for q in range(8):
                if mode == "regular":
                    s, st, dd = material._reg_step(props, stress[e, q], state[e, q], strain[e, q], dt)
                else:
                    s, st, dd, ds, dst = material._oti_step(
                        props, stress[e, q], state[e, q],
                        dstress[e, q] if sensitivities else zero_ds,
                        dstate[e, q] if sensitivities else zero_dst, strain[e, q], dt)
                    part_s[e, q] = ds; part_st[e, q] = dst
                new_s[e, q], new_st[e, q], tangent[e, q] = s, st, dd
        return strain, new_s, new_st, tangent, part_s, part_st

    def assemble(tangent, field_values):
        k_local = batch_element_tangent(ops, tangent, weights)
        f_local = batch_element_force(ops, field_values, weights)
        stiffness = np.zeros((problem.ndof, problem.ndof))
        force = np.zeros((problem.ndof,) + field_values.shape[3:])
        for e in range(ne):
            stiffness[np.ix_(edofs[e], edofs[e])] += k_local[e]
            np.add.at(force, edofs[e], f_local[e])
        return stiffness, force

    for n, dt in enumerate(problem.dts):
        u = u_prev.copy()
        u[problem.prescribed] = problem.values[n]
        for iteration in range(1, max_iter + 1):
            strain, new_s, new_st, tangent, part_s, part_st = evaluate(u)
            stiffness, force = assemble(tangent, new_s)
            scale = max(float(np.max(np.abs(force))), 1.0)
            error = float(np.max(np.abs(force[free]))) / scale if free.size else 0.0
            if error < tol:
                break
            u[free] += np.linalg.solve(stiffness[np.ix_(free, free)], -force[free])
        else:
            raise RuntimeError(f"increment {n + 1}: Newton did not converge (scaled R={error:.3e})")
        out.iterations.append(iteration); out.residual.append(error)
        if sensitivities:
            d_bu_prev = np.einsum("eqai,eim->eqam", ops, du_prev[edofs])
            history = part_s - np.einsum("eqab,eqbm->eqam", tangent, d_bu_prev)
            rhs = assemble(tangent, history)[1]
            du = np.zeros((problem.ndof, npar))
            if free.size:
                du[free] = np.linalg.solve(stiffness[np.ix_(free, free)], -rhs[free])
            d_strain = np.einsum("eqai,eim->eqam", ops, (du - du_prev)[edofs])
            dstress = part_s + np.einsum("eqab,eqbm->eqam", tangent, d_strain)
            dstate = part_st
            out.max_strain_sensitivity_increment = max(out.max_strain_sensitivity_increment,
                                                       float(np.max(np.abs(d_strain))))
            out.max_du_dp = max(out.max_du_dp, float(np.max(np.abs(du))))
            out.du_dp.append(du.copy()); out.dstress_dp.append(dstress.copy())
            du_prev = du
        stress, state, u_prev = new_s, new_st, u
        out.u.append(u.copy()); out.stress.append(new_s.copy()); out.state.append(new_st.copy())
    return out


def volume_average(problem: Problem, values: np.ndarray) -> np.ndarray:
    """Volume average over all integration points of a per-IP field (ne, 8, ...)."""
    return np.einsum("eq,eq...->...", problem.weights, values) / problem.volume


def mises_field(stress: np.ndarray) -> np.ndarray:
    s = stress
    return np.sqrt(0.5 * ((s[..., 0] - s[..., 1]) ** 2 + (s[..., 1] - s[..., 2]) ** 2
                          + (s[..., 2] - s[..., 0]) ** 2)
                   + 3.0 * (s[..., 3] ** 2 + s[..., 4] ** 2 + s[..., 5] ** 2))


def dmises_field(stress: np.ndarray, dstress: np.ndarray) -> np.ndarray:
    """d sigma_vM / dp per IP from (.., 6) stress and (.., 6, P) derivative."""
    s = stress
    vm = mises_field(s)
    hydro = (s[..., 0] + s[..., 1] + s[..., 2]) / 3.0
    g = np.empty_like(s)
    safe = np.where(vm > 1e-30, vm, 1.0)
    for i in range(3):
        g[..., i] = 1.5 * (s[..., i] - hydro) / safe
    for i in range(3, 6):
        g[..., i] = 3.0 * s[..., i] / safe
    g[vm <= 1e-30] = 0.0
    return np.einsum("...i,...ip->...p", g, dstress)
