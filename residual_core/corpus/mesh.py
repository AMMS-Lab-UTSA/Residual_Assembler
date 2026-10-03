"""Small C3D8 meshes and displacement-controlled boundary-value problems.

Everything a problem prescribes is a displacement proportional to one load
factor ``lambda``; no external force acts, so at a converged state the
internal force vector R(u) vanishes on the free DOFs and equals the reaction
on the prescribed ones. Prescribed values do not depend on the material
parameters, so ``du_c/dp = 0`` exactly.

Node numbering of a brick: ``i + (nx+1) * (j + (ny+1) * k)``; element nodes in
Abaqus C3D8 order (bottom face counter-clockwise, then top face). DOF of node
``a`` and direction ``d`` is ``3 a + d``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

__all__ = ["Mesh", "Problem", "brick", "uniaxial", "clamped_tension", "clamped_shear",
           "patch_test", "DEFAULT_PATH", "rotation_matrix", "with_superposed_rotation"]

#: load, partial unload, reload past the first peak: (lambda at segment end, increments)
DEFAULT_PATH: Tuple[Tuple[float, int], ...] = ((1.0, 4), (0.4, 3), (1.25, 3))


@dataclass
class Mesh:
    coords: np.ndarray                 # (nn, 3)
    conn: np.ndarray                   # (ne, 8) zero-based
    size: Tuple[float, float, float]
    divisions: Tuple[int, int, int]

    @property
    def ndof(self) -> int:
        return 3 * len(self.coords)

    def nodes_where(self, axis: int, value: float, tol: float = 1e-9) -> List[int]:
        return [int(a) for a in np.flatnonzero(np.abs(self.coords[:, axis] - value) < tol)]

    def node_at(self, point: Sequence[float], tol: float = 1e-9) -> int:
        hit = np.flatnonzero(np.linalg.norm(self.coords - np.asarray(point), axis=1) < tol)
        if hit.size != 1:
            raise ValueError("no unique node at %s" % (list(point),))
        return int(hit[0])


def brick(divisions=(1, 1, 1), size=(1.0, 1.0, 1.0), distort: float = 0.0) -> Mesh:
    """A structured brick; ``distort`` moves interior nodes (patch-test geometry)."""
    nx, ny, nz = divisions
    lx, ly, lz = size
    coords = []
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                coords.append((lx * i / nx, ly * j / ny, lz * k / nz))
    coords = np.array(coords, dtype=float)

    def nid(i, j, k):
        return i + (nx + 1) * (j + (ny + 1) * k)

    conn = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                conn.append([nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
                             nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1),
                             nid(i, j + 1, k + 1)])
    if distort:
        h = np.array([lx / nx, ly / ny, lz / nz])
        pattern = np.array([[0.31, -0.17, 0.23], [-0.29, 0.21, -0.13], [0.19, 0.27, -0.25]])
        for a, x in enumerate(coords):
            interior = all(1e-9 < x[d] < (lx, ly, lz)[d] - 1e-9 for d in range(3))
            if interior:
                coords[a] = x + distort * h * pattern[a % 3]
    return Mesh(coords=coords, conn=np.array(conn, dtype=int), size=(lx, ly, lz),
                divisions=(nx, ny, nz))


@dataclass
class Problem:
    name: str
    mesh: Mesh
    prescribed: Dict[int, float]       # dof -> value at lambda = 1
    path: Tuple[Tuple[float, int], ...] = DEFAULT_PATH
    segment_period: float = 1.0
    qois: List[dict] = field(default_factory=list)
    description: str = ""
    #: optional proper orthogonal Q per increment (superposed rigid rotation,
    #: x = Q (X + V)); None = no rotation. See engine.Frame.
    rotations: Optional[List[np.ndarray]] = None
    #: optional (8, 3) corner coordinates (Abaqus C3D8 node order) of one
    #: element of the AUTHOR's mesh: COORDS handed to the routine are the probe
    #: points mapped into that element (``material_coordinates``); the
    #: mechanics stay on ``mesh``. None = COORDS are the probe's own positions.
    material_element: Optional[np.ndarray] = None

    @property
    def constrained(self) -> np.ndarray:
        return np.array(sorted(self.prescribed), dtype=int)

    @property
    def free(self) -> np.ndarray:
        return np.setdiff1d(np.arange(self.mesh.ndof), self.constrained)

    def amplitudes(self) -> np.ndarray:
        return np.array([self.prescribed[d] for d in self.constrained], dtype=float)

    def schedule(self) -> List[dict]:
        """Increments: lambda at end, step time, total time at start, dtime."""
        out, start, total = [], 0.0, 0.0
        for segment, (end, count) in enumerate(self.path):
            dtime = self.segment_period / count
            for m in range(1, count + 1):
                out.append({"lambda_start": start + (end - start) * (m - 1) / count,
                            "lambda": start + (end - start) * m / count,
                            "dtime": dtime, "time_start": total,
                            "segment": segment, "segment_end": m == count})
                total += dtime
            start = end
        return out


def _face(mesh: Mesh, axis: int, at_max: bool) -> List[int]:
    return mesh.nodes_where(axis, mesh.size[axis] if at_max else 0.0)


def uniaxial(mesh: Mesh, strain: float, path=DEFAULT_PATH) -> Problem:
    """Symmetry planes x=0, y=0, z=0; x = Lx face pulled by ``strain * Lx``.

    Homogeneous for a homogeneous material. QoIs: the x-reaction on the pulled
    face, and the y-displacement of the far corner (lateral contraction).
    """
    lx, ly, lz = mesh.size
    prescribed = {}
    for d in range(3):
        for a in _face(mesh, d, False):
            prescribed[3 * a + d] = 0.0
    pulled = _face(mesh, 0, True)
    for a in pulled:
        prescribed[3 * a + 0] = strain * lx
    corner = mesh.node_at((lx, ly, lz))
    return Problem(
        name="uniaxial_%dx%dx%d" % mesh.divisions, mesh=mesh, prescribed=prescribed, path=path,
        qois=[{"name": "reaction_x_pulled_face", "kind": "reaction", "dofs": [3 * a for a in pulled]},
              {"name": "uy_far_corner", "kind": "displacement", "dof": 3 * corner + 1}],
        description="symmetry planes x=0,y=0,z=0; face x=Lx displaced by lambda*%g*Lx" % strain)


def clamped_tension(mesh: Mesh, strain: float, path=DEFAULT_PATH) -> Problem:
    """z=0 clamped; z=Lz face displaced by ``strain*Lz`` in z with ux=uy=0.

    Inhomogeneous: the clamps restrain lateral contraction at both ends.
    """
    lx, ly, lz = mesh.size
    prescribed = {}
    for a in _face(mesh, 2, False):
        for d in range(3):
            prescribed[3 * a + d] = 0.0
    top = _face(mesh, 2, True)
    for a in top:
        prescribed[3 * a + 0] = 0.0
        prescribed[3 * a + 1] = 0.0
        prescribed[3 * a + 2] = strain * lz
    side = mesh.node_at((lx, ly / 2.0, lz / 2.0)) if mesh.divisions[1] % 2 == 0 and \
        mesh.divisions[2] % 2 == 0 else mesh.node_at((lx, ly, lz))
    return Problem(
        name="clamped_tension_%dx%dx%d" % mesh.divisions, mesh=mesh, prescribed=prescribed,
        path=path,
        qois=[{"name": "reaction_z_top", "kind": "reaction", "dofs": [3 * a + 2 for a in top]},
              {"name": "ux_side_midheight", "kind": "displacement", "dof": 3 * side + 0}],
        description="z=0 clamped; z=Lz face uz=lambda*%g*Lz, ux=uy=0" % strain)


def clamped_shear(mesh: Mesh, shear: float, path=DEFAULT_PATH) -> Problem:
    """z=0 clamped; z=Lz face displaced by ``shear*Lz`` in x, uy=uz=0."""
    lx, ly, lz = mesh.size
    prescribed = {}
    for a in _face(mesh, 2, False):
        for d in range(3):
            prescribed[3 * a + d] = 0.0
    top = _face(mesh, 2, True)
    for a in top:
        prescribed[3 * a + 0] = shear * lz
        prescribed[3 * a + 1] = 0.0
        prescribed[3 * a + 2] = 0.0
    side = mesh.node_at((lx, ly / 2.0, lz / 2.0)) if mesh.divisions[1] % 2 == 0 and \
        mesh.divisions[2] % 2 == 0 else mesh.node_at((lx, ly, lz))
    return Problem(
        name="clamped_shear_%dx%dx%d" % mesh.divisions, mesh=mesh, prescribed=prescribed,
        path=path,
        qois=[{"name": "reaction_x_top", "kind": "reaction", "dofs": [3 * a + 0 for a in top]},
              {"name": "uz_side_midheight", "kind": "displacement", "dof": 3 * side + 2}],
        description="z=0 clamped; z=Lz face ux=lambda*%g*Lz, uy=uz=0" % shear)


def patch_test(mesh: Mesh, gradient: np.ndarray, increments: int = 2) -> Problem:
    """Every boundary node carries u = gradient . X; interior nodes are free.

    The exact solution for ANY homogeneous material is the linear field, so the
    interior nodes must land on ``gradient . X`` and every point must see the
    same strain -- on a distorted mesh as well.
    """
    gradient = np.asarray(gradient, dtype=float)
    prescribed = {}
    for a, x in enumerate(mesh.coords):
        on_boundary = any(abs(x[d]) < 1e-9 or abs(x[d] - mesh.size[d]) < 1e-9 for d in range(3))
        if on_boundary:
            u = gradient @ x
            for d in range(3):
                prescribed[3 * a + d] = float(u[d])
    return Problem(name="patch_%dx%dx%d" % mesh.divisions, mesh=mesh, prescribed=prescribed,
                   path=((1.0, increments),), qois=[],
                   description="linear displacement field on every boundary node")


def rotation_matrix(axis, angle: float) -> np.ndarray:
    """Proper orthogonal Q = exp(angle [n]x), Rodrigues."""
    n = np.asarray(axis, dtype=float)
    n = n / np.linalg.norm(n)
    K = np.array([[0.0, -n[2], n[1]], [n[2], 0.0, -n[0]], [-n[1], n[0], 0.0]])
    return np.eye(3) + np.sin(angle) * K + (1.0 - np.cos(angle)) * K @ K


def with_superposed_rotation(problem: Problem, axis=(1.0, 2.0, 3.0), angle: float = np.pi / 3,
                             after_segment: int = 1, increments: int = 3):
    """Two problems for an objectivity test: the original path with a HOLD
    segment (load factor constant, ``increments`` increments) inserted after
    segment ``after_segment``, unrotated; and the same with a rigid rotation
    growing from I to Q(axis, angle) over the hold and kept to the end.

    Rotation-only increments: the Hughes-Winget rotation of a pure rotation
    increment is that rotation exactly and DSTRAN is exactly zero, so an
    objective material and update must give the same frame displacements and
    reactions, and sigma' = Q sigma Q^T, to round-off.
    """
    path = list(problem.path)
    hold_value = path[after_segment - 1][0]
    path.insert(after_segment, (hold_value, increments))
    count = sum(c for _, c in path)
    before = sum(c for _, c in path[:after_segment])
    rotations = []
    for index in range(count):
        k = min(max(index - before + 1, 0), increments)
        rotations.append(rotation_matrix(axis, angle * k / increments))
    import dataclasses
    plain = dataclasses.replace(problem, path=tuple(path), rotations=None,
                                name=problem.name + "_hold")
    rotated = dataclasses.replace(problem, path=tuple(path), rotations=rotations,
                                  name=problem.name + "_rotated")
    return plain, rotated


#: Abaqus C3D8 corner order in the parent cube [-1, 1]^3
_CORNERS = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                     [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], dtype=float)


def material_coordinates(points: np.ndarray, mesh: Mesh, element: np.ndarray) -> np.ndarray:
    """Map probe positions ``points`` (reference configuration of ``mesh``)
    into the hexahedron ``element`` (8 corners, Abaqus order): the probe box
    [lo, lo + size] is the parent cube, the map the trilinear isoparametric one.
    """
    points = np.asarray(points, dtype=float)
    element = np.asarray(element, dtype=float).reshape(8, 3)
    lo = mesh.coords.min(axis=0)
    size = np.asarray(mesh.coords.max(axis=0) - lo, dtype=float)
    xi = 2.0 * (points - lo) / size - 1.0                              # (..., 3) in [-1, 1]
    N = np.prod(1.0 + xi[..., None, :] * _CORNERS, axis=-1) / 8.0      # (..., 8)
    return N @ element


def author_element(case) -> Optional[np.ndarray]:
    """(8, 3) corners of the author's element recorded by the verification run
    (``manifest.node_coordinates``: [label, x, y, z] rows), or None."""
    rows = (getattr(case, "extra", None) or {}).get("node_coordinates") or []
    if len(rows) != 8:
        return None
    return np.array([[float(v) for v in row[-3:]] for row in rows])
