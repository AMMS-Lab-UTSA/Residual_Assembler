"""
c3d8_residual.py  --  C3D8 element internal-force / tangent residual assembler.

Pure numpy. NO dependency on the parser (operates on plain arrays), so it is
independently testable against analytic / divergence-theorem checks that need
no Abaqus.

This module implements the public API pinned in `residual_core/CONTRACT.md`
section 4, obeying the conventions in sections 1, 2 and 6.

=====================================================================
CONVENTIONS (binding -- see CONTRACT sections 1 & 2; do not guess silently)
=====================================================================
* Voigt order (stress & strain):  (11, 22, 33, 12, 13, 23)  -- Abaqus order.
* Strain uses ENGINEERING shear:   eps_voigt = (e11,e22,e33, 2e12, 2e13, 2e23).
  Internal virtual work is then  delta_eps_voigt . sigma_voigt  with the true
  (not doubled) shear stresses. The B matrix is built to produce engineering
  shear strains (rows for shear = dN/dxj on dof i + dN/dxi on dof j).
* Element DOF ordering: node-major, 24 dofs per C3D8:
  [u1x,u1y,u1z, u2x,u2y,u2z, ... , u8x,u8y,u8z], node k = k-th entry of the
  element connectivity list.
* Global DOF ordering: nodes sorted ascending by id -> 0-based index i;
  node id -> dofs [3i, 3i+1, 3i+2].

CONFIGURATION / STRESS MEASURE per function (CONTRACT section 0):
* `b_matrix_reference`  : shape-fn derivatives w.r.t. REFERENCE coords Xe.
  Used for the small-strain (F = I) reference-config path -- patch tests only.
* `b_matrix_spatial`    : shape-fn derivatives w.r.t. CURRENT coords xe = Xe+Ue.
  Used for the finite-strain (nlgeom=YES) path that matches Abaqus:
      f_int,e = integral_v  B_spatial(x)^T . sigma_cauchy  dv
  with dv = detJ_current * w. `sigma_ip` there are CAUCHY (true) stresses in
  Abaqus Voigt order -- exactly what the ngrilli UMAT stores (kmat.f:1023).
"""

import numpy as np

# ---------------------------------------------------------------------------
# Section 2 geometry constants
# ---------------------------------------------------------------------------

# Node natural coordinates in Abaqus C3D8 connectivity order (node 1..8).
ABAQUS_C3D8_NODES = np.array([
    [-1.0, -1.0, -1.0],   # n1
    [+1.0, -1.0, -1.0],   # n2
    [+1.0, +1.0, -1.0],   # n3
    [-1.0, +1.0, -1.0],   # n4
    [-1.0, -1.0, +1.0],   # n5
    [+1.0, -1.0, +1.0],   # n6
    [+1.0, +1.0, +1.0],   # n7
    [-1.0, +1.0, +1.0],   # n8
], dtype=float)


def _abaqus_c3d8_gauss():
    """Build the 8 full-integration Gauss points in ABAQUS C3D8 integration-
    point (SDV / S output) order: fastest index is xi1, then xi2, then xi3.

        IP1 (-g,-g,-g)  IP2 (+g,-g,-g)  IP3 (-g,+g,-g)  IP4 (+g,+g,-g)
        IP5 (-g,-g,+g)  IP6 (+g,-g,+g)  IP7 (-g,+g,+g)  IP8 (+g,+g,+g)
    with g = 1/sqrt(3), all weights = 1.

    Note this "raster" order (xi1 fastest) is deliberately DIFFERENT from the
    C3D8 *node* order (which walks the face CCW): IP3/IP4 are NOT in node order.
    Getting that distinction right is a classic pitfall.

    !! ORDERING CAVEAT (verified by audit 2026-07-10) !!
    This per-IP ORDERING is NOT validated by any offline test in this repo: every
    offline check uses a spatially-uniform stress, which is *provably blind* to IP
    ordering (all IPs share the same stress). The ordering here is the standard
    Abaqus/CalculiX lexicographic convention (corroborated against CalculiX
    `gauss3d2` / `e_c3d.f`) -- high confidence, but it can only be *proven* against
    a real ODB. When the first `fields.json` from an Abaqus run exists, verify the
    ordering by comparing a single-C3D8 job with a known spatially-varying stress
    (see tests/cp_c3d8_umat/stress_driven_residual/README.md). The INTERNAL pairing
    "sigma_ip[k] <-> points[k]" is verified; only its match to Abaqus's export
    index is pending an ODB.
    """
    g = 1.0 / np.sqrt(3.0)
    pts = []
    for x3 in (-g, +g):
        for x2 in (-g, +g):
            for x1 in (-g, +g):
                pts.append((x1, x2, x3))
    points = np.array(pts, dtype=float)
    weights = np.ones(8, dtype=float)
    return points, weights


class _GaussRule:
    """Small container exposing `.points` (8,3) and `.weights` (8,).

    Also iterable/indexable as (points, weights) so callers can unpack it or
    access the arrays by attribute -- points and weights are exposed clearly.
    """

    __slots__ = ("points", "weights")

    def __init__(self, points, weights):
        self.points = points
        self.weights = weights

    def __iter__(self):
        yield self.points
        yield self.weights

    def __getitem__(self, i):
        return (self.points, self.weights)[i]

    def __len__(self):
        return 2


_pts, _wts = _abaqus_c3d8_gauss()
# Named constant: (points:(8,3), weights:(8,)) in Abaqus IP order.
ABAQUS_C3D8_GAUSS = _GaussRule(_pts, _wts)


# ---------------------------------------------------------------------------
# Shape functions
# ---------------------------------------------------------------------------

def shape_functions(xi):
    """Trilinear C3D8 shape functions N_a(xi), a = 1..8  ->  (8,).

    N_a = 1/8 (1 + xi1*xi1a)(1 + xi2*xi2a)(1 + xi3*xi3a).
    """
    xi = np.asarray(xi, dtype=float)
    n = ABAQUS_C3D8_NODES  # (8,3) of +-1
    return 0.125 * (1.0 + n[:, 0] * xi[0]) * (1.0 + n[:, 1] * xi[1]) * (1.0 + n[:, 2] * xi[2])


def shape_grad_natural(xi):
    """Gradients dN_a/dxi_j w.r.t. natural coords  ->  (8,3).

    dN_a/dxi1 = 1/8 * xi1a * (1+xi2*xi2a)(1+xi3*xi3a)   etc.
    """
    xi = np.asarray(xi, dtype=float)
    n = ABAQUS_C3D8_NODES
    a1, a2, a3 = n[:, 0], n[:, 1], n[:, 2]
    f1 = 1.0 + a1 * xi[0]
    f2 = 1.0 + a2 * xi[1]
    f3 = 1.0 + a3 * xi[2]
    dN = np.empty((8, 3), dtype=float)
    dN[:, 0] = 0.125 * a1 * f2 * f3
    dN[:, 1] = 0.125 * a2 * f1 * f3
    dN[:, 2] = 0.125 * a3 * f1 * f2
    return dN


# ---------------------------------------------------------------------------
# B matrices
# ---------------------------------------------------------------------------

def _b_from_coords(coords, xi):
    """Core B-matrix construction from a set of nodal coords (reference OR
    current) and a natural point. Returns (B (6,24), detJ, dNdx (8,3)).

    Jacobian  J_ij = d x_i / d xi_j = sum_a coords[a,i] * dN_a/dxi_j
              => J = coords^T @ dNdxi.
    Spatial gradient  dNdx = dNdxi @ inv(J).
    """
    coords = np.asarray(coords, dtype=float)
    dNdxi = shape_grad_natural(xi)             # (8,3)
    J = coords.T @ dNdxi                        # (3,3)
    detJ = np.linalg.det(J)
    invJ = np.linalg.inv(J)
    dNdx = dNdxi @ invJ                         # (8,3): dN_a/dx_i

    B = np.zeros((6, 24), dtype=float)
    for a in range(8):
        bx, by, bz = dNdx[a, 0], dNdx[a, 1], dNdx[a, 2]
        c = 3 * a
        # Voigt rows: 0:11, 1:22, 2:33, 3:12(eng), 4:13(eng), 5:23(eng)
        B[0, c + 0] = bx
        B[1, c + 1] = by
        B[2, c + 2] = bz
        B[3, c + 0] = by
        B[3, c + 1] = bx
        B[4, c + 0] = bz
        B[4, c + 2] = bx
        B[5, c + 1] = bz
        B[5, c + 2] = by
    return B, detJ, dNdx


def b_matrix_reference(Xe, xi):
    """Small-strain B in the REFERENCE configuration.

    Xe : (8,3) reference nodal coords.
    Returns (B0 (6,24), detJ0). detJ0 is the reference-config Jacobian det.
    Voigt order (11,22,33,12,13,23), engineering shear.
    """
    B, detJ, _ = _b_from_coords(Xe, xi)
    return B, detJ


def b_matrix_spatial(xe, xi):
    """B in the CURRENT (deformed) configuration.

    xe : (8,3) CURRENT nodal coords (xe = Xe + Ue).
    Returns (B (6,24), detJ) where detJ is the current-config Jacobian det.
    Same construction as reference but derivatives are w.r.t. current coords.
    """
    B, detJ, _ = _b_from_coords(xe, xi)
    return B, detJ


# ---------------------------------------------------------------------------
# Voigt helpers
# ---------------------------------------------------------------------------

def voigt_to_tensor(s):
    """(6,) Voigt (11,22,33,12,13,23) symmetric stress -> (3,3) tensor."""
    s = np.asarray(s, dtype=float)
    return np.array([[s[0], s[3], s[4]],
                     [s[3], s[1], s[5]],
                     [s[4], s[5], s[2]]], dtype=float)


def _as_ue_2d(Ue):
    """Accept (8,3) or (24,) nodal displacements -> (8,3)."""
    Ue = np.asarray(Ue, dtype=float)
    if Ue.shape == (24,):
        return Ue.reshape(8, 3)
    if Ue.shape == (8, 3):
        return Ue
    raise ValueError("Ue must have shape (8,3) or (24,), got %r" % (Ue.shape,))


# ---------------------------------------------------------------------------
# Internal forces
# ---------------------------------------------------------------------------

def element_internal_force_small_strain(Xe, sigma_ip, gauss=ABAQUS_C3D8_GAUSS):
    """Small-strain (reference-config, F = I) internal force. Patch-test path.

        f = sum_k  B0(Xe, gp_k)^T @ sigma_ip[k] * detJ0_k * w_k

    Xe       : (8,3) reference nodal coords.
    sigma_ip : (8,6) stress in Voigt order at the 8 IPs (Cauchy ~ stress).
    Returns  : (24,) node-major internal force.
    """
    Xe = np.asarray(Xe, dtype=float)
    sigma_ip = np.asarray(sigma_ip, dtype=float)
    points, weights = gauss.points, gauss.weights
    f = np.zeros(24, dtype=float)
    for k in range(len(weights)):
        B0, detJ0 = b_matrix_reference(Xe, points[k])
        f += (B0.T @ sigma_ip[k]) * detJ0 * weights[k]
    return f


def element_internal_force_finite_strain(Xe, Ue, sigma_ip, gauss=ABAQUS_C3D8_GAUSS):
    """Finite-strain (current-config) internal force -- the Abaqus-matching
    path for nlgeom=YES C3D8.

        xe = Xe + Ue
        f  = sum_k  B_spatial(xe, gp_k)^T @ sigma_ip[k] * detJ_current_k * w_k

    Xe       : (8,3) reference nodal coords.
    Ue       : (8,3) or (24,) nodal displacements.
    sigma_ip : (8,6) CAUCHY stress in Abaqus Voigt order at the 8 IPs.
    Returns  : (24,) node-major internal force.
    """
    Xe = np.asarray(Xe, dtype=float)
    xe = Xe + _as_ue_2d(Ue)
    sigma_ip = np.asarray(sigma_ip, dtype=float)
    points, weights = gauss.points, gauss.weights
    f = np.zeros(24, dtype=float)
    for k in range(len(weights)):
        B, detJ = b_matrix_spatial(xe, points[k])
        f += (B.T @ sigma_ip[k]) * detJ * weights[k]
    return f


# ---------------------------------------------------------------------------
# Element tangent
# ---------------------------------------------------------------------------

def element_tangent(Xe, Ue, Dmat_ip, sigma_ip=None, mode='finite',
                    gauss=ABAQUS_C3D8_GAUSS):
    """Element tangent stiffness (24,24).

    Material term (always):     K_mat = sum_k  B^T D_k B * detJ * w
    Geometric/initial-stress term (only when mode='finite' AND sigma_ip given):

        K_geo[(a,m),(b,n)] = delta_mn * sum_k ( dNdx[a] . sigma_k . dNdx[b] )
                                         * detJ_k * w_k

    i.e. the standard current-configuration initial-stress stiffness with the
    Cauchy stress. In block form for the 8 nodes this is  kron(G, I3)  with
    G[a,b] = dNdx[a] @ sigma_tensor @ dNdx[b]. It is SYMMETRIC and is the
    conventional geometric term used together with an objective-rate material
    tangent D to form the consistent nlgeom stiffness.

    IT IS NOT the directional derivative of  int B_spatial^T sigma dv  at FIXED
    sigma. That fixed-sigma force Jacobian is generally NON-symmetric and equals
    outer(sigma.g_a, g_b) - outer(sigma.g_b, g_a) per pair; it is provided
    separately as force_tangent_fixed_sigma() and differs from K_geo by exactly
    the objective-rate terms. A finite-difference of element_internal_force_
    finite_strain at fixed sigma reproduces force_tangent_fixed_sigma(), NOT this
    K_geo (verified numerically; see tests/cp_c3d8_umat/tangent_fd_check).

    Xe       : (8,3) reference coords.
    Ue       : (8,3)/(24,) displacements (used only for mode='finite').
    Dmat_ip  : (8,6,6) consistent material tangent per IP (or (6,6) broadcast).
    sigma_ip : (8,6) Cauchy stress per IP; required for the geometric term.
    mode     : 'finite' (current config, B_spatial + K_geo) or
               'small'  (reference config, B0, no geometric term).

    NOTE on measures: for the true nlgeom C3D8 the material term needs D to be
    the Cauchy/Jaumann (or chosen objective-rate) spatial tangent consistent
    with the UMAT DDSDDE mapping; that mapping is the error-prone item to be
    closed against Abaqus later (see the tests/ README). Here the geometric
    term is verified independently.
    """
    Xe = np.asarray(Xe, dtype=float)
    Dmat_ip = np.asarray(Dmat_ip, dtype=float)
    if Dmat_ip.shape == (6, 6):
        Dmat_ip = np.broadcast_to(Dmat_ip, (8, 6, 6))
    points, weights = gauss.points, gauss.weights

    if mode == 'finite':
        coords = Xe + _as_ue_2d(Ue if Ue is not None else np.zeros((8, 3)))
    elif mode == 'small':
        coords = Xe
    else:
        raise ValueError("mode must be 'finite' or 'small', got %r" % (mode,))

    K = np.zeros((24, 24), dtype=float)
    I3 = np.eye(3)
    for k in range(len(weights)):
        B, detJ, dNdx = _b_from_coords(coords, points[k])
        D = Dmat_ip[k]
        K += (B.T @ D @ B) * detJ * weights[k]
        if mode == 'finite' and sigma_ip is not None:
            S = voigt_to_tensor(np.asarray(sigma_ip)[k])   # (3,3) Cauchy
            G = dNdx @ S @ dNdx.T                           # (8,8) scalars
            K += np.kron(G, I3) * detJ * weights[k]
    return K


# ---------------------------------------------------------------------------
# Global assembly
# ---------------------------------------------------------------------------

def force_tangent_fixed_sigma(Xe, Ue, sigma_ip, gauss=ABAQUS_C3D8_GAUSS):
    """EXACT Jacobian d/dU of `element_internal_force_finite_strain` evaluated
    with the Cauchy stress array HELD FIXED (dsigma/dU = 0).

    This is the analytic derivative of the actual code path
        f_a,i(U) = integral_v  (dN_a/dx_j) sigma_ji  dv ,   x = X + U
    accounting for BOTH the change of the spatial gradient dN/dx and of the
    current volume dv with U. Writing g_a = grad N_a (current config):

        dF_a,i / dU_b,k = integral_v [ (sigma g_a)_i g_b,k
                                       - (sigma g_b)_i g_a,k ] dv

    so the (a,b) 3x3 block is  outer(sigma g_a, g_b) - outer(sigma g_b, g_a),
    integrated over the element.

    NOTE: this is NOT the standard geometric/initial-stress stiffness returned
    inside `element_tangent` (that is delta_ik * g_a.sigma.g_b, the physically
    conventional form used with an objective stress rate). The two differ by
    exactly the objective-rate terms. This helper exists only so the FD spot-
    check has the correct analytic target for the force routine as written.
    """
    Xe = np.asarray(Xe, dtype=float)
    xe = Xe + _as_ue_2d(Ue)
    sigma_ip = np.asarray(sigma_ip, dtype=float)
    points, weights = gauss.points, gauss.weights
    K = np.zeros((24, 24), dtype=float)
    for k in range(len(weights)):
        _, detJ, dNdx = _b_from_coords(xe, points[k])   # dNdx rows = g_a
        S = voigt_to_tensor(sigma_ip[k])
        Sg = dNdx @ S                                    # rows = sigma g_a
        scale = detJ * weights[k]
        for a in range(8):
            for b in range(8):
                block = np.outer(Sg[a], dNdx[b]) - np.outer(Sg[b], dNdx[a])
                K[3 * a:3 * a + 3, 3 * b:3 * b + 3] += block * scale
    return K


def assemble_global_internal_force(node_ids, coords, connectivity, U, sigma_all,
                                   mode='finite'):
    """Scatter element internal forces into a global vector.

    node_ids     : list[int] of all node ids present.
    coords       : (nnode,3) reference coords, row i corresponds to node_ids[i].
    connectivity : list of (eid, [8 node ids]).
    U            : (ndof,) global displacement (global DOF order) OR
                   dict[node_id] -> (3,) displacement.
    sigma_all    : dict[eid] -> (8,6) Cauchy stress at the 8 IPs (Voigt order).
    mode         : 'finite' (current config) or 'small' (reference config).

    Returns (F_global (ndof,), node_id_to_index) where node_id_to_index maps
    node id -> 0-based global node index (nodes sorted ascending), and
    ndof = 3 * nnode with node dofs [3i, 3i+1, 3i+2].
    """
    node_ids = list(node_ids)
    sorted_ids = sorted(node_ids)
    node_id_to_index = {nid: i for i, nid in enumerate(sorted_ids)}
    nnode = len(sorted_ids)
    ndof = 3 * nnode
    coords = np.asarray(coords, dtype=float)
    id_to_row = {nid: i for i, nid in enumerate(node_ids)}

    is_dict_U = isinstance(U, dict)
    U_arr = None if is_dict_U else np.asarray(U, dtype=float)

    def get_u(nid):
        if is_dict_U:
            return np.asarray(U[nid], dtype=float)
        gi = node_id_to_index[nid]
        return U_arr[3 * gi:3 * gi + 3]

    F = np.zeros(ndof, dtype=float)
    for eid, conn in connectivity:
        Xe = np.array([coords[id_to_row[nid]] for nid in conn], dtype=float)
        Ue = np.array([get_u(nid) for nid in conn], dtype=float)   # (8,3)
        sig = np.asarray(sigma_all[eid], dtype=float)              # (8,6)
        if mode == 'finite':
            fe = element_internal_force_finite_strain(Xe, Ue, sig)
        elif mode == 'small':
            fe = element_internal_force_small_strain(Xe, sig)
        else:
            raise ValueError("mode must be 'finite' or 'small'")
        for a, nid in enumerate(conn):
            gi = node_id_to_index[nid]
            F[3 * gi:3 * gi + 3] += fe[3 * a:3 * a + 3]
    return F, node_id_to_index


# ===========================================================================
# Verification helpers (Abaqus-independent) -- exercised by __main__ below and
# reused by tests/cp_c3d8_umat/tangent_fd_check/run_checks.py
# ===========================================================================

# Abaqus C3D8 face node lists (0-based within element connectivity), ordered
# around each face perimeter (S1..S6).
_C3D8_FACES = [
    [0, 1, 2, 3],   # S1  xi3 = -1
    [4, 7, 6, 5],   # S2  xi3 = +1
    [0, 4, 5, 1],   # S3  xi2 = -1
    [1, 5, 6, 2],   # S4  xi1 = +1
    [2, 6, 7, 3],   # S5  xi2 = +1
    [3, 7, 4, 0],   # S6  xi1 = -1
]

# Element 1 of Compression111.inp (connectivity 37,38,44,43, 1,2,8,7), verified
# against the .inp on 2026-07-10. Rows are in connectivity (node-major) order.
COMPRESSION111_ELEM1_XE = np.array([
    [8.0, 10.0, 10.0],   # node 37
    [8.0,  8.0, 10.0],   # node 38
    [8.0,  8.0,  8.0],   # node 44
    [8.0, 10.0,  8.0],   # node 43
    [10.0, 10.0, 10.0],  # node 1
    [10.0,  8.0, 10.0],  # node 2
    [10.0,  8.0,  8.0],  # node 8
    [10.0, 10.0,  8.0],  # node 7
], dtype=float)


def _face_bilinear(u, v):
    """4-node bilinear face shape functions and their (u,v) derivatives.

    Node local order matches the perimeter traversal of _C3D8_FACES.
    Returns (M (4,), dMdu (4,), dMdv (4,)).
    """
    M = 0.25 * np.array([(1 - u) * (1 - v),
                         (1 + u) * (1 - v),
                         (1 + u) * (1 + v),
                         (1 - u) * (1 + v)])
    dMdu = 0.25 * np.array([-(1 - v), (1 - v), (1 + v), -(1 + v)])
    dMdv = 0.25 * np.array([-(1 - u), -(1 + u), (1 + u), (1 - u)])
    return M, dMdu, dMdv


def surface_traction_nodal_forces(Xe, sigma0_voigt):
    """Consistent nodal forces of surface tractions t = sigma0 . n for a
    CONSTANT stress sigma0, integrated over the 6 C3D8 faces:

        f_a = integral_{dOmega} N_a t dA ,   t = sigma0 . n_outward

    Each face is a 4-node bilinear quad, integrated with 2x2 Gauss. The outward
    normal is guaranteed by orienting the surface area vector away from the
    element centroid. Returns (24,) node-major.

    By the divergence theorem (constant sigma0, div sigma0 = 0) this MUST equal
    element_internal_force_small_strain(Xe, [sigma0]*8) to numerical precision.
    """
    Xe = np.asarray(Xe, dtype=float)
    S = voigt_to_tensor(sigma0_voigt)
    centroid = Xe.mean(axis=0)
    g = 1.0 / np.sqrt(3.0)
    gp2 = [(-g, -g), (g, -g), (g, g), (-g, g)]   # 2x2 face Gauss, weights 1
    f = np.zeros(24, dtype=float)
    for face in _C3D8_FACES:
        Xf = Xe[face]                    # (4,3) face node coords
        face_center = Xf.mean(axis=0)
        for (u, v) in gp2:
            M, dMdu, dMdv = _face_bilinear(u, v)
            xg = M @ Xf                  # physical point
            t_u = dMdu @ Xf              # dx/du  (3,)
            t_v = dMdv @ Xf              # dx/dv  (3,)
            area_vec = np.cross(t_u, t_v)  # n * dA per unit du dv
            # orient outward (away from element centroid)
            if np.dot(area_vec, xg - centroid) < 0.0:
                area_vec = -area_vec
            traction_dA = S @ area_vec   # sigma0 . (n dA)   (weight = 1)
            for j, a in enumerate(face):
                f[3 * a:3 * a + 3] += M[j] * traction_dA
    return f


def linear_stress_target(Xe, sigma_fn, gauss=ABAQUS_C3D8_GAUSS):
    """Independent nodal-force target for a spatially-VARYING stress field.

    For any differentiable Cauchy stress sigma(x), integration by parts gives the
    exact identity for the element internal force

        f_a = int_Omega  grad(N_a) . sigma dV
            = int_dOmega  N_a (sigma . n) dA  -  int_Omega  N_a (div sigma) dV .

    This target is computed WITHOUT using the volume Gauss-point stress pairing
    (surface term from face quadrature; volume term from int N_a dV times the
    constant-for-linear div sigma), so comparing it to
    element_internal_force_small_strain (which DOES pair sigma_ip[k] with
    points[k]) exercises the per-IP weights, detJ, and the sigma<->point pairing
    for a NON-uniform field -- unlike the uniform-stress patch test, which is
    blind to all of those. For an affine element map (parallelepiped, e.g. the
    Compression111 cubes) and a linear sigma, both sides are integrated exactly
    by 2x2(x2) Gauss, so they must agree to machine precision.

    sigma_fn : callable x(3,) -> sigma Voigt(6,) (Abaqus order), linear in x.
    Returns (24,) node-major target force.
    """
    Xe = np.asarray(Xe, dtype=float)

    def sig_tensor(x):
        return voigt_to_tensor(np.asarray(sigma_fn(x), dtype=float))

    # constant divergence of a linear field, via central differences (exact)
    h = 1e-3
    b = np.zeros(3)
    x0 = Xe.mean(axis=0)
    for j in range(3):
        ep = x0.copy(); ep[j] += h
        em = x0.copy(); em[j] -= h
        b += (sig_tensor(ep)[:, j] - sig_tensor(em)[:, j]) / (2 * h)   # d sigma_ij/dx_j

    points, weights = gauss.points, gauss.weights

    # volume term: - int N_a dV * b
    f = np.zeros(24, dtype=float)
    for k in range(len(weights)):
        N = shape_functions(points[k])
        _, detJ, _ = _b_from_coords(Xe, points[k])
        vol_w = detJ * weights[k]
        for a in range(8):
            f[3 * a:3 * a + 3] += -(N[a] * vol_w) * b

    # surface term: int N_a (sigma(x).n) dA over the 6 faces
    centroid = Xe.mean(axis=0)
    g = 1.0 / np.sqrt(3.0)
    gp2 = [(-g, -g), (g, -g), (g, g), (-g, g)]
    for face in _C3D8_FACES:
        Xf = Xe[face]
        for (u, v) in gp2:
            M, dMdu, dMdv = _face_bilinear(u, v)
            xg = M @ Xf
            area_vec = np.cross(dMdu @ Xf, dMdv @ Xf)
            if np.dot(area_vec, xg - centroid) < 0.0:
                area_vec = -area_vec
            traction_dA = sig_tensor(xg) @ area_vec
            for j, a in enumerate(face):
                f[3 * a:3 * a + 3] += M[j] * traction_dA
    return f


def isotropic_D(E, nu):
    """Isotropic linear-elastic 6x6 tangent in Abaqus Voigt order with
    ENGINEERING shear (so it pairs with the engineering-shear B matrix).

    Private duplicate of the canonical residual_core/core/voigt.py::isotropic_D,
    kept here only so this kernel stays runnable as a bare script. Materials
    import the core copy (avoids a materials->formulations dependency)."""
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    D = np.zeros((6, 6), dtype=float)
    D[:3, :3] = lam
    for i in range(3):
        D[i, i] += 2 * mu
    for i in range(3, 6):
        D[i, i] = mu          # engineering shear: tau = mu * gamma
    return D


def unit_cube_Xe():
    """C3D8 unit cube in Abaqus node order (node-major (8,3))."""
    return ABAQUS_C3D8_NODES.copy() * 0.5 + 0.5   # maps [-1,1] -> [0,1]


# ===========================================================================
# __main__: run the three Abaqus-independent checks and print PASS/FAIL
# ===========================================================================

def _check_divergence_theorem(verbose=True):
    """Check 1: divergence-theorem patch test on the real Compression111
    element 1, and on a distorted version of it. Returns (ok, details)."""
    rng = np.random.default_rng(12345)
    results = []
    geoms = {
        "Compression111 elem1 (real)": COMPRESSION111_ELEM1_XE.copy(),
        "distorted hex": COMPRESSION111_ELEM1_XE
        + rng.uniform(-0.6, 0.6, size=(8, 3)),
    }
    # a random symmetric stress and a uniaxial stress
    srand = rng.uniform(-100, 100, size=6)
    stresses = {
        "random symmetric": srand,
        "uniaxial sigma33=250": np.array([0, 0, 250.0, 0, 0, 0]),
    }
    ok_all = True
    for gname, Xe in geoms.items():
        for sname, s0 in stresses.items():
            f_vol = element_internal_force_small_strain(Xe, np.tile(s0, (8, 1)))
            f_surf = surface_traction_nodal_forces(Xe, s0)
            denom = max(np.linalg.norm(f_surf), 1e-30)
            rel = np.linalg.norm(f_vol - f_surf) / denom
            equil = np.linalg.norm(f_vol.reshape(8, 3).sum(axis=0))
            ok = (rel < 1e-8) and (equil < 1e-8 * max(denom, 1.0))
            ok_all = ok_all and ok
            results.append((gname, sname, rel, equil, ok))
            if verbose:
                print("   [%s] %-22s rel=%.3e  |sum F|=%.3e  %s"
                      % ("OK " if ok else "BAD", "%s / %s" % (gname, sname),
                         rel, equil, "PASS" if ok else "FAIL"))
    return ok_all, results


def _check_uniaxial(verbose=True):
    """Check 2: uniform uniaxial sanity on a unit cube."""
    Xe = unit_cube_Xe()
    S33 = 300.0
    s0 = np.array([0, 0, S33, 0, 0, 0.0])
    f = element_internal_force_small_strain(Xe, np.tile(s0, (8, 1))).reshape(8, 3)
    # top face (z=1): nodes 5,6,7,8 (0-based 4..7). bottom (z=0): 0..3.
    top = [4, 5, 6, 7]
    bot = [0, 1, 2, 3]
    A = 1.0
    expect_top = S33 * A / 4.0
    err = 0.0
    err = max(err, np.max(np.abs(f[top, 2] - expect_top)))
    err = max(err, np.max(np.abs(f[bot, 2] + expect_top)))
    err = max(err, np.max(np.abs(f[:, 0])))   # x comps zero
    err = max(err, np.max(np.abs(f[:, 1])))   # y comps zero
    ok = err < 1e-10
    if verbose:
        print("   [%s] +z nodes carry +S33*A/4=%.4f, -z carry -, sides 0; "
              "max err=%.3e  %s" % ("OK " if ok else "BAD", expect_top, err,
                                    "PASS" if ok else "FAIL"))
    return ok, err


def _tangent_fd(Xe, D, h_scale=1e-6):
    """Central-difference the small-strain residual r(U)=int B0^T sigma(U) dV,
    sigma(U)=D @ (B0 U) per IP, and compare to K = sum B0^T D B0 detJ0 w."""
    points, weights = ABAQUS_C3D8_GAUSS.points, ABAQUS_C3D8_GAUSS.weights
    # precompute per-IP B0, detJ0
    B0s, dets = [], []
    for k in range(8):
        B0, d0 = b_matrix_reference(Xe, points[k])
        B0s.append(B0)
        dets.append(d0)

    def residual(U):
        sig = np.array([D @ (B0s[k] @ U) for k in range(8)])   # (8,6)
        return element_internal_force_small_strain(Xe, sig)

    # analytic K via element_tangent (mode='small')
    K = element_tangent(Xe, np.zeros(24), D, mode='small')

    # finite-difference Jacobian
    U0 = np.zeros(24)
    scale = np.linalg.norm(Xe) / np.sqrt(8)
    h = h_scale * max(scale, 1.0)
    Kfd = np.zeros((24, 24))
    for j in range(24):
        Up = U0.copy(); Up[j] += h
        Um = U0.copy(); Um[j] -= h
        Kfd[:, j] = (residual(Up) - residual(Um)) / (2 * h)
    abs_err = np.max(np.abs(Kfd - K))
    fro = np.linalg.norm(Kfd - K) / max(np.linalg.norm(K), 1e-30)
    max_entry = abs_err
    return K, Kfd, abs_err, fro, max_entry


def _check_tangent(verbose=True):
    """Check 3a: linear tangent finite-difference check."""
    Xe = COMPRESSION111_ELEM1_XE.copy()
    D = isotropic_D(E=200000.0, nu=0.3)
    K, Kfd, abs_err, fro, max_entry = _tangent_fd(Xe, D)
    ok = fro < 1e-6
    if verbose:
        print("   [%s] linear tangent FD: abs=%.3e rel(Fro)=%.3e max=%.3e  %s"
              % ("OK " if ok else "BAD", abs_err, fro, max_entry,
                 "PASS" if ok else "FAIL"))
    return ok, (abs_err, fro, max_entry)


def _check_finite_strain_fd(verbose=True, U_amp=0.0):
    """Check 3b: finite-strain FD spot-check.

    Central-difference d/dU of `element_internal_force_finite_strain` at a FIXED
    Cauchy-stress field and compare to the EXACT analytic Jacobian-at-fixed-sigma
    `force_tangent_fixed_sigma`.

    What it PROVES: the finite-strain internal-force routine is smooth in U and
    its current-config kinematics (spatial gradient dN/dx and current volume
    detJ, both varying with U) are differentiated correctly -- the assembly
    machinery is consistent.

    What it does NOT prove: the physical finite-strain element tangent. That
    additionally needs the material derivative dsigma/dU (the UMAT DDSDDE mapped
    to a spatial/objective-rate tangent) plus the conventional geometric split.
    That mapping is the error-prone item to be closed later against Abaqus.
    """
    rng = np.random.default_rng(7)
    Xe = COMPRESSION111_ELEM1_XE.copy()
    sigma0 = rng.uniform(-50, 50, size=6)
    sigma_ip = np.tile(sigma0, (8, 1))
    # baseline displacement: exercise a genuinely deformed config if U_amp>0
    U0 = rng.uniform(-U_amp, U_amp, size=24) if U_amp > 0 else np.zeros(24)

    def force(U):
        return element_internal_force_finite_strain(Xe, U, sigma_ip)

    scale = np.linalg.norm(Xe) / np.sqrt(8)
    h = 1e-6 * max(scale, 1.0)
    Kfd = np.zeros((24, 24))
    for j in range(24):
        Up = U0.copy(); Up[j] += h
        Um = U0.copy(); Um[j] -= h
        Kfd[:, j] = (force(Up) - force(Um)) / (2 * h)

    Kan = force_tangent_fixed_sigma(Xe, U0, sigma_ip)
    abs_err = np.max(np.abs(Kfd - Kan))
    fro = np.linalg.norm(Kfd - Kan) / max(np.linalg.norm(Kan), 1e-30)
    ok = fro < 1e-6
    if verbose:
        print("   [%s] finite-strain force FD vs exact fixed-sigma J (U_amp=%.2f): "
              "abs=%.3e rel(Fro)=%.3e  %s"
              % ("OK " if ok else "BAD", U_amp, abs_err, fro,
                 "PASS" if ok else "FAIL"))
    return ok, (abs_err, fro)


def _check_linear_stress(verbose=True):
    """Check 1b: LINEAR (spatially varying) stress patch test.

    Unlike the uniform patch test, a linear stress makes the per-IP weights and
    the sigma_ip[k] <-> points[k] pairing matter, and gives node-specific target
    forces (so it also catches intra-element force redistribution). Run on the
    real Compression111 cube element (affine map -> quadrature exact -> machine
    precision). Still does NOT probe the Abaqus IP *ordering* (we sample sigma at
    our own points), which requires an ODB.
    """
    rng = np.random.default_rng(2024)
    # linear symmetric stress:  sigma(x) = S0 + Sx*x + Sy*y + Sz*z
    S0 = rng.uniform(-50, 50, size=(3, 3)); S0 = 0.5 * (S0 + S0.T)
    grads = [0.5 * (G + G.T) for G in rng.uniform(-8, 8, size=(3, 3, 3))]

    def sigma_fn(x):
        M = S0 + grads[0] * x[0] + grads[1] * x[1] + grads[2] * x[2]
        return np.array([M[0, 0], M[1, 1], M[2, 2], M[0, 1], M[0, 2], M[1, 2]])

    ok_all = True
    for gname, Xe in (("Compression111 elem1 (cube)", COMPRESSION111_ELEM1_XE.copy()),):
        # sample sigma at each physical Gauss point and assemble
        pts = ABAQUS_C3D8_GAUSS.points
        sigma_ip = np.array([sigma_fn(shape_functions(pts[k]) @ Xe) for k in range(8)])
        f_asm = element_internal_force_small_strain(Xe, sigma_ip)
        f_tgt = linear_stress_target(Xe, sigma_fn)
        denom = max(np.linalg.norm(f_tgt), 1e-30)
        rel = np.linalg.norm(f_asm - f_tgt) / denom
        ok = rel < 1e-9
        ok_all = ok_all and ok
        if verbose:
            print("   [%s] %-30s rel=%.3e  (tol 1e-9)  %s"
                  % ("OK " if ok else "BAD", gname, rel, "PASS" if ok else "FAIL"))
    return ok_all, None


def run_all(verbose=True):
    if verbose:
        print("=" * 72)
        print("C3D8 residual assembler -- Abaqus-independent verification")
        print("=" * 72)
        print("Check 1: divergence-theorem patch test (real + distorted hex)")
    ok1, _ = _check_divergence_theorem(verbose)
    if verbose:
        print("Check 1b: LINEAR-stress patch test (per-IP weights + pairing)")
    ok1b, _ = _check_linear_stress(verbose)
    if verbose:
        print("Check 2: uniform uniaxial sanity (unit cube)")
    ok2, _ = _check_uniaxial(verbose)
    if verbose:
        print("Check 3a: linear tangent finite-difference")
    ok3, _ = _check_tangent(verbose)
    if verbose:
        print("Check 3b: finite-strain force finite-difference (undeformed + deformed)")
    ok4a, _ = _check_finite_strain_fd(verbose, U_amp=0.0)
    ok4b, _ = _check_finite_strain_fd(verbose, U_amp=0.5)
    ok4 = ok4a and ok4b
    all_ok = ok1 and ok1b and ok2 and ok3 and ok4
    if verbose:
        print("-" * 72)
        print("OVERALL: %s" % ("ALL PASS" if all_ok else "FAILURE"))
        print("=" * 72)
    return all_ok


if __name__ == "__main__":
    import sys
    sys.exit(0 if run_all(verbose=True) else 1)
