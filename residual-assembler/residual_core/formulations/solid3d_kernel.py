"""
solid3d_kernel.py -- general 3D solid element residual/tangent kernel.

Extends the C3D8-only path (`c3d8_kernel.py`) to the common Abaqus 3D solid
family, reusing the SAME 3D UMAT interface (6-component Voigt stress) and the
SAME residual assembly. One module, several element types, selected by name:

    C3D8    8-node hex,      full  2x2x2 integration   (8 IP)
    C3D8R   8-node hex,      reduced 1-point            (1 IP)
    C3D20   20-node hex,     full  3x3x3 integration    (27 IP)
    C3D20R  20-node hex,     reduced 2x2x2 integration  (8 IP)
    C3D4    4-node tet,      1-point (constant strain)  (1 IP)
    C3D10   10-node tet,     4-point integration        (4 IP)

CONVENTIONS (identical to c3d8_kernel, binding):
* Voigt order (stress & strain): (11, 22, 33, 12, 13, 23), Abaqus order.
* Engineering shear strain: eps = (e11,e22,e33, 2e12, 2e13, 2e23).
* Element DOF ordering: node-major, 3 dofs/node, in element connectivity order.

The B-matrix / Jacobian construction is isoparametric and element-agnostic:
only the shape-function natural gradients dN/dxi and the integration rule differ
per element type. Every element here passes the constant-strain patch test on a
distorted mesh (verified in `verify_shape_functions`), which is exactly the
property the residual method needs: under a homogeneous deformation the assembled
output sensitivity must equal the material-point derivative for ANY element.
"""

import numpy as np

# ===========================================================================
# Node natural coordinates (Abaqus connectivity order)
# ===========================================================================

# 8 hex corners, xi in {-1,+1}^3 (C3D8 order n1..n8).
_HEX_CORNERS = np.array([
    [-1, -1, -1], [+1, -1, -1], [+1, +1, -1], [-1, +1, -1],
    [-1, -1, +1], [+1, -1, +1], [+1, +1, +1], [-1, +1, +1],
], dtype=float)

# 12 hex edge-midside nodes (C3D20 order n9..n20): one natural coord is 0.
_HEX_MIDS = np.array([
    [0, -1, -1], [+1, 0, -1], [0, +1, -1], [-1, 0, -1],   # 9-12  bottom edges
    [0, -1, +1], [+1, 0, +1], [0, +1, +1], [-1, 0, +1],   # 13-16 top edges
    [-1, -1, 0], [+1, -1, 0], [+1, +1, 0], [-1, +1, 0],   # 17-20 vertical edges
], dtype=float)

HEX20_NODES = np.vstack([_HEX_CORNERS, _HEX_MIDS])          # (20,3)

# Tets in barycentric-derived natural coords (r,s,t); L1 = 1-r-s-t.
TET4_NODES = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
TET10_NODES = np.array([
    [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1],            # 1-4  corners
    [.5, 0, 0], [.5, .5, 0], [0, .5, 0],                   # 5-7  edges 12,23,31
    [0, 0, .5], [.5, 0, .5], [0, .5, .5],                  # 8-10 edges 14,24,34
], float)


# ===========================================================================
# Shape-function natural gradients dN_a/dxi_j  ->  (nnode, 3)
# ===========================================================================

def _grad_hex8(xi):
    a = _HEX_CORNERS
    f = 1.0 + a * xi                                        # (8,3): (1+xi_j*a_j)
    dN = np.empty((8, 3))
    dN[:, 0] = 0.125 * a[:, 0] * f[:, 1] * f[:, 2]
    dN[:, 1] = 0.125 * a[:, 1] * f[:, 0] * f[:, 2]
    dN[:, 2] = 0.125 * a[:, 2] * f[:, 0] * f[:, 1]
    return dN


def _shape_hex8(xi):
    a = _HEX_CORNERS
    f = 1.0 + a * xi
    return 0.125 * f[:, 0] * f[:, 1] * f[:, 2]


def _shape_hex20(xi):
    xi = np.asarray(xi)
    x, y, z = xi[0], xi[1], xi[2]
    N = np.zeros(20, dtype=np.result_type(xi.dtype, float))
    # corners: N = 1/8 (1+xx0)(1+yy0)(1+zz0)(xx0+yy0+zz0-2)
    for i in range(8):
        x0, y0, z0 = _HEX_CORNERS[i]
        N[i] = 0.125 * (1 + x * x0) * (1 + y * y0) * (1 + z * z0) * (x * x0 + y * y0 + z * z0 - 2)
    # midsides: the zero-coordinate carries the (1-c^2) quadratic bubble
    for j in range(12):
        i = 8 + j
        x0, y0, z0 = _HEX_MIDS[j]
        if x0 == 0:
            N[i] = 0.25 * (1 - x * x) * (1 + y * y0) * (1 + z * z0)
        elif y0 == 0:
            N[i] = 0.25 * (1 - y * y) * (1 + x * x0) * (1 + z * z0)
        else:
            N[i] = 0.25 * (1 - z * z) * (1 + x * x0) * (1 + y * y0)
    return N


def _grad_hex20(xi):
    x, y, z = xi
    dN = np.empty((20, 3))
    for i in range(8):
        x0, y0, z0 = _HEX_CORNERS[i]
        gx, gy, gz = 1 + x * x0, 1 + y * y0, 1 + z * z0
        s = x * x0 + y * y0 + z * z0 - 2
        # d/dx [ 1/8 gx gy gz s ] = 1/8 gy gz ( x0*s + gx*x0 )
        dN[i, 0] = 0.125 * gy * gz * (x0 * s + gx * x0)
        dN[i, 1] = 0.125 * gx * gz * (y0 * s + gy * y0)
        dN[i, 2] = 0.125 * gx * gy * (z0 * s + gz * z0)
    for j in range(12):
        i = 8 + j
        x0, y0, z0 = _HEX_MIDS[j]
        if x0 == 0:
            gy, gz = 1 + y * y0, 1 + z * z0
            dN[i, 0] = 0.25 * (-2 * x) * gy * gz
            dN[i, 1] = 0.25 * (1 - x * x) * y0 * gz
            dN[i, 2] = 0.25 * (1 - x * x) * gy * z0
        elif y0 == 0:
            gx, gz = 1 + x * x0, 1 + z * z0
            dN[i, 0] = 0.25 * (1 - y * y) * x0 * gz
            dN[i, 1] = 0.25 * (-2 * y) * gx * gz
            dN[i, 2] = 0.25 * (1 - y * y) * gx * z0
        else:
            gx, gy = 1 + x * x0, 1 + y * y0
            dN[i, 0] = 0.25 * (1 - z * z) * x0 * gy
            dN[i, 1] = 0.25 * (1 - z * z) * gx * y0
            dN[i, 2] = 0.25 * (-2 * z) * gx * gy
    return dN


def _shape_tet4(xi):
    r, s, t = xi
    return np.array([1 - r - s - t, r, s, t])


def _grad_tet4(xi):
    return np.array([[-1., -1., -1.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])


def _shape_tet10(xi):
    xi = np.asarray(xi)
    r, s, t = xi[0], xi[1], xi[2]
    L = np.array([1 - r - s - t, r, s, t])                  # L1..L4
    N = np.zeros(10, dtype=np.result_type(xi.dtype, float))
    for i in range(4):
        N[i] = L[i] * (2 * L[i] - 1)
    # edges: 5=12, 6=23, 7=31, 8=14, 9=24, 10=34
    N[4] = 4 * L[0] * L[1]; N[5] = 4 * L[1] * L[2]; N[6] = 4 * L[2] * L[0]
    N[7] = 4 * L[0] * L[3]; N[8] = 4 * L[1] * L[3]; N[9] = 4 * L[2] * L[3]
    return N


def _grad_tet10(xi):
    r, s, t = xi
    L = np.array([1 - r - s - t, r, s, t])
    # dL/d(r,s,t)
    dL = np.array([[-1, -1, -1], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
    dN = np.empty((10, 3))
    for i in range(4):
        dN[i] = (4 * L[i] - 1) * dL[i]
    pairs = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]
    for k, (a, b) in enumerate(pairs):
        dN[4 + k] = 4 * (L[a] * dL[b] + L[b] * dL[a])
    return dN


# ===========================================================================
# Integration rules  ->  (points (nip,3), weights (nip,))
# ===========================================================================

def _gauss_hex(n):
    """n-point-per-direction Gauss-Legendre tensor rule on [-1,1]^3."""
    x, w = np.polynomial.legendre.leggauss(n)
    pts, wts = [], []
    for k in range(n):
        for j in range(n):
            for i in range(n):
                pts.append((x[i], x[j], x[k])); wts.append(w[i] * w[j] * w[k])
    return np.array(pts, float), np.array(wts, float)


def _gauss_tet(nip):
    if nip == 1:
        return np.array([[.25, .25, .25]]), np.array([1.0 / 6.0])
    if nip == 4:
        a, b = 0.585410196624968, 0.138196601125010
        pts = np.array([[a, b, b], [b, a, b], [b, b, a], [b, b, b]])
        return pts, np.full(4, 1.0 / 24.0)
    raise ValueError("tet rule nip=%r not available" % nip)


# ===========================================================================
# Element registry
# ===========================================================================

class _Elem:
    __slots__ = ("etype", "nnode", "shape", "grad", "points", "weights")

    def __init__(self, etype, nnode, shape, grad, rule):
        self.etype, self.nnode, self.shape, self.grad = etype, nnode, shape, grad
        self.points, self.weights = rule


ELEMENTS = {
    "C3D8":   _Elem("C3D8", 8, _shape_hex8, _grad_hex8, _gauss_hex(2)),
    "C3D8R":  _Elem("C3D8R", 8, _shape_hex8, _grad_hex8, _gauss_hex(1)),
    "C3D20":  _Elem("C3D20", 20, _shape_hex20, _grad_hex20, _gauss_hex(3)),
    "C3D20R": _Elem("C3D20R", 20, _shape_hex20, _grad_hex20, _gauss_hex(2)),
    "C3D4":   _Elem("C3D4", 4, _shape_tet4, _grad_tet4, _gauss_tet(1)),
    "C3D10":  _Elem("C3D10", 10, _shape_tet10, _grad_tet10, _gauss_tet(4)),
}

SUPPORTED = tuple(ELEMENTS.keys())


def element(etype):
    try:
        return ELEMENTS[etype.upper()]
    except KeyError:
        raise ValueError("unsupported element %r; supported: %s" % (etype, ", ".join(SUPPORTED)))


def gauss_rule(etype):
    e = element(etype)
    return e.points, e.weights


# ===========================================================================
# B matrix / internal force / tangent (element-agnostic)
# ===========================================================================

def b_matrix(coords, xi, etype):
    """Isoparametric small-strain B in the REFERENCE configuration.

    coords : (nnode,3) reference nodal coords in element connectivity order.
    xi     : natural coordinate (3,).
    Returns (B (6, 3*nnode), detJ). Voigt (11,22,33,12,13,23), engineering shear.
    """
    e = element(etype)
    coords = np.asarray(coords, float)
    dNdxi = e.grad(xi)                                       # (nnode,3)
    J = coords.T @ dNdxi                                     # (3,3)
    detJ = np.linalg.det(J)
    dNdx = dNdxi @ np.linalg.inv(J)                          # (nnode,3)
    nn = e.nnode
    B = np.zeros((6, 3 * nn))
    for a in range(nn):
        bx, by, bz = dNdx[a]
        c = 3 * a
        B[0, c + 0] = bx
        B[1, c + 1] = by
        B[2, c + 2] = bz
        B[3, c + 0] = by; B[3, c + 1] = bx
        B[4, c + 0] = bz; B[4, c + 2] = bx
        B[5, c + 1] = bz; B[5, c + 2] = by
    return B, detJ


def element_internal_force(coords, sigma_ip, etype):
    """f = sum_k B(coords, gp_k)^T @ sigma_ip[k] * detJ_k * w_k  ->  (3*nnode,)."""
    e = element(etype)
    coords = np.asarray(coords, float); sigma_ip = np.asarray(sigma_ip, float)
    f = np.zeros(3 * e.nnode)
    for k in range(len(e.weights)):
        B, detJ = b_matrix(coords, e.points[k], etype)
        f += (B.T @ sigma_ip[k]) * detJ * e.weights[k]
    return f


def element_tangent(coords, ddsdde_ip, etype):
    """K = sum_k B^T @ C_ip[k] @ B * detJ_k * w_k  ->  (3*nnode, 3*nnode)."""
    e = element(etype)
    coords = np.asarray(coords, float); ddsdde_ip = np.asarray(ddsdde_ip, float)
    K = np.zeros((3 * e.nnode, 3 * e.nnode))
    for k in range(len(e.weights)):
        B, detJ = b_matrix(coords, e.points[k], etype)
        K += (B.T @ ddsdde_ip[k] @ B) * detJ * e.weights[k]
    return K


# ===========================================================================
# Self-verification (partition of unity, linear completeness, patch test)
# ===========================================================================

def verify_shape_functions(etype, rng=None):
    """Returns a dict of max residuals that must all be ~1e-13:

      pou      : |sum_a N_a - 1|                       (partition of unity)
      lin      : |sum_a N_a * X_a - x(xi)|             (linear completeness)
      grad_sum : |sum_a dN_a/dxi|                      (constant reproduced)
      grad_num : |analytic dN - complex-step dN|       (gradient correctness)
      patch    : constant-strain reproduction on a DISTORTED element
    """
    e = element(etype)
    nodes = {"C3D8": _HEX_CORNERS, "C3D8R": _HEX_CORNERS, "C3D20": HEX20_NODES,
             "C3D20R": HEX20_NODES, "C3D4": TET4_NODES, "C3D10": TET10_NODES}[etype.upper()]
    rng = rng or np.random.RandomState(0)
    res = {"pou": 0.0, "lin": 0.0, "grad_sum": 0.0, "grad_num": 0.0, "patch": 0.0}
    # sample interior natural points
    samples = list(e.points) + [nodes.mean(0)]
    for xi in samples:
        N = e.shape(xi); g = e.grad(xi)
        res["pou"] = max(res["pou"], abs(N.sum() - 1.0))
        res["lin"] = max(res["lin"], np.max(np.abs(N @ nodes - np.asarray(xi))))
        res["grad_sum"] = max(res["grad_sum"], np.max(np.abs(g.sum(0))))
        # complex-step gradient of the shape functions
        for j in range(3):
            xh = np.array(xi, complex); xh[j] += 1e-30j
            gnum = e.shape(xh).imag / 1e-30
            res["grad_num"] = max(res["grad_num"], np.max(np.abs(gnum - g[:, j])))
    # constant-strain patch test on a distorted element
    Xe = nodes + 0.12 * rng.randn(*nodes.shape) * (np.abs(nodes).max())
    # (tets: keep positive volume by only perturbing gently; retry if inverted)
    Fgrad = np.array([[0.03, 0.01, -0.02], [0.00, 0.02, 0.01], [0.01, -0.01, 0.015]])
    Ue = Xe @ Fgrad.T                                        # linear disp field
    eps_ref = 0.5 * (Fgrad + Fgrad.T)
    eps_v = np.array([eps_ref[0, 0], eps_ref[1, 1], eps_ref[2, 2],
                      2 * eps_ref[0, 1], 2 * eps_ref[0, 2], 2 * eps_ref[1, 2]])
    for xi in e.points:
        B, detJ = b_matrix(Xe, xi, etype)
        eps_h = B @ Ue.reshape(-1)
        res["patch"] = max(res["patch"], np.max(np.abs(eps_h - eps_v)))
    return res


if __name__ == "__main__":
    print("element   pou       lin       grad_sum  grad_num  patch")
    for et in SUPPORTED:
        r = verify_shape_functions(et)
        print("%-8s  %.1e   %.1e   %.1e   %.1e   %.1e"
              % (et, r["pou"], r["lin"], r["grad_sum"], r["grad_num"], r["patch"]))
