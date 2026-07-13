"""Stress-driven residual driver  (Verification Mode 1, project Step 5).

Ties the pieces together:

    abaqus_inp_parser.parse_inp   (mesh, sections, materials, BCs)
        +
    c3d8_residual.assemble_global_internal_force   (f_int = integral B^T sigma dv)
        +
    fields.json   (nodal U, integration-point Cauchy stress S, reactions RF,
                   exported from an Abaqus ODB by extract_abaqus_fields.py)

and answers the Step-5 verification question WITHOUT trusting our UMAT replay:

    At a converged Abaqus increment, using Abaqus' OWN exported integration-point
    stresses and displacements,

      * for FREE dofs           :  r = f_int - f_ext  ~  0
      * for PRESCRIBED-disp dofs :  f_int  ~  Abaqus reaction RF   (up to sign)

If that holds, our finite-element residual assembly reproduces what Abaqus does
internally -- the first gold-standard check, independent of the UMAT adapter.

--------------------------------------------------------------------------------
Two modes
--------------------------------------------------------------------------------
1. ``--fields fields.json``  : the real comparison against an Abaqus ODB export.
   (Requires an Abaqus run; Abaqus is not installed in this dev environment, so
   this path is delivered ready-to-run and is exercised the moment a fields.json
   exists.)

2. ``--selftest`` (default when no --fields given) : an ABAQUS-INDEPENDENT
   manufactured-equilibrium test of the whole pipeline.  Impose a *uniform*
   Cauchy stress field over the parsed mesh.  For a uniform stress:
       - divergence is zero  =>  assembled internal force at every INTERIOR node
         must vanish (this is the global analogue of the free-dof residual ~ 0),
       - boundary nodes carry exactly the surface traction sigma.n (the analogue
         of "internal force = reaction"),
       - the total force sums to zero.
   This verifies parse -> global scatter -> interior/boundary dof bookkeeping to
   machine precision without needing Abaqus.  (The single-element proof of
   integral B^T sigma dV itself lives in tests/cp_c3d8_umat/tangent_fd_check.)

Convention (see residual_core/CONTRACT.md): these CP examples are nlgeom=YES,
the UMAT returns Cauchy stress (Abaqus Voigt order 11,22,33,12,13,23), and the
internal force is assembled in the CURRENT configuration f_int = int B_spatial^T
sigma dv.  The self-test uses U=0 (current == reference) purely to exercise the
plumbing; the --fields path uses the real exported U with mode='finite'.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

# Post-refactor locations: parser under io/, C3D8 numerics under formulations/.
from residual_core.io import abaqus_inp_parser as parser        # noqa: E402
from residual_core.formulations import c3d8_kernel as c3d8       # noqa: E402


# --------------------------------------------------------------------------- #
# Mesh / DOF helpers derived from a parsed AbaqusModel
# --------------------------------------------------------------------------- #
def build_mesh(model, elem_type="C3D8"):
    """Return (node_ids, coords, connectivity) for the C3D8 part of the model.

    node_ids     : sorted list[int]
    coords       : (nnode,3) float, row i == node_ids[i]
    connectivity : list of (eid, [8 node ids]) for every C3D8 element
    """
    node_ids = sorted(model.nodes.keys())
    coords = np.array([model.nodes[nid] for nid in node_ids], dtype=float)
    connectivity = []
    for eid, el in sorted(model.elements.items()):
        if el.etype.upper() == elem_type.upper():
            connectivity.append((eid, list(el.connectivity)))
    return node_ids, coords, connectivity


def _resolve_target_nodes(model, target):
    """A boundary/load target -> list of node ids (int passthrough or set lookup)."""
    if isinstance(target, int):
        return [target]
    # case-insensitive set lookup
    if target in model.node_sets:
        return list(model.node_sets[target])
    for name, ids in model.node_sets.items():
        if name.upper() == str(target).upper():
            return list(ids)
    # a bare integer written as string
    try:
        return [int(target)]
    except (TypeError, ValueError):
        return []


# Abaqus symmetry/antisymmetry BC -> constrained TRANSLATIONAL dofs (1=x,2=y,3=z).
# XSYMM fixes U1 (+ rotations, absent for C3D8); XASYMM fixes U2,U3; etc.
_SYMM_DOF = {"XSYMM": [1], "YSYMM": [2], "ZSYMM": [3],
             "XASYMM": [2, 3], "YASYMM": [1, 3], "ZASYMM": [1, 2],
             "ENCASTRE": [1, 2, 3], "PINNED": [1, 2, 3]}


def prescribed_dofs(model, node_id_to_index):
    """Global 0-based dof indices constrained by *Boundary, with their values.

    Returns dict {global_dof_index: prescribed_value}.  Symmetry BCs contribute
    value 0.0.  Abaqus dof d in {1,2,3} maps to global dof 3*node_index+(d-1);
    rotational dofs (d>3) are ignored for C3D8.
    """
    out = {}
    for b in model.boundaries:
        nids = _resolve_target_nodes(model, b.target)
        if b.kind in _SYMM_DOF:
            dofs = _SYMM_DOF[b.kind]
            value = 0.0
        else:  # numeric value form
            dofs = list(range(b.dof_start, b.dof_end + 1))
            value = float(b.value)
        for nid in nids:
            if nid not in node_id_to_index:
                continue
            gi = node_id_to_index[nid]
            for d in dofs:
                if 1 <= d <= 3:
                    out[3 * gi + (d - 1)] = value
    return out


# --------------------------------------------------------------------------- #
# Real comparison against an Abaqus ODB export (fields.json, CONTRACT sec 5)
# --------------------------------------------------------------------------- #
def _pick_frame(fields, which="last"):
    frames = fields["frames"]
    if not frames:
        raise ValueError("fields.json contains no frames")
    if which in ("last", None):
        return frames[-1]
    if which == "first":
        return frames[0]
    fid = int(which)
    for fr in frames:
        if int(fr.get("frame", -1)) == fid:
            return fr
    raise ValueError("frame %s not found; available: %s"
                     % (which, [f.get("frame") for f in frames]))


def residual_from_fields(model, fields, which="last", mode="finite"):
    """Assemble f_int from an exported frame and compare to RF / equilibrium.

    Returns a dict of metrics and the assembled vectors.
    """
    node_ids, coords, connectivity = build_mesh(model)
    _, n2i = c3d8.assemble_global_internal_force(
        node_ids, coords, [], np.zeros(3 * len(node_ids)), {}, mode="small")

    frame = _pick_frame(fields, which)

    # displacements (dict node_id -> (3,))
    U = {int(k): np.asarray(v, float) for k, v in frame["U"].items()}
    for nid in node_ids:
        U.setdefault(nid, np.zeros(3))

    # integration-point Cauchy stress: dict eid -> (8,6)
    sigma_all = {}
    for k, v in frame["S"].items():
        sigma_all[int(k)] = np.asarray(v, float).reshape(-1, 6)

    # only keep elements we have stress for
    conn = [(eid, c) for (eid, c) in connectivity if eid in sigma_all]

    F_int, n2i = c3d8.assemble_global_internal_force(
        node_ids, coords, conn, U, sigma_all, mode=mode)

    ndof = 3 * len(node_ids)

    # reactions RF (may be sparse -- only nonzero at constrained/loaded nodes)
    RF = np.zeros(ndof)
    for k, v in frame.get("RF", {}).items():
        gi = n2i[int(k)]
        RF[3 * gi:3 * gi + 3] = np.asarray(v, float)

    # external concentrated loads, if any
    F_ext = np.zeros(ndof)
    for cl in model.cloads:
        for nid in _resolve_target_nodes(model, cl.target):
            if nid in n2i and 1 <= cl.dof <= 3:
                F_ext[3 * n2i[nid] + (cl.dof - 1)] += float(cl.value)

    pres = prescribed_dofs(model, n2i)
    pres_idx = np.array(sorted(pres.keys()), dtype=int)
    free_mask = np.ones(ndof, dtype=bool)
    free_mask[pres_idx] = False

    # free-dof residual r = f_int - f_ext  (should be ~0 at converged increment)
    r = F_int - F_ext
    r_free = r[free_mask]

    # prescribed-dof: f_int vs reaction (try both signs)
    fi_p = F_int[pres_idx]
    rf_p = RF[pres_idx]
    err_plus = np.linalg.norm(fi_p - rf_p)
    err_minus = np.linalg.norm(fi_p + rf_p)
    sign = "+" if err_plus <= err_minus else "-"
    rf_err = min(err_plus, err_minus)
    rf_scale = np.linalg.norm(rf_p) + 1e-30

    metrics = {
        "frame": frame.get("frame"),
        "time": frame.get("time"),
        "n_dof": ndof,
        "n_prescribed": int(pres_idx.size),
        "n_free": int(free_mask.sum()),
        "free_residual_Linf": float(np.max(np.abs(r_free))) if r_free.size else 0.0,
        "free_residual_L2": float(np.linalg.norm(r_free)),
        "internal_force_scale": float(np.linalg.norm(F_int)),
        "reaction_sign_match": sign,
        "reaction_abs_err_L2": float(rf_err),
        "reaction_rel_err": float(rf_err / rf_scale),
        "sum_F_int": float(np.linalg.norm(F_int.reshape(-1, 3).sum(axis=0))),
    }
    return metrics, {"F_int": F_int, "RF": RF, "residual": r,
                     "free_mask": free_mask, "prescribed": pres}


# --------------------------------------------------------------------------- #
# Abaqus-independent manufactured-equilibrium self-test
# --------------------------------------------------------------------------- #
def manufactured_equilibrium(model, sigma0_voigt=None, verbose=True):
    """Impose a UNIFORM Cauchy stress on the mesh and check global equilibrium.

    For a spatially-constant stress field, the assembled internal force must
    vanish at every interior node and equal the surface traction at boundary
    nodes; the total must sum to zero.  This exercises parse -> assemble ->
    interior/boundary bookkeeping with no Abaqus dependence.
    """
    if sigma0_voigt is None:
        # a generic non-degenerate symmetric stress
        sigma0_voigt = np.array([120.0, -60.0, 30.0, 25.0, -15.0, 10.0])
    sigma0_voigt = np.asarray(sigma0_voigt, float)

    node_ids, coords, connectivity = build_mesh(model)
    nnode = len(node_ids)
    ndof = 3 * nnode
    U = np.zeros(ndof)                       # reference == current
    sigma_all = {eid: np.tile(sigma0_voigt, (8, 1)) for eid, _ in connectivity}

    F_int, n2i = c3d8.assemble_global_internal_force(
        node_ids, coords, connectivity, U, sigma_all, mode="small")

    # classify boundary vs interior nodes by face membership: a node is on the
    # mesh boundary iff it belongs to a facet not shared by two elements.
    boundary_nodes = _mesh_boundary_nodes(connectivity)
    interior_mask = np.ones(ndof, dtype=bool)
    for nid in boundary_nodes:
        gi = n2i[nid]
        interior_mask[3 * gi:3 * gi + 3] = False

    interior_res = F_int[interior_mask]
    total = F_int.reshape(-1, 3).sum(axis=0)

    fscale = np.linalg.norm(F_int) + 1e-30
    interior_Linf = float(np.max(np.abs(interior_res))) if interior_res.size else 0.0
    metrics = {
        "n_nodes": nnode,
        "n_elements": len(connectivity),
        "n_boundary_nodes": len(boundary_nodes),
        "n_interior_dofs": int(interior_mask.sum()),
        "interior_residual_Linf": interior_Linf,
        "interior_residual_rel": float(np.linalg.norm(interior_res) / fscale),
        "sum_F_int_Linf": float(np.max(np.abs(total))),
        "internal_force_scale": float(fscale),
    }
    ok = (metrics["interior_residual_rel"] < 1e-10
          and metrics["sum_F_int_Linf"] < 1e-6 * fscale + 1e-9)
    if verbose:
        print("Manufactured uniform-stress equilibrium test")
        print("  mesh: %d nodes, %d C3D8 elements, %d boundary nodes"
              % (nnode, len(connectivity), len(boundary_nodes)))
        print("  uniform sigma (Voigt) = %s" % (sigma0_voigt.tolist(),))
        print("  interior-node residual : Linf=%.3e  rel=%.3e   (must ~0)"
              % (interior_Linf, metrics["interior_residual_rel"]))
        print("  sum of all nodal forces: Linf=%.3e   (must ~0)"
              % metrics["sum_F_int_Linf"])
        print("  internal-force scale   : %.3e" % fscale)
        print("  ==> %s" % ("PASS" if ok else "FAIL"))
    return ok, metrics


def manufactured_linear(model, verbose=True):
    """Stronger global test: impose a LINEAR (spatially varying) Cauchy stress.

    A uniform stress makes each element's internal force self-equilibrated, so
    the uniform test is blind to per-IP weights and to intra-element force
    redistribution (audit 2026-07-10). A linear stress removes that blindness:
    at an INTERIOR node a the assembled force must equal the consistent nodal
    force of the constant body load b = div(sigma):

        F_a = int grad(N_a).sigma dV = - int N_a (div sigma) dV = - b * V_a ,
        V_a = int N_a dV   (interior node: the surface term vanishes).

    The target -b*V_a is computed independently of the Gauss stress pairing, so a
    match validates the global scatter + per-IP weights + detJ for a non-uniform
    field. Compression111 is all affine cubes -> exact quadrature -> machine
    precision. (Still does NOT probe the Abaqus IP *ordering*; that needs an ODB.)
    """
    rng = np.random.default_rng(2024)
    S0 = rng.uniform(-50, 50, (3, 3)); S0 = 0.5 * (S0 + S0.T)
    grads = [0.5 * (G + G.T) for G in rng.uniform(-8, 8, (3, 3, 3))]

    def sigma_vec(x):
        M = S0 + grads[0] * x[0] + grads[1] * x[1] + grads[2] * x[2]
        return np.array([M[0, 0], M[1, 1], M[2, 2], M[0, 1], M[0, 2], M[1, 2]])

    # constant divergence b_i = d sigma_ij / dx_j  (central diff, exact for linear)
    h = 1e-3
    b = np.zeros(3)
    for j in range(3):
        ep = np.zeros(3); ep[j] += h
        em = np.zeros(3); em[j] -= h
        b += (c3d8.voigt_to_tensor(sigma_vec(ep))[:, j]
              - c3d8.voigt_to_tensor(sigma_vec(em))[:, j]) / (2 * h)

    node_ids, coords, connectivity = build_mesh(model)
    n2i = {nid: i for i, nid in enumerate(sorted(node_ids))}
    id_to_row = {nid: i for i, nid in enumerate(node_ids)}
    pts, wts = c3d8.ABAQUS_C3D8_GAUSS.points, c3d8.ABAQUS_C3D8_GAUSS.weights

    sigma_all = {}
    V = np.zeros(len(node_ids))          # int N_a dV per node
    for eid, conn in connectivity:
        Xe = np.array([coords[id_to_row[nid]] for nid in conn], dtype=float)
        sig = np.zeros((8, 6))
        for k in range(8):
            N = c3d8.shape_functions(pts[k])
            sig[k] = sigma_vec(N @ Xe)
            _, detJ = c3d8.b_matrix_reference(Xe, pts[k])
            for a, nid in enumerate(conn):
                V[n2i[nid]] += N[a] * detJ * wts[k]
        sigma_all[eid] = sig

    U = np.zeros(3 * len(node_ids))
    F_int, _ = c3d8.assemble_global_internal_force(
        node_ids, coords, connectivity, U, sigma_all, mode="small")

    boundary = _mesh_boundary_nodes(connectivity)
    interior_ids = [nid for nid in node_ids if nid not in boundary]
    if not interior_ids:
        if verbose:
            print("Manufactured LINEAR-stress test: mesh has no interior nodes "
                  "(single-layer); skipping interior check.")
        return True, {"n_interior_nodes": 0, "skipped": True}

    fi = np.array([F_int[3 * n2i[nid]:3 * n2i[nid] + 3] for nid in interior_ids])
    tgt = np.array([-b * V[n2i[nid]] for nid in interior_ids])
    denom = max(np.linalg.norm(tgt), np.linalg.norm(fi), 1e-30)
    rel = float(np.linalg.norm(fi - tgt) / denom)
    ok = rel < 1e-9
    metrics = {"n_interior_nodes": len(interior_ids),
               "div_sigma": b.tolist(),
               "interior_force_vs_bodyload_rel": rel,
               "target_scale": float(np.linalg.norm(tgt))}
    if verbose:
        print("Manufactured LINEAR-stress equilibrium test")
        print("  interior nodes: %d   div(sigma)=[%.3f %.3f %.3f]"
              % (len(interior_ids), b[0], b[1], b[2]))
        print("  interior force vs -b*V_a (consistent body load): rel=%.3e (tol 1e-9)"
              % rel)
        print("  ==> %s" % ("PASS" if ok else "FAIL"))
    return ok, metrics


def _mesh_boundary_nodes(connectivity):
    """Nodes lying on the exterior surface of a C3D8 mesh.

    A quad facet on exactly one element is exterior; its 4 nodes are boundary
    nodes.  Facets are keyed by their frozenset of node ids (orientation-free).
    """
    faces = c3d8._C3D8_FACES
    count = {}
    facenodes = {}
    for eid, conn in connectivity:
        for f in faces:
            key = frozenset(conn[i] for i in f)
            count[key] = count.get(key, 0) + 1
            facenodes[key] = key
    boundary = set()
    for key, c in count.items():
        if c == 1:
            boundary |= set(key)
    return boundary


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _default_inp():
    return os.path.normpath(os.path.join(
        _HERE, "..", "sources", "permissive", "ngrilli_Oxford_Crystal_Plasticity",
        "ExampleInputFiles", "HCPnoTwin", "Compression111.inp"))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inp", default=None,
                    help="Abaqus .inp (default: Compression111.inp)")
    ap.add_argument("--fields", default=None,
                    help="fields.json exported from an Abaqus ODB (real Mode-1 check)")
    ap.add_argument("--frame", default="last", help="frame selector: last|first|<id>")
    ap.add_argument("--mode", default="finite", choices=["finite", "small"])
    ap.add_argument("--selftest", action="store_true",
                    help="run the Abaqus-independent manufactured-equilibrium tests")
    ap.add_argument("--report", default=None,
                    help="write the metrics dict as JSON to this path")
    args = ap.parse_args(argv)

    def _write_report(payload):
        if args.report:
            with open(args.report, "w") as fh:
                json.dump(payload, fh, indent=2)
            print("wrote report:", args.report)

    inp = args.inp or _default_inp()
    print("Parsing:", inp)
    model = parser.parse_inp(inp)
    ne = sum(1 for e in model.elements.values() if e.etype.upper() == "C3D8")
    print("  %d nodes, %d C3D8 elements, %d materials"
          % (len(model.nodes), ne, len(model.materials)))

    if args.fields:
        with open(args.fields) as fh:
            fields = json.load(fh)
        metrics, _ = residual_from_fields(model, fields, args.frame, args.mode)
        print("\nStress-driven residual vs Abaqus (frame %s, t=%s):"
              % (metrics["frame"], metrics["time"]))
        for k in ("n_dof", "n_free", "n_prescribed",
                  "free_residual_Linf", "free_residual_L2",
                  "internal_force_scale", "reaction_sign_match",
                  "reaction_abs_err_L2", "reaction_rel_err", "sum_F_int"):
            print("  %-22s : %s" % (k, metrics[k]))
        rel = metrics["reaction_rel_err"]
        freerel = metrics["free_residual_Linf"] / (metrics["internal_force_scale"] + 1e-30)
        ok = rel < 1e-3 and freerel < 1e-3
        print("  ==> %s (reaction rel err %.2e, free-residual rel %.2e; tol 1e-3)"
              % ("PASS" if ok else "REVIEW", rel, freerel))
        metrics["verdict"] = "PASS" if ok else "REVIEW"
        _write_report(metrics)
        return 0 if ok else 1

    # default: Abaqus-independent manufactured self-tests (uniform + linear)
    ok_u, m_u = manufactured_equilibrium(model)
    print()
    ok_l, m_l = manufactured_linear(model)
    # also confirm the finite-strain global-assembly path matches small-strain at U=0
    ok_fs, m_fs = _check_finite_mode_consistency(model)
    all_ok = ok_u and ok_l and ok_fs
    print("\nOVERALL manufactured self-test: %s" % ("PASS" if all_ok else "FAIL"))
    _write_report({"uniform": m_u, "linear": m_l, "finite_mode_check": m_fs,
                   "verdict": "PASS" if all_ok else "FAIL"})
    return 0 if all_ok else 1


def _check_finite_mode_consistency(model, verbose=True):
    """Exercise the mode='finite' global-assembly branch (untested otherwise):
    at U=0 the current config equals the reference config, so finite and small
    assembly of the same stress field must be identical."""
    node_ids, coords, connectivity = build_mesh(model)
    rng = np.random.default_rng(7)
    s0 = rng.uniform(-40, 40, 6)
    sigma_all = {eid: np.tile(s0, (8, 1)) for eid, _ in connectivity}
    U = np.zeros(3 * len(node_ids))
    F_small, _ = c3d8.assemble_global_internal_force(
        node_ids, coords, connectivity, U, sigma_all, mode="small")
    F_fin, _ = c3d8.assemble_global_internal_force(
        node_ids, coords, connectivity, U, sigma_all, mode="finite")
    rel = float(np.linalg.norm(F_fin - F_small) / max(np.linalg.norm(F_small), 1e-30))
    ok = rel < 1e-12
    if verbose:
        print("Finite-mode global-assembly consistency (U=0: finite==small): "
              "rel=%.3e  ==> %s" % (rel, "PASS" if ok else "FAIL"))
    return ok, {"finite_vs_small_at_U0_rel": rel}


if __name__ == "__main__":
    raise SystemExit(main())
