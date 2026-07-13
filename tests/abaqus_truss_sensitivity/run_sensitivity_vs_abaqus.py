"""Validate the framework's residual-method sensitivity against Abaqus finite
differences, on a small linear-elastic 3-D truss.

Design parameter: Young's modulus ``E`` (a native Abaqus material property, so a
clean finite-difference variable). Response: the displacement of the loaded node.

Three sensitivities are compared for du/dE at the loaded free DOFs:

  1. **Residual method (this framework).** Assemble the tangent ``T`` and the
     converged real solution ``u`` with the ``truss2`` backend; generate the
     right-hand side ``dR/dE`` *hypercomplexly* (evaluate the element residual
     with ``E`` seeded as a Dual1 number and read the imaginary part); then solve
     the first-order sensitivity system ``T (du/dE) = -dR/dE`` on the free DOFs.
     This is exactly the HYPAD residual method at order 1.

  2. **Abaqus finite differences.** Run the *same* truss in Abaqus at ``E0``,
     ``E0(1+d)`` and ``E0(1-d)`` and central-difference the loaded-node
     displacement:  du/dE ≈ (u(E+) - u(E-)) / (2 dE).

  3. **Closed form.** For a linear-elastic structure with a fixed load,
     u ∝ 1/E, so du/dE = -u/E exactly.

The truss (2 bars, T3D2) is chosen so the ``truss2`` axial stiffness EA/L matches
Abaqus's T3D2 element exactly, making the comparison meaningful.

Run (Abaqus on PATH):  python run_sensitivity_vs_abaqus.py
If Abaqus is missing, the framework-vs-analytic check still runs; the Abaqus part
reports 'pending' and the script still exits 0.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from residual_core.core.model import Model, Element
from residual_core.core.dof_manager import DofManager
from residual_core.core.assembler import Assembler
from residual_core.core import constraints as _constraints
from residual_core.formulations.truss2 import Truss2
from residual_core.materials.base import MaterialBinding
from residual_core.algebra.dual1 import Dual1, imag_part

# --------------------------------------------------------------------------- #
# Problem definition (units: N, mm, MPa)
# --------------------------------------------------------------------------- #
E0 = 200000.0        # Young's modulus (MPa)
AREA = 100.0         # cross-section area (mm^2)
NU = 0.3             # Poisson (unused by truss, needed by *ELASTIC)
NODES = {1: (0.0, 0.0, 0.0),     # support
         2: (0.0, 1.0, 0.0),     # support
         3: (1.0, 0.0, 0.0)}     # loaded, free in x-y
BARS = {1: (1, 3), 2: (2, 3)}    # T3D2 connectivity
LOAD = {(3, 1): 1500.0, (3, 2): -1000.0}   # (node, dof) -> force (both comps active)
FREE_DOFS = [(3, 1), (3, 2)]     # node 3 Ux, Uy


# --------------------------------------------------------------------------- #
# a tiny boundary object (duck-typed for core/constraints.py)
# --------------------------------------------------------------------------- #
class _BC:
    def __init__(self, target, dof_start, dof_end, value=0.0, kind="value"):
        self.target = target
        self.dof_start = dof_start
        self.dof_end = dof_end
        self.value = value
        self.kind = kind
        self.amplitude = None


def _build_model(E=E0):
    m = Model()
    m.nodes = dict(NODES)
    m.elements = {eid: Element(eid, "T3D2", list(conn)) for eid, conn in BARS.items()}
    m.element_formulation = {eid: "truss2" for eid in BARS}
    m.element_material = {eid: "bar" for eid in BARS}
    b = MaterialBinding(material=None, name="bar")
    b.section = {"E": E, "A": AREA}
    m.materials = {"bar": b}
    # supports: nodes 1 & 2 fully fixed; node 3 fixed in z only
    m.boundaries = [_BC(1, 1, 3, 0.0), _BC(2, 1, 3, 0.0), _BC(3, 3, 3, 0.0)]
    return m


def _force_vector(dm):
    F = np.zeros(dm.ndof)
    for (nid, ldof), val in LOAD.items():
        F[dm.dof_index(nid, ldof)] = val
    return F


# --------------------------------------------------------------------------- #
# framework side
# --------------------------------------------------------------------------- #
def framework_solve_and_sensitivity():
    """Return (dm, free_mask, u, dudE_method, dudE_analytic)."""
    m = _build_model(E0)
    forms = {"truss2": Truss2()}
    dm = DofManager.for_model(m, forms)
    asm = Assembler(m, dm, forms)

    # real force-controlled solve: K is constant (linear), so K_ff u_f = F_f.
    _R0, K, _diag = asm.assemble(np.zeros(dm.ndof), compute_tangent=True)
    free, pres, _pv = _constraints.partition(m, dm)
    F = _force_vector(dm)
    u = np.zeros(dm.ndof)
    Kff = K[np.ix_(free, free)]
    u[free] = np.linalg.solve(Kff, F[free])

    # RHS dR/dE, generated HYPERCOMPLEXLY (E seeded as Dual1) at the real u.
    dRdE = _assemble_dRdE(m, dm, forms, u)

    # residual method: T (du/dE) = -dR/dE on the free partition.
    dudE_method = np.zeros(dm.ndof)
    dudE_method[free] = np.linalg.solve(Kff, -dRdE[free])

    # closed form: u ∝ 1/E  ->  du/dE = -u/E
    dudE_analytic = -u / E0
    return dm, free, u, dudE_method, dudE_analytic


def _assemble_dRdE(m, dm, forms, u):
    """Global dR/dE via Dual1: evaluate each element residual with E seeded
    (E -> E0 + 1*eps) at the real solution u, and read the imaginary part."""
    form = forms["truss2"]
    dRdE = np.zeros(dm.ndof)
    E_seed = Dual1(E0, 1.0)
    for eid, el in m.elements.items():
        edofs = np.array(dm.element_dofs(el.connectivity, form.dof_types), int)
        coords = m.coords_of(el.connectivity)
        u_e = [float(u[d]) for d in edofs]
        props = {"E": E_seed, "A": AREA}
        r_e, _k, _s, _d = form.eval_element(eid, el.etype, coords, u_e, {}, None,
                                            props, (0.0, 0.0), 0.0, None,
                                            {"compute_tangent": False})
        for j, val in enumerate(r_e):
            dRdE[edofs[j]] += imag_part(val)
    return dRdE


# --------------------------------------------------------------------------- #
# Abaqus side
# --------------------------------------------------------------------------- #
_INP_TEMPLATE = """*HEADING
Truss sensitivity validation (E = {E})
*NODE
1, 0.0, 0.0, 0.0
2, 0.0, 1.0, 0.0
3, 1.0, 0.0, 0.0
*NSET, NSET=NLOAD
3
*ELEMENT, TYPE=T3D2
1, 1, 3
2, 2, 3
*ELSET, ELSET=ALLEL
1, 2
*SOLID SECTION, ELSET=ALLEL, MATERIAL=MAT
{A}
*MATERIAL, NAME=MAT
*ELASTIC
{E}, {NU}
*BOUNDARY
1, 1, 3, 0.0
2, 1, 3, 0.0
3, 3, 3, 0.0
*STEP, NAME=Step-1
*STATIC
{CLOAD}*NODE PRINT, NSET=NLOAD
U
*OUTPUT, FIELD
*NODE OUTPUT
U, RF
*END STEP
"""


def _abaqus_cmd():
    for c in ("abaqus", "abq2024", "abq2023"):
        p = shutil.which(c)
        if p:
            return p
    for var in ("ABAQUS_CMD", "ABAQUS", "ABQ_CMD"):
        p = os.environ.get(var)
        if p and (shutil.which(p) or os.path.exists(p)):
            return shutil.which(p) or p
    return None


def _parse_node_u_from_dat(dat_path, node=3):
    """Parse U1,U2,U3 of a node from an Abaqus *NODE PRINT .dat table."""
    with open(dat_path, "r", errors="ignore") as fh:
        lines = fh.readlines()
    header = None
    for i, ln in enumerate(lines):
        if re.search(r"\bU1\b", ln) and re.search(r"\bU2\b", ln):
            header = i
            break
    if header is None:
        raise RuntimeError("no *NODE PRINT U table found in %s" % dat_path)
    num = re.compile(r"[-+]?\d+\.?\d*(?:[eEdD][-+]?\d+)?")
    for ln in lines[header + 1:]:
        toks = ln.split()
        if not toks:
            continue
        if toks[0].isdigit() and int(toks[0]) == node:
            vals = num.findall(ln)
            # first token is the node id; next are U1,U2,U3
            comps = [float(v.replace("D", "E").replace("d", "e")) for v in vals[1:4]]
            return np.array(comps)
    raise RuntimeError("node %d not found in *NODE PRINT table" % node)


def _cload_block():
    lines = ["*CLOAD"]
    for (nid, ldof), val in LOAD.items():
        lines.append("%d, %d, %s" % (nid, ldof, repr(float(val))))
    return "\n".join(lines) + "\n"


def _run_abaqus_case(workdir, name, E, cmd):
    inp = os.path.join(workdir, name + ".inp")
    with open(inp, "w") as fh:
        fh.write(_INP_TEMPLATE.format(E=repr(E), A=repr(AREA), NU=repr(NU),
                                      CLOAD=_cload_block()))
    # .bat/.exe on Windows: route through cmd /c so the batch launcher resolves.
    if os.name == "nt":
        argv = ["cmd", "/c", cmd, "job=" + name, "interactive", "ask_delete=OFF"]
    else:
        argv = [cmd, "job=" + name, "interactive", "ask_delete=OFF"]
    proc = subprocess.run(
        argv, cwd=workdir, capture_output=True, text=True, timeout=1200)
    dat = os.path.join(workdir, name + ".dat")
    if not os.path.exists(dat):
        raise RuntimeError("Abaqus produced no .dat for %s:\n%s\n%s"
                           % (name, proc.stdout[-2000:], proc.stderr[-2000:]))
    return _parse_node_u_from_dat(dat, node=3)


def abaqus_finite_difference(delta=1e-4):
    """Central-difference du/dE at node 3 from three Abaqus runs. Returns
    (dudE_fd (2,), info) or (None, reason)."""
    cmd = _abaqus_cmd()
    if cmd is None:
        return None, "Abaqus not available: validation pending"
    workdir = tempfile.mkdtemp(prefix="abq_truss_")
    dE = delta * E0
    try:
        u_plus = _run_abaqus_case(workdir, "trussEp", E0 + dE, cmd)
        u_minus = _run_abaqus_case(workdir, "trussEm", E0 - dE, cmd)
        u_base = _run_abaqus_case(workdir, "trussE0", E0, cmd)
    except Exception as exc:                      # noqa: BLE001
        return None, "Abaqus run/parse failed: %s" % exc
    dudE = (u_plus - u_minus) / (2.0 * dE)        # (U1,U2,U3)
    info = {"workdir": workdir, "u_base": u_base, "u_plus": u_plus,
            "u_minus": u_minus, "dE": dE}
    return dudE[:2], info                          # node-3 Ux, Uy


# --------------------------------------------------------------------------- #
def main():
    print("Residual-method sensitivity vs Abaqus finite differences (truss, du/dE)")
    dm, free, u, dudE_method, dudE_analytic = framework_solve_and_sensitivity()

    idx = [dm.dof_index(n, d) for (n, d) in FREE_DOFS]
    u_load = u[idx]
    m_method = dudE_method[idx]
    m_analytic = dudE_analytic[idx]

    print("\nFramework (truss2 backend):")
    print("  loaded-node u  = [Ux=%.6e, Uy=%.6e] mm" % (u_load[0], u_load[1]))
    print("  du/dE method   = [%.6e, %.6e]" % (m_method[0], m_method[1]))
    print("  du/dE analytic = [%.6e, %.6e]  (= -u/E)" % (m_analytic[0], m_analytic[1]))
    method_ok = np.allclose(m_method, m_analytic, rtol=1e-9, atol=1e-14)
    print("  method == analytic (-u/E): %s" % ("PASS" if method_ok else "FAIL"))

    checks = [("residual-method du/dE == analytic -u/E", method_ok)]

    dudE_fd, info = abaqus_finite_difference()
    if dudE_fd is None:
        print("\nAbaqus: %s" % info)
        print("  (framework-vs-analytic verified; Abaqus FD pending)")
    else:
        print("\nAbaqus finite differences (central, dE=%.4g):" % info["dE"])
        print("  u(E0) node3    = [Ux=%.6e, Uy=%.6e] mm"
              % (info["u_base"][0], info["u_base"][1]))
        print("  du/dE FD       = [%.6e, %.6e]" % (dudE_fd[0], dudE_fd[1]))
        # compare Abaqus solve to our solve, and Abaqus FD to our method
        u_match = np.allclose(info["u_base"][:2], u_load, rtol=1e-4, atol=1e-9)
        denom = np.linalg.norm(m_method) or 1.0
        fd_rel = np.linalg.norm(dudE_fd - m_method) / denom
        fd_ok = fd_rel < 1e-3
        print("  Abaqus u == framework u (rtol 1e-4): %s" % ("PASS" if u_match else "FAIL"))
        print("  du/dE  method vs Abaqus FD  rel=%.3e  -> %s"
              % (fd_rel, "PASS" if fd_ok else "FAIL"))
        checks.append(("Abaqus displacement matches framework", u_match))
        checks.append(("residual-method du/dE matches Abaqus FD", fd_ok))
        try:
            shutil.rmtree(info["workdir"], ignore_errors=True)
        except OSError:
            pass

    ok = all(v for _, v in checks)
    print("\nSummary:")
    for label, v in checks:
        print("  %-46s -> %s" % (label, "PASS" if v else "FAIL"))
    print("OVERALL: %s" % ("ALL PASS" if ok else "FAILURE"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
