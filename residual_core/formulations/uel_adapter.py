"""UEL / direct-residual adapter (Mode 3).

Some elements do NOT expose a material update or an integration-point field: they
hand you the *element* right-hand side and tangent directly, the way an Abaqus
user element (UEL) subroutine does. This adapter wraps such a routine behind the
standard ``Formulation.eval_element(...)`` contract so the formulation-agnostic
assembler (core/assembler.py) can drive it unchanged, with NO knowledge of
continuum kinematics or materials. Nothing here imports a solid formulation or a
UMAT path -- the element math lives entirely inside the user-supplied callable.

The wrapped callable is Abaqus-UEL-like::

    uel_fn(element_id, element_type, coords, dofs, svars, props,
           time, dtime, fields, options)
        -> (RHS, AMATRX, SVARS_new, diagnostics)

``RHS``      : (ndof,)        element right-hand side
``AMATRX``   : (ndof, ndof)   element operator matrix (or None)
``SVARS_new``: solution-dependent state variables after the increment (or None)
``diagnostics``: dict (or None)

======================================================================
SIGN CONVENTION  (the critical part -- must-document-and-verify)
======================================================================

Abaqus defines its user-element output relative to the *negative* residual:

      AMATRX = dR/du            (the element stiffness / Jacobian)
      RHS    = -R = Fext - Fint (the negative residual)

i.e. Abaqus assembles the Newton system ``AMATRX . du = RHS`` and RHS is what is
left of the external load after subtracting the internal force.

The rest of THIS framework uses ``sign_convention = "residual"`` (see
formulations/base.py): ``eval_element`` returns the element's *positive*
``F_internal`` contribution to the global residual

      R(U) = F_internal(U, state) - F_external

and the assembler adds ``-F_external`` itself (assembler.py: ``R = R - F_ext``).
So the vector this adapter must return is ``+R`` (= +F_internal contribution),
NOT the Abaqus RHS. The conversion is therefore:

      element_residual = -RHS = -(-R) = +R          (default Abaqus case)
      element_tangent  =  AMATRX = dR/du

Because ``element_residual == R`` and ``AMATRX == dR/du``, the returned tangent is
exactly ``d(element_residual)/d(dofs)`` -- which is what the Level-4
finite-difference check expects, and what the self-test at the bottom proves
numerically.

Construction flag ``rhs_is_negative_residual`` (default True) documents and
selects this convention:

  * True  (Abaqus): RHS = -R          -> element_residual = -RHS
  * False (rare)  : RHS = +R directly -> element_residual = +RHS

In the ``False`` case the UEL is assumed to also report ``AMATRX`` as the tangent
of the *positive* residual it returns (``AMATRX = dR/du``), so ``element_tangent =
AMATRX`` holds unchanged and ``d(element_residual)/d(dofs) == element_tangent``
remains true. Only the residual sign is flag-dependent; the tangent is passed
through as-is in both cases.

Because the adapter normalises the output to the framework's ``+R`` convention,
the class attribute ``sign_convention`` is set to ``"residual"`` (NOT ``"rhs"``):
downstream code sees a standard residual formulation.
"""

from __future__ import annotations

from typing import Any, Callable, Dict, Optional, Tuple, Union

import numpy as np

try:
    from .base import Formulation
except ImportError:  # executed directly as a script (python uel_adapter.py)
    import os as _os
    import sys as _sys
    _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
    from base import Formulation  # type: ignore


class UelAdapter(Formulation):
    """Wrap an Abaqus-UEL-like callable as a standard residual Formulation.

    Parameters
    ----------
    uel_fn : callable or None
        The user element routine (signature documented in the module docstring).
        May be None at construction (so the adapter registers as an available but
        unconfigured backend); ``eval_element`` then raises a clear error rather
        than silently returning zeros.
    dofs_per_node : int, optional
        DOFs per node for this user element (a UEL may have arbitrary DOFs per
        node). Default 3 (displacement DOFs). Ignored if ``dof_types`` is given.
    dof_types : tuple of str, optional
        Explicit per-node DOF labels; its length overrides ``dofs_per_node``.
    n_node : int or dict, optional
        Nodes per element. Either a single int (all element types) or a mapping
        ``{element_type: n_nodes}``. If None, ``n_nodes`` raises when queried.
    element_types : tuple of str, optional
        Abaqus element type keys this adapter handles (e.g. ("U1",)).
    rhs_is_negative_residual : bool, optional
        See the module docstring's SIGN CONVENTION section. Default True
        (Abaqus: RHS = -R).
    name : str, optional
        Registry key. Default "uel_direct".
    """

    name = "uel_direct"
    element_types: Tuple[str, ...] = ("U1",)
    dof_types: Tuple[str, ...] = ("UX", "UY", "UZ")
    # Normalised to the framework's +R convention by eval_element, so downstream
    # code treats this exactly like any other residual formulation.
    sign_convention = "residual"
    # A direct-residual element has no material/IP field to replay, so the
    # material-replay levels (6, and the material half of 7) do not apply.
    #   1 zero-field residual, 4 FD tangent, 5 solver comparison,
    #   7 full residual replay (assembly half).
    verification_levels: Tuple[int, ...] = (1, 4, 5, 7)
    supported_modes = ("direct-residual",)
    required_inputs = ("coords", "connectivity", "dofs", "uel_routine")
    optional_inputs = ("svars", "props")
    limitations = ("element math lives entirely in the user routine",
                   "DOFs/nodes declared by the user element",
                   "sign convention (RHS = -R) must be verified per routine")
    verification_tests = ("zero-field", "FD-tangent", "sign-convention")
    required_inputs_by_mode = {
        "direct-residual": ("coords", "connectivity", "dofs", "uel_routine")}
    optional_inputs_by_mode = {"direct-residual": ("svars", "props")}
    material_interface_needed = False
    state_requirements = "user-defined (SVARS)"
    tangent_support = "user-provided (AMATRX)"
    verification_status = "skeleton"

    def __init__(self,
                 uel_fn: Optional[Callable] = None,
                 dofs_per_node: int = 3,
                 dof_types: Optional[Tuple[str, ...]] = None,
                 n_node: Optional[Union[int, Dict[str, int]]] = None,
                 element_types: Optional[Tuple[str, ...]] = None,
                 rhs_is_negative_residual: bool = True,
                 name: Optional[str] = None):
        self.uel_fn = uel_fn
        self.rhs_is_negative_residual = bool(rhs_is_negative_residual)

        if dof_types is not None:
            self.dof_types = tuple(dof_types)
        else:
            n = int(dofs_per_node)
            if n < 1:
                raise ValueError("dofs_per_node must be >= 1")
            # Reuse the standard displacement labels for the common 3-DOF case;
            # otherwise synthesise generic labels DOF1..DOFn.
            if n == 3:
                self.dof_types = ("UX", "UY", "UZ")
            elif n == 2:
                self.dof_types = ("UX", "UY")
            elif n == 1:
                self.dof_types = ("UX",)
            else:
                self.dof_types = tuple("DOF%d" % (i + 1) for i in range(n))

        if element_types is not None:
            self.element_types = tuple(element_types)
        if name is not None:
            self.name = str(name)
        self._n_node = n_node

    # ------------------------------------------------------------------ #
    # base defines dofs_per_node as len(self.dof_types); setting the instance
    # attribute self.dof_types above makes that property return the configured
    # value, so no override is required.

    def n_nodes(self, element_type: str) -> int:
        """Nodes per element, from the configured value or per-type mapping."""
        n = self._n_node
        if isinstance(n, dict):
            if element_type not in n:
                raise KeyError("UelAdapter %r: no n_node configured for element "
                               "type %r (have %r)" % (self.name, element_type,
                                                      tuple(n.keys())))
            return int(n[element_type])
        if n is None:
            raise ValueError("UelAdapter %r: n_node not configured; pass "
                             "n_node=<int> or n_node={element_type: n} at "
                             "construction." % self.name)
        return int(n)

    # ------------------------------------------------------------------ #
    def eval_element(self, element_id, element_type, coords, dofs,
                     solution_state, material_state, properties,
                     time, dtime, fields, options):
        """Call the wrapped UEL routine and normalise it to the +R convention.

        Returns the standard assembler 4-tuple
        ``(element_residual, element_tangent, updated_state, diagnostics)`` where
        ``element_residual`` is the +F_internal contribution to the global
        residual (see the module SIGN CONVENTION section) and
        ``element_tangent`` is ``AMATRX = d(element_residual)/d(dofs)``.
        """
        if self.uel_fn is None:
            raise RuntimeError(
                "UelAdapter %r has no uel_fn: it was constructed without a user "
                "element routine and cannot evaluate an element. Construct it as "
                "UelAdapter(uel_fn=<callable>, dofs_per_node=..., n_node=...). "
                "(No silent zeros are returned.)" % self.name)

        coords = np.asarray(coords, dtype=float)
        u_e = np.asarray(dofs, dtype=float)
        # svars = solution-dependent state variables carried by the UEL.
        svars = material_state

        RHS, AMATRX, SVARS_new, uel_diag = self.uel_fn(
            element_id, element_type, coords, u_e, svars, properties,
            time, dtime, fields, options)

        RHS = np.asarray(RHS, dtype=float)
        # --- SIGN CONVERSION (see module docstring) ------------------------- #
        if self.rhs_is_negative_residual:
            element_residual = -RHS          # RHS = -R  ->  R = -RHS
        else:
            element_residual = np.array(RHS, dtype=float)   # RHS = +R already

        element_tangent = None if AMATRX is None else np.asarray(AMATRX, dtype=float)

        diagnostics: Dict[str, Any] = {
            "formulation": self.name,
            "rhs_is_negative_residual": self.rhs_is_negative_residual,
            "sign_convention": self.sign_convention,
            "ndof": int(element_residual.size),
        }
        if isinstance(uel_diag, dict):
            diagnostics.update(uel_diag)
        elif uel_diag is not None:
            diagnostics["uel_diagnostics"] = uel_diag

        return element_residual, element_tangent, SVARS_new, diagnostics


# ====================================================================== #
# SELF-TEST: trivial mock UEL for a 2-node, 1-DOF/node linear spring.
# ====================================================================== #
if __name__ == "__main__":
    import sys

    k_spring = 137.0
    Kmat = k_spring * np.array([[1.0, -1.0],
                                [-1.0, 1.0]], dtype=float)

    def mock_spring_uel(element_id, element_type, coords, dofs, svars, props,
                        time, dtime, fields, options):
        """R = k*[[1,-1],[-1,1]] @ u ; Abaqus: RHS = -R, AMATRX = dR/du = Kmat."""
        u = np.asarray(dofs, dtype=float)
        R = Kmat @ u
        RHS = -R                       # Abaqus convention
        AMATRX = Kmat.copy()           # dR/du
        return RHS, AMATRX, svars, {"k": k_spring}

    adapter = UelAdapter(uel_fn=mock_spring_uel, dofs_per_node=1, n_node=2,
                         element_types=("U1",))

    u_test = np.array([0.3, -0.2], dtype=float)
    R_true = Kmat @ u_test                     # +R  (internal force)
    RHS_true = -R_true                          # what the UEL returns

    res, tan, state, diag = adapter.eval_element(
        1, "U1", coords=np.array([[0.0], [1.0]]), dofs=u_test,
        solution_state={}, material_state=None, properties=None,
        time=(0.0, 0.0), dtime=0.0, fields={}, options={"compute_tangent": True})

    ok = True

    # (a) returned element_tangent == AMATRX
    err_tan = float(np.max(np.abs(tan - Kmat)))
    a_ok = err_tan == 0.0
    ok &= a_ok
    print("[a] element_tangent == AMATRX          : max|dK| = %.3e  -> %s"
          % (err_tan, "PASS" if a_ok else "FAIL"))

    # (b) returned element_residual == +R == -RHS  (numeric sign check, u != 0)
    err_res = float(np.max(np.abs(res - R_true)))
    err_vs_rhs = float(np.max(np.abs(res - (-RHS_true))))
    b_ok = (err_res < 1e-12) and (err_vs_rhs < 1e-12) and np.any(res != 0.0)
    ok &= b_ok
    print("[b] element_residual == +R == -RHS     : R=%s RHS=%s res=%s "
          "max|res-R|=%.3e max|res-(-RHS)|=%.3e -> %s"
          % (R_true.tolist(), RHS_true.tolist(), res.tolist(),
             err_res, err_vs_rhs, "PASS" if b_ok else "FAIL"))

    # (c) finite-difference d(residual)/d(dofs) matches returned tangent ~1e-6
    #     (proves the sign+tangent are mutually consistent).
    h = 1e-6
    n = u_test.size
    Kfd = np.zeros((n, n))
    for j in range(n):
        up = u_test.copy(); up[j] += h
        um = u_test.copy(); um[j] -= h
        rp = adapter.eval_element(1, "U1", np.array([[0.0], [1.0]]), up, {}, None,
                                  None, (0.0, 0.0), 0.0, {}, {})[0]
        rm = adapter.eval_element(1, "U1", np.array([[0.0], [1.0]]), um, {}, None,
                                  None, (0.0, 0.0), 0.0, {}, {})[0]
        Kfd[:, j] = (rp - rm) / (2.0 * h)
    err_fd = float(np.max(np.abs(Kfd - tan)))
    rel_fd = err_fd / max(float(np.max(np.abs(tan))), 1e-30)
    c_ok = rel_fd < 1e-6
    ok &= c_ok
    print("[c] FD d(residual)/d(dofs) == tangent  : max|Kfd-K|=%.3e rel=%.3e "
          "-> %s" % (err_fd, rel_fd, "PASS" if c_ok else "FAIL"))

    # Bonus: the unconfigured (no uel_fn) adapter must raise, never zeros.
    d_ok = False
    try:
        UelAdapter().eval_element(1, "U1", [[0.0], [1.0]], [0.0, 0.0], {}, None,
                                  None, (0.0, 0.0), 0.0, {}, {})
    except RuntimeError:
        d_ok = True
    ok &= d_ok
    print("[d] no-uel_fn eval_element raises       : -> %s"
          % ("PASS" if d_ok else "FAIL"))

    print("diagnostics:", diag)
    print("OVERALL:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)
