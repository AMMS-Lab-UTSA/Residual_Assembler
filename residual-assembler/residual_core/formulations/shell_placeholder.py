"""Shell backend placeholder — declares support, refuses to fake results.

Registered so the model inspector can *recognise* shell elements (S3/S4/…) and
report them as "contract defined, backend not implemented yet" — instead of
silently skipping them or pretending a residual exists. ``eval_element`` raises a
clear ``NotImplementedError`` listing exactly what a real shell must provide.

This is the honest way to show the framework can host shells without claiming the
implementation is finished.
"""

from __future__ import annotations

from .shell_base import ShellFormulation, SHELL_CONTRACT


class ShellPlaceholder(ShellFormulation):
    name = "shell_placeholder"
    element_types = ("S3", "S4", "S4R", "S8R", "STRI3", "SHELL")
    dof_types = ("UX", "UY", "UZ", "RX", "RY", "RZ")
    sign_convention = "residual"
    verification_levels = ()
    supported_modes = ()               # nothing runnable yet
    required_inputs = SHELL_CONTRACT
    limitations = ("NOT IMPLEMENTED — placeholder only",
                   "recognises shell elements but cannot assemble them yet")
    verification_tests = ()
    notes = ("Contract defined in formulations/shell_base.py; a real shell must "
             "implement membrane strains, bending curvatures, through-thickness "
             "integration and section resultants.")
    required_inputs_by_mode = {}
    optional_inputs_by_mode = {}
    material_interface_needed = False
    state_requirements = "unknown (not implemented)"
    tangent_support = "none"
    verification_status = "contract-only"

    def n_nodes(self, element_type):
        et = (element_type or "").upper()
        if et in ("S3", "STRI3"):
            return 3
        if et in ("S8R",):
            return 8
        return 4

    def membrane_strain(self, coords, dofs, xi):
        raise NotImplementedError("shell_placeholder: no membrane kinematics yet")

    def bending_curvature(self, coords, dofs, xi):
        raise NotImplementedError("shell_placeholder: no bending kinematics yet")

    def eval_element(self, element_id, element_type, coords, dofs,
                     solution_state, material_state, properties,
                     time, dtime, fields, options):
        raise NotImplementedError(
            "Shell backend is a contract placeholder, not implemented. A real "
            "shell must provide: %s. See formulations/shell_base.py and "
            "docs/adding_a_formulation.md." % ", ".join(SHELL_CONTRACT))
