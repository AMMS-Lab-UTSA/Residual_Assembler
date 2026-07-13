"""Shell formulation contract (base) — the interface a real shell must satisfy.

Production shell elements are intentionally NOT implemented yet. This base makes
the *requirements* explicit so the framework can advertise shell support as
"contract defined, implementation pending" rather than pretending it is done.

A concrete shell backend must provide, on top of the generic Formulation
contract, the shell-specific pieces enumerated in ``SHELL_CONTRACT`` below:

    translational DOFs (UX, UY, UZ)
    rotational DOFs    (RX, RY, RZ)     (drilling DOF handling if applicable)
    mid-surface interpolation           (shape functions on the reference surface)
    through-thickness integration       (layered / Gauss points across thickness)
    membrane strains                    (in-plane stretching)
    bending curvatures                  (change of surface curvature)
    transverse shear treatment          (or Kirchhoff constraint if thin)
    section resultants or stress integration
                                        (N, M, Q resultants, or 3D stress * thickness)

See docs/limitations.md and the ``adding_a_formulation.md`` guide.
"""

from __future__ import annotations

from abc import abstractmethod

from .base import Formulation


SHELL_CONTRACT = (
    "translational DOFs (UX, UY, UZ)",
    "rotational DOFs (RX, RY, RZ)",
    "mid-surface interpolation (reference-surface shape functions)",
    "through-thickness integration (layered or Gauss points)",
    "membrane strains (in-plane stretching)",
    "bending curvatures",
    "transverse shear treatment (or Kirchhoff thin-shell constraint)",
    "section resultants (N, M, Q) or through-thickness stress integration",
    "drilling DOF handling if the element has one",
)


class ShellFormulation(Formulation):
    """Abstract shell base. A real shell subclasses this and implements the
    membrane/bending/shear kinematics and section integration."""

    dof_types = ("UX", "UY", "UZ", "RX", "RY", "RZ")
    supported_modes = ("formulation", "stress-driven")

    @abstractmethod
    def membrane_strain(self, coords, dofs, xi):
        """In-plane membrane strain at a surface point (concrete shell only)."""
        raise NotImplementedError

    @abstractmethod
    def bending_curvature(self, coords, dofs, xi):
        """Bending curvature at a surface point (concrete shell only)."""
        raise NotImplementedError
