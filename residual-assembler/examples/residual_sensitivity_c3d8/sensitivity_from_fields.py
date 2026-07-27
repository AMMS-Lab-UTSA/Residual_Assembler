# -*- coding: utf-8 -*-
"""The residual assembler consuming a material formulation's OUTPUTS (thin caller).

"Our new thing": form the residual and its parameter sensitivities from the
fields a user's material formulation produces, without re-running the material.
An OTI-overloaded UMAT, seeded with the strain AND the parameter, writes per IP:

    S     : integration-point stress  sigma
    SDV   : DDSDDE = d sigma / d epsilon  (consistent tangent)  -> builds K
          : d sigma / d a  per parameter  (parameter sensitivity) -> builds R^(p)

This module is now a THIN CALLER. All of the assembly + solve logic
(assemble K from per-IP DDSDDE, assemble R^(p) from d sigma/d a, partition the
DOFs, solve K U^(1) = -R^(1) for every parameter at once) lives in the reusable
core engine ``residual_core.core.field_sensitivity.solve_field_sensitivities``,
which is multi-element and multi-parameter. Here we only:

  1. synthesize the per-IP fields an OTI-overloaded UMAT would export (for a real
     UMAT these come from the ODB; the assembler code is identical),
  2. describe the SDV layout as metadata and hand the fields to the core engine,
  3. chain-rule the displacement sensitivity to a stress output for display.

Verified against the analytic single-cube engine in ``sensitivity_engine.py``.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from residual_core.core.model import Element, Model
from residual_core.core.field_sensitivity import (
    solve_field_sensitivities, fields_from_statev)
from residual_core.formulations import c3d8_kernel as k
from examples.residual_sensitivity_c3d8 import sensitivity_engine as eng

# --- SDV layout (metadata, not hard-coded slices in the engine) ------------- #
#   Abaqus 1-based inclusive SDV numbers:
#     SDV1..36  : DDSDDE, 6x6 row-major (Abaqus Voigt 11,22,33,12,13,23)
#     SDV37..42 : d sigma / d a   for one parameter, Voigt
SDV_LAYOUT = {"indexing": "abaqus_1_based", "ddsdde": [1, 36],
              "parameters": {"E": [37, 42]}}


@dataclass
class _BC:
    """Duck-typed *Boundary object understood by core.constraints."""
    target: int
    dof_start: int = 0
    dof_end: int = 0
    value: float = 0.0
    kind: str = "value"
    amplitude: Optional[str] = None


def _one_cube_model(fixed):
    """A one-element C3D8 model on eng.XE, with `fixed` (0-based global dof ->
    value) turned into single-component Dirichlet BCs. Node ids are 1..8 so the
    core's node-major global dof(node-1, comp) matches eng.dof(node0, comp)."""
    nodes = {i + 1: tuple(eng.XE[i]) for i in range(8)}
    model = Model(nodes=nodes,
                  elements={1: Element(1, "C3D8", [1, 2, 3, 4, 5, 6, 7, 8])})
    model.element_formulation = {1: "solid_c3d8_small_strain"}
    model.element_material = {1: "MAT"}
    model.boundaries = [_BC(target=d // 3 + 1, dof_start=d % 3 + 1,
                            dof_end=d % 3 + 1, value=v) for d, v in fixed.items()]
    return model


def pack_statev(D, dsig_da):
    """Pack the per-IP SDV vector an overloaded UMAT would write (DDSDDE | dσ/da)."""
    return list(np.asarray(D, float).ravel()) + list(np.asarray(dsig_da, float))


def synthetic_overloaded_umat_fields(fixed, fext, E=210000.0):
    """Stand-in for an OTI-overloaded UMAT's ODB export: solve the elastic state
    and write DDSDDE and d sigma/dE into the SDV slots exactly as a transformed
    UMAT would. (For a real UMAT these come from Abaqus; the code below is
    identical.)"""
    u, sig, _K, _free = eng.solve(E, fixed, fext)
    D = k.isotropic_D(E, eng.NU)
    statev = [pack_statev(D, sig[i] / E) for i in range(8)]   # d sigma/dE = sigma/E
    disp = {str(n + 1): [u[eng.dof(n, c)] for c in range(3)] for n in range(8)}
    return {"stress_ip": {"1": sig.tolist()}, "statev": {"1": statev},
            "displacements": disp}


def assemble_sensitivity(fields, fixed):
    """Thin caller: parse the SDV layout, delegate assembly+solve to the core
    engine, then chain-rule to a stress output.

    fields : dict like the ODB export, with fields['statev']['1'] = per-IP SDV
             vectors packed as pack_statev(DDSDDE, d sigma/da).
    fixed  : {global_dof (0-based) -> value} boundary conditions.
    """
    model = _one_cube_model(fixed)
    statev = {1: fields["statev"]["1"]}
    tangent_fields, stress_fields = fields_from_statev(
        statev, SDV_LAYOUT, parameters=["E"])

    res = solve_field_sensitivities(
        model=model, tangent_fields=tangent_fields,
        stress_derivative_fields=stress_fields, parameters=["E"])
    U1 = res.du_da("E")                                  # (24,) du/dE

    # chain rule to a stress output: total d sigma / dE = d sigma/dE|_u + D B (du/dE)
    D_ip = tangent_fields[1]
    dsig_da = stress_fields["E"][1]
    B = [k.b_matrix_reference(eng.XE, xi)[0] for xi in eng.G.points]
    dsig_tot = np.array([dsig_da[i] + D_ip[i] @ (B[i] @ U1) for i in range(8)])
    return dict(dudp=U1, dsig_dp_total=dsig_tot)


def _verify():
    E0 = 210000.0
    # Case A: load control
    fixed = eng.symmetry_bcs()
    fext = np.zeros(24)
    for n in eng.X1:
        fext[eng.dof(n, 0)] = 100.0 / 4
    fields = synthetic_overloaded_umat_fields(fixed, fext, E0)
    got = assemble_sensitivity(fields, fixed)
    ref = eng.sensitivity_dE(E0, fixed, fext)          # the direct analytic engine
    d = eng.dof(1, 0)
    rel = abs(got["dudp"][d] - ref["dudE"][d]) / abs(ref["dudE"][d])
    print("field-consuming assembler (core engine) vs direct engine (Case A):")
    print("  du/dE (loaded node): from-fields %+.6e   engine %+.6e   rel %.2e"
          % (got["dudp"][d], ref["dudE"][d], rel))

    # Case B: displacement control (stress output)
    fixed = eng.symmetry_bcs()
    for n in eng.X1:
        fixed[eng.dof(n, 0)] = 1.0e-3
    fields = synthetic_overloaded_umat_fields(fixed, np.zeros(24), E0)
    got = assemble_sensitivity(fields, fixed)
    ref = eng.sensitivity_dE(E0, fixed, np.zeros(24))
    rel = abs(got["dsig_dp_total"][0][0] - ref["dsig_dE_total"][0][0]) \
        / abs(ref["dsig_dE_total"][0][0])
    print("  dSigma11/dE (IP0):  from-fields %+.6e   engine %+.6e   rel %.2e"
          % (got["dsig_dp_total"][0][0], ref["dsig_dE_total"][0][0], rel))
    print("  => the core engine reproduces the sensitivities purely from the "
          "exported SDV fields (DDSDDE -> K, d sigma/dE -> R^(p)).")


if __name__ == "__main__":
    _verify()
