"""Independent ORIGINAL and archived Abaqus checks for the bounded workflow."""

import copy
import json
from pathlib import Path

import numpy as np

from ..core.dof_manager import DofManager
from ..core.field_sensitivity import fields_from_statev, solve_field_sensitivities
from ..core.model import from_abaqus
from ..formulations.c3d8_kernel import isotropic_D
from ..io.abaqus_inp_parser import parse_inp
from ..io.derivative_fields import load_derivative_fields
from .connected import solve_history
from .kinematics import batch_operators, batch_strain, element_internal_force_bbar
from .record import ReplayRecord


STEPS = (1e-4, 3e-5, 1e-5)
TOLERANCE = 2e-6


def scaled_error(actual, expected):
    actual, expected = np.asarray(actual), np.asarray(expected)
    if not np.all(np.isfinite(actual)) or not np.all(np.isfinite(expected)):
        raise ValueError("nonfinite numerical verification data")
    return float(np.max(np.abs(actual - expected)) / max(np.max(np.abs(expected)), 1e-12))


def check_error(name, actual, expected, tolerance=TOLERANCE):
    error = scaled_error(actual, expected)
    if error >= tolerance:
        raise ValueError(f"{name}: scaled error {error:.6e} >= {tolerance:.6e}")
    return error


def branches(states):
    states = np.asarray(states)
    return np.diff(states, axis=0, prepend=np.zeros_like(states[:1])) > 1e-12


def verify_connected(record, material, result):
    """FD restarts ORIGINAL from virgin state, never the transformed evaluator."""
    base = copy.deepcopy(record.raw)
    for increment in base["increments"]:
        for field in ("u", "stress_ip", "state_ip"):
            increment.pop(field, None)
    base_record = ReplayRecord(base)
    model = base_record.build_model()
    manager = DofManager(model.nodes)
    element_ids = sorted(model.elements)
    indices = np.array([manager.element_dofs(model.elements[eid].connectivity, ("UX", "UY", "UZ"))
                        for eid in element_ids])
    operators, _ = batch_operators([model.coords_of(model.elements[eid].connectivity)
                                    for eid in element_ids], result["integration"])
    displacements = np.array([row["u"] for row in result["increments"]])
    strains = np.array([batch_strain(operators, displacement[indices]) for displacement in displacements])
    paths = np.diff(strains, axis=0, prepend=np.zeros_like(strains[:1]))
    times = [row["dt"] for row in base["increments"]]
    nominal_states = np.array([row["state"] for row in result["increments"]])
    nominal_branch = branches(nominal_states)
    props = np.array(base_record.props)
    original = solve_history(base_record, material, original=True)
    primal_error = check_error("ORIGINAL primal displacement", displacements,
                               [row["u"] for row in original["increments"]], 1e-9)
    rows = []
    for column, name in enumerate(material.params):
        previous = {}
        for relative in STEPS:
            step = relative * max(abs(props[column]), 1.)
            perturbed = []
            material_paths = []
            for sign in (1, -1):
                data = copy.deepcopy(base)
                data["material"]["props"][column] += sign * step
                solved = solve_history(ReplayRecord(data), material, original=True)
                if not np.array_equal(branches([row["state"] for row in solved["increments"]]), nominal_branch):
                    raise ValueError(f"{name}: whole-model FD crosses a material branch")
                perturbed.append(solved)
                fields = {"stress": np.empty_like(strains), "state": np.empty_like(nominal_states)}
                for element in range(len(element_ids)):
                    for point in range(8):
                        marched = material.march_regular(data["material"]["props"], paths[:, element, point], times)
                        fields["stress"][:, element, point] = [entry["stress"] for entry in marched]
                        fields["state"][:, element, point] = [entry["statev"] for entry in marched]
                if not np.array_equal(branches(fields["state"]), nominal_branch):
                    raise ValueError(f"{name}: fixed-path FD crosses a material branch")
                material_paths.append(fields)
            errors = {"parameter": name, "relative_step": relative, "absolute_step": float(step)}
            checks = [("u", "du_dp", "total_u"), ("stress", "dsigma_dp", "total_stress"),
                      ("state", "dstate_dp", "total_state")]
            references = {}
            for value, derivative, label in checks:
                reference = (np.array([row[value] for row in perturbed[0]["increments"]]) -
                             np.array([row[value] for row in perturbed[1]["increments"]])) / (2 * step)
                actual = np.array([row[derivative] for row in result["increments"]])[..., column]
                errors[label] = check_error(f"{name} {label} h={relative}", actual, reference)
                references[label] = reference
            for value, derivative in (("stress", "fixed_path_dsigma_dp"), ("state", "fixed_path_dstate_dp")):
                reference = (material_paths[0][value] - material_paths[1][value]) / (2 * step)
                actual = np.array([row[derivative] for row in result["increments"]])[..., column]
                label = "fixed_path_" + value
                errors[label] = check_error(f"{name} {label} h={relative}", actual, reference)
                references[label] = reference
            if relative == STEPS[-1]:
                for label, reference in references.items():
                    errors[label + "_plateau"] = check_error(f"{name} {label} FD plateau", previous[label], reference)
            previous = references
            rows.append(errors)
    actual = np.array([row["du_dp"] for row in result["increments"]])
    naive = np.array([row["fixed_path_local_solve"] for row in result["increments"]])
    return {"passed": True, "reference": "ORIGINAL whole paths and independently re-equilibrated ORIGINAL model",
            "tolerance": TOLERANCE, "steps": list(STEPS), "rows": rows,
            "primal_u_error": primal_error,
            "fixed_path_local_solve_is_not_du_dp": float(np.max(np.abs(actual - naive))),
            "branch_checks": "all increments and all IPs at every perturbation"}


def verify_abaqus_fixture(asset_directory):
    base = Path(asset_directory)
    fixtures = base / "fixtures"
    model = from_abaqus(parse_inp(str(base / "nonuniform_c3d8.inp")))
    manager = DofManager(model.nodes)
    fields = load_derivative_fields(str(fixtures / "derivative_fields.json"))
    if fields.metadata.get("provenance", {}).get("odb") != "nonuniform_c3d8.odb":
        raise ValueError("genuine nonuniform_c3d8 Abaqus export provenance is required")
    tangent, derivatives = fields_from_statev(fields.statev, fields.sdv_layout)
    result = solve_field_sensitivities(model, tangent, derivatives, integration="selective_reduced")

    def displacement(filename):
        values = json.loads((fixtures / filename).read_text())
        vector = np.zeros(manager.ndof)
        for node, components in values.items():
            vector[manager.node_dofs(int(node))[:3]] = components
        return vector

    errors = {}
    for name, upper, lower, step, tolerance in (
            ("E", "U_Ep.json", "U_Em.json", 840., 3e-5),
            ("nu", "U_nup.json", "U_num.json", .004, 1e-5)):
        reference = (displacement(upper) - displacement(lower)) / (2 * step)
        errors["du_d" + name] = check_error("Abaqus FD " + name, result.du_da(name), reference, tolerance)
    coordinates = model.coords_of(model.elements[1].connectivity)
    operators, _ = batch_operators([coordinates])
    stress = json.loads((fixtures / "S_ip.json").read_text())
    stress = np.array([stress[str(point)] for point in range(1, 9)])
    reconstructed = batch_strain(operators, displacement("U_base.json")[None])[0] @ isotropic_D(210000., .3).T
    errors["stress_ip"] = check_error("Abaqus IP stress", reconstructed, stress, 1e-4)
    external = np.zeros(24)
    for node, fx, fz in ((5, 15, 30), (6, 25, 10), (7, 5, 40), (8, 35, 20)):
        external[3 * (node - 1)] = fx
        external[3 * (node - 1) + 2] = fz
    internal = element_internal_force_bbar(coordinates, stress)
    errors["equilibrium"] = check_error("Abaqus equilibrium", internal[12:], external[12:], 1e-5)
    permutation = [0, 6, 2, 3, 4, 5, 1, 7]
    errors["permuted_ip_equilibrium"] = scaled_error(element_internal_force_bbar(coordinates, stress[permutation])[12:], external[12:])
    if errors["permuted_ip_equilibrium"] <= .1:
        raise ValueError("Abaqus fixture no longer discriminates incorrect integration-point order")
    return {"passed": True, "kind": "real_archived_Abaqus_2021_elastic_export",
            "provenance": fields.metadata["provenance"], "errors": errors,
            "limits": "elastic fixture validates FE assembly, not a J2 Abaqus run; float32 ODB limits FD accuracy"}