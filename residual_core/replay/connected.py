"""Bounded nonlinear C3D8 J2 equilibrium and total-history sensitivity replay."""

import copy

import numpy as np

from ..core.constraints import partition
from ..core.dof_manager import DofManager
from ..core.field_sensitivity import solve_field_sensitivities
from .j2_history import real_array, require_j2, step_j2
from .kinematics import batch_element_force, batch_element_tangent, batch_operators, batch_strain
from .record import ReplayRecord


def solve_history(record, material, *, original=False, replay=False, odb_tolerances=False):
    """Solve/replay every increment from virgin state; no state commits in Newton.

    For sensitivities, the RHS includes total incoming history derivatives and
    -B*du_previous/dp, holding current u fixed. After the constrained solve the
    update is differentiated again with B*(du_current-du_previous)/dp. The
    independent ORIGINAL mode never evaluates the OTI update.
    """
    require_j2(material)
    if odb_tolerances and not replay:
        raise ValueError("ODB tolerances are only valid for recorded replay")
    raw = record.raw
    if raw.get("schema") != "resasm_replay_record_v1" or record.kinematics != "small_strain":
        raise ValueError("connected replay requires resasm_replay_record_v1 small_strain")
    unsupported = set(raw) - {"schema", "kinematics", "integration", "mesh", "material",
                             "boundaries", "loads", "increments", "provenance"}
    if unsupported:
        raise ValueError(f"unsupported record fields: {sorted(unsupported)}")
    if set(raw.get("loads", {})) - {"cload"}:
        raise ValueError("only parameter-independent concentrated loads (cload) are supported")
    if set(raw.get("material", {})) != {"props"}:
        raise ValueError("one homogeneous material.props array is required")
    if raw.get("provenance", {}).get("regular_source_hash") != material.contract["regular_source_hash"]:
        raise ValueError("record provenance.regular_source_hash must match the provider ORIGINAL source")
    props = real_array(raw["material"]["props"], (4,), "props")
    if props[0] <= 0 or not -1 < props[1] < .5 or props[2] <= 0 or props[3] < 0:
        raise ValueError("invalid J2 properties: E>0, -1<nu<0.5, SIGY0>0, H>=0 required")
    if not record.increments:
        raise ValueError("at least one increment is required")
    if not raw.get("boundaries"):
        raise ValueError("boundaries are required")
    for boundary in raw["boundaries"]:
        if set(boundary) - {"target", "kind", "dof", "value"}:
            raise ValueError("unsupported boundary metadata or parameter dependence")
        if boundary.get("kind", "value") not in ("value", "ENCASTRE"):
            raise ValueError("only value and ENCASTRE boundaries are supported")
        if boundary.get("value", 0) != 0:
            raise ValueError("bounded replay supports homogeneous, parameter-independent boundaries only")
        if boundary.get("kind", "value") == "value" and boundary.get("dof") not in (1, 2, 3):
            raise ValueError("boundary dof must be 1, 2, or 3")
    for load in raw.get("loads", {}).get("cload", []):
        if set(load) != {"node", "dof", "value"} or load["dof"] not in (1, 2, 3):
            raise ValueError("cload requires node, dof (1..3), value only")
    model = record.build_model()
    element_ids = sorted(model.elements)
    if not element_ids or any(model.elements[eid].etype != "C3D8" for eid in element_ids):
        raise ValueError("connected replay supports C3D8 elements only")
    if any(len(model.elements[eid].connectivity) != 8 for eid in element_ids):
        raise ValueError("C3D8 requires eight nodes per element")
    manager = DofManager(model.nodes)
    free, prescribed, _ = partition(model, manager)
    if not free.any():
        raise ValueError("no free displacement degrees of freedom")
    integration = raw.get("integration", "selective_reduced")
    operators, weights = batch_operators([
        model.coords_of(model.elements[eid].connectivity) for eid in element_ids], integration)
    indices = np.array([manager.element_dofs(model.elements[eid].connectivity, ("UX", "UY", "UZ"))
                        for eid in element_ids])
    force = real_array(record.external_load(manager), (manager.ndof,), "external force")
    shape = (len(element_ids), 8)
    stress, state = np.zeros((*shape, 6)), np.zeros((*shape, 1))
    dstress, dstate = np.zeros((*shape, 6, 4)), np.zeros((*shape, 1, 4))
    fixed_stress, fixed_state = dstress.copy(), dstate.copy()
    previous_u = np.zeros(manager.ndof)
    previous_du = np.zeros((manager.ndof, 4))
    results = []
    output_record = copy.deepcopy(raw)
    if not replay:
        output_record["provenance"]["kind"] = "synthetic_converged_fe"

    def assemble(tangent, values):
        stiffness = np.zeros((manager.ndof, manager.ndof))
        residual = np.zeros((manager.ndof, *values.shape[3:]))
        local_stiffness = batch_element_tangent(operators, tangent, weights)
        local_force = batch_element_force(operators, values, weights)
        for position, edofs in enumerate(indices):
            stiffness[np.ix_(edofs, edofs)] += local_stiffness[position]
            np.add.at(residual, edofs, local_force[position])
        return stiffness, residual

    for number, increment in enumerate(record.increments):
        if set(increment) - {"dt", "load_factor", "u", "stress_ip", "state_ip"}:
            raise ValueError(f"increment {number}: unsupported fields")
        dt = increment.get("dt")
        factor = increment.get("load_factor")
        if not isinstance(dt, (float, int)) or not np.isfinite(dt) or dt <= 0:
            raise ValueError(f"increment {number}: dt must be finite and positive")
        if not isinstance(factor, (float, int)) or not np.isfinite(factor):
            raise ValueError(f"increment {number}: finite load_factor is required")
        displacement = (real_array(increment.get("u"), (manager.ndof,), "recorded u").copy()
                        if replay else previous_u.copy())
        boundary_tolerance = (1e-12 * np.max(np.ptp(np.asarray(list(model.nodes.values())), axis=0))
                              if odb_tolerances else 0.)
        if np.any(np.abs(displacement[prescribed]) > boundary_tolerance):
            raise ValueError("recorded u violates prescribed boundaries")
        displacement[prescribed] = 0.
        external = force * factor
        def evaluate(trial_displacement):
            strain_increment = batch_strain(operators, trial_displacement[indices] - previous_u[indices])
            new_stress, new_state = np.empty_like(stress), np.empty_like(state)
            tangent = np.empty((*shape, 6, 6))
            for element in range(shape[0]):
                for point in range(8):
                    arguments = (props, stress[element, point], state[element, point])
                    if original:
                        updated = material._reg_step(*arguments, strain_increment[element, point], dt)
                    else:
                        updated = material._oti_step(*arguments, np.zeros((6, 4)), np.zeros((1, 4)),
                                                     strain_increment[element, point], dt)
                    new_stress[element, point], new_state[element, point], tangent[element, point] = updated[:3]
            stiffness, internal = assemble(tangent, new_stress)
            residual = internal - external
            return strain_increment, new_stress, new_state, tangent, stiffness, residual

        for iteration in range(40):
            strain_increment, new_stress, new_state, tangent, stiffness, residual = evaluate(displacement)
            equilibrium_error = np.max(np.abs(residual[free])) / max(np.max(np.abs(external)), 1.)
            if np.isfinite(equilibrium_error) and equilibrium_error < (1e-5 if odb_tolerances else 1e-14):
                break
            if replay:
                raise ValueError(f"increment {number}: recorded path is not equilibrated (scaled R={equilibrium_error:.3e})")
            correction = np.linalg.solve(stiffness[np.ix_(free, free)], -residual[free])
            baseline = np.linalg.norm(residual[free])
            for backtrack in range(20):
                scale = 0.5 ** backtrack
                candidate = displacement.copy()
                candidate[free] += scale * correction
                trial_residual = evaluate(candidate)[-1]
                if np.linalg.norm(trial_residual[free]) <= (1 - 1e-4 * scale) * baseline:
                    displacement = candidate
                    break
            else:
                raise ValueError(f"increment {number}: Newton line search failed (scaled R={equilibrium_error:.3e})")
        else:
            raise ValueError(f"increment {number}: Newton did not converge in 40 iterations")
        if "stress_ip" in increment:
            expected = np.array([increment["stress_ip"][str(eid)] for eid in element_ids])
            np.testing.assert_allclose(new_stress, expected,
                                       rtol=2e-5 if odb_tolerances else 2e-9,
                                       atol=2e-5 * max(np.max(np.abs(expected)), 1.) if odb_tolerances else 1e-9)
        if "state_ip" in increment:
            expected = np.array([increment["state_ip"][str(eid)] for eid in element_ids])
            np.testing.assert_allclose(new_state, expected,
                                       rtol=2e-5 if odb_tolerances else 2e-9,
                                       atol=1e-8 if odb_tolerances else 1e-12)
        row = {"u": displacement.tolist(), "stress": new_stress.tolist(), "state": new_state.tolist(),
               "K": stiffness.tolist(), "R": residual.tolist(), "equilibrium_error": float(equilibrium_error),
               "iterations": iteration + 1}
        if not original:
            incoming_strain = -np.einsum("eqai,eim->eqam", operators, previous_du[indices])
            partial = np.empty_like(dstress)
            next_fixed_stress, next_fixed_state = np.empty_like(dstress), np.empty_like(dstate)
            for element in range(shape[0]):
                for point in range(8):
                    common = (material, props, stress[element, point], state[element, point])
                    updated = step_j2(*common, dstress[element, point], dstate[element, point],
                                      strain_increment[element, point], incoming_strain[element, point], dt)
                    partial[element, point] = updated[3]
                    fixed = step_j2(*common, fixed_stress[element, point], fixed_state[element, point],
                                    strain_increment[element, point], np.zeros((6, 4)), dt)
                    next_fixed_stress[element, point], next_fixed_state[element, point] = fixed[3:]
            fields = {name: {eid: partial[position, :, :, column]
                             for position, eid in enumerate(element_ids)}
                      for column, name in enumerate(material.params)}
            sensitivity = solve_field_sensitivities(model,
                {eid: tangent[position] for position, eid in enumerate(element_ids)}, fields,
                solution=displacement, dof_manager=manager, integration=integration)
            derivative = sensitivity.displacement_sensitivities
            total_strain = np.einsum("eqai,eim->eqam", operators, (derivative - previous_du)[indices])
            for element in range(shape[0]):
                for point in range(8):
                    updated = step_j2(material, props, stress[element, point], state[element, point],
                                      dstress[element, point], dstate[element, point],
                                      strain_increment[element, point], total_strain[element, point], dt)
                    dstress[element, point], dstate[element, point] = updated[3:]
            fixed_rhs = assemble(tangent, next_fixed_stress)[1]
            naive = np.zeros_like(derivative)
            naive[free] = np.linalg.solve(stiffness[np.ix_(free, free)], -fixed_rhs[free])
            row.update(du_dp=derivative.tolist(), Rp_history=sensitivity.residual_derivatives.tolist(),
                       Rp_fixed_path=fixed_rhs.tolist(), fixed_path_local_solve=naive.tolist(),
                       dsigma_dp=dstress.tolist(), dstate_dp=dstate.tolist(),
                       fixed_path_dsigma_dp=next_fixed_stress.tolist(),
                       fixed_path_dstate_dp=next_fixed_state.tolist())
            previous_du = derivative.copy()
            fixed_stress, fixed_state = next_fixed_stress, next_fixed_state
        stress, state, previous_u = new_stress.copy(), new_state.copy(), displacement.copy()
        output_record["increments"][number].update(u=displacement.tolist(),
            stress_ip={str(eid): stress[position].tolist() for position, eid in enumerate(element_ids)},
            state_ip={str(eid): state[position].tolist() for position, eid in enumerate(element_ids)})
        results.append(row)
    return {"schema": "resasm_connected_j2_v1", "parameters": material.params,
            "sensitivity_semantics": "total_equilibrated_history" if not original else "original_primal_only",
            "provenance": output_record["provenance"], "integration": integration,
            "free_dofs": np.flatnonzero(free).tolist(), "element_ids": element_ids,
            "increments": results, "record": output_record}