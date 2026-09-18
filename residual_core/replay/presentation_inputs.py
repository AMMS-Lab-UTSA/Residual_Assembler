"""Strict INP/ODB exchange conversion for the bounded presentation interface."""

import numpy as np

from ..io.abaqus_inp_parser import _iter_blocks, parse_inp
from .j2_history import real_array
from .record import ReplayRecord


OPTIONS = {
    "heading": set(), "preprint": {"echo", "model", "history", "contact"},
    "node": set(), "element": {"type", "elset"},
    "nset": {"nset", "generate", "instance"}, "elset": {"elset", "generate", "instance"},
    "part": {"name"}, "end part": set(), "assembly": {"name"},
    "end assembly": set(), "instance": {"name", "part"}, "end instance": set(),
    "material": {"name"}, "user material": {"constants"}, "depvar": set(),
    "solid section": {"elset", "material"}, "boundary": set(), "cload": set(),
    "step": {"name", "nlgeom", "inc"}, "static": set(), "end step": set(),
    "output": {"field", "history", "frequency"},
    "node output": {"nset"}, "element output": {"elset", "directions"},
}


def read_model(path):
    blocks = list(_iter_blocks(path.read_text().splitlines()))
    counts = {}
    for keyword, spelling, options, data in blocks:
        counts[keyword] = counts.get(keyword, 0) + 1
        if keyword not in OPTIONS:
            raise ValueError("unsupported INP keyword *%s" % spelling)
        unknown = set(options) - OPTIONS[keyword]
        if unknown:
            raise ValueError("unsupported *%s options: %s" % (spelling, sorted(unknown)))
        if keyword == "solid section" and any(line.replace(",", "").strip() for line in data):
            raise ValueError("unsupported solid section data")
        if keyword == "depvar" and len(data) != 1:
            raise ValueError("only one physical SDV without custom labels is supported")
        if keyword in ("node", "cload", "boundary"):
            for line in data:
                tokens = [token.strip() for token in line.rstrip().rstrip(",").split(",")]
                allowed_lengths = {"node": (4,), "cload": (3,), "boundary": (2, 3, 4)}[keyword]
                if len(tokens) not in allowed_lengths or not all(tokens):
                    raise ValueError("malformed *%s data: %s" % (spelling, line))
                integer_columns = [0] if keyword == "node" else [1] if keyword == "cload" else [1, 2][:len(tokens) - 1]
                if keyword == "boundary" and not tokens[1].isdigit():
                    if len(tokens) != 2 or tokens[1].upper() not in ("ENCASTRE", "PINNED", "XSYMM", "YSYMM", "ZSYMM"):
                        raise ValueError("unsupported symbolic boundary or malformed DOF: %s" % line)
                    integer_columns = []
                if any(not tokens[column].isdigit() or int(tokens[column]) < 1 for column in integer_columns):
                    raise ValueError("node ids and DOFs must be positive integers: %s" % line)
    for keyword in ("step", "static", "end step", "material", "user material", "solid section", "depvar"):
        if counts.get(keyword) != 1:
            raise ValueError("exactly one *%s block is required" % keyword)
    if counts.get("part", 0) > 1 or counts.get("instance", 0) > 1:
        raise ValueError("only one untransformed part/instance is supported")
    model = parse_inp(str(path))
    if model.warnings:
        raise ValueError("INP parser diagnostics: " + "; ".join(model.warnings))
    step = model.steps[0]
    if step.nlgeom or step.procedure != "static" or not step.name:
        raise ValueError("one named *Static step with NLGEOM=NO is required")
    if not step.time_period or not np.isfinite(step.time_period) or step.time_period <= 0:
        raise ValueError("positive finite static time period is required")
    if not model.nodes or not model.elements or any(element.etype != "C3D8" for element in model.elements.values()):
        raise ValueError("nonempty C3D8 mesh is required")
    material = next(iter(model.materials.values()))
    if not material.user_material or material.depvar != 1 or len(material.constants) != 4:
        raise ValueError("pinned J2 requires four USER MATERIAL constants and one DEPVAR")
    if set(model.element_material) != set(model.elements) or any(
            name.upper() != material.name.upper() for name in model.element_material.values()):
        raise ValueError("every element must use the same J2 solid section")
    return model


def target_nodes(model, target):
    nodes = [target] if isinstance(target, int) else next(
        (members for name, members in model.node_sets.items() if name.upper() == target.upper()), [])
    if not nodes or any(node not in model.nodes for node in nodes):
        raise ValueError("unknown or empty node target %r" % target)
    return sorted(set(nodes))


def convert_fields(model, fields, contract):
    if fields.get("export_mode") != "strict_request":
        raise ValueError("complete strict_request ODB export is required")
    step = model.steps[0]
    if fields.get("step_names") != [step.name] or len(fields.get("instance_names", [])) != 1:
        raise ValueError("INP/ODB step or instance scope mismatch")
    if fields.get("element_type") != "C3D8":
        raise ValueError("ODB element type must be C3D8")
    if set(fields["nodes"]) != {str(node) for node in model.nodes}:
        raise ValueError("INP/ODB node ids mismatch")
    if set(fields["elements"]) != {str(element) for element in model.elements}:
        raise ValueError("INP/ODB element ids mismatch")
    node_ids = sorted(model.nodes)
    element_ids = sorted(model.elements)
    for node in node_ids:
        np.testing.assert_allclose(real_array(fields["nodes"][str(node)], (3,), "ODB coordinates"),
                                   model.nodes[node], rtol=2e-7, atol=1e-10,
                                   err_msg="INP/ODB coordinates mismatch")
    for element in element_ids:
        if fields["elements"][str(element)] != model.elements[element].connectivity:
            raise ValueError("INP/ODB connectivity mismatch for element %s" % element)
    boundaries = []
    for boundary in model.boundaries:
        if boundary.value != 0 or boundary.amplitude or boundary.op != "MOD":
            raise ValueError("only zero, parameter-independent boundaries are supported")
        end = 3 if boundary.kind == "ENCASTRE" else boundary.dof_end
        if boundary.dof_start < 1 or end > 3 or end < boundary.dof_start:
            raise ValueError("unsupported boundary DOFs")
        for node in target_nodes(model, boundary.target):
            for dof in range(boundary.dof_start, end + 1):
                boundaries.append({"target": node, "dof": dof, "value": 0.})
    loads = []
    force = np.zeros((len(node_ids), 3))
    for load in model.cloads:
        if load.step != 1 or load.amplitude or load.op != "MOD" or load.dof not in (1, 2, 3) or not np.isfinite(load.value):
            raise ValueError("only single-step ramped concentrated loads with DOFs 1..3 are supported")
        for node in target_nodes(model, load.target):
            loads.append({"node": node, "dof": load.dof, "value": load.value})
            force[node_ids.index(node), load.dof - 1] += load.value
    if not loads or not boundaries:
        raise ValueError("nonempty CLOAD and boundary conditions are required")
    frames = fields.get("frames", [])
    if len(frames) < 2 or fields.get("sdv_labels") != [1]:
        raise ValueError("ODB needs virgin frame and every converged frame with SDV1")
    increments = []
    previous_time = 0.
    for number, frame in enumerate(frames):
        if frame.get("increment") != number or frame.get("step") != step.name:
            raise ValueError("ODB must save every increment consecutively from virgin increment 0")
        time = frame["time"]
        if not np.isfinite(time) or (number == 0 and time != 0) or (number and time <= previous_time):
            raise ValueError("invalid or nonmonotone ODB time")
        for field in ("U", "RF"):
            if set(frame.get(field, {})) != set(fields["nodes"]):
                raise ValueError("ODB %s requires every node" % field)
        displacement = real_array([frame["U"][str(node)] for node in node_ids], (len(node_ids), 3), "ODB U")
        real_array([frame["RF"][str(node)] for node in node_ids], (len(node_ids), 3), "ODB RF")
        if "CF" not in frame or set(frame["CF"]) - set(fields["nodes"]):
            raise ValueError("ODB CF output is required to validate INP concentrated loads")
        concentrated = real_array([frame["CF"].get(str(node), [0., 0., 0.]) for node in node_ids],
                                  (len(node_ids), 3), "ODB CF")
        np.testing.assert_allclose(concentrated, force * time / step.time_period, rtol=2e-6, atol=1e-7,
                                   err_msg="INP/ODB loads or ramp time mismatch")
        for field, width in (("S", 6), ("SDV", 1)):
            if set(frame.get(field, {})) != set(fields["elements"]):
                raise ValueError("ODB %s requires every element" % field)
            real_array([frame[field][str(element)] for element in element_ids],
                       (len(element_ids), 8, width), "ODB " + field)
        if number == 0:
            if np.any(displacement) or any(np.any(frame[field][str(element)])
                    for field in ("S", "SDV") for element in element_ids):
                raise ValueError("nonzero initial ODB state is unsupported")
        else:
            increments.append({"dt": time - previous_time, "load_factor": time / step.time_period,
                               "u": displacement.ravel().tolist(), "stress_ip": frame["S"], "state_ip": frame["SDV"]})
        previous_time = time
    if not np.isclose(previous_time, step.time_period, rtol=1e-10, atol=1e-12):
        raise ValueError("ODB does not reach the INP step time period")
    return ReplayRecord({"schema": "resasm_replay_record_v1", "kinematics": "small_strain",
        "integration": "selective_reduced", "mesh": {
            "nodes": fields["nodes"], "elements": {str(element): {
                "type": "C3D8", "connectivity": model.elements[element].connectivity} for element in element_ids}},
        "material": {"props": next(iter(model.materials.values())).constants},
        "boundaries": boundaries, "loads": {"cload": loads}, "increments": increments,
        "provenance": {"kind": "abaqus_odb", "regular_source_hash": contract["regular_source_hash"],
                       "step": step.name, "instance": fields["instance_names"][0]}})