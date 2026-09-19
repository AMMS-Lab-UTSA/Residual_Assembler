"""Application service shared by the four-input CLI and GUI."""

import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import traceback

import numpy as np

from .connected import solve_history
from .j2_history import J2_SOURCE_HASH, require_j2
from .path_material import PathMaterial
from .presentation_inputs import convert_fields, read_model
from .request_reductions import reduce_scalar


PUBLIC_FILES = ("sensitivity_results.json", "sensitivity_tables.csv", "run_report.txt")


class RequestFailure(ValueError):
    """Public diagnostic containing only trusted category/action text."""

    def __init__(self, category, action):
        self.category = category
        self.action = action
        super().__init__("Category: %s\nAction: %s\nPrivate diagnostics: private/error_report.txt" % (category, action))


class NotThePinnedProvider(ValueError):
    """The completed mapping is valid and belongs to the object, but the object
    is another material than the pinned m3_j2 provider: the model is outside
    this engine's scope, not broken (``resasm request`` routes it to the
    history engine)."""


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def mapping_for(material, mapping=None):
    candidates = [Path(mapping)] if mapping else list(dict.fromkeys(
        candidate for candidate in (material.with_suffix(".json"), material.parent / "Mapping.json")
        if candidate.is_file()))
    if not candidates:
        raise ValueError("missing completed mapping: place the generated object-stem .json or Mapping.json beside OTI_UMAT.obj, or use --mapping")
    contracts = [(path, json.loads(path.read_text())) for path in candidates]
    if len(contracts) > 1 and any(contract != contracts[0][1] for _, contract in contracts[1:]):
        raise ValueError("ambiguous mapping sidecars; choose --mapping explicitly")
    path, contract = contracts[0]
    if contract.get("schema") != "resasm_umat_oti_contract_v1":
        raise ValueError("Mapping.json must be the generated completed resasm_umat_oti_contract_v1, not a transformation request")
    if contract.get("object", {}).get("sha256_full") != digest(material):
        raise ValueError("mapping object.sha256_full is missing or does not match OTI_UMAT.obj")
    if contract.get("regular_source_hash") != J2_SOURCE_HASH:
        raise NotThePinnedProvider(
            "mapping regular_source_hash does not match the pinned m3_j2 source: the material "
            "is %s, not the fingerprint-pinned m3_j2 J2 provider"
            % (contract.get("model_id") or "another provider"))
    if contract.get("layouts") != {
            "DSIGMA_DP": "fortran(NTENS,NPARAM)", "DSTATEV_DP": "fortran(NSTATV,NPARAM)",
            "DDSDDE": "fortran(NTENS,NTENS)", "voigt": ["11", "22", "33", "12", "13", "23"]}:
        raise ValueError("unsupported mapping derivative layouts or Voigt order")
    if [parameter.get("oti_direction") for parameter in contract.get("parameters", [])] != [1, 2, 3, 4]:
        raise ValueError("mapping parameter OTI directions must be 1,2,3,4 in generated order")
    return path, contract


def export_odb(odb, destination, abaqus="abaqus"):
    executable = shutil.which(str(abaqus))
    if executable is None:
        raise ValueError("Abaqus executable %r not found on PATH; licensed Abaqus Python with odbAccess is required to read Analysis.odb" % str(abaqus))
    exporter = Path(__file__).resolve().parents[1] / "io/abaqus_odb_export.py"
    command = [executable, "python", str(exporter), "--", "--odb", str(odb.resolve()),
               "--frames", "all", "--strict", "yes", "--out", str(destination.resolve())]
    completed = subprocess.run(command, capture_output=True, text=True)
    (destination.parent / "odb_export.log").write_text(completed.stdout + completed.stderr)
    if completed.returncode:
        raise ValueError("Abaqus ODB extraction failed (exit %s): %s\n%s" % (
            completed.returncode, " ".join(command), completed.stdout + completed.stderr))
    return json.loads(destination.read_text()), command


def validate_request(request, parameters):
    if not isinstance(request, dict) or set(request) != {"outputs", "parameters", "domain", "increments"}:
        raise ValueError("request requires exactly outputs, parameters, domain, increments")
    selected = request["parameters"]
    if not isinstance(selected, list) or not selected or any(name not in parameters for name in selected) or len(set(selected)) != len(selected):
        raise ValueError("parameters must be a nonempty unique list drawn from %s" % parameters)
    if not isinstance(request["domain"], dict) or set(request["domain"]) - {"nodes", "elements"}:
        raise ValueError("domain supports nodes and elements only")
    outputs = request["outputs"]
    if not isinstance(outputs, list) or not outputs:
        raise ValueError("outputs must be a nonempty list")
    names = []
    for output in outputs:
        if not isinstance(output, dict) or set(output) - {"name", "field", "component", "reduction", "domain"} or not {"name", "field", "component", "reduction"} <= set(output):
            raise ValueError("each output needs name, field, component, reduction; optional domain")
        if not isinstance(output["name"], str) or not output["name"]:
            raise ValueError("output name must be a nonempty string")
        names.append(output["name"])
        width = {"U": 3, "RF": 3, "S": 6, "SDV": 1}.get(output["field"])
        if width is None:
            raise ValueError("unsupported output field %r; use U, RF, S, SDV" % output["field"])
        component = output["component"]
        if component != "ALL" and (type(component) is not int or not 1 <= component <= width):
            raise ValueError("component must be ALL or a one-based field component")
        if output["reduction"] not in ("component", "sum", "mean", "L2", "max"):
            raise ValueError("unsupported reduction %r" % output["reduction"])
    if len(set(names)) != len(names):
        raise ValueError("output names must be unique")
    increments = request["increments"]
    if increments not in ("ALL", "LAST") and (not isinstance(increments, list) or not increments
            or any(type(number) is not int or number < 1 for number in increments)
            or len(set(increments)) != len(increments)):
        raise ValueError("increments must be ALL, LAST, or a unique list of one-based increment numbers")


def selected_ids(domain, key, available):
    if not isinstance(domain, dict) or set(domain) - {"nodes", "elements"} or key not in domain:
        raise ValueError("output domain must specify %s; no implicit node/element conversion" % key)
    selected = domain[key]
    if selected == "ALL":
        return available
    if not isinstance(selected, list) or not selected or any(type(value) is not int or value not in available for value in selected):
        raise ValueError("domain %s must be ALL or a nonempty list of existing ids" % key)
    if len(set(selected)) != len(selected):
        raise ValueError("duplicate domain ids")
    return selected


def scalar_results(result, request):
    validate_request(request, result["parameters"])
    count = len(result["increments"])
    selected = request["increments"]
    numbers = list(range(1, count + 1)) if selected == "ALL" else [count] if selected == "LAST" else selected
    if any(number > count for number in numbers):
        raise ValueError("requested increment exceeds saved history length %d" % count)
    node_ids = sorted(int(node) for node in result["record"]["mesh"]["nodes"])
    element_ids = result["element_ids"]
    for key in request["domain"]:
        selected_ids(request["domain"], key, node_ids if key == "nodes" else element_ids)
    columns = [result["parameters"].index(name) for name in request["parameters"]]
    rows = []
    times = np.cumsum([increment["dt"] for increment in result["record"]["increments"]])
    for number in numbers:
        increment = result["increments"][number - 1]
        for output in request["outputs"]:
            field = output["field"]
            key = "nodes" if field in ("U", "RF") else "elements"
            available = node_ids if key == "nodes" else element_ids
            domain = output.get("domain", request["domain"])
            for domain_key in domain:
                selected_ids(domain, domain_key, node_ids if domain_key == "nodes" else element_ids)
            ids = selected_ids(domain, key, available)
            positions = [available.index(identifier) for identifier in ids]
            if field in ("U", "RF"):
                values = np.asarray(increment["u" if field == "U" else "R"]).reshape(len(node_ids), 3)
                derivatives = np.asarray(increment["du_dp"]) if field == "U" else (
                    np.asarray(increment["Rp_history"]) + np.asarray(increment["K"]) @ np.asarray(increment["du_dp"]))
                derivatives = derivatives.reshape(len(node_ids), 3, len(result["parameters"]))
            else:
                values = np.asarray(increment["stress" if field == "S" else "state"])
                derivatives = np.asarray(increment["dsigma_dp" if field == "S" else "dstate_dp"])
            values, derivatives = values[positions], derivatives[positions][..., columns]
            component = output["component"]
            if component != "ALL":
                values, derivatives = values[..., component - 1], derivatives[..., component - 1, :]
            value, gradient = reduce_scalar(values, derivatives, output["reduction"])
            rows.append({"output": output["name"], "field": field, "component": component,
                         "reduction": output["reduction"], "domain": {key: ids}, "increment": number,
                         "time": float(times[number - 1]), "value": value,
                         "derivatives": dict(zip(request["parameters"], gradient.tolist()))})
    return rows


def run_request(*, model, odb, material, request, out, mapping=None, abaqus="abaqus", validate=False):
    model, odb, material, request, out = map(Path, (model, odb, material, request, out))
    if out.exists() and any(out.iterdir()):
        raise ValueError("output directory must be empty or new; refusing stale public results")
    private = out / "private"
    private.mkdir(parents=True)
    category, action = "input_path", "Check the four required input paths and file suffixes."
    try:
        for role, path, suffix in (("model", model, ".inp"), ("odb", odb, ".odb"),
                                   ("material", material, ".obj"), ("request", request, ".json")):
            action = "Check --%s: provide an existing readable %s file." % (role, suffix)
            if path.suffix.lower() != suffix or not path.is_file():
                raise ValueError("required %s input does not exist or has wrong suffix: %s" % (suffix, path))
        category, action = "material_mapping", "Check the completed mapping sidecar or --mapping and its object hash, source fingerprint and layouts."
        mapping_path, contract = mapping_for(material, mapping)
        category, action = "request_schema", "Check the request JSON: outputs, parameters, domain and increments must match the supported schema."
        specification = json.loads(request.read_text())
        validate_request(specification, [parameter["name"] for parameter in contract["parameters"]])
        category, action = "model_input", "Check --model for a supported single-step C3D8 INP with pinned J2 material, loads and boundaries."
        parsed = read_model(model)
        category, action = "material_provider", "Check the compiled material object, completed mapping and compiler/linker availability."
        provider = PathMaterial(str(material.resolve()), contract, str(private / "link"))
        require_j2(provider)
        category, action = "odb_export", "Check --odb and --abaqus; licensed Abaqus Python with odbAccess and a complete field history are required."
        fields, command = export_odb(odb, private / "fields.json", abaqus)
        category, action = "odb_fields", "Check ODB U, RF, CF, S and SDV1 fields, complete increment history and consistency with --model."
        record = convert_fields(parsed, fields, contract)
        category, action = "history_replay", "Check the saved history, material mapping and equilibrium consistency in the private diagnostics."
        result = solve_history(record, provider, replay=True, odb_tolerances=True)
        category, action = "reaction_mismatch", "Check ODB RF output and INP/ODB load, boundary and equilibrium consistency."
        for number, increment in enumerate(result["increments"], 1):
            reactions = np.asarray([fields["frames"][number]["RF"][str(node)] for node in sorted(parsed.nodes)]).ravel()
            np.testing.assert_allclose(increment["R"], reactions, rtol=2e-5,
                                       atol=2e-5 * max(np.max(np.abs(reactions)), 1.),
                                       err_msg="INP/ODB reaction or equilibrium mismatch")
        category, action = "output_selection", "Check requested fields, components, reductions, domains and increments."
        rows = scalar_results(result, specification)
        verification = {"status": "not_run", "passed": None}
        if validate:
            category, action = "derivative_verification", "Inspect private diagnostics for the independent derivative verification failure."
            from .verification import verify_connected
            verification = verify_connected(record, provider, result)
        category, action = "output_write", "Check output directory permissions and available space; inspect private diagnostics."
        summary = {"schema": "resasm_presentation_results_v1", "status": "completed", "request": specification,
            "scope": {"history_increments_replayed": len(result["increments"]),
                      "output_increments": sorted(set(row["increment"] for row in rows)),
                      "element_type": "C3D8", "material": "fingerprint-pinned m3_j2", "integration": result["integration"]},
            "metadata": {"sensitivity_semantics": result["sensitivity_semantics"], "verified": bool(validate),
                "verification": verification, "max_equilibrium_error": max(row["equilibrium_error"] for row in result["increments"]),
                "input_sha256": {"Analysis.inp": digest(model), "Analysis.odb": digest(odb),
                                 "OTI_UMAT.obj": digest(material), "Mapping.json": digest(mapping_path),
                                 "sensitivity_request.json": digest(request)},
                "regular_source_hash": contract["regular_source_hash"],
                "checks": ["mesh ids/coordinates/connectivity", "complete virgin-to-final history", "time/ramped CF loads",
                           "prescribed displacements", "free equilibrium", "IP stress and state", "reactions"],
                "odb_tolerances": {"scaled_equilibrium": 1e-5, "stress_rtol": 2e-5,
                                   "stress_reaction_atol": "2e-5 * max(max_abs_recorded_field, 1)",
                                   "zero_bc_atol": "1e-12 * max_mesh_extent",
                                   "state_rtol": 2e-5, "state_atol": 1e-8},
                "limitations": "one small-strain static step; zero BCs; homogeneous pinned J2; no load/geometry/BC derivatives"},
            "results": rows}
        (private / "result.json").write_text(json.dumps(result, allow_nan=False) + "\n")
        (private / "export_command.json").write_text(json.dumps(command) + "\n")
        (out / PUBLIC_FILES[0]).write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
        with (out / PUBLIC_FILES[1]).open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["output", "field", "component", "reduction", "domain", "increment", "time", "value", "parameter", "derivative"])
            for row in rows:
                for parameter, derivative in row["derivatives"].items():
                    writer.writerow([row["output"], row["field"], row["component"], row["reduction"],
                                     json.dumps(row["domain"], sort_keys=True), row["increment"], row["time"], row["value"], parameter, derivative])
        (out / PUBLIC_FILES[2]).write_text(
            "Status: executed successfully\nSensitivity semantics: total equilibrated history\n"
            "Production analysis rerun: no\nMaterial source read/transformed: no\n"
            "Independent derivative verification: %s\nHistory increments replayed: %d\n"
            "Maximum scaled free residual: %.9g (limit 1e-5)\nChecks: %s\n"
            "Solver: constrained dense tangent solve; reactions Rp_history + K du/dp\n"
            "ODB stress/reaction: rtol 2e-5, atol 2e-5 * max(max_abs_recorded_field, 1); state: rtol 2e-5, atol 1e-8\n"
            "Zero BC roundoff: at most 1e-12 * max_mesh_extent, then enforced as zero\n"
            "Limits: %s\n" % ("passed" if validate else "NOT RUN (execution is not independent verification)",
                len(result["increments"]), summary["metadata"]["max_equilibrium_error"],
                "; ".join(summary["metadata"]["checks"]), summary["metadata"]["limitations"]))
        return summary
    except Exception as error:
        (private / "error_report.txt").write_text(traceback.format_exc())
        if category == "odb_fields":
            for field in ("U", "RF", "CF", "S", "SDV"):
                if str(error) in ("ODB %s requires every node" % field,
                                  "ODB %s requires every element" % field,
                                  "ODB %s output is required to validate INP concentrated loads" % field) or (
                                      isinstance(error, KeyError) and error.args == (field,)):
                    category = "odb_missing_field"
                    action = "Request complete ODB %s output for every required entity and saved increment." % ("SDV1" if field == "SDV" else field)
                    break
        failure = RequestFailure(category, action)
        (out / PUBLIC_FILES[2]).write_text("Status: failed\nIndependent derivative verification: not established\n" + str(failure) + "\n")
        raise failure from error