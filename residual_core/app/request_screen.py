"""The parts of the Solve screen the presentation shows (slides 18 and 42).

"Point at the compiled OTI object from the previous build. Point at the saved
analysis .odb and its setup .inp. Tick the parameters, then choose the output
and the region. Click Solve. Then the residual equation is assembled and
solved, letting the user see the results of the sensitivity in the defined
domain."

The request screen in :mod:`residual_core.app.streamlit_app` already takes the
object, the .inp and the .odb and runs ``resasm request``. This module adds the
three things it took from a hand-written ``sensitivity_request.json`` instead:

* the parameter tick boxes, read from the ``Mapping.json`` the provider build
  wrote beside the object;
* the output and the region, the region offered from the sets in the .inp;
* the full-field result, read back from what the request wrote.

It builds the same ``sensitivity_request.json`` a user would write, so the run
is still exactly ``resasm request --request <that file>``: nothing is solved,
reduced or differentiated here, and nothing from the engine is imported (see
tests/framework/test_gui_is_a_thin_cli_front_end.py). The full-field table shows
the arrays the request itself stored in ``private/result.json`` -- ``u`` and
``du_dp`` at the nodes, ``stress``/``dsigma_dp`` and ``state``/``dstate_dp`` at
every integration point. Reaction derivatives are not stored per node there, so
for a reaction output the screen shows the reduced value the request computed.
The one repository module used is the .inp reader, to offer the .inp's own sets
in the region dropdown.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

#: What ``resasm request`` can return (``validate_request`` in
#: residual_core/replay/presentation.py): field, one-based component, label.
OUTPUTS: dict[str, tuple[str, int, str]] = {
    "displacement U1": ("U", 1, "U1"), "displacement U2": ("U", 2, "U2"),
    "displacement U3": ("U", 3, "U3"),
    "reaction force RF1": ("RF", 1, "RF1"), "reaction force RF2": ("RF", 2, "RF2"),
    "reaction force RF3": ("RF", 3, "RF3"),
    "stress S11": ("S", 1, "S11"), "stress S22": ("S", 2, "S22"), "stress S33": ("S", 3, "S33"),
    "stress S12": ("S", 4, "S12"), "stress S13": ("S", 5, "S13"), "stress S23": ("S", 6, "S23"),
    "state variable SDV1": ("SDV", 1, "SDV1"),
}
REDUCTIONS = ("mean", "sum", "max", "L2", "component")
WHOLE_MESH = "whole mesh"
INCREMENTS = {"last increment": "LAST", "every increment": "ALL"}


def mapping_path(material: str | None, mapping: str | None) -> Path | None:
    """The mapping ``resasm request`` will use: explicit, else beside the object.

    Discovery follows :func:`residual_core.replay.presentation.mapping_for`:
    ``<object-stem>.json`` or ``Mapping.json`` in the object's directory.
    """
    if mapping and Path(mapping).is_file():
        return Path(mapping)
    if material:
        obj = Path(material)
        for candidate in (obj.with_suffix(".json"), obj.parent / "Mapping.json"):
            if candidate.is_file():
                return candidate
    return None


def mapping_parameters(mapping: dict[str, Any]) -> list[str]:
    """Parameter names in the mapping's own order (its OTI directions)."""
    if mapping.get("schema") != "resasm_umat_oti_contract_v1":
        raise ValueError("Mapping.json is not a completed resasm_umat_oti_contract_v1")
    return [str(entry["name"]) for entry in mapping.get("parameters") or []]


def model_sets(inp_path: Path) -> dict[str, dict[str, list[int]]]:
    """Node and element sets of the .inp, by the repository's own INP parser."""
    from residual_core.io.abaqus_inp_parser import parse_inp

    model = parse_inp(str(inp_path))
    return {"nodes": {name: sorted(set(ids)) for name, ids in model.node_sets.items()},
            "elements": {name: sorted(set(ids)) for name, ids in model.element_sets.items()}}


def region_choices(sets: dict[str, dict[str, list[int]]], field: str) -> list[str]:
    kind = "nodes" if field in ("U", "RF") else "elements"
    return [WHOLE_MESH] + [f"{'node' if kind == 'nodes' else 'element'} set {name}"
                           for name in sorted(sets.get(kind) or {})]


def build_request(*, parameters: list[str], output: str, region: str, reduction: str,
                  increments: str, sets: dict[str, dict[str, list[int]]]) -> dict[str, Any]:
    """The ``sensitivity_request.json`` the ticks and choices stand for."""
    if not parameters:
        raise ValueError("tick at least one parameter")
    field, component, label = OUTPUTS[output]
    kind = "nodes" if field in ("U", "RF") else "elements"
    if region == WHOLE_MESH:
        ids: Any = "ALL"
        where = "all"
    else:
        name = region.split(" set ", 1)[1]
        ids = list((sets.get(kind) or {}).get(name) or [])
        if not ids:
            raise ValueError(f"{region} is empty in the .inp")
        where = name
    return {"outputs": [{"name": f"{reduction}_{label}_{where}", "field": field,
                         "component": component, "reduction": reduction,
                         "domain": {kind: ids}}],
            "parameters": list(parameters),
            "domain": {"nodes": "ALL", "elements": "ALL"},
            "increments": INCREMENTS[increments]}


def write_request(request: dict[str, Any], directory: Path) -> Path:
    path = Path(directory) / "sensitivity_request.json"
    path.write_text(json.dumps(request, indent=2) + "\n", encoding="utf-8")
    return path


def full_field(out_dir: Path, props_index: dict[str, int] | None = None) -> dict[str, Any]:
    """Every location's value and derivatives for the solved request.

    Rows are read from ``private/result.json``; nothing is re-solved. Returns
    the table rows plus two headline numbers: how well the replay reproduced
    the .odb stress (both arrays as the request stored them), and the size of
    the solution sensitivity du/dp for the ticked parameters. ``props_index``
    (from Mapping.json) locates each parameter's value among the .inp's USER
    MATERIAL constants, for the governing-parameter column; without it the
    column is left out.
    """
    out_dir = Path(out_dir)
    summary = json.loads((out_dir / "sensitivity_results.json").read_text())
    result = json.loads((out_dir / "private" / "result.json").read_text())
    fields = json.loads((out_dir / "private" / "fields.json").read_text())
    request = summary["request"]
    names = list(result["parameters"])
    ticked = list(request["parameters"])
    columns = [names.index(name) for name in ticked]
    props = np.asarray(result["record"]["material"]["props"], dtype=float)
    values_of = {name: float(props[index - 1]) for name, index in (props_index or {}).items()
                 if name in ticked and 0 < int(index) <= len(props)}
    node_ids = sorted(int(node) for node in result["record"]["mesh"]["nodes"])
    element_ids = list(result["element_ids"])
    rows: list[dict[str, Any]] = []
    for increment_number in summary["scope"]["output_increments"]:
        increment = result["increments"][increment_number - 1]
        for output in request["outputs"]:
            field, component = output["field"], output["component"]
            kind = "nodes" if field in ("U", "RF") else "elements"
            domain = output.get("domain", request["domain"])[kind]
            selected = (node_ids if kind == "nodes" else element_ids) if domain == "ALL" else list(domain)
            components = range(1, {"U": 3, "RF": 3, "S": 6, "SDV": 1}[field] + 1) \
                if component == "ALL" else [component]
            for comp in components:
                if field == "U":
                    values = np.asarray(increment["u"]).reshape(len(node_ids), 3)
                    derivatives = np.asarray(increment["du_dp"]).reshape(len(node_ids), 3, len(names))
                    for node in selected:
                        at = node_ids.index(node)
                        rows.append(_row(output["name"], increment_number, f"node {node}", comp,
                                         values[at, comp - 1], derivatives[at, comp - 1, columns], ticked))
                elif field == "RF":
                    continue  # not stored per node; the reduced result is shown above
                else:
                    key, dkey = ("stress", "dsigma_dp") if field == "S" else ("state", "dstate_dp")
                    values = np.asarray(increment[key])
                    derivatives = np.asarray(increment[dkey])
                    for element in selected:
                        at = element_ids.index(element)
                        for point in range(values.shape[1]):
                            rows.append(_row(output["name"], increment_number,
                                             f"element {element} IP {point + 1}", comp,
                                             values[at, point, comp - 1],
                                             derivatives[at, point, comp - 1, columns], ticked))
    if len(values_of) == len(ticked):
        for row in rows:
            scaled = {name: abs(row[f"d/d{name}"] * values_of[name]) for name in ticked}
            row["governing parameter"] = max(scaled, key=scaled.get) if max(scaled.values()) > 0 else "none"
    last = result["increments"][summary["scope"]["output_increments"][-1] - 1]
    frame = fields["frames"][summary["scope"]["output_increments"][-1]]
    replayed = np.asarray(last["stress"])
    recorded = np.asarray([frame["S"][str(element)] for element in element_ids], dtype=float)
    scale = max(float(np.max(np.abs(recorded))), 1.0)
    du_dp = np.asarray(last["du_dp"])[:, columns]
    return {"rows": rows, "parameters": ticked,
            "stress_reproduced": float(np.max(np.abs(replayed - recorded)) / scale),
            "du_dp_norm": float(np.linalg.norm(du_dp)),
            "locations": len({row["location"] for row in rows}),
            "increments": summary["scope"]["output_increments"]}


def _row(output: str, increment: int, location: str, component: int, value: float,
         derivatives: Any, names: list[str]) -> dict[str, Any]:
    row = {"output": output, "increment": increment, "location": location,
           "component": component, "value": float(value)}
    row.update({f"d/d{name}": float(derivative) for name, derivative in zip(names, derivatives)})
    return row


def render_outputs_and_parameters(st, material: str | None, mapping: str | None,
                                  model: str | None, *, mapping_upload: Any = None,
                                  request_supplied: bool = False) -> dict[str, Any] | None:
    """Section 3 of the slide: tick parameters, choose output and region.

    Returns the request those choices stand for, or None (with the reason on
    the screen) when they do not make one yet.
    """
    st.subheader("3. Outputs and parameters")
    found = mapping_path(material, mapping)
    try:
        if found is not None:
            payload, origin = json.loads(found.read_text(encoding="utf-8")), str(found)
        elif mapping_upload is not None:
            payload, origin = json.loads(mapping_upload.getvalue()), mapping_upload.name
        else:
            st.caption("The parameters are read from the Mapping.json the provider build "
                       "wrote beside the compiled object (or the one chosen under Advanced). "
                       "Point at the object to tick them.")
            return None
        names = mapping_parameters(payload)
        st.session_state["request_mapping_indices"] = {
            str(entry["name"]): int(entry["props_index"]) for entry in payload["parameters"]}
    except (OSError, ValueError, KeyError, TypeError) as error:
        st.error(f"Mapping.json could not be read: {error}")
        return None
    st.caption(f"parameters from `{origin}`")
    ticks = st.columns(max(len(names), 1))
    parameters = [name for column, name in zip(ticks, names)
                  if column.checkbox(name, value=True, key=f"request_param_{name}")]
    sets: dict[str, dict[str, list[int]]] = {"nodes": {}, "elements": {}}
    if model and Path(model).is_file():
        try:
            sets = model_sets(Path(model))
        except (OSError, ValueError) as error:
            st.warning(f"could not read the sets in the .inp: {error}")
    left, middle, right = st.columns(3)
    output = left.selectbox("output", list(OUTPUTS), key="request_output_choice")
    region = middle.selectbox("region", region_choices(sets, OUTPUTS[output][0]),
                              key=f"request_region_{OUTPUTS[output][0] in ('U', 'RF')}")
    reduction = right.selectbox("summary over the region", REDUCTIONS, key="request_reduction")
    increments = st.radio("increments", list(INCREMENTS), horizontal=True, key="request_increments")
    if request_supplied:
        st.caption("A sensitivity_request.json was supplied above; Solve uses it instead of "
                   "these choices.")
    try:
        return build_request(parameters=parameters, output=output, region=region,
                             reduction=reduction, increments=increments, sets=sets)
    except ValueError as error:
        st.warning(str(error))
        return None


def render_full_field(st, completed: str | None) -> None:
    """The full-field sensitivities of the last Solve, in the chosen region."""
    if not completed or not (Path(completed) / "private" / "result.json").is_file():
        return
    try:
        view = full_field(Path(completed), st.session_state.get("request_mapping_indices"))
    except (OSError, ValueError, KeyError, IndexError) as error:
        st.warning(f"full-field view unavailable: {error}")
        return
    rows = view["rows"]
    increments = ", ".join(map(str, view["increments"]))
    if rows:
        st.success(f"Solved — d({rows[0]['output']})/dp at all {view['locations']:,} locations "
                   f"of the region (increment {increments})")
    else:
        st.success(f"Solved (increment {increments}). Reaction derivatives are not stored per "
                   "node by the request; the reduced value is in sensitivity_results.json above.")
    left, right = st.columns(2)
    left.metric("stress reproduced vs. the .odb", f"{view['stress_reproduced']:.1e}")
    right.metric("solution sensitivity ‖du/dp‖", f"{view['du_dp_norm']:.3g}")
    st.caption("Stress: max |replayed − .odb| / max(|.odb|, 1) at the last output increment. "
               "‖du/dp‖: Frobenius norm over every DOF and the ticked parameters. "
               + ("Governing parameter: the largest |p·∂y/∂p| at that location, p read "
                  "from the .inp's USER MATERIAL constants at the Mapping.json PROPS index."
                  if rows and "governing parameter" in rows[0] else ""))
    # The table's own toolbar downloads it as CSV; the public downloads above
    # stay the three the request writes.
    if not rows:
        return
    numbers = [key for key in rows[0] if key == "value" or key.startswith("d/d")]
    st.dataframe(rows, hide_index=True, width="stretch",
                 column_config={key: st.column_config.NumberColumn(key, format="%.6e") for key in numbers})
