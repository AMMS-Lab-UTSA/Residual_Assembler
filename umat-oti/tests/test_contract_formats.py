"""The compact transformation contract and the separate verification case.

Pins the user-facing contract: the developer writes only source, interface
sizes, parameters (name -> PROPS index), and requested derivatives; nparam and
everything downstream is inferred; property values live in a separate
verification case. All accepted input shapes normalize to one internal rep.
"""
import json
import os
import sys

import pytest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "oti_provider"))
sys.path.insert(0, os.path.join(_ROOT, "src"))

from umat_transform import (  # noqa: E402
    ContractError,
    find_verification_case,
    normalize_contract,
    resolve_props,
)

# the documented compact contract (J2), and the equivalent legacy compact-v2
NEW = {
    "material": "j2_linear_hardening",
    "source": "umat_j2.for",
    "kinematics": "small_strain",
    "dimensions": {"stress_components": 6, "material_properties": 4, "state_variables": 1},
    "parameters": {"E": 1, "nu": 2, "SIGY0": 3, "H": 4},
    "derivatives": {"stress": True, "statev": "auto"},
}
V2 = {
    "schema": "resasm_umat_transform_v2",
    "source": {"main_file": "umat_j2.for"},
    "kinematics": "small_strain",
    "dimensions": {"ntens": 6, "nprops": 4, "nstatev": 1},
    "parameters": [{"name": "E", "props_index": 1}, {"name": "nu", "props_index": 2},
                   {"name": "SIGY0", "props_index": 3}, {"name": "H", "props_index": 4}],
}


def test_new_and_v2_normalize_identically():
    a, fa = normalize_contract(NEW, "new.json")
    b, fb = normalize_contract(V2, "v2.json")
    assert (fa, fb) == ("compact", "compact-v2")
    assert a["parameters"] == b["parameters"]
    assert a["dimensions"] == b["dimensions"] == {"ntens": 6, "nprops": 4, "nstatev": 1}
    assert a["source"] == b["source"] == {"main_file": "umat_j2.for", "additional_files": []}


def test_nparam_is_inferred_not_declared():
    c, _ = normalize_contract(NEW, "p.json")
    assert len(c["parameters"]) == 4
    assert [p["oti_direction"] for p in c["parameters"]] == [1, 2, 3, 4]


def test_manual_nparam_is_rejected():
    with pytest.raises(ContractError) as e:
        normalize_contract({**NEW, "nparam": 9}, "bad.json")
    assert "nparam is inferred" in str(e.value)


def test_dimensions_forbid_nparam():
    bad = {**NEW, "dimensions": {**NEW["dimensions"], "nparam": 4}}
    with pytest.raises(ContractError):
        normalize_contract(bad, "bad.json")


def test_clear_and_raw_dimension_names_both_work():
    raw_names = {**NEW, "dimensions": {"ntens": 6, "nprops": 4, "nstatev": 1}}
    c, _ = normalize_contract(raw_names, "raw.json")
    assert c["dimensions"] == {"ntens": 6, "nprops": 4, "nstatev": 1}


def test_string_and_object_source():
    s, _ = normalize_contract({**NEW, "source": "u.for"}, "s.json")
    assert s["source"] == {"main_file": "u.for", "additional_files": []}
    o, _ = normalize_contract(
        {**NEW, "source": {"main": "u.for", "additional_files": ["dep.for"]}}, "o.json")
    assert o["source"] == {"main_file": "u.for", "additional_files": ["dep.for"]}


def test_statev_auto_true_false():
    auto, _ = normalize_contract(NEW, "a.json")
    assert auto["_derivatives"]["path_dependent"] is True         # nstatev=1
    elastic = {**NEW, "dimensions": {"stress_components": 6, "material_properties": 2,
                                     "state_variables": 0},
               "parameters": {"E": 1, "nu": 2}, "derivatives": {"stress": True, "statev": False}}
    e, _ = normalize_contract(elastic, "e.json")
    assert e["_derivatives"]["path_dependent"] is False


def test_statev_false_with_state_is_rejected():
    with pytest.raises(ContractError) as e:
        normalize_contract({**NEW, "derivatives": {"stress": True, "statev": False}}, "b.json")
    assert "path-dependent" in str(e.value)


def test_stress_must_be_requested():
    with pytest.raises(ContractError):
        normalize_contract({**NEW, "derivatives": {"stress": False}}, "b.json")


def test_out_of_range_and_duplicate_parameters_rejected():
    with pytest.raises(ContractError):
        normalize_contract({**NEW, "parameters": {"E": 1, "X": 9}}, "b.json")
    with pytest.raises(ContractError):
        normalize_contract({**NEW, "parameters": {"E": 1, "F": 1}}, "b.json")  # dup index


def test_verbose_still_accepted():
    verbose = {
        "schema": "resasm_umat_transform_v2",
        "source": {"entry_point": "UMAT", "main_file": "umat.for"},
        "interface": {"ntens": 6, "nprops": 2, "nstatev": 0, "kinematics": "small_strain"},
        "derivative_requests": [{"seed": {"components": [
            {"index": 1, "name": "E", "oti_direction": 1},
            {"index": 2, "name": "nu", "oti_direction": 2}]}}],
        "resasm_provider": {"props_values": [210000.0, 0.3]},
    }
    c, fmt = normalize_contract(verbose, "v.json")
    assert fmt == "verbose"
    assert [p["name"] for p in c["parameters"]] == ["E", "nu"]
    assert c["_verification"]["props_values"] == [210000.0, 0.3]


def test_completed_contract_is_not_an_input():
    with pytest.raises(ContractError) as e:
        normalize_contract({"schema": "resasm_umat_oti_contract_v1"}, "done.json")
    assert "COMPLETED" in str(e.value)


# --- verification case: property values live outside the contract ---

def test_props_come_from_verification_case(tmp_path):
    (tmp_path / "verification.json").write_text(json.dumps({"props": [210000.0, 0.3, 250.0, 2000.0],
                                                            "loading_path": "uniaxial"}))
    c, _ = normalize_contract(NEW, "c.json")
    values, origin, loading = resolve_props(c, 4, None, cdir=str(tmp_path))
    assert values == [210000.0, 0.3, 250.0, 2000.0]
    assert "verification case" in origin
    assert loading == "uniaxial"


def test_props_cli_overrides_verification(tmp_path):
    (tmp_path / "verification.json").write_text(json.dumps({"props": [1, 2, 3, 4]}))
    c, _ = normalize_contract(NEW, "c.json")
    values, origin, _ = resolve_props(c, 4, [9, 8, 7, 6], cdir=str(tmp_path))
    assert values == [9, 8, 7, 6] and origin == "--props"


def test_missing_props_is_an_explicit_error(tmp_path):
    c, _ = normalize_contract(NEW, "c.json")
    with pytest.raises(ContractError) as e:
        resolve_props(c, 4, None, cdir=str(tmp_path))
    assert "verification case" in str(e.value) and "operating point" in str(e.value)


def test_wrong_length_props_rejected(tmp_path):
    (tmp_path / "verification.json").write_text(json.dumps({"props": [1, 2, 3]}))
    c, _ = normalize_contract(NEW, "c.json")
    with pytest.raises(ContractError) as e:
        resolve_props(c, 4, None, cdir=str(tmp_path))
    assert "material_properties=4" in str(e.value)


def test_find_verification_case(tmp_path):
    assert find_verification_case(str(tmp_path)) is None
    (tmp_path / "verify.json").write_text("{}")
    assert os.path.basename(find_verification_case(str(tmp_path))) == "verify.json"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
