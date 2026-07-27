"""Single source of truth for every resasm interface version.

Each JSON-contract schema tag and the shared C-ABI / material-package contract
version is declared here exactly once. Modules import their schema constant from
this module rather than re-declaring the literal, and `check()` gives every
loader the same actionable drift message. `contract_version()` returns the
combined hash of the C-ABI + material-package contract (see
replay/contract/CONTRACT_VERSION.json) so a stale binary is caught early.
"""
from __future__ import annotations

import json
import os
from typing import Dict

# --- Program 2: job / result / request contracts ---------------------------
SENSITIVITY_JOB = "resasm_sensitivity_job_v1"
SENSITIVITY_RESULT = "resasm_sensitivity_result_v1"
SENSITIVITY_REQUEST = "resasm_sensitivity_request_v1"
VERIFICATION_JOB = "resasm_verification_job_v1"
CONFIG = "resasm_config_v1"

# --- replay / record / material --------------------------------------------
REPLAY_RECORD = "resasm_replay_record_v1"
MATERIAL_PACKAGE = "resasm_material_package_v1"
UMAT_OTI_CONTRACT = "resasm_umat_oti_contract_v1"
DERIVATIVE_FIELDS = "resasm_derivative_fields_v1"
SDV_LAYOUT = "resasm_sdv_layout_v1"

# --- shared C-ABI (matches replay/contract/resasm_mat_abi_v1.h) -------------
MAT_ABI = "resasm_mat_abi_v1"

# --- transform contract authored on the Program 1 side (kept in sync) ------
UMAT_TRANSFORM = "resasm_umat_transform_v2"

# name -> tag, for enumeration / drift guards
SCHEMAS: Dict[str, str] = {
    "SENSITIVITY_JOB": SENSITIVITY_JOB,
    "SENSITIVITY_RESULT": SENSITIVITY_RESULT,
    "SENSITIVITY_REQUEST": SENSITIVITY_REQUEST,
    "VERIFICATION_JOB": VERIFICATION_JOB,
    "CONFIG": CONFIG,
    "REPLAY_RECORD": REPLAY_RECORD,
    "MATERIAL_PACKAGE": MATERIAL_PACKAGE,
    "UMAT_OTI_CONTRACT": UMAT_OTI_CONTRACT,
    "DERIVATIVE_FIELDS": DERIVATIVE_FIELDS,
    "SDV_LAYOUT": SDV_LAYOUT,
    "MAT_ABI": MAT_ABI,
    "UMAT_TRANSFORM": UMAT_TRANSFORM,
}

_CONTRACT_JSON = os.path.join(os.path.dirname(__file__), os.pardir, "replay",
                              "contract", "CONTRACT_VERSION.json")


def check(expected: str, got, what: str = "schema") -> None:
    """Raise a uniform, actionable error when a loaded tag does not match."""
    if got != expected:
        raise ValueError("%s must be %r, got %r (interface drift -- regenerate or "
                         "pin to this resasm version)" % (what, expected, got))


def contract_version() -> str:
    """Combined hash of the C-ABI + material-package contract (empty if absent)."""
    try:
        return json.load(open(_CONTRACT_JSON)).get("combined_hash", "")
    except (OSError, ValueError):
        return ""


if __name__ == "__main__":
    print("resasm interface versions")
    for name, tag in SCHEMAS.items():
        print("  %-22s %s" % (name, tag))
    print("  %-22s %s" % ("CONTRACT_HASH", contract_version()))
