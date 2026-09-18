"""Strain-history chain rule for the fingerprint-pinned m3_j2 provider only."""

import numpy as np


J2_SOURCE_HASH = "9b779f0c6cadf9c4"


def require_j2(material):
    contract = material.contract
    if (contract.get("regular_source_hash") != J2_SOURCE_HASH
            or contract.get("kinematics") != "small_strain"
            or (material.ntens, material.nprops, material.nstatev, material.nparam) != (6, 4, 1, 4)
            or material.props_index != [1, 2, 3, 4]
            or material.params != ["E", "nu", "SIGY0", "H"]):
        raise ValueError("total-history replay requires the fingerprint-pinned m3_j2 small-strain provider")
    signature = contract.get("symbols", {}).get("oti_eval_signature", [])
    if len(signature) != 18 or len(contract.get("march", {}).get("signature", [])) != 11:
        raise ValueError("m3_j2 requires the 18-argument carry EVAL and 11-argument MARCH ABI")
    if not material.has_march:
        raise ValueError("declared MARCH ABI is not available in the linked provider")


def real_array(values, shape, name):
    array = np.asarray(values)
    if array.dtype.kind not in "fiu" or array.shape != shape or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} requires finite extracted real coefficients of shape {shape}; live OTI is not accepted")
    return array.astype(float, copy=False)


def step_j2(material, props, stress, state, dstress, dstate, increment, dincrement, dt):
    """Differentiate the pinned update with respect to parameters AND strain.

    The original update depends on strain only via stress + C_elastic*increment.
    Thus seeding dstress + C_elastic*dincrement before EVAL is exactly the
    multivariate chain rule, including the physical STATEV derivative. C_elastic
    is obtained from the ORIGINAL UMAT at virgin zero strain. EVAL itself adds
    the direct parameter derivatives, including dC/dp*increment. No constitutive
    return mapping or finite differences are implemented here.
    """
    require_j2(material)
    props = real_array(props, (4,), "props")
    if props[0] <= 0 or not -1 < props[1] < 0.5 or props[2] <= 0 or props[3] < 0:
        raise ValueError("m3_j2 requires E>0, -1<nu<0.5, SIGY0>0, H>=0")
    stress = real_array(stress, (6,), "stress")
    state = real_array(state, (1,), "state")
    dstress = real_array(dstress, (6, 4), "dstress")
    dstate = real_array(dstate, (1, 4), "dstate")
    increment = real_array(increment, (6,), "increment")
    dincrement = real_array(dincrement, (6, 4), "dincrement")
    if not np.isfinite(dt) or dt <= 0:
        raise ValueError("dt must be finite and positive")
    elastic = material._reg_step(props, np.zeros(6), np.zeros(1), np.zeros(6), dt)[2]
    return material._oti_step(props, stress, state, dstress + elastic @ dincrement,
                              dstate, increment, dt)