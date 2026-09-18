"""Shared fixtures for the history-replay tests.

Providers are built from the UMAT-OTI repository with the real compiler
(``umat_oti.provider.build_provider``); nothing constitutive is mocked. The
UMAT-OTI checkout is ``$UMAT_OTI_REPO`` or the repository that the imported
``umat_oti`` package lives in.
"""
import json
import os
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
RA_ROOT = HERE.parents[1]
EXAMPLE = RA_ROOT / "examples" / "replay_history"


def pytest_configure(config):
    config.addinivalue_line("markers", "abaqus: needs a licensed Abaqus installation")


def umat_repository() -> Path:
    candidates = []
    if os.environ.get("UMAT_OTI_REPO"):
        candidates.append(Path(os.environ["UMAT_OTI_REPO"]))
    import umat_oti
    candidates.append(Path(umat_oti.__file__).resolve().parents[2])
    for candidate in candidates:
        if (candidate / "parameter_sensitivity" / "models" / "m3_j2" / "contract_v2.json").is_file():
            return candidate
    raise RuntimeError("UMAT-OTI repository with parameter_sensitivity/models not found "
                       "(set UMAT_OTI_REPO)")


CONTRACTS = {
    "m3_j2": ("parameter_sensitivity/models/m3_j2/contract_v2.json", "umat_m3_j2_oti"),
    "m6_fcc": ("parameter_sensitivity/models/m6_fcc/contract_v2.json", "umat_m6_fcc_oti"),
    "damage": ("tests/fixtures/provider_total_strain/contract_v2.json", "umat_total_strain_damage_oti"),
}


@pytest.fixture(scope="session")
def provider_factory(tmp_path_factory):
    from umat_oti.provider import build_provider
    from residual_core.replay.history_material import HistoryMaterial
    built = {}

    def get(name):
        if name not in built:
            relative, _ = CONTRACTS[name]
            directory = tmp_path_factory.mktemp("provider-" + name)
            result = build_provider(umat_repository() / relative, directory / "out")
            contract = json.loads(Path(result["contract"]).read_text())
            built[name] = (Path(result["object"]), contract,
                           HistoryMaterial(result["object"], contract, str(directory / "link")))
        return built[name]
    return get


def beam_model(tmp_path, *, n=(6, 2, 2), push=0.03, steps=6, props=(200000.0, 0.3, 250.0, 2000.0),
               depvar=1, cload=None, name="beam.inp", extra=None):
    """Write a structured beam deck with the example generator and read it."""
    import sys
    sys.path.insert(0, str(EXAMPLE))
    from make_beam_deck import deck
    from residual_core.replay.history_inputs import read_history_model
    text = deck(n, push, steps, list(props), depvar, cload)
    if extra:
        text = extra(text)
    path = tmp_path / name
    path.write_text(text)
    return path, read_history_model(path)


@pytest.fixture
def times():
    return lambda steps: np.linspace(0.0, 1.0, steps + 1)
