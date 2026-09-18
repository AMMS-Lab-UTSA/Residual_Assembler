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

from history_support import CONTRACTS, EXAMPLE, HERE, RA_ROOT, beam_model, umat_repository  # noqa: F401


def pytest_configure(config):
    config.addinivalue_line("markers", "abaqus: needs a licensed Abaqus installation")


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


@pytest.fixture
def times():
    return lambda steps: np.linspace(0.0, 1.0, steps + 1)
