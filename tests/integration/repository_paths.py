"""Explicit source inputs for cross-repository integration tests."""

import os
from pathlib import Path


def umat_repo_root() -> Path:
    root = Path(os.environ.get(
        "UMAT_OTI_REPO",
        Path(__file__).resolve().parents[3] / "UMAT_source_transformation",
    )).expanduser().resolve()
    if not (root / "src/umat_oti").is_dir():
        raise FileNotFoundError(f"UMAT source checkout not found: {root}; set UMAT_OTI_REPO")
    return root