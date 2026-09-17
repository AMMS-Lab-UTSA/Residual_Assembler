"""Launcher for the Residual_Assembler Streamlit GUI.

    streamlit run scripts/app.py

Adds the repository root to ``sys.path`` so the app runs straight from a
checkout, with or without ``pip install -e .``.
"""

from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from residual_core.app.streamlit_app import main  # noqa: E402


if __name__ == "__main__":
    main()
