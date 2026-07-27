"""User-facing layer: high-level API, CLI, config and runnable examples.

Most users only need ``residual_core.ResidualProblem`` (re-exported at the
package root) or the ``resasm`` command line. Everything here orchestrates the
agnostic core + backend registries; no physics lives in this layer.
"""

from .wizard import ResidualProblem
from .config import Config, write_config_template

__all__ = ["ResidualProblem", "Config", "write_config_template"]
