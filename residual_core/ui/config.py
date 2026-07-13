"""Configuration for advanced/developer users (progressive disclosure).

Basic users never touch this — auto-detection covers them. A config is only
needed when auto-detection cannot infer something (custom formulation policy,
field export paths, material parameter overrides). ``resasm doctor --write-config-
template`` emits one pre-filled from what was detected.

YAML is used if PyYAML is installed; otherwise JSON is used (same schema). The
emitted template is human-readable YAML-style text regardless.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from typing import Any, Dict, Optional

try:
    import yaml   # type: ignore
    _HAVE_YAML = True
except Exception:
    _HAVE_YAML = False


@dataclass
class Config:
    """User-facing knobs. Everything optional; absent -> auto-detect."""
    mode: Optional[str] = None                 # stress-driven | material-replay | ...
    odb: Optional[str] = None                  # exported field source (ODB/CSV/JSON)
    subroutine: Optional[str] = None           # UMAT/UEL source file
    formulation_policy: Dict[str, str] = field(default_factory=dict)  # etype -> backend
    material_backend: Dict[str, str] = field(default_factory=dict)    # matname -> backend
    material_parameters: Dict[str, list] = field(default_factory=dict)  # matname -> PROPS
    options: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "mode": self.mode,
            "odb": self.odb,
            "subroutine": self.subroutine,
            "formulation_policy": dict(self.formulation_policy),
            "material_backend": dict(self.material_backend),
            "material_parameters": {k: list(v) for k, v in self.material_parameters.items()},
            "options": dict(self.options),
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "Config":
        d = d or {}
        return cls(
            mode=d.get("mode"),
            odb=d.get("odb"),
            subroutine=d.get("subroutine"),
            formulation_policy=dict(d.get("formulation_policy", {})),
            material_backend=dict(d.get("material_backend", {})),
            material_parameters={k: list(v) for k, v in d.get("material_parameters", {}).items()},
            options=dict(d.get("options", {})),
        )

    # ------------------------------------------------------------------ #
    def save(self, path: str) -> None:
        data = self.to_dict()
        if path.lower().endswith((".yml", ".yaml")) and _HAVE_YAML:
            with open(path, "w", encoding="utf-8") as fh:
                yaml.safe_dump(data, fh, sort_keys=False)
        elif path.lower().endswith((".yml", ".yaml")):
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(_yaml_text(data))
        else:
            with open(path, "w", encoding="utf-8") as fh:
                json.dump(data, fh, indent=2)

    @classmethod
    def load(cls, path: str) -> "Config":
        with open(path, "r", encoding="utf-8") as fh:
            text = fh.read()
        if path.lower().endswith((".yml", ".yaml")) and _HAVE_YAML:
            return cls.from_dict(yaml.safe_load(text) or {})
        try:
            return cls.from_dict(json.loads(text))
        except json.JSONDecodeError:
            if _HAVE_YAML:
                return cls.from_dict(yaml.safe_load(text) or {})
            raise ValueError(
                "cannot parse %s: install PyYAML or use a .json config" % path)


def _yaml_text(data: Dict[str, Any]) -> str:
    """Minimal YAML emitter (no dependency) for the config template."""
    lines = []

    def emit(obj, indent=0):
        pad = "  " * indent
        if isinstance(obj, dict):
            if not obj:
                return " {}"
            out = ""
            for k, v in obj.items():
                rendered = emit(v, indent + 1)
                if rendered.startswith("\n") or not rendered.startswith(" "):
                    out += "\n%s%s:%s" % (pad, k, rendered)
                else:
                    out += "\n%s%s:%s" % (pad, k, rendered)
            return out
        if isinstance(obj, list):
            if not obj:
                return " []"
            out = ""
            for v in obj:
                out += "\n%s- %s" % (pad, json.dumps(v))
            return out
        if obj is None:
            return " null"
        return " %s" % json.dumps(obj)

    for k, v in data.items():
        lines.append("%s:%s" % (k, emit(v, 1)))
    return "\n".join(lines) + "\n"


def write_config_template(path: str, detected: Optional[Config] = None) -> str:
    """Write a commented starter config, pre-filled from detection if given."""
    cfg = detected or Config()
    header = (
        "# resasm configuration (auto-generated template)\n"
        "# Only fill what auto-detection could NOT infer. Delete the rest.\n"
        "#   mode                : stress-driven | material-replay | direct-residual | formulation\n"
        "#   odb                 : path to exported field (ODB/CSV/JSON) for stress-driven\n"
        "#   subroutine          : UMAT/UEL source for material-replay / direct-residual\n"
        "#   formulation_policy  : {ELEMENT_TYPE: backend_name}\n"
        "#   material_backend    : {material_name: backend_name}\n"
        "#   material_parameters : {material_name: [PROPS...]}\n")
    if path.lower().endswith((".yml", ".yaml")):
        text = header + _yaml_text(cfg.to_dict())
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(text)
    else:
        cfg.save(path)
        text = header
    return text
