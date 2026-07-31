"""Tiny YAML-subset loader (zero dependency).

The user config (`resasm.yml`) uses only a small YAML subset: nested mappings by
indentation, scalar values, inline flow lists ``[a, b]`` and ``- `` block lists.
If PyYAML is installed it is used instead (full YAML). Otherwise this minimal
reader handles the documented config shape.
"""

from __future__ import annotations

import json
from typing import Any, List


def load_str(text: str) -> Any:
    try:
        import yaml  # type: ignore
        return yaml.safe_load(text)
    except Exception:
        pass
    return _mini_load(text)


def load_file(path: str) -> Any:
    with open(path, "r", encoding="utf-8") as fh:
        text = fh.read()
    # allow a .json config too
    stripped = text.lstrip()
    if stripped[:1] in "{[":
        return json.loads(text)
    return load_str(text)


# --------------------------------------------------------------------------- #
def _indent_of(line: str) -> int:
    return len(line) - len(line.lstrip(" "))


def _strip_inline_comment(s: str) -> str:
    if not s:
        return s
    q = s[0] in "\"'"
    if q:
        return s
    # cut a comment that follows whitespace
    for i in range(1, len(s)):
        if s[i] == "#" and s[i - 1] == " ":
            return s[:i].strip()
    return s.strip()


def _scalar(s: str) -> Any:
    s = _strip_inline_comment(s.strip())
    if s == "" or s in ("null", "~", "None"):
        return None
    if len(s) >= 2 and s[0] in "\"'" and s[-1] == s[0]:
        return s[1:-1]
    low = s.lower()
    if low in ("true", "yes", "on"):
        return True
    if low in ("false", "no", "off"):
        return False
    if s.startswith("[") and s.endswith("]"):
        inner = s[1:-1].strip()
        if not inner:
            return []
        return [_scalar(x) for x in _split_flow(inner)]
    try:
        return int(s)
    except ValueError:
        pass
    try:
        return float(s)
    except ValueError:
        pass
    return s


def _split_flow(inner: str) -> List[str]:
    out, depth, cur = [], 0, ""
    for ch in inner:
        if ch == "[":
            depth += 1
        elif ch == "]":
            depth -= 1
        if ch == "," and depth == 0:
            out.append(cur)
            cur = ""
        else:
            cur += ch
    if cur.strip():
        out.append(cur)
    return out


def _mini_load(text: str) -> Any:
    lines = []
    for raw in text.splitlines():
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        lines.append(raw.rstrip())
    pos = [0]

    def parse_block(indent: int) -> Any:
        result: Any = None
        while pos[0] < len(lines):
            line = lines[pos[0]]
            cur = _indent_of(line)
            if cur < indent:
                break
            if cur > indent:
                break
            stripped = line.strip()
            pos[0] += 1
            if stripped.startswith("- "):
                if result is None:
                    result = []
                item = stripped[2:].strip()
                if item == "":
                    nxt = _indent_of(lines[pos[0]]) if pos[0] < len(lines) else indent
                    result.append(parse_block(nxt) if nxt > indent else None)
                else:
                    result.append(_scalar(item))
                continue
            key, sep, val = stripped.partition(":")
            key = key.strip()
            val = val.strip()
            if result is None:
                result = {}
            if val == "":
                nxt = _indent_of(lines[pos[0]]) if pos[0] < len(lines) else indent
                if nxt > indent:
                    result[key] = parse_block(nxt)
                else:
                    result[key] = None
            else:
                result[key] = _scalar(val)
        return result

    parsed = parse_block(0)
    return parsed if parsed is not None else {}
