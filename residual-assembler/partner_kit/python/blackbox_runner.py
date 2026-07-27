"""Black-box executable runner — the no-linking integration path.

Some partners will not link against our code. Instead they expose a command:

    ./their_solver --input request.json --output response.json

The kit writes ``request.json``, runs the command, and reads ``response.json``
(and an optional ``response.npz`` for large arrays). Everything proprietary stays
behind the executable boundary — the kit performs no overloading and never sees
the mesh or source.

Request schema  (resasm-partner-request/1):
    { order, solution[], parameters{name:val}, seed_directions{name:basis},
      time[2], dtime }

Response schema (resasm-partner-response/1):
    { status: "ok"|"error",
      residual_coefficients: [[...]],        # ndof x m  (R^(p)); or a ref to .npz
      residual_real: [...],                  # optional R (order 0)
      tangent: [[...]] | null,               # optional dense T
      arrays_npz: "response.npz" | null,      # optional array sidecar
      diagnostics: {...}, message: "" }
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from typing import Any, Dict, Optional

import numpy as np


class BlackBoxRunner:
    def __init__(self, command, workdir: Optional[str] = None, timeout: Optional[float] = None):
        """``command`` is a list or string. ``{input}`` / ``{output}`` placeholders
        are substituted; if absent, ``--input``/``--output`` are appended."""
        self.command = command
        self.workdir = workdir or tempfile.mkdtemp(prefix="resasm_partner_")
        self.timeout = timeout

    def _argv(self, in_path, out_path):
        cmd = self.command
        if isinstance(cmd, str):
            if "{input}" in cmd or "{output}" in cmd:
                cmd = cmd.format(input=in_path, output=out_path)
                return cmd, True     # shell string
            return "%s --input %s --output %s" % (cmd, in_path, out_path), True
        # list form
        argv = list(cmd)
        if any("{input}" in a or "{output}" in a for a in argv):
            argv = [a.format(input=in_path, output=out_path) for a in argv]
        else:
            argv += ["--input", in_path, "--output", out_path]
        return argv, False

    def run(self, request: Dict[str, Any]) -> Dict[str, Any]:
        os.makedirs(self.workdir, exist_ok=True)
        in_path = os.path.join(self.workdir, "request.json")
        out_path = os.path.join(self.workdir, "response.json")
        with open(in_path, "w", encoding="utf-8") as fh:
            json.dump(request, fh, indent=2)
        argv, shell = self._argv(in_path, out_path)
        proc = subprocess.run(argv, shell=shell, cwd=self.workdir,
                              capture_output=True, text=True, timeout=self.timeout)
        if proc.returncode != 0:
            raise RuntimeError("black-box solver failed (rc=%d): %s"
                               % (proc.returncode, proc.stderr[-500:]))
        if not os.path.exists(out_path):
            raise RuntimeError("black-box solver produced no response.json")
        with open(out_path, "r", encoding="utf-8") as fh:
            resp = json.load(fh)
        # optional array sidecar
        npz_ref = resp.get("arrays_npz")
        if npz_ref:
            npz_path = npz_ref if os.path.isabs(npz_ref) else os.path.join(self.workdir, npz_ref)
            if os.path.exists(npz_path):
                data = np.load(npz_path, allow_pickle=True)
                resp["_arrays"] = {k: data[k] for k in data.files}
        if resp.get("status", "ok") == "error":
            raise RuntimeError("black-box solver reported error: %s"
                               % resp.get("message", ""))
        return resp

    def residual_coefficients(self, response: Dict[str, Any]) -> np.ndarray:
        """Extract R^(p) (ndof x m) from a response (inline or from the npz)."""
        if "residual_coefficients" in response and response["residual_coefficients"] is not None:
            return np.asarray(response["residual_coefficients"], float)
        arrays = response.get("_arrays", {})
        for key in ("residual_coefficients", "R", "R_order_1"):
            if key in arrays:
                return np.asarray(arrays[key], float)
        raise KeyError("response has no residual_coefficients")
