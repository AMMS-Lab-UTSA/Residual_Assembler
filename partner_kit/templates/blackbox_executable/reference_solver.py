#!/usr/bin/env python
"""reference_solver.py — dependency-free black-box solver template.

Speaks the kit's contract:  reference_solver.py --input request.json --output response.json

Computes R^(1) by CENTRAL FINITE DIFFERENCES of a private residual, so it works
for any model you drop into ``residual(u, params)`` without needing AD. (Finite
differences are for getting started / validation; prefer a dual/OTI or
source-transformed derivative for production accuracy.)

Replace ``residual`` with a call into your own code. Nothing about your model
leaves the process except the coefficients you write.
"""

import argparse
import json


def residual(u, params):
    """PRIVATE model residual. Return a length-ndof list. Edit this.
    Example: 1-DOF cubic spring R = k u^3 - f."""
    k = params.get("k", 1.0)
    f = params.get("f", 0.0)
    return [k * u[0] ** 3 - f]


def tangent(u, params, h=1e-7):
    """Optional dR/du by finite differences (dense). Return ndof x ndof or None."""
    n = len(u)
    R0 = residual(u, params)
    ndof = len(R0)
    T = [[0.0] * n for _ in range(ndof)]
    for j in range(n):
        up = list(u); up[j] += h
        um = list(u); um[j] -= h
        Rp = residual(up, params); Rm = residual(um, params)
        for i in range(ndof):
            T[i][j] = (Rp[i] - Rm[i]) / (2 * h)
    return T


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    with open(args.input, "r", encoding="utf-8") as fh:
        req = json.load(fh)

    u = list(req["solution"])
    params = dict(req["parameters"])
    seeds = req["seed_directions"]                 # {param: basis_index}
    cols = sorted(seeds.items(), key=lambda kv: kv[1])
    R0 = residual(u, params)
    ndof = len(R0)
    m = len(cols)

    # R^(1)[:, j] = dR/dp_j by central finite differences
    R1 = [[0.0] * m for _ in range(ndof)]
    for j, (pname, _basis) in enumerate(cols):
        p0 = float(params.get(pname, 0.0))
        dp = 1e-6 * max(1.0, abs(p0))
        pp = dict(params); pp[pname] = p0 + dp
        pm = dict(params); pm[pname] = p0 - dp
        Rp = residual(u, pp); Rm = residual(u, pm)
        for i in range(ndof):
            R1[i][j] = (Rp[i] - Rm[i]) / (2 * dp)

    resp = {
        "schema": "resasm-partner-response/1",
        "status": "ok",
        "order": int(req.get("order", 1)),
        "residual_real": R0,
        "residual_coefficients": R1,
        "tangent": tangent(u, params),
        "diagnostics": {"method": "central-finite-difference", "ndof": ndof, "m": m},
        "message": "",
    }
    with open(args.output, "w", encoding="utf-8") as fh:
        json.dump(resp, fh, indent=2)


if __name__ == "__main__":
    main()
