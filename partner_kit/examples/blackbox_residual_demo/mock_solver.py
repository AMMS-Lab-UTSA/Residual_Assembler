#!/usr/bin/env python
"""Mock black-box solver — stands in for a partner's proprietary executable.

Implements the kit's black-box contract:

    python mock_solver.py --input request.json --output response.json

It reads the requested solution, parameter values and seeded directions, computes
the residual and its parameter derivatives *with its own internal method* (here,
closed form for the cubic spring R = k u^3 - f), and returns the residual
coefficients R^(1). The kit never sees how this is done — the model stays behind
the executable boundary.

Request  (resasm-partner-request/1):  order, solution[], parameters{}, seed_directions{}, time[], dtime
Response (resasm-partner-response/1): status, residual_coefficients[[...]], residual_real[], tangent[[...]], diagnostics{}
"""

import argparse
import json


def solve(request):
    u = request["solution"][0]
    k = request["parameters"]["k"]
    f = request["parameters"]["f"] if "f" in request["parameters"] else 16.0
    seeds = request["seed_directions"]              # {param: basis_index}
    m = len(seeds)

    # residual (order 0)
    R0 = [k * u ** 3 - f]

    # dR/dparam for each seeded parameter (the partner's own AD/analytic method)
    deriv = {"k": u ** 3, "f": -1.0}
    # order columns by basis index
    cols = sorted(seeds.items(), key=lambda kv: kv[1])
    R1 = [[0.0] * m]                                 # ndof(=1) x m
    for j, (pname, _basis) in enumerate(cols):
        R1[0][j] = float(deriv.get(pname, 0.0))

    tangent = [[3.0 * k * u * u]]                    # optional dR/du
    return {
        "schema": "resasm-partner-response/1",
        "status": "ok",
        "residual_real": R0,
        "residual_coefficients": R1,
        "tangent": tangent,
        "diagnostics": {"method": "internal-closed-form", "ndof": 1, "m": m},
        "message": "",
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    with open(args.input, "r", encoding="utf-8") as fh:
        request = json.load(fh)
    response = solve(request)
    with open(args.output, "w", encoding="utf-8") as fh:
        json.dump(response, fh, indent=2)


if __name__ == "__main__":
    main()
