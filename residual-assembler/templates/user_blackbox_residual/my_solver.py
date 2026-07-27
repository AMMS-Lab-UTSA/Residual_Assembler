#!/usr/bin/env python
"""Black-box residual solver (template).

Your private model/code stays entirely inside this executable. The framework
never sees it -- it only exchanges JSON/NPZ files:

    read   request.json   (u, parameters, order, direction_map, ...)
    write  response.npz    (R_order_<p> arrays, optional `tangent`, diagnostics)

This template implements a 1-DOF cubic spring R = k u^3 - f and produces the
order-1 sensitivity right-hand-side R^(1) and the tangent using its OWN finite
differences -- so it needs nothing from us (not even OTILib). Replace `_residual`
with your model; for order >= 2, also consume `u_star_coefficients` to rebuild u*.
"""

import argparse
import json

import numpy as np


def _residual(u, params):
    k = params["k"]
    f = params["f"]
    return np.array([k * u[0] ** 3 - f], float)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--request", required=True)
    ap.add_argument("--response", required=True)
    args = ap.parse_args()

    with open(args.request, "r", encoding="utf-8") as fh:
        req = json.load(fh)

    u = np.array(req["u"], float)
    params = dict(req["parameters"])
    names = list(params.keys())
    order = int(req["order"])
    direction_map = req["direction_map"]           # {"1": [[1,0],[0,1]], ...}
    ndof = u.size

    out = {}

    # tangent dR/du by central finite differences on our own residual
    T = np.zeros((ndof, ndof))
    hu = 1e-6
    for j in range(ndof):
        up = u.copy(); up[j] += hu
        um = u.copy(); um[j] -= hu
        T[:, j] = (_residual(up, params) - _residual(um, params)) / (2 * hu)
    out["tangent"] = T

    # order-1 RHS: column per parameter direction, R^(1)[:,c] = dR/da_i
    dirs1 = direction_map.get("1", [])
    if order >= 1 and dirs1:
        R1 = np.zeros((ndof, len(dirs1)))
        for c, exps in enumerate(dirs1):
            i = int(np.argmax(exps))               # which parameter this direction differentiates
            v0 = params[names[i]]
            h = 1e-6 * max(1.0, abs(v0))
            pp = dict(params); pp[names[i]] = v0 + h
            pm = dict(params); pm[names[i]] = v0 - h
            R1[:, c] = (_residual(u, pp) - _residual(u, pm)) / (2 * h)
        out["R_order_1"] = R1

    out["diagnostics"] = np.array(json.dumps(
        {"solver": "template_blackbox", "method": "finite-difference",
         "ndof": ndof, "order": order}))

    np.savez(args.response, **out)


if __name__ == "__main__":
    main()
