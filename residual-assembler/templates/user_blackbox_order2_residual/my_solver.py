#!/usr/bin/env python
"""Black-box residual solver — ORDER-2 REFERENCE.

Your private model stays inside this executable. The framework only exchanges
files with it:

    read   request.json     (u, parameters, order, direction_map, u_star_coefficients, ...)
    write  response.npz     (R_order_<p>, optional tangent, optional diagnostics)

===============================================================================
 THE ONE RULE THAT MATTERS AT ORDER >= 2
===============================================================================
 The arrays you return are OTI **TAYLOR COEFFICIENTS**, *not* partial derivatives.

     coefficient = derivative / prod_i factorial(kappa_i)

 For the direction [2, 0] (i.e. d2/dk2) you must return

     (1/2!) * d2R/dk2        NOT   d2R/dk2

 The framework applies the recovery factor (prod_i kappa_i!) itself when it
 exports `*_derivatives` and writes the public report. If you pre-multiply, your
 second derivatives come out 2x too large (6x at order 3, ...).

 The safest way to get this right is what this file does: DO NOT hand-derive the
 coefficients. Write your residual once in ordinary arithmetic and evaluate it in
 a truncated Taylor algebra, then read the coefficients straight off. That is
 exactly what the framework does internally with OTILib.
===============================================================================

Model (edit `residual()` and `dR_du()` only):

    R(u, k, f) = k^2 * u^3 - f

    chosen so the order-2 residual coefficients are VISIBLY NON-ZERO (unlike
    R = k*u^3 - f, whose second parameter-derivatives all vanish).

Reference values at k=2, f=32 (=> u=2, T=dR/du=48):

    R^(1) coefficients : [dR/dk, dR/df]         = [32, -1]
    R^(2) coefficients : [k^2,   k*f,   f^2]    = [-40/3, 1/3, 1/96]
    => U^(2) coefficients                       = [5/18, -1/144, -1/4608]
    => recovered derivatives (x 2!, 1!, 2!)     = [5/9,  -1/144, -1/2304]
       which are exactly d2u/dk2, d2u/dkdf, d2u/df2 of u(k,f) = (f/k^2)^(1/3).
"""

import argparse
import itertools
import json
import math

import numpy as np


# =============================================================================
# 1. A minimal truncated Taylor algebra (order <= 2, m bases).
#    Coefficients are TAYLOR coefficients: the coefficient of the monomial
#    e_1^k1 ... e_m^km is  (1/prod ki!) * d^p R / da_1^k1 ... da_m^km.
#    Products convolve:  (A*B)_kappa = sum_{alpha+beta=kappa} A_alpha B_beta.
#    You should not need to touch this.
# =============================================================================
class T2(object):
    """Truncated Taylor number: real part + order-1 + order-2 coefficients."""

    __slots__ = ("m", "re", "o1", "o2")

    def __init__(self, m, re=0.0, o1=None, o2=None):
        self.m = m
        self.re = float(re)
        self.o1 = list(o1) if o1 is not None else [0.0] * m
        # o2 keyed by (i, j) with i <= j
        self.o2 = dict(o2) if o2 is not None else {}

    @staticmethod
    def const(m, v):
        return T2(m, v)

    @staticmethod
    def seed(m, v, i):
        """Real value v, perturbed along basis i (0-based): v + 1*e_i."""
        o1 = [0.0] * m
        o1[i] = 1.0
        return T2(m, v, o1)

    def get2(self, i, j):
        return self.o2.get((min(i, j), max(i, j)), 0.0)

    def __add__(self, other):
        b = other if isinstance(other, T2) else T2.const(self.m, other)
        o2 = dict(self.o2)
        for kk, vv in b.o2.items():
            o2[kk] = o2.get(kk, 0.0) + vv
        return T2(self.m, self.re + b.re,
                  [x + y for x, y in zip(self.o1, b.o1)], o2)

    __radd__ = __add__

    def __neg__(self):
        return T2(self.m, -self.re, [-x for x in self.o1],
                  {k: -v for k, v in self.o2.items()})

    def __sub__(self, other):
        return self + (-(other if isinstance(other, T2)
                         else T2.const(self.m, other)))

    def __rsub__(self, other):
        return (-self) + other

    def __mul__(self, other):
        b = other if isinstance(other, T2) else T2.const(self.m, other)
        m = self.m
        re = self.re * b.re
        o1 = [self.re * b.o1[i] + self.o1[i] * b.re for i in range(m)]
        o2 = {}
        for i in range(m):
            for j in range(i, m):
                # real x order2  +  order2 x real  +  order1 x order1
                v = self.re * b.get2(i, j) + self.get2(i, j) * b.re
                if i == j:
                    v += self.o1[i] * b.o1[i]
                else:
                    v += self.o1[i] * b.o1[j] + self.o1[j] * b.o1[i]
                if v:
                    o2[(i, j)] = v
        return T2(m, re, o1, o2)

    __rmul__ = __mul__

    def __pow__(self, n):
        out = T2.const(self.m, 1.0)
        for _ in range(int(n)):
            out = out * self
        return out


def directions(m, p):
    """Canonical order-p exponent vectors — IDENTICAL ordering to the framework
    (`combinations_with_replacement` over the 1-based basis indices)."""
    out = []
    for combo in itertools.combinations_with_replacement(range(m), p):
        e = [0] * m
        for b in combo:
            e[b] += 1
        out.append(e)
    return out


# =============================================================================
# 2. YOUR MODEL — this is the only part you edit.
#    Write it in ORDINARY arithmetic. It must work for plain floats AND for T2
#    numbers, so use only + - * ** (no float() casts, no numpy ufuncs).
# =============================================================================
def residual(u, k, f):
    """R(u, k, f) = k^2 u^3 - f    (u is a length-ndof list)."""
    return [k * k * u[0] ** 3 - f]


def dR_du(u, k, f):
    """Tangent T = dR/du at the converged real solution (plain floats)."""
    return [[3.0 * k * k * u[0] ** 2]]


# =============================================================================
# 3. Plumbing — evaluate the model in the Taylor algebra and emit COEFFICIENTS.
# =============================================================================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--request", required=True)
    ap.add_argument("--response", required=True)
    args = ap.parse_args()

    with open(args.request, "r") as fh:
        req = json.load(fh)

    u = [float(x) for x in req["u"]]
    ndof = len(u)
    params = req["parameters"]
    # seed_directions: name -> 1-based basis index. Defines the basis order.
    seeds = req.get("seed_directions") or {n: i + 1 for i, n in enumerate(params)}
    names = sorted(seeds, key=lambda n: seeds[n])          # basis order
    m = len(names)
    p = int(req["order"])                                  # the order asked for NOW
    dmap = req["direction_map"][str(p)]                    # exponent vectors, in column order

    # --- build u*: real part u, plus the LOWER-order sensitivities the framework
    #     already solved and handed back (u_star_coefficients). The order-p
    #     coefficients of u* are ZERO -- they are exactly what is being solved for.
    usc = req.get("u_star_coefficients") or {}
    ustar = [T2(m, u[j]) for j in range(ndof)]
    if p >= 2 and "1" in usc:
        U1 = usc["1"]                                      # (ndof, m), canonical order-1 columns
        for j in range(ndof):
            # order-1 column i corresponds to basis i (see directions(m, 1))
            ustar[j].o1 = [float(U1[j][i]) for i in range(m)]

    # --- seed the parameters: a_i* = a_i + 1*e_i
    pstar = {n: T2.seed(m, float(params[n]), seeds[n] - 1) for n in names}

    # --- evaluate the residual ONCE in the Taylor algebra
    R = residual(ustar, pstar[names[0]], pstar[names[1]]) if m == 2 else None
    if R is None:                                          # generic fallback
        R = residual(ustar, *[pstar[n] for n in names])

    # --- read the order-p TAYLOR COEFFICIENTS straight off. No factorials here:
    #     the framework multiplies by prod(kappa_i!) when it recovers derivatives.
    R_p = np.zeros((ndof, len(dmap)))
    for col, exps in enumerate(dmap):
        idx = [i for i, e in enumerate(exps) for _ in range(e)]   # e.g. [2,0] -> [0,0]
        for j in range(ndof):
            if p == 1:
                R_p[j, col] = R[j].o1[idx[0]]
            elif p == 2:
                R_p[j, col] = R[j].get2(idx[0], idx[1])
            else:
                raise SystemExit("this reference template implements order <= 2")

    out = {"R_order_%d" % p: R_p,
           "tangent": np.asarray(dR_du(u, *[float(params[n]) for n in names]), float),
           "diagnostics": np.array(json.dumps({
               "solver": "template_blackbox_order2",
               "method": "truncated Taylor algebra (order 2)",
               "convention": "TAYLOR COEFFICIENTS (not derivatives)",
               "ndof": ndof, "order": p, "basis": names}))}
    np.savez(args.response, **out)


if __name__ == "__main__":
    main()
