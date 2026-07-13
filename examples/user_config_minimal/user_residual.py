def residual(u, params, state=None, time=None):
    k = params["k"]
    f = params["f"]
    return [k * u[0] ** 3 - f]


def tangent(u, params, state=None, time=None):
    k = params["k"]
    return [[3.0 * k * u[0] ** 2]]
