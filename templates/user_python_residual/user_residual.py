# Edit this file: implement your residual (and, optionally, your tangent).
#
# The residual is evaluated two ways by the framework:
#   * with plain floats (real check), and
#   * with OTILib scalars (to extract sensitivities) -- so use ordinary
#     arithmetic and it will carry the derivatives automatically.
#
# u      : the unknowns (indexable; u[i] may be a float or an OTI scalar)
# params : a dict of design parameters (values may be float or OTI scalars)
# Return : the residual R as a length-`unknowns` list/array.


def residual(u, params, state=None, time=None):
    k = params["k"]
    f = params["f"]
    return [k * u[0] ** 3 - f]        # cubic spring:  R = k u^3 - f


def tangent(u, params, state=None, time=None):
    # optional. dR/du at the real solution (plain floats here).
    k = params["k"]
    return [[3.0 * k * u[0] ** 2]]
