"""Scalar reductions of extracted real fields and their total derivatives."""

import numpy as np


def reduce_scalar(values, derivatives, reduction):
    values = np.asarray(values, dtype=float).reshape(-1)
    derivatives = np.asarray(derivatives, dtype=float)
    if derivatives.ndim < 2 or derivatives.size == 0:
        raise ValueError("derivatives must have a final parameter axis")
    derivatives = derivatives.reshape(-1, derivatives.shape[-1])
    if not values.size or derivatives.shape[0] != values.size:
        raise ValueError("empty domain or field/derivative shape mismatch")
    if not np.isfinite(values).all() or not np.isfinite(derivatives).all():
        raise ValueError("nonfinite field or derivative")
    if reduction == "component":
        if values.size != 1:
            raise ValueError("component reduction requires exactly one scalar location")
        return float(values[0]), derivatives[0]
    if reduction == "sum":
        return float(values.sum()), derivatives.sum(axis=0)
    if reduction == "mean":
        return float(values.mean()), derivatives.mean(axis=0)
    if reduction == "L2":
        norm = float(np.linalg.norm(values))
        if norm == 0:
            raise ValueError("L2 at zero is nondifferentiable; select a component instead")
        return norm, values @ derivatives / norm
    if reduction == "max":
        maximum = float(values.max())
        active = np.flatnonzero(np.isclose(values, maximum, rtol=1e-12, atol=1e-14))
        if len(active) != 1:
            raise ValueError("max has tied active values; derivative is not uniquely established")
        return maximum, derivatives[active[0]]
    raise ValueError("unsupported reduction %r; use component, sum, mean, L2, max" % reduction)