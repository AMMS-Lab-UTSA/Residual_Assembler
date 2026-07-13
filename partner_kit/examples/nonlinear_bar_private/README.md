# Example — private bar (element provider)

A two-element nonlinear axial bar chain (nodes 0-1-2), node 0 fixed and node 2
given a prescribed displacement. Each element has its own stiffness parameter
(`k1`, `k2`). Implemented in [provider.py](provider.py) as an
`ElementResidualProvider` — the kit scatters element residuals it never has to
understand, and seeds `k1`/`k2` independently to build `R^(1)`.

## Run

```bash
python -m partner_kit.python.partner_cli --config config.json
```

## Expected output

```
provider        : private_bar (PrivateBar)
parameters      : ['k1', 'k2']
tangent source  : backend-assembled
R^(1) shape     : (3, 2)
hypercomplex ok : True
  dU/d(k1        ) = +1.666667e-01
  dU/d(k2        ) = +1.666667e-01
validation      : True
```

- The tangent is **assembled from the element tangents** the provider returns.
- Two design parameters ⟹ `R^(1)` has 2 columns; `dU/dk1` and `dU/dk2` are the
  free-DOF sensitivities (node 1), verified against finite differences locally.

## What stays private

The mesh (3 nodes, 2 elements), the connectivity, and the cubic law all live in
`provider.py` on the partner machine. The public report shares only the
per-parameter sensitivity norms and the validation summary.
