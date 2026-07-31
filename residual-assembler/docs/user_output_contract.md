# User output contract

Every run writes one directory (default `resasm_output/`) split into **private**
(full numerical arrays, keep local) and **public** (safe to share).

```
resasm_output/
    private/
        metadata.json                       run info (name, order, tangent source)
        parameter_map.json                  parameter -> imaginary basis index
        dof_map.json                        ndof, free/prescribed DOF indices
        residual_real.npz                   R at the converged u
        tangent.npz                         T = dR/du
        rhs_order1.npz                      R^(1) and rhs^(1) = -R^(1)
        rhs_order2.npz                      R^(2) ... (per solved order)
        solution_sensitivities_order1.npz   U^(1) = du/da
        solution_sensitivities_order2.npz   U^(2) ... (per solved order)
        validation_full.json                residual norm, FD check, diagnostics

    public/
        summary.md                          human-readable overview
        timing.json                         wall-clock timing
        parameter_ranking.csv               parameters ranked by influence
        sensitivity_norms.csv               per-direction solution/RHS norms
        validation_summary.json             machine-readable status
```

## Privacy

- **private/** may contain proprietary numerical arrays (full residual, tangent,
  RHS, solution sensitivities). It never leaves your machine unless you share it.
- **public/** is safe to share by default. It contains only **norms, rankings,
  status, timing** — no arrays that reveal your model.

The public report **must not** include:

- source code, mesh, or element/material data,
- the full residual vector,
- the full tangent matrix,
- state variables,
- private parameter *values* (only parameter **names** and sensitivity **norms**
  appear) unless you choose to add them.

## Reading it back

```python
from resasm_user import read_report
rep = read_report("resasm_output")
print(rep["validation_summary"])   # orders solved, tangent source, status
```

or `resasm report resasm_output/`.
