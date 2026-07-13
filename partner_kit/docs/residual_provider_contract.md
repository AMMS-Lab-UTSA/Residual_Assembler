# Residual Provider Contract

Implement **one** of three provider levels locally. All are scalar-generic: a
parameter value may be a plain float or a hypercomplex seed (`Dual1`) — write the
residual with ordinary arithmetic and the derivative rides along in the imaginary
part. Tangent is optional and can be supplied several ways.

The kit never requires you to expose implementation details. Declared metadata is
limited to: a provider `name`, the ordered `parameters`, and `ndof`.

## Level 1 — Element residual provider

```python
class ElementResidualProvider:
    name: str
    parameters: tuple[str, ...]
    ndof: int
    def parameter_values(self) -> dict[str, float]: ...
    def free_mask(self) -> np.ndarray | None: ...          # True = free DOF
    def elements(self) -> Iterable[Any]: ...
    def element_dof_map(self, element_id) -> list[int]: ...  # global DOF scatter
    def element_state(self, element_id): ...                 # or None
    def eval_element_residual(self, element_id, element_dofs, element_state,
                              parameters, time, dtime):
        # -> (element_residual, element_tangent_optional, updated_state_optional)
```

The kit scatters element residuals into the global system using
`element_dof_map`. Return an element tangent to enable Newton solves and the
`backend-assembled` tangent.

## Level 2 — Global residual provider

```python
class GlobalResidualProvider:
    name, parameters, ndof, parameter_values, free_mask ...
    def eval_global_residual(self, global_solution, state, parameters, time, dtime):
        # -> residual vector (length ndof)
```

## Tangent (optional, any one)

```python
def get_tangent(self, U, state, parameters, time, dtime) -> np.ndarray        # dense T
def apply_tangent(self, U, state, parameters, x, time, dtime) -> np.ndarray    # matrix-free T·x
# ...or pass an external file:   config: { "tangent_file": "tangent.npz" }
```

If no tangent is available the kit still generates and exports $\mathbf{R}^{(1)}$;
you (or we) can solve $\mathbf{T}\,\mathbf{U}^{(1)} = -\mathbf{R}^{(1)}$ later.

## Level 3 — Black-box executable (no linking)

You provide a command; the kit writes `request.json`, runs it, and reads
`response.json`. Your model stays entirely behind the executable.

```
your_solver --input request.json --output response.json
```

**Request** (`resasm-partner-request/1`):

```json
{ "order": 1,
  "solution": [ ... ],
  "parameters": { "k": 2.0, "f": 16.0 },
  "seed_directions": { "k": 1 },
  "time": [0.0, 0.0], "dtime": 0.0 }
```

**Response** (`resasm-partner-response/1`):

```json
{ "status": "ok",
  "residual_coefficients": [[ ... ]],   // ndof x m  == R^(p) columns per seeded param
  "residual_real": [ ... ],             // optional R (order 0)
  "tangent": [[ ... ]],                 // optional dense T (or null)
  "arrays_npz": "response.npz",         // optional sidecar for large arrays
  "diagnostics": { ... }, "message": "" }
```

`residual_coefficients` column `j` is the derivative of the residual with respect
to the parameter whose `seed_directions` basis index is `j+1` (order 1). The kit
negates to form `rhs = -R^(p)` and solves with `tangent` if provided.

Config to use it:

```json
{ "blackbox": { "command": ["./your_solver"], "ndof": 1,
                "parameters": ["k"], "parameter_values": {"k": 2.0, "f": 16.0},
                "free_mask": [true] },
  "parameters": ["k"], "order": 1, "solution": [2.0], "output_dir": "out" }
```

## Scalar-type seam (OTI/HYPAD later)

Order 1 uses a first-order dual number. To go higher order without rewriting your
residual, keep it **generic in the scalar type**:

- **C++**: `template<class Scalar> void eval_residual(..., const Scalar* params, Scalar* R)`
  — instantiate with `double`, `Dual1`, `OTI`, or a HYPAD scalar. See
  [../include/residual_provider.hpp](../include/residual_provider.hpp).
- **Fortran**: use a source-transformation tool, an operator-overloaded derived
  type, a generated wrapper, or (for validation only) a finite-difference
  fallback. See [../templates/fortran_uel_wrapper](../templates/fortran_uel_wrapper).
- **Black-box**: return `R^(p)` for the requested `order`; your internal method is
  your choice.
