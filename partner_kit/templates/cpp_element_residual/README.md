# Template — C++ element residual (templated Scalar)

Implement your element residual once, generic in the scalar type, so it evaluates
with `double` (real), `resasm::Dual1` (order-1 sensitivity), or a future
OTI/HYPAD scalar — no rewrite.

```cpp
#include "residual_provider.hpp"     // from ../../include

template <class Scalar>
void my_element_residual(int /*eid*/,
                         const Scalar* u, std::size_t /*ndofs*/,
                         const Scalar* p, std::size_t /*np*/,
                         double /*t*/, double /*dt*/,
                         Scalar* r) {
  // Example: cubic axial bar between local dofs 0,1 with stiffness p[0]
  Scalar d = u[1] - u[0];
  Scalar N = p[0] * d * d * d;       // ordinary arithmetic -> derivative rides along
  r[0] = -N;
  r[1] =  N;
}
```

- Seed a parameter with `resasm::seed(value)` and read `resasm::imag_part(r[i])`
  to get `dR/dp` — that is a column of `R^(1)`.
- To expose this to the kit, either (a) assemble the global `R^(1)` yourself and
  write `response.json` with `resasm::write_response_json`
  (see [../../include/sensitivity_output.hpp](../../include/sensitivity_output.hpp)),
  or (b) drive it from Python via the element provider contract.

Nothing here reveals your mesh or material model. Only the numbers you choose to
write leave the process.

## Build

```bash
c++ -std=c++17 -I../../include your_program.cpp -o your_solver
```
