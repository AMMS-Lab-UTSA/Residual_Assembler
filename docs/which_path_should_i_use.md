# Which path should I use?

Pick one — they all produce the same private/ + public/ sensitivity package.

| Your situation | Use | `residual.type` |
|---|---|---|
| Your residual is easy to expose in Python | **Python provider** | `python` |
| Your code already exists in compiled C++/Fortran | **C++/Fortran provider** | `executable` (wrap it) |
| You cannot link or expose internals at all | **Black-box executable** | `executable` |
| You want *us* to assemble your residual from a mesh | **Element backend** | `element` |

Guidance:

- **Python provider** — write `residual(u, params, state, time)` with ordinary
  arithmetic; the framework seeds parameters with OTILib and extracts `R^(p)`
  automatically. Simplest if you can express the residual in Python.
  Template: `templates/user_python_residual/`.

- **C++/Fortran provider** — you already have compiled code. Write a
  `template<class Scalar> eval_residual(...)` (C++) or an OTI-typed subroutine
  (Fortran), and wrap it as an executable behind the request/response contract.
  Templates: `templates/user_cpp_residual/`, `templates/user_fortran_residual/`.

- **Black-box executable** — your model stays entirely private. Your program
  reads `request.json` and writes `response.npz` with the `R^(p)` coefficients
  (and optionally the tangent). The framework never sees your code, mesh, or
  parameters' meaning. Template: `templates/user_blackbox_residual/`.

- **Element backend** — only if you *want* the framework to assemble the global
  residual from your mesh/element data. Most private-code users do **not** need
  this — a global residual (Python) or an executable that returns `R^(p)`
  (black-box) keeps everything local.

See [user_input_contract.md](user_input_contract.md) for the minimum inputs and
[minimal_user_config.md](minimal_user_config.md) for the config.
