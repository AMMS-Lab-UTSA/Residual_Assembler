# Which path should I use?

> **You should not have to provide R. You provide the ingredients, and
> Residual_Assembler builds R.**

Start at the top. Only fall through if the row above does not apply. All three
paths produce the same `private/` + `public/` sensitivity package.

| Your situation | Path | How you configure it |
|---|---|---|
| **Your solver hides the global residual** (Abaqus and friends) — but you have a mesh, a converged solution, a material, and BCs | **Path A — assembly from ingredients** *(the main path)* | a **recipe**: `mesh:` + `solution:` + `material:` + `parameters:` |
| Your model is **private** and you will not expose the code, but your solver can return residual coefficients | **Path B — black-box** | `residual.type: executable` |
| You can **already write `R(u, params)`** yourself (toy, prototype, small custom model) | **Path C — direct residual** *(shortcut)* | `residual.type: python` |

```
Path A   mesh + formulation + material + solution + stimuli  ->  WE assemble R
Path B   your private executable                             ->  it returns R^(p) coefficients
Path C   your own residual(u, params)                        ->  you hand us R
```

`resasm.yml` has **two dialects**, decided by one key:

- names a **`mesh:`**     → an **assembly recipe** (Path A)
- names a **`residual:`** → a provider config (Path B / Path C)

## Not sure? Ask the tool

```
resasm inspect-model model.inp
```

It reports what it found, what it **inferred for you** (element type, backend, DOF
map, integration rule, BCs), the **minimum missing ingredient**, and — crucially —
whether the model can be OTI-differentiated at all.

## Guidance per path

- **Path A — assembly (main).** Give the ingredients; we build
  `R = F_int(u,a,q) - F_ext(a,t) + F_constraints(u,t)` by integrating element
  residuals. Almost everything is inferred from the mesh. See
  [residual_assembly_recipe.md](residual_assembly_recipe.md) and
  [abaqus_user_path.md](abaqus_user_path.md).

- **Path B — black-box.** Your model stays entirely private: your program reads
  `request.json` and writes the `R^(p)` **coefficients** (and optionally the
  tangent). We never see your code, mesh, or what the parameters mean.
  Templates: `templates/user_blackbox_residual/` (order 1),
  `templates/user_blackbox_order2_residual/` (order ≥ 2 — **read
  [blackbox_order2_contract.md](blackbox_order2_contract.md) first: you must return
  Taylor coefficients, not derivatives**). Compiled C++/Fortran code plugs in the
  same way: `templates/user_cpp_residual/`, `templates/user_fortran_residual/`.

- **Path C — direct residual (shortcut).** Write
  `residual(u, params, state, time)` in ordinary arithmetic; we seed the parameters
  with OTILib and extract `R^(p)`. Fine for toys and prototypes. **An Abaqus user
  generally cannot use this** — the solver never exposes R. Template:
  `templates/user_python_residual/`.

## ⚠ Honest routing caveat — read before choosing Path A *for sensitivities*

Assembling `R` and **differentiating** `R` are different capabilities.

| backend | assemble R | OTI-differentiate R |
|---|---|---|
| `solid_c3d8_*` (C3D8), `stress_driven_c3d8`, `truss2`, `beam2` | yes | **no** |
| `nonlinear_spring1`, `nonlinear_bar1` | yes | yes |

The solid kernels store into numpy **float** arrays (e.g. `core/voigt.py::isotropic_D`
allocates `np.zeros((6,6), dtype=float)`), which cannot hold a hypercomplex number —
it either raises or silently truncates to the real part, destroying the imaginary
directions.

So **today**, for a C3D8 model you can *assemble and verify* `R` through Path A, but
you cannot yet get OTI sensitivities through it. `resasm run` **refuses and says so**
rather than emitting a plausible-looking wrong number. For sensitivities on such a
model use **Path B**, or an **OTI-transformed UMAT** (the companion source-transformation
project) once wired in. Making the solid kernels OTI-safe is a known engineering task,
not a physics limit — see
[assembly_minimum_information.md](assembly_minimum_information.md).

## ⚠ `residual.type: element` is NOT the assembly path

Despite its name, `residual.type: element` is currently just an **alias for
`python`**: it loads `residual.module` and calls a global
`residual(u, params, state, time)`. No element-wise assembly happens.

If you want us to assemble from a mesh, use a **recipe** (`mesh:`), not
`residual.type: element`.

## Path-letter note

Some older pages use "Path A/B/C" for the three *provider kinds*
(Python / executable / compiled). The canonical lettering is the one on this page:
**A = assembly, B = black-box, C = direct residual.** When in doubt, the unambiguous
identifiers are the config keys: `mesh:` (assembly) vs
`residual.type: python | executable`.

See [assembly_minimum_information.md](assembly_minimum_information.md) for the
minimum inputs and [minimal_user_config.md](minimal_user_config.md) for the
provider-config reference.
