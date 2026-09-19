# Which path should I use?

This page helps you choose how to use Residual_Assembler for your model. Read
the table from the top and stop at the first row that describes your
situation.

> **You do not have to provide R. You provide the ingredients, and
> Residual_Assembler builds R.**

| Your situation | Path | How you run it |
|---|---|---|
| You have a **finished Abaqus analysis** (`.inp` + `.odb`) with a UMAT, and want sensitivities of its results | **Analysis replay** *(the main workflow)* | `resasm request` with the compiled OTI provider built by UMAT-OTI ([REQUEST_INTERFACE.md](REQUEST_INTERFACE.md)) |
| Your solver hides the global residual, and you want **`R` itself**: a mesh, a converged solution, a material and boundary conditions | **Path A: assembly from ingredients** | a **recipe**: `mesh:` + `solution:` + `material:` + `parameters:` |
| Your model is **private** and you will not expose the code, but your solver can return residual coefficients | **Path B: black-box** | `residual.type: executable` |
| You can **already write `R(u, params)`** yourself (a prototype, a small custom model) | **Path C: direct residual** | `residual.type: python` |

```
Analysis replay   .inp + .odb + OTI provider + request  ->  WE replay the material and assemble R, K, dR/dp
Path A            mesh + formulation + material + fields ->  WE assemble R
Path B            your private executable                ->  it returns R^(p) coefficients
Path C            your own residual(u, params)           ->  you hand us R
```

Paths A, B and C are configured with `resasm.yml`, which has **two
dialects**, decided by one key:

- it names a **`mesh:`**: an **assembly recipe** (Path A);
- it names a **`residual:`**: a provider configuration (Path B or Path C).

Paths B and C write the same `private/` + `public/` sensitivity package. The
analysis replay writes `sensitivity_results.json`, `sensitivity_tables.csv` and
`run_report.txt`, with full fields under `private/`.

## Not sure? Ask the tool

```
resasm inspect-model model.inp
```

It reports what it found, what it **inferred for you** (element type, backend,
DOF map, integration rule, boundary conditions), the **minimum missing
ingredient**, and whether the recipe path can OTI-differentiate the model.

## Guidance per path

- **Analysis replay (main).** The material developer builds the provider once
  (`umat-oti-provider build`) and shares `OTI_UMAT.obj` with its
  `Mapping.json`; the person running the replay needs no material source.
  `resasm request` handles one bounded J2 model itself and hands every other
  readable model to `resasm history`, which covers small-strain C3D8 analyses
  with any UMAT-OTI provider, prescribed displacements and many increments.
  See [REQUEST_INTERFACE.md](REQUEST_INTERFACE.md),
  [REPLAY_HISTORY.md](REPLAY_HISTORY.md) and
  [examples/cantilevers](../examples/cantilevers/README.md).

- **Path A: assembly.** Give the ingredients; we build
  `R = F_int(u,a,q) - F_ext(a,t) + F_constraints(u,t)` by integrating element
  residuals. Almost everything is inferred from the mesh. See
  [residual_assembly_recipe.md](residual_assembly_recipe.md) and
  [abaqus_user_path.md](abaqus_user_path.md).

- **Path B: black-box.** Your model stays entirely private: your program reads
  `request.json` and writes the `R^(p)` **coefficients** (and optionally the
  tangent). We never see your code, mesh, or what the parameters mean.
  Templates: `templates/user_blackbox_residual/` (order 1) and
  `templates/user_blackbox_order2_residual/` (order 2 and above; **read
  [blackbox_order2_contract.md](blackbox_order2_contract.md) first: you must
  return Taylor coefficients, not derivatives**). Compiled C++ and Fortran code
  plugs in the same way: `templates/user_cpp_residual/`,
  `templates/user_fortran_residual/`.

- **Path C: direct residual.** Write `residual(u, params, state, time)` in
  ordinary arithmetic; we seed the parameters with OTILib and extract `R^(p)`.
  Suited to prototypes and small models. **An Abaqus user generally cannot use
  this**, because the solver never exposes `R`. Template:
  `templates/user_python_residual/`.

## Assembling R is not the same as differentiating it

For Path A, assembling `R` and **OTI-differentiating** `R` are different
capabilities:

| backend | assemble R | OTI-differentiate R through the recipe |
|---|---|---|
| `solid_c3d8_*` (C3D8), `stress_driven_c3d8`, `truss2`, `beam2` | yes | **no** |
| `nonlinear_spring1`, `nonlinear_bar1` | yes | yes |

The solid kernels store into NumPy **float** arrays (for example
`core/voigt.py::isotropic_D` allocates `np.zeros((6,6), dtype=float)`), which
cannot hold a hypercomplex number: they either raise or silently truncate to
the real part, destroying the imaginary directions. So for a C3D8 model the
recipe path can *assemble and verify* `R`, and `resasm run` **refuses** a
sensitivity it cannot compute rather than emitting a plausible-looking wrong
number.

Sensitivities of C3D8 models come from elsewhere:

- **with a UMAT:** the analysis replay (`resasm request`, `resasm history`),
  which links the OTI-transformed provider and carries the derivative arrays
  alongside the real ones, so the float kernels never have to hold an OTI
  number;
- **bounded finite-strain neo-Hookean:** `resasm sensitivity` with the example
  configuration in [examples/finite_strain_c3d8](../examples/finite_strain_c3d8/README.md)
  (first-order material parameters);
- **a private solver:** Path B.

See [assembly_minimum_information.md](assembly_minimum_information.md) for the
detail.

## `residual.type: element` is not the assembly path

Despite its name, `residual.type: element` is an **alias for `python`**: it
loads `residual.module` and calls a global `residual(u, params, state, time)`.
No element-wise assembly happens. To have the mesh assembled for you, use a
**recipe** (`mesh:`), not `residual.type: element`.

## Path letters

The documentation uses one lettering throughout: **A = assembly, B =
black-box, C = direct residual**; the analysis replay has no letter. Compiled
C++ and Fortran providers are black-box providers (Path B). When in doubt, the
unambiguous identifiers are the commands (`resasm request`, `resasm history`)
and the configuration keys: `mesh:` (assembly) versus
`residual.type: python | executable`.

See [assembly_minimum_information.md](assembly_minimum_information.md) for the
minimum inputs and [minimal_user_config.md](minimal_user_config.md) for the
provider configuration reference.
