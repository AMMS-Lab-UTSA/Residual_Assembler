# Which layer did it come from

A residual assembled from a converted UMAT passes through nine places where it
can go wrong, and they need nine different fixes. "The residual is wrong by
12%" tells nobody which of nine people to ask, so this repository answers the
question directly: `residual_core.diagnose(fixture, ...)` returns a
`Diagnosis` whose `blame` is the layer the failure came from and whose
`findings` carry the numbers each verdict was reached on.

```python
from residual_core import diagnose
from residual_core.materials.verified_fixture import load

found = diagnose(load("tests/fixtures/verified/neohookean_umat--10759f1ffd.json"),
                 assemble=my_assembly, reference=an_independent_one)
print(found.blame)        # None, or one of the nine layer names
print(found.report())     # every verdict, with its measurement
found.as_dict()           # the same thing for a program
```

## The nine layers

| Layer | Who fixes it | What the check measures |
|---|---|---|
| `umat` | the author's routine, as published | every carried value is a number, and the stress moves when the strain does |
| `transformation` | the OTI source transformation | the converted build's stress and state against the original's, over the whole window |
| `constitutive_derivative` | the derivative the converted build extracts | whether the reported tangent differentiates the reported stress, in either of the two readings of DDSDDE |
| `state_sensitivity` | the internal-state derivative and how it is carried | a supplied `dR/dq` against a centred difference of `R`, over a plateau of step sizes |
| `mapping` | the convention this assembler reads the UMAT's arrays in | `ndi + nshr == ntens`, the tangent's shape and shear diagonal, and how much shear the window actually exercises |
| `element_integration` | the quadrature and integration-point order here | a constant stress against its surface tractions, a linearly varying one against an independent target, a rigid translation against zero |
| `residual_assembly` | the integration of `B^T sigma` in this repository | the assembled force against a second, independently written integration |
| `boundary_conditions` | the prescribed/free partition here | the partition covers every DOF without overlap, and how many DOFs a step-blind read would get wrong |
| `global_dof_mapping` | the element-to-global scatter here | the map is a dense bijection, and every element gathers the rows its own connectivity resolves to |

They run in that order and the run **stops at the first failure**: everything
downstream of a NaN is meaningless, and reporting nine verdicts would bury the
one that matters.

## Three answers, not two

A layer reports `holds`, `fails`, or `not_established`. The third is not a
softer version of the first. It means this diagnosis was given nothing that
could decide the layer, and the finding says what would:

```
mapping: NOT ESTABLISHED -- 6 components as 3 direct and 3 engineering shear,
  and the tangent is shaped and signed accordingly. But the carried window
  applies effectively NO shear strain -- the largest shear component is
  2.03e-20 of the largest extension -- so these numbers cannot tell
  engineering shear from tensorial: the two differ by a factor of two in
  components that are all zero here.
  what would establish it: a fixture window inside the deck's shear step,
  where a factor of two in the shear convention changes the stress
```

`Diagnosis.ok` is true when nothing that was checked failed.
`Diagnosis.complete` is true only when every layer was checked AND held. Use
`complete` when you need a clean bill; use `ok` when you need to know whether
to stop.

## The two readings of DDSDDE

Abaqus's finite-strain material Jacobian is the tangent of the Jaumann rate of
**Kirchhoff** stress over `J`, so the Cauchy increment it predicts is

    Dsigma = D : Deps  -  sigma tr(Deps)

and a small-strain routine's tangent is just `D : Deps`. **Which one a UMAT
returns is a property of the routine, not of the step's `NLGEOM` flag.**
Measured on the committed fixtures, at the mid-point of each increment so that
a chord across a nonlinear increment is not mistaken for a wrong tangent:

| fixture | NLGEOM | plain reading | Jaumann reading | verdict |
|---|---|---|---|---|
| `irfancn umat_elastic` | YES | 3.4e-16 | 5.5e-03 | `material` |
| `AlexanderJFDR NeoHookean` | YES | 8.5e-03 | 6.9e-08 | `jaumann` |
| `Sina-Taghizadeh CompresibleNeoHookean` | YES | 8.5e-03 | 6.9e-08 | `jaumann` |
| `keisuke58 biofilm visco` | YES | 8.4e-03 | 1.0e-08 | `jaumann` |
| `CAEAssistant isotropic elasticity` | NO | 4.0e-16 | 5.5e-03 | `material` |

An assembler that picks by the flag builds the wrong stiffness for the first
row. `residual_core.materials.verified_fixture.tangent_convention(fixture)`
measures it and says how decisively; a window whose strain increments are
isochoric cannot separate the two and reports so rather than choosing.

## Provenance comes before the layers

A fixture is evidence about the transformation that produced it. `load()`
refuses one whose `transform_fingerprint` is not the current store's
(`b0d27ee53c630500`, read from `corpus_run/pass11/results/store_verification.jsonl`,
237 entries) — because verifying today's assembler against a fixture frozen
under an older transformation is verifying it against somebody else's run.
Reading one anyway is possible and deliberate: `load(path, fingerprint=None)`.

## What the frozen set can and cannot decide

Every committed fixture's window sits inside its deck's **uniaxial** step. That
window decides the Voigt ORDER and the direct components; it cannot decide the
shear convention, and the `mapping` layer says so. The shear question is
settled separately, against the deck's second step, whose affine motion has a
strain known in closed form — see
`tests/framework/test_a_verified_deck_drives_the_global_assembly.py`.

The integration-point ORDERING is internally consistent and is not verified
against Abaqus's export index: every offline check uses a uniform field, which
is provably blind to it. Filed as
`abaqus_queue/requests/A3_integration_point_ordering.json`.
