# History replay: the collaborator workflow for any provider

`resasm history` is the slide-12 collaborator program for models outside the
bounded presentation scope: it turns

    OTI_UMAT.obj (+ its completed contract) + Analysis.inp + Analysis.odb + sensitivity_request.json

into `sensitivity_results.json`, `sensitivity_tables.csv` and `run_report.txt`
(plus `sensitivity_shares.csv` and a full-field `fields.npz` when requested)
without rerunning the production analysis and without the UMAT source.

    resasm history --model Analysis.inp --odb Analysis.odb \
        --material OTI_UMAT.obj --request sensitivity_request.json --out results

`--fields fields.npz` replaces `--odb` when the ODB has been exported already
(`abaqus python residual_core/replay/odb_export_npz.py -- Analysis.odb fields.npz`).
Options: `--mapping` (the completed contract; default `<object>.json` or
`Mapping.json` beside the object), `--reequilibrate`, `--verify tangent|fd`,
`--fd-steps`, `--abaqus`.

Code: `residual_core/replay/history.py` (engine), `history_material.py`
(provider link), `history_inputs.py` (deck and ODB export),
`history_outputs.py` (requests, reductions, shares, files),
`history_verify.py` (independent references), `odb_export_npz.py` (Abaqus
Python exporter), `residual_core/ui/cmd_history.py` (command and routing).
Tests: `tests/replay_history/`. Example: `examples/replay_history/j2_beam/`.

## What it supports, and what it refuses

Any provider built by `umat-oti-provider build <contract_v2.json> --out DIR`
from a UMAT-OTI version that emits `UMAT_OTI_EVAL_TOTAL` (the completed
contract lists `symbols.oti_eval_total`): any NPARAM, any NSTATV, the
provider's parameter order and PROPS slots. Older objects are refused with a
rebuild message.

Model: one `*Static` step, NLGEOM=NO (small strain), C3D8 with Abaqus's
default selective-reduced (B-bar) volumetric integration, one user material on
every element, zero-valued model-data boundaries (ENCASTRE, PINNED, XSYMM, ...
or DOF ranges), step boundaries with any value ramped over the step (the
Abaqus default for a static step), step `*Cload` ramped from zero, `*Controls`
(solver settings only). Everything else - amplitudes, `OP=NEW`, other element
types, several steps or materials, NLGEOM, distributed loads, nonzero initial
state - is refused and named in `run_report.txt` ("Unsupported feature
detected"). The UMAT receives STRESS, STATEV, STRAN (B-bar strain at the start
of the increment), DSTRAN, TIME, DTIME, PROPS, COORDS (integration point),
CELENT (element volume^(1/3)), NOEL, NPT, KSTEP=1, KINC; DROT, DFGRD0 and
DFGRD1 are the identity and TEMP/PREDEF zero, so UMATs that read those are
outside the scope.

## The mathematics

Increment n, displacements u_n, prescribed values u_c(t_n) and loads F(t_n)
independent of the parameters p:

    R_n(u_n, p) = sum_e sum_q w_q B_q^T sigma_{n,q} - F(t_n),
    sigma_{n,q} = UMAT(sigma_{n-1,q}, xi_{n-1,q}, STRAN = B_q u_{n-1}, DSTRAN = B_q (u_n - u_{n-1}), p).

Differentiating R_n = 0 on the free DOFs, with du_c/dp = 0:

    K_ff du_f/dp = - dR_f/dp |_(u_n fixed),     K = sum w B^T DDSDDE B,
    dR/dp |_(u_n fixed) = sum w B^T dsigma^A/dp,

where dsigma^A/dp comes from ONE call of `UMAT_OTI_EVAL_TOTAL` per point whose
parameter direction j is seeded with dsigma_{n-1}/dp_j, dxi_{n-1}/dp_j,
dSTRAN/dp_j = B du_{n-1}/dp_j and dDSTRAN/dp_j = -B du_{n-1}/dp_j (u_n held
fixed) plus the unit PROPS seed. First-order OTI is linear in its seeds, so
after the solve

    dsigma_n/dp = dsigma^A/dp + DDSDDE B du_n/dp,
    dxi_n/dp    = dxi^A/dp    + (dxi/dDSTRAN) B du_n/dp,
    dRF_c/dp    = dR_c/dp|_(u fixed) + K_cf du_f/dp,
    dsigma_vM/dp = (d sigma_vM / d sigma) . dsigma_n/dp.

DDSDDE and dxi/dDSTRAN come from the provider's unit DSTRAN directions in the
same call. Why a new entry point: the older `UMAT_OTI_EVAL` carries dsigma/dp
and dxi/dp but neither seeds dSTRAN/dp nor returns dxi/dDSTRAN, so the
equilibrium displacement sensitivities cannot be chained through a general
UMAT's state (the bounded J2 engine works around this for m3_j2 only).

Assembly is sparse (`scipy.sparse`, COO -> CSR with the element DOF pattern
computed once); K_ff is factorised once per increment (SuperLU) and reused for
all parameter columns; no dense K is formed or stored. The provider is called
once per increment for all integration points through a Fortran loop.

## Replay of an ODB and its checks

The recorded (single-precision) displacements drive the kinematics; the
prescribed DOFs take the deck's ramped values (the export must agree with them
to single precision). At every increment and integration point the replayed
stress, state and reactions are compared with the ODB, and the free-DOF
residual with its admissible size. Any excess is a failure (`ReplayMismatch`),
never a warning. With eps32 = 2^-23 and u32 = 8 eps32 max|U_n|:

| check | limit |
| --- | --- |
| stress | 8 eps32 max\|S_n\| + \|\|D\|\|_inf \|\|B\|\|_inf u32 |
| state | 8 eps32 max\|SDV_n\| + \|\|dxi/dDSTRAN\|\|_inf \|\|B\|\|_inf u32 |
| reactions | 8 eps32 max\|RF_n\| + \|\|K\|\|_inf u32 |
| free residual | 5e-3 mean\|nodal force\| (Abaqus's default R_n) + \|\|K_ff\|\|_inf u32 |

The first term is the storage rounding of the ODB value, the second the effect
of the rounded displacements; the 5e-3 is the force tolerance Abaqus itself
accepted. The worst error/limit ratios are printed in `run_report.txt`.

`--reequilibrate` Newton-polishes every recorded increment (starting from the
recorded state, provider tangent) to max|R_free| <= 1e-11 of the reaction
scale; the stress/state limits then widen by the size of the correction. This
matters when the ODB was converged loosely: on the FCC slide-15 cantilever
(Abaqus with the UMAT's elastic DDSDDE, default tolerances) the recorded-state
sensitivities differ from the re-equilibrated ones by 0.1-4.5 % (dsigma_vM/dm
largest); on the J2 cantilever by <= 6.6e-5. The re-equilibrated result is the
exact derivative of the discrete problem that Abaqus approximates.

## Request

The four keys of the bounded request (`outputs`, `parameters`, `domain`,
`increments`) with extensions:

```json
{
  "outputs": [
    {"name": "tip_RF2", "field": "RF", "component": 2, "reduction": "sum", "domain": {"nset": "TIP"}},
    {"name": "mises_mean", "field": "MISES", "component": 1, "reduction": "volume_mean", "domain": {"elements": "ALL"}},
    {"name": "mises_root_max", "field": "MISES", "component": 1, "reduction": "max", "domain": {"elset": "ROOTEL"}},
    {"name": "S11_e1_ip1", "field": "S", "component": 1, "reduction": "component", "domain": {"elements": [1], "points": [1]}}
  ],
  "parameters": "ALL",
  "domain": {"nodes": "ALL", "elements": "ALL"},
  "increments": "ALL",
  "weighted_shares": {"field": "MISES", "domain": {"elements": "ALL"}},
  "full_field": true
}
```

Fields U, RF (nodal), S, SDV, MISES (integration points). Reductions
component, sum, mean, volume_mean (integration-point volumes det(J) w_q), max,
min, L2. A max/min attained at several locations (mirror points of a symmetric
mesh) is accepted only when their derivatives agree; otherwise it is refused
as non-differentiable.

## Weighted sensitivity shares (slide 33)

The slide's "weighted sensitivity, so parameters with different units can be
compared" is the definition of `results/cp_residual_sensitivities.py` and
`results/fcc_crystal_results.py` (branch `cross-platform-hardening`):
W_j = |p_j dq/dp_j|, share_j = 100 W_j / sum_k W_k at every increment. For a
field q_i (von Mises stress at integration point i) the weights are aggregated
with the integration-point volumes V_i:

    W_j(n) = sum_i V_i |p_j dq_i(n)/dp_j| / sum_i V_i,      field_share_j = W_j / sum_k W_k.

`scalar_share` applies the scalar form to the volume-averaged von Mises stress.
Both are written to `sensitivity_shares.csv` and `sensitivity_results.json`.

## Outputs and report

Public: `sensitivity_results.json` (request, scope, metadata with parity,
tolerances, timings, input hashes, and the scalar results with derivatives and
weighted derivatives p dq/dp), `sensitivity_tables.csv`, `run_report.txt`,
optional `sensitivity_shares.csv` and `fields.npz` (U, dU, RF, dRF, S, dS, SDV,
dSDV, MISES, dMISES at every increment; parameter axis last). Private
(`results/private/`): the linked library, the ODB export and `run_details.json`
(per-increment residuals, parity, Newton corrections, sensitivity-solve
residuals). `run_report.txt` states separately: command executed, residual
assembled, equilibrium checked, equilibrium passed, tangent available, tangent
verified, derivative calculated, derivative verified, reference resolved,
Abaqus comparison available, unsupported feature detected, public and private
outputs separated.

`--verify tangent`: provider DDSDDE vs central FD of the ORIGINAL UMAT (the
regular path compiled into the same object) at sampled points.
`--verify fd`: additionally, the whole model is re-equilibrated in Python with
the ORIGINAL UMAT at p(1 +/- h) for a ladder of h, and the OTI derivatives are
compared with the adjacent pair of steps that forms the plateau. It costs
2 x (number of steps) x NPARAM solves - minutes for the reduced meshes,
longer at cantilever size.

## `resasm request` routing (one edit in `cmd_request.py`)

`residual_core/ui/cmd_history.py:route_request` keeps requests inside the
bounded presentation scope (pinned m3_j2 mapping, zero-valued boundaries,
concentrated loads, the four-key request with U/RF/S/SDV and
component/sum/mean/L2/max) on `presentation.run_request` and sends everything
else to this engine (`--validate` maps to `--verify fd`). To make
`resasm request` use it, replace in `residual_core/ui/cmd_request.py`, function
`register`:

```python
    parser.set_defaults(func=run)
```

with

```python
    from .cmd_history import route_request
    parser.set_defaults(func=route_request)
```

`route_request` calls `cmd_request.run` itself for bounded models, so nothing
else changes for them. (Not applied here: `cmd_request.py` belongs to the
recovery branch; the lead merges.)

## Performance (measured, Python 3.11, 24-core workstation)

| model | increments x points x parameters | replay (recorded) | re-equilibrated |
| --- | --- | --- | --- |
| J2 cantilever (slide 39) | 40 x 12,288 x 4 | 9.9 s engine, 11.1 s wall | 25.2 s engine, 26.3 s wall |
| FCC cantilever (slide 15) | 25 x 3,072 x 10 | 18.8 s engine, 19.4 s wall | 71.5 s engine, 72.1 s wall |

One cProfile of the J2 replay: SuperLU factorisation 5.3 s (0.13 s per
increment), provider 1.4 s (491,520 point evaluations), einsum kernels about
1.5 s, B-bar operators 0.9 s. Details: `docs/evidence/claude_C.md`.

## Limitations

Small strain, one static step, C3D8 (B-bar) only, ramped boundaries and
concentrated loads only, homogeneous material, parameter-independent loads,
boundaries and geometry (no load/BC/shape sensitivities), identity DROT/DFGRD,
no temperature fields. The ODB's single precision bounds how closely a replay
can reproduce it; `--reequilibrate` removes the equilibrium part of that error
but not the recorded state's own tolerance.
