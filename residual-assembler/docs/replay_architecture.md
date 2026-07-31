# Replay sensitivity architecture — TWO independent programs

Two **separate programs**, in **separate repositories**, shared with **different
audiences**, connected only by a versioned public contract. This is not one
application with two modules.

| | Program 1 — provider | Program 2 — consumer (**this repo**) |
|---|---|---|
| repo | `UMAT_source_transformation/oti_provider` | `Residual_Assembler/residual_core/replay` |
| audience | JHU material-model developers (hold the source) | NASA / Vanderbilt / NG / … collaborators (no source) |
| does | promote source to OTI, build the matched regular+OTI binaries, emit the manifest, verify real-part parity, package `dist/` | replay the saved analysis through the opaque binary, assemble R/K/R_,p, solve `K du/dp = -R_,p`, propagate `dq/dp` |
| never does | FE mesh, ODB, residual assembly, sensitivity solves | source transformation, OTI promotion, wrapper generation, compiling private source |

The **JHU source never leaves Program 1.** The only thing that crosses is the
compiled package (`dist/`: opaque `.so` + `material_manifest.json`) plus the
shared **contract** (`contract/`: the ABI header + manifest schema + a version
hash to catch drift). This document describes Program 2, the collaborator tool.

It is fundamentally an **offline analysis-replay and sensitivity system**, not
merely a residual assembler.

```
JHU material source ──┬── regular compilation ── material_regular.obj ─┐
                      │                                                 │  Stage 1 (collaborator)
                      └── OTI-enabled compilation ─ material_oti.obj    │  production FE analysis
                                                          │             ▼
                                material_manifest.json ───┤        saved replay record
                                                          │             │
                                                          ▼             ▼   Stage 2 (collaborator, offline)
                                              ┌───────────────────────────────────┐
   requested parameters + outputs ──────────▶│  offline residual reconstruction  │──▶ sensitivity results
                                              │  replay ▸ R,K,R_p ▸ du/dp ▸ dq/dp  │     du/dp, dq/dp
                                              └───────────────────────────────────┘
```

## Stage 1 — production analysis (unchanged)

The collaborator runs the real analysis with **`material_regular.obj`**, solving
`R(u,p)=0`. No OTI, ordinary cost, existing workflow untouched. The run saves a
**replay record** (`schemas/resasm_replay_record_v1`) — defined by what the
residual method mathematically needs, not by what an ODB happens to contain:
mesh, connectivity, element type, IP ordering, DOF map, BCs and loads, the
material parameter values used, and the converged displacement per required
increment. Path-dependent models additionally need per-increment kinematics and
the initial state.

## Stage 2 — offline replay (this package, `residual_core/replay/`)

The collaborator picks parameters and outputs and runs the offline tool with the
record + **`material_oti.obj`** + its manifest. The tool:

1. **preflights** the record (reports any mathematically-missing field) and
   **rejects a mismatched twin** (the OTI binary must be the matched pair of the
   model that produced the record — same `model_id` and regular-binary hash);
2. **replays** every integration point through the OTI binary via the C ABI
   (`abi/resasm_mat_abi.h`), obtaining real stress, the consistent tangent, and
   `d(sigma)/d(p)` for each seeded parameter — OTI arithmetic stays *inside* the
   binary;
3. assembles `K` and `R_,p` and **reuses the existing tested solver**
   (`core/field_sensitivity.py`) to solve `K du/dp = -R_,p`;
4. **propagates** to the requested outputs: `dq/dp = @q/@p + (@q/@u) du/dp`
   (`outputs.py`), and reports the replay-equilibrium and stress-parity checks.

## The math

```
R(u,p) = 0                    (converged primal, from Stage 1)
K (du/dp) = -R_,p             K = @R/@u ,   R_,p = @R/@p  (from the OTI binary)
dq/dp = @q/@p + (@q/@u) du/dp (per requested response q)
```

The **binary** provides *local* constitutive derivatives; the **framework**
turns them into *global* residual derivatives; the **solver** turns those into
*solution* sensitivities; the **output manager** turns those into the
*quantities the collaborator asked for*.

## What each side owns

| Concern | Owner |
|---|---|
| constitutive equations, OTI promotion/seeding, building the matched pair + manifest | material developer (JHU) — separate build tool |
| the C ABI spec, the record/manifest/request/result schemas | shared contract (this repo) |
| replay, residual/tangent assembly, sensitivity solve, output propagation, verification | collaborator tool (this repo) |

## Status (elastic v1)

Implemented and verified **offline, end-to-end, through the real compiled C
ABI** using a clearly-labelled REFERENCE elastic provider:

- replay stress parity vs production ~1e-16; equilibrium `||R_free||` ~1e-13;
- `du/dp` and `dq/dp` (displacement, reaction, stress component, von Mises) vs
  central finite differences of the regular model to ~1e-9;
- manifest validation, ABI dimension/seed errors, matched-twin rejection,
  record preflight, and multi-direction seeding all tested.

See `tests/framework/test_replay_elastic.py` and
`examples/replay_elastic_c3d8/`. The reference provider is **not** the real OTI
binary; see "Remaining work" below.

## Remaining work

- **Real OTI binary:** a real `material_oti.obj` built by JHU's OTI toolchain,
  exposing `mat_eval_v1`. The reference provider computes elastic derivatives
  analytically; it exists only to exercise the pipeline.
- **Finite strain:** the record must carry `F0,F1` per IP; the solver's finite-
  strain assembly (spatial B, geometric term, objective-rate tangent) is a
  separate milestone.
- **Path-dependent (crystal plasticity):** replay must march every increment
  carrying `state_in`/`dstate_dseed_in` → `state_out`/`dstate_dseed_out`; the
  record must store per-increment kinematics and the initial state
  (preflight already demands these when `nstatev>0`).
