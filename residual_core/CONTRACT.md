# Residual Assembler — Shared Contract (C3D8 + CP conventions)

> **Post-refactor note.** This document pins the numerical/data conventions for
> the **C3D8 solid + crystal-plasticity** path. After the framework refactor the
> files moved: `abaqus_inp_parser.py` → `io/abaqus_inp_parser.py`;
> `c3d8_residual.py` (the numerics) → `formulations/c3d8_kernel.py`;
> `extract_abaqus_fields.py` → `io/abaqus_odb_export.py`; the CP tests →
> `tests/cp_c3d8_umat/`. The conventions below are unchanged and are still
> authoritative for that backend. For the framework-level contracts (the
> `Formulation`/`Material` interfaces and how to add a backend) see
> [docs/formulation_contract.md](docs/formulation_contract.md) and
> [docs/architecture.md](docs/architecture.md).

This file pins down the **data model, APIs, and mechanical/numerical conventions**
that every module in the C3D8/CP path must agree on. All parallel work builds
against this contract so the pieces interlock. Do not silently deviate — if a
convention here is wrong, flag it, don't work around it.

Scope of the first milestone: **standard continuum C3D8 element + crystal-plasticity
UMAT stress update + residual assembly**, verified against Abaqus. No OTI/HYPAD,
no UEL path, no new CP physics.

---

## 0. Ground truth for the first target

Primary example:
`sources/permissive/ngrilli_Oxford_Crystal_Plasticity/ExampleInputFiles/HCPnoTwin/Compression111.inp`

- 216 nodes, 125 `C3D8` elements, one material `CPuranium`.
- `*User Material, constants=11`, `*Depvar 125`.
- `*Step, nlgeom=YES` → **finite strain**.
- UMAT: `sources/permissive/ngrilli_Oxford_Crystal_Plasticity/umat.for` (+ `kmat.f` etc.).

### The stress convention (verified from source, do not re-guess)
- `kmat.f:1023`  `usvars(47+i) = xstressdef(i)` stores the **Cauchy (true) stress**.
- `umat.for:388` `stress(K) = STATEV(47+K)` copies it to Abaqus `STRESS`.
- Voigt order is **Abaqus order**: `(σ11, σ22, σ33, σ12, σ13, σ23)` (`kmat.f:78`).
- Therefore, for these `nlgeom=YES` C3D8 elements, Abaqus assembles the internal
  force in the **current (deformed) configuration**:

  ```
  f_int,e = ∫_v  B_spatial(x)^T · σ_cauchy  dv           (x = X + u,  dv = det(J_current)·w)
  ```

  This is the convention the residual assembler must reproduce to match Abaqus
  reaction forces. A separate small-strain path (reference config, `F = I`) exists
  only for analytic patch-test verification.

---

## 1. Voigt / DOF conventions (MANDATORY, identical everywhere)

- **Stress/strain Voigt order:** `(11, 22, 33, 12, 13, 23)` — Abaqus order.
- **Strain uses engineering shear:** `ε_voigt = (ε11, ε22, ε33, 2ε12, 2ε13, 2ε23)`.
  With this pairing, internal virtual work is `δε_voigt · σ_voigt` with the true
  (not doubled) shear stresses `σ12, σ13, σ23`. `B` is built to produce engineering
  shear strains.
- **Element DOF ordering:** node-major, 24 DOFs per C3D8:
  `[u1x,u1y,u1z, u2x,u2y,u2z, …, u8x,u8y,u8z]`, where node k is the k-th entry of
  the element connectivity list.
- **Global DOF ordering:** nodes sorted ascending by node id → index `i` (0-based);
  node id `nid` → DOFs `[3*i, 3*i+1, 3*i+2]`. A `node_id_to_index: dict[int,int]`
  map is part of the assembled-model output.

## 2. C3D8 element geometry conventions (MANDATORY)

Isoparametric natural coords ξ ∈ [-1,1]^3. **Node natural coordinates, in Abaqus
C3D8 connectivity order (node 1..8):**

```
n1 (-1,-1,-1)   n2 (+1,-1,-1)   n3 (+1,+1,-1)   n4 (-1,+1,-1)
n5 (-1,-1,+1)   n6 (+1,-1,+1)   n7 (+1,+1,+1)   n8 (-1,+1,+1)
```

Trilinear shape functions `N_a(ξ) = 1/8 (1+ξ1·ξ1a)(1+ξ2·ξ2a)(1+ξ3·ξ3a)`.

**Gauss integration:** full 2×2×2, points at ±1/√3, all weights = 1. The 8 Gauss
points MUST be ordered to match **Abaqus C3D8 integration-point (SDV/S) output
order**, which is (fastest index first = ξ1):

```
IP1 (-g,-g,-g)  IP2 (+g,-g,-g)  IP3 (-g,+g,-g)  IP4 (+g,+g,-g)
IP5 (-g,-g,+g)  IP6 (+g,-g,+g)  IP7 (-g,+g,+g)  IP8 (+g,+g,+g)      g = 1/√3
```

This ordering matters only when pairing an assembler Gauss point with an
Abaqus-exported integration-point stress. Expose it as a named constant
`ABAQUS_C3D8_GAUSS` (8×3 coords + weights) so it can be audited/changed in one place.

## 3. Parser output — `abaqus_inp_parser.py`

`parse_inp(path) -> AbaqusModel`. `AbaqusModel` is a dataclass (JSON-serializable
via a `.to_dict()`), fields:

- `nodes: dict[int, tuple[float,float,float]]` — reference coords (assembly/global).
- `elements: dict[int, Element]`, `Element(eid:int, etype:str, connectivity:list[int])`.
- `node_sets: dict[str, list[int]]`  (resolve `generate`, instance-qualified names OK).
- `element_sets: dict[str, list[int]]`.
- `materials: dict[str, Material]`,
  `Material(name, user_material:bool, constants:list[float], depvar:int|None)`.
- `sections: list[Section]`, `Section(elset:str, material:str, kind='solid')`.
- `element_material: dict[int, str]` — resolved element id → material name.
- `boundaries: list[Boundary]`,
  `Boundary(target:str|int, dof_start:int, dof_end:int, value:float, kind:str, amplitude:str|None)`
  where `kind ∈ {'value','XSYMM','YSYMM','ZSYMM','ENCASTRE',...}`. Encode symmetry
  BCs into the affected DOFs (XSYMM ⇒ dof 1, YSYMM ⇒ dof 2, ZSYMM ⇒ dof 3 fixed).
- `cloads: list[Cload]`, `Cload(target, dof, value, amplitude)`.
- `dsloads: list` (parse if present, may be minimally structured).
- `equations: list` (parse rows if present; not fully applied yet).
- `includes: list[str]` — files pulled in via `*Include, input=...` (parser must
  follow `*Include` relative to the parent file's directory).

Parser rules:
- Keyword lines start with `*`; comments start with `**`; keywords are
  case-insensitive; parameters are `key=value`, comma-separated.
- Support `*Node, *Element, *Nset, *Elset (incl. generate), *Solid Section,
  *Material, *User Material, *Depvar, *Boundary, *Cload, *Dsload, *Equation,
  *Include`. Ignore unknown keywords gracefully (record them in
  `model.unsupported_keywords: list[str]`), never crash.
- `*Part`/`*Instance`/`*Assembly`: for these single-instance, no-transform files,
  flatten to one global model. Record instance name; if a real transform appears,
  raise a clear NotImplementedError rather than silently ignoring it.

Self-test (`python abaqus_inp_parser.py`): parse `Compression111.inp` and assert
216 nodes, 125 C3D8 elements, material `CPURANIUM` user_material=True depvar=125,
11 constants.

## 4. Residual assembler — `c3d8_residual.py`

Pure numpy, **no dependency on the parser** (takes plain arrays), so it is
independently testable. Required public API:

```python
ABAQUS_C3D8_NODES      # (8,3) node natural coords, section 2
ABAQUS_C3D8_GAUSS      # (points:(8,3), weights:(8,)) Abaqus IP order, section 2

shape_functions(xi)        -> (8,)      # N_a(xi)
shape_grad_natural(xi)     -> (8,3)     # dN_a/dξ
b_matrix_reference(Xe, xi) -> (B0:(6,24), detJ0)   # small-strain B in ref config
b_matrix_spatial(xe, xi)   -> (B :(6,24), detJ )   # B in CURRENT config (xe deformed coords)

# Small-strain internal force (reference config), for patch tests:
element_internal_force_small_strain(Xe, sigma_ip, gauss=ABAQUS_C3D8_GAUSS) -> (24,)
#   Xe: (8,3) ref coords; sigma_ip: (8,6) Voigt Cauchy≈stress at the 8 IPs.

# Finite-strain internal force (current config) — the Abaqus-matching path:
element_internal_force_finite_strain(Xe, Ue, sigma_ip, gauss=ABAQUS_C3D8_GAUSS) -> (24,)
#   Ue: (8,3) or (24,) nodal displacements; xe = Xe + Ue; integrate B_spatial^T σ dv.

# Material tangent (used by tangent FD test); geometric term included for finite strain:
element_tangent(Xe, Ue, Dmat_ip, sigma_ip=None, mode='finite') -> (24,24)
#   Dmat_ip: (8,6,6) consistent material tangent per IP.

# Global scatter:
assemble_global_internal_force(node_ids, coords, connectivity, U, sigma_all,
                               mode='finite') -> (F_global:(ndof,), node_id_to_index)
#   node_ids: list[int]; coords: (nnode,3); connectivity: list[(eid, [8 node ids])];
#   U: (ndof,) or dict; sigma_all: dict[eid]->(8,6). mode in {'finite','small'}.
```

Verification (`python c3d8_residual.py`, all Abaqus-independent — these are the
Step-4/Step-7 gold checks that do NOT need Abaqus):

1. **Divergence-theorem patch test:** for a constant stress field σ0 over the real
   element geometry (pull a hex out of `Compression111.inp`), the assembled nodal
   forces `∫ B^T σ0 dV` must equal the consistent nodal forces of the surface
   tractions `t = σ0·n` (compute via `∫_∂Ω N^T t dA` on the 6 faces). Assert match
   to ~1e-8 (relative). Also assert Σ over nodes = 0 (self-equilibrium).
2. **Uniform uniaxial sanity:** σ = diag(0,0,S33) on a unit cube ⇒ +z face carries
   +S33·A split equally over its 4 nodes; −z face the opposite; side faces zero.
3. **Tangent finite-difference check** (writes results to
   `tests/c3d8_tangent_fd_check/`): pick one element; use a constant isotropic
   linear-elastic `D`; define `r(U) = element_internal_force_small_strain(Xe, σ(U))`
   with `σ(U) = D·(B0·U)` per IP; build `K = Σ B0^T D B0 detJ0 w`; FD-perturb each of
   the 24 DOFs (central difference, h≈1e-6·scale) and compare `dr/dU` to `K`.
   Report absolute error, relative error (Frobenius), and max entry error.

## 5. Field-exchange format (Abaqus ↔ Python), written by `extract_abaqus_fields.py`

Portable JSON (`fields.json`), Abaqus-Python-2.7 friendly (json only, numpy optional):

```json
{
  "odb": "Job-1.odb", "instance": "PART-1-1", "element_type": "C3D8",
  "nodes":    { "1": [x,y,z], ... },                         // undeformed ref coords
  "elements": { "1": [n1,...,n8], ... },                     // connectivity
  "frames": [
    { "frame": 5, "step": "Step-1", "time": 1.0,
      "U":  { "1": [ux,uy,uz], ... },
      "RF": { "1": [rx,ry,rz], ... },
      "S":  { "1": [[s11,s22,s33,s12,s13,s23]×8 IPs], ... }, // per element, 8 IPs, Abaqus IP order
      "SDV":{ "1": [[...]×8], ... }                          // optional, all SDVs per IP
    }
  ]
}
```

`stress_driven_residual.py` (integration driver, written separately) consumes this
+ the parsed model to assemble the residual and compare to `RF`.

## 6. Hard rules

- Never modify anything under `sources/`. Never run `git` commands. Do not commit.
- Never copy code from `sources/copyleft/` or `sources/license-unknown/` — those are
  reference-only. The permissive core builds only against `sources/permissive/`.
- Do not claim a comparison against Abaqus passed unless an Abaqus run actually
  produced the numbers. Abaqus is not installed in this environment; the Abaqus
  scripts are delivered ready-to-run and clearly marked as awaiting an Abaqus run.
- No OTI/HYPAD. No 150-parameter model. No final-step parameter overloading —
  plasticity is history-dependent; any replay must march through increments.
