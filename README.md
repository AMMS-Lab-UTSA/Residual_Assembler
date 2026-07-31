# HYPAD material-sensitivity framework

**Exact material Jacobians and parameter sensitivities for Abaqus UMATs — automatically, from the constitutive code you already have, with the source kept private.**

Take a normal Abaqus **UMAT** and get back, *without hand-coding a tangent and without finite differences*:

- the consistent material Jacobian **DDSDDE** = ∂σ/∂ε,
- the parameter sensitivities **DSIGMA_DP** = ∂σ/∂p and **DSTATEV_DP** = ∂(state)/∂p,

all to **machine precision**, and then assemble them into **full-field structural sensitivities** ∂(response)/∂p by replaying a single saved analysis — never re-solving.

```
  real UMAT  ─►  umat-oti  ─►  OTI object (.obj) + contract (.json)  ─►  residual-assembler  ─►  d(response)/d(parameter)
 (STRESS,STATEV,DDSDDE)        (… + DSIGMA_DP, DSTATEV_DP)                (one replay, no re-solve)
```

It works by **HYPAD / OTI** (hypercomplex, order-truncated-imaginary automatic differentiation): the constitutive routine is evaluated once over an enriched number system whose imaginary parts carry the exact derivatives. No perturbation step to tune, no approximation.

---

## Contents

- [Why](#why)
- [The two programs](#the-two-programs)
- [Requirements](#requirements)
- [Install](#install)
- [5-minute smoke test](#5-minute-smoke-test) — confirm it works on your machine
- [Quick start 1 — transform a UMAT (Program 1)](#quick-start-1--transform-a-umat-program-1)
- [Quick start 2 — assemble sensitivities (Program 2)](#quick-start-2--assemble-sensitivities-program-2)
- [Bring your own UMAT: the contract](#bring-your-own-umat-the-contract)
- [What you get back](#what-you-get-back)
- [A worked example (plate with a hole, J2 plasticity)](#a-worked-example-plate-with-a-hole-j2-plasticity)
- [Repository layout](#repository-layout)
- [How it is tested](#how-it-is-tested)
- [Platform notes (Windows & Linux)](#platform-notes-windows--linux)
- [Troubleshooting](#troubleshooting)
- [Deeper documentation](#deeper-documentation)
- [License](#license)

---

## Why

Calibration, optimization, and uncertainty quantification all need **derivatives of the FE response with respect to material parameters**. The usual options are bad:

- **Finite differences** re-run the whole analysis once per parameter (2N+1 runs for N parameters), are only approximate, and need a step size you have to babysit.
- **Hand-coded analytical tangents/sensitivities** are error-prone and have to be re-derived for every new constitutive model.

This framework gives you the **exact** derivatives **automatically** from the stress-update code, and — crucially for industry — **the material source never has to leave the developer's machine**. The material developer ships a compiled object + a small JSON contract; the collaborator computes sensitivities against it without ever seeing the equations.

---

## The two programs

Each subfolder is an independent, separately-distributable repository with its own `README.md`, `INTERFACE.md`, and tests.

| repo | who runs it | what it does |
|---|---|---|
| [**`umat-oti/`**](umat-oti) | material developer | Transforms a normal UMAT into an **OTI-enabled** object that also returns `DSIGMA_DP` / `DSTATEV_DP`, by hypercomplex AD. Emits the object **plus** a completed JSON interface contract. Original source stays private. |
| [**`residual-assembler/`**](residual-assembler) | collaborator / analyst | Consumes that object + a saved FE record and computes `d(response)/d(parameter)` by the **residual method**: assemble `R = ∫Bᵀσ dV − f_ext`, then `∂R/∂p = ∫Bᵀ(∂σ/∂p) dV`, and solve `K (∂u/∂p) = −∂R/∂p` **once** on the already-factored stiffness — no re-analysis. |

The hand-off is a compiled object + JSON contract over a shared C-ABI; **no OTI object crosses back**. It is documented identically in each repo's `INTERFACE.md`.

---

## Requirements

| | Program 1 (`umat-oti`) | Program 2 (`residual-assembler`) |
|---|---|---|
| Python | 3.10+ | 3.10+ |
| Python deps | `numpy>=1.26`, `sympy>=1.12` | `numpy` (`scipy` for the FE demo) |
| Compiler | **`gfortran` on `PATH`** (or `FC`) — the transform emits and links Fortran | none for assembly |
| Abaqus | only for optional end-to-end verification | only for `verify-job` against a real `.odb` |
| OS | **Windows and Linux**, same repo & CLI | Windows and Linux |

`gfortran` check: `gfortran --version` should print a version. On Windows, MSYS2/MinGW-w64 or the equivalent works; the gfortran runtime directory is discovered automatically at load time (see [Platform notes](#platform-notes-windows--linux)).

---

## Install

Each program installs independently. From the repo root:

```bash
# Program 1 — the transformer
cd umat-oti
pip install -e .

# Program 2 — the residual assembler (exposes the `resasm` CLI)
cd ../residual-assembler
pip install -e .
```

---

## 5-minute smoke test

Confirm the whole toolchain works on your box before doing anything with your own material. This builds a bundled material end-to-end and checks the derivatives against a converged finite-difference reference of the original UMAT.

```bash
cd umat-oti
python oti_provider/umat_transform.py build oti_provider/materials/m3_j2/contract.json
```

**What success looks like:** the build compiles the OTI object, runs the original vs. transformed **parity checks**, and prints per-parameter agreement of `DSIGMA_DP` against the finite-difference reference — you should see relative errors around **1e-9 or smaller**. It writes `umat_m3_j2_oti.obj` and `umat_m3_j2_oti.json` next to the material.

Try a few more (they exercise different physics):

```bash
python oti_provider/umat_transform.py build oti_provider/materials/m1_elastic/contract.json   # linear elastic
python oti_provider/umat_transform.py build oti_provider/materials/m6_fcc/contract.json        # FCC crystal plasticity (10 params)
```

Then confirm the assembler imports and runs:

```bash
cd ../residual-assembler
python -c "import residual_core; print('residual-assembler OK')"
resasm --help
```

If both of those work, you're ready to test with your own model.

---

## Quick start 1 — transform a UMAT (Program 1)

A material lives in `umat-oti/oti_provider/materials/<name>/` as a `umat.for` (your stress update) + a `contract.json` (what to differentiate). Build:

```bash
cd umat-oti
python oti_provider/umat_transform.py build oti_provider/materials/<name>/contract.json
```

You get two artifacts next to the material:

- **`umat_<name>_oti.obj`** — relocatable object exposing three Fortran symbols: `UMAT` (original, unchanged), `UMAT_OTI_INTERNAL` (the OTI-lifted UMAT), `UMAT_OTI_EVAL` (offline entry returning the derivative arrays), plus the bundled OTI runtime.
- **`umat_<name>_oti.json`** — the completed interface contract (symbols, array layouts, dimensions, binary metadata, verification evidence).

For **your own** UMAT source (not staged as a material folder), drive the transform straight from a compact JSON:

```bash
python transform_from_json.py path/to/your_contract.json --out path/to/output_dir
```

> The transform never runs Abaqus and never needs the ODB. Optional end-to-end Abaqus verification lives in `verify_abaqus_local.py` / `validate_all_local.py` and is only used on a machine with Abaqus configured.

---

## Quick start 2 — assemble sensitivities (Program 2)

`residual-assembler` takes a **converged FE record** (mesh + material + solution field + loads/BCs) and the OTI object, and returns `d(response)/d(parameter)`. You never hand it the global residual `R` — it assembles `R` for you. Your model stays on your machine.

**Verify the assembler reproduces a solved job** (one file, one command):

```jsonc
// job.json
{ "name": "my_check", "model": "my_model.inp", "odb": "my_model.odb",
  "mode": "stress-driven", "compare_reactions": true, "tol": 1e-6 }
```

```bash
resasm verify-job job.json
```

It exports the ODB fields for you, assembles the residual, and checks the free-DOF residual is ~0 and the assembled reaction equals Abaqus `RF`. A runnable elastic C3D8 example that passes to ~1e-16 ships in [`residual-assembler/examples/abaqus_elastic_c3d8/`](residual-assembler/examples/abaqus_elastic_c3d8/). No Abaqus on the box? Point `"fields": "fields.json"` at an export made elsewhere.

**The assembly path** (each step tells you exactly what it still needs):

```bash
resasm inspect model.inp                            # what do I have? what's missing?
resasm requirements model.inp --mode stress-driven  # the single next thing to provide
resasm init-assembly --model model.inp --solution U.npy
resasm check resasm.yml
resasm run   resasm.yml
```

---

## Bring your own UMAT: the contract

The contract is deliberately small — you declare only genuine decisions; everything else (stress-update line, DDSDDE block, promoted variables, OTI directions) is inferred. This is the bundled J2 example, [`oti_provider/materials/m3_j2/contract.json`](umat-oti/oti_provider/materials/m3_j2/contract.json):

```jsonc
{
  "material":   "j2_linear_hardening",   // a name
  "source":     "umat.for",              // your UMAT, next to this file
  "kinematics": "small_strain",
  "dimensions": {
    "stress_components":   6,            // NTENS
    "material_properties": 4,            // NPROPS
    "state_variables":     1             // NSTATV
  },
  "parameters": {                        // which PROPS to differentiate, by name -> PROPS index
    "E": 1, "nu": 2, "SIGY0": 3, "H": 4  // parameter order fixes the derivative-column order
  },
  "derivatives": { "stress": true, "statev": "auto" }
}
```

- **`parameters`** is the only thing that fixes what sensitivities you get, and in what column order. Drop the block entirely and you get a Jacobian-only build (DDSDDE, no `DSIGMA_DP`).
- **`derivatives.statev: "auto"`** lets the tool decide which state variables need derivative tracking for path-dependent materials.

---

## What you get back

One enriched evaluation of the OTI UMAT returns, **synchronized through the same nonlinear execution** (original equations untouched):

| output | meaning | shape |
|---|---|---|
| `STRESS` | stress | `(NTENS)` |
| `STATEV` | physical state (history only) | `(NSTATV)` |
| `DDSDDE` | consistent tangent ∂σ/∂ε | `(NTENS, NTENS)` |
| `DSIGMA_DP` | stress–parameter sensitivity ∂σ/∂p | `(NTENS, NPARAM)` |
| `DSTATEV_DP` | state–parameter sensitivity ∂s/∂p | `(NSTATV, NPARAM)` |

`DSIGMA_DP` / `DSTATEV_DP` are **dedicated output arguments** — they are never packed into `STATEV`; `STATEV` stays the physical material state only.

---

## A worked example (plate with a hole, J2 plasticity)

A quarter-symmetry plate with a central hole is pulled into plasticity (E, ν, σ_y0, H). Running the material point through the analysis and comparing HYPAD against finite differences:

- **DDSDDE** (all independent tangent components) agrees with FD to **≈ 10⁻¹¹** across the whole elastic→plastic transition.
- **DSIGMA_DP** (∂σ/∂p for all four parameters, per increment) agrees with FD to **10⁻¹² – 10⁻¹⁰**.
- The residual method then produces **full-field** ∂σ_vM/∂p across the mesh from **one** replay — showing, for free, *which material property controls the stress at every point* (elastic regions governed by E, the plastic band by σ_y0), and even the **stress-uncertainty field** from parameter scatter, with no Monte-Carlo.

The physics comes out exactly right: E and ν sensitivities collapse after yield, σ_y0 switches on at yield, H grows with accumulated plastic strain.

---

## Repository layout

```
.
├── README.md                    <- you are here
├── umat-oti/                    <- Program 1: UMAT -> OTI transformer
│   ├── oti_provider/
│   │   ├── umat_transform.py    <- build CLI:  python oti_provider/umat_transform.py build <contract>
│   │   └── materials/           <- bundled examples (m1_elastic, m3_j2, m6_fcc, ...)
│   ├── transform_from_json.py   <- transform your own UMAT from a compact JSON
│   ├── src/umat_oti/            <- the transformer package
│   ├── UMATs/OTI/               <- bundled OTILib / MultiZ hypercomplex runtime
│   ├── examples/  templates/  tests/
│   ├── README.md · START_HERE.md · INTERFACE.md
│   └── pyproject.toml
└── residual-assembler/          <- Program 2: residual sensitivity solver
    ├── residual_core/           <- assembly + replay + CLI (`resasm`)
    ├── examples/abaqus_elastic_c3d8/   <- runnable verification example
    ├── schemas/  tests/  validation/
    ├── README.md · QUICKSTART_USER.md · METHODOLOGY.md · INTERFACE.md
    └── pyproject.toml
```

---

## How it is tested

- **Program 1** validates every transform against a **converged central-finite-difference reference of the original UMAT** (adaptive step ladder, plateau-selected, Richardson-checked) and reports per-parameter evidence. `validate_all_local.py` runs the whole material suite; `tests/` holds unit tests (`pip install -e ".[test]" && pytest`).
- **Program 2** ships runnable verification jobs (`resasm verify-job`) and a self-contained elastic C3D8 that reproduces the Abaqus reaction to ~1e-16. `pytest` covers assembly and replay.
- The framework has been exercised across a spread of physics families (linear/anisotropic elasticity, J2 plasticity with several hardening laws, viscoplasticity, crystal plasticity, generalized continua, damage, hyperelasticity).

A tester's fastest confidence check is the [5-minute smoke test](#5-minute-smoke-test): if the bundled builds print ~1e-9 parity against finite differences, the toolchain is healthy on that machine.

---

## Platform notes (Windows & Linux)

Program 1 runs natively on **both** from the same repo and CLI. The build compiles two throwaway libraries and loads them through `ctypes`:

- the library gets the platform suffix automatically (`.dll` / `.so`);
- on **Windows**, the active gfortran runtime directory is discovered dynamically and registered with `os.add_dll_directory` (Python 3.8+ ignores `PATH` for *dependent* DLLs, so gfortran on `PATH` is necessary but not sufficient) — **no compiler path is hard-coded**;
- if a load ever fails you get a diagnostic naming the platform, compiler, library, and the missing dependency;
- `UMAT_OTI_STATIC_FORTRAN=1` links the Fortran runtime statically as a fallback.

Both platforms are exercised in CI.

---

## Troubleshooting

| symptom | fix |
|---|---|
| `gfortran: command not found` during build | install gfortran and put it on `PATH` (or set `FC=/path/to/gfortran`). |
| Windows: DLL load fails after a successful compile | a gfortran runtime isn't discoverable — ensure the gfortran that compiled is on `PATH`; or set `UMAT_OTI_STATIC_FORTRAN=1` to link statically. |
| `ModuleNotFoundError: sympy` | `pip install -e .` inside `umat-oti` (sympy powers the OTI module generator). |
| parity check prints large errors | the contract's `parameters` indices or `dimensions` don't match the UMAT's `PROPS`/`NTENS`/`NSTATV` — recheck them. |
| `resasm: command not found` | `pip install -e .` inside `residual-assembler` (it registers the `resasm` entry point). |

---

## Deeper documentation

- **Program 1:** [`umat-oti/README.md`](umat-oti/README.md), [`START_HERE.md`](umat-oti/START_HERE.md), [`INTERFACE.md`](umat-oti/INTERFACE.md)
- **Program 2:** [`residual-assembler/README.md`](residual-assembler/README.md), [`QUICKSTART_USER.md`](residual-assembler/QUICKSTART_USER.md), [`METHODOLOGY.md`](residual-assembler/METHODOLOGY.md), [`INTERFACE.md`](residual-assembler/INTERFACE.md)
- The **interface contract** (the hand-off) is specified identically in both `INTERFACE.md` files.

---

## License

Both programs are **MIT-licensed** (matching OTILib / pyoti). See each subfolder's `LICENSE`.
