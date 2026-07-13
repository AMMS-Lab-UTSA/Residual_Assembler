# UMAT material-point replay adapter

Call the Oxford / Grilli crystal-plasticity **UMAT** for a single integration
point **outside Abaqus**, marching a deformation-gradient history
**increment by increment**. This is the material-point half of Verification
Mode 2 (`tests/umat_replay_vs_abaqus/`).

> **History dependence is mandatory.** The UMAT state (`STATEV` + the `/UMPS/`
> common block: `kFp`, `kgausscoords`, `kcurlFp`, `kSigma0`, ...) evolves from
> one increment to the next. You must replay the **whole loading path**; you may
> never jump straight to the final increment. The driver therefore marches the
> entire history inside **one process** so both `STATEV` and the common block
> persist across increments.

The UMAT source lives in
`sources/permissive/ngrilli_Oxford_Crystal_Plasticity/` and is **MIT-licensed**
(Copyright (c) 2023 Nicolò Grilli; see that folder's `LICENSE`). Nothing under
`sources/` is copied or modified by this adapter — the build compiles the
original files **by path**. If you instead copy `umat.for`, `kmat.f`,
`mycommon.f`, `kMaterialParam.f`, etc. into a build directory, they remain under
the MIT license and must keep the copyright/attribution notice.

---

## Contents

| File | Role |
|------|------|
| `umat_driver.f90` | Standalone `program umat_driver`: reads a history, calls `UMAT`, writes `STRESS`/`STATEV`/`DDSDDE` per increment. |
| `aba_stubs/ABA_PARAM.INC` | The Abaqus double-precision include (not part of the MIT sources) so the UMAT compiles standalone. |
| `aba_stubs/CoheleBulkMap.f` | Empty shim so `UEL.for`'s `include 'CoheleBulkMap.f'` resolves (UEL/cohesive path is unused). |
| `aba_stubs/aba_stubs.f` | Link-time shims for Abaqus solver utilities: `MutexInit/Lock/Unlock`, `SMAIntArrayCreate/Access`, `XIT`. |
| `aba_stubs/lapack_stub.f` | Portable reference `DGETRF`/`DGETRI` (MKL replacement) used by `utils.f::lapinverse`. |
| `aba_stubs/umat_mock.f90` | A **mock** isotropic-elastic UMAT with the exact Abaqus signature, to exercise the plumbing without ifort/Abaqus. |
| `build.sh` / `build.bat` | gfortran build recipes (bash / cmd). |
| `umat_replay.py` | Python 3 orchestrator: builds the per-increment input, runs the driver, compares to ODB `S`/`SDV`. |

---

## The Abaqus UMAT argument list (matches `umat.for` exactly)

```fortran
SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
 1 RPL,DDSDDT,DRPLDE,DRPLDT,
 2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
 3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
 4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
```

Dimensions/kinds (double precision via `ABA_PARAM.INC`, `IMPLICIT REAL*8(A-H,O-Z)`):
`STRESS(NTENS)`, `STATEV(NSTATV)`, `DDSDDE(NTENS,NTENS)`, `DDSDDT(NTENS)`,
`DRPLDE(NTENS)`, `STRAN(NTENS)`, `DSTRAN(NTENS)`, `TIME(2)`, `PREDEF(1)`,
`DPRED(1)`, `PROPS(NPROPS)`, `COORDS(3)`, `DROT(3,3)`, `DFGRD0(3,3)`,
`DFGRD1(3,3)`, `CMNAME` = `CHARACTER*80`; the rest are scalars
(`NDI,NSHR,NTENS,NSTATV,NPROPS,NOEL,NPT,LAYER,KSPT,KSTEP,KINC` integer).

What this UMAT actually consumes (from `umat.for`): it forms `F = DFGRD0`,
`Fdot = (DFGRD1-DFGRD0)/DTIME`, `L = Fdot·F⁻¹`, then calls `kmat` for the stress
update. It **ignores** `STRAN`/`DSTRAN` (kinematics come from `DFGRD0/DFGRD1`),
so the driver leaves those zero. `PROPS` for this model:

* `PROPS(1)` = crystal type: **0=HCP**, 1/2=BCC, 3=Carbide, 4=Olivine, 5=α-U.
  (The `Compression111.inp` example sets `PROPS(1)=0` → HCP branch, `nSys=12`.)
* `PROPS(2:10)` = rotation matrix `R11,R12,R13,R21,R22,R23,R31,R32,R33`
  (crystal → sample).
* `PROPS(11)` = grain index.

**First increment must be `KSTEP=1, KINC=1`** — the UMAT's one-time init block
(`if (kinc<=1 .and. kstep==1)`) seeds `STATEV` (rotation into 1–9, `Fp=I`,
cohesive `kSigma0`, ...). Skipping it corrupts the whole replay.

---

## STATEV (SDV) layout, inferred from `kmat.f` / `umat.for`

The stress convention is pinned by CONTRACT §0: `kmat.f` stores the **Cauchy**
(true) stress in `usvars(48:53)` (Abaqus order `σ11,σ22,σ33,σ12,σ13,σ23`) and
`umat.for:388` copies it into `STRESS`.

| STATEV index | Meaning (source) |
|---|---|
| `1–9`  | rotation matrix `gmatinv` (crystal→sample); `usvars(j+(i-1)*3)=gmatinv(i,j)` |
| `10`   | equivalent plastic strain `p` |
| `11–16`| total plastic strain (Voigt) |
| `17–22`| total strain (Voigt) |
| `26`   | total GND density `gndtot` |
| `32`   | max plastic strain rate |
| `33`   | von Mises plastic strain rate `pdot` |
| `34`   | `xtau` |
| **`35`** | **cumulative plastic slip** `slip` |
| `36`   | solute/irradiation hardening `tauSolute` |
| `38–46`| curl of `Fp` (9) |
| **`48–53`** | **Cauchy stress** `σ11,σ22,σ33,σ12,σ13,σ23` → copied to `STRESS` |
| `54`   | sessile SSD density `rhossd` |
| `55`   | von Mises stress `vms` |
| `56`   | `max(tauc)` (then `56+i` reused for `gndcut`/`rhofor`) |
| `65`   | `rhosub` |
| `71`   | temperature |
| `72–80`| elastic strain in crystal frame `EECrys` (9) |
| **`81–89`** | **plastic deformation gradient** `Fp` (9); init to `I` at `KINC=1` |
| `90–101` | accumulated slip per system `usvars(89+i)` (phase-dependent) |
| `99–101` | critical resolved shear stress `tauc` (these are the `SDV99–101` the example requests) |
| `106+i`| twin volume fraction (twin runs only) |
| `124`  | cohesive effective opening `kDeltaEff` |
| `125`  | cohesive max stress `kSigma0` (or `max|tau|`, branch-dependent) |

Indices in the `90+` range are reused differently per crystal type / twin
branch; treat 1–9, 11–22, 35, 48–53, 81–89 as the stable, model-independent
ones.

---

## `mycommon.f` and the `/UMPS/` common block

`umat.for` does `include 'mycommon.f'`, which declares the `/UMPS/` **static**
common block sized by compile-time parameters:

```
nElements = 18315      ! must be >= the real element count in your mesh
nintpts   = 8          ! C3D8 integration points
```

* In this source the arrays are **static COMMON** (not allocatable), so
  `UEXTERNALDB` only initialises mutexes — the driver does **not** need to
  allocate anything; zero-initialised BSS is correct at start.
* For `Compression111` (125 elements) the shipped `nElements=18315` is a valid
  upper bound (it just over-allocates ~30 MB). For a lean single-IP replay you
  may copy `mycommon.f` into a build dir and set `nElements = 1` (keep the MIT
  notice). Element/IP indices used in a single-point replay are then `NOEL=1`,
  `NPT=1`.
* The correct `kMaterialParam.f` / `mycommon.f` for your target mesh must be on
  the include path. `build.sh` points `-I` at the sources folder, which uses the
  repo-root copies.

---

## Driver I/O contract

All plain text, **list-directed** (whitespace/newline separated) so Fortran and
Python agree trivially. `umat_replay.py` writes the input and reads the output.

### Input  (`umat_driver <in> <out>`)

```
NSTATV NPROPS NTENS NINC
PROPS(1) ... PROPS(NPROPS)
STATEV0(1) ... STATEV0(NSTATV)        # initial state; usually all 0
# then NINC increment records, each:
NOEL NPT KSTEP KINC
DTIME TEMP DTEMP TIME1 TIME2
F0_11 F0_12 F0_13 F0_21 F0_22 F0_23 F0_31 F0_32 F0_33     # DFGRD0, row-major
F1_11 F1_12 F1_13 F1_21 F1_22 F1_23 F1_31 F1_32 F1_33     # DFGRD1, row-major
```

### Output (one record per increment)

```
# two comment lines
NSTATV NTENS NINC
KINC
STRESS(1..NTENS)
STATEV(1..NSTATV)
DDSDDE(1..NTENS*NTENS)     # row-major
PNEWDT
```

---

## Build status (gfortran 15.2, this environment) — honest

`build.sh` / `build.bat` do two things:

1. **Mock path — builds and runs.** `umat_driver.f90` + `aba_stubs/umat_mock.f90`
   link into `build/umat_driver_mock`. Verified: a 2-increment uniaxial history
   produces `σ33 = 134.6 → 269.2`, `STATEV(35)` accumulates `9.99e-6 → 3.0e-5`
   (proving the common block / state persists), `STATEV(1:9)` holds the rotation
   matrix, `STATEV(48:53)` mirrors `STRESS`. Identity deformation → zero stress.

2. **Real UMAT path — does NOT compile with gfortran (expected).** gfortran
   preprocesses `#include <SMAAspUserSubroutines.hdr>`, then compiles cleanly
   through `mycommon.f`, `kmat.f`, and the whole slip/hardening include stack,
   stopping only on **two ifort-specific constructs in the untouchable sources**:

   * `Error: CRAY POINTEE attribute conflicts with TARGET attribute`
     — `umat.for:150–154` declares the twin arrays both `integer, target ::`
     **and** `pointer(ptr…, …)` (Cray pointer). ifort allows a symbol to be
     both; gfortran forbids it. (Twin machinery only; unused at `twinon=0`, but
     it is a declaration so it blocks the whole translation unit.)
   * `Error: Return type mismatch of function 'trace' (REAL(4)/REAL(8))`
     — `utils.f::trace` is `real*8`, but `kmat.f` (which intentionally comments
     out `ABA_PARAM.INC`) uses it under default `real(4)` implicit typing. ifort
     tolerates this; gfortran errors.

   `-std=legacy -fallow-argument-mismatch` downgrades the other differences
   (`TEMP/DTEMP/PNEWDT REAL8→REAL4`, `cubicslip REAL8→INTEGER4`, `gndmob` rank)
   to warnings, but the two above remain hard errors and cannot be fixed without
   editing `sources/` (forbidden).

   **Conclusion:** the real crystal-plasticity object must be built with the
   toolchain the model targets — **Intel `ifort` + Abaqus (MKL)**, per
   `sources/.../abaqus_v6.env` (`/Qmkl /extend-source /fpp /iface:cref
   /recursive`). Typical route:

   ```
   # inside the Abaqus environment (provides ABA_PARAM.INC, MKL, SMA*, Mutex*):
   abaqus make library=umat.for              # or ifort -c with the env flags
   # then link our driver against the resulting object:
   ifort umat_driver.o umat.o -o umat_driver   # MKL supplies DGETRF/DGETRI
   ```

   When linking the real object **outside** Abaqus (e.g. a patched gfortran or a
   bare ifort build), also link `aba_stubs/aba_stubs.f` (Mutex/SMA/XIT) and
   `aba_stubs/lapack_stub.f` (if no MKL/LAPACK). **Do not** link `lapack_stub.f`
   when MKL/LAPACK is present (duplicate `DGETRF`/`DGETRI`).

### What `lapinverse` / `kdeter` depend on
`lapinverse` and `KDETER` are **both defined in the repo's own `utils.f`** — not
Abaqus-provided. `KDETER` is self-contained. `lapinverse`'s only external
dependency is LAPACK `DGETRF`/`DGETRI` (MKL under Abaqus). `aba_stubs/lapack_stub.f`
supplies correct portable versions for a non-MKL build.

---

## Replay procedure

```bash
# 1. build
./build.sh                       # or build.bat on Windows cmd

# 2a. plumbing / known-answer check (mock material, no Abaqus needed)
python umat_replay.py --dry-run --mock

# 2b. replay one IP from an extracted ODB (Verification Mode 2)
#     (needs the REAL umat_driver; otherwise pass --mock to test the wiring)
python umat_replay.py --fields fields.json --element 63 --ip 1
```

`umat_replay.py`:
* reconstructs `DFGRD1` at the chosen element/IP from the ODB nodal `U` and the
  C3D8 shape-function gradients (CONTRACT §2/§4; imports
  `residual_core/c3d8_residual.py` if available, else uses an inline copy),
* marches frame-to-frame (`DFGRD0` = previous frame's F, `DFGRD1` = this frame's,
  `DTIME` = Δt), starting from the undeformed reference `F = I` at `t = 0`,
* runs the driver and prints per-increment `STRESS` vs ODB `S` and `STATEV` vs
  ODB `SDV` (abs/rel error).

> **Accuracy caveat.** Exact replay needs `DFGRD1` at **every** solver
> increment. ODB field output written every N increments (the example uses
> `frequency=200`) is too coarse for an exact plasticity replay — frame-to-frame
> marching is then a sub-stepping approximation. For an exact Mode-2 check,
> request field output **every increment** (`*Output, field, frequency=1`) or
> use fixed time increments so each frame equals one increment.
```
