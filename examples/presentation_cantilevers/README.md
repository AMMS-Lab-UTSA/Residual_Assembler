# The presentation cantilevers (IMQCAM Annual Meeting, slides 15, 33, 39)

Two full-size models the presentation shows, rebuilt from the numbers the
slides print (the original decks and ODBs were not preserved):

| model | mesh | DOF | integration points | steps | parameters |
| --- | --- | ---: | ---: | ---: | --- |
| `j2` (slide 39) | 48 x 16 x 2 C3D8 | 7,497 | 12,288 | 40 | E, nu, SIGY0, H = 200000, 0.3, 250, 2000 |
| `fcc` (slide 15) | 24 x 8 x 2 C3D8 | 2,025 | 3,072 | 25 | C11, C12, C44, g0, gsat, h0, a, q, gd0, m (slide 15 values) |

Both are clamped at x = 0 with every tip node pushed in -y (0.7 mm and
0.25 mm), small strain, one increment per load step, every increment written
to the ODB. Units mm, N, MPa. The UMATs are the ORIGINAL sources in
UMAT_source_transformation: `parameter_sensitivity/models/m3_j2/umat.for` and
`parameter_sensitivity/models/m6_fcc/umat.for`.

Nothing here is large: the decks are generated, and the Abaqus output (about
650 MB for the J2 ODB and its reruns) belongs in a work directory outside the
repository.

## Files

| file | what it does |
| --- | --- |
| `gen_cantilever.py` | writes a deck: `python gen_cantilever.py j2 --out j2/claude_j2_nominal.inp` (`--props` overrides PROPS) |
| `export_odb.py` | every frame of an ODB to `.npz` (run with `abaqus python`) |
| `run_fd.sh` | perturbed-parameter Abaqus reruns at h = 0.02, 0.01, 0.005, one job at a time |
| `fd_reference.py` | central differences from those reruns, with the step-size plateau and the single-precision floor of ODB field output |
| `j2_request.json`, `fcc_request.json` | the collaborator's requests (tip reaction, displacement, stress, von Mises mean / root max, plastic strain, full fields, weighted shares) |

## Run it

With Abaqus (tested: 2021.HF5) and both packages installed, in a work
directory outside the repositories (`$EX` is this folder, `$UMAT` the UMAT
repository):

```sh
mkdir -p j2 && python $EX/gen_cantilever.py j2 --out j2/claude_j2_nominal.inp
cd j2 && abaqus job=claude_j2_nominal input=claude_j2_nominal.inp \
    user=$UMAT/parameter_sensitivity/models/m3_j2/umat.for double=both interactive && cd ..
grep "COMPLETED SUCCESSFULLY" j2/claude_j2_nominal.sta     # the verdict; see below

umat-oti-provider build $UMAT/parameter_sensitivity/models/m3_j2/contract_v2.json --out provider_j2
resasm request --model j2/claude_j2_nominal.inp --odb j2/claude_j2_nominal.odb \
    --material provider_j2/umat_m3_j2_oti.obj --request $EX/j2_request.json --out j2_results
```

`resasm request` sees a model outside the bounded presentation scope and hands
it to the history replay engine (`resasm history`), which says so. Add
`--reequilibrate` to `resasm history` to Newton-polish every recorded increment
to double-precision equilibrium before the sensitivities are taken.

Abaqus 2021.HF5 on a machine whose process IDs exceed 999,999 aborts with
signal 6 during teardown, after writing a complete ODB (a fixed-size buffer in
its bundled Intel runtime; see UMAT_source_transformation `docs/GUI.md`). Judge
a job by "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" in the `.sta` file, not by
the exit code.

Independent reference from Abaqus itself:

```sh
UMAT=$UMAT/parameter_sensitivity/models/m3_j2/umat.for bash $EX/run_fd.sh j2
abaqus python $EX/export_odb.py -- j2/claude_j2_nominal.odb j2/claude_j2_nominal_fields.npz
python $EX/fd_reference.py j2 .
```

The full comparison (both cantilevers, recorded and re-equilibrated, Abaqus
FD, figures and the slide-33 shares) is
`scripts/replay_history_cantilevers.py`; its results and the measured numbers
are in [docs/evidence/claude_C.md](../../docs/evidence/claude_C.md).

## Checks that need no finite differences

Both models are homogeneous of degree one in their stress-dimensioned
parameters: J2 with linear hardening in (E, SIGY0, H) at fixed nu, the FCC
crystal in (C11, C12, C44, g0, gsat, h0), because its slip rate depends only
on the ratio of resolved shear stress to slip resistance. Under prescribed
displacements, scaling those parameters together leaves every displacement
and plastic strain unchanged and scales every stress and reaction. So at every
increment, elastic or plastic,

    sum_p p dQ/dp = Q   for reactions, stresses and von Mises
    sum_p p dQ/dp = 0   for displacements and plastic strain

Measured on 2026-09-18 with `--reequilibrate`: largest residual 1.0e-12 of the
largest term for J2 and 1.2e-13 for FCC. On the recorded (single-precision)
state the residual is the size of the recorded equilibrium residual: 8.4e-5
(J2) and 5.4e-4 (FCC). `tests/replay_history/test_history_example.py` checks
the same identity on the committed small beam, and the clean-install gate
(`scripts/clean_install_gate.py --cantilever`) checks it on the J2 cantilever
from installed wheels.
