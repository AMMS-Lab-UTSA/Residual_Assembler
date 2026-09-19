# Clean-clone reproduction of the published `main` branches (2026-09-19)

Both repositories were cloned from GitHub into a new directory, and everything
below ran from those clones by the documented steps only. Nothing was copied
from a development checkout. The run is one command,
[scripts/reproduce_from_clean_clones.sh](../../scripts/reproduce_from_clean_clones.sh)
([docs/INSTALL.md](../INSTALL.md), section 7):

```sh
export PYOTI_PATH=<OTILib build> OTILIB_ROOT=<OTILib build>
bash scripts/reproduce_from_clean_clones.sh --python python3.11 \
    --odb <Analysis.odb of examples/presentation_request> \
    --cantilever <folder with j2/cantilever_j2_nominal.inp and .odb> NEW_DIR
```

| repository | commit tested |
| --- | --- |
| Residual_Assembler | `bf3600cb58d8f267020573476dc26c826c5048f5` |
| UMAT_source_transformation | `5dcdd8d6041f970c315b999e3e0ace4fbf376bd3` |

Environment:

- Linux; Python 3.11.7 in a new virtual environment made by the gate; gfortran
  9.4.0.
- Abaqus 2021.HF5, used only to read the two ODBs; no analysis was started.
- OTILib built as INSTALL.md section 4 documents, at the pinned commit
  `a4b7a05ca275e8d441b0b717b7b272e96728ffcc`.
- Installed versions: numpy 2.4.6, scipy 1.17.1, streamlit 1.64.0, pytest
  9.1.1.

Total time: 26 min 40 s. Every step exited 0, and both clones were left
unmodified (empty `git status` after the run).

| step | exit | time | result |
| --- | ---: | ---: | --- |
| clone both repositories | 0 | 5 s | both on `main`, each commit the published head |
| clean-install gate | 0 | 109 s | all 22 commands exited 0; the gate passed (`final_branch_clean_clone: true`) |
| Residual_Assembler offline suite | 0 | 462 s | **647 passed, 20 skipped, 0 failed** |
| UMAT_source_transformation suite | 0 | 820 s | **3399 passed, 160 skipped, 0 failed** |
| worked examples (`scripts/audit_recovery_usage.py --phase examples`) | 0 | 204 s | 31 of 31 commands exited 0 and passed their checks (16 Residual_Assembler, 15 UMAT) |

The suites ran with the documented commands, in the environment the gate made:
Residual_Assembler `pytest -q -m "not abaqus and not arc and not network"` and
UMAT_source_transformation `python -m pytest -q`, each with `--junitxml` added.

Every skip names a missing prerequisite. In Residual_Assembler that is the
permissive source submodule, the frozen corpus store or Playwright. In
UMAT_source_transformation it is Abaqus-only data, the corpus cache or optional
corpora.

## What the gate checked

| check | result |
| --- | --- |
| wheels built from the clones | `residual_assembler-0.1.0` sha256 `1a1c7d49…0b05`, `umat_oti-1.1.0` sha256 `47126212…6135` |
| the installed packages | every module imported from the new environment's site-packages, nothing editable, the packaged data files present |
| provider build with the installed `umat-oti-provider` | built and verified (J2) |
| `resasm request` on the genuine ODB, against the analytic uniaxial J2 derivatives | relative error: E 1.7e-7, initial yield stress 8.7e-16, H 8.8e-9 |
| the same request with every read of a Fortran source denied | byte-identical outputs |
| both GUIs (installed Streamlit apps) | rendered headlessly, reached HTTP readiness and stopped cleanly |
| full-size J2 cantilever (1,536 C3D8, 40 increments): `resasm request` on its ODB | completed in 21.7 s |
| the same cantilever re-equilibrated (`resasm history --reequilibrate`) | 30.8 s; homogeneity identity at every increment 1.4e-12 of the largest term (bound 1e-10); 34 plastic increments, 240 outputs |

## The completion ledger

`tools/update_completion_ledger.py` set the status of every row of
[docs/COMPLETION_LEDGER.md](../COMPLETION_LEDGER.md) from this run alone.

- **265 of 274 rows** are `PASS: reproduced from a clean installation`: every
  test, gate check and example step mapped to them succeeded here.
- **8 rows are `BLOCKED`**, each naming the external resource this run does not
  exercise: an Abaqus job that compiles the transformed UMAT, the third-party
  corpus cache, the network-fetched source submodule, or a browser.
- **1 row (finite-strain material state) is `NOT STARTED`.**

## Earlier run

The first clean-clone run of the final branches, on 2026-09-18 (Residual_Assembler
`3504a02`, UMAT_source_transformation `1352114`), passed the gate and the UMAT
suite. It found six Residual_Assembler failures, all with one cause: the
verification-record scripts required `umat_oti` to be imported from the UMAT
checkout's `src/`, not from an identical installed copy. That was corrected
(`verification/common.py`), and this run has no failures.
