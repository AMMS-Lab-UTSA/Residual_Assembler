# Clean-clone reproduction of the published `main` branches (2026-09-18)

Both repositories were cloned from GitHub (`git clone --branch main`) into a
new directory, and everything below ran from those clones: the joint
clean-install gate, both offline test suites and the examples. Nothing was
copied from a development checkout.

| repository | commit tested |
| --- | --- |
| Residual_Assembler | `3504a0294107e17d32de451f8dd3b7a94abccd1c` |
| UMAT_source_transformation | `1352114164a936973441c30b8afe7520e582f106` |

Environment: Linux, Python 3.11.7 (a new venv), gfortran 9.4.0, Abaqus
2021.HF5 (ODB export only), the external OTILib build for the source suites.
Installed versions: numpy 2.4.6, scipy 1.17.1, streamlit 1.64.0, pytest 9.1.1.

## The clean-install gate

```sh
python Residual_Assembler/scripts/clean_install_gate.py \
  --umat-repo UMAT_source_transformation --branch main \
  --python python3.11 --odb <presentation ODB> --cantilever <cantilever runs> \
  --work <new directory>
```

All 22 commands exited 0 and the gate passed. Both trees were clean and each
commit was its branch's published head (`final_branch_clean_clone: true`).

| check | result |
| --- | --- |
| wheels built from the clones | `residual_assembler-0.1.0` sha256 `825ffeb8…a994`, `umat_oti-1.1.0` sha256 `215d32cd…1a14` |
| provider build with the installed `umat-oti-provider` | built and verified (J2) |
| `resasm request` on the genuine ODB, against the analytic uniaxial J2 derivatives | relative error E 1.7e-7, initial yield stress 8.7e-16, H 8.8e-9 |
| the same request with every read of a Fortran source denied | byte-identical outputs |
| both GUIs (installed Streamlit apps) | reached HTTP readiness and stopped cleanly |
| full-size J2 cantilever (1,536 C3D8, 40 increments): `resasm request` on its ODB | completed, 22.7 s |
| same cantilever re-equilibrated (`resasm history --reequilibrate`) | 33.0 s; homogeneity identity at every increment 1.4e-12 of the largest term (bound 1e-10), 34 plastic increments, 240 outputs |

## Test suites and examples (from the clones, with the gate's venv)

| step | result |
| --- | --- |
| UMAT_source_transformation `python -m pytest -q` | **3370 passed, 158 skipped, 0 failed** (skips name a missing prerequisite: Abaqus-only data, the corpus cache, optional corpora) |
| Residual_Assembler `pytest -q -m "not abaqus and not arc and not network"` | 494 passed, 20 skipped, 6 failed |
| the five examples of each repository (`scripts/audit_recovery_usage.py --phase examples`) | every command exited 0 |

The six Residual_Assembler failures were one cause, in the verification-record
scripts: they required `umat_oti` to be imported from the UMAT checkout's
`src/`, and the clean-clone venv held the identical copy installed from the
wheel instead. The check now accepts an installed package whose files are
byte-identical to the checkout's (`verification/common.py`); a different
version is still refused. The development-environment failure of
`tests/integration/test_connected_j2.py::test_reproducer_fresh_build_and_nonzero_failure`
does not occur here: it came from an editable install of an older checkout,
and it passes from a clean clone.
