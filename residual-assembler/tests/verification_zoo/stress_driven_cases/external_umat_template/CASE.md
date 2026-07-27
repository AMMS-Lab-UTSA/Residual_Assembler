# CASE — External UMAT via stress-driven (generic template)

Applies to any external UMAT candidate (viscoelastic, J2 plasticity, Mazars
damage, Neo-Hookean hyperelastic). The material model is irrelevant to the
residual once its stress field is exported.

- **Source**: any external UMAT in `sources/external_manifest.md` (respect its license)
- **License**: per source (permissive → may run; copyleft/unknown → export data only, never reuse code)
- **Physics**: any stress-producing constitutive law
- **Element type**: continuum solid (C3D8 today)
- **User subroutine type**: UMAT (not called in this mode)
- **Required residual backend**: `stress_driven_c3d8` (Mode 1)
- **Residual class**: F — stress-driven verification candidate
- **Required inputs**: mesh + exported integration-point stress `S` (Voigt) + DOF field
- **Expected outputs**: `f_int = sum_k B^T sigma_k dV`; free-DOF `R ~ 0`; reactions vs `RF`
- **Offline tests possible**: **yes** — once a stress export (JSON/CSV) exists, `resasm assemble MODEL --mode stress-driven --fields fields.json` runs with no Abaqus and no external code
- **Abaqus tests required**: only to *produce* the export (run the job, extract the ODB)
- **Current status**: planned per source (blocked only on the export)
- **Next missing item**: run the Abaqus job for the chosen UMAT and export the stress field with `scripts/extract_odb_fields.py`, then `scripts/compare_residuals.py`

## Procedure (with Abaqus)
1. `scripts/run_abaqus_validation.py --job <case>` → runs the `.inp` + UMAT
2. `scripts/extract_odb_fields.py --odb <case>.odb --out fields.json` → S at IPs
3. `resasm assemble <model> --mode stress-driven --fields fields.json`
4. `scripts/compare_residuals.py` → free-DOF residual ~ 0, reactions vs `RF`
