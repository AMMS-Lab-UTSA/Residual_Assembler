# Authors

Residual_Assembler is developed at the University of Texas at San Antonio,
alongside [UMAT-OTI](https://github.com/AMMS-Lab-UTSA/UMAT_source_transformation).
The two are separate products connected by a versioned contract.

| Author | Affiliation |
|---|---|
| Santiago García Botero | The University of Texas at San Antonio |
| Harry Millwater | UTSA, Department of Mechanical Engineering |
| Arturo Montoya | UTSA, Department of Civil Engineering |
| David Restrepo | UTSA, Department of Mechanical Engineering |

## Licence status

**This repository does not yet declare a licence for its own code.** Without one,
default copyright applies and no redistribution rights are granted, which is at
odds with the design intent recorded in
[`sources/LICENSING.md`](sources/LICENSING.md) of building "a clean permissive
core". Choosing the licence is an authors' decision and is deliberately not made
here; `tools/audit_repository_standards.py` reports the missing `LICENSE` file
until it is.

## Third-party material

External Abaqus sources are pinned git submodules organised by licence tier and
are never vendored. See [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md) and
[`sources/LICENSING.md`](sources/LICENSING.md). OTILib is GPL-3.0 and stays an
external dependency.
