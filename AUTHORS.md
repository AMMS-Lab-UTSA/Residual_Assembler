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

## Licence

Residual_Assembler is licensed **GPL-3.0-only**; see [`LICENSE`](LICENSE). This
matches its companion product UMAT-OTI and the OTILib algebra both build on.

The `sources/` licence tiers predate this decision and still apply, but their
purpose changed with it: they are a distribution and obligation boundary rather
than a firewall around a non-copyleft core. See
[`sources/LICENSING.md`](sources/LICENSING.md).

## Third-party material

External Abaqus sources are pinned git submodules organised by licence tier and
are never vendored. See [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md) and
[`sources/LICENSING.md`](sources/LICENSING.md). OTILib is GPL-3.0 and stays an
external dependency.
