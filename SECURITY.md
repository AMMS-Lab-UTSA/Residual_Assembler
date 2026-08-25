# Security policy

## Scope

Residual_Assembler parses Abaqus input decks and exported solution data, and
invokes a Fortran compiler on material subroutines you supply. The
security-relevant consequences are that it **parses untrusted input** and
**compiles and runs code you point it at**.

Treat a `.inp` file or a UMAT source the way you would treat any code you are
about to run. The framework does not sandbox them.

Your model never leaves your machine: the assembler makes no network calls
during assembly or sensitivity computation. The only network access in this
repository is the explicit submodule bootstrap in `scripts/`.

## Reporting a vulnerability

Use GitHub's private advisory workflow:

- <https://github.com/AMMS-Lab-UTSA/Residual_Assembler/security/advisories/new>

Please do not open a public issue for a vulnerability before it is addressed.
Include the input that triggers it and the command you ran. Expect an
acknowledgement within ten working days; this is academic research software
maintained alongside other duties.

## Out of scope

- Wrong residuals or sensitivities are correctness bugs. Open a normal issue with
  the recipe and solution data that reproduce them.
- Crashes in Abaqus or `gfortran` belong upstream.
- Anything requiring an attacker to already run code as you.

## Licence-tier isolation

Fetching a copyleft or licence-unknown submodule is possible but never automatic
(`update = none`). If you believe ordinary setup fetched a restricted source,
report it — that is a policy failure worth fixing, and
`scripts/verify_source_submodules.py` exists to catch it.
