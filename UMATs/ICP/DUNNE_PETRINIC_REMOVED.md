# Dunne & Petrinic book code: removed from the current tree

The folders `elasticity/`, `plasticity_exp/`, `plasticity_imp/`, `spin/` and
`visco/` in this directory, with their `.inp` decks, and `index.html`, held the
companion Fortran code of

> F. Dunne and N. Petrinic, *Introduction to Computational Plasticity*,
> Oxford University Press, 2005 (ISBN 0 19 856826 6), Appendix B,
> "Fortran coding available via the OUP website".

Neither the files nor the publisher's page state a licence, so they cannot be
redistributed. UMAT-OTI removed the same files in its 1.2.0 release (its
`THIRD_PARTY_NOTICES.md`, section 2b); this repository removed its copies on
2026-11-06, for use with UMAT-OTI 1.2.1. The `UMAT_*.for` files in this
directory are copies of jgomezc1/ABAQUS-US (MIT) and are not affected.

Removed (sha256 prefix of each source file):

| File | sha256 |
|---|---|
| `elasticity/elastic.f` | `7c7217e4` |
| `plasticity_exp/code_exp.f` | `efea003a` |
| `plasticity_imp/code_imp.f` | `82defecd` |
| `spin/spin_elas_def.f` | `03458588` |
| `spin/spin_elastic.f` | `d59b9596` |
| `visco/uni_visco_exp.f` | `b0dde38b` |
| `visco/uni_visco_imp.f` | `1bbec9f5` |
| `visco/visco_beam.f`, `visco/visco_imp.f` (identical) | `33663a26` |

with the 24 `.inp` decks of those folders and `index.html`.

**They had already been distributed**: they remain in this repository's git
history. Removing them from the current tree does not withdraw those copies;
rewriting the history is a separate decision.

Nothing in this repository reads these files. The claim-4 benchmark check
(`verification/claim4_benchmark_ddsdde.py`) runs UMAT-OTI's benchmark
contracts, seven of which name these sources at UMAT-OTI's old paths; it
reports them as "source not available" unless a user places the files there.
