# Local-model analysis of this repository

Repetitive analysis can run on a local model on this machine rather than on a
metered reviewer. The orchestration stack is **not here** — it lives once, in
The Council — and this file is the thin project-specific part.

- Infrastructure, tiers, benchmark and safety model:
  `../The-council/docs/local/README.md`
- Start and stop: `../The-council/scripts/local/council-local up|status|down`
- What a model may read here: [`.councilignore`](../.councilignore)

## What is delegated from this repository

Reading and explaining, not deciding:

| task | what it reads |
| --- | --- |
| convention audit | index ranges, Voigt ordering, transpose and sign conventions in assembly code |
| interface contracts | what a provider must supply and what the assembler assumes |
| failure triage | a solver or verification log, to name the first cause |
| documentation first pass | an explanation of an existing, already-verified result |

A convention audit is a reading task with a right answer in the source, which
is why it is delegable. Everything downstream of it is not.

## What is not delegated

- **The homogeneity (Euler) identity.** It is this repository's independent
  check — it held to 1.4e-12 of the largest term against a 1e-10 bound on the
  full cantilever — and an independent check that a language model has
  touched is no longer independent.
- **Any residual or Jacobian assembly change.** A local model may propose a
  patch and must cite file and line; Vera reviews every material change to
  mathematical assembly, and that review is not a local-model task.
- **Any sensitivity number.** dR/dp, dσ/dp and dstatev/dp are established
  against an analytical or finite-difference reference. The reference is the
  authority, never a model's reading of the code.

## The rule that matters here

A local model is asked to *read*, never to *recall*. All three models on this
machine, asked from memory for the 37 arguments of the UMAT interface,
answered confidently and wrongly. Prompts quote the source with line numbers;
every claim carries a `path:line` a reviewer can open; an answer whose
citations do not resolve escalates itself rather than being returned.
