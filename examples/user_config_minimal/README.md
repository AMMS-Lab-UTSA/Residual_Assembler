# Minimal user config example

The smallest complete job: a 1-DOF cubic spring `R = k·u³ − f`, sensitivities of
`u` w.r.t. `k` and `f`.

```
resasm check resasm.yml
resasm run   resasm.yml
resasm report resasm_output/
```

Order 1: `du/dk = −1/3`, `du/df = 1/24`. Order 2 needs OTILib
(see [QUICKSTART_USER.md](../../QUICKSTART_USER.md)); the finite-difference
cross-check confirms order 1 without it.
