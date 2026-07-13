# Template — C++ global residual as a black-box executable

Build a `your_solver` that speaks the kit's black-box JSON contract:

```
your_solver --input request.json --output response.json
```

`main.cpp` shows the shape: read the request, seed each parameter with
`resasm::Dual1`, evaluate your **templated** global residual, extract `R^(1)`, and
write `response.json` with `resasm::write_response_json`. Your model stays inside
the binary; only the coefficients you write are exposed.

> The parsing in `main.cpp` is a **minimal, illustrative** number extractor so the
> template is dependency-free. For production, use a real JSON library
> (e.g. nlohmann/json) — the schema is in
> [../../docs/residual_provider_contract.md](../../docs/residual_provider_contract.md).

## Build

```bash
c++ -std=c++17 -I../../include main.cpp -o your_solver
```

## Use from the kit

```json
{ "blackbox": { "command": ["./your_solver"], "ndof": 1,
                "parameters": ["k"], "parameter_values": {"k": 2.0, "f": 16.0} },
  "parameters": ["k"], "order": 1, "solution": [2.0], "output_dir": "out" }
```
