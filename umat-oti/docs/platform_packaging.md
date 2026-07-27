# Platform-specific deliverables

**The source and the CLI are cross-platform. A compiled binary is not.** One
object file cannot serve both Windows and Linux — they use different object
formats (COFF/PE vs ELF), different runtime libraries, and different linkers.
This page states exactly what crosses the boundary and what is platform-bound.

## What crosses the boundary

```
   ┌───────────────────────────── per platform ─────────────────────────────┐
   │  Windows deliverable                    Linux deliverable               │
   │  ─────────────────────                  ────────────────────            │
   │  umat_<name>_oti.obj   (COFF)           umat_<name>_oti.obj   (ELF .o)   │
   │    links into a .dll                      links into a .so              │
   │  built with: gfortran (MinGW-w64)       built with: gfortran (GNU/Linux)│
   └─────────────────────────────────────────────────────────────────────────┘
                                    +
   ┌───────────────────── platform-independent ──────────────────────────────┐
   │  umat_<name>_oti.json   (the completed contract)                         │
   │    dimensions, parameters, symbols, layouts, ABI version, validation,    │
   │    AND a `binary` metadata block describing the object it was built with  │
   └──────────────────────────────────────────────────────────────────────────┘
```

- **The JSON contract is shared.** It is the same bytes on every platform and
  carries no machine code. It *describes* the binary (including a `binary`
  metadata block, below) but is not itself platform-bound.
- **The compiled object is per-platform.** A collaborator on Linux needs the
  object produced by the Linux toolchain; the Windows object will not link there,
  and vice-versa. Ship the object built for the collaborator's platform (or ship
  both, alongside the one shared JSON).

Note the file *extension* is not the discriminator — the *format* is. Program 2
sometimes names a linked library `libmat.so` on both platforms; on Windows that
file is a PE. Compatibility is decided from the object's actual magic bytes, not
its name.

## Distributable object vs temporary validation library

Two different things are easy to confuse:

| | distributable material object | temporary validation library |
|---|---|---|
| **file** | `umat_<name>_oti.obj` (COFF/ELF, relocatable) | `orig.{dll,so}`, `oti.{dll,so}` under `build/` |
| **who makes it** | Program 1, shipped to collaborators | Program 1 **and** Program 2, internally |
| **purpose** | the product; Program 2 links it into a loadable lib | to call the material through `ctypes` for a check |
| **lifetime** | permanent, versioned, hashed in the contract | throwaway; regenerated every run, git-ignored |
| **crosses the boundary?** | **yes** (per platform) | **never** |

The temporary `.dll`/`.so` files are an implementation detail of the in-process
validation drivers. They are **not** the deliverable and are never shipped; only
the relocatable `.obj` and the JSON contract are.

## The `binary` metadata block

Program 1 records, in the completed contract (and in the material-package
manifest), exactly what the object is and how it was built:

```json
"binary": {
  "schema": "resasm_binary_metadata_v1",
  "os": "windows",              "arch": "amd64",        "pointer_bits": 64,
  "compiler": "gfortran",       "compiler_version": "15.2.0",
  "binary_format": "coff/pe",   "object_format_probed": "COFF",
  "loadable_suffix": ".dll",
  "abi_version": "<contract hash>",
  "build_id": "<deterministic id>",
  "source_hash": "<umat.for sha>",   "transform_hash": "<transform sha>",
  "object_file": "umat_<name>_oti.obj",  "object_sha256": "<obj sha>"
}
```

Fields: **operating system**, **architecture**, **compiler + version**, **binary
format** (and the format actually probed from the file), **ABI/interface
version**, a deterministic **build identifier**, and the **source / transform
hashes**.

## Program 2 rejects an incompatible package before linking

Before Program 2 invokes the linker on a distributed object it calls
`residual_core.runtime.check_binary_compatibility`, which:

1. probes the object's **actual format** from its magic bytes (ELF / COFF / PE);
2. checks it belongs to **this platform's format family** (Windows: COFF/PE;
   Linux: ELF);
3. cross-checks the declared `binary` metadata (os / arch / pointer width) when
   present.

A mismatch raises a clear diagnostic naming the produced platform, the current
platform, and the format conflict — e.g.

```
material object is not compatible with this platform.
  object format     : ELF (from magic bytes)
  this platform     : os=windows arch=amd64 pointer=64-bit, uses COFF/PE
  problem(s)        : object format is ELF but this platform (windows) uses COFF/PE
  -> obtain the material object built by the linux toolchain, or rebuild it there.
     The JSON contract is platform-independent; the compiled object is not.
```

— instead of an opaque `ld` error partway through the link.
