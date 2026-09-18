#!/usr/bin/env python3
"""Write the two example cantilever decks as Abaqus input files.

Two models:

  j2  : 1,536 C3D8, 2,499 nodes (7,497 DOF), 12,288 integration points,
        J2 with linear hardening E=200 GPa, nu=0.3, sigy=250 MPa, H=2 GPa,
        clamped root, tip pushed 0.7 mm over 40 load steps.
  fcc : 384 C3D8, 675 nodes (2,025 DOF), 3,072 integration points,
        FCC single crystal, 12 {111}<110> slip systems (m6_fcc UMAT),
        clamped root, tip pushed 0.25 mm over 25 load steps.

Both counts are met exactly by a structured nx x ny x nz grid of unit cubes:
48 x 16 x 2 for j2 and 24 x 8 x 2 for fcc. Units: mm, N, MPa. Length along x,
height along y (the push direction), thickness along z. Small strain
(NLGEOM=NO). The root face x=0 is clamped (ENCASTRE); every node of the tip
face x=L is displaced in -y by the prescribed amount, ramped linearly with
exactly one increment per load step (DIRECT). Every increment is written to
the ODB. PROPS can be overridden for the finite-difference reruns.

    python gen_cantilever.py j2  --out j2/cantilever_j2_nominal.inp
    python gen_cantilever.py j2  --props 200000,0.3,250,2000 --out ...
"""
import argparse
import sys

MODELS = {
    "j2": {"n": (48, 16, 2), "push": 0.7, "steps": 40, "depvar": 1,
           "props": [200000.0, 0.3, 250.0, 2000.0],
           "prop_names": ["E", "nu", "SIGY0", "H"], "material": "J2"},
    "fcc": {"n": (24, 8, 2), "push": 0.25, "steps": 25, "depvar": 12,
            # C11-C12-C44 = 168-121-75 GPa, g0 = 13 MPa, gsat = 55,
            # h0 = 800, a = 2, q = 1.4, gamma0 = 0.001 /s, m = 0.05 (MPa, s).
            # (The model's contract_v2.json validation point differs: g0 16,
            # gsat 40, h0 300.)
            "props": [168000.0, 121000.0, 75000.0, 13.0, 55.0, 800.0, 2.0,
                      1.4, 0.001, 0.05],
            "prop_names": ["C11", "C12", "C44", "g0", "gsat", "h0", "a", "q",
                           "gd0", "m"], "material": "FCC"},
}


def node_id(i, j, k, nx, ny):
    return 1 + i + (nx + 1) * (j + (ny + 1) * k)


def deck(model, props):
    spec = MODELS[model]
    nx, ny, nz = spec["n"]
    lines = ["*Heading",
             f"{model} cantilever: {nx}x{ny}x{nz} C3D8, "
             f"tip -{spec['push']} mm in {spec['steps']} steps",
             "*Node"]
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                lines.append(f"{node_id(i, j, k, nx, ny)}, {float(i)}, {float(j)}, {float(k)}")
    lines.append("*Element, type=C3D8, elset=ALL")
    eid = 0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                eid += 1
                n = lambda a, b, c: node_id(a, b, c, nx, ny)
                conn = [n(i, j, k), n(i + 1, j, k), n(i + 1, j + 1, k), n(i, j + 1, k),
                        n(i, j, k + 1), n(i + 1, j, k + 1), n(i + 1, j + 1, k + 1), n(i, j + 1, k + 1)]
                lines.append(f"{eid}, " + ", ".join(str(c) for c in conn))
    root = [node_id(0, j, k, nx, ny) for k in range(nz + 1) for j in range(ny + 1)]
    tip = [node_id(nx, j, k, nx, ny) for k in range(nz + 1) for j in range(ny + 1)]

    def nset(name, ids):
        out = [f"*Nset, nset={name}"]
        for s in range(0, len(ids), 16):
            out.append(", ".join(str(x) for x in ids[s:s + 16]))
        return out

    lines += nset("ROOT", root) + nset("TIP", tip)
    # the tip node at mid-height of the z=0 face: the scalar output q = U2 there
    lines += nset("TIPMID", [node_id(nx, ny // 2, 0, nx, ny)])
    lines += [f"*Elset, elset=ROOTEL, generate", f"1, {1 + nx * (ny * nz - 1)}, {nx}"]
    lines += [f"*Material, name={spec['material']}", "*Depvar", str(spec["depvar"]),
              f"*User Material, constants={len(props)}"]
    for s in range(0, len(props), 8):
        lines.append(", ".join(repr(float(p)) for p in props[s:s + 8]))
    lines += [f"*Solid Section, elset=ALL, material={spec['material']}", ",",
              "*Boundary", "ROOT, ENCASTRE",
              f"*Step, name=Push, nlgeom=NO, inc={spec['steps'] + 10}",
              "*Static, direct",
              f"{1.0 / spec['steps']!r}, 1.0",
              # The m6_fcc UMAT returns the ELASTIC stiffness as DDSDDE (its
              # header says the OTI transform is what replaces it), so Newton
              # converges linearly once slip starts; with fixed increments
              # Abaqus cannot cut back, so the iteration limits are raised
              # instead of letting the increment be subdivided (the replay
              # needs exactly the stated load steps).
              *( ["*Controls, parameters=time incrementation",
                  "50, 60, , 400, , , , , , "] if model == "fcc" else [] ),
              "*Boundary", f"TIP, 2, 2, {-spec['push']!r}",
              "*Output, field, frequency=1",
              "*Node Output", "U, RF",
              "*Element Output, position=INTEGRATION POINTS", "S, E, SDV",
              "*Output, history, frequency=1", "*Node Output, nset=TIPMID", "U2",
              "*End Step"]
    return "\n".join(lines) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("model", choices=sorted(MODELS))
    ap.add_argument("--props", help="comma-separated PROPS overriding the nominal values")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    props = MODELS[a.model]["props"]
    if a.props:
        props = [float(x) for x in a.props.split(",")]
        if len(props) != len(MODELS[a.model]["props"]):
            sys.exit(f"{a.model} needs {len(MODELS[a.model]['props'])} PROPS, got {len(props)}")
    with open(a.out, "w") as fh:
        fh.write(deck(a.model, props))
    nx, ny, nz = MODELS[a.model]["n"]
    print(f"wrote {a.out}: {nx*ny*nz} C3D8, {(nx+1)*(ny+1)*(nz+1)} nodes, "
          f"{3*(nx+1)*(ny+1)*(nz+1)} DOF, {8*nx*ny*nz} IPs, PROPS={props}")


if __name__ == "__main__":
    main()
