#!/usr/bin/env python3
"""Write a structured C3D8 cantilever deck for the history replay.

The mesh is nx x ny x nz unit cubes (x length, y height = push direction,
z thickness; mm, N, MPa). The root face x=0 is clamped (ENCASTRE) in the model
data; every node of the tip face x=nx gets U2 = -push inside one static step
(NLGEOM=NO), ramped over ``steps`` fixed increments (*Static, direct). Every
increment writes U, RF, S and SDV. Optional ``--cload`` replaces the pushed
tip by a concentrated U2 force on the TIPMID node (then CF is written too).

    python make_beam_deck.py --n 12 4 2 --push 0.08 --steps 10 \
        --props 200000,0.3,250,2000 --depvar 1 --out beam.inp
"""
import argparse


def node_id(i, j, k, nx, ny):
    return 1 + i + (nx + 1) * (j + (ny + 1) * k)


def deck(n, push, steps, props, depvar, cload=None, controls=False):
    nx, ny, nz = n
    lines = ["*Heading", "history replay beam %dx%dx%d C3D8, %d steps" % (nx, ny, nz, steps), "*Node"]
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                lines.append("%d, %r, %r, %r" % (node_id(i, j, k, nx, ny), float(i), float(j), float(k)))
    lines.append("*Element, type=C3D8, elset=ALL")
    eid = 0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                eid += 1
                c = [node_id(a, b, d, nx, ny) for a, b, d in (
                    (i, j, k), (i + 1, j, k), (i + 1, j + 1, k), (i, j + 1, k),
                    (i, j, k + 1), (i + 1, j, k + 1), (i + 1, j + 1, k + 1), (i, j + 1, k + 1))]
                lines.append("%d, %s" % (eid, ", ".join(str(x) for x in c)))

    def nset(name, ids):
        out = ["*Nset, nset=%s" % name]
        out += [", ".join(str(x) for x in ids[s:s + 16]) for s in range(0, len(ids), 16)]
        return out

    root = [node_id(0, j, k, nx, ny) for k in range(nz + 1) for j in range(ny + 1)]
    tip = [node_id(nx, j, k, nx, ny) for k in range(nz + 1) for j in range(ny + 1)]
    lines += nset("ROOT", root) + nset("TIP", tip) + nset("TIPMID", [node_id(nx, ny // 2, 0, nx, ny)])
    lines += ["*Elset, elset=ROOTEL, generate", "1, %d, %d" % (1 + nx * (ny * nz - 1), nx)]
    lines += ["*Material, name=UMATMAT", "*Depvar", str(depvar),
              "*User Material, constants=%d" % len(props)]
    lines += [", ".join(repr(float(p)) for p in props[s:s + 8]) for s in range(0, len(props), 8)]
    lines += ["*Solid Section, elset=ALL, material=UMATMAT", ",", "*Boundary", "ROOT, ENCASTRE",
              "*Step, name=Push, nlgeom=NO, inc=%d" % (steps + 10), "*Static, direct",
              "%r, 1.0" % (1.0 / steps)]
    if controls:
        lines += ["*Controls, parameters=time incrementation", "50, 60, , 400, , , , , , "]
    if cload is None:
        lines += ["*Boundary", "TIP, 2, 2, %r" % (-push)]
    else:
        lines += ["*Cload", "TIPMID, 2, %r" % float(cload)]
    lines += ["*Output, field, frequency=1", "*Node Output", "U, RF" + (", CF" if cload is not None else ""),
              "*Element Output, position=INTEGRATION POINTS", "S, SDV", "*End Step"]
    return "\n".join(lines) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", nargs=3, type=int, required=True, metavar=("NX", "NY", "NZ"))
    ap.add_argument("--push", type=float, default=0.0)
    ap.add_argument("--cload", type=float)
    ap.add_argument("--steps", type=int, required=True)
    ap.add_argument("--props", required=True)
    ap.add_argument("--depvar", type=int, required=True)
    ap.add_argument("--controls", action="store_true", help="raise Abaqus iteration limits")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    props = [float(x) for x in a.props.split(",")]
    with open(a.out, "w") as stream:
        stream.write(deck(a.n, a.push, a.steps, props, a.depvar, a.cload, a.controls))


if __name__ == "__main__":
    main()
