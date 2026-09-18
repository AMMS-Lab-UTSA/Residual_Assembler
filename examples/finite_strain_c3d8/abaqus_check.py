"""Prepare or compare one seven-step Abaqus job; never launches a solver."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from residual_core.formulations.c3d8_kernel import unit_cube_Xe
from residual_core.formulations.solid_c3d8_finite_strain import SolidC3D8FiniteStrain
from residual_core.materials.base import MaterialBinding
from residual_core.materials.neo_hookean import CompressibleNeoHookean


def prepare(work):
    work.mkdir(parents=True, exist_ok=True)
    path = work / "imqr6_hyper.inp"
    if path.exists():
        raise FileExistsError("Refusing to overwrite an existing Abaqus job: " + str(path))
    coords = unit_cube_Xe()
    angle = 0.8
    rotation = np.array([[np.cos(angle), -np.sin(angle), 0],
                         [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
    deformation = rotation @ np.array([[1.35, 0.25, -0.1], [0, 0.85, 0.15], [0, 0, 1.12]])
    displacement = (coords @ (deformation - np.eye(3)).T).ravel()
    direction = (coords @ np.array([[0.6, 0.2, 0.1], [-0.1, 0.3, 0], [0.05, 0, -0.2]]).T).ravel()
    binding = MaterialBinding(CompressibleNeoHookean(), [2.3, 4.1])
    force, tangent, _, _ = SolidC3D8FiniteStrain().eval_element(
        1, "C3D8", coords, displacement, {}, None, binding, (0, 1), 1, None, {})
    steps = [("base", displacement)]
    sizes = [1e-3, 3e-4, 1e-4]
    for index, size in enumerate(sizes):
        steps.extend([("plus%d" % index, displacement + size * direction),
                      ("minus%d" % index, displacement - size * direction)])
    lines = ["*Heading", "Bounded log neo-Hookean energy / homogeneous C3D8", "*Node"]
    lines += ["%d, %.16g, %.16g, %.16g" % (index + 1, *point)
              for index, point in enumerate(coords)]
    lines += ["*Element, type=C3D8, elset=solid", "1,1,2,3,4,5,6,7,8",
              "*Nset, nset=allnodes, generate", "1,8,1",
              "*Solid Section, elset=solid, material=nh", ",", "*Material, name=nh",
              "*Hyperelastic, user, type=compressible, properties=2", "2.3,4.1"]
    for name, values in steps:
        lines += ["*Step, name=%s, nlgeom=YES, inc=100" % name, "*Static", "0.1,1.,1e-8,0.1",
                  "*Boundary"]
        lines += ["%d,%d,%d,%.16g" % (index // 3 + 1, index % 3 + 1, index % 3 + 1, value)
                  for index, value in enumerate(values)]
        lines += ["*Output, field, frequency=1", "*Node Output, nset=allnodes", "U,RF", "*End Step"]
    path.write_text("\n".join(lines) + "\n")
    reference = {"force": force.tolist(), "directional_tangent": (tangent @ direction).tolist(),
                 "step_sizes": sizes, "steps": {name: values.tolist() for name, values in steps},
                 "scope": "homogeneous affine full-integration C3D8; energy-based UHYPER"}
    (work / "reference.json").write_text(json.dumps(reference, indent=2) + "\n")
    return {"prepared": str(path), "steps": len(steps)}


def compare(work):
    reference = json.loads((work / "reference.json").read_text())
    fields = json.loads((work / "fields.json").read_text())
    def values(step, field):
        return np.array([fields[step][field][str(node)] for node in range(1, 9)]).ravel()
    def error(actual, expected):
        return float(np.linalg.norm(actual - expected) / np.linalg.norm(expected))
    displacement_error = max(float(np.max(np.abs(values(name, "U") - expected)))
                             for name, expected in reference["steps"].items())
    force_error = error(values("base", "RF"), reference["force"])
    sweep = []
    for index, size in enumerate(reference["step_sizes"]):
        numerical = (values("plus%d" % index, "RF") - values("minus%d" % index, "RF")) / (2 * size)
        sweep.append({"step": size, "relative_error": error(numerical, reference["directional_tangent"])})
    completed = all(abs(fields[name]["time"] - 1) < 1e-10 for name in reference["steps"])
    report = {"job": "imqr6_hyper", "scope": reference["scope"], "completed": completed,
        "displacement_max_error": displacement_error, "reaction_relative_error": force_error,
        "reaction_derivative_fd": sweep, "AMATRX_exported": False,
        "passed": completed and displacement_error < 1e-7 and force_error < 2e-6
                  and max(entry["relative_error"] for entry in sweep) < 2e-5}
    report["artifacts_sha256"] = {name: hashlib.sha256((work / name).read_bytes()).hexdigest()
        for name in ("imqr6_hyper.inp", "imqr6_hyper.odb", "imqr6_hyper.sta", "fields.json", "reference.json")}
    source = Path(__file__).with_name("neo_hookean_uhyper.for")
    report["uhyper_sha256"] = hashlib.sha256(source.read_bytes()).hexdigest()
    (work / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "compare"))
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    report = prepare(args.work) if args.action == "prepare" else compare(args.work)
    print(json.dumps(report, indent=2))
    raise SystemExit(0 if report.get("passed", True) else 1)