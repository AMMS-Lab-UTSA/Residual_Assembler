"""``resasm`` command line interface.

Progressive disclosure: the common path is ``resasm inspect model.inp`` and
``resasm assemble model.inp --mode stress-driven --fields fields.json``. A config
file is optional and only needed for overrides; ``resasm doctor`` can emit a
pre-filled template.

Commands
--------
  inspect  MODEL                         model inspection summary (elements,
                                         materials, possible modes, missing data)
  requirements MODEL --mode M            what MODE needs vs. what is present
  assemble MODEL --mode M [--fields F]   assemble the global residual
           [--subroutine S] [--tangent]
  verify   MODEL --fields F              stress-driven residual + reaction check
  doctor   MODEL [--write-config-template PATH]
                                         diagnose + optionally emit a config
  template --formulation NAME            print a backend's declared contract
  backends                               user-readable registry audit (all backends)
  sensitivity MODEL --params P           generate R^(1) and solve T U^(1) = -R^(1)
  modes                                  list assembly modes
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

from . import wizard as _wiz
from .config import Config, write_config_template
from ..core import requirements as _req
from ..formulations.registry import build_formulation_registry
from ..materials.registry import build_material_registry


def _load_problem(model_path: str, config_path=None):
    cfg = Config.load(config_path) if config_path else None
    if model_path.lower().endswith((".json",)):
        prob = _wiz.ResidualProblem.from_neutral(model_path, config=cfg)
    else:
        prob = _wiz.ResidualProblem.from_abaqus(model_path, config=cfg)
    if cfg:
        if cfg.subroutine:
            prob.attach_subroutine(cfg.subroutine)
        if cfg.odb:
            prob.attach_results(cfg.odb)
    return prob, cfg


def _attach_common(prob, args):
    """Attach a field export (--fields / --odb) and/or subroutine if given."""
    fields = getattr(args, "fields", None) or getattr(args, "odb", None)
    if fields:
        prob.attach_results(fields)
    if getattr(args, "subroutine", None):
        prob.attach_subroutine(args.subroutine)
    return prob


def _cmd_inspect(args):
    prob, _ = _load_problem(args.model, args.config)
    rep = prob.inspect(print_it=True)
    if getattr(args, "detail", False):
        print()
        print(rep.render_backends())
    return 0


def _cmd_requirements(args):
    prob, _ = _load_problem(args.model, args.config)
    _attach_common(prob, args)
    print(prob.requirements(args.mode).render())
    return 0


def _cmd_assemble(args):
    prob, _ = _load_problem(args.model, args.config)
    _attach_common(prob, args)
    # progressive disclosure: if data is missing, print the clean requirements
    # report (minimum missing input only) instead of crashing.
    report = prob.requirements(args.mode)
    if not report.runnable:
        print(report.render(), file=sys.stderr)
        return 2
    try:
        out = prob.assemble(mode=args.mode, compute_tangent=args.tangent)
    except (RuntimeError, NotImplementedError) as exc:
        print("Cannot assemble in %s mode." % _req.canonical_mode(args.mode),
              file=sys.stderr)
        print("  %s" % exc, file=sys.stderr)
        return 2
    R = out[0] if isinstance(out, tuple) else out
    print("assembled mode '%s': ndof=%d  ||R||=%.6e  max|R|=%.6e"
          % (_req.canonical_mode(args.mode), R.size,
             float(np.linalg.norm(R)), float(np.max(np.abs(R))) if R.size else 0.0))
    if args.out:
        np.save(args.out, R)
        print("residual saved to %s" % args.out)
    return 0


def _cmd_verify(args):
    prob, _ = _load_problem(args.model, args.config)
    fields = getattr(args, "fields", None) or getattr(args, "odb", None)
    if fields:
        prob.attach_results(fields)
    report = prob.requirements("stress-driven")
    if not report.runnable:
        print(report.render(), file=sys.stderr)
        return 2
    R = prob.assemble(mode="stress-driven")
    free_R, pres_idx, reac = _split_safe(prob, R)
    print("stress-driven verification:")
    print("  ndof                 = %d" % R.size)
    print("  ||R_free||           = %.6e   (should be ~0 at equilibrium)"
          % float(np.linalg.norm(free_R)))
    print("  ||reaction (at BC)|| = %.6e" % float(np.linalg.norm(reac)))
    return 0


def _split_safe(prob, R):
    try:
        from ..core import constraints as _c
        free_mask, pres_idx, _ = _c.partition(prob.model, prob.dof_manager)
        return R[free_mask], pres_idx, (R[pres_idx] if pres_idx.size else np.zeros(0))
    except Exception:
        return R, np.zeros(0, int), np.zeros(0)


def _cmd_doctor(args):
    prob, _ = _load_problem(args.model, args.config)
    rep = prob.inspect(print_it=True)
    print("\nper-mode readiness:")
    for mode in _req.known_modes():
        r = prob.requirements(mode)
        status = "ready" if r.runnable else ("needs: %s" % r.minimum_next.display())
        print("  %-16s %s" % (mode, status))
    if args.write_config_template:
        detected = Config(subroutine=(prob._subroutine))
        write_config_template(args.write_config_template, detected)
        print("\nconfig template written to %s" % args.write_config_template)
    return 0


def _print_spec(d, name):
    print("backend: %s" % d.get("backend_name", name))
    for key in ("kind", "verification_status", "supported_element_types",
                "dof_types", "supported_modes", "material_interface_needed",
                "state_requirements", "tangent_support", "required_inputs",
                "optional_inputs", "verification_tests", "limitations", "notes"):
        if key in d and (d[key] or d[key] is False):
            print("  %-26s %s" % (key + ":", d[key]))
    rbm = d.get("required_inputs_by_mode") or {}
    if rbm:
        print("  required_inputs_by_mode:")
        for mode, items in rbm.items():
            print("      %-18s %s" % (mode + ":", items))


def _cmd_template(args):
    forms = build_formulation_registry()
    mats = build_material_registry()
    name = args.formulation or args.material
    reg = forms if args.formulation else mats
    if name not in reg.names():
        print("unknown backend %r; known formulations: %s; materials: %s"
              % (name, ", ".join(forms.names()), ", ".join(mats.names())),
              file=sys.stderr)
        return 2
    _print_spec(reg.get(name).spec.as_dict(), name)
    return 0


def _cmd_backends(args):
    """User-readable registry audit: every formulation backend's declaration."""
    forms = build_formulation_registry()
    print("Registered formulation backends")
    print("===============================")
    for name in forms.names():
        spec = forms.get(name).spec
        modes = ", ".join(spec.supported_modes) or "none"
        print("\n%s" % name)
        print("  element types : %s" % ", ".join(spec.supported_element_types))
        print("  dof types     : %s" % ", ".join(spec.dof_types))
        print("  status        : %s" % spec.verification_status)
        print("  available modes: %s" % modes)
        print("  material iface : %s" % spec.material_interface_needed)
        print("  state         : %s" % spec.state_requirements)
        print("  tangent       : %s" % spec.tangent_support)
        if spec.limitations:
            print("  limitations   : %s" % "; ".join(spec.limitations))
    print("\nMaterials: %s"
          % ", ".join(build_material_registry().names()))
    return 0


def _cmd_modes(args):
    print("assembly modes:")
    for mode in _req.known_modes():
        items = _req.MODE_REQUIREMENTS[mode]
        print("  %-16s min inputs: %s" % (mode, ", ".join(i.key for i in items)))
    return 0


def _load_params(path):
    """Load a params file: JSON always, YAML if PyYAML is installed. Returns a
    dict with keys: parameters (list), mode, order, expected (optional map)."""
    import json as _json
    with open(path, "r", encoding="utf-8") as fh:
        text = fh.read()
    if path.lower().endswith((".yml", ".yaml")):
        try:
            import yaml
            return yaml.safe_load(text) or {}
        except Exception:
            pass
    return _json.loads(text)


def _cmd_sensitivity(args):
    prob, _ = _load_problem(args.model, args.config)
    spec = _load_params(args.params) if args.params else {}
    params = list(args.param) if args.param else spec.get("parameters")
    mode = args.mode or spec.get("mode") or "formulation"
    order = args.order or int(spec.get("order", 1))
    backend = args.backend or spec.get("backend") or "otilib"
    expected = spec.get("expected", {}) or {}

    try:
        pkg = prob.sensitivity_package(mode=mode, parameters=params, max_order=order,
                                       generate_rhs=True, fd_check=not args.no_fd,
                                       backend=backend)
    except RuntimeError as exc:
        # OTILib missing (or backend error) -> report cleanly, do NOT fall back
        print("hypercomplex backend: %s" % backend, file=sys.stderr)
        print("ERROR: %s" % exc, file=sys.stderr)
        if backend.lower().startswith(("oti", "hypad")):
            print("\nInstall OTILib from https://github.com/mauriaristi/otilib.git "
                  "(GPLv3) — run scripts/setup_otilib.sh or see "
                  "docs/otilib_integration.md. For a first-order-only smoke test "
                  "you may pass --backend dual1 (NOT the production path).",
                  file=sys.stderr)
        return 3
    if not pkg.runnable:
        print("Cannot assemble in %s mode." % mode, file=sys.stderr)
        print("  minimum missing input: %s" % pkg.minimum_missing, file=sys.stderr)
        return 2

    s = pkg.sensitivity
    print("hypercomplex backend     : %s" % s.algebra.algebra)
    print("basis count (m)          : %d" % s.algebra.n_bases)
    print("truncation order (nt)    : %d" % s.algebra.truncation_order)
    print("total coefficients (N)   : %d" % s.algebra.n_coefficients())
    print("sensitivity analysis (mode=%s, order=%d)" % (mode, order))
    print("  residual norm ||R_free|| = %.6e" % pkg.residual.free_residual_norm)
    print("  tangent source           = %s" % pkg.residual.tangent.source)
    print("  parameters               = %s" % list(s.parameter_map))
    for p in range(1, order + 1):
        print("  R^(%d) shape              = %s  (%d directions)"
              % (p, s.expected_shape(p), s.algebra.n_directions(p)))
    print("  hypercomplex ready       = %s"
          % s.diagnostics.get("hypercomplex_ready", "n/a"))
    if not pkg.system_ready(1):
        print("  system NOT ready: tangent unavailable or R^(1) not generated "
              "(stress-driven residual assembly alone cannot produce parameter "
              "sensitivities — use a parameterized formulation/material backend).",
              file=sys.stderr)
        return 2

    fd = None
    if not args.no_fd:
        try:
            fd = prob.finite_difference_sensitivity(mode, list(s.parameter_map))
        except Exception:
            fd = None
    print("  solved derivative orders : %s"
          % [p for p in range(1, order + 1) if s.R(p) is not None])
    for p in range(1, order + 1):
        if s.R(p) is None:
            continue
        Up = pkg.solve(p)
        dirs = s.direction_map(p)
        print("  order %d:" % p)
        for col, d in enumerate(dirs):
            up = Up[:, col]
            raw = float(up[0]) if up.size == 1 else float(np.linalg.norm(up))
            deriv = d["recovery_factor"] * raw     # true partial derivative
            label = d["label"]
            line = "    d^%d/%-10s = %+.6e" % (p, label, deriv)
            key = _expected_key(label, d, expected)
            if key is not None:
                ex = float(expected[key])
                rel = abs(deriv - ex) / max(abs(ex), 1e-30)
                line += "   [analytic %+.6e, rel %.2e]" % (ex, rel)
            if p == 1 and fd is not None:
                fv = fd[:, col]
                fs = float(fv[0]) if fv.size == 1 else float(np.linalg.norm(fv))
                rel = abs(raw - fs) / max(abs(fs), 1e-30)
                line += "   [FD %+.6e, rel %.2e]" % (fs, rel)
            print(line)
    if args.out:
        pkg.save(args.out)
        print("  saved package -> %s.{json,npz,md}" % args.out)
    return 0


def _expected_key(label, direction, expected):
    """Match a direction to an entry in the params-file 'expected' map. Accepts
    either the OTI label (e.g. 'e1^2') or a parameter-derivative spelling."""
    if label in expected:
        return label
    params = direction.get("parameters") if isinstance(direction, dict) else None
    if params:
        joined = "*".join(params)
        if joined in expected:
            return joined
    return None


# --------------------------------------------------------------------------- #
# user-facing commands (resasm init / check / run / report) -> resasm_user
# --------------------------------------------------------------------------- #
# Short user-facing name -> template folder under templates/
_TEMPLATES = {
    "python":  "user_python_residual",
    "blackbox": "user_blackbox_residual",
    "cpp":     "user_cpp_residual",
    "fortran": "user_fortran_residual",
}


def _templates_root():
    """Locate the templates/ directory (repo root; works for `pip install -e .`)."""
    here = os.path.dirname(os.path.abspath(__file__))            # residual_core/ui
    root = os.path.dirname(os.path.dirname(here))                # repo root
    for cand in (os.path.join(root, "templates"),
                 os.path.join(os.getcwd(), "templates")):
        if os.path.isdir(cand):
            return cand
    return None


def _cmd_init(args):
    # Non-interactive: copy a ready-to-run template folder. This is the path real
    # users want -- `resasm init --template python --out my_job`.
    if getattr(args, "template", None):
        import shutil
        root = _templates_root()
        if root is None:
            print("ERROR: could not locate the templates/ directory.\n"
                  "Run from a checkout, or `pip install -e .` from the repo root.",
                  file=sys.stderr)
            return 2
        src = os.path.join(root, _TEMPLATES[args.template])
        if not os.path.isdir(src):
            print("ERROR: template not found: %s" % src, file=sys.stderr)
            return 2
        dest = args.out or ("resasm_%s" % args.template)
        if os.path.exists(dest):
            if not args.force:
                print("ERROR: %s already exists (use --force to overwrite)" % dest,
                      file=sys.stderr)
                return 2
            shutil.rmtree(dest)
        # never copy build/run artifacts a previous run may have left in the
        # template folder (a stale resasm_output/ would look like real results)
        shutil.copytree(src, dest,
                        ignore=shutil.ignore_patterns(
                            "__pycache__", "*.pyc", "resasm_output", "build",
                            "*.o", "*.mod", "*.exe"))
        print("created %s (from template '%s')" % (dest, args.template))
        for f in sorted(os.listdir(dest)):
            print("  %s" % f)
        print("\nNext:\n  cd %s\n  resasm check resasm.yml\n  resasm run resasm.yml"
              % dest)
        return 0

    # Interactive wizard (still supported).
    from resasm_user.wizard import init_wizard
    init_wizard(out_path=args.out or "resasm.yml")
    print("\nNext:\n  resasm check %s\n  resasm run %s"
          % (args.out or "resasm.yml", args.out or "resasm.yml"))
    return 0


def _cmd_check(args):
    from resasm_user import check_config
    rep = check_config(args.config_file)
    print(rep.render())
    return 0 if rep.ok else 1


def _cmd_run(args):
    from resasm_user import run_from_config, ConfigError
    try:
        res = run_from_config(args.config_file)
    except ConfigError as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        return 2
    print("run complete: %s" % res.summary)
    print("  private outputs -> %s" % res.private_dir)
    print("  public  outputs -> %s" % res.public_dir)
    print("  read: %s" % os.path.join(res.public_dir, "summary.md"))
    return 0


def _cmd_report(args):
    from resasm_user import read_report
    rep = read_report(args.output_dir)
    meta = rep.get("metadata", {})
    vs = rep.get("validation_summary", {})
    print("Sensitivity run report: %s" % args.output_dir)
    print("  residual norm (free) : %s" % vs.get("residual_free_norm"))
    print("  tangent source       : %s" % vs.get("tangent_source"))
    print("  parameters           : %s" % ", ".join(meta.get("parameters", [])))
    print("  derivative orders    : %s" % vs.get("orders_solved"))
    print("  validation status    : %s" % vs.get("status"))
    print("  private outputs      : %s" % rep["private_dir"])
    print("  public  outputs      : %s" % rep["public_dir"])
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="resasm",
                                description="Model-agnostic residual assembly framework")
    p.add_argument("--config", help="optional config file (.yml/.json)")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("inspect", help="model inspection summary")
    s.add_argument("model")
    s.add_argument("--detail", action="store_true",
                   help="show per-element backend selection detail")
    s.set_defaults(func=_cmd_inspect)

    s = sub.add_parser("requirements", help="what a mode needs vs. what is present")
    s.add_argument("model")
    s.add_argument("--mode", required=True)
    s.add_argument("--fields", help="exported element field (JSON) for stress-driven")
    s.add_argument("--odb", help="alias for --fields")
    s.add_argument("--subroutine")
    s.set_defaults(func=_cmd_requirements)

    s = sub.add_parser("assemble", help="assemble the global residual")
    s.add_argument("model")
    s.add_argument("--mode", default="stress-driven")
    s.add_argument("--fields", help="exported element field (JSON) for stress-driven")
    s.add_argument("--odb", help="alias for --fields")
    s.add_argument("--subroutine", help="UMAT/UEL source for material-replay")
    s.add_argument("--tangent", action="store_true")
    s.add_argument("--out", help="save residual to .npy")
    s.set_defaults(func=_cmd_assemble)

    s = sub.add_parser("verify", help="stress-driven residual + reaction check")
    s.add_argument("model")
    s.add_argument("--fields", help="exported element field (JSON)")
    s.add_argument("--odb", help="alias for --fields")
    s.set_defaults(func=_cmd_verify)

    s = sub.add_parser("doctor", help="diagnose readiness; optionally emit config")
    s.add_argument("model")
    s.add_argument("--write-config-template", dest="write_config_template")
    s.set_defaults(func=_cmd_doctor)

    s = sub.add_parser("template", help="print a backend's declared contract")
    g = s.add_mutually_exclusive_group(required=True)
    g.add_argument("--formulation")
    g.add_argument("--material")
    s.set_defaults(func=_cmd_template)

    s = sub.add_parser("backends", help="user-readable registry audit (all backends)")
    s.set_defaults(func=_cmd_backends)

    s = sub.add_parser("sensitivity",
                       help="generate sensitivity RHS and solve T U^(1) = -R^(1)")
    s.add_argument("model")
    s.add_argument("--params", help="params file (.json/.yml): parameters, mode, order, expected")
    s.add_argument("--param", action="append", help="a design parameter 'material.key' (repeatable)")
    s.add_argument("--mode", help="assembly mode (default: formulation or params file)")
    s.add_argument("--order", type=int, help="derivative order (OTILib supports >=1)")
    s.add_argument("--backend", help="hypercomplex backend: otilib (default, production) "
                   "or dual1 (legacy order-1 smoke test)")
    s.add_argument("--no-fd", action="store_true", help="skip finite-difference cross-check")
    s.add_argument("--out", help="save the sensitivity package to PREFIX.{json,npz,md}")
    s.set_defaults(func=_cmd_sensitivity)

    s = sub.add_parser("modes", help="list assembly modes")
    s.set_defaults(func=_cmd_modes)

    # ---- user-facing minimal workflow (resasm_user) ----
    s = sub.add_parser("init",
                       help="start a job: copy a template (--template) or run the "
                            "interactive wizard")
    s.add_argument("--template", choices=sorted(_TEMPLATES),
                   help="copy a ready-to-run template folder (non-interactive)")
    s.add_argument("--out", help="output dir for --template (e.g. my_job); "
                                 "otherwise the wizard's config path (default resasm.yml)")
    s.add_argument("--force", action="store_true",
                   help="overwrite an existing --out directory")
    s.set_defaults(func=_cmd_init)

    s = sub.add_parser("check", help="check a resasm.yml is ready to run")
    s.add_argument("config_file")
    s.set_defaults(func=_cmd_check)

    s = sub.add_parser("run", help="run a sensitivity job from resasm.yml")
    s.add_argument("config_file")
    s.set_defaults(func=_cmd_run)

    s = sub.add_parser("report", help="summarize a completed run's output dir")
    s.add_argument("output_dir")
    s.set_defaults(func=_cmd_report)
    return p


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
