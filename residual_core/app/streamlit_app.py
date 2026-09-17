"""Residual_Assembler Streamlit GUI.

Six-tab workflow over exactly the same backend the ``resasm`` CLI uses. Every
button in this app builds an argv list and calls ``residual_core.ui.cli.main``
in-process, then prints the captured stdout/stderr and the real exit code. The
GUI therefore cannot drift from the CLI, and it cannot show you a result the
CLI would not show you:

  1. Model         - pick / upload a model, ``inspect`` it, list ``modes``
  2. Requirements  - ``requirements --mode``, ``doctor``, config template
  3. Assemble      - ``assemble`` and ``verify`` the global residual R
  4. Sensitivity   - ``sensitivity`` (T U^(p) = -R^(p)) with the FD cross-check
  5. Job           - the config workflow: ``init`` / ``init-assembly`` /
                     ``check`` / ``run`` / ``report``
  6. Backends      - ``backends`` registry audit and ``template`` contracts

Each panel shows the command it is about to run, so the GUI teaches the CLI
rather than hiding it. Exit codes are reported verbatim:

  0  the command succeeded
  1  the command ran and the answer is "not ready" / "checks failed"
  2  the command could not run with the inputs given (missing ingredient)
  3  ``sensitivity`` only: OTILib was requested and is not installed

Start it with::

    pip install -e ".[gui]"
    streamlit run scripts/app.py
"""

from __future__ import annotations

import contextlib
import io
import os
import time
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

import streamlit as st


# ----------------------------------------------------------------------------
# Path discovery
# ----------------------------------------------------------------------------

def _repo_root() -> Path:
    """The repository checkout: .../residual_core/app/streamlit_app.py -> ..."""
    return Path(__file__).resolve().parents[2]


REPO_ROOT = _repo_root()
EXAMPLES_DIR = REPO_ROOT / "residual_core" / "examples"
TEMPLATES_DIR = REPO_ROOT / "templates"
DEFAULT_WORKDIR = REPO_ROOT / "resasm_gui_workspace"


# ----------------------------------------------------------------------------
# CLI bridge
# ----------------------------------------------------------------------------

@dataclass
class CliResult:
    """One invocation of ``resasm``: what was asked, and everything it said."""

    argv: list[str]
    code: int
    stdout: str
    stderr: str
    seconds: float
    cwd: str

    @property
    def command(self) -> str:
        return "resasm " + " ".join(_quote(a) for a in self.argv)


def _quote(token: str) -> str:
    return f'"{token}"' if (" " in token or not token) else token


def _run(argv: Sequence[str], cwd: str | os.PathLike[str] | None = None) -> CliResult:
    """Call the real CLI entry point in-process and capture everything.

    ``cwd`` matters: ``init``, ``init-assembly``, ``check`` and ``run`` resolve
    relative paths against the working directory, exactly as they would in a
    terminal. We chdir around the call and always restore, so a failing command
    cannot leave the app pointing somewhere else.
    """
    from residual_core.ui import cli

    argv = [str(a) for a in argv]
    out, err = io.StringIO(), io.StringIO()
    previous = os.getcwd()
    target = str(cwd or previous)
    code = 0
    started = time.time()
    try:
        os.chdir(target)
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            code = int(cli.main(argv) or 0)
    except SystemExit as exc:  # argparse usage errors call sys.exit()
        code = int(exc.code or 0)
    except BaseException:  # noqa: BLE001 - the traceback is the useful output
        err.write(traceback.format_exc())
        code = 1
    finally:
        os.chdir(previous)
    return CliResult(argv=argv, code=code, stdout=out.getvalue(),
                     stderr=err.getvalue(), seconds=time.time() - started,
                     cwd=target)


_EXIT_MEANING = {
    0: "success",
    1: "ran, and the answer is negative (not ready / checks failed)",
    2: "could not run with these inputs (missing ingredient or usage error)",
    3: "OTILib was requested and is not installed (no silent fallback)",
}


# ----------------------------------------------------------------------------
# Session state
# ----------------------------------------------------------------------------

_DEFAULTS: dict[str, Any] = {
    "model_path": "",
    "model_label": "",
    "fields_path": "",
    "subroutine_path": "",
    "workdir": str(DEFAULT_WORKDIR),
    "results": {},
}


def _init_state() -> None:
    for key, value in _DEFAULTS.items():
        st.session_state.setdefault(key, value.copy() if isinstance(value, dict) else value)


def _remember(key: str, result: CliResult) -> None:
    st.session_state.results[key] = result


def _recall(key: str) -> CliResult | None:
    return st.session_state.results.get(key)


# ----------------------------------------------------------------------------
# Rendering helpers
# ----------------------------------------------------------------------------

def _preview(argv: Sequence[str]) -> None:
    """Show the command this button will run, before it runs."""
    st.code("resasm " + " ".join(_quote(str(a)) for a in argv), language="bash")


def _render(key: str) -> None:
    """Render the stored result for ``key``, if the user has run it."""
    res = _recall(key)
    if res is None:
        return
    meaning = _EXIT_MEANING.get(res.code, "unrecognised exit code")
    line = f"exit code {res.code} - {meaning}  ({res.seconds:.2f} s)"
    if res.code == 0:
        st.success(line)
    elif res.code == 1:
        st.warning(line)
    else:
        st.error(line)
    st.caption(f"$ {res.command}    (in {res.cwd})")
    if res.stdout.strip():
        st.code(res.stdout, language="text")
    if res.stderr.strip():
        st.markdown("**stderr**")
        st.code(res.stderr, language="text")
    if not res.stdout.strip() and not res.stderr.strip():
        st.caption("(the command printed nothing)")


def _action(label: str, key: str, argv: Sequence[str], *, cwd=None,
            disabled: bool = False, spinner: str | None = None,
            help: str | None = None) -> None:
    """A preview + button + captured-output block for one CLI invocation."""
    _preview(argv)
    if st.button(label, key=f"btn_{key}", disabled=disabled, help=help):
        with st.spinner(spinner or f"running {label.lower()} ..."):
            _remember(key, _run(argv, cwd=cwd))
    _render(key)


def _need_model() -> bool:
    if not st.session_state.model_path:
        st.info("Load a model on tab **1. Model** first.")
        return False
    return True


def _optional_args(flag: str, value: str) -> list[str]:
    return [flag, value] if value else []


# ----------------------------------------------------------------------------
# Discovery of shipped examples
# ----------------------------------------------------------------------------

def _shipped_examples() -> dict[str, Path]:
    """Neutral-format models shipped in the repository, newest name order."""
    found: dict[str, Path] = {}
    for model in sorted(EXAMPLES_DIR.glob("*/model.json")):
        found[f"examples/{model.parent.name}"] = model
    return found


def _sibling(model: str, name: str) -> str:
    """A file next to the model, if it exists (e.g. fields.json, params.json)."""
    if not model:
        return ""
    candidate = Path(model).resolve().parent / name
    return str(candidate) if candidate.exists() else ""


def _templates() -> list[str]:
    from residual_core.ui.cli import _TEMPLATES

    return sorted(_TEMPLATES)


def _modes() -> list[str]:
    from residual_core.core import requirements as req

    return list(req.known_modes())


def _formulations() -> list[str]:
    from residual_core.formulations.registry import build_formulation_registry

    return list(build_formulation_registry().names())


def _materials() -> list[str]:
    from residual_core.materials.registry import build_material_registry

    return list(build_material_registry().names())


def _otilib_status() -> dict[str, Any]:
    try:
        from residual_core.algebra.otilib_adapter import otilib_status

        return otilib_status()
    except Exception as exc:  # noqa: BLE001 - reported, never raised at the user
        return {"available": False, "error": f"{type(exc).__name__}: {exc}"}


# ----------------------------------------------------------------------------
# 1. Model
# ----------------------------------------------------------------------------

def _tab_model() -> None:
    st.header("1. Model")
    st.markdown(
        "A *model* is a mesh plus whatever ingredients came with it: an Abaqus "
        "`.inp` deck, or a neutral `.json` model. Everything downstream needs "
        "one. Pick a shipped example, type a path, or upload a file."
    )

    examples = _shipped_examples()
    left, right = st.columns(2)

    with left:
        st.subheader("Shipped examples")
        if examples:
            labels = ["-"] + list(examples)
            chosen = st.selectbox("Example model", labels, key="example_pick")
            if chosen != "-" and st.button("Use this example", key="btn_use_example"):
                st.session_state.model_path = str(examples[chosen])
                st.session_state.model_label = chosen
                st.session_state.fields_path = _sibling(str(examples[chosen]), "fields.json")
        else:
            st.warning(
                f"No `model.json` found under {EXAMPLES_DIR}. Generate them with "
                "`python -m residual_core.examples.generate_minimal`."
            )

        st.subheader("Path on this machine")
        typed = st.text_input("Model file (.inp or .json)", value="",
                              key="typed_model",
                              placeholder="/absolute/path/to/model.inp")
        if st.button("Use this path", key="btn_use_path", disabled=not typed):
            resolved = Path(typed).expanduser()
            if resolved.exists():
                st.session_state.model_path = str(resolved.resolve())
                st.session_state.model_label = resolved.name
                st.session_state.fields_path = _sibling(str(resolved.resolve()), "fields.json")
            else:
                st.error(f"No such file: {resolved}")

    with right:
        st.subheader("Upload")
        st.caption(
            "Uploaded files are written into the working directory shown in the "
            "sidebar so that the CLI can open them by path. An `.inp` deck that "
            "`*INCLUDE`s other files will only work if you upload it together "
            "with them, or point at the original folder instead."
        )
        upload = st.file_uploader("Model file", type=["inp", "json"], key="model_upload")
        if upload is not None and st.button("Use uploaded file", key="btn_use_upload"):
            workdir = Path(st.session_state.workdir)
            workdir.mkdir(parents=True, exist_ok=True)
            target = workdir / upload.name
            target.write_bytes(upload.getbuffer())
            st.session_state.model_path = str(target)
            st.session_state.model_label = upload.name
            st.session_state.fields_path = _sibling(str(target), "fields.json")
            st.success(f"saved to {target}")

    st.markdown("---")
    if not st.session_state.model_path:
        st.info("No model loaded yet.")
        return

    st.markdown(f"**Loaded model:** `{st.session_state.model_path}`")
    if st.session_state.fields_path:
        st.caption(f"field export found next to it: {st.session_state.fields_path}")

    st.subheader("Inspect")
    st.markdown(
        "`inspect` reports the element types it recognised, the materials it "
        "found, which ingredients are still missing, and which assembly modes "
        "are therefore possible. It never guesses a material constant."
    )
    detail = st.checkbox("--detail (per-element backend selection)", key="inspect_detail")
    argv = ["inspect", st.session_state.model_path] + (["--detail"] if detail else [])
    _action("Run inspect", "inspect", argv)

    st.markdown("---")
    st.subheader("Assembly modes")
    st.markdown(
        "The four modes differ only in *where the stress comes from*: "
        "`stress-driven` reads an exported field, `material-replay` re-runs a "
        "UMAT, `formulation` evaluates a self-contained element, and "
        "`direct-residual` takes R from a UEL."
    )
    _action("List modes", "modes", ["modes"])


# ----------------------------------------------------------------------------
# 2. Requirements
# ----------------------------------------------------------------------------

def _tab_requirements() -> None:
    st.header("2. Requirements & doctor")
    st.markdown(
        "Before assembling anything, ask the framework what a given mode needs "
        "and what it already has. This is the panel that answers *why can't I "
        "run yet* - with a named missing ingredient rather than a stack trace."
    )
    if not _need_model():
        return
    model = st.session_state.model_path

    st.subheader("requirements --mode")
    modes = _modes()
    mode = st.selectbox("Mode", modes,
                        index=modes.index("stress-driven") if "stress-driven" in modes else 0,
                        key="req_mode")
    fields = st.text_input("--fields (exported element field, JSON)",
                           value=st.session_state.fields_path, key="req_fields")
    subroutine = st.text_input("--subroutine (UMAT/UEL source, for material-replay)",
                               value=st.session_state.subroutine_path, key="req_sub")
    argv = (["requirements", model, "--mode", mode]
            + _optional_args("--fields", fields)
            + _optional_args("--subroutine", subroutine))
    _action("Run requirements", "requirements", argv)
    st.caption(
        "Exit code 1 here is not a bug: it is the tool saying *this mode is not "
        "runnable yet*, and the report above names what to supply."
    )

    st.markdown("---")
    st.subheader("doctor")
    st.markdown(
        "`doctor` is `inspect` plus a per-mode readiness table, and it can write "
        "a pre-filled config template you then edit by hand."
    )
    write_template = st.checkbox("--write-config-template", key="doctor_write")
    template_path = st.text_input(
        "Template path", value=str(Path(st.session_state.workdir) / "resasm_config.yml"),
        key="doctor_template", disabled=not write_template)
    argv = ["doctor", model]
    if write_template and template_path:
        argv += ["--write-config-template", template_path]
    _action("Run doctor", "doctor", argv, cwd=st.session_state.workdir)


# ----------------------------------------------------------------------------
# 3. Assemble & verify
# ----------------------------------------------------------------------------

def _tab_assemble() -> None:
    st.header("3. Assemble & verify R")
    st.markdown(
        "Abaqus never hands you the global residual, so this assembles it from "
        "ingredients:\n\n"
        "$$R(u,a) = F_\\mathrm{internal} - F_\\mathrm{external} + F_\\mathrm{constraints}$$\n\n"
        "The printed `||R||` is the norm over all degrees of freedom; `max|R|` "
        "is the worst single DOF."
    )
    if not _need_model():
        return
    model = st.session_state.model_path

    st.subheader("assemble")
    modes = _modes()
    mode = st.selectbox("--mode", modes,
                        index=modes.index("stress-driven") if "stress-driven" in modes else 0,
                        key="asm_mode")
    col_a, col_b = st.columns(2)
    with col_a:
        fields = st.text_input("--fields (JSON field export)",
                               value=st.session_state.fields_path, key="asm_fields")
        tangent = st.checkbox("--tangent (also assemble T)", key="asm_tangent")
    with col_b:
        subroutine = st.text_input("--subroutine (UMAT/UEL source)",
                                   value=st.session_state.subroutine_path, key="asm_sub")
        out = st.text_input("--out (save R to .npy)", value="", key="asm_out",
                            placeholder="R.npy")

    argv = (["assemble", model, "--mode", mode]
            + _optional_args("--fields", fields)
            + _optional_args("--subroutine", subroutine)
            + (["--tangent"] if tangent else [])
            + _optional_args("--out", out))
    _action("Run assemble", "assemble", argv, cwd=st.session_state.workdir)
    st.caption(
        "If the mode is not runnable, `assemble` returns 2 and prints the same "
        "requirements report as tab 2 rather than assembling something wrong."
    )

    st.markdown("---")
    st.subheader("verify")
    st.markdown(
        "`verify` is the sanity check for a stress-driven assembly: on a "
        "converged solution the residual on the *free* DOFs should be small, "
        "and what is left should sit on the constrained DOFs as reactions. A "
        "large `||R_free||` means the ingredients do not describe the same "
        "state - not that the assembler is broken."
    )
    vfields = st.text_input("--fields (JSON field export)",
                            value=st.session_state.fields_path, key="ver_fields")
    _action("Run verify", "verify",
            ["verify", model] + _optional_args("--fields", vfields),
            cwd=st.session_state.workdir)


# ----------------------------------------------------------------------------
# 4. Sensitivity
# ----------------------------------------------------------------------------

def _tab_sensitivity() -> None:
    st.header("4. Sensitivity")
    st.markdown(
        "Differentiating the residual with respect to design parameters gives "
        "the sensitivity system\n\n"
        "$$T\\,U^{(p)} = -R^{(p)}$$\n\n"
        "solved order by order. `--backend otilib` is the production path; "
        "`--backend dual1` is a legacy order-1 backend that is useful as a "
        "smoke test when OTILib is not installed."
    )
    status = _otilib_status()
    if status.get("available"):
        st.success(f"OTILib is available (api module: {status.get('api_module')})")
    else:
        st.warning(
            "OTILib is **not** installed in this environment, so "
            "`--backend otilib` will stop with exit code 3 instead of quietly "
            "falling back to a less accurate method. `--backend dual1` still "
            "works for order 1."
        )
        if status.get("error"):
            st.code(str(status["error"]), language="text")

    if not _need_model():
        return
    model = st.session_state.model_path

    st.subheader("sensitivity")
    params_file = st.text_input(
        "--params (params file: parameters, mode, order, expected)",
        value=_sibling(model, "params.json"), key="sens_params")
    raw_params = st.text_area(
        "--param (one 'group.key' per line; repeatable flag)",
        value="", key="sens_param_list",
        placeholder="material.E\nspring.k")
    col_a, col_b, col_c = st.columns(3)
    with col_a:
        modes = ["(from params file / formulation)"] + _modes()
        mode = st.selectbox("--mode", modes, key="sens_mode")
    with col_b:
        order = st.number_input("--order", min_value=0, max_value=6, value=0,
                                key="sens_order",
                                help="0 means: do not pass --order, let the "
                                     "params file or the default decide")
    with col_c:
        backend = st.selectbox("--backend", ["(default: otilib)", "otilib", "dual1"],
                               key="sens_backend")
    no_fd = st.checkbox("--no-fd (skip the finite-difference cross-check)",
                        key="sens_nofd")
    out = st.text_input("--out (save PREFIX.{json,npz,md})", value="",
                        key="sens_out", placeholder="sens_run")

    argv = ["sensitivity", model]
    argv += _optional_args("--params", params_file)
    for name in [p.strip() for p in raw_params.splitlines() if p.strip()]:
        argv += ["--param", name]
    if not mode.startswith("("):
        argv += ["--mode", mode]
    if int(order) > 0:
        argv += ["--order", str(int(order))]
    if not backend.startswith("("):
        argv += ["--backend", backend]
    if no_fd:
        argv += ["--no-fd"]
    argv += _optional_args("--out", out)

    _action("Run sensitivity", "sensitivity", argv, cwd=st.session_state.workdir,
            spinner="solving T U^(p) = -R^(p) ...")
    st.caption(
        "Keep the finite-difference cross-check on while you are learning: it "
        "is the only thing in the report that was computed by a completely "
        "different method."
    )


# ----------------------------------------------------------------------------
# 5. Job (the config workflow)
# ----------------------------------------------------------------------------

def _tab_job() -> None:
    st.header("5. Job: init -> check -> run -> report")
    workdir = st.session_state.workdir
    st.markdown(
        "The four-step user workflow. A *job* is a folder with a residual (or a "
        "model), a small `resasm.yml`, and an output package. All paths on this "
        f"tab are resolved relative to the working directory: `{workdir}`."
    )
    Path(workdir).mkdir(parents=True, exist_ok=True)

    st.subheader("a. Start from a template")
    st.markdown(
        "Each template is a complete, runnable job. `python` implements R in "
        "Python and differentiates it with OTILib; the `blackbox` variants get "
        "derivative coefficients from an executable you own, which is the route "
        "that works without OTILib installed here; `cpp` and `fortran` are the "
        "same idea in those languages."
    )
    col_a, col_b = st.columns(2)
    with col_a:
        template = st.selectbox("--template", _templates(), key="job_template")
    with col_b:
        out_dir = st.text_input("--out (job folder)", value="my_job", key="job_out")
    force = st.checkbox("--force (overwrite an existing folder)", key="job_force")
    argv = ["init", "--template", template, "--out", out_dir] + (["--force"] if force else [])
    _action("Copy template", "job_init", argv, cwd=workdir)
    st.caption(
        "`resasm init` without `--template` starts an interactive wizard that "
        "reads answers from the terminal. A web page has no terminal, so the "
        "wizard is CLI-only; use the template here and edit its `resasm.yml`."
    )

    st.markdown("---")
    st.subheader("b. Or write a recipe from a model (assembly path)")
    st.markdown(
        "`inspect-model` says what is inferable and what is missing; "
        "`init-assembly` writes the `resasm.yml`. **Pass an absolute model "
        "path**: the recipe records `mesh:` exactly as you typed it, so a "
        "relative path written into a folder elsewhere will not resolve when "
        "`check` reads it back."
    )
    model = st.session_state.model_path
    if not model:
        st.info("Load a model on tab 1 to use this section.")
    else:
        solution = st.text_input("--solution (converged U.npy)", value="", key="job_sol")
        material = st.text_input("--material (material evaluator, e.g. umat.f)",
                                 value="", key="job_mat")
        raw = st.text_area("--param (one 'group.key' per line)", value="", key="job_params")
        name = st.text_input("--name (job name)", value="assembly_job", key="job_name")
        recipe_out = st.text_input("--out (config path)", value="resasm.yml",
                                   key="job_recipe_out")
        params = [p.strip() for p in raw.splitlines() if p.strip()]

        common = (_optional_args("--solution", solution)
                  + _optional_args("--material", material))
        inspect_argv = ["inspect-model", model] + common
        for p in params:
            inspect_argv += ["--param", p]
        _action("Run inspect-model", "job_inspect_model", inspect_argv, cwd=workdir)

        assembly_argv = ["init-assembly", "--model", model] + common
        for p in params:
            assembly_argv += ["--param", p]
        assembly_argv += ["--name", name, "--out", recipe_out]
        _action("Write recipe", "job_init_assembly", assembly_argv, cwd=workdir)

    st.markdown("---")
    st.subheader("c. Check the config")
    st.markdown(
        "`check` reads the config and reports one `[ok]` or `[fail]` line per "
        "requirement. Fix every `[fail]` before running: a failing check is a "
        "refusal to produce a number, and that is the point."
    )
    config_path = st.text_input("Config file", value="my_job/resasm.yml",
                                key="job_config")
    _action("Run check", "job_check", ["check", config_path], cwd=workdir)

    st.markdown("---")
    st.subheader("d. Run the job")
    st.markdown(
        "`run` produces two output folders. `private/` holds the full arrays; "
        "`public/` holds only norms, rankings, status and timing, so it can be "
        "shared when the model itself cannot be."
    )
    _action("Run job", "job_run", ["run", config_path], cwd=workdir,
            spinner="running the sensitivity job ...")

    st.markdown("---")
    st.subheader("e. Read the report")
    report_dir = st.text_input("Output directory", value="my_job/out", key="job_report")
    _action("Read report", "job_report", ["report", report_dir], cwd=workdir)

    summary = Path(workdir) / report_dir / "public" / "summary.md"
    if summary.exists():
        with st.expander(f"public/summary.md ({summary})"):
            st.markdown(summary.read_text(encoding="utf-8"))


# ----------------------------------------------------------------------------
# 6. Backends & contracts
# ----------------------------------------------------------------------------

def _tab_backends() -> None:
    st.header("6. Backends & contracts")
    st.markdown(
        "What is actually registered in this checkout, and what each backend "
        "declares about itself: which element types it takes, which modes it "
        "supports, whether it can produce a tangent, how far it has been "
        "verified, and - importantly - its stated limitations."
    )
    _action("Audit all backends", "backends", ["backends"])

    st.markdown("---")
    st.subheader("template: one backend's declared contract")
    st.markdown(
        "Print the exact inputs a backend expects. Read this before writing a "
        "config against it; the contract is the specification."
    )
    kind = st.radio("Contract kind", ["--formulation", "--material"],
                    horizontal=True, key="tpl_kind")
    names = _formulations() if kind == "--formulation" else _materials()
    if not names:
        st.warning("The registry reported no entries of that kind.")
        return
    name = st.selectbox("Name", names, key="tpl_name")
    _action("Print contract", "template", ["template", kind, name])


# ----------------------------------------------------------------------------
# Sidebar & entry point
# ----------------------------------------------------------------------------

def _sidebar() -> None:
    with st.sidebar:
        st.title("Residual_Assembler")
        st.caption("Assemble R, then differentiate it")
        st.markdown("---")

        if st.session_state.model_path:
            st.markdown(f"**Model:** `{Path(st.session_state.model_path).name}`")
            st.caption(st.session_state.model_path)
        else:
            st.markdown(":grey[No model loaded]")

        st.markdown("---")
        st.markdown("**Working directory**")
        st.caption(
            "Where relative paths resolve and where generated jobs are written, "
            "exactly as if you had `cd`-ed there in a terminal."
        )
        workdir = st.text_input("Path", value=st.session_state.workdir,
                                key="sidebar_workdir", label_visibility="collapsed")
        if workdir != st.session_state.workdir:
            st.session_state.workdir = workdir
        if st.button("Create it", key="btn_mkdir"):
            Path(st.session_state.workdir).mkdir(parents=True, exist_ok=True)
            st.success("created")

        st.markdown("---")
        status = _otilib_status()
        if status.get("available"):
            st.success("OTILib: available")
        else:
            st.error("OTILib: not installed")
            st.caption(
                "OTILib is an external GPLv3 package and is not vendored here. "
                "Without it, order >= 1 derivatives are available through the "
                "black-box path and through `--backend dual1`; the Python "
                "OTILib path will refuse to run rather than approximate."
            )
        st.caption(f"repo: {status.get('repo', 'n/a')}")

        st.markdown("---")
        st.caption(f"Repo root: {REPO_ROOT}")
        st.caption("Every button runs the real `resasm` CLI in-process.")


def main() -> None:
    st.set_page_config(page_title="Residual_Assembler", layout="wide")
    _init_state()
    _sidebar()

    tabs = st.tabs([
        "1. Model",
        "2. Requirements",
        "3. Assemble",
        "4. Sensitivity",
        "5. Job",
        "6. Backends",
    ])
    with tabs[0]:
        _tab_model()
    with tabs[1]:
        _tab_requirements()
    with tabs[2]:
        _tab_assemble()
    with tabs[3]:
        _tab_sensitivity()
    with tabs[4]:
        _tab_job()
    with tabs[5]:
        _tab_backends()


if __name__ == "__main__":
    main()
