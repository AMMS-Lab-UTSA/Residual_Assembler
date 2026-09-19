"""Residual_Assembler Streamlit GUI.

A guided workflow over exactly the same backend the ``resasm`` CLI uses. Every
button in this app builds an argv list and calls ``residual_core.ui.cli.main``
in-process, then prints the captured stdout/stderr and the real exit code. The
GUI therefore cannot drift from the CLI, and it cannot show you a result the
CLI would not show you:

  0. Start here    - what this tool is, and a demo that runs end to end
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
import re
import time
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Sequence

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

#: The example the "Start here" tab drives, and the one every locked panel
#: offers to load for you. It is the only shipped example that ships with the
#: exported stress field a stress-driven assembly needs.
DEMO_EXAMPLE = "examples/minimal_c3d8_stress_driven"

_NONE = "- select -"


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

#: What to *do* about each exit code. The meaning above is what happened; this
#: is the next move, which is the part people actually get stuck on.
_EXIT_ADVICE = {
    1: "This is an answer, not a crash. The report above names the check that "
       "failed or the ingredient that is missing - fix that and run again.",
    2: "The command refused to start. Usually a path that does not exist, a "
       "missing `--fields` / `--subroutine`, or a flag typo. The text above "
       "says which.",
    3: "OTILib is an external GPLv3 package and is not vendored here. Either "
       "build it (`bash scripts/setup_otilib.sh`), or take a route that does "
       "not need it: `--backend dual1` for order 1, or the `blackbox` job "
       "templates for any order.",
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


def _succeeded(key: str) -> bool:
    res = _recall(key)
    return res is not None and res.code == 0


# ----------------------------------------------------------------------------
# Rendering helpers
# ----------------------------------------------------------------------------

def _preview(argv: Sequence[str]) -> None:
    """Show the command this button will run, before it runs."""
    st.code("resasm " + " ".join(_quote(str(a)) for a in argv), language="bash")


def _render(key: str, next_hint: str | None = None) -> None:
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
        _long_text(res.stdout, "output")
    if res.stderr.strip():
        st.markdown("**stderr**")
        _long_text(res.stderr, "stderr")
    if not res.stdout.strip() and not res.stderr.strip():
        st.caption("(the command printed nothing)")

    advice = _EXIT_ADVICE.get(res.code)
    if advice:
        st.info(advice)
    elif res.code == 0 and next_hint:
        st.caption(f"Next: {next_hint}")


def _long_text(text: str, what: str, head: int = 18) -> None:
    """Print short output inline; fold long output so the page stays readable."""
    lines = text.rstrip("\n").splitlines()
    if len(lines) <= head + 4:
        st.code(text, language="text")
        return
    st.code("\n".join(lines[:head]) + "\n...", language="text")
    with st.expander(f"full {what} ({len(lines)} lines)"):
        st.code(text, language="text")


def _action(label: str, key: str, argv: Sequence[str], *, cwd=None,
            disabled: bool = False, spinner: str | None = None,
            help: str | None = None, next_hint: str | None = None,
            primary: bool = False) -> None:
    """A preview + button + captured-output block for one CLI invocation."""
    _preview(argv)
    if st.button(label, key=f"btn_{key}", disabled=disabled, help=help,
                 type="primary" if primary else "secondary"):
        with st.spinner(spinner or f"running {label.lower()} ..."):
            _remember(key, _run(argv, cwd=cwd))
    _render(key, next_hint=next_hint)


def _reported_public_dir(result: CliResult | None) -> str | None:
    """The public folder a successful ``resasm report`` named in its output."""
    if result is None or result.code != 0:
        return None
    for line in result.stdout.splitlines():
        if line.strip().startswith("public  outputs"):
            return line.split(":", 1)[1].strip()
    return None


def _download(path: Path, label: str, key: str) -> None:
    """Offer a produced artefact for download, only once it actually exists."""
    if not path.is_file():
        return
    size = path.stat().st_size
    st.download_button(f"Download {label} ({size:,} bytes)",
                       data=path.read_bytes(), file_name=path.name,
                       key=f"dl_{key}")


def _need_model(key: str, what: str) -> bool:
    """Guard a panel that cannot do anything without a model - and unblock it.

    A locked panel that only says "go somewhere else" wastes a click. This one
    loads the demo model for you, so every tab is reachable from every tab.
    """
    if st.session_state.model_path:
        return True
    st.warning(f"**Locked.** {what} needs a model.")
    left, right = st.columns([1, 2])
    with left:
        if st.button("Load the demo model", key=f"unlock_{key}", type="primary"):
            _load_example(DEMO_EXAMPLE)
            st.rerun()
    with right:
        st.caption(
            f"Loads `{DEMO_EXAMPLE}` - a single unit-cube C3D8 solid that ships "
            "with its stress field, so every panel here becomes runnable. Tab "
            "**1. Model** has the other examples, your own path, and upload."
        )
    return False


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


def _blurb(readme: Path) -> tuple[str, str]:
    """The title and opening paragraph of an example's README, de-marked-down.

    The examples already document themselves; a dropdown that shows only the
    folder name throws that away and makes the user open five files to find
    out which example does what.
    """
    if not readme.is_file():
        return "", ""
    lines = readme.read_text(encoding="utf-8", errors="replace").splitlines()
    title, para = "", []
    for index, line in enumerate(lines):
        if not line.startswith("# "):
            continue
        title = line[2:].strip()
        for rest in lines[index + 1:]:
            if rest.strip():
                para.append(rest.strip())
            elif para:
                break
        break
    text = re.sub(r"[*`$\\]", "", " ".join(para))
    return title, re.sub(r"\s+", " ", text).strip()


def _catalogue() -> dict[str, dict[str, Any]]:
    """Every shipped example with enough context to choose between them."""
    out: dict[str, dict[str, Any]] = {}
    for label, model in _shipped_examples().items():
        folder = model.parent
        title, text = _blurb(folder / "README.md")
        extras = sorted(p.name for p in folder.iterdir()
                        if p.is_file() and p.name not in {"model.json", "README.md"})
        out[label] = {"path": model, "title": title or folder.name,
                      "blurb": text, "extras": extras, "folder": folder}
    return out


def _load_example(label: str) -> None:
    """Point every tab at a shipped example, including its side files."""
    model = _shipped_examples().get(label)
    if model is None:
        return
    st.session_state.model_path = str(model)
    st.session_state.model_label = label
    st.session_state.fields_path = _sibling(str(model), "fields.json")


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


def _mode_index(modes: Sequence[str], preferred: str = "stress-driven") -> int:
    return modes.index(preferred) if preferred in modes else 0


# ----------------------------------------------------------------------------
# Progress: which steps of the workflow have actually succeeded
# ----------------------------------------------------------------------------

_STEPS: list[tuple[str, Callable[[], bool]]] = [
    ("Load a model", lambda: bool(st.session_state.model_path)),
    ("Inspect it", lambda: _succeeded("inspect") or _succeeded("demo_inspect")),
    ("Assemble R", lambda: _succeeded("assemble") or _succeeded("demo_assemble")),
    ("Verify R", lambda: _succeeded("verify") or _succeeded("demo_verify")),
    ("Differentiate", lambda: _succeeded("sensitivity")),
    ("Package a job", lambda: _succeeded("job_run")),
]


def _progress() -> list[tuple[bool, str]]:
    return [(bool(done()), label) for label, done in _STEPS]


# ----------------------------------------------------------------------------
# 0. Start here
# ----------------------------------------------------------------------------

def _tab_start() -> None:
    st.header("Start here")
    st.markdown(
        "**Residual_Assembler assembles the thing your solver never hands you.**"
        " A commercial FE solver reports displacements and stresses, but not the"
        " global residual\n\n"
        "$$R(u,a) = F_\\mathrm{internal} - F_\\mathrm{external} + F_\\mathrm{constraints}$$\n\n"
        "and without $R$ you cannot differentiate the *model* with respect to a "
        "design parameter $a$. This tool builds $R$ from whatever ingredients "
        "you actually have, checks it, and then solves the sensitivity system "
        "$T\\,U^{(p)} = -R^{(p)}$ order by order."
    )

    col_a, col_b, col_c = st.columns(3)
    with col_a:
        st.markdown("**What you give it**")
        st.caption(
            "A mesh (Abaqus `.inp` or neutral `.json`), plus *one* source of "
            "stress: an exported field, a UMAT, a self-contained element, or a "
            "UEL residual."
        )
    with col_b:
        st.markdown("**What you get back**")
        st.caption(
            "`||R||` and `max|R|`, the residual vector as `.npy`, a reaction "
            "check, and - with parameters - derivative coefficients with a "
            "finite-difference cross-check."
        )
    with col_c:
        st.markdown("**What it refuses to do**")
        st.caption(
            "Guess a material constant, approximate a derivative you asked to "
            "be exact, or assemble a mode whose ingredients are missing. Those "
            "are exit codes 1, 3 and 2 - never a wrong number."
        )

    st.markdown("---")
    st.subheader("Run the whole thing once, right now")
    catalogue = _catalogue()
    demo = catalogue.get(DEMO_EXAMPLE)
    if demo is None:
        st.error(
            f"The demo example is missing from {EXAMPLES_DIR}. Regenerate it "
            "with `python -m residual_core.examples.generate_minimal`."
        )
        return

    st.markdown(
        f"This runs four real CLI commands against `{DEMO_EXAMPLE}` - "
        "*inspect*, *requirements*, *assemble*, *verify* - in order, and shows "
        "you each command and its verbatim output. Nothing is mocked and "
        "nothing is cached; it takes a couple of seconds and needs no Abaqus "
        "and no OTILib."
    )
    st.caption(f"{demo['title']} — {demo['blurb']}")

    workdir = Path(st.session_state.workdir)
    if st.button("Run the demo", key="btn_demo", type="primary"):
        workdir.mkdir(parents=True, exist_ok=True)
        _load_example(DEMO_EXAMPLE)
        model = st.session_state.model_path
        fields = st.session_state.fields_path
        steps = [
            ("demo_inspect", ["inspect", model], None),
            ("demo_requirements",
             ["requirements", model, "--mode", "stress-driven", "--fields", fields],
             None),
            ("demo_assemble",
             ["assemble", model, "--mode", "stress-driven", "--fields", fields,
              "--out", "demo_R.npy"], str(workdir)),
            ("demo_verify", ["verify", model, "--fields", fields], str(workdir)),
        ]
        with st.spinner("running inspect -> requirements -> assemble -> verify ..."):
            for key, argv, cwd in steps:
                _remember(key, _run(argv, cwd=cwd))

    titles = {
        "demo_inspect": "1. inspect - what is in this model?",
        "demo_requirements": "2. requirements - can stress-driven run?",
        "demo_assemble": "3. assemble - build R and save it",
        "demo_verify": "4. verify - is R small where it should be?",
    }
    ran = [k for k in titles if _recall(k) is not None]
    if ran:
        st.markdown("---")
        for key, title in titles.items():
            if _recall(key) is None:
                continue
            st.markdown(f"##### {title}")
            _render(key)
        _download(workdir / "demo_R.npy", "demo_R.npy (the residual vector)",
                  "demo_R")
        if all(_succeeded(k) for k in titles):
            st.success(
                "All four succeeded, and the model is now loaded, so every "
                "other tab is unlocked. Tab **3. Assemble** is the same "
                "assembly with the knobs exposed; tab **4. Sensitivity** is "
                "the part that needs parameters."
            )

    st.markdown("---")
    st.subheader("What the exit codes mean")
    st.markdown(
        "| code | meaning | what to do |\n"
        "|---|---|---|\n"
        "| `0` | success | read the numbers |\n"
        "| `1` | it ran, and the answer is negative | the report names the "
        "failing check or missing ingredient |\n"
        "| `2` | it could not run with these inputs | fix the path or the "
        "missing flag |\n"
        "| `3` | OTILib was requested and is not installed | build it, or use "
        "`--backend dual1` / a `blackbox` template |\n"
    )
    st.caption(
        "Every panel in this app reports the code it really got. A yellow box "
        "is the tool telling you something true, not the GUI failing."
    )


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

    catalogue = _catalogue()
    pick, own = st.tabs(["Shipped examples", "Your own model"])

    with pick:
        if not catalogue:
            st.warning(
                f"No `model.json` found under {EXAMPLES_DIR}. Generate them with "
                "`python -m residual_core.examples.generate_minimal`."
            )
        else:
            labels = [_NONE] + list(catalogue)
            index = labels.index(st.session_state.model_label) \
                if st.session_state.model_label in labels else 0
            chosen = st.selectbox(
                "Example model", labels, index=index, key="example_pick",
                format_func=lambda name: name if name == _NONE
                else f"{name}  —  {catalogue[name]['title']}",
                help="Choosing an example loads it immediately; there is no "
                     "second button to press.")
            if chosen != _NONE and chosen != st.session_state.model_label:
                _load_example(chosen)
                st.rerun()

            if chosen != _NONE:
                entry = catalogue[chosen]
                st.info(entry["blurb"] or "(this example has no README summary)")
                extras = entry["extras"]
                st.caption(
                    "Ships alongside the model: " + ", ".join(f"`{e}`" for e in extras)
                    if extras else
                    "Ships with the mesh only - a stress-driven assembly would "
                    "need a `--fields` export you supply.")
                with st.expander("This example's README"):
                    readme = entry["folder"] / "README.md"
                    st.markdown(readme.read_text(encoding="utf-8")
                                if readme.is_file() else "(no README)")

    with own:
        left, right = st.columns(2)
        with left:
            st.subheader("Path on this machine")
            typed = st.text_input("Model file (.inp or .json)", value="",
                                  key="typed_model",
                                  placeholder="/absolute/path/to/model.inp")
            if st.button("Use this path", key="btn_use_path", disabled=not typed,
                         help="Give an absolute path: recipes record the path "
                              "exactly as typed."):
                resolved = Path(typed).expanduser()
                if resolved.exists():
                    st.session_state.model_path = str(resolved.resolve())
                    st.session_state.model_label = resolved.name
                    st.session_state.fields_path = _sibling(str(resolved.resolve()),
                                                            "fields.json")
                    st.rerun()
                else:
                    st.error(f"No such file: {resolved}")

        with right:
            st.subheader("Upload")
            st.caption(
                "Uploaded files are written into the working directory shown in "
                "the sidebar so that the CLI can open them by path. An `.inp` "
                "deck that `*INCLUDE`s other files will only work if you upload "
                "it together with them, or point at the original folder instead."
            )
            upload = st.file_uploader("Model file", type=["inp", "json"],
                                      key="model_upload")
            if upload is not None and st.button("Use uploaded file",
                                                key="btn_use_upload"):
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
        st.info("No model loaded yet - pick an example above to unlock the rest "
                "of the app.")
        return

    st.success(f"Loaded: `{st.session_state.model_path}`")
    if st.session_state.fields_path:
        st.caption(f"field export found next to it, and filled in for you: "
                   f"{st.session_state.fields_path}")
    else:
        st.caption("No `fields.json` next to it, so `--fields` starts empty on "
                   "tabs 2 and 3.")

    st.subheader("Inspect")
    st.markdown(
        "`inspect` reports the element types it recognised, the materials it "
        "found, which ingredients are still missing, and which assembly modes "
        "are therefore possible. It never guesses a material constant."
    )
    detail = st.checkbox("--detail (per-element backend selection)", key="inspect_detail")
    argv = ["inspect", st.session_state.model_path] + (["--detail"] if detail else [])
    _action("Run inspect", "inspect", argv, primary=True,
            next_hint="tab **2. Requirements** turns the 'missing' list into a "
                      "per-mode readiness answer.")

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
    if not _need_model("req", "Checking requirements"):
        return
    model = st.session_state.model_path

    st.subheader("requirements --mode")
    modes = _modes()
    mode = st.selectbox("Mode", modes, index=_mode_index(modes), key="req_mode")
    fields = st.text_input("--fields (exported element field, JSON)",
                           value=st.session_state.fields_path, key="req_fields",
                           help="Needed by stress-driven. Left empty for the "
                                "other modes.")
    subroutine = st.text_input("--subroutine (UMAT/UEL source, for material-replay)",
                               value=st.session_state.subroutine_path, key="req_sub")
    argv = (["requirements", model, "--mode", mode]
            + _optional_args("--fields", fields)
            + _optional_args("--subroutine", subroutine))
    _action("Run requirements", "requirements", argv, primary=True,
            next_hint="if this said ready, go to tab **3. Assemble** with the "
                      "same mode and the same flags.")
    st.caption(
        "Exit code 1 here is not a bug: it is the tool saying *this mode is not "
        "runnable yet*, and the report above names what to supply."
    )

    st.markdown("---")
    st.subheader("doctor")
    st.markdown(
        "`doctor` is `inspect` plus a per-mode readiness table, and it can write "
        "a pre-filled config template you then edit by hand. Run this one when "
        "you do not yet know *which* mode your model can support."
    )
    write_template = st.checkbox("--write-config-template", key="doctor_write")
    template_path = st.text_input(
        "Template path", value=str(Path(st.session_state.workdir) / "resasm_config.yml"),
        key="doctor_template", disabled=not write_template,
        help=None if write_template else "Tick the box above to enable this.")
    argv = ["doctor", model]
    if write_template and template_path:
        argv += ["--write-config-template", template_path]
    _action("Run doctor", "doctor", argv, cwd=st.session_state.workdir)
    if write_template and template_path:
        _download(Path(template_path), "the config template", "doctor_tpl")


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
    if not _need_model("asm", "Assembling a residual"):
        return
    model = st.session_state.model_path
    workdir = Path(st.session_state.workdir)

    st.subheader("assemble")
    modes = _modes()
    mode = st.selectbox("--mode", modes, index=_mode_index(modes), key="asm_mode")
    col_a, col_b = st.columns(2)
    with col_a:
        fields = st.text_input("--fields (JSON field export)",
                               value=st.session_state.fields_path, key="asm_fields")
        tangent = st.checkbox("--tangent (also assemble T)", key="asm_tangent",
                              help="T is the tangent the sensitivity solve "
                                   "needs; assembling it costs more time.")
    with col_b:
        subroutine = st.text_input("--subroutine (UMAT/UEL source)",
                                   value=st.session_state.subroutine_path, key="asm_sub")
        out = st.text_input("--out (save R to .npy)", value="R.npy", key="asm_out",
                            placeholder="R.npy",
                            help=f"Written inside the working directory: {workdir}")

    argv = (["assemble", model, "--mode", mode]
            + _optional_args("--fields", fields)
            + _optional_args("--subroutine", subroutine)
            + (["--tangent"] if tangent else [])
            + _optional_args("--out", out))
    _action("Run assemble", "assemble", argv, cwd=str(workdir), primary=True,
            next_hint="`verify` below is the cheap sanity check; tab "
                      "**4. Sensitivity** is what R was built for.")
    if out:
        _download(workdir / out, f"{out} (the residual vector)", "asm_out")
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
            cwd=str(workdir))


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
            "works for order 1, and the `blackbox` job templates on tab 5 work "
            "at any order."
        )
        if status.get("error"):
            with st.expander("why the adapter says it is unavailable"):
                st.code(str(status["error"]), language="text")

    if not _need_model("sens", "Differentiating a residual"):
        return
    model = st.session_state.model_path
    workdir = Path(st.session_state.workdir)

    params_file_default = _sibling(model, "params.json")
    if not params_file_default:
        st.info(
            "This model ships no `params.json`, so there is nothing to "
            "differentiate *with respect to* yet. Either name parameters by "
            "hand below, or load `examples/minimal_nonlinear_spring_sensitivity` "
            "on tab 1 - it is the example built for this tab."
        )

    st.subheader("sensitivity")
    params_file = st.text_input(
        "--params (params file: parameters, mode, order, expected)",
        value=params_file_default, key="sens_params")
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

    _action("Run sensitivity", "sensitivity", argv, cwd=str(workdir), primary=True,
            spinner="solving T U^(p) = -R^(p) ...",
            next_hint="tab **5. Job** turns this one-off command into a "
                      "reproducible, shareable job folder.")
    if out:
        for suffix, what in ((".md", "the report"), (".json", "the coefficients"),
                             (".npz", "the arrays")):
            _download(workdir / f"{out}{suffix}", f"{out}{suffix} ({what})",
                      f"sens_{suffix.strip('.')}")
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
        "model), a small `resasm.yml`, and an output package. Use this instead "
        "of tabs 3-4 when the run has to be reproducible or handed to someone "
        f"else. All paths here resolve relative to the working directory: "
        f"`{workdir}`."
    )
    Path(workdir).mkdir(parents=True, exist_ok=True)

    job = st.text_input(
        "Job folder", value="my_job", key="job_folder",
        help="One name drives the whole tab: the template is copied here, the "
             "config is read from here, and the output lands inside it.")
    job = job.strip() or "my_job"
    config_path = f"{job}/resasm.yml"
    st.caption(f"config: `{config_path}`   ·   output: where its `output: dir:` "
               f"points (default `{job}/resasm_output`)   ·   "
               f"absolute: `{Path(workdir) / job}`")

    st.markdown("---")
    st.subheader("a. Start from a template")
    st.markdown(
        "Each template is a complete, runnable job. `python` implements R in "
        "Python and differentiates it with OTILib; the `blackbox` variants get "
        "derivative coefficients from an executable you own, which is the route "
        "that works without OTILib installed here; `cpp` and `fortran` are the "
        "same idea in those languages."
    )
    col_a, col_b = st.columns([2, 1])
    with col_a:
        template = st.selectbox("--template", _templates(), key="job_template")
    with col_b:
        force = st.checkbox("--force (overwrite)", key="job_force",
                            help="Without this, `init` refuses to write into a "
                                 "folder that already exists.")
    argv = ["init", "--template", template, "--out", job] + (["--force"] if force else [])
    _action("Copy template", "job_init", argv, cwd=workdir, primary=True,
            next_hint=f"open `{Path(workdir) / job}` in an editor, adjust "
                      "`resasm.yml`, then run the check below.")
    st.caption(
        "`resasm init` without `--template` starts an interactive wizard that "
        "reads answers from the terminal. A web page has no terminal, so the "
        "wizard is CLI-only; use the template here and edit its `resasm.yml`."
    )

    live_config = Path(workdir) / config_path
    if live_config.is_file():
        with st.expander(f"{config_path} as it stands now"):
            st.code(live_config.read_text(encoding="utf-8"), language="yaml")

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
        st.info("This section needs a model; the template route above does not. "
                "Load one on tab **1. Model** to use it.")
    else:
        solution = st.text_input("--solution (converged U.npy)", value="", key="job_sol")
        material = st.text_input("--material (material evaluator, e.g. umat.f)",
                                 value="", key="job_mat")
        raw = st.text_area("--param (one 'group.key' per line)", value="", key="job_params")
        name = st.text_input("--name (job name)", value="assembly_job", key="job_name")
        recipe_out = st.text_input("--out (config path)", value=config_path,
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
    _action("Run check", "job_check", ["check", config_path], cwd=workdir,
            disabled=not live_config.is_file(),
            help=None if live_config.is_file()
            else f"No {config_path} yet - copy a template in step (a) first.",
            next_hint="all `[ok]`? run the job below.")

    st.markdown("---")
    st.subheader("d. Run the job")
    st.markdown(
        "`run` produces two output folders. `private/` holds the full arrays; "
        "`public/` holds only norms, rankings, status and timing, so it can be "
        "shared when the model itself cannot be."
    )
    _action("Run job", "job_run", ["run", config_path], cwd=workdir,
            disabled=not live_config.is_file(),
            help=None if live_config.is_file()
            else f"No {config_path} yet - copy a template in step (a) first.",
            primary=True, spinner="running the sensitivity job ...")

    st.markdown("---")
    st.subheader("e. Read the report")
    st.markdown(
        "`report` is given the job's `resasm.yml` and reads the folder that job "
        "writes to, so it follows an `output: dir:` in the file."
    )
    _action("Read report", "job_report", ["report", config_path], cwd=workdir,
            disabled=not live_config.is_file(),
            help=None if live_config.is_file()
            else f"No {config_path} yet - copy a template in step (a) first.")

    public = _reported_public_dir(_recall("job_report"))
    summary = Path(workdir) / public / "summary.md" if public else None
    if summary is not None and summary.is_file():
        with st.expander(f"public/summary.md ({summary})", expanded=True):
            st.markdown(summary.read_text(encoding="utf-8"))
        _download(summary, "public/summary.md", "job_summary")


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
        st.markdown("**Where you are**")
        done = _progress()
        for finished, label in done:
            st.markdown(f"{'✅' if finished else '⬜'} {label}")
        st.progress(sum(1 for finished, _ in done if finished) / len(done))
        st.caption("A step ticks only when the command really exited 0.")

        st.markdown("---")
        if st.session_state.model_path:
            st.markdown(f"**Model:** `{Path(st.session_state.model_path).name}`")
            st.caption(st.session_state.model_path)
            # Rendered before the tabs, so resetting the example widget here is
            # safe: the selectbox has not been instantiated yet this run.
            if st.button("Clear model", key="btn_clear_model"):
                st.session_state.model_path = ""
                st.session_state.model_label = ""
                st.session_state.fields_path = ""
                st.session_state["example_pick"] = _NONE
                st.rerun()
        else:
            st.markdown(":grey[No model loaded]")
            if st.button("Load the demo model", key="btn_sidebar_demo"):
                _load_example(DEMO_EXAMPLE)
                st.rerun()

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
        exists = Path(st.session_state.workdir).is_dir()
        st.caption("exists" if exists else ":red[does not exist yet]")
        if st.button("Create it", key="btn_mkdir", disabled=exists):
            Path(st.session_state.workdir).mkdir(parents=True, exist_ok=True)
            st.rerun()

        st.markdown("---")
        status = _otilib_status()
        if status.get("available"):
            st.success("OTILib: available")
        else:
            st.error("OTILib: not installed")
        with st.expander("What that means"):
            st.caption(
                "OTILib is an external GPLv3 package and is not vendored here. "
                "Without it, order >= 1 derivatives are available through the "
                "black-box path and through `--backend dual1`; the Python "
                "OTILib path will refuse to run rather than approximate. Build "
                "it with `bash scripts/setup_otilib.sh`. Do not "
                "`pip install pyoti` - that name is an unrelated package."
            )
            st.caption(f"repo: {status.get('repo', 'n/a')}")

        with st.expander("Exit codes"):
            for code, meaning in _EXIT_MEANING.items():
                st.caption(f"`{code}` — {meaning}")

        st.markdown("---")
        st.caption(f"Repo root: {REPO_ROOT}")
        st.caption("Every button runs the real `resasm` CLI in-process.")


def _tab_request() -> None:
    import json
    import tempfile

    from residual_core.app import request_screen

    st.header("Residual Sensitivity Solver")
    inputs = {}
    uploads = {}
    for key, label, extension in (("material", "OTI_UMAT.obj", "obj"),
                                  ("model", "Analysis.inp", "inp"),
                                  ("odb", "Analysis.odb", "odb"),
                                  ("request", "sensitivity_request.json", "json")):
        left, right = st.columns(2)
        with left:
            uploads[key] = st.file_uploader(label, type=[extension], key="request_upload_" + key)
        with right:
            inputs[key] = st.text_input(label + " path", key="request_path_" + key)
    generated = request_screen.render_outputs_and_parameters(
        st, inputs["material"], st.session_state.get("request_mapping"), inputs["model"],
        mapping_upload=st.session_state.get("request_upload_mapping"),
        request_supplied=bool(inputs["request"] or uploads["request"] is not None))
    output = st.text_input("Output directory", str(DEFAULT_WORKDIR / "request"), key="request_output")
    with st.expander("Advanced"):
        mapping = st.text_input("Mapping.json path (optional with adjacent generated sidecar)", key="request_mapping")
        mapping_upload = st.file_uploader("Mapping.json", type=["json"], key="request_upload_mapping")
        abaqus = st.text_input("Abaqus executable", "abaqus", key="request_abaqus")
        validate = st.checkbox("Independent finite-difference validation", value=False, key="request_validate")
    if st.button("Solve", key="btn_request_run", type="primary",
                 disabled=not output or not all(inputs[key] or uploads[key] is not None for key in inputs
                                                if key != "request")
                 or not (inputs["request"] or uploads["request"] is not None or generated)):
        st.session_state.pop("request_completed_output", None)
        st.session_state.pop("request_completed_message", None)
        try:
            with tempfile.TemporaryDirectory(prefix="resasm_request_upload_") as temporary:
                for key, upload in uploads.items():
                    if upload is not None:
                        path = Path(temporary) / {"material": "OTI_UMAT.obj", "model": "Analysis.inp",
                                                 "odb": "Analysis.odb", "request": "sensitivity_request.json"}[key]
                        path.write_bytes(upload.getvalue())
                        inputs[key] = str(path)
                if not inputs["request"] and generated:
                    inputs["request"] = str(request_screen.write_request(generated, Path(temporary)))
                if mapping_upload is not None:
                    path = Path(temporary) / "Mapping.json"
                    path.write_bytes(mapping_upload.getvalue())
                    mapping = str(path)
                argv = ["request", "--out", output, "--abaqus", abaqus]
                for key, value in inputs.items():
                    argv.extend(["--" + key, value])
                if mapping:
                    argv.extend(["--mapping", mapping])
                if validate:
                    argv.append("--validate")
                execution = _run(argv)
                if execution.code:
                    raise ValueError(execution.stderr or execution.stdout)
                result = json.loads((Path(output) / "sensitivity_results.json").read_text())
            st.session_state["request_completed_output"] = str(Path(output).resolve())
            st.session_state["request_completed_message"] = "Executed: %d scalar results. Independent validation: %s." % (
                len(result["results"]), "passed" if result["metadata"]["verified"]
                else "run, NOT passed (see run_report.txt)" if validate else "not run")
        except (ValueError, OSError) as error:
            st.error(str(error))
    completed = st.session_state.get("request_completed_output")
    if completed:
        st.success(st.session_state["request_completed_message"])
        for filename in ("sensitivity_results.json", "sensitivity_tables.csv", "run_report.txt"):
            path = Path(completed) / filename
            if path.is_file():
                st.download_button(filename, path.read_bytes(), file_name=filename, key="request_download_" + filename)
    request_screen.render_full_field(st, completed)


def _tab_replay() -> None:
    st.header("7. Connected J2 Replay")
    record = st.text_input("Record or model JSON", str(REPO_ROOT / "examples/bounded_j2_c3d8/model.json"), key="replay_record")
    material_object = st.text_input("Compiled provider object", key="replay_object")
    contract = st.text_input("Provider contract JSON", key="replay_contract")
    output = st.text_input("Replay output directory", str(DEFAULT_WORKDIR / "replay"), key="replay_output")
    solve = st.checkbox("Solve model", value=True, key="replay_solve")
    verify = st.checkbox("Verify with ORIGINAL finite differences", value=True, key="replay_verify")
    argv = ["replay", record, "--object", material_object, "--contract", contract, "--out", output]
    if solve:
        argv.append("--solve")
    if verify:
        argv.append("--verify")
    _action("Run replay", "replay_run", argv, cwd=REPO_ROOT,
            disabled=not (material_object and contract and record and output))


def main() -> None:
    st.set_page_config(page_title="Residual_Assembler", layout="wide",
                       page_icon="🧮", initial_sidebar_state="collapsed")
    _init_state()
    _sidebar()

    tabs = st.tabs([
        "Sensitivity Request",
        "Start here",
        "1. Model",
        "2. Requirements",
        "3. Assemble",
        "4. Sensitivity",
        "5. Job",
        "6. Backends",
        "Advanced Replay",
    ])
    with tabs[0]:
        _tab_request()
    with tabs[1]:
        _tab_start()
    with tabs[2]:
        _tab_model()
    with tabs[3]:
        _tab_requirements()
    with tabs[4]:
        _tab_assemble()
    with tabs[5]:
        _tab_sensitivity()
    with tabs[6]:
        _tab_job()
    with tabs[7]:
        _tab_backends()
    with tabs[8]:
        _tab_replay()


if __name__ == "__main__":
    main()
