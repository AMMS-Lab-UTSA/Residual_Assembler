#!/usr/bin/env python3
"""Evidence for the presentation cantilevers (slides 15, 33, 39) with the history replay.

Runs ``resasm history`` on the Abaqus data of the J2 cantilever (slide 39) and
the FCC crystal-plasticity cantilever (slide 15) -- recorded-state and
re-equilibrated -- then

* compares the OTI sensitivities with central differences of the perturbed
  Abaqus reruns (read directly from their field exports; single-precision
  floor and step-size plateau respected);
* writes the weighted per-increment sensitivity shares (slide 33) as CSV and
  Markdown tables and as a stacked chart;
* plots p * d(sigma_vM)/dp at the final increment for every parameter.

Large outputs (full fields) go to ``--work``; the small tables and figures go
to ``--evidence`` (committed). Nothing here is used by the library.

    python scripts/replay_history_cantilevers.py \\
        --abaqus-data /path/imq_abaqus/claude_cantilevers \\
        --providers /path/imq_abaqus/claude_C/providers \\
        --work /path/imq_abaqus/claude_C/cantilevers --evidence docs/evidence/claude_C
"""
import argparse
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from residual_core.replay.history import von_mises, von_mises_gradient  # noqa: E402
from residual_core.replay.history_inputs import read_history_model  # noqa: E402
from residual_core.ui.cmd_history import run_history_request  # noqa: E402

EPS32 = float(np.finfo(np.float32).eps)
MODELS = {
    "j2": {"deck": "j2/claude_j2_nominal.inp", "fields": "j2/claude_j2_nominal_fields.npz",
           "provider": "prov_j2/umat_m3_j2_oti.obj", "fd_dir": "j2",
           "props_order": ["E", "nu", "SIGY0", "H"],
           "labels": {"E": "E", "nu": "ν", "SIGY0": "σy0", "H": "H"}},
    "fcc": {"deck": "fcc/claude_fcc_nominal.inp", "fields": "fcc/claude_fcc_nominal_fields.npz",
            "provider": "prov_fcc/umat_m6_fcc_oti.obj", "fd_dir": "fcc",
            "props_order": ["C11", "C12", "C44", "g0", "gsat", "h0", "a", "q", "gd0", "m"],
            "labels": {"g0": "g0", "h0": "h0", "q": "q", "gd0": "γ̇0", "m": "m", "gsat": "gsat",
                       "C11": "C11", "C12": "C12", "C44": "C44", "a": "a"}},
}
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
OTHER = "#8a8984"
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"


def request_for(model):
    return {
        "outputs": [
            {"name": "tip_RF2", "field": "RF", "component": 2, "reduction": "sum", "domain": {"nset": "TIP"}},
            {"name": "mises_mean", "field": "MISES", "component": 1, "reduction": "volume_mean",
             "domain": {"elements": "ALL"}},
            {"name": "mises_root_mean", "field": "MISES", "component": 1, "reduction": "volume_mean",
             "domain": {"elset": "ROOTEL"}},
            {"name": "e1_ip1_S11", "field": "S", "component": 1, "reduction": "component",
             "domain": {"elements": [1], "points": [1]}}],
        "parameters": "ALL", "domain": {"nodes": "ALL", "elements": "ALL"}, "increments": "ALL",
        "weighted_shares": {"field": "MISES", "domain": {"elements": "ALL"}}, "full_field": True}


def run(model, data, providers, work, reequilibrate):
    spec = MODELS[model]
    out = work / ("%s_%s" % (model, "reequilibrated" if reequilibrate else "recorded"))
    if (out / "sensitivity_results.json").is_file():
        return out, json.loads((out / "sensitivity_results.json").read_text())
    work.mkdir(parents=True, exist_ok=True)
    request = work / ("%s_request.json" % model)
    request.write_text(json.dumps(request_for(model), indent=1))
    report = run_history_request(model=data / spec["deck"], fields=data / spec["fields"],
                                 material=providers / spec["provider"], request=request, out=out,
                                 reequilibrate=reequilibrate,
                                 command="resasm history --model %s --fields %s --material %s --request %s "
                                         "--out %s%s" % (spec["deck"], spec["fields"], spec["provider"],
                                                         request.name, out.name,
                                                         " --reequilibrate" if reequilibrate else ""))
    return out, json.loads((out / "sensitivity_results.json").read_text())


# ------------------------------------------------------------------ Abaqus FD
def fd_outputs(npz, model):
    d = np.load(npz)
    coords, conn = d["coords"], d["conn"]
    L, H = coords[:, 0].max(), coords[:, 1].max()
    tip = np.isclose(coords[:, 0], L)
    labels = list(d["node_labels"])
    index = {int(label): i for i, label in enumerate(labels)}
    X = coords[np.vectorize(index.get)(conn)]
    xs = X[:, :, 0].mean(axis=1)
    region = coords[:, 0] >= L - 4 - 1e-9                      # tip region: last four element columns
    root_col = np.isclose(xs, xs.min())
    s = d["S"]
    vm = von_mises(s)
    return {"tip_RF2": d["RF"][:, tip, 1].sum(axis=1),
            "tip_region_U": d["U"][:, region, :].reshape(len(d["time"]), -1),
            "U_all": d["U"].reshape(len(d["time"]), -1),
            "e1_ip1_S11": s[:, 0, 0, 0],
            "root_mises": vm[:, root_col, :].reshape(len(d["time"]), -1),
            "_region": region, "_root": root_col}


def fd_reference(data, model, names, values):
    manifest = data / MODELS[model]["fd_dir"] / "fd_manifest.tsv"
    if not manifest.is_file():
        return None
    rows = [line.split("\t") for line in manifest.read_text().splitlines()[1:] if line.strip()]
    done = {(r[1], r[2], r[3]): r[0] for r in rows if len(r) > 5 and r[5] == "completed"}
    reference = {}
    for j, name in enumerate(names):
        per_h = {}
        for h in ("0.02", "0.01", "0.005"):
            if (name, "p", h) in done and (name, "m", h) in done:
                plus = fd_outputs(data / MODELS[model]["fd_dir"] / (done[(name, "p", h)] + "_fields.npz"), model)
                minus = fd_outputs(data / MODELS[model]["fd_dir"] / (done[(name, "m", h)] + "_fields.npz"), model)
                per_h[h] = {k: (plus[k] - minus[k]) / (2 * float(h) * values[j])
                            for k in plus if not k.startswith("_")}
        if len(per_h) >= 2:
            reference[name] = per_h
    return reference


def compare_with_abaqus(result_dir, data, model, names, values):
    reference = fd_reference(data, model, names, values)
    if not reference:
        return None
    fields = np.load(result_dir / "fields.npz")
    base = fd_outputs(data / MODELS[model]["fields"], model)
    region, root = base["_region"], base["_root"]
    dU, dRF, dS = fields["dU"], fields["dRF"], fields["dS"]
    tip = np.isclose(np.load(data / MODELS[model]["fields"])["coords"][:, 0],
                     np.load(data / MODELS[model]["fields"])["coords"][:, 0].max())
    S = fields["S"]
    dvm = np.einsum("neqa,neqam->neqm", von_mises_gradient(S), dS)
    count = dU.shape[0]
    ours_all = {
        "tip_RF2": dRF[:, tip, 1, :].sum(axis=1),
        "tip_region_U": dU[:, region, :, :].reshape(count, -1, dU.shape[-1]),
        "e1_ip1_S11": dS[:, 0, 0, 0, :],
        "root_mises": dvm[:, root, :, :].reshape(count, -1, dU.shape[-1]),
    }
    table = []
    for j, name in enumerate(names):
        if name not in reference:
            continue
        steps = sorted(reference[name], key=float, reverse=True)
        for key, ours in ours_all.items():
            q = base[key]
            for n in range(1, count):
                estimates = [np.atleast_1d(reference[name][h][key][n]) for h in steps]
                mine = np.atleast_1d(ours[n][..., j])
                fine, finer = estimates[-2], estimates[-1]
                scale = max(np.abs(finer).max(), 1e-300)
                natural = max(np.abs(np.atleast_1d(q[n])).max(), 1e-300) / abs(values[j])
                floor = EPS32 * np.abs(np.atleast_1d(q[n])).max() / (2 * float(steps[-1]) * abs(values[j]))
                table.append({
                    "output": key, "parameter": name, "increment": n,
                    "fd": float(finer.flat[np.abs(finer).argmax()]),
                    "oti": float(mine.flat[np.abs(finer).argmax()]),
                    "rel_error": float(np.abs(mine - finer).max() / scale),
                    "plateau_spread": float(np.abs(fine - finer).max() / scale),
                    "truncation_spread": float(np.abs(estimates[0] - finer).max() / scale),
                    "weighted_fd": float(scale / natural), "float32_floor_rel": float(floor / scale),
                })
    return table


def summarize_abaqus(table, resolved_spread=1e-3):
    summary = {}
    for row in table:
        entry = summary.setdefault((row["output"], row["parameter"]), {
            "resolved": 0, "unresolved": 0, "zero": 0, "max_rel_error": 0.0, "max_spread": 0.0,
            "max_error_over_spread": 0.0, "zero_max_weighted_oti": 0.0})
        if row["weighted_fd"] < 1e-5:
            entry["zero"] += 1
            continue
        if row["plateau_spread"] <= resolved_spread and row["float32_floor_rel"] < 1e-4:
            entry["resolved"] += 1
            entry["max_rel_error"] = max(entry["max_rel_error"], row["rel_error"])
            entry["max_spread"] = max(entry["max_spread"], row["plateau_spread"])
            entry["max_error_over_spread"] = max(entry["max_error_over_spread"],
                                                 row["rel_error"] / max(row["plateau_spread"], row["float32_floor_rel"]))
        else:
            entry["unresolved"] += 1
    return summary


# ------------------------------------------------------------------ figures
def _style(plt):
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK2,
                         "xtick.color": INK2, "ytick.color": INK2, "axes.facecolor": SURFACE,
                         "figure.facecolor": SURFACE, "savefig.facecolor": SURFACE,
                         "axes.spines.top": False, "axes.spines.right": False})


def shares_figure(results, model, path, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _style(plt)
    shares = results["weighted_shares"]
    names = shares["parameters"]
    rows = shares["rows"]
    steps = np.array([r["increment"] for r in rows])
    data = np.array([[r["field_share_percent"][n] for n in names] for r in rows])
    order = np.argsort(-data[-1])
    labels = MODELS[model]["labels"]
    if len(names) > 8:
        keep = list(order[:7])
        other = [i for i in range(len(names)) if i not in keep]
        series = [(labels[names[i]], data[:, i], SERIES[k]) for k, i in enumerate(keep)]
        series.append(("other (%s)" % ", ".join(labels[names[i]] for i in other), data[:, other].sum(axis=1), OTHER))
    else:
        series = [(labels[names[i]], data[:, i], SERIES[k]) for k, i in enumerate(order)]
    fig, ax = plt.subplots(figsize=(7.4, 3.6), dpi=130)
    bottom = np.zeros(len(steps))
    for label, values, color in series:
        ax.bar(steps, values, bottom=bottom, width=0.82, color=color, edgecolor=SURFACE, linewidth=0.6,
               label=label)
        bottom += values
    ax.set_xlim(steps[0] - 0.6, steps[-1] + 0.6)
    ax.set_ylim(0, 100)
    ax.set_xlabel("load step (increment)")
    ax.set_ylabel("weighted share of |p ∂σvM/∂p| (%)")
    ax.set_title(title, loc="left", color=INK, fontsize=10)
    ax.grid(axis="y", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    handles, texts = ax.get_legend_handles_labels()            # legend top = stack top
    ax.legend(handles[::-1], texts[::-1], loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False,
              fontsize=8, labelcolor=INK)
    last = 0.0
    for label, values, color in series:                          # direct labels for segments >= 6 %
        if values[-1] >= 6:
            ax.text(steps[-1] + 0.55, last + values[-1] / 2, "%s %.0f%%" % (label.split(" (")[0], values[-1]),
                    va="center", ha="left", fontsize=7, color=INK)
        last += values[-1]
    fig.tight_layout()
    fig.savefig(path, dpi=130, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)


def field_figure(result_dir, deck, model, path, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
    _style(plt)
    fields = np.load(result_dir / "fields.npz")
    names = [str(n) for n in fields["parameters"]]
    values = fields["parameter_values"]
    dvm = fields["dMISES"][-1]                                 # (ne, 8, npar)
    weighted = dvm * values                                    # p d sigma_vM / dp  (MPa)
    X = deck.coords[deck.connectivity]
    centre = X.mean(axis=1)
    nx = int(round(deck.coords[:, 0].max())); ny = int(round(deck.coords[:, 1].max()))
    front = centre[:, 2] < 1.0                                 # z = 0 element layer
    order = MODELS[model]["props_order"]
    labels = MODELS[model]["labels"]
    maps = {}
    for j, name in enumerate(names):
        grid = np.full((ny, nx), np.nan)
        cells = np.flatnonzero(front)
        grid[np.floor(centre[cells, 1]).astype(int), np.floor(centre[cells, 0]).astype(int)] = \
            weighted[cells, :, j].mean(axis=1)
        maps[name] = grid
    limit = max(np.nanmax(np.abs(g)) for g in maps.values())
    cmap = LinearSegmentedColormap.from_list("div", ["#1c5cab", "#86b6ef", "#f0efec", "#f19a8e", "#b3261e"])
    columns = 2
    rows = int(np.ceil(len(order) / columns))
    fig, axes = plt.subplots(rows, columns, figsize=(7.4, 0.2 + rows * (7.4 / columns) * ny / nx * 1.25 + 0.3),
                             dpi=130, squeeze=False)
    norm = TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)
    for ax, name in zip(axes.flat, order):
        image = ax.imshow(maps[name], origin="lower", cmap=cmap, norm=norm, extent=(0, nx, 0, ny),
                          interpolation="nearest")
        ax.set_title("p ∂σvM/∂p,  p = %s   (max |.| %.3g MPa)" % (labels[name], np.nanmax(np.abs(maps[name]))),
                     fontsize=8, color=INK, loc="left")
        ax.set_xticks([]); ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
    for ax in axes.flat[len(order):]:
        ax.axis("off")
    fig.suptitle(title, x=0.01, ha="left", fontsize=10, color=INK)
    bar = fig.colorbar(image, ax=axes, orientation="horizontal", fraction=0.04, pad=0.04, aspect=50)
    bar.set_label("MPa  (element average over the 8 IPs, z = 0 layer; root at left, tip at right)",
                  color=INK2, fontsize=8)
    fig.savefig(path, dpi=130, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)


def steps_figure(result_dir, deck, model, path, parameter, steps, title):
    """One panel per load step: p d(sigma_vM)/dp of one parameter (slide 33 layout)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
    _style(plt)
    fields = np.load(result_dir / "fields.npz")
    names = [str(n) for n in fields["parameters"]]
    j = names.index(parameter)
    value = fields["parameter_values"][j]
    X = deck.coords[deck.connectivity]
    centre = X.mean(axis=1)
    nx = int(round(deck.coords[:, 0].max())); ny = int(round(deck.coords[:, 1].max()))
    cells = np.flatnonzero(centre[:, 2] < 1.0)
    maps = []
    for n in steps:
        grid = np.full((ny, nx), np.nan)
        grid[np.floor(centre[cells, 1]).astype(int), np.floor(centre[cells, 0]).astype(int)] = \
            (fields["dMISES"][n][cells, :, j] * value).mean(axis=1)
        maps.append(grid)
    limit = max(np.nanmax(np.abs(g)) for g in maps)
    cmap = LinearSegmentedColormap.from_list("div", ["#1c5cab", "#86b6ef", "#f0efec", "#f19a8e", "#b3261e"])
    columns = 3
    rows = int(np.ceil(len(steps) / columns))
    fig, axes = plt.subplots(rows, columns, figsize=(7.4, 0.5 + rows * (7.4 / columns) * ny / nx * 1.35 + 0.3),
                             dpi=130, squeeze=False)
    norm = TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)
    for ax, n, grid in zip(axes.flat, steps, maps):
        image = ax.imshow(grid, origin="lower", cmap=cmap, norm=norm, extent=(0, nx, 0, ny),
                          interpolation="nearest")
        ax.set_title("step %d" % n, fontsize=8, color=INK, loc="left")
        ax.set_xticks([]); ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
    for ax in axes.flat[len(steps):]:
        ax.axis("off")
    fig.suptitle(title, x=0.01, ha="left", fontsize=10, color=INK)
    bar = fig.colorbar(image, ax=axes, orientation="horizontal", fraction=0.05, pad=0.05, aspect=50)
    bar.set_label("MPa  (element average, z = 0 layer; root at left)", color=INK2, fontsize=8)
    fig.savefig(path, dpi=130, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)


def shares_tables(results, model, evidence, suffix):
    shares = results["weighted_shares"]
    names = shares["parameters"]
    with open(evidence / ("%s_shares_%s.csv" % (model, suffix)), "w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["increment", "time", "volume_mean_mises"] + ["field_share_%s" % n for n in names]
                        + ["scalar_share_%s" % n for n in names])
        for r in shares["rows"]:
            writer.writerow([r["increment"], r["time"], r["volume_mean"]]
                            + [r["field_share_percent"][n] for n in names]
                            + [r["scalar_share_percent"][n] for n in names])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--abaqus-data", required=True, type=Path)
    ap.add_argument("--providers", required=True, type=Path)
    ap.add_argument("--work", required=True, type=Path)
    ap.add_argument("--evidence", required=True, type=Path)
    ap.add_argument("--models", default="j2,fcc")
    a = ap.parse_args(argv)
    a.evidence.mkdir(parents=True, exist_ok=True)
    summary = {}
    for model in a.models.split(","):
        deck = read_history_model(a.abaqus_data / MODELS[model]["deck"])
        entry = summary[model] = {}
        for reequilibrate in (False, True):
            tag = "reequilibrated" if reequilibrate else "recorded"
            started = time.perf_counter()
            out, results = run(model, a.abaqus_data, a.providers, a.work, reequilibrate)
            meta = results["metadata"]
            names = results["weighted_shares"]["parameters"]
            values = np.array([meta["parameter_values"][n] for n in names])
            shares_tables(results, model, a.evidence, tag)
            rows = results["weighted_shares"]["rows"]
            plastic = []
            fields = np.load(out / "fields.npz")
            for n in range(1, fields["SDV"].shape[0]):
                plastic.append(int((np.abs(fields["SDV"][n]).max(axis=-1) > 0).sum()))
            entry[tag] = {
                "timings_s": meta["timings_s"], "max_scaled_free_residual": meta["max_scaled_free_residual"],
                "parity": meta["parity"], "points_with_nonzero_state": plastic,
                "shares_field": {r["increment"]: r["field_share_percent"] for r in rows},
                "shares_scalar": {r["increment"]: r["scalar_share_percent"] for r in rows},
                "wall_s": round(time.perf_counter() - started, 1)}
            abaqus = compare_with_abaqus(out, a.abaqus_data, model, names, values)
            if abaqus:
                with open(a.evidence / ("%s_abaqus_fd_%s.csv" % (model, tag)), "w", newline="") as stream:
                    writer = csv.DictWriter(stream, fieldnames=list(abaqus[0]))
                    writer.writeheader()
                    writer.writerows(abaqus)
                entry[tag]["abaqus_fd"] = {"%s/%s" % key: value
                                           for key, value in summarize_abaqus(abaqus).items()}
            if reequilibrate:
                title = {"j2": "J2 cantilever (slide 39): weighted σvM sensitivity shares per load step",
                         "fcc": "FCC cantilever (slide 15): weighted σvM sensitivity shares per load step"}[model]
                shares_figure(results, model, a.evidence / ("%s_weighted_shares.png" % model), title)
                field_figure(out, deck, model, a.evidence / ("%s_final_step_dmises.png" % model),
                             {"j2": "J2 cantilever, step 40: p ∂σvM/∂p for every parameter",
                              "fcc": "FCC cantilever, step 25: p ∂σvM/∂p for every parameter"}[model])
                plastic_steps = [n for n, count in enumerate(plastic, 1) if count]
                first = plastic_steps[0] if plastic_steps else 1
                last = len(plastic)
                picks = sorted(set([first, first + 3, (first + last) // 3, (first + last) // 2,
                                    (2 * last + first) // 3, last]))
                parameter = {"j2": "SIGY0", "fcc": "g0"}[model]
                steps_figure(out, deck, model, a.evidence / ("%s_steps_%s.png" % (model, parameter)),
                             parameter, picks,
                             "%s cantilever: p ∂σvM/∂p for p = %s, one panel per load step"
                             % (model.upper(), MODELS[model]["labels"][parameter]))
        recorded = np.load(a.work / ("%s_recorded" % model) / "fields.npz")
        polished = np.load(a.work / ("%s_reequilibrated" % model) / "fields.npz")
        entry["recorded_vs_reequilibrated_rel"] = {
            key: (np.abs(recorded[key] - polished[key]).max(axis=tuple(range(recorded[key].ndim - 1)))
                  / np.maximum(np.abs(polished[key]).max(axis=tuple(range(polished[key].ndim - 1))), 1e-300)).tolist()
            for key in ("dU", "dRF", "dS", "dMISES")}
    (a.evidence / "cantilever_summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({m: {t: {k: v for k, v in e.items() if k in ("timings_s", "max_scaled_free_residual")}
                          for t, e in summary[m].items() if isinstance(e, dict) and "timings_s" in e}
                      for m in summary}, indent=1))


if __name__ == "__main__":
    main()
