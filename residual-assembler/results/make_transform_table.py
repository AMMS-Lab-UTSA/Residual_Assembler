#!/usr/bin/env python3
"""Render the UMAT source-transformation verification table (Program 1) as a
figure in the study's style, from the Abaqus-based batch results.

Reads umat_oti_workspace/validate_all/summary.json (per-UMAT original-vs-OTI
comparison: STRESS / STATEV / DDSDDE against the original UMAT run in Abaqus) and
renders results/figures/umat_transform_table.png plus a CSV.

    UMAT_OTI_ROOT=~/Documents/UMAT_source_transformation \
      python results/make_transform_table.py
"""
import csv, json, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RA = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
P1 = os.path.expanduser(os.environ.get("UMAT_OTI_ROOT", "~/Documents/UMAT_source_transformation"))
FIGDIR = os.path.join(RA, "results", "figures"); os.makedirs(FIGDIR, exist_ok=True)
SUMMARY = os.path.join(P1, "umat_oti_workspace", "validate_all", "summary.json")


def _fmt(a, r):
    def g(x):
        if x in (0, 0.0):
            return "0"
        return ("%.3g" % x)
    return "%s / %s" % (g(a), g(r))


def main():
    rows = json.load(open(SUMMARY))
    rows = sorted(rows, key=lambda r: r["case"])
    header = ["Example UMAT", "DDSDDE max diff\n(abs / rel)", "Difference explanation", "Physics model"]
    data, colors, npass = [], [], 0
    for r in rows:
        ok = bool(r.get("overall_pass"))
        transform_ok = bool(r.get("transform_ok"))
        npass += int(ok and transform_ok)
        reason = r.get("difference_reason", "")
        if len(reason) > 52:
            reason = reason[:49] + "..."
        data.append([r["case"], _fmt(r.get("ddsdde_max_abs", 0), r.get("ddsdde_max_rel", 0)),
                     reason, r.get("physics_model", "")])
        # green = verified, amber = documented difference, red = transform/job failure
        if not transform_ok or r.get("transform_exit", 0) not in (0, None):
            c = "#f4b8b8"
        elif r.get("ddsdde_max_rel", 0) and r.get("ddsdde_max_rel", 0) > 1e-6:
            c = "#fde6b8"
        else:
            c = "#c9e7c0"
        colors.append([c] * 4)

    # CSV
    with open(os.path.join(FIGDIR, "umat_transform_table.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["UMAT", "ddsdde_abs", "ddsdde_rel", "stress_rel",
                                        "statev_pass", "physics_model", "reason", "pass"])
        for r in rows:
            w.writerow([r["case"], r.get("ddsdde_max_abs"), r.get("ddsdde_max_rel"),
                        r.get("stress_max_rel"), r.get("statev_pass"),
                        r.get("physics_model"), r.get("difference_reason"),
                        r.get("overall_pass")])

    fig, ax = plt.subplots(figsize=(12.5, 7.6)); ax.axis("off")
    tbl = ax.table(cellText=data, colLabels=header, cellColours=colors,
                   colColours=["#4a6d8c"] * 4, loc="center", cellLoc="left")
    tbl.auto_set_font_size(False); tbl.set_fontsize(8.5); tbl.scale(1, 1.5)
    for (row, col), cell in tbl.get_celld().items():
        if row == 0:
            cell.set_text_props(color="white", fontweight="bold")
        cell.set_edgecolor("white")
    tbl.auto_set_column_width([0, 1, 2, 3])
    ax.set_title("UMAT source-transformation verification  (%d/%d UMATs verified;\n"
                 "OTI STRESS, STATEV and DDSDDE compared with the original UMAT in Abaqus)"
                 % (npass, len(rows)), fontsize=12, pad=14)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "umat_transform_table.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("verified %d/%d UMATs; figure + CSV -> %s" % (npass, len(rows), FIGDIR))
    # a couple of noteworthy rows
    for r in rows:
        if r.get("ddsdde_max_rel", 0) and r.get("ddsdde_max_rel") > 1e-6:
            print("  %-16s DDSDDE %s  (%s)" % (r["case"],
                  _fmt(r.get("ddsdde_max_abs"), r.get("ddsdde_max_rel")), r.get("difference_reason")[:60]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
