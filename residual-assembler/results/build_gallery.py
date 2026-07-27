#!/usr/bin/env python3
"""Scrapeable results gallery: every figure as its own labelled, downloadable
block, and all the error data as REAL HTML tables (selectable/copyable text).

    UMAT_OTI_ROOT=~/Documents/UMAT_source_transformation python results/build_gallery.py
"""
import base64, json, os, sys, html

RA = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
FIG = os.path.join(RA, "results", "figures")


def datauri(name):
    return "data:image/png;base64," + base64.b64encode(open(os.path.join(FIG, name), "rb").read()).decode()


def jload(name):
    p = os.path.join(FIG, name)
    return json.load(open(p)) if os.path.exists(p) else {}


# figure registry: (file, title, one-line description)
GROUPS = [
    ("20-material framework sweep (both programs, per parameter)", [
        ("table_sweep_summary.png", "20 material models — P1 & P2 worst error per material",
         "The whole framework (transform + residual) validated across 20 UMATs spanning elasticity, plasticity, viscoplasticity, viscoelasticity and crystal plasticity. 18/20 exact to <1e-5 on both programs."),
        ("table_program1_sweep.png", "Program 1: DSIGMA_DP per parameter (84 parameters)",
         "OTI stress derivative vs central FD of the original UMAT, for every parameter of all 20 materials."),
        ("table_program2_sweep.png", "Program 2: residual d(sigma_vM)/dp per parameter (84 parameters)",
         "Residual-method output sensitivity vs full-analysis FD, for every parameter of all 20 materials."),
        ("sweep_oti_vs_fd_j2kin.png", "OTI vs FD along the path — J2 kinematic hardening (new material)",
         "A newly authored sweep material: OTI (lines) on central-FD (points) across the loading path, worst RMSE ~6e-10."),
    ]),
    ("Element formulations (residual method across element types)", [
        ("element_validation_table.png", "6 element types validated: C3D8/C3D8R, C3D20/C3D20R, C3D4, C3D10",
         "Residual method validated per element vs full-analysis finite differences: shape-function patch test ~1e-16, "
         "assembly ∂R/∂p ≤5e-9, full solve ∂u/∂p ≤9e-8. FCC crystal plasticity (10 params) flows through every element at ~3e-9."),
    ]),
    ("FCC single crystal (headline: 10 sensitivities, one run)", [
        ("fcc_weighted_regimes.png", "Weighted σ_vM sensitivities with regimes",
         "FCC [100], 12 slip systems. Stacked weighted sensitivities across elastic / yield / hardening / saturation regimes; σ_vM on the right axis."),
        ("fcc_oti_vs_fd.png", "OTI vs finite differences (9 parameters)",
         "(p/σ)∂σ/∂p along the load path, OTI (line) vs central FD (dashed). Relative RMSE ~1e-9."),
        ("fcc_timing.png", "Analysis-level timing (OTI ~15x faster than FD)",
         "Nominal analysis + one enriched replay for all 10 sensitivities, vs 2x10+1 = 21 full nonlinear analyses for finite differences."),
    ]),
    ("Viscoplastic flow model (M5) — residual method", [
        ("weighted_sigvm_sensitivities.png", "Weighted σ_vM sensitivities (single element)",
         "Kocks thermally-activated viscoplastic, 6 flow parameters, simple shear."),
        ("weighted_sigvm_sensitivities_4x4x4.png", "Weighted σ_vM sensitivities (4×4×4 mesh)",
         "Full residual solve on a 64-element mesh — identical to the single element (mesh independent)."),
        ("dsigvm_dp_oti_vs_fd.png", "OTI vs finite differences (6 parameters)",
         "dσ_vM/dp along the shear path, OTI vs central FD. Relative RMSE ~1e-9."),
        ("normalized_cost.png", "Normalized cost (M5)",
         "One enriched OTI run vs 13 finite-difference runs."),
    ]),
    ("Program 1 — source-transformation verification", [
        ("umat_transform_table.png", "18 / 19 UMAT verification suite",
         "OTI STRESS / STATEV / DDSDDE vs the original UMAT in Abaqus, across elasticity / plasticity / viscoplasticity / finite strain."),
        ("program1_dsigma_dstatev_vs_fd.png", "DSIGMA_DP & DSTATEV_DP vs finite differences (J2)",
         "Both transformer outputs — stress derivative and state-variable derivative — vs central FD along the loading path. OTI lines on FD points, worst rel. RMSE 8.9e-09."),
        ("table_program1_errors.png", "Program 1 per-parameter error table (image)",
         "Same data as the HTML table below, rendered as a figure for slides."),
        ("table_program2_errors.png", "Program 2 per-parameter error table (image)",
         "Same data as the HTML table below, rendered as a figure for slides."),
    ]),
]

DESC = {"umat_m1_elastic": "Isotropic linear elasticity",
        "umat_m2_cubic": "Cubic anisotropic elasticity",
        "umat_m3_j2": "J2 plasticity, linear hardening",
        "umat_m5_cpflow": "Thermally-activated viscoplastic flow",
        "umat_m6_fcc": "FCC single crystal, 12 slip systems"}


def err_rows(prog):
    et = jload("error_tables.json")
    rows = []
    for mid, blocks in et.items():
        for param, err in blocks.get(prog, {}).items():
            cls = "ok" if err < 1e-5 else ("warn" if err < 1e-3 else "bad")
            rows.append((mid, DESC.get(mid, ""), param, "%.2e" % err if err else "0", cls))
    return rows


def html_table(prog, caption):
    rows = err_rows(prog)
    body = "\n".join(
        f'<tr><td class="mono">{html.escape(m)}</td><td>{html.escape(d)}</td>'
        f'<td class="mono">{html.escape(p)}</td><td class="num {c}">{e}</td></tr>'
        for m, d, p, e, c in rows)
    return f"""<figure class="card"><figcaption class="tcap">{caption}</figcaption>
    <div class="scroll"><table><thead><tr><th>Material (file)</th><th>Physics</th>
    <th>Parameter</th><th class="num">OTI vs FD (rel.)</th></tr></thead>
    <tbody>{body}</tbody></table></div></figure>"""


def fig_card(fname, title, desc):
    uri = datauri(fname)
    return f"""<figure class="card fig">
    <div class="fighead"><h3>{html.escape(title)}</h3>
      <a class="dl" href="{uri}" download="{fname}">↓ {fname}</a></div>
    <img alt="{html.escape(title)}" src="{uri}">
    <figcaption>{html.escape(desc)}</figcaption></figure>"""


def main():
    fcc = jload("fcc_results.json")
    fcc_rmse = fcc.get("rmse_rel_oti_vs_fd", {})
    fcc_rows = "\n".join(
        f'<tr><td class="mono">{html.escape(k)}</td><td class="num ok">{v:.2e}</td></tr>'
        for k, v in fcc_rmse.items())

    groups_html = ""
    for gtitle, figs in GROUPS:
        cards = "\n".join(fig_card(f, t, d) for f, t, d in figs)
        groups_html += f'<section><h2>{html.escape(gtitle)}</h2><div class="grid">{cards}</div></section>\n'

    page = f"""<title>OTI residual-method results — gallery</title>
<style>
:root {{ --bg:#f6f7f9; --panel:#fff; --ink:#17212b; --muted:#5b6b7a; --line:#e2e6ea;
  --accent:#1f63c4; --ok:#2e8b57; --warn:#b6822a; --bad:#c0392b; }}
@media (prefers-color-scheme:dark){{ :root{{ --bg:#0f151b; --panel:#161f28; --ink:#e7eef4;
  --muted:#9fb0bf; --line:#26313d; --accent:#5aa0ff; --ok:#5cc98a; --warn:#e0b25a; --bad:#e8776b; }} }}
:root[data-theme=dark]{{ --bg:#0f151b; --panel:#161f28; --ink:#e7eef4; --muted:#9fb0bf; --line:#26313d;
  --accent:#5aa0ff; --ok:#5cc98a; --warn:#e0b25a; --bad:#e8776b; }}
:root[data-theme=light]{{ --bg:#f6f7f9; --panel:#fff; --ink:#17212b; --muted:#5b6b7a; --line:#e2e6ea;
  --accent:#1f63c4; --ok:#2e8b57; --warn:#b6822a; --bad:#c0392b; }}
*{{box-sizing:border-box}}
body{{margin:0;background:var(--bg);color:var(--ink);line-height:1.55;
  font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}}
.wrap{{max-width:70rem;margin:0 auto;padding:2.4rem 1.2rem 5rem}}
h1{{font-size:1.9rem;margin:0 0 .3rem;letter-spacing:-.02em}}
h2{{font-size:1.25rem;margin:2.4rem 0 .8rem;padding-top:1.2rem;border-top:1px solid var(--line)}}
h3{{font-size:1rem;margin:0}}
.lede{{color:var(--muted);margin:.2rem 0 1rem;max-width:52rem}}
.mono,.num{{font-variant-numeric:tabular-nums;font-family:ui-monospace,Menlo,Consolas,monospace}}
.grid{{display:grid;gap:1.1rem}}
@media(min-width:820px){{.grid{{grid-template-columns:1fr 1fr}}}}
.card{{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:.9rem 1rem 1rem;
  box-shadow:0 1px 2px rgba(20,30,40,.04),0 6px 18px rgba(20,30,40,.05)}}
.fig img{{width:100%;height:auto;display:block;border-radius:6px;margin:.5rem 0}}
.fighead{{display:flex;justify-content:space-between;align-items:baseline;gap:.6rem}}
.dl{{font-size:.76rem;font-weight:600;color:var(--accent);text-decoration:none;white-space:nowrap;
  border:1px solid var(--line);border-radius:6px;padding:.15rem .5rem}}
.dl:hover{{background:color-mix(in srgb,var(--accent) 12%,transparent)}}
figcaption{{font-size:.82rem;color:var(--muted)}}
.tcap{{font-size:.9rem;font-weight:600;color:var(--ink);margin-bottom:.5rem}}
.scroll{{overflow-x:auto}}
table{{border-collapse:collapse;width:100%;font-size:.82rem}}
th,td{{text-align:left;padding:.32rem .55rem;border-bottom:1px solid var(--line)}}
th{{font-size:.72rem;text-transform:uppercase;letter-spacing:.04em;color:var(--muted)}}
td.num,th.num{{text-align:right}}
td.ok{{color:var(--ok)}} td.warn{{color:var(--warn)}} td.bad{{color:var(--bad)}}
.tables{{display:grid;gap:1.1rem}}
@media(min-width:820px){{.tables{{grid-template-columns:1fr 1fr}}}}
.hint{{font-size:.8rem;color:var(--muted);background:color-mix(in srgb,var(--accent) 8%,transparent);
  border:1px solid var(--line);border-radius:8px;padding:.5rem .8rem;margin:.6rem 0 0}}
</style>

<div class="wrap">
  <h1>OTI residual-method results — gallery</h1>
  <p class="lede">Every figure and data table from the two-program validation, laid out to grab individually:
  each plot has a download link, and the error data is in real HTML tables you can select and copy.</p>
  <p class="hint">Tip: use the ↓ link on any figure to save the PNG, or select a table's cells to copy the numbers.</p>

  {groups_html}

  <section><h2>Per-parameter error data (copyable tables)</h2>
  <div class="tables">
    {html_table('program1', 'Program 1 — DSIGMA_DP per parameter vs FD of the original UMAT')}
    {html_table('program2', 'Program 2 — ∂σ_vM/∂p per parameter vs full-analysis FD')}
  </div></section>

  <section><h2>FCC per-parameter relative RMSE (OTI vs FD)</h2>
  <figure class="card"><div class="scroll"><table>
    <thead><tr><th>Parameter</th><th class="num">rel. RMSE</th></tr></thead>
    <tbody>{fcc_rows}</tbody></table></div>
    <figcaption>FCC single crystal [100], all ten sensitivities from one enriched OTI run;
    σ_vM at 5% strain = {fcc.get('sigma_vM_final_MPa', 0):.0f} MPa. C44 ≈ 0 (no shear in [100]).</figcaption>
  </figure></section>
</div>
"""
    out = os.path.join(RA, "results", "gallery.html")
    open(out, "w").write(page)
    print("wrote", out, "(%d KB)" % (len(page) // 1024))
    return 0


if __name__ == "__main__":
    sys.exit(main())
