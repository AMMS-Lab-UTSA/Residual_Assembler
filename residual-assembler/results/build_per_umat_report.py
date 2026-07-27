#!/usr/bin/env python3
"""Comprehensive per-UMAT HTML report: for every material in the 20-UMAT sweep,
the Program-1 transform contract, the Program-2 residual request, the weighted
sensitivity plot, OTI-vs-FD comparisons for both programs, and a per-parameter
error table. Everything is embedded (base64 images, inline JSON) -> one
self-contained, scrapeable HTML.

    python results/build_per_umat_report.py        # -> results/per_umat_report.html
"""
import base64
import glob
import html
import json
import os
import re
import sys

RA = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
sys.path.insert(0, os.path.join(RA, "results"))
import build_sweep_tables as B                                       # noqa: E402

PU = os.path.join(RA, "results", "figures", "per_umat")


def datauri(path):
    return "data:image/png;base64," + base64.b64encode(open(path, "rb").read()).decode()


def hl_json(obj):
    """Minimal JSON syntax highlight -> HTML."""
    s = json.dumps(obj, indent=2)
    s = html.escape(s)
    s = re.sub(r'(&quot;.*?&quot;)(\s*:)', r'<span class="jk">\1</span>\2', s)   # keys
    s = re.sub(r':\s(&quot;.*?&quot;)', r': <span class="js">\1</span>', s)      # string values
    s = re.sub(r'(?<![\w&])(-?\d+\.?\d*(?:[eE][-+]?\d+)?)', r'<span class="jn">\1</span>', s)
    s = re.sub(r'\b(true|false|null)\b', r'<span class="jb">\1</span>', s)
    return s


def status_class(v):
    v = float(v)
    return "ok" if v < 1e-5 else ("warn" if v < 1e-3 else "bad")


def order_key(mdir):
    return (0 if not mdir.startswith("sweep_") else 1, mdir)


def load():
    mats = []
    for dj in sorted(glob.glob(os.path.join(PU, "*_data.json"))):
        d = json.load(open(dj))
        mats.append(d)
    mats.sort(key=lambda d: order_key(d["mdir"]))
    return mats


def err_table(d):
    p1, p2 = d.get("p1_errors", {}), d.get("p2_errors", {})
    rows = ""
    for nm in d["params"]:
        e1 = p1.get(nm); e2 = p2.get(nm)
        c1 = ('<td class="num %s">%s</td>' % (status_class(e1), "%.1e" % e1)) if e1 is not None else '<td class="num">-</td>'
        c2 = ('<td class="num %s">%s</td>' % (status_class(e2), "%.1e" % e2)) if e2 is not None else '<td class="num">-</td>'
        rows += '<tr><td class="mono">%s</td>%s%s</tr>' % (html.escape(nm), c1, c2)
    return ('<table class="err"><thead><tr><th>parameter</th><th class="num">Program&nbsp;1<br>DSIGMA_DP vs FD</th>'
            '<th class="num">Program&nbsp;2<br>residual vs FD</th></tr></thead><tbody>%s</tbody></table>' % rows)


def section(d):
    mdir = d["mdir"]
    badge_type = "path-dependent" if d["path_dependent"] else "elastic"
    p1w, p2w = d["p1_worst"], d["p2_worst"]
    # contracts
    p1c = json.load(open(d["transform_contract"])) if os.path.exists(d["transform_contract"]) else {}
    p2r_path = os.path.join(PU, "%s_p2_request.json" % mdir)
    p2c = json.load(open(p2r_path)) if os.path.exists(p2r_path) else {}
    # plots
    def img(sfx, alt):
        p = os.path.join(PU, "%s_%s.png" % (mdir, sfx))
        return ('<figure><img alt="%s" src="%s"><figcaption>%s</figcaption></figure>'
                % (html.escape(alt), datauri(p), html.escape(alt))) if os.path.exists(p) else ""
    props = ", ".join("%g" % v for v in (d.get("props_values") or []))
    return f"""
  <section id="{html.escape(mdir)}" class="mat">
    <div class="mathead">
      <h2>{html.escape(d['physics'])}</h2>
      <div class="badges">
        <span class="badge">{html.escape(mdir)}</span>
        <span class="badge">{len(d['params'])} parameters</span>
        <span class="badge {'pd' if d['path_dependent'] else 'el'}">{badge_type}</span>
        <span class="badge {status_class(p1w)}">P1 {p1w:.1e}</span>
        <span class="badge {status_class(p2w)}">P2 {p2w:.1e}</span>
      </div>
    </div>
    <p class="props"><b>parameters:</b> {html.escape(", ".join(d['params']))} &nbsp;·&nbsp; <b>PROPS:</b> [{html.escape(props)}]</p>
    <div class="grid2">
      <div class="card"><div class="ct">Program&nbsp;1 &mdash; transform contract <span class="tag">resasm_umat_transform_v2</span></div>
        <pre class="code">{hl_json(p1c)}</pre></div>
      <div class="card"><div class="ct">Program&nbsp;2 &mdash; residual request <span class="tag">resasm_sensitivity_request_v1</span></div>
        <pre class="code">{hl_json(p2c)}</pre>
        {err_table(d)}</div>
    </div>
    <div class="plots">
      {img('weighted', 'Weighted sigma_vM parameter sensitivities')}
      {img('p1', 'Program 1: DSIGMA_DP vs FD of the original UMAT')}
      {img('p2', 'Program 2: residual d(sigma_vM)/dp vs full-analysis FD')}
    </div>
  </section>"""


def main():
    mats = load()
    nboth = sum(1 for d in mats if float(d["p1_worst"]) < 1e-5 and float(d["p2_worst"]) < 1e-5)
    nav = " ".join('<a href="#%s">%s</a>' % (html.escape(d["mdir"]), html.escape(d["physics"].split("(")[0].strip()))
                   for d in mats)
    summ = datauri(os.path.join(RA, "results", "figures", "table_sweep_summary.png"))
    sections = "\n".join(section(d) for d in mats)
    page = f"""<title>Residual-method validation — 20 UMATs, per material</title>
<style>
:root {{
  --bg:#f5f7fa; --panel:#ffffff; --ink:#141d2b; --muted:#5a6675; --line:#e3e8ef;
  --accent:#25538f; --accent2:#3b7dd8; --code:#f2f5f9; --codeink:#22303f;
  --ok:#2e8b57; --warn:#a9791f; --bad:#c0392b;
  --jk:#25538f; --js:#2e7d4f; --jn:#a9791f; --jb:#8e44ad;
}}
@media (prefers-color-scheme:dark) {{
  :root {{ --bg:#0d131b; --panel:#151e29; --ink:#e6edf5; --muted:#93a2b5; --line:#24303e;
    --accent:#6aa3ff; --accent2:#8dbaff; --code:#0c1826; --codeink:#c7d6e6;
    --ok:#5cc98a; --warn:#e0b25a; --bad:#e8776b; --jk:#7bb0ff; --js:#6fd39a; --jn:#e0b25a; --jb:#c99bff; }}
}}
:root[data-theme=dark] {{ --bg:#0d131b; --panel:#151e29; --ink:#e6edf5; --muted:#93a2b5; --line:#24303e;
  --accent:#6aa3ff; --accent2:#8dbaff; --code:#0c1826; --codeink:#c7d6e6; --ok:#5cc98a; --warn:#e0b25a; --bad:#e8776b;
  --jk:#7bb0ff; --js:#6fd39a; --jn:#e0b25a; --jb:#c99bff; }}
:root[data-theme=light] {{ --bg:#f5f7fa; --panel:#ffffff; --ink:#141d2b; --muted:#5a6675; --line:#e3e8ef;
  --accent:#25538f; --accent2:#3b7dd8; --code:#f2f5f9; --codeink:#22303f; --ok:#2e8b57; --warn:#a9791f; --bad:#c0392b;
  --jk:#25538f; --js:#2e7d4f; --jn:#a9791f; --jb:#8e44ad; }}
*{{box-sizing:border-box}}
body{{margin:0;background:var(--bg);color:var(--ink);line-height:1.55;
  font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
  font-feature-settings:"tnum";}}
.mono,.num,.code{{font-family:ui-monospace,"SF Mono",Menlo,Consolas,monospace;font-variant-numeric:tabular-nums}}
header.top{{background:linear-gradient(180deg,var(--panel),transparent);border-bottom:1px solid var(--line);
  padding:2.2rem 1.4rem 1.2rem}}
.wrap{{max-width:76rem;margin:0 auto;padding:0 1.2rem}}
h1{{font-size:1.75rem;margin:0 0 .3rem;letter-spacing:-.02em}}
.lede{{color:var(--muted);max-width:56rem;margin:.2rem 0 1rem}}
.stats{{display:flex;gap:1.4rem;flex-wrap:wrap;margin:.6rem 0}}
.stat b{{font-size:1.5rem;color:var(--accent);font-variant-numeric:tabular-nums}}
.stat span{{color:var(--muted);font-size:.82rem;display:block}}
nav.idx{{position:sticky;top:0;z-index:5;background:color-mix(in srgb,var(--bg) 88%,transparent);
  backdrop-filter:blur(8px);border-bottom:1px solid var(--line);padding:.5rem 0}}
nav.idx .wrap{{display:flex;gap:.4rem .8rem;flex-wrap:wrap;font-size:.78rem}}
nav.idx a{{color:var(--muted);text-decoration:none;white-space:nowrap}}
nav.idx a:hover{{color:var(--accent)}}
section.mat{{padding:2rem 0;border-top:1px solid var(--line)}}
.mathead{{display:flex;justify-content:space-between;align-items:baseline;gap:1rem;flex-wrap:wrap}}
h2{{font-size:1.25rem;margin:0;letter-spacing:-.01em}}
.badges{{display:flex;gap:.4rem;flex-wrap:wrap}}
.badge{{font-size:.72rem;font-weight:600;padding:.16rem .5rem;border-radius:999px;border:1px solid var(--line);
  color:var(--muted);background:var(--panel);font-variant-numeric:tabular-nums}}
.badge.ok{{color:var(--ok);border-color:color-mix(in srgb,var(--ok) 40%,var(--line))}}
.badge.warn{{color:var(--warn);border-color:color-mix(in srgb,var(--warn) 40%,var(--line))}}
.badge.bad{{color:#fff;background:var(--bad);border-color:var(--bad)}}
.badge.pd{{color:var(--accent2)}} .badge.el{{color:var(--muted)}}
.props{{color:var(--muted);font-size:.85rem;margin:.5rem 0 1rem}}
.grid2{{display:grid;gap:1rem;grid-template-columns:1fr 1fr}}
@media(max-width:820px){{.grid2{{grid-template-columns:1fr}}}}
.card{{background:var(--panel);border:1px solid var(--line);border-radius:12px;overflow:hidden}}
.ct{{font-size:.82rem;font-weight:600;padding:.55rem .8rem;border-bottom:1px solid var(--line);
  display:flex;justify-content:space-between;align-items:center;gap:.5rem}}
.tag{{font-family:ui-monospace,monospace;font-size:.68rem;color:var(--muted);font-weight:500}}
pre.code{{margin:0;padding:.7rem .8rem;background:var(--code);color:var(--codeink);overflow-x:auto;
  font-size:.74rem;line-height:1.5;max-height:22rem}}
.jk{{color:var(--jk)}} .js{{color:var(--js)}} .jn{{color:var(--jn)}} .jb{{color:var(--jb)}}
table.err{{width:100%;border-collapse:collapse;font-size:.8rem}}
table.err th,table.err td{{text-align:left;padding:.3rem .8rem;border-top:1px solid var(--line)}}
table.err th{{font-size:.68rem;text-transform:uppercase;letter-spacing:.03em;color:var(--muted);font-weight:600}}
td.num,th.num{{text-align:right;font-family:ui-monospace,monospace}}
td.ok{{color:var(--ok)}} td.warn{{color:var(--warn)}} td.bad{{color:var(--bad);font-weight:700}}
.plots{{display:grid;gap:1rem;grid-template-columns:1fr;margin-top:1.1rem}}
.plots figure{{margin:0;background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:.6rem}}
.plots img{{width:100%;height:auto;display:block;border-radius:6px}}
.plots figcaption{{font-size:.76rem;color:var(--muted);padding:.4rem .2rem 0}}
.summ{{margin:1.4rem 0}}
.summ img{{width:100%;max-width:60rem;height:auto;border:1px solid var(--line);border-radius:12px;background:var(--panel)}}
footer{{color:var(--muted);font-size:.8rem;padding:2.5rem 0 4rem;border-top:1px solid var(--line);margin-top:1rem}}
</style>

<header class="top"><div class="wrap">
  <h1>Residual-method validation &mdash; 20 UMATs, per material</h1>
  <p class="lede">For every material model: the Program&nbsp;1 transform contract, a Program&nbsp;2 residual request,
  the weighted &sigma;<sub>vM</sub> parameter sensitivities, OTI-vs-finite-difference comparisons for <em>both</em>
  programs across the loading path, and a per-parameter error table. Everything is embedded &mdash; this page is
  self-contained.</p>
  <div class="stats">
    <div class="stat"><b>{len(mats)}</b><span>material models</span></div>
    <div class="stat"><b>{sum(len(d['params']) for d in mats)}</b><span>parameters validated</span></div>
    <div class="stat"><b>{nboth}/{len(mats)}</b><span>exact to &lt;1e-5 on both programs</span></div>
    <div class="stat"><b>2</b><span>programs (transform + residual)</span></div>
  </div>
  <div class="summ"><img alt="20-material summary" src="{summ}"></div>
</div></header>
<nav class="idx"><div class="wrap">{nav}</div></nav>
<main class="wrap">
{sections}
</main>
<footer><div class="wrap">Program&nbsp;1 = OTI DSIGMA_DP vs central finite differences of the original UMAT.
Program&nbsp;2 = residual-method d(&sigma;<sub>vM</sub>)/dp vs full-analysis FD (single C3D8, controlled deformation).
Drucker&ndash;Prager and Perzyna are the conditioning-sensitive nonlinear return maps. The three real ICP entries
(ECO, ECL_TEMP, PCO) are self-contained adaptations of production UMATs. Reproducible via
<span class="mono">results/plot_one_umat.py</span> + <span class="mono">results/build_per_umat_report.py</span>.</div></footer>
"""
    out = os.path.join(RA, "results", "per_umat_report.html")
    open(out, "w").write(page)
    print("wrote", out, "(%d KB, %d materials)" % (len(page) // 1024, len(mats)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
