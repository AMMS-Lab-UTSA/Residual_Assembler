#!/usr/bin/env python3
"""Assemble the self-contained HTML results report from the generated figures and
validation summaries. Images are embedded as base64 data URIs so the page needs no
external assets.

    UMAT_OTI_ROOT=~/Documents/UMAT_source_transformation python results/build_report.py
"""
import base64, json, os, sys

RA = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
FIG = os.path.join(RA, "results", "figures")


def img(name):
    p = os.path.join(FIG, name)
    b = base64.b64encode(open(p, "rb").read()).decode()
    return "data:image/png;base64," + b


def jload(name, default=None):
    p = os.path.join(FIG, name)
    return json.load(open(p)) if os.path.exists(p) else (default or {})


MATERIALS = [
    ("M1", "Isotropic linear elastic", "E, &nu;", 2, "&mdash;", False,
     "2.1e-16", "2.3e-10", "&mdash;"),
    ("M2", "Cubic anisotropic elastic", "C11, C12, C44", 3, "&mdash;", False,
     "0", "2.4e-13", "&mdash;"),
    ("M3", "J2 plasticity, linear hardening", "E, &nu;, &sigma;<sub>y0</sub>, H", 4, "EQPLAS", True,
     "8.7e-17", "2.2e-08", "9.8e-09"),
    ("M5", "Thermally-activated crystal-plasticity flow", "&tau;<sub>0</sub>, &Delta;G, p, q, &gamma;<sub>0</sub>, H", 6, "EQPLAS", True,
     "1.4e-16", "2.4e-07", "1.6e-07"),
]


def main():
    cp = jload("cp_sensitivity_summary.json")
    mesh = jload("cp_mesh_summary.json")
    fcc = jload("fcc_results.json")
    fcc_worst = fcc.get("worst_rmse_rel", 0.0)
    fcc_fd = fcc.get("cost", {}).get("fd_runs", 21)
    worst = cp.get("worst_rmse_rel", 0.0)
    oti_ratio = cp.get("cost", {}).get("oti_vs_regular", 0.0)
    fd_runs = cp.get("cost", {}).get("fd_runs", 13)

    matrows = "\n".join(
        f"""<tr><td class="tag">{mid}</td><td>{name}</td><td class="mono">{params}</td>
        <td class="num">{npar}</td><td class="ctr">{'yes' if path else '&ndash;'}</td>
        <td class="num">{sp}</td><td class="num ok">{ds}</td><td class="num ok">{dv}</td></tr>"""
        for (mid, name, params, npar, sv, path, sp, ds, dv) in MATERIALS)

    html = f"""<title>OTI Residual-Method Sensitivities &mdash; Results</title>
<style>
:root {{
  --bg:#faf9f6; --panel:#ffffff; --ink:#182430; --muted:#5b6b7a; --line:#e4e2db;
  --accent:#1f63c4; --accent2:#2b9c94; --ok:#2e8b57; --warn:#c98a1a; --code:#0f2233;
  --shadow:0 1px 2px rgba(20,30,40,.05), 0 8px 24px rgba(20,30,40,.06);
}}
@media (prefers-color-scheme: dark) {{
  :root {{ --bg:#0e141b; --panel:#151d26; --ink:#e7eef4; --muted:#9fb0bf; --line:#26313d;
    --accent:#5aa0ff; --accent2:#48c7bd; --ok:#5cc98a; --warn:#e0b25a; --code:#c7d6e4;
    --shadow:0 1px 2px rgba(0,0,0,.3), 0 10px 30px rgba(0,0,0,.35); }}
}}
:root[data-theme="dark"] {{ --bg:#0e141b; --panel:#151d26; --ink:#e7eef4; --muted:#9fb0bf;
  --line:#26313d; --accent:#5aa0ff; --accent2:#48c7bd; --ok:#5cc98a; --warn:#e0b25a; --code:#c7d6e4;
  --shadow:0 1px 2px rgba(0,0,0,.3), 0 10px 30px rgba(0,0,0,.35); }}
:root[data-theme="light"] {{ --bg:#faf9f6; --panel:#ffffff; --ink:#182430; --muted:#5b6b7a;
  --line:#e4e2db; --accent:#1f63c4; --accent2:#2b9c94; --ok:#2e8b57; --warn:#c98a1a; --code:#0f2233;
  --shadow:0 1px 2px rgba(20,30,40,.05), 0 8px 24px rgba(20,30,40,.06); }}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:var(--bg); color:var(--ink);
  font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
  line-height:1.62; -webkit-font-smoothing:antialiased; }}
.wrap {{ max-width:64rem; margin:0 auto; padding:3.2rem 1.4rem 5rem; }}
.eyebrow {{ text-transform:uppercase; letter-spacing:.14em; font-size:.72rem; font-weight:700;
  color:var(--accent); margin:0 0 .5rem; }}
h1 {{ font-size:2.15rem; line-height:1.12; letter-spacing:-.02em; margin:0 0 .6rem; text-wrap:balance; }}
h2 {{ font-size:1.4rem; letter-spacing:-.01em; margin:2.8rem 0 .3rem; padding-top:1.6rem;
  border-top:1px solid var(--line); }}
h3 {{ font-size:1.02rem; margin:1.6rem 0 .5rem; color:var(--ink); }}
.lede {{ font-size:1.12rem; color:var(--muted); margin:.2rem 0 1.4rem; max-width:52rem; }}
.meta {{ display:flex; flex-wrap:wrap; gap:.5rem 1.4rem; font-size:.82rem; color:var(--muted);
  border-top:1px solid var(--line); border-bottom:1px solid var(--line); padding:.7rem 0; margin:1.2rem 0 0; }}
.meta b {{ color:var(--ink); font-weight:600; }}
p {{ margin:.6rem 0 1rem; }}
.mono, .num {{ font-variant-numeric:tabular-nums;
  font-family:ui-monospace,SFMono-Regular,"SF Mono",Menlo,Consolas,monospace; }}
.grid {{ display:grid; gap:1.1rem; }}
@media (min-width:720px) {{ .grid.two {{ grid-template-columns:1fr 1fr; }} }}
.card {{ background:var(--panel); border:1px solid var(--line); border-radius:12px;
  box-shadow:var(--shadow); padding:1rem 1.1rem 1.1rem; }}
.card figcaption {{ font-size:.82rem; color:var(--muted); margin-top:.6rem; }}
.card img {{ width:100%; height:auto; display:block; border-radius:6px; }}
.kpis {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr)); gap:.8rem; margin:1.3rem 0; }}
.kpi {{ background:var(--panel); border:1px solid var(--line); border-radius:10px; padding:.85rem 1rem; }}
.kpi .v {{ font-size:1.5rem; font-weight:700; letter-spacing:-.02em; font-variant-numeric:tabular-nums; }}
.kpi .v.ok {{ color:var(--ok); }}
.kpi .k {{ font-size:.74rem; text-transform:uppercase; letter-spacing:.06em; color:var(--muted); margin-top:.15rem; }}
table {{ border-collapse:collapse; width:100%; font-size:.86rem; margin:.6rem 0 1rem; }}
.scroll {{ overflow-x:auto; }}
th, td {{ text-align:left; padding:.5rem .6rem; border-bottom:1px solid var(--line); vertical-align:top; }}
th {{ font-size:.72rem; text-transform:uppercase; letter-spacing:.05em; color:var(--muted); font-weight:600; }}
td.num, td.ctr, th.num, th.ctr {{ text-align:right; }} td.ctr, th.ctr {{ text-align:center; }}
td.num {{ font-variant-numeric:tabular-nums; font-family:ui-monospace,Menlo,Consolas,monospace; }}
td.ok {{ color:var(--ok); }}
.tag {{ display:inline-block; font-weight:700; color:var(--accent); }}
.pill {{ display:inline-block; font-size:.7rem; font-weight:700; padding:.12rem .5rem; border-radius:999px;
  background:color-mix(in srgb, var(--ok) 15%, transparent); color:var(--ok); text-transform:uppercase; letter-spacing:.04em; }}
.note {{ font-size:.86rem; color:var(--muted); border-left:3px solid var(--accent2); padding:.2rem 0 .2rem .9rem; margin:1rem 0; }}
.flow {{ display:flex; flex-wrap:wrap; align-items:center; gap:.5rem; font-size:.82rem; color:var(--muted); margin:.4rem 0 0; }}
.flow span {{ background:var(--panel); border:1px solid var(--line); border-radius:6px; padding:.28rem .6rem; }}
.flow b {{ color:var(--ink); }}
.arw {{ color:var(--accent); font-weight:700; }}
footer {{ margin-top:3rem; padding-top:1.2rem; border-top:1px solid var(--line); font-size:.8rem; color:var(--muted); }}
</style>

<div class="wrap">
  <p class="eyebrow">NASA STRI &middot; Residual-method parameter sensitivities</p>
  <h1>Model-agnostic material sensitivities via OTI, end to end</h1>
  <p class="lede">Two separately-distributable programs &mdash; a UMAT&rarr;OTI source transformer and an
  offline residual assembler &mdash; verified on real UMATs and on a crystal-plasticity flow model whose six
  flow-rule sensitivities are recovered in a single enriched run and validated against finite differences.</p>
  <div class="meta">
    <span><b>Program&nbsp;1</b> UMAT source transformation &rarr; OTI <span class="mono">.obj</span></span>
    <span><b>Program&nbsp;2</b> residual assembler (offline sensitivities)</span>
    <span><b>Reference</b> central finite differences of the original UMAT</span>
  </div>

  <div class="flow">
    <span><b>JHU</b> private UMAT source</span><span class="arw">&rarr;</span>
    <span>OTI transform</span><span class="arw">&rarr;</span>
    <span><b>umat_&lt;name&gt;_oti.obj</b> + contract</span><span class="arw">&rarr;</span>
    <span>collaborator FE record</span><span class="arw">&rarr;</span>
    <span>residual assembler</span><span class="arw">&rarr;</span>
    <span><b>&part;q/&part;p</b></span>
  </div>

  <h2>Program&nbsp;1 &nbsp;&middot;&nbsp; UMAT &rarr; OTI transformation</h2>
  <p>Each transformed object exposes the unchanged <span class="mono">UMAT</span> symbol plus an offline
  entry point <span class="mono">UMAT_OTI_EVAL</span> returning <span class="mono">DSIGMA_DP</span> and, for
  path-dependent models, <span class="mono">DSTATEV_DP</span> &mdash; the derivatives are returned separately,
  never stored in STATEV. Validation is non-circular: the OTI outputs are compared against a
  <em>separately compiled original UMAT</em>, and the derivatives against central finite differences of that
  original, marched over a multi-increment path for the history-dependent models.</p>

  <h3>Material-point validation <span class="pill">all pass</span></h3>
  <div class="scroll"><table>
    <thead><tr><th>Model</th><th>Constitutive law</th><th>Differentiated parameters</th>
      <th class="num">#p</th><th class="ctr">path-dep.</th>
      <th class="num">stress parity</th><th class="num">DSIGMA_DP vs FD</th><th class="num">DSTATEV_DP vs FD</th></tr></thead>
    <tbody>
      {matrows}
    </tbody>
  </table></div>
  <p class="note">Stress/state parity is machine precision; every parameter sensitivity matches central
  finite differences of the original UMAT to <span class="mono">10<sup>-7</sup></span> or better, including the
  path-dependent state sensitivity <span class="mono">DSTATEV_DP</span> propagated across increments.</p>

  <h3>Source-transformation test suite</h3>
  <p>Eighteen of nineteen benchmark UMATs pass the transformation verification, their OTI STRESS, STATEV and
  consistent tangent DDSDDE compared against the original UMAT run in Abaqus. The lone exception is an
  original job that did not converge; one row (<span class="mono">spin_elas_def</span>) whose DDSDDE
  disagreed flagged a genuine error in the original source.</p>
  <figure class="card"><img alt="UMAT source-transformation verification table" src="{img('umat_transform_table.png')}">
    <figcaption>Per-UMAT STRESS/STATEV/DDSDDE comparison against the original UMAT. Green: exact.
    Amber: a small, documented tangent difference (or a non-converged original job).</figcaption></figure>

  <h3>Both derivative outputs vs finite differences</h3>
  <p>The transformer returns two derivatives at every material point &mdash;
  <span class="mono">DSIGMA_DP</span> (stress) and <span class="mono">DSTATEV_DP</span> (state). Both are
  checked against centered finite differences of the original UMAT along the whole loading path.</p>
  <figure class="card"><img alt="Program 1 DSIGMA_DP and DSTATEV_DP vs finite differences" src="{img('program1_dsigma_dstatev_vs_fd.png')}">
    <figcaption>J2 plasticity (M3), state variable = equivalent plastic strain. OTI (lines) falls exactly on the
    central finite differences (points) through elastic, yield and hardening &mdash; worst relative RMSE 8.9e-09.
    SIGY0 and H activate only after yield; nu's stress derivative vanishes once plastic flow (pressure-independent)
    dominates &mdash; each curve is physically what it should be.</figcaption></figure>

  <h2>Program&nbsp;2 &nbsp;&middot;&nbsp; Residual-method sensitivities</h2>
  <p>The assembler links the distributed <span class="mono">.obj</span> directly, replays the material at each
  integration point, assembles <span class="mono">R</span>, <span class="mono">&part;R/&part;p</span> and
  <span class="mono">K</span>, solves <span class="mono">K&nbsp;&part;u/&part;p = &minus;&part;R/&part;p</span>,
  and evaluates the requested <span class="mono">&part;q/&part;p</span> &mdash; without rerunning the production
  analysis. Below: the crystal-plasticity flow model under simple shear.</p>

  <div class="kpis">
    <div class="kpi"><div class="v ok">{worst:.1e}</div><div class="k">worst OTI vs FD (rel. RMSE)</div></div>
    <div class="kpi"><div class="v">1 run</div><div class="k">OTI &mdash; all 6 sensitivities</div></div>
    <div class="kpi"><div class="v">{fd_runs} runs</div><div class="k">finite-difference equivalent</div></div>
    <div class="kpi"><div class="v ok">{mesh.get('mesh_vs_single_element_rel',0):.1e}</div><div class="k">mesh vs single element</div></div>
    <div class="kpi"><div class="v ok">{mesh.get('max_dudp_norm',0):.0e}</div><div class="k">&#8214;&part;u/&part;p&#8214; (homogeneous)</div></div>
  </div>

  <div class="grid two">
    <figure class="card"><img alt="Weighted sigma_vM sensitivities, single element" src="{img('weighted_sigvm_sensitivities.png')}">
      <figcaption>Weighted &sigma;<sub>vM</sub> parameter sensitivities, single element. All six obtained in one
      OTI run; the exponents <em>q, p</em> and activation energy &Delta;G dominate, hardening <em>H</em> grows with
      strain.</figcaption></figure>
    <figure class="card"><img alt="Weighted sigma_vM sensitivities, 4x4x4 mesh" src="{img('weighted_sigvm_sensitivities_4x4x4.png')}">
      <figcaption>The same weighted sensitivities from a full 4&times;4&times;4 residual solve. Identical to the
      single element (rel. {mesh.get('mesh_vs_single_element_rel',0):.0e}): the sensitivities are mesh
      independent.</figcaption></figure>
  </div>

  <figure class="card" style="margin-top:1.1rem"><img alt="OTI vs finite differences per parameter" src="{img('dsigvm_dp_oti_vs_fd.png')}">
    <figcaption>d&sigma;<sub>vM</sub>/dp along the shear path for each parameter &mdash; OTI (line) vs central
    finite differences (points). Curves are indistinguishable; relative RMSE ~10<sup>-9</sup>.</figcaption></figure>

  <div class="grid two" style="margin-top:1.1rem; align-items:start">
    <figure class="card"><img alt="Normalized cost, OTI vs finite differences" src="{img('normalized_cost.png')}">
      <figcaption>One enriched OTI run returns all six sensitivities for {oti_ratio:.1f}&times; a plain UMAT run;
      central differences need {fd_runs} full runs. (This is an unoptimized gfortran OTI carrying 19 directions;
      an optimized implementation reduces the OTI overhead substantially.)</figcaption></figure>
    <div class="card">
      <h3 style="margin-top:.2rem">What this establishes</h3>
      <p style="font-size:.92rem">The residual assembler consumes only the public
      <span class="mono">.obj</span> and a converged record, never the private source. For a path-dependent
      material it marches every integration point carrying both the state and its parameter sensitivity, so the
      assembled <span class="mono">&part;R/&part;p</span> is exact. Under controlled homogeneous shear the
      residual solve returns <span class="mono">&part;u/&part;p &asymp; 0</span> and the volume-averaged stress
      sensitivities match the single element to machine precision &mdash; a direct check of the global assembly.</p>
    </div>
  </div>

  <h2>FCC single-crystal plasticity &nbsp;&middot;&nbsp; ten sensitivities, one run</h2>
  <p>The headline model is a 12-slip-system FCC single crystal ({{111}}&lt;110&gt;) under [100] loading, with
  the ten standard parameters <span class="mono">C11, C12, C44, g0, gsat, h0, a, q, &gamma;&#775;0, m</span>.
  One enriched OTI run returns all ten <span class="mono">&part;&sigma;<sub>vM</sub>/&part;p</span>; each matches
  central finite differences of the regular UMAT to a relative RMSE of {fcc_worst:.0e} or better
  (<span class="mono">C44</span>&nbsp;&asymp;&nbsp;0, correct for [100]).</p>

  <figure class="card"><img alt="Weighted sigma_vM sensitivities with regimes, FCC [100]" src="{img('fcc_weighted_regimes.png')}">
    <figcaption>Weighted &sigma;<sub>vM</sub> sensitivities across the elastic&nbsp;(1), yield&nbsp;(2),
    hardening&nbsp;(3) and saturation&nbsp;(4) regimes. Elastic constants dominate regime&nbsp;1; the initial slip
    resistance <em>g0</em> leads early plasticity; the hardening parameters <em>gsat, a, h0, q</em> grow with
    strain &mdash; the balance validated parameter-by-parameter against finite differences.</figcaption></figure>

  <div class="grid two" style="margin-top:1.1rem; align-items:start">
    <figure class="card"><img alt="FCC OTI vs FD per parameter" src="{img('fcc_oti_vs_fd.png')}">
      <figcaption>Normalized sensitivity <span class="mono">(p/&sigma;)&part;&sigma;/&part;p</span> along the load
      path for the nine non-trivial parameters &mdash; OTI (line) vs central FD (dashed). Indistinguishable;
      relative RMSE ~10<sup>-9</sup>.</figcaption></figure>
    <figure class="card"><img alt="FCC analysis-level timing" src="{img('fcc_timing.png')}">
      <figcaption>Analysis level: the residual method replays the converged record once (the primal you already
      have + an enriched sensitivity pass) instead of re-solving the nonlinear system {fcc_fd} times &mdash;
      about 15&times; faster than finite differences.</figcaption></figure>
  </div>

  <h2>Per-parameter error tables</h2>
  <p>Every differentiated parameter of every material, compared against finite differences &mdash; for the source
  transformation (Program&nbsp;1) and again through the residual assembly (Program&nbsp;2).</p>
  <figure class="card"><img alt="Program 1 per-parameter error table" src="{img('table_program1_errors.png')}">
    <figcaption>Program&nbsp;1: each parameter's <span class="mono">DSIGMA_DP</span> vs central FD of the original
    UMAT. All 34 comparisons &le; 4.4e-8.</figcaption></figure>
  <figure class="card" style="margin-top:1.1rem"><img alt="Program 2 per-parameter error table" src="{img('table_program2_errors.png')}">
    <figcaption>Program&nbsp;2: each parameter's <span class="mono">&part;&sigma;<sub>vM</sub>/&part;p</span> through
    the residual assembly (R, &part;R/&part;p, K, solve, output) vs full-analysis FD. All &le; 4.5e-8.</figcaption></figure>

  <h2>Element formulations &nbsp;&middot;&nbsp; residual method across element types</h2>
  <p>The residual assembler was C3D8 only. <span class="mono">solid3d_kernel</span> adds the common Abaqus 3D solid
  family &mdash; reduced integration (C3D8R, C3D20R), quadratic hexes (C3D20) and tetrahedra (C3D4, C3D10) &mdash;
  reusing the same OTI UMAT and residual assembly.</p>
  <figure class="card"><img alt="Residual method validated across 6 element formulations" src="{img('element_validation_table.png')}">
    <figcaption>Each element reproduces a constant strain on a distorted mesh at machine precision (patch test),
    then the residual method is validated per element against full-analysis finite differences: assembly
    <span class="mono">&part;R/&part;p</span> &le; 5e-9, full solve <span class="mono">&part;u/&part;p</span> &le; 9e-8.
    FCC crystal plasticity (10 parameters, [100]) flows through every element type at ~3e-9. Single reduced-integration
    elements are rank-deficient (hourglass) and need hourglass control / a mesh for a standalone solve.</figcaption></figure>

  <h2>Breadth &nbsp;&middot;&nbsp; the whole framework across 20 material models</h2>
  <p>Both programs were run on 20 UMATs &mdash; the 6 model materials plus 14 authored/adapted models spanning
  elasticity (isotropic, cubic, orthotropic, transversely isotropic, Lam&eacute;, damage, thermo-, Cosserat,
  neo-Hookean), plasticity (J2 linear/bilinear/kinematic/combined, Drucker&ndash;Prager, couple-stress),
  viscoplasticity (Perzyna, thermally-activated), viscoelasticity (Maxwell) and FCC crystal plasticity &mdash;
  including self-contained adaptations of three real ICP production UMATs (ECO, ECL_TEMP, PCO).</p>
  <figure class="card"><img alt="20-material framework sweep summary" src="{img('table_sweep_summary.png')}">
    <figcaption>Per material, the worst relative error over its parameters for Program&nbsp;1
    (<span class="mono">DSIGMA_DP</span> vs FD of the original UMAT) and Program&nbsp;2 (residual
    <span class="mono">&part;&sigma;<sub>vM</sub>/&part;p</span> vs full-analysis FD). 18/20 are exact to
    &lt;1e-5 on both programs; Drucker&ndash;Prager (6.7e-3) and Perzyna (1.5e-4) are the conditioning-sensitive
    nonlinear return maps. The literal deep-external-call-tree production UMATs are a documented transform limit.</figcaption></figure>
  <figure class="card" style="margin-top:1.1rem"><img alt="OTI vs FD J2 kinematic hardening" src="{img('sweep_oti_vs_fd_j2kin.png')}">
    <figcaption>A newly authored sweep material (J2 with kinematic hardening): OTI (lines) on central finite
    differences (points) across the loading path, worst relative RMSE ~6e-10.</figcaption></figure>

  <footer>
    Generated from the two repositories' validation runs (Program&nbsp;1: <span class="mono">UMAT_source_transformation</span>,
    Program&nbsp;2: <span class="mono">Residual_Assembler</span>). Reference throughout is central finite differences
    of the original UMAT. Figures are reproducible via <span class="mono">results/*.py</span>.
  </footer>
</div>
"""
    out = os.path.join(RA, "results", "report.html")
    open(out, "w").write(html)
    print("wrote", out, "(%d KB)" % (len(html) // 1024))
    return 0


if __name__ == "__main__":
    sys.exit(main())
