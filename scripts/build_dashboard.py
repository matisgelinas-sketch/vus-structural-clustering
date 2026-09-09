"""
Builds the standalone HTML performance dashboard for the VUS structural
clustering project (Phase 1 + Phase 2, 24 genes, dual Ca/Cb distance).
Reads dashboard_data.json, embeds fonts as base64, writes a
self-contained artifact-ready HTML file.
"""
import json
from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parent.parent  # project root

DATA = json.loads((PROJECT_ROOT / "output/phase2/dashboard_data.json").read_text())
FONT_DIR = Path("/tmp/fonts_dashboard")
SANS_B64 = (FONT_DIR / "sans.b64").read_text()
MONO400_B64 = (FONT_DIR / "mono400.b64").read_text()
MONO500_B64 = (FONT_DIR / "mono500.b64").read_text()

# ---- derived numbers for KPI row ----
n_genes = len(DATA["panelA"])
total_variants = sum(r["pathogenic"] + r["benign"] + r["vus"] for r in DATA["panelC"])
total_candidates = sum(r["both"] for r in DATA["panelD"])
rf_bal_acc = [r["rf"] for r in DATA["panelA"]]
mean_rf_bal_acc = sum(rf_bal_acc) / len(rf_bal_acc)

# ---- Panel A: paired bars, sorted by RF descending ----
panelA_sorted = sorted(DATA["panelA"], key=lambda r: -r["rf"])
maxA = 1.0  # balanced accuracy scale 0-1, fixed axis

def bar_row_pair(r):
    gene = r["gene"]
    lr_pct = r["logreg"] * 100
    rf_pct = r["rf"] * 100
    reliable = r["reliable"]
    dim = "" if reliable else " unreliable"
    note = "" if reliable else '<span class="tag-warn">n=3 Benign</span>'
    return f'''
    <div class="pair-row{dim}">
      <div class="pair-label"><span class="gene-name">{gene}</span>{note}</div>
      <div class="pair-bars">
        <div class="pair-bar-track">
          <div class="pair-bar s-logreg" style="width:{lr_pct:.1f}%"></div>
          <span class="pair-val">{r['logreg']:.2f}</span>
        </div>
        <div class="pair-bar-track">
          <div class="pair-bar s-rf" style="width:{rf_pct:.1f}%"></div>
          <span class="pair-val">{r['rf']:.2f}</span>
        </div>
      </div>
    </div>'''

panelA_html = "\n".join(bar_row_pair(r) for r in panelA_sorted)

# ---- Panel B: dumbbell, feature importance shift ----
FEATURE_LABELS = {
    "dist3d_A": "Backbone (Ca) distance to nearest pathogenic residue",
    "dist3d_cb_A": "Side-chain (Cb) distance to nearest pathogenic residue",
    "n_pathogenic_within_threshold": "Pathogenic residues within 6Å",
    "plddt": "pLDDT (structural confidence)",
    "seqdist_nearest3d": "Sequence distance to nearest pathogenic residue (Ca)",
    "seqdist_cb": "Sequence distance to nearest pathogenic residue (Cb)",
    "seq_minus_3d_diff": "Sequence − 3D distance",
    "seq_over_3d_ratio": "Sequence / 3D distance ratio",
    "disordered_region": "Disordered region flag",
}
feat_order = sorted(DATA["panelB_after"].keys(), key=lambda k: -DATA["panelB_after"][k])
maxB = max(max(DATA["panelB_after"].values()), max(DATA["panelB_before"].values()))

def dumbbell_row(feat):
    before = DATA["panelB_before"][feat]
    after = DATA["panelB_after"][feat]
    b_pct = before / maxB * 100
    a_pct = after / maxB * 100
    lo, hi = (b_pct, a_pct) if b_pct <= a_pct else (a_pct, b_pct)
    is_clustering = feat in ("dist3d_A", "n_pathogenic_within_threshold")
    highlight = " highlight" if is_clustering else ""
    return f'''
    <div class="dumb-row{highlight}">
      <div class="dumb-label">{FEATURE_LABELS[feat]}</div>
      <div class="dumb-track">
        <div class="dumb-connector" style="left:{lo:.1f}%; width:{(hi-lo):.1f}%"></div>
        <div class="dumb-dot before" style="left:{b_pct:.1f}%"></div>
        <div class="dumb-dot after" style="left:{a_pct:.1f}%"></div>
      </div>
    </div>'''

panelB_html = "\n".join(dumbbell_row(f) for f in feat_order)

# ---- Panel C: stacked bar, Pathogenic/Benign per gene, VUS as annotation ----
panelC_sorted = sorted(DATA["panelC"], key=lambda r: -(r["pathogenic"] + r["benign"]))
maxC = max(r["pathogenic"] + r["benign"] for r in DATA["panelC"])

def stacked_row(r):
    total = r["pathogenic"] + r["benign"]
    p_pct = r["pathogenic"] / maxC * 100
    b_pct = r["benign"] / maxC * 100
    return f'''
    <div class="stack-row">
      <div class="stack-label"><span class="gene-name">{r['gene']}</span></div>
      <div class="stack-track">
        <div class="stack-seg s-patho" style="width:{p_pct:.1f}%"></div>
        <div class="stack-seg s-benign" style="width:{b_pct:.1f}%"></div>
      </div>
      <div class="stack-meta">{r['pathogenic']}P / {r['benign']}B <span class="muted">+{r['vus']:,} VUS</span></div>
    </div>'''

panelC_html = "\n".join(stacked_row(r) for r in panelC_sorted)

# ---- Panel D: Phase 1 flagged rate per gene ----
panelD_sorted = sorted(DATA["panelE"], key=lambda r: -r["pct_flagged"])
maxD = max(r["pct_flagged"] for r in DATA["panelE"])

def single_bar_row(r, max_val, unit_fmt):
    pct = r["pct_flagged"] / max_val * 100
    return f'''
    <div class="single-row">
      <div class="single-label"><span class="gene-name">{r['gene']}</span></div>
      <div class="single-track">
        <div class="single-bar" style="width:{pct:.1f}%"></div>
      </div>
      <div class="single-val">{unit_fmt(r)}</div>
    </div>'''

panelD_html = "\n".join(single_bar_row(r, maxD, lambda r: f"{r['pct_flagged']:.1f}% <span class=\"muted\">({r['n_flagged']})</span>") for r in panelD_sorted)

# ---- Panel E: cross-validated high-confidence candidates per gene ----
panelE_sorted = sorted(DATA["panelD"], key=lambda r: -r["both"])
maxE = max(r["both"] for r in DATA["panelD"])

def candidate_row(r):
    pct = r["both"] / maxE * 100 if maxE else 0
    return f'''
    <div class="single-row">
      <div class="single-label"><span class="gene-name">{r['gene']}</span></div>
      <div class="single-track">
        <div class="single-bar s-candidate" style="width:{pct:.1f}%"></div>
      </div>
      <div class="single-val">{r['both']}</div>
    </div>'''

panelE_html = "\n".join(candidate_row(r) for r in panelE_sorted)

HTML = f'''<title>VUS Clustering Signal</title>
<style>
@font-face {{
  font-family: 'Plex Sans';
  font-weight: 100 700;
  font-style: normal;
  font-display: swap;
  src: url(data:font/woff2;base64,{SANS_B64}) format('woff2');
}}
@font-face {{
  font-family: 'Plex Mono';
  font-weight: 400;
  font-style: normal;
  font-display: swap;
  src: url(data:font/woff2;base64,{MONO400_B64}) format('woff2');
}}
@font-face {{
  font-family: 'Plex Mono';
  font-weight: 500;
  font-style: normal;
  font-display: swap;
  src: url(data:font/woff2;base64,{MONO500_B64}) format('woff2');
}}

:root {{
  color-scheme: light;
  --paper:    #f7f8fa;
  --surface:  #ffffff;
  --ink:      #12161c;
  --ink-2:    #4b5563;
  --ink-3:    #8892a0;
  --line:     #e3e7ed;
  --accent:       #2a78d6;
  --accent-soft:  #eaf2fc;
  --s-logreg:   #2a78d6;
  --s-rf:       #eb6834;
  --s-patho:    #2a78d6;
  --s-benign:   #eb6834;
  --s-candidate:#1baf7a;
  --dumb-before:#9ec5f4;
  --dumb-after: #184f95;
  --good:     #0ca30c;
  --warning:  #fab219;
  --critical: #d03b3b;
}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
    color-scheme: dark;
    --paper:    #0e1116;
    --surface:  #171b21;
    --ink:      #f4f6f8;
    --ink-2:    #b9c0cc;
    --ink-3:    #737d8c;
    --line:     #2a303a;
    --accent:       #3987e5;
    --accent-soft:  rgba(57,135,229,0.16);
    --s-logreg:   #3987e5;
    --s-rf:       #d95926;
    --s-patho:    #3987e5;
    --s-benign:   #d95926;
    --s-candidate:#199e70;
    --dumb-before:#3a5a86;
    --dumb-after: #7fb3f7;
  }}
}}
:root[data-theme="dark"] {{
  color-scheme: dark;
  --paper:    #0e1116;
  --surface:  #171b21;
  --ink:      #f4f6f8;
  --ink-2:    #b9c0cc;
  --ink-3:    #737d8c;
  --line:     #2a303a;
  --accent:       #3987e5;
  --accent-soft:  rgba(57,135,229,0.16);
  --s-logreg:   #3987e5;
  --s-rf:       #d95926;
  --s-patho:    #3987e5;
  --s-benign:   #d95926;
  --s-candidate:#199e70;
  --dumb-before:#3a5a86;
  --dumb-after: #7fb3f7;
}}

* {{ box-sizing: border-box; }}
body {{
  background: var(--paper);
  color: var(--ink);
  font-family: 'Plex Sans', system-ui, -apple-system, sans-serif;
  font-size: 15px;
  line-height: 1.5;
  margin: 0;
  padding: 0;
}}
.wrap {{
  max-width: 980px;
  margin: 0 auto;
  padding: 48px 24px 96px;
  display: flex;
  flex-direction: column;
  gap: 56px;
}}

/* ---- header ---- */
.header {{
  display: flex;
  flex-direction: column;
  gap: 14px;
}}
.eyebrow {{
  font-family: 'Plex Mono', monospace;
  font-size: 12px;
  letter-spacing: 0.08em;
  text-transform: uppercase;
  color: var(--accent);
  font-weight: 500;
}}
h1 {{
  font-size: 32px;
  font-weight: 600;
  letter-spacing: -0.01em;
  margin: 0;
  text-wrap: balance;
}}
.dek {{
  color: var(--ink-2);
  font-size: 16px;
  max-width: 62ch;
  margin: 0;
}}
.caveat {{
  display: flex;
  gap: 10px;
  align-items: flex-start;
  border: 1px solid var(--line);
  background: var(--surface);
  border-radius: 6px;
  padding: 12px 16px;
  font-size: 13.5px;
  color: var(--ink-2);
  max-width: 62ch;
}}
.caveat b {{ color: var(--ink); }}

/* ---- KPI row ---- */
.kpi-row {{
  display: grid;
  grid-template-columns: repeat(4, 1fr);
  gap: 1px;
  background: var(--line);
  border: 1px solid var(--line);
  border-radius: 8px;
  overflow: hidden;
}}
.kpi {{
  background: var(--surface);
  padding: 20px 20px 18px;
  display: flex;
  flex-direction: column;
  gap: 6px;
}}
.kpi-label {{
  font-size: 12px;
  color: var(--ink-3);
  text-transform: uppercase;
  letter-spacing: 0.04em;
}}
.kpi-value {{
  font-family: 'Plex Mono', monospace;
  font-weight: 500;
  font-size: 28px;
  letter-spacing: -0.01em;
}}
.kpi-sub {{
  font-size: 12.5px;
  color: var(--ink-3);
}}

/* ---- sections ---- */
.section {{
  display: flex;
  flex-direction: column;
  gap: 4px;
  border-top: 1px solid var(--line);
  padding-top: 32px;
}}
.section-head {{
  display: flex;
  flex-direction: column;
  gap: 6px;
  margin-bottom: 20px;
}}
.section-num {{
  font-family: 'Plex Mono', monospace;
  font-size: 12px;
  color: var(--ink-3);
}}
h2 {{
  font-size: 19px;
  font-weight: 600;
  margin: 0;
  letter-spacing: -0.005em;
}}
.section-desc {{
  font-size: 14px;
  color: var(--ink-2);
  max-width: 68ch;
  margin: 0;
}}
.finding {{
  font-size: 13.5px;
  color: var(--ink-2);
  background: var(--accent-soft);
  border-radius: 6px;
  padding: 10px 14px;
  margin-top: 16px;
  max-width: 68ch;
}}
.finding b {{ color: var(--ink); }}

/* ---- legend ---- */
.legend {{
  display: flex;
  gap: 18px;
  flex-wrap: wrap;
  margin-bottom: 14px;
}}
.legend-item {{
  display: flex;
  align-items: center;
  gap: 7px;
  font-size: 12.5px;
  color: var(--ink-2);
}}
.swatch {{
  width: 10px; height: 10px;
  border-radius: 2px;
  flex: none;
}}

/* ---- paired bar (Panel A) ---- */
.pair-row {{
  display: grid;
  grid-template-columns: 92px 1fr;
  align-items: center;
  gap: 14px;
  padding: 6px 0;
}}
.pair-row.unreliable {{ opacity: 0.55; }}
.pair-label {{
  display: flex;
  align-items: baseline;
  gap: 6px;
  min-width: 0;
}}
.gene-name {{
  font-family: 'Plex Mono', monospace;
  font-weight: 500;
  font-size: 13px;
}}
.tag-warn {{
  font-size: 10px;
  color: var(--warning);
  border: 1px solid var(--warning);
  border-radius: 3px;
  padding: 1px 4px;
  white-space: nowrap;
}}
.pair-bars {{
  display: flex;
  flex-direction: column;
  gap: 3px;
}}
.pair-bar-track {{
  position: relative;
  height: 12px;
  background: var(--line);
  border-radius: 3px;
  display: flex;
  align-items: center;
}}
.pair-bar {{
  height: 100%;
  border-radius: 3px 0 3px 3px;
}}
.pair-bar.s-logreg {{ background: var(--s-logreg); }}
.pair-bar.s-rf {{ background: var(--s-rf); }}
.pair-val {{
  font-family: 'Plex Mono', monospace;
  font-size: 11px;
  color: var(--ink-2);
  margin-left: 6px;
}}
.baseline-note {{
  font-size: 11.5px;
  color: var(--ink-3);
  margin-top: 10px;
  font-family: 'Plex Mono', monospace;
}}

/* ---- dumbbell (Panel B) ---- */
.dumb-row {{
  display: grid;
  grid-template-columns: 260px 1fr;
  align-items: center;
  gap: 16px;
  padding: 10px 0;
  border-bottom: 1px solid var(--line);
}}
.dumb-row:last-child {{ border-bottom: none; }}
.dumb-row.highlight .dumb-label {{ color: var(--ink); font-weight: 500; }}
.dumb-label {{
  font-size: 13px;
  color: var(--ink-2);
}}
.dumb-track {{
  position: relative;
  height: 18px;
}}
.dumb-connector {{
  position: absolute;
  top: 50%;
  height: 2px;
  background: var(--line);
  transform: translateY(-50%);
}}
.dumb-dot {{
  position: absolute;
  top: 50%;
  width: 11px; height: 11px;
  border-radius: 50%;
  transform: translate(-50%, -50%);
  border: 2px solid var(--surface);
}}
.dumb-dot.before {{ background: var(--dumb-before); }}
.dumb-dot.after {{ background: var(--dumb-after); z-index: 1; }}

/* ---- stacked bar (Panel C) ---- */
.stack-row {{
  display: grid;
  grid-template-columns: 72px 1fr 170px;
  align-items: center;
  gap: 12px;
  padding: 5px 0;
}}
.stack-label {{ min-width: 0; }}
.stack-track {{
  display: flex;
  height: 14px;
  border-radius: 3px;
  overflow: hidden;
  background: var(--line);
}}
.stack-seg {{ height: 100%; }}
.stack-seg.s-patho {{ background: var(--s-patho); }}
.stack-seg.s-benign {{ background: var(--s-benign); }}
.stack-meta {{
  font-family: 'Plex Mono', monospace;
  font-size: 11.5px;
  color: var(--ink-2);
  white-space: nowrap;
}}
.stack-meta .muted {{ color: var(--ink-3); }}

/* ---- single bar (Panels D/E) ---- */
.single-row {{
  display: grid;
  grid-template-columns: 72px 1fr 110px;
  align-items: center;
  gap: 12px;
  padding: 5px 0;
}}
.single-label {{ min-width: 0; }}
.single-track {{
  height: 12px;
  border-radius: 3px;
  background: var(--line);
}}
.single-bar {{
  height: 100%;
  border-radius: 3px;
  background: var(--accent);
}}
.single-bar.s-candidate {{ background: var(--s-candidate); }}
.single-val {{
  font-family: 'Plex Mono', monospace;
  font-size: 12px;
  color: var(--ink-2);
  white-space: nowrap;
}}
.single-val .muted {{ color: var(--ink-3); }}

.two-col {{
  display: grid;
  grid-template-columns: 1fr 1fr;
  gap: 40px;
}}
@media (max-width: 760px) {{
  .two-col {{ grid-template-columns: 1fr; }}
  .kpi-row {{ grid-template-columns: repeat(2, 1fr); }}
}}

footer {{
  border-top: 1px solid var(--line);
  padding-top: 24px;
  font-size: 12px;
  color: var(--ink-3);
  font-family: 'Plex Mono', monospace;
}}
</style>

<div class="wrap">

  <div class="header">
    <div class="eyebrow">VUS structural clustering &mdash; performance dashboard</div>
    <h1>Does spatial clustering to known pathogenic residues predict variant pathogenicity?</h1>
    <p class="dek">Phase 1 flags candidate variants of uncertain significance (VUS) by 3D proximity to known pathogenic residues on AlphaFold structures. Phase 2 pools that signal across genes and tests whether it generalizes with gene-held-out cross-validation.</p>
    <div class="caveat">
      <span>&#9888;&#65039;</span>
      <span><b>Hypothesis-generation only.</b> Every number on this page describes a computational clustering signal, not a clinical classification. Nothing here is "likely pathogenic," "reclassified," or "confirmed" &mdash; the ceiling is <b>candidate for further investigation</b>.</span>
    </div>
  </div>

  <div class="kpi-row">
    <div class="kpi">
      <div class="kpi-label">Genes covered</div>
      <div class="kpi-value">{n_genes}</div>
      <div class="kpi-sub">up from 3 at project start</div>
    </div>
    <div class="kpi">
      <div class="kpi-label">Classified missense variants</div>
      <div class="kpi-value">{total_variants:,}</div>
      <div class="kpi-sub">germline, structure- &amp; numbering-validated</div>
    </div>
    <div class="kpi">
      <div class="kpi-label">Mean balanced accuracy</div>
      <div class="kpi-value">{mean_rf_bal_acc:.2f}</div>
      <div class="kpi-sub">random forest, {n_genes} gene-held-out folds</div>
    </div>
    <div class="kpi">
      <div class="kpi-label">Cross-validated candidates</div>
      <div class="kpi-value">{total_candidates}</div>
      <div class="kpi-sub">flagged by geometry <i>and</i> the model</div>
    </div>
  </div>

  <div class="section">
    <div class="section-head">
      <div class="section-num">01 &mdash; Phase 2, task 4</div>
      <h2>Gene-held-out validation: how well does it generalize?</h2>
      <p class="section-desc">Each gene is held out in turn, the model trained on the other {n_genes - 1}, then tested on the held-out gene. Balanced accuracy averages Pathogenic recall and Benign recall, so per-gene class imbalance doesn't distort the picture. Sorted by random forest score.</p>
    </div>
    <div class="legend">
      <div class="legend-item"><span class="swatch" style="background:var(--s-logreg)"></span>Logistic regression</div>
      <div class="legend-item"><span class="swatch" style="background:var(--s-rf)"></span>Random forest</div>
    </div>
    {panelA_html}
    <div class="baseline-note">0.50 = random baseline &nbsp;&middot;&nbsp; 1.00 = perfect &nbsp;&middot;&nbsp; PTEN dimmed: only 3 Benign variants exist for this gene, so its Benign-class recall isn't statistically meaningful</div>
    <div class="finding"><b>The spread is the finding.</b> TP53, LDLR and VHL generalize well (recall &gt; 0.9 for the majority class); MYBPC3, PMS2 and RB1 barely beat chance. A hypothesis-generation tool should report that range honestly rather than average it into one confident-sounding number.</div>
  </div>

  <div class="section">
    <div class="section-head">
      <div class="section-num">02 &mdash; Phase 2, task 5 + ablation</div>
      <h2>What the model actually learned, before and after the Cα/Cβ dual-atom fix</h2>
      <p class="section-desc">Random-forest feature importance, averaged across held-out folds. Light dot = the earlier 15-gene, Cα-only run (also affected by a nearest-neighbor masking bug, since fixed); dark dot = the current 24-gene, dual Cα/Cβ run. Sorted by current importance.</p>
    </div>
    <div class="legend">
      <div class="legend-item"><span class="swatch" style="background:var(--dumb-before);border-radius:50%"></span>Earlier run (15 genes, Cα-only)</div>
      <div class="legend-item"><span class="swatch" style="background:var(--dumb-after);border-radius:50%"></span>Current run ({n_genes} genes, Cα + Cβ)</div>
    </div>
    {panelB_html}
    <div class="finding"><b>Side-chain (Cβ) distance is now the single most important feature in both models</b> &mdash; ahead of backbone (Cα) distance, which itself outranks pLDDT. Removing pLDDT entirely from the model costs less than 0.001 balanced accuracy: the clustering signal is real and largely independent of structural confidence, not a proxy for it.</div>
  </div>

  <div class="section">
    <div class="section-head">
      <div class="section-num">03 &mdash; Phase 2, task 1&ndash;2</div>
      <h2>The training data: scale and class balance, gene by gene</h2>
      <p class="section-desc">Germline Pathogenic vs. Benign missense variants used for training (VUS count shown alongside, not stacked &mdash; VUS are scored afterward, not trained on). Sorted by labeled-variant count.</p>
    </div>
    <div class="legend">
      <div class="legend-item"><span class="swatch" style="background:var(--s-patho)"></span>Pathogenic / Likely pathogenic</div>
      <div class="legend-item"><span class="swatch" style="background:var(--s-benign)"></span>Benign / Likely benign</div>
    </div>
    {panelC_html}
    <div class="finding"><b>Balance varies by an order of magnitude across genes</b> &mdash; from LDLR's 18:1 Pathogenic skew to MSH2's roughly 1:2 Benign skew. Class weighting was recomputed independently within every training fold specifically to handle this, rather than relying on the pooled 62:38 average.</div>
  </div>

  <div class="two-col">
    <div class="section" style="border-top:none; padding-top:0;">
      <div class="section-head">
        <div class="section-num">04 &mdash; Phase 1</div>
        <h2>Geometric flagging rate</h2>
        <p class="section-desc">Share of each gene's VUS flagged by the geometric rule alone (&le;6&#8491; in 3D <i>and</i> &gt;10 residues apart in sequence &mdash; excludes trivial chain-adjacent cases).</p>
      </div>
      {panelD_html}
    </div>

    <div class="section" style="border-top:none; padding-top:0;">
      <div class="section-head">
        <div class="section-num">05 &mdash; Phase 2, task 6</div>
        <h2>Cross-validated candidates</h2>
        <p class="section-desc">VUS flagged by <i>both</i> the geometric rule and the statistical model &mdash; the strongest candidate list this project currently produces, {total_candidates} variants across {n_genes} genes.</p>
      </div>
      {panelE_html}
    </div>
  </div>

  <footer>
    Data: ClinVar (NCBI, germline classifications only) &middot; AlphaFold DB structures &middot; UniProt canonical sequences &middot; {n_genes}/125 ACMG SF v3.2 candidate genes cleared the 30 Pathogenic / 30 Benign data bar and structure validation. Variants with "Conflicting classifications of pathogenicity" or no classification are excluded throughout &mdash; for the original 3 genes this was done at the ClinVar search stage; for genes sourced via the NCBI API it happens in the same parsing step, same end result. 3D distance is computed by both backbone (Cα) and side-chain (Cβ) atom position; a VUS is flagged if either is within 6&#8491; of a known pathogenic residue and &gt;10 residues away in sequence.
  </footer>

</div>
'''

out_path = (PROJECT_ROOT / "output/phase2/dashboard.html")
out_path.write_text(HTML)
print(f"Wrote {out_path} ({len(HTML):,} chars)")
