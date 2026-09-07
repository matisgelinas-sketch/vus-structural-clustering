"""
External benchmark: compares this project's structural-clustering model
against AlphaMissense (Google DeepMind's genome-wide missense pathogenicity
predictor) on the exact same variants.

AlphaMissense scores are fetched per-protein from AlphaFold DB's
"aa-substitutions" CSV (https://alphafold.ebi.ac.uk/files/AF-{accession}
-F1-aa-substitutions.csv), which covers every possible amino-acid
substitution at every position of the canonical AlphaFold model -- the
same canonical UniProt numbering this project's own structures use, so
matching is exact (wt_aa + position + mt_aa), no genomic-coordinate or
transcript-mapping ambiguity.

Three comparisons, all on data AlphaMissense was never involved in
producing:
  1. On labeled (Pathogenic/Benign) variants: how well does AlphaMissense's
     OWN classification match ClinVar truth, per gene and overall --
     directly comparable to this project's own balanced-accuracy numbers.
  2. Correlation between this project's model probability and
     AlphaMissense's pathogenicity score, across all VUS.
  3. Agreement/disagreement specifically on the project's own
     high-confidence candidate list (flagged by both the geometric rule
     and the trained model) -- the money number: does independent
     external validation back up the project's own top candidates?
"""
import sys
import json
import time
import urllib.request
import io
from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parent.parent  # project root

sys.path.insert(0, str(Path(__file__).parent))
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr

OUT = (PROJECT_ROOT / "output/phase2")
LOG_PATH = (PROJECT_ROOT / "scripts/phase2b/expansion_log.json")

# --- gene -> UniProt accession map ---
log = json.loads(LOG_PATH.read_text()) if LOG_PATH.exists() else {}
ACCESSIONS = {g: e["accession"] for g, e in log.items() if e.get("structure_passed")}
ACCESSIONS.update({"TP53": "P04637", "BRCA1": "P38398", "PTEN": "P60484"})

AA3TO1 = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Gln': 'Q',
    'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
    'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W',
    'Tyr': 'Y', 'Val': 'V',
}


def fetch_alphamissense(accession: str) -> dict:
    """Returns {(wt_aa, position, mt_aa): (am_score, am_class)} for one protein."""
    url = f"https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-aa-substitutions.csv"
    with urllib.request.urlopen(url, timeout=30) as r:
        text = r.read().decode()
    df = pd.read_csv(io.StringIO(text))
    out = {}
    for _, row in df.iterrows():
        pv = row['protein_variant']
        wt, mt = pv[0], pv[-1]
        pos = int(pv[1:-1])
        out[(wt, pos, mt)] = (row['am_pathogenicity'], row['am_class'])
    return out


# --- fetch AlphaMissense for every gene in the project ---
print(f"Fetching AlphaMissense per-protein scores for {len(ACCESSIONS)} genes...\n")
am_lookup = {}
for i, (gene, acc) in enumerate(sorted(ACCESSIONS.items()), 1):
    try:
        am_lookup[gene] = fetch_alphamissense(acc)
        print(f"[{i}/{len(ACCESSIONS)}] {gene} ({acc}): {len(am_lookup[gene])} substitutions")
    except Exception as e:
        print(f"[{i}/{len(ACCESSIONS)}] {gene} ({acc}): FAILED -- {e}")
    time.sleep(0.15)

# --- load the full pooled feature set (all buckets, all genes) ---
all_features = pd.read_csv(OUT / "all_genes_features.csv")


def lookup_am(row):
    d = am_lookup.get(row['gene'])
    if d is None:
        return pd.Series([np.nan, None])
    key = (row['wt_aa'], int(row['position']), row['mt_aa'])
    res = d.get(key)
    if res is None:
        return pd.Series([np.nan, None])
    return pd.Series(res)


all_features[['am_score', 'am_class']] = all_features.apply(lookup_am, axis=1)
coverage = all_features['am_score'].notna().mean()
print(f"\nAlphaMissense coverage: {coverage:.1%} of {len(all_features)} variants matched.\n")

all_features.to_csv(OUT / "all_genes_features_with_alphamissense.csv", index=False)

# ============================================================
# Comparison 1: labeled variants -- AlphaMissense's own accuracy on
# ClinVar truth, per gene, directly comparable to this project's numbers.
# ============================================================
print("=" * 70)
print("COMPARISON 1: AlphaMissense classification vs. ClinVar truth")
print("(am_class: Likely benign / Ambiguous / Likely pathogenic)")
print("=" * 70)

labeled = all_features[all_features['bucket'].isin(['Pathogenic', 'Benign'])].dropna(subset=['am_score']).copy()
labeled['label'] = (labeled['bucket'] == 'Pathogenic').astype(int)
# AlphaMissense's own 0.564 pathogenic / 0.34 benign thresholds (their paper's cutoffs);
# score >= 0.564 predicted pathogenic-like, for symmetry with our 0.5-threshold reporting.
labeled['am_pred'] = (labeled['am_score'] >= 0.564).astype(int)

per_gene_am = []
for gene, sub in labeled.groupby('gene'):
    tp = ((sub['label'] == 1) & (sub['am_pred'] == 1)).sum()
    tn = ((sub['label'] == 0) & (sub['am_pred'] == 0)).sum()
    n_patho = (sub['label'] == 1).sum()
    n_benign = (sub['label'] == 0).sum()
    patho_recall = tp / n_patho if n_patho else np.nan
    benign_recall = tn / n_benign if n_benign else np.nan
    bal_acc = np.nanmean([patho_recall, benign_recall])
    per_gene_am.append(dict(gene=gene, n_patho=n_patho, n_benign=n_benign,
                             am_patho_recall=patho_recall, am_benign_recall=benign_recall,
                             am_balanced_accuracy=bal_acc))

am_summary = pd.DataFrame(per_gene_am).sort_values('am_balanced_accuracy', ascending=False)
print(am_summary.to_string(index=False))
print(f"\nAlphaMissense mean balanced accuracy across {len(am_summary)} genes: "
      f"{am_summary['am_balanced_accuracy'].mean():.3f}")
am_summary.to_csv(OUT / "alphamissense_per_gene_accuracy.csv", index=False)

# ============================================================
# Comparison 2: correlation between this project's model probability
# and AlphaMissense score, on VUS (never labeled by either source's
# ground truth -- this is agreement between two independent predictors).
# ============================================================
print("\n" + "=" * 70)
print("COMPARISON 2: correlation with this project's model, on VUS")
print("=" * 70)

vus_scored = pd.read_csv(OUT / "vus_scored_phase1_vs_phase2.csv")
vus_scored['key'] = list(zip(vus_scored['gene'], vus_scored['variation']))
all_features['key'] = list(zip(all_features['gene'], all_features['variation']))
vus_merged = vus_scored.merge(
    all_features[['key', 'am_score', 'am_class']], on='key', how='left'
).dropna(subset=['am_score', 'phase2_pathogenic_probability'])

r_pearson, p_pearson = pearsonr(vus_merged['phase2_pathogenic_probability'], vus_merged['am_score'])
r_spearman, p_spearman = spearmanr(vus_merged['phase2_pathogenic_probability'], vus_merged['am_score'])
print(f"VUS compared: {len(vus_merged)}")
print(f"Pearson r  = {r_pearson:.3f} (p={p_pearson:.1e})")
print(f"Spearman r = {r_spearman:.3f} (p={p_spearman:.1e})")

# ============================================================
# Comparison 3: the money number -- does AlphaMissense back up this
# project's own highest-confidence candidates (flagged by BOTH the
# geometric rule and the trained model)?
# ============================================================
print("\n" + "=" * 70)
print("COMPARISON 3: external validation of the high-confidence candidate list")
print("=" * 70)

both_flagged = vus_merged[vus_merged['agreement'] == 'agree: both flag'].copy()
both_flagged['am_agrees'] = both_flagged['am_score'] >= 0.564
n_agree = both_flagged['am_agrees'].sum()
n_total = len(both_flagged)
print(f"Of the {n_total} VUS flagged by BOTH this project's geometric rule AND its "
      f"trained model,\nAlphaMissense independently rates {n_agree} ({n_agree/n_total:.1%}) "
      f"as pathogenic-like (score >= 0.564).")

both_flagged.to_csv(OUT / "high_confidence_candidates_vs_alphamissense.csv", index=False)

disagree = both_flagged[~both_flagged['am_agrees']].sort_values('phase2_pathogenic_probability', ascending=False)
print(f"\n{len(disagree)} candidates where this project's structural model is confident "
      f"but AlphaMissense is NOT -- worth flagging as genuinely divergent cases:")
print(disagree[['gene', 'variation', 'phase2_pathogenic_probability', 'am_score', 'am_class']]
      .head(15).to_string(index=False))

print(f"\nSaved: all_genes_features_with_alphamissense.csv, alphamissense_per_gene_accuracy.csv, "
      f"high_confidence_candidates_vs_alphamissense.csv")
