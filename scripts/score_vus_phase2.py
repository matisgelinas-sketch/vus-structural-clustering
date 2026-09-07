"""
Phase 2, task 6: train a final model on ALL pooled Pathogenic/Benign labels
(all 3 genes -- this is the "deployed" model, distinct from the
gene-held-out folds used only to estimate generalization performance) and
apply it to every gene's VUS. Compare the model's probability score against
Phase 1's geometric flag; report agreement/disagreement.
"""
import sys
from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parent.parent  # project root
sys.path.insert(0, str(Path(__file__).parent))
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

OUT = (PROJECT_ROOT / "output/phase2")
train_df = pd.read_csv(OUT / "training_pooled.csv")
train_df['disordered_region'] = train_df['domain'].str.contains('Disordered', na=False).astype(int)

FEATURES = [
    'dist3d_A', 'seqdist_nearest3d', 'seq_over_3d_ratio', 'seq_minus_3d_diff',
    'plddt', 'n_pathogenic_within_threshold', 'disordered_region',
]

X = train_df[FEATURES].values
y = train_df['label'].values
scaler = StandardScaler().fit(X)
final_model = LogisticRegression(class_weight='balanced', max_iter=2000, random_state=0)
final_model.fit(scaler.transform(X), y)

# --- score every gene's VUS ---
all_features = pd.read_csv(OUT / "all_genes_features.csv")
vus = all_features[all_features['bucket'] == 'VUS'].copy()
vus['disordered_region'] = vus['domain'].str.contains('Disordered', na=False).astype(int)
X_vus = vus[FEATURES].values
vus['phase2_pathogenic_probability'] = final_model.predict_proba(scaler.transform(X_vus))[:, 1]

# --- bring in Phase 1's geometric flag per gene (auto-discovered from the pooled data) ---
ALL_GENES = sorted(train_df['gene'].unique())
OUTPUT_ROOT = (PROJECT_ROOT / "output")
GENE_DIRS = {g: OUTPUT_ROOT / g / f"{g}_all_VUS_annotated.csv" for g in ALL_GENES}
flag_frames = []
for gene, path in GENE_DIRS.items():
    d = pd.read_csv(path)[['variation', 'flagged_candidate']]
    d['gene'] = gene
    flag_frames.append(d)
flags = pd.concat(flag_frames, ignore_index=True)

merged = vus.merge(flags, on=['gene', 'variation'], how='left')
merged['flagged_candidate'] = merged['flagged_candidate'].fillna(False)

# A model-side "positive" call needs a threshold -- use 0.5 as the
# standard default (not tuned/cherry-picked).
merged['model_flags_pathogenic_like'] = merged['phase2_pathogenic_probability'] >= 0.5

def agreement_label(row):
    g, m = row['flagged_candidate'], row['model_flags_pathogenic_like']
    if g and m:
        return 'agree: both flag'
    if not g and not m:
        return 'agree: neither flags'
    if g and not m:
        return 'disagree: geometric only'
    return 'disagree: model only'

merged['agreement'] = merged.apply(agreement_label, axis=1)

OUT_COLS = ['gene', 'variation', 'position', 'domain', 'plddt',
            'dist3d_A', 'seqdist_nearest3d', 'n_pathogenic_within_threshold',
            'flagged_candidate', 'phase2_pathogenic_probability',
            'model_flags_pathogenic_like', 'agreement', 'condition']
merged[OUT_COLS].sort_values('phase2_pathogenic_probability', ascending=False).to_csv(
    OUT / "vus_scored_phase1_vs_phase2.csv", index=False)

print("=== Task 6: Phase 1 geometric flag vs Phase 2 model probability, per gene ===\n")
for gene in ALL_GENES:
    sub = merged[merged['gene'] == gene]
    print(f"--- {gene} ({len(sub)} VUS scored) ---")
    print(sub['agreement'].value_counts())
    print()

print("Overall agreement breakdown:")
print(merged['agreement'].value_counts())
print()
print(f"Saved: {OUT / 'vus_scored_phase1_vs_phase2.csv'} ({len(merged)} rows)")

print("\n=== Disagreement worth extra scrutiny: geometric flag says candidate, model says low probability ===")
disagree_geo_only = merged[merged['agreement'] == 'disagree: geometric only'].sort_values('phase2_pathogenic_probability')
print(f"{len(disagree_geo_only)} cases")
if len(disagree_geo_only):
    print(disagree_geo_only[['gene', 'variation', 'phase2_pathogenic_probability', 'dist3d_A', 'plddt']].head(20).to_string(index=False))

print("\n=== Disagreement worth extra scrutiny: model says high probability, geometric flag does NOT flag ===")
disagree_model_only = merged[merged['agreement'] == 'disagree: model only'].sort_values('phase2_pathogenic_probability', ascending=False)
print(f"{len(disagree_model_only)} cases")
if len(disagree_model_only):
    print(disagree_model_only[['gene', 'variation', 'phase2_pathogenic_probability', 'dist3d_A', 'seqdist_nearest3d', 'plddt']].head(20).to_string(index=False))
