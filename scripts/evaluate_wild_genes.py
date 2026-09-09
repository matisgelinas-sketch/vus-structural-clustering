"""
Tests the already-trained Phase 2 model on genes that were REJECTED from
training (failed the 30/30 Pathogenic/Benign data bar) -- the opposite of
cherry-picking. Every gene here was excluded specifically because its
ClinVar data isn't well-curated/balanced, which makes this a genuinely
unbiased test of whether the model says anything useful in the wild,
rather than only on genes hand-selected for good statistics.

Genes already known to fail structure/numbering validation (MECP2,
COL4A5, CACNA1A, NF1, SOS1, FBN1, MUTYH, RPGR, PTCH1, STXBP1) are
skipped -- re-fetching them would just reproduce a known failure, not
test anything new.

Nothing here is saved to Genes/ or folded into training -- this is
inference-only, using the final model trained on the 24 curated genes.
"""
import sys
import io
import json
from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parent.parent  # project root
sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent / "phase2b"))
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import precision_recall_fscore_support

from ncbi_clinvar import fetch_clinvar_germline_missense, screen_and_fetch_gene, fetch_domain_features
from clustering_lib import cross_check_numbering, compute_distances_leave_one_out

OUT = (PROJECT_ROOT / "output/phase2")
LOW_PLDDT_THRESHOLD = 70
DIST_THRESHOLD_A = 6.0

FEATURES = [
    'dist3d_A', 'seqdist_nearest3d', 'seq_over_3d_ratio', 'seq_minus_3d_diff',
    'dist3d_cb_A', 'seqdist_cb',
    'plddt', 'n_pathogenic_within_threshold', 'disordered_region',
]

# gene, accession, (n_pathogenic, n_benign) already known from expansion_log.json
# -- picked for >=10 of each class, diverse disease categories, and NOT
# already known to fail structure validation.
TEST_GENES = [
    ("FLNC", "Q14315"), ("SCN5A", "Q14524"), ("ATP7B", "P35670"), ("BAP1", "Q92560"),
    ("RAF1", "P04049"), ("TSC1", "Q92574"), ("DSP", "P15924"), ("APC", "P25054"),
    ("CDH1", "P12830"), ("WT1", "P19544"), ("MYH7", "P12883"), ("GAA", "P10253"),
    ("SERPINA1", "P01009"), ("F9", "P00740"), ("CACNA1S", "Q13698"), ("COL3A1", "P02461"),
    ("POLG", "P54098"), ("RAD51C", "O43502"), ("CDKN2A", "P42771"), ("PKD2", "Q13563"),
]

# --- train the final model once, on the 24 curated genes (same recipe as score_vus_phase2.py) ---
train_df = pd.read_csv(OUT / "training_pooled.csv")
train_df['disordered_region'] = train_df['domain'].str.contains('Disordered', na=False).astype(int)
X = train_df[FEATURES].values
y = train_df['label'].values
scaler = StandardScaler().fit(X)
logreg = LogisticRegression(class_weight='balanced', max_iter=2000, random_state=0)
logreg.fit(scaler.transform(X), y)
rf = RandomForestClassifier(n_estimators=200, max_depth=4, class_weight='balanced', random_state=0)
rf.fit(X, y)
print(f"Trained final model on {len(train_df)} labeled variants from "
      f"{train_df['gene'].nunique()} curated genes.\n")

results = []
for gene, accession in TEST_GENES:
    print(f"--- {gene} ---")
    clinvar_df, dropped = fetch_clinvar_germline_missense(gene, gene)
    n_patho_raw = int((clinvar_df['bucket'] == 'Pathogenic').sum())
    n_benign_raw = int((clinvar_df['bucket'] == 'Benign').sum())

    struct_result = screen_and_fetch_gene(gene, gene, accession, None)
    if not struct_result.get("passed"):
        print(f"  SKIPPED -- structure validation failed: {struct_result.get('reason')}\n")
        results.append(dict(gene=gene, status="structure_failed", reason=struct_result.get("reason")))
        continue

    seq = struct_result["seq"]
    struct_df = struct_result["struct_df"]
    mismatches = cross_check_numbering(clinvar_df, seq, struct_df)
    if len(mismatches):
        print(f"  SKIPPED -- {len(mismatches)} numbering mismatches\n")
        results.append(dict(gene=gene, status="numbering_mismatch", n_mismatches=len(mismatches)))
        continue

    domain_feats = fetch_domain_features(accession)

    def annotate(pos, feats=domain_feats):
        for s, e, label in feats:
            if s <= pos <= e:
                return label
        return "Unannotated / linker"

    dist_df = compute_distances_leave_one_out(clinvar_df, struct_df, count_threshold_A=DIST_THRESHOLD_A)
    dist_df = dist_df.merge(struct_df[['position', 'plddt']], on='position', how='left')
    dist_df['domain'] = dist_df['position'].apply(annotate)
    dist_df['disordered_region'] = dist_df['domain'].str.contains('Disordered', na=False).astype(int)
    dist_df = dist_df.dropna(subset=['dist3d_A'])
    safe_dist3d = dist_df['dist3d_A'].replace(0, 0.01)
    dist_df['seq_over_3d_ratio'] = dist_df['seqdist_nearest3d'] / safe_dist3d
    dist_df['seq_minus_3d_diff'] = dist_df['seqdist_nearest3d'] - dist_df['dist3d_A']

    labeled = dist_df[dist_df['bucket'].isin(['Pathogenic', 'Benign'])].copy()
    labeled['label'] = (labeled['bucket'] == 'Pathogenic').astype(int)

    if labeled['label'].nunique() < 2 or (labeled['label'] == 0).sum() < 5 or (labeled['label'] == 1).sum() < 5:
        print(f"  SKIPPED -- too few labeled examples with valid structure coords "
              f"(P={int((labeled['label']==1).sum())} B={int((labeled['label']==0).sum())})\n")
        results.append(dict(gene=gene, status="insufficient_after_structure_join"))
        continue

    Xt = labeled[FEATURES].values
    yt = labeled['label'].values

    pred_lr = logreg.predict(scaler.transform(Xt))
    pred_rf = rf.predict(Xt)
    prec_lr, rec_lr, f1_lr, _ = precision_recall_fscore_support(yt, pred_lr, labels=[0, 1], zero_division=0)
    prec_rf, rec_rf, f1_rf, _ = precision_recall_fscore_support(yt, pred_rf, labels=[0, 1], zero_division=0)
    bal_acc_lr = (rec_lr[0] + rec_lr[1]) / 2
    bal_acc_rf = (rec_rf[0] + rec_rf[1]) / 2

    n_p, n_b = int((yt == 1).sum()), int((yt == 0).sum())
    print(f"  n=P{n_p}/B{n_b} (raw ClinVar counts were P{n_patho_raw}/B{n_benign_raw})")
    print(f"  RF:     balanced_acc={bal_acc_rf:.3f}  Patho recall={rec_rf[1]:.2f}  Benign recall={rec_rf[0]:.2f}")
    print(f"  LogReg: balanced_acc={bal_acc_lr:.3f}  Patho recall={rec_lr[1]:.2f}  Benign recall={rec_lr[0]:.2f}\n")

    results.append(dict(
        gene=gene, status="tested", n_pathogenic=n_p, n_benign=n_b,
        bal_acc_rf=bal_acc_rf, patho_recall_rf=rec_rf[1], benign_recall_rf=rec_rf[0],
        bal_acc_logreg=bal_acc_lr, patho_recall_logreg=rec_lr[1], benign_recall_logreg=rec_lr[0],
    ))

results_df = pd.DataFrame(results)
results_df.to_csv(OUT / "wild_gene_evaluation.csv", index=False)

tested = results_df[results_df['status'] == 'tested']
print("=" * 60)
print(f"Tested {len(tested)} / {len(TEST_GENES)} genes (rest skipped for data-quality reasons -- see CSV)")
if len(tested):
    print(f"\nMean balanced accuracy (RF):     {tested['bal_acc_rf'].mean():.3f}")
    print(f"Mean balanced accuracy (LogReg): {tested['bal_acc_logreg'].mean():.3f}")
    print(f"\nFor comparison, mean balanced accuracy on the 24 CURATED (bar-passing) genes was ~0.71-0.73 (RF).")
    print(f"\nPer-gene results:")
    print(tested[['gene', 'n_pathogenic', 'n_benign', 'bal_acc_rf', 'bal_acc_logreg']]
          .sort_values('bal_acc_rf', ascending=False).to_string(index=False))
print(f"\nSaved full results to {OUT / 'wild_gene_evaluation.csv'}")
