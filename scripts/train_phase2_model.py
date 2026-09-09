"""
Phase 2, tasks 3-5: train simple interpretable models (logistic regression,
small random forest) with GENE-based held-out validation (never random
splitting), report precision/recall + feature importance.

Feature design note: domain is encoded as a single gene-agnostic
'disordered_region' binary flag (derived from the curated domain label
containing "Disordered"), NOT as one-hot of the literal per-gene domain
name (e.g. "BRCT domain 1", "DNA-binding domain"). Those literal domain
names are almost entirely gene-unique, so one-hot-encoding them directly
would let the model key off gene identity rather than a generalizable
structural property -- which would undermine the entire point of
gene-held-out validation ("do held-out genes test whether the pattern
generalizes"). pLDDT (continuous) + disordered_region (binary) is the
gene-agnostic stand-in for "domain/region" called for in the spec.
"""
import sys
from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parent.parent  # project root
sys.path.insert(0, str(Path(__file__).parent))
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import precision_recall_fscore_support, confusion_matrix, classification_report

OUT = (PROJECT_ROOT / "output/phase2")
df = pd.read_csv(OUT / "training_pooled.csv")
df['disordered_region'] = df['domain'].str.contains('Disordered', na=False).astype(int)

EXCLUDE_PLDDT = '--no-plddt' in sys.argv
FEATURES = [
    'dist3d_A', 'seqdist_nearest3d', 'seq_over_3d_ratio', 'seq_minus_3d_diff',
    'dist3d_cb_A', 'seqdist_cb',
    'n_pathogenic_within_threshold', 'disordered_region',
]
if not EXCLUDE_PLDDT:
    FEATURES.append('plddt')
SUFFIX = '_no_plddt' if EXCLUDE_PLDDT else ''
print(f"Features used: {FEATURES}\n")

# Auto-discovered: the 3 original genes + every gene that's passed
# gene_pipeline.py's validation, per training_pooled.csv's own gene column
# (which build_phase2_dataset.py already derives from expansion_log.json).
GENES = sorted(df['gene'].unique())

print("=== Task 3-4: gene-held-out cross-validation ===\n")

results = {'logreg': [], 'rf': []}
coef_records = []
importance_records = []

for held_out in GENES:
    train_df = df[df['gene'] != held_out]
    test_df = df[df['gene'] == held_out]

    X_train, y_train = train_df[FEATURES].values, train_df['label'].values
    X_test, y_test = test_df[FEATURES].values, test_df['label'].values

    scaler = StandardScaler().fit(X_train)
    X_train_s, X_test_s = scaler.transform(X_train), scaler.transform(X_test)

    logreg = LogisticRegression(class_weight='balanced', max_iter=2000, random_state=0)
    logreg.fit(X_train_s, y_train)
    pred_lr = logreg.predict(X_test_s)

    rf = RandomForestClassifier(n_estimators=200, max_depth=4, class_weight='balanced', random_state=0)
    rf.fit(X_train, y_train)  # tree models don't need scaling
    pred_rf = rf.predict(X_test)

    n_benign_test = int((y_test == 0).sum())
    n_patho_test = int((y_test == 1).sum())
    caveat = ""
    if n_benign_test < 10:
        caveat = f"  !! Benign class n={n_benign_test} in this test fold -- metrics below are NOT statistically meaningful for that class."

    print(f"--- Held-out gene: {held_out} (test set: {n_patho_test} Pathogenic, {n_benign_test} Benign) ---")
    if caveat:
        print(caveat)

    for name, pred in [('LogisticRegression', pred_lr), ('RandomForest', pred_rf)]:
        prec, rec, f1, support = precision_recall_fscore_support(
            y_test, pred, labels=[0, 1], zero_division=0)
        print(f"  {name}: Benign  precision={prec[0]:.2f} recall={rec[0]:.2f} f1={f1[0]:.2f} (n={support[0]})")
        print(f"  {name}: Pathogenic precision={prec[1]:.2f} recall={rec[1]:.2f} f1={f1[1]:.2f} (n={support[1]})")
        cm = confusion_matrix(y_test, pred, labels=[0, 1])
        print(f"  {name}: confusion matrix [[TN FP],[FN TP]] =\n{cm}")
        key = 'logreg' if name == 'LogisticRegression' else 'rf'
        results[key].append(dict(held_out=held_out, precision_benign=prec[0], recall_benign=rec[0],
                                  precision_patho=prec[1], recall_patho=rec[1],
                                  n_benign_test=n_benign_test, n_patho_test=n_patho_test))
    print()

    coef_records.append(pd.Series(logreg.coef_[0], index=FEATURES, name=held_out))
    importance_records.append(pd.Series(rf.feature_importances_, index=FEATURES, name=held_out))

for key in ['logreg', 'rf']:
    res_df = pd.DataFrame(results[key])
    res_df['balanced_accuracy'] = (res_df['recall_benign'] + res_df['recall_patho']) / 2
    res_df.to_csv(OUT / f"cv_results_{key}{SUFFIX}.csv", index=False)
    reliable = res_df[res_df['n_benign_test'] >= 10]
    print(f"\n{key}: mean balanced accuracy across all {len(res_df)} folds = {res_df['balanced_accuracy'].mean():.3f}")
    print(f"{key}: mean balanced accuracy across the {len(reliable)} folds with a statistically "
          f"usable Benign class (n>=10) = {reliable['balanced_accuracy'].mean():.3f}")

print(f"\n=== Task 5: feature importance (averaged across the {len(GENES)} held-out folds) ===\n")
coef_df = pd.DataFrame(coef_records)
imp_df = pd.DataFrame(importance_records)

print("Logistic regression standardized coefficients (sign = direction, magnitude = strength):")
print(coef_df.mean().sort_values(key=abs, ascending=False))
print()
print("Random forest feature importances (mean decrease in impurity):")
print(imp_df.mean().sort_values(ascending=False))

coef_df.to_csv(OUT / f"logreg_coefficients_per_fold{SUFFIX}.csv")
imp_df.to_csv(OUT / f"rf_importances_per_fold{SUFFIX}.csv")

print("\n=== Dataset size caveat ===")
print(f"Total labeled training examples: {len(df)} "
      f"({int((df['label']==1).sum())} Pathogenic, {int((df['label']==0).sum())} Benign "
      f"across {len(GENES)} genes).")
print("For comparison, AlphaMissense was trained on ~71M variants genome-wide using deep "
      f"learning on evolutionary + structural features; this dataset is still roughly 4 orders "
      f"of magnitude smaller and covers only {len(GENES)} genes. Results here characterize "
      "whether a simple structural-clustering signal carries ANY generalizable predictive "
      "information across genes -- they are not competitive with, and are not intended to "
      "compete with, genome-scale variant effect predictors.")
