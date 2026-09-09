"""
Phase 2, task 1: pool labeled variants (+ features) across all 15 genes
(TP53, BRCA1, PTEN + the 12 sourced this round) into one table. Also
computes features for VUS rows (used later for scoring in task 6), but
training will only use Pathogenic/Benign labels.
"""
import sys
import glob
from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parent.parent  # project root
sys.path.insert(0, str(Path(__file__).parent))
from clustering_lib import extract_fasta, parse_clinvar, parse_structure, compute_distances_leave_one_out
from domains import annotate_domain
import pandas as pd

GENES_ROOT = (PROJECT_ROOT / "Genes")

# Original 3: raw ClinVar TSV export, parsed via parse_clinvar().
RAW_TSV_GENES = {
    "TP53": dict(seq_file="TP53 sequence.rtf", struct_file="AF-P04637-F1-model_v6.pdb",
                 clinvar_file="Clinvar data 2.tsv"),
    "BRCA1": dict(seq_file="BRCA1 sequence.rtf", struct_file="AF-P38398-F1-model_v6.pdb",
                  clinvar_file="Clinvar Data.tsv"),
    "PTEN": dict(seq_file="PTEN sequence.rtf", struct_file="AF-P60484-F1-model_v6.pdb",
                 clinvar_file="clinvar_download_20260817_143437.tsv"),
}
# Every gene sourced via NCBI E-utilities (auto-discovered from
# expansion_log.json -- adding a gene via gene_pipeline.py is enough to
# have it picked up here automatically, no edits needed): already-parsed
# CSV, same schema parse_clinvar() returns, loaded directly.
import json as _json
_log_path = (PROJECT_ROOT / "scripts/phase2b/expansion_log.json")
_log = _json.loads(_log_path.read_text()) if _log_path.exists() else {}
PREPARSED_CSV_GENES = sorted(g for g, e in _log.items() if e.get("structure_passed"))

LOW_PLDDT_THRESHOLD = 70
DIST_THRESHOLD_A = 6.0

OUT = (PROJECT_ROOT / "output/phase2")
OUT.mkdir(parents=True, exist_ok=True)

all_frames = []

for gene, cfg in RAW_TSV_GENES.items():
    base = GENES_ROOT / gene
    header, seq = extract_fasta(base / cfg["seq_file"])
    struct_df = parse_structure(base / cfg["struct_file"])
    clinvar_df, dropped = parse_clinvar(base / cfg["clinvar_file"], gene)

    dist_df = compute_distances_leave_one_out(clinvar_df, struct_df, count_threshold_A=DIST_THRESHOLD_A)
    dist_df = dist_df.merge(struct_df[['position', 'plddt']], on='position', how='left')
    dist_df['low_confidence'] = dist_df['plddt'] < LOW_PLDDT_THRESHOLD
    dist_df['domain'] = dist_df['position'].apply(lambda p: annotate_domain(gene, p))
    dist_df['gene'] = gene

    safe_dist3d = dist_df['dist3d_A'].replace(0, 0.01)
    dist_df['seq_over_3d_ratio'] = dist_df['seqdist_nearest3d'] / safe_dist3d
    dist_df['seq_minus_3d_diff'] = dist_df['seqdist_nearest3d'] - dist_df['dist3d_A']
    safe_dist3d_cb = dist_df['dist3d_cb_A'].replace(0, 0.01)
    dist_df['seq_over_3d_ratio_cb'] = dist_df['seqdist_cb'] / safe_dist3d_cb
    dist_df['seq_minus_3d_diff_cb'] = dist_df['seqdist_cb'] - dist_df['dist3d_cb_A']

    all_frames.append(dist_df)
    print(f"{gene}: {len(dist_df)} rows with features "
          f"(Pathogenic={len(dist_df[dist_df.bucket=='Pathogenic'])}, "
          f"Benign={len(dist_df[dist_df.bucket=='Benign'])}, "
          f"VUS={len(dist_df[dist_df.bucket=='VUS'])})")

for gene in PREPARSED_CSV_GENES:
    base = GENES_ROOT / gene
    fasta_path = glob.glob(str(base / "*_sequence.fasta"))[0]
    pdb_path = glob.glob(str(base / "*.pdb"))[0]
    csv_path = glob.glob(str(base / "*_clinvar_germline_missense.csv"))[0]

    header, seq = extract_fasta(fasta_path)
    struct_df = parse_structure(pdb_path)
    clinvar_df = pd.read_csv(csv_path)

    dist_df = compute_distances_leave_one_out(clinvar_df, struct_df, count_threshold_A=DIST_THRESHOLD_A)
    dist_df = dist_df.merge(struct_df[['position', 'plddt']], on='position', how='left')
    dist_df['low_confidence'] = dist_df['plddt'] < LOW_PLDDT_THRESHOLD
    dist_df['domain'] = dist_df['position'].apply(lambda p: annotate_domain(gene, p))
    dist_df['gene'] = gene

    # Phase 2 features: ratio and difference between sequence and 3D distance
    safe_dist3d = dist_df['dist3d_A'].replace(0, 0.01)
    dist_df['seq_over_3d_ratio'] = dist_df['seqdist_nearest3d'] / safe_dist3d
    dist_df['seq_minus_3d_diff'] = dist_df['seqdist_nearest3d'] - dist_df['dist3d_A']
    safe_dist3d_cb = dist_df['dist3d_cb_A'].replace(0, 0.01)
    dist_df['seq_over_3d_ratio_cb'] = dist_df['seqdist_cb'] / safe_dist3d_cb
    dist_df['seq_minus_3d_diff_cb'] = dist_df['seqdist_cb'] - dist_df['dist3d_cb_A']

    all_frames.append(dist_df)
    print(f"{gene}: {len(dist_df)} rows with features "
          f"(Pathogenic={len(dist_df[dist_df.bucket=='Pathogenic'])}, "
          f"Benign={len(dist_df[dist_df.bucket=='Benign'])}, "
          f"VUS={len(dist_df[dist_df.bucket=='VUS'])})")

pooled_all = pd.concat(all_frames, ignore_index=True)
pooled_all = pooled_all.dropna(subset=['dist3d_A'])  # drop rows with no structure coordinate

FEATURE_COLS = [
    'gene', 'variation', 'position', 'wt_aa', 'mt_aa', 'bucket', 'domain',
    'plddt', 'low_confidence', 'dist3d_A', 'seqdist_nearest3d',
    'seq_over_3d_ratio', 'seq_minus_3d_diff', 'n_pathogenic_within_threshold',
    'nearest_pathogenic_pos', 'dist3d_cb_A', 'seqdist_cb',
    'seq_over_3d_ratio_cb', 'seq_minus_3d_diff_cb',
    'n_pathogenic_within_threshold_cb', 'nearest_pathogenic_pos_cb',
    'condition',
]
pooled_all = pooled_all[FEATURE_COLS]
pooled_all.to_csv(OUT / "all_genes_features.csv", index=False)
print(f"\nSaved full pooled feature table (all buckets, all genes): "
      f"{OUT / 'all_genes_features.csv'} ({len(pooled_all)} rows)")

training = pooled_all[pooled_all['bucket'].isin(['Pathogenic', 'Benign'])].copy()
training['label'] = (training['bucket'] == 'Pathogenic').astype(int)
training.to_csv(OUT / "training_pooled.csv", index=False)
print(f"Saved training table (Pathogenic/Benign only): "
      f"{OUT / 'training_pooled.csv'} ({len(training)} rows)")

print("\n=== Phase 2, task 2: class balance (before training) ===")
print(training.groupby(['gene', 'bucket']).size().unstack(fill_value=0))
print()
print("Overall:")
print(training['bucket'].value_counts())
print(f"Pathogenic fraction: {training['label'].mean():.3f}")
