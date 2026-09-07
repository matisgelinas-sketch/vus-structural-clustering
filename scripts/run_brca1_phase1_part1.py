import sys
from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parent.parent  # project root
sys.path.insert(0, str(Path(__file__).parent))
from clustering_lib import (
    extract_fasta, parse_clinvar, parse_structure,
    cross_check_numbering, compute_distances,
)

BASE = (PROJECT_ROOT / "Genes/BRCA1")
OUT = (PROJECT_ROOT / "output/BRCA1")
OUT.mkdir(parents=True, exist_ok=True)

GENE_TAG = "BRCA1"
LOW_PLDDT_THRESHOLD = 70

print("=== Task 2/3: sequence + structure ===")
header, seq = extract_fasta(BASE / "BRCA1 sequence.rtf")
print(f"UniProt header: {header}")
print(f"Sequence length: {len(seq)}")

struct_df = parse_structure(BASE / "AF-P38398-F1-model_v6.pdb")
print(f"Structure residues (Ca atoms): {len(struct_df)}")
if len(struct_df) != len(seq):
    print(f"!! LENGTH MISMATCH: structure={len(struct_df)} vs UniProt seq={len(seq)}")

print()
print("=== Task 1: ClinVar parsing ===")
clinvar_df, dropped = parse_clinvar(BASE / "Clinvar Data.tsv", GENE_TAG)
print(f"Rows kept after transcript/protein-change filter: {len(clinvar_df)}")
print(f"Dropped breakdown: {dropped}")
print()
bucket_counts = clinvar_df['bucket'].value_counts(dropna=False)
print("Classification bucket counts (germline only):")
print(bucket_counts)

print()
print("=== Task 3: cross-check ClinVar wt residue vs UniProt seq vs structure ===")
mismatches = cross_check_numbering(clinvar_df, seq, struct_df)
print(f"Mismatches found: {len(mismatches)}")
if len(mismatches):
    print(mismatches.head(30).to_string(index=False))
    mismatches.to_csv(OUT / "numbering_mismatches.csv", index=False)

print()
print(f"=== Task 4: low-confidence (pLDDT < {LOW_PLDDT_THRESHOLD}) residues among Pathogenic/VUS ===")
merged = clinvar_df.merge(struct_df[['position', 'plddt']], on='position', how='left')
low_conf = merged[(merged['plddt'] < LOW_PLDDT_THRESHOLD) & (merged['bucket'].isin(['Pathogenic', 'VUS']))]
denom = len(merged[merged['bucket'].isin(['Pathogenic','VUS'])])
print(f"Pathogenic/VUS residues in low-confidence regions: {len(low_conf)} / {denom}")
if len(low_conf):
    low_conf.to_csv(OUT / "low_confidence_flagged.csv", index=False)
    print(f"(full list saved; {len(low_conf)} rows)")

print()
print("=== Task 5: 3D + sequence distance for every VUS to nearest Pathogenic residue ===")
dist_df = compute_distances(clinvar_df, struct_df)
vus_dist = dist_df[dist_df['bucket'] == 'VUS'].dropna(subset=['dist3d_A'])
print(f"VUS with computable distances: {len(vus_dist)} / {len(dist_df[dist_df['bucket']=='VUS'])}")
print()
print("Distance distribution (Angstrom) summary:")
print(vus_dist['dist3d_A'].describe())
print()
print("Counts of VUS at various 3D-distance cutoffs (for threshold discussion):")
for cutoff in [5, 6, 7, 8, 9, 10, 12, 15]:
    n = (vus_dist['dist3d_A'] <= cutoff).sum()
    print(f"  <= {cutoff:>2} A: {n}")

dist_df.to_csv(OUT / "brca1_vus_distances_full.csv", index=False)
print()
print(f"Saved full distance table to {OUT / 'brca1_vus_distances_full.csv'}")
