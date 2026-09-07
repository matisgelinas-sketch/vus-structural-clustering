import sys
from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parent.parent  # project root
sys.path.insert(0, str(Path(__file__).parent))
from clustering_lib import extract_fasta, parse_clinvar, parse_structure, compute_distances
from domains import annotate_domain
import py3Dmol

BASE = (PROJECT_ROOT / "Genes/BRCA1")
OUT = (PROJECT_ROOT / "output/BRCA1")
GENE = "BRCA1"
GENE_TAG = "BRCA1"
LOW_PLDDT_THRESHOLD = 70
DIST_THRESHOLD_A = 6.0  # same methodology/threshold as TP53, per user instruction
MIN_SEQDIST_RESIDUES = 10  # confirmed with user: excludes trivial chain-neighbor "clustering"

# --- rebuild task 1-5 state ---
header, seq = extract_fasta(BASE / "BRCA1 sequence.rtf")
struct_df = parse_structure(BASE / "AF-P38398-F1-model_v6.pdb")
clinvar_df, dropped = parse_clinvar(BASE / "Clinvar Data.tsv", GENE_TAG)
dist_df = compute_distances(clinvar_df, struct_df)
dist_df = dist_df.merge(struct_df[['position', 'plddt']], on='position', how='left')
dist_df['low_confidence'] = dist_df['plddt'] < LOW_PLDDT_THRESHOLD
dist_df['domain'] = dist_df['position'].apply(lambda p: annotate_domain(GENE, p))

# --- Task 6: flag VUS within distance threshold, excluding trivial chain-adjacent cases ---
vus = dist_df[dist_df['bucket'] == 'VUS'].dropna(subset=['dist3d_A']).copy()
vus['flagged_candidate'] = (vus['dist3d_A'] <= DIST_THRESHOLD_A) & (vus['seqdist_nearest3d'] > MIN_SEQDIST_RESIDUES)
flagged = vus[vus['flagged_candidate']].copy()
flagged['seq_over_3d'] = flagged['seqdist_nearest3d'] / flagged['dist3d_A'].replace(0, 0.01)
flagged = flagged.sort_values('seq_over_3d', ascending=False)

print(f"=== Task 6: VUS flagged as candidates (3D distance <= {DIST_THRESHOLD_A} A AND sequence distance > {MIN_SEQDIST_RESIDUES} residues) ===")
print(f"{len(flagged)} / {len(vus)} VUS flagged")
print(f"  of which low-confidence (pLDDT<{LOW_PLDDT_THRESHOLD}): {flagged['low_confidence'].sum()}")

# --- Task 8: summary table ---
summary_cols = ['variation', 'position', 'wt_aa', 'mt_aa', 'domain', 'plddt',
                 'low_confidence', 'nearest_pathogenic_pos', 'dist3d_A',
                 'seqdist_nearest3d', 'seq_over_3d', 'condition']
summary = flagged[summary_cols].rename(columns={
    'dist3d_A': 'dist3d_angstrom', 'plddt': 'pLDDT',
    'nearest_pathogenic_pos': 'nearest_pathogenic_residue',
})
summary.to_csv(OUT / "BRCA1_phase1_flagged_candidates.csv", index=False)
print(f"Saved summary table: {OUT / 'BRCA1_phase1_flagged_candidates.csv'} ({len(summary)} rows)")

vus_full_cols = ['variation', 'position', 'wt_aa', 'mt_aa', 'domain', 'plddt',
                  'low_confidence', 'nearest_pathogenic_pos', 'dist3d_A',
                  'seqdist_nearest3d', 'flagged_candidate', 'condition']
vus[vus_full_cols].sort_values('dist3d_A').to_csv(OUT / "BRCA1_all_VUS_annotated.csv", index=False)
print(f"Saved full annotated VUS table: {OUT / 'BRCA1_all_VUS_annotated.csv'} ({len(vus)} rows)")

print()
print("Flagged candidates by domain:")
print(flagged['domain'].value_counts())
print()
print("Flagged candidates: high-confidence vs low-confidence split")
print(flagged['low_confidence'].value_counts())

# --- Task 7: interactive py3Dmol visualization ---
pdb_text = (BASE / "AF-P38398-F1-model_v6.pdb").read_text()
view = py3Dmol.view(width=1000, height=750)
view.addModel(pdb_text, 'pdb')
view.setStyle({}, {'cartoon': {'color': 'lightgrey'}})

pathogenic_positions = sorted(dist_df.loc[dist_df['bucket'] == 'Pathogenic', 'position'].unique().tolist())
flagged_high_conf = flagged[~flagged['low_confidence']]['position'].unique().tolist()
flagged_low_conf = flagged[flagged['low_confidence']]['position'].unique().tolist()

view.addStyle({'resi': [int(p) for p in pathogenic_positions]},
               {'sphere': {'color': 'red', 'radius': 0.9}, 'stick': {'color': 'red'}})
view.addStyle({'resi': [int(p) for p in flagged_high_conf]},
               {'sphere': {'color': 'green', 'radius': 0.9}, 'stick': {'color': 'green'}})
view.addStyle({'resi': [int(p) for p in flagged_low_conf]},
               {'sphere': {'color': 'blue', 'radius': 0.9}, 'stick': {'color': 'blue'}})
view.zoomTo()

html = view._make_html()
legend = """
<div style="font-family: sans-serif; padding: 12px; max-width: 1000px;">
  <h3>BRCA1 &mdash; structural clustering candidates (hypothesis-generation only)</h3>
  <p style="color:#a00; font-weight:bold;">
    This view shows computational/structural clustering ONLY. It is NOT a
    classification or reclassification. Flagged residues are
    "candidates for further investigation," nothing more.
  </p>
  <p>
    Note: BRCA1's AlphaFold model is largely low-confidence outside the
    RING domain (residues 24-65) and the two BRCT domains (1642-1736,
    1756-1855) &mdash; the long central region is intrinsically disordered,
    not just poorly predicted. 74.5% of Pathogenic/VUS residues fall in
    low-confidence regions, far more than TP53. Interpret clustering in
    the central region with real caution.
  </p>
  <ul>
    <li><span style="color:red;">&#9679;</span> Red = known Pathogenic/Likely pathogenic (germline, ClinVar)</li>
    <li><span style="color:green;">&#9679;</span> Green = VUS flagged as candidate (&le;6&Aring; in 3D AND &gt;10 residues away in sequence from a pathogenic residue &mdash; excludes trivial chain-neighbor cases, high-confidence structure region)</li>
    <li><span style="color:blue;">&#9679;</span> Blue = VUS flagged as candidate but in a LOW-CONFIDENCE (pLDDT&lt;70) structure region &mdash; interpret with extra caution</li>
    <li>Grey cartoon = rest of the modeled protein</li>
  </ul>
</div>
"""
(OUT / "BRCA1_structure_viz.html").write_text(legend + html)
print(f"Saved visualization: {OUT / 'BRCA1_structure_viz.html'}")
