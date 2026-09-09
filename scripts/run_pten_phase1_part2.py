import sys
from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parent.parent  # project root
sys.path.insert(0, str(Path(__file__).parent))
from clustering_lib import extract_fasta, parse_clinvar, parse_structure, compute_distances
from domains import annotate_domain
import numpy as np
import py3Dmol

BASE = (PROJECT_ROOT / "Genes/PTEN")
OUT = (PROJECT_ROOT / "output/PTEN")
GENE = "PTEN"
GENE_TAG = "PTEN"
LOW_PLDDT_THRESHOLD = 70
DIST_THRESHOLD_A = 6.0  # same methodology/threshold as TP53 and BRCA1
MIN_SEQDIST_RESIDUES = 10  # confirmed with user: excludes trivial chain-neighbor "clustering"

# --- rebuild task 1-5 state ---
header, seq = extract_fasta(BASE / "PTEN sequence.rtf")
struct_df = parse_structure(BASE / "AF-P60484-F1-model_v6.pdb")
clinvar_df, dropped = parse_clinvar(BASE / "clinvar_download_20260817_143437.tsv", GENE_TAG)
dist_df = compute_distances(clinvar_df, struct_df)
dist_df = dist_df.merge(struct_df[['position', 'plddt']], on='position', how='left')
dist_df['low_confidence'] = dist_df['plddt'] < LOW_PLDDT_THRESHOLD
dist_df['domain'] = dist_df['position'].apply(lambda p: annotate_domain(GENE, p))

# --- Task 6: flag VUS within distance threshold (Ca OR Cb qualifies; the
# nearest-QUALIFYING logic in compute_distances already excludes trivial
# chain-adjacent cases, and already fixes the nearest-neighbor masking bug
# found via the Cb sensitivity check) ---
vus = dist_df[dist_df['bucket'] == 'VUS'].copy()
vus['flagged_ca'] = vus['dist3d_A'] <= DIST_THRESHOLD_A
vus['flagged_cb'] = vus['dist3d_cb_A'] <= DIST_THRESHOLD_A
vus['flagged_candidate'] = vus['flagged_ca'] | vus['flagged_cb']
vus['confirmed_by_both_atoms'] = vus['flagged_ca'] & vus['flagged_cb']
flagged = vus[vus['flagged_candidate']].copy()
flagged['seq_over_3d'] = np.where(
    flagged['flagged_ca'],
    flagged['seqdist_nearest3d'] / flagged['dist3d_A'].replace(0, 0.01),
    flagged['seqdist_cb'] / flagged['dist3d_cb_A'].replace(0, 0.01),
)
flagged = flagged.sort_values('seq_over_3d', ascending=False)

print(f"=== Task 6: VUS flagged as candidates (Ca<={DIST_THRESHOLD_A}A OR Cb<={DIST_THRESHOLD_A}A, "
      f"each already requiring >{MIN_SEQDIST_RESIDUES} residues apart in sequence) ===")
print(f"{len(flagged)} / {len(vus)} VUS flagged "
      f"({flagged['confirmed_by_both_atoms'].sum()} confirmed by both atom types)")
print(f"  of which low-confidence (pLDDT<{LOW_PLDDT_THRESHOLD}): {flagged['low_confidence'].sum()}")

# --- Task 8: summary table ---
summary_cols = ['variation', 'position', 'wt_aa', 'mt_aa', 'domain', 'plddt',
                 'low_confidence', 'nearest_pathogenic_pos', 'dist3d_A', 'seqdist_nearest3d',
                 'nearest_pathogenic_pos_cb', 'dist3d_cb_A', 'seqdist_cb',
                 'confirmed_by_both_atoms', 'seq_over_3d', 'condition']
summary = flagged[summary_cols].rename(columns={
    'dist3d_A': 'dist3d_ca_angstrom', 'plddt': 'pLDDT',
    'nearest_pathogenic_pos': 'nearest_pathogenic_residue_ca',
    'dist3d_cb_A': 'dist3d_cb_angstrom',
    'nearest_pathogenic_pos_cb': 'nearest_pathogenic_residue_cb',
})
summary.to_csv(OUT / "PTEN_phase1_flagged_candidates.csv", index=False)
print(f"Saved summary table: {OUT / 'PTEN_phase1_flagged_candidates.csv'} ({len(summary)} rows)")

vus_full_cols = ['variation', 'position', 'wt_aa', 'mt_aa', 'domain', 'plddt',
                  'low_confidence', 'nearest_pathogenic_pos', 'dist3d_A', 'seqdist_nearest3d',
                  'nearest_pathogenic_pos_cb', 'dist3d_cb_A', 'seqdist_cb',
                  'flagged_candidate', 'confirmed_by_both_atoms', 'condition']
vus[vus_full_cols].sort_values('dist3d_A').to_csv(OUT / "PTEN_all_VUS_annotated.csv", index=False)
print(f"Saved full annotated VUS table: {OUT / 'PTEN_all_VUS_annotated.csv'} ({len(vus)} rows)")

print()
print("Flagged candidates by domain:")
print(flagged['domain'].value_counts())
print()
print("Flagged candidates: high-confidence vs low-confidence split")
print(flagged['low_confidence'].value_counts())
print()
print("Reminder: only 3 germline Benign variants in the whole PTEN dataset (class imbalance note for Phase 2).")

# --- Task 7: interactive py3Dmol visualization ---
pdb_text = (BASE / "AF-P60484-F1-model_v6.pdb").read_text()
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
  <h3>PTEN &mdash; structural clustering candidates (hypothesis-generation only)</h3>
  <p style="color:#a00; font-weight:bold;">
    This view shows computational/structural clustering ONLY. It is NOT a
    classification or reclassification. Flagged residues are
    "candidates for further investigation," nothing more.
  </p>
  <p>
    Note: PTEN is a small, mostly single-domain folded protein
    (phosphatase + C2 domain, residues ~14-350) with only a short
    disordered C-terminal tail, so a distance-only rule over-flags here.
    Requiring &gt;10 residues of sequence separation (see legend) removes
    trivial chain-neighbor cases and keeps the list to genuine
    tertiary-structure clustering.
  </p>
  <ul>
    <li><span style="color:red;">&#9679;</span> Red = known Pathogenic/Likely pathogenic (germline, ClinVar)</li>
    <li><span style="color:green;">&#9679;</span> Green = VUS flagged as candidate (&le;6&Aring; by backbone Ca <em>or</em> side-chain Cb distance to a pathogenic residue, AND &gt;10 residues away in sequence &mdash; excludes trivial chain-neighbor cases, high-confidence structure region)</li>
    <li><span style="color:blue;">&#9679;</span> Blue = VUS flagged as candidate but in a LOW-CONFIDENCE (pLDDT&lt;70) structure region &mdash; interpret with extra caution</li>
    <li>Grey cartoon = rest of the modeled protein</li>
  </ul>
  <p style="font-size:0.9em; color:#555;">Distance is computed both ways (backbone Ca and side-chain Cb); a candidate confirmed by both is stronger evidence than one found by only one atom type &mdash; see the flagged-candidates CSV for the <code>confirmed_by_both_atoms</code> column.</p>
</div>
"""
(OUT / "PTEN_structure_viz.html").write_text(legend + html)
print(f"Saved visualization: {OUT / 'PTEN_structure_viz.html'}")
