"""
Phase 1, generalized across every auto-sourced gene (originally 12:
VHL, MLH1, MSH2, MSH6, PMS2, RB1, RET, TSC2, LDLR, MYBPC3, HNF1A, ENG --
now whichever have been added since via gene_pipeline.py). Same
methodology as TP53/BRCA1/PTEN: dual-criterion flagging (<=6A in 3D AND
>10 residues apart in sequence), low-confidence (pLDDT<70) flagged not
excluded, domain annotation from UniProt features, interactive py3Dmol
viz, plain-language per-gene summary.

Unlike the original 3 genes, these already have a pre-parsed ClinVar
CSV (from scripts/phase2b/ncbi_clinvar.py) in the same schema
clustering_lib.parse_clinvar() returns, so it's loaded directly rather
than re-parsed from a raw ClinVar TSV.

Gene list is auto-discovered from expansion_log.json (every gene that
fully passed the data bar + structure validation) -- adding a gene via
gene_pipeline.py is enough to have it picked up here, no edits needed.
Genes that already have Phase 1 output are skipped (idempotent/resumable).
"""
import sys
import glob
import json
from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parent.parent  # project root
sys.path.insert(0, str(Path(__file__).parent))
import pandas as pd
from clustering_lib import extract_fasta, parse_structure, cross_check_numbering, compute_distances
from domains import annotate_domain
import py3Dmol

LOG_PATH = (PROJECT_ROOT / "scripts/phase2b/expansion_log.json")
_log = json.loads(LOG_PATH.read_text()) if LOG_PATH.exists() else {}
ALL_VALIDATED = sorted(g for g, e in _log.items() if e.get("structure_passed"))

BASE_DIR = (PROJECT_ROOT / "Genes")
OUT_DIR_CHECK = (PROJECT_ROOT / "output")
GENES = [g for g in ALL_VALIDATED
         if not (OUT_DIR_CHECK / g / f"{g}_phase1_flagged_candidates.csv").exists()]
skipped = [g for g in ALL_VALIDATED if g not in GENES]
if skipped:
    print(f"Skipping {len(skipped)} genes with existing Phase 1 output: {', '.join(skipped)}")
print(f"Running Phase 1 on {len(GENES)} genes: {', '.join(GENES) if GENES else '(none -- all up to date)'}\n")
OUT_DIR = (PROJECT_ROOT / "output")

LOW_PLDDT_THRESHOLD = 70
DIST_THRESHOLD_A = 6.0
MIN_SEQDIST_RESIDUES = 10

results_summary = []

for gene in GENES:
    gdir = BASE_DIR / gene
    odir = OUT_DIR / gene
    odir.mkdir(parents=True, exist_ok=True)

    fasta_path = glob.glob(str(gdir / "*_sequence.fasta"))[0]
    pdb_path = glob.glob(str(gdir / "*.pdb"))[0]
    csv_path = glob.glob(str(gdir / "*_clinvar_germline_missense.csv"))[0]

    header, seq = extract_fasta(fasta_path)
    struct_df = parse_structure(pdb_path)
    clinvar_df = pd.read_csv(csv_path)

    length_ok = len(struct_df) == len(seq)

    mismatches = cross_check_numbering(clinvar_df, seq, struct_df)
    mismatches.to_csv(odir / "numbering_mismatches.csv", index=False)

    merged = clinvar_df.merge(struct_df[['position', 'plddt']], on='position', how='left')
    merged['low_confidence'] = merged['plddt'] < LOW_PLDDT_THRESHOLD
    merged['domain'] = merged['position'].apply(lambda p: annotate_domain(gene, p))

    dist_df = compute_distances(clinvar_df, struct_df)
    dist_df = dist_df.merge(struct_df[['position', 'plddt']], on='position', how='left')
    dist_df['low_confidence'] = dist_df['plddt'] < LOW_PLDDT_THRESHOLD
    dist_df['domain'] = dist_df['position'].apply(lambda p: annotate_domain(gene, p))

    vus = dist_df[dist_df['bucket'] == 'VUS'].dropna(subset=['dist3d_A']).copy()
    vus['flagged_candidate'] = (vus['dist3d_A'] <= DIST_THRESHOLD_A) & (vus['seqdist_nearest3d'] > MIN_SEQDIST_RESIDUES)
    flagged = vus[vus['flagged_candidate']].copy()
    flagged['seq_over_3d'] = flagged['seqdist_nearest3d'] / flagged['dist3d_A'].replace(0, 0.01)
    flagged = flagged.sort_values('seq_over_3d', ascending=False)

    n_patho = (clinvar_df['bucket'] == 'Pathogenic').sum()
    n_benign = (clinvar_df['bucket'] == 'Benign').sum()
    n_vus = (clinvar_df['bucket'] == 'VUS').sum()
    low_conf_pb = merged[(merged['low_confidence']) & (merged['bucket'].isin(['Pathogenic', 'VUS']))]

    summary_cols = ['variation', 'position', 'wt_aa', 'mt_aa', 'domain', 'plddt',
                     'low_confidence', 'nearest_pathogenic_pos', 'dist3d_A',
                     'seqdist_nearest3d', 'seq_over_3d', 'condition']
    flagged[summary_cols].rename(columns={
        'dist3d_A': 'dist3d_angstrom', 'plddt': 'pLDDT',
        'nearest_pathogenic_pos': 'nearest_pathogenic_residue',
    }).to_csv(odir / f"{gene}_phase1_flagged_candidates.csv", index=False)

    vus_full_cols = ['variation', 'position', 'wt_aa', 'mt_aa', 'domain', 'plddt',
                      'low_confidence', 'nearest_pathogenic_pos', 'dist3d_A',
                      'seqdist_nearest3d', 'flagged_candidate', 'condition']
    vus[vus_full_cols].sort_values('dist3d_A').to_csv(odir / f"{gene}_all_VUS_annotated.csv", index=False)

    low_conf_pb.to_csv(odir / "low_confidence_flagged.csv", index=False)

    # --- visualization ---
    pdb_text = Path(pdb_path).read_text()
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
    legend = f"""
<div style="font-family: sans-serif; padding: 12px; max-width: 1000px;">
  <h3>{gene} &mdash; structural clustering candidates (hypothesis-generation only)</h3>
  <p style="color:#a00; font-weight:bold;">
    This view shows computational/structural clustering ONLY. It is NOT a
    classification or reclassification. Flagged residues are
    "candidates for further investigation," nothing more.
  </p>
  <ul>
    <li><span style="color:red;">&#9679;</span> Red = known Pathogenic/Likely pathogenic (germline, ClinVar)</li>
    <li><span style="color:green;">&#9679;</span> Green = VUS flagged as candidate (&le;6&Aring; in 3D AND &gt;10 residues away in sequence, high-confidence structure region)</li>
    <li><span style="color:blue;">&#9679;</span> Blue = VUS flagged as candidate but in a LOW-CONFIDENCE (pLDDT&lt;70) structure region</li>
    <li>Grey cartoon = rest of the modeled protein</li>
  </ul>
</div>
"""
    (odir / f"{gene}_structure_viz.html").write_text(legend + html)

    domain_counts = flagged['domain'].value_counts().to_dict()
    n_low_conf_flagged = int(flagged['low_confidence'].sum())

    # --- summary.md ---
    summary_md = f"""# {gene} — structural clustering of VUS against known pathogenic residues

*Draft working notes, not submission-ready.*

## What this is, and isn't

Same framing as TP53/BRCA1/PTEN: **hypothesis-generation only**. Every
flagged variant is a **candidate for further investigation**, not a
verdict.

## Data and methodology

- Source: NCBI ClinVar via E-utilities (germline-only, transcript-matched
  to {gene}, protein change parsed from the variant title) — same
  filtering methodology as the original three genes, sourced
  programmatically this round.
- {n_patho + n_benign + n_vus} germline missense variants: {n_patho}
  Pathogenic/Likely pathogenic, {n_benign} Benign/Likely benign,
  {n_vus} VUS.
- Structure: AlphaFold model, {len(struct_df)} residues
  ({'matches' if length_ok else 'DOES NOT MATCH'} UniProt canonical
  length {len(seq)}).
- Residue numbering cross-checked across ClinVar, UniProt, and
  structure: **{len(mismatches)} mismatches** found.
- Flagging threshold: **≤ 6 Å in 3D AND > 10 residues apart in
  sequence** (same dual criterion validated on TP53, to exclude trivial
  chain-adjacent cases).

## Structural confidence

{len(low_conf_pb)} of {len(merged[merged['bucket'].isin(['Pathogenic','VUS'])])}
Pathogenic+VUS residues fall in low-confidence regions (pLDDT<70). Kept
in every table, not excluded. Of the {len(flagged)} flagged candidates,
{n_low_conf_flagged} are low-confidence.

## Results

**{len(flagged)} of {len(vus)} VUS ({100*len(flagged)/max(len(vus),1):.1f}%)**
flagged as candidates.

By domain: {domain_counts}

## Outputs

- `{gene}_phase1_flagged_candidates.csv`
- `{gene}_all_VUS_annotated.csv`
- `numbering_mismatches.csv`
- `low_confidence_flagged.csv`
- `{gene}_structure_viz.html`
"""
    (odir / f"{gene}_summary.md").write_text(summary_md)

    results_summary.append(dict(
        gene=gene, length=len(seq), length_ok=length_ok, mismatches=len(mismatches),
        n_patho=n_patho, n_benign=n_benign, n_vus=n_vus,
        n_flagged=len(flagged), pct_flagged=100*len(flagged)/max(len(vus),1),
        n_low_conf_flagged=n_low_conf_flagged,
    ))
    print(f"{gene}: {len(struct_df)}res, {len(mismatches)} mismatches, "
          f"{n_patho}P/{n_benign}B/{n_vus}V, flagged={len(flagged)} ({100*len(flagged)/max(len(vus),1):.1f}%), "
          f"low_conf_flagged={n_low_conf_flagged}")

summary_path = OUT_DIR / "new_genes_phase1_summary.csv"
new_df = pd.DataFrame(results_summary)
if summary_path.exists() and len(new_df):
    prior = pd.read_csv(summary_path)
    summary_df = pd.concat([prior[~prior['gene'].isin(new_df['gene'])], new_df], ignore_index=True)
elif summary_path.exists():
    summary_df = pd.read_csv(summary_path)
else:
    summary_df = new_df
summary_df = summary_df.sort_values('gene').reset_index(drop=True)
print(f"\n=== Summary across all {len(summary_df)} auto-sourced genes ===")
print(summary_df.to_string(index=False))
summary_df.to_csv(summary_path, index=False)
