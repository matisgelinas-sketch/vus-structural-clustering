# RET — structural clustering of VUS against known pathogenic residues

*Draft working notes, not submission-ready.*

## What this is, and isn't

Same framing as TP53/BRCA1/PTEN: **hypothesis-generation only**. Every
flagged variant is a **candidate for further investigation**, not a
verdict.

## Data and methodology

- Source: NCBI ClinVar via E-utilities (germline-only, transcript-matched
  to RET, protein change parsed from the variant title) — same
  filtering methodology as the original three genes, sourced
  programmatically this round.
- 1989 germline missense variants: 81
  Pathogenic/Likely pathogenic, 47 Benign/Likely benign,
  1861 VUS.
- Structure: AlphaFold model, 1114 residues
  (matches UniProt canonical
  length 1114).
- Residue numbering cross-checked across ClinVar, UniProt, and
  structure: **0 mismatches** found.
- Flagging threshold: **≤ 6 Å in 3D AND > 10 residues apart in
  sequence** (same dual criterion validated on TP53, to exclude trivial
  chain-adjacent cases).

## Structural confidence

519 of 1942
Pathogenic+VUS residues fall in low-confidence regions (pLDDT<70). Kept
in every table, not excluded. Of the 73 flagged candidates,
0 are low-confidence.

## Results

**73 of 1861 VUS (3.9%)**
flagged as candidates.

By domain: {'Unannotated / linker': 24, 'Protein kinase': 20, 'Cadherin-like region 1 (CLD1)': 13, 'Cadherin': 12, 'Cadherin-like region 4 (CLD4)': 2, 'Cadherin-like region 3 (CLD3)': 2}

## Outputs

- `RET_phase1_flagged_candidates.csv`
- `RET_all_VUS_annotated.csv`
- `numbering_mismatches.csv`
- `low_confidence_flagged.csv`
- `RET_structure_viz.html`
