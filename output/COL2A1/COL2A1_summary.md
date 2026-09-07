# COL2A1 — structural clustering of VUS against known pathogenic residues

*Draft working notes, not submission-ready.*

## What this is, and isn't

Same framing as TP53/BRCA1/PTEN: **hypothesis-generation only**. Every
flagged variant is a **candidate for further investigation**, not a
verdict.

## Data and methodology

- Source: NCBI ClinVar via E-utilities (germline-only, transcript-matched
  to COL2A1, protein change parsed from the variant title) — same
  filtering methodology as the original three genes, sourced
  programmatically this round.
- 1284 germline missense variants: 489
  Pathogenic/Likely pathogenic, 129 Benign/Likely benign,
  666 VUS.
- Structure: AlphaFold model, 1487 residues
  (matches UniProt canonical
  length 1487).
- Residue numbering cross-checked across ClinVar, UniProt, and
  structure: **0 mismatches** found.
- Flagging threshold: **≤ 6 Å in 3D AND > 10 residues apart in
  sequence** (same dual criterion validated on TP53, to exclude trivial
  chain-adjacent cases).

## Structural confidence

975 of 1155
Pathogenic+VUS residues fall in low-confidence regions (pLDDT<70). Kept
in every table, not excluded. Of the 11 flagged candidates,
1 are low-confidence.

## Results

**11 of 666 VUS (1.7%)**
flagged as candidates.

By domain: {'Fibrillar collagen NC1': 6, 'VWFC': 4, 'Disordered': 1}

## Outputs

- `COL2A1_phase1_flagged_candidates.csv`
- `COL2A1_all_VUS_annotated.csv`
- `numbering_mismatches.csv`
- `low_confidence_flagged.csv`
- `COL2A1_structure_viz.html`
