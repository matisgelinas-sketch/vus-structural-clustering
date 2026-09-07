# COL5A1 — structural clustering of VUS against known pathogenic residues

*Draft working notes, not submission-ready.*

## What this is, and isn't

Same framing as TP53/BRCA1/PTEN: **hypothesis-generation only**. Every
flagged variant is a **candidate for further investigation**, not a
verdict.

## Data and methodology

- Source: NCBI ClinVar via E-utilities (germline-only, transcript-matched
  to COL5A1, protein change parsed from the variant title) — same
  filtering methodology as the original three genes, sourced
  programmatically this round.
- 1160 germline missense variants: 41
  Pathogenic/Likely pathogenic, 205 Benign/Likely benign,
  914 VUS.
- Structure: AlphaFold model, 1838 residues
  (matches UniProt canonical
  length 1838).
- Residue numbering cross-checked across ClinVar, UniProt, and
  structure: **0 mismatches** found.
- Flagging threshold: **≤ 6 Å in 3D AND > 10 residues apart in
  sequence** (same dual criterion validated on TP53, to exclude trivial
  chain-adjacent cases).

## Structural confidence

710 of 955
Pathogenic+VUS residues fall in low-confidence regions (pLDDT<70). Kept
in every table, not excluded. Of the 22 flagged candidates,
2 are low-confidence.

## Results

**22 of 914 VUS (2.4%)**
flagged as candidates.

By domain: {'Fibrillar collagen NC1': 12, 'Laminin G-like': 8, 'Disordered': 2}

## Outputs

- `COL5A1_phase1_flagged_candidates.csv`
- `COL5A1_all_VUS_annotated.csv`
- `numbering_mismatches.csv`
- `low_confidence_flagged.csv`
- `COL5A1_structure_viz.html`
