# VHL — structural clustering of VUS against known pathogenic residues

*Draft working notes, not submission-ready.*

## What this is, and isn't

Same framing as TP53/BRCA1/PTEN: **hypothesis-generation only**. Every
flagged variant is a **candidate for further investigation**, not a
verdict.

## Data and methodology

- Source: NCBI ClinVar via E-utilities (germline-only, transcript-matched
  to VHL, protein change parsed from the variant title) — same
  filtering methodology as the original three genes, sourced
  programmatically this round.
- 597 germline missense variants: 180
  Pathogenic/Likely pathogenic, 35 Benign/Likely benign,
  382 VUS.
- Structure: AlphaFold model, 213 residues
  (matches UniProt canonical
  length 213).
- Residue numbering cross-checked across ClinVar, UniProt, and
  structure: **0 mismatches** found.
- Flagging threshold: **≤ 6 Å in 3D AND > 10 residues apart in
  sequence** (same dual criterion validated on TP53, to exclude trivial
  chain-adjacent cases).

## Structural confidence

166 of 562
Pathogenic+VUS residues fall in low-confidence regions (pLDDT<70). Kept
in every table, not excluded. Of the 11 flagged candidates,
0 are low-confidence.

## Results

**11 of 382 VUS (2.9%)**
flagged as candidates.

By domain: {'Involved in binding to CCT complex': 11}

## Outputs

- `VHL_phase1_flagged_candidates.csv`
- `VHL_all_VUS_annotated.csv`
- `numbering_mismatches.csv`
- `low_confidence_flagged.csv`
- `VHL_structure_viz.html`
