# SCN1A — structural clustering of VUS against known pathogenic residues

*Draft working notes, not submission-ready.*

## What this is, and isn't

Same framing as TP53/BRCA1/PTEN: **hypothesis-generation only**. Every
flagged variant is a **candidate for further investigation**, not a
verdict.

## Data and methodology

- Source: NCBI ClinVar via E-utilities (germline-only, transcript-matched
  to SCN1A, protein change parsed from the variant title) — same
  filtering methodology as the original three genes, sourced
  programmatically this round.
- 2420 germline missense variants: 1037
  Pathogenic/Likely pathogenic, 63 Benign/Likely benign,
  1320 VUS.
- Structure: AlphaFold model, 2009 residues
  (matches UniProt canonical
  length 2009).
- Residue numbering cross-checked across ClinVar, UniProt, and
  structure: **0 mismatches** found.
- Flagging threshold: **≤ 6 Å in 3D AND > 10 residues apart in
  sequence** (same dual criterion validated on TP53, to exclude trivial
  chain-adjacent cases).

## Structural confidence

642 of 2357
Pathogenic+VUS residues fall in low-confidence regions (pLDDT<70). Kept
in every table, not excluded. Of the 25 flagged candidates,
2 are low-confidence.

## Results

**25 of 1320 VUS (1.9%)**
flagged as candidates.

By domain: {'IV': 8, 'III': 8, 'II': 7, 'Unannotated / linker': 2}

## Outputs

- `SCN1A_phase1_flagged_candidates.csv`
- `SCN1A_all_VUS_annotated.csv`
- `numbering_mismatches.csv`
- `low_confidence_flagged.csv`
- `SCN1A_structure_viz.html`
