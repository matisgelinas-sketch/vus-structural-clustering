# RB1 — structural clustering of VUS against known pathogenic residues

*Draft working notes, not submission-ready.*

## What this is, and isn't

Same framing as TP53/BRCA1/PTEN: **hypothesis-generation only**. Every
flagged variant is a **candidate for further investigation**, not a
verdict.

## Data and methodology

- Source: NCBI ClinVar via E-utilities (germline-only, transcript-matched
  to RB1, protein change parsed from the variant title) — same
  filtering methodology as the original three genes, sourced
  programmatically this round.
- 1580 germline missense variants: 41
  Pathogenic/Likely pathogenic, 64 Benign/Likely benign,
  1475 VUS.
- Structure: AlphaFold model, 928 residues
  (matches UniProt canonical
  length 928).
- Residue numbering cross-checked across ClinVar, UniProt, and
  structure: **0 mismatches** found.
- Flagging threshold: **≤ 6 Å in 3D AND > 10 residues apart in
  sequence** (same dual criterion validated on TP53, to exclude trivial
  chain-adjacent cases).

## Structural confidence

562 of 1516
Pathogenic+VUS residues fall in low-confidence regions (pLDDT<70). Kept
in every table, not excluded. Of the 22 flagged candidates,
0 are low-confidence.

## Results

**22 of 1475 VUS (1.5%)**
flagged as candidates.

By domain: {'Pocket; binds T and E1A': 20, 'Unannotated / linker': 2}

## Outputs

- `RB1_phase1_flagged_candidates.csv`
- `RB1_all_VUS_annotated.csv`
- `numbering_mismatches.csv`
- `low_confidence_flagged.csv`
- `RB1_structure_viz.html`
