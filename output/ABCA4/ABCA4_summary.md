# ABCA4 — structural clustering of VUS against known pathogenic residues

*Draft working notes, not submission-ready.*

## What this is, and isn't

Same framing as TP53/BRCA1/PTEN: **hypothesis-generation only**. Every
flagged variant is a **candidate for further investigation**, not a
verdict.

## Data and methodology

- Source: NCBI ClinVar via E-utilities (germline-only, transcript-matched
  to ABCA4, protein change parsed from the variant title) — same
  filtering methodology as the original three genes, sourced
  programmatically this round.
- 1630 germline missense variants: 659
  Pathogenic/Likely pathogenic, 33 Benign/Likely benign,
  938 VUS.
- Structure: AlphaFold model, 2273 residues
  (matches UniProt canonical
  length 2273).
- Residue numbering cross-checked across ClinVar, UniProt, and
  structure: **0 mismatches** found.
- Flagging threshold: **≤ 6 Å in 3D AND > 10 residues apart in
  sequence** (same dual criterion validated on TP53, to exclude trivial
  chain-adjacent cases).

## Structural confidence

305 of 1597
Pathogenic+VUS residues fall in low-confidence regions (pLDDT<70). Kept
in every table, not excluded. Of the 64 flagged candidates,
2 are low-confidence.

## Results

**64 of 938 VUS (6.8%)**
flagged as candidates.

By domain: {'Unannotated / linker': 42, 'ABC transporter 1': 13, 'ABC transporter 2': 8, 'Disordered': 1}

## Outputs

- `ABCA4_phase1_flagged_candidates.csv`
- `ABCA4_all_VUS_annotated.csv`
- `numbering_mismatches.csv`
- `low_confidence_flagged.csv`
- `ABCA4_structure_viz.html`
