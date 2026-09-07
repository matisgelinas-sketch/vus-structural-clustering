# MSH2 — structural clustering of VUS against known pathogenic residues

*Draft working notes, not submission-ready.*

## What this is, and isn't

Same framing as TP53/BRCA1/PTEN: **hypothesis-generation only**. Every
flagged variant is a **candidate for further investigation**, not a
verdict.

## Data and methodology

- Source: NCBI ClinVar via E-utilities (germline-only, transcript-matched
  to MSH2, protein change parsed from the variant title) — same
  filtering methodology as the original three genes, sourced
  programmatically this round.
- 2181 germline missense variants: 153
  Pathogenic/Likely pathogenic, 304 Benign/Likely benign,
  1724 VUS.
- Structure: AlphaFold model, 934 residues
  (matches UniProt canonical
  length 934).
- Residue numbering cross-checked across ClinVar, UniProt, and
  structure: **0 mismatches** found.
- Flagging threshold: **≤ 6 Å in 3D AND > 10 residues apart in
  sequence** (same dual criterion validated on TP53, to exclude trivial
  chain-adjacent cases).

## Structural confidence

197 of 1877
Pathogenic+VUS residues fall in low-confidence regions (pLDDT<70). Kept
in every table, not excluded. Of the 90 flagged candidates,
0 are low-confidence.

## Results

**90 of 1724 VUS (5.2%)**
flagged as candidates.

By domain: {'Unannotated / linker': 83, 'Interaction with EXO1': 7}

## Outputs

- `MSH2_phase1_flagged_candidates.csv`
- `MSH2_all_VUS_annotated.csv`
- `numbering_mismatches.csv`
- `low_confidence_flagged.csv`
- `MSH2_structure_viz.html`
